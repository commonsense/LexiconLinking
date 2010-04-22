/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This file contains routines which handle the pivot rejection list. In
  addition to routines to add and remove entries, there's a routine to
  examine the list and decide on an appropriate action when a simplex punts.

  The basic mechanism used for pivot rejection is a qualifier flag,
  vstatNOPIVOT, which is added to a variable's status.  Flagged variables
  are not considered for entering the basis (primal) or leaving the basis
  (dual).

  For management purposes, information about flagged variables is kept in an
  array, pivrejlst. There is also a control structure, pivrej_ctl.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_pivreject.c	4.4	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_pivreject.c 269 2009-04-02 05:38:19Z lou $" ;



/*
  pivrej_struct
  
  Each pivrej_struct holds the index of the variable, the value of basis.pivs
  at the time of rejection, and an indication of why the variable was
  rejected:
    * dyrSINGULAR indicates that the basis became singular when this
      variable tried to enter the basis.
    * dyrMADPIV indicates that abar<ij> was rejected by dy_chkpiv; in this
      case the rejection ratio is stored.

  Keeping a list saves scanning the entire nonbasic (primal) or
  basic (dual) partition when we come to clear marked variables after a
  successful pivot. Keeping the acceptance ratio allows us to make an
  intelligent decision about whether to reduce the pivot tolerance or break
  out of simplex and attempt to modify the constraint system.

  Field		Definition
  -----		----------
  ndx		index of rejected variable
  iter		value of basis.pivs when the pivot was rejected
  why		the reason the pivot was rejected (dyrSINGULAR or dyrMADPIV)
  ratio         the ratio returned by dy_chkpiv when the pivot was rejected,
		multiplied by the value of dy_tols.pivot at time of rejection
		(for independence from possible changes to tols.pivot)
*/

typedef struct
{ int ndx ;
  int iter ;
  dyret_enum why ;
  double ratio ; } pivrej_struct ;

static pivrej_struct *pivrejlst = NULL ;

/*
  pivrejctl_struct

  Control structure for the pivot rejection mechanism.

  Field		Definition
  -----		----------
  sze		allocated capacity of pivrejlst
  cnt		number of entries in pivrejlst
  mad		number of entries rejected due to small abar<ij>
  sing		number of entries rejected due to singular basis
  iter_reduced	value of basis.pivs when the pivot tolerance was reduced;
		-1 if we're running with the default tolerance.
  savedtol	saved copy of the default pivot tolerance; captured by
		initpivrej
  pivmul	reduction factor for pivot tolerance; currently hardcoded in
		initpivrej; when reduction is required, it's in steps of
		1/pivmul.
*/

typedef struct
{ int sze ;
  int cnt ;
  int mad ;
  int sing ;
  int iter_reduced ;
  double savedtol ;
  double pivmul ; } pivrejctl_struct ;

static pivrejctl_struct pivrej_ctl ;


void dy_initpivrej (int sze)
/*
  Allocate the pivrej control structure and an initial pivrejlst array.

  Parameter:
    sze:	allocated capacity of the pivot rejection list
*/
{ pivrej_ctl.sze = maxx(sze,5) ;
  pivrejlst =
    (pivrej_struct *) MALLOC(pivrej_ctl.sze*sizeof(pivrej_struct)) ;
  pivrej_ctl.cnt = 0 ;
  pivrej_ctl.mad = 0 ;
  pivrej_ctl.sing = 0 ;
  pivrej_ctl.iter_reduced = -1 ;
  pivrej_ctl.savedtol = dy_tols->pivot ;
  pivrej_ctl.pivmul = 1000.0 ;

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->pivrej.min_pivtol = dy_tols->pivot ;
# endif

  return ; }

void dy_freepivrej (void)
/*
  Free the pivot rejection list. Solely for information hiding.
*/

{ if (pivrejlst != NULL)
  { FREE(pivrejlst) ;
    pivrejlst = NULL ; }
  
  return ; }



void dy_checkpivtol (void)
/*
  A quick little routine to see if we've been running long enough on a reduced
  pivot multiplier. If so, get back to the default. The primary purpose of this
  is information hiding.
*/

{ if (pivrej_ctl.iter_reduced > 0 &&
      dy_lp->basis.pivs-pivrej_ctl.iter_reduced > dy_opts->factor)
  { dy_tols->pivot = pivrej_ctl.savedtol ;
    pivrej_ctl.iter_reduced = -1 ; }
    
  return ; }




static int int_nonincreasing (const void *p_i, const void *p_j)
/*
  Reverse integer comparison so we can sort arrays of indices in nonincreasing
  order.

  Returns:	< 0	if i > j
		  0	if i = j
		> 0	if i < j
*/
{ int i = *((const int *) p_i) ;
  int j = *((const int *) p_j) ;

  return ((j)-(i)) ; }


bool dy_clrpivrej (int *entries)

/*
  This routine removes variables from rejected pivot list. As far as the rest
  of dylp is concerned, all that needs to be done is to clear the NOPIVOT
  qualifier from the variable's status entry.

  Internally, there are really two modes: clear specified entries, and clear
  the entire list. If the client supplies an array of indices in the entries
  parameter, selective removal is performed, otherwise the entire list is
  cleared.

  Parameters:
    entries:	an array with indices of entries to be removed
		entries[0] is expected to contain the number of entries

  Returns: TRUE if the clearing operation is successful, FALSE otherwise.
*/

{ int n,j,ndx,last,endx,elast ;
  flags statj ;

  const char *rtnnme = "dy_clrpivrej" ;

# ifdef DYLP_PARANOIA
  flags chkflgs ;

/*
  For dual simplex, only out-of-bound basic variables are considered for
  pivoting, but subsequent dual pivots could change that to pretty much
  any basic status. For primal simplex, any nonbasic status is ok, including
  the exotic ones.
*/
  if (dy_lp->phase == dyDUAL)
  { chkflgs = vstatBASIC ; }
  else
  { chkflgs = vstatNONBASIC|vstatEXOTIC ; }
# endif

/*
  Are we clearing the entire list? If so, also restore the default pivot
  tolerance. If we're being selective about clearing, leave the tolerance
  unchanged and assume the client will take care of it. If there are no
  entries in pivrejlst, that's all we need to do.
*/
  if (entries == NULL)
  { dy_tols->pivot = pivrej_ctl.savedtol ;
    pivrej_ctl.iter_reduced = -1 ; }
  if (pivrej_ctl.cnt == 0) return (TRUE) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivreject >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    %s pivot reject list ... ",
	        (entries == NULL)?"clearing":"winnowing") ; }
# endif

  n = dy_sys->varcnt ;
  last = pivrej_ctl.cnt-1 ;

/*
  If the client hasn't supplied entries, we're clearing the entire list.
*/
  if (entries == NULL)
  { 
    for (ndx = 0 ; ndx <= last ; ndx++)
    { j = pivrejlst[ndx].ndx ;
      statj = dy_status[j] ;
#     ifdef DYLP_PARANOIA
      if (j < 1 || j > n)
      { errmsg(102,rtnnme,dy_sys->nme,"rejected variable",j,1,n) ;
	return (FALSE) ; }
      if (flgoff(statj,vstatNOPIVOT) || flgoff(statj,chkflgs))
      { errmsg(329,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       consys_nme(dy_sys,'v',j,FALSE,NULL),j,ndx,dy_prtvstat(statj),
	       (dy_lp->phase == dyDUAL)?"basic":"nonbasic") ;
	return (FALSE) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivreject >= 2)
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\trestoring %s (%d) as eligible for pivoting.",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
#     endif
      clrflg(dy_status[j],vstatNOPIVOT) ; }

    last = -1 ;
    pivrej_ctl.mad = 0 ;
    pivrej_ctl.sing = 0 ; }
/*
  The more complicated case: Remove the set of entries specified by the
  client. The sort is necessary so that we can compress in place, moving the
  last entry to replace the deleted entry.
*/
  else
  { elast = entries[0] ;
    if (elast > 1)
    { qsort(&entries[1],elast,sizeof(int),int_nonincreasing) ; }
    for (endx = 1 ; endx <= elast ; endx++)
    { ndx = entries[endx] ;
#     ifdef DYLP_PARANOIA
      if (ndx < 0 || ndx >= pivrej_ctl.cnt)
      { errmsg(102,rtnnme,dy_sys->nme,"pivrej list index",ndx,
	       0,pivrej_ctl.cnt-1) ;
	return (FALSE) ; }
#     endif
      j = pivrejlst[ndx].ndx ;
      statj = dy_status[j] ;
#     ifdef DYLP_PARANOIA
      if (j < 1 || j > n)
      { errmsg(102,rtnnme,dy_sys->nme,"rejected variable",j,1,n) ;
	return (FALSE) ; }
      if (flgoff(statj,vstatNOPIVOT) || flgoff(statj,chkflgs))
      { errmsg(329,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       consys_nme(dy_sys,'v',j,FALSE,NULL),j,ndx,dy_prtvstat(statj),
	       (dy_lp->phase == dyDUAL)?"basic":"nonbasic") ;
	return (FALSE) ; }
#     endif
      clrflg(dy_status[j],vstatNOPIVOT) ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivreject >= 2)
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\trestoring %s (%d) as eligible for pivoting.",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
#     endif
      if (ndx < last)
      { pivrejlst[ndx] = pivrejlst[last] ;
	switch (pivrejlst[ndx].why)
	{ case dyrSINGULAR:
	  { pivrej_ctl.sing-- ;
	    break ; }
	  case dyrMADPIV:
	  { pivrej_ctl.mad-- ;
	    break ; }
	  default:
	  { errmsg(1,rtnnme,__LINE__) ;
	    return (FALSE) ; } } }
      last-- ; } }

    last++ ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivreject >= 1)
    { if (dy_opts->print.pivreject >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,"restored %d variables.",
		  pivrej_ctl.cnt-last) ; }
#   endif
    pivrej_ctl.cnt = last ;

  return (TRUE) ; }



dyret_enum dy_addtopivrej (int j, dyret_enum why,
			   double abarij, double maxabarij)

/*
  This routine adds x<j> to the rejected pivot list by adding an entry to
  pivrejlst and adding the NOPIVOT qualifier to x<j>'s status.
  If necessary, it expands the size of the list.

  Parameter:
    j:		the variable x<j> 
    why:	the reason it's going on the pivot reject list; one of
		dyrSINGULAR or dyrMADPIV
    abarij:	(why == dyrMADPIV) the pivot element
    maxabarij:	(why == dyrMADPIV) the maximum pivot element in the pivot
		column (primal) or row (dual).
  
  Returns: dyrOK if the entry is added without error, dyrFATAL if we can't
	   get more space, or if a paranoid check fails.
*/

{ int n,ndx,newsze ;
  double ratio ;
  const char *rtnnme = "dy_addtopivrej" ;

# ifndef DYLP_NDEBUG
  int saveprint ;

  saveprint = dy_opts->print.pivoting ;
  dy_opts->print.pivoting = 0 ;
# endif

/*
  We don't actually need the pivot ratio until further down, but it's handy
  to do it here where we can easily suppress the internal print, then restore
  the print level.
*/
  ratio = dy_chkpiv(abarij,maxabarij) ;
  n = dy_sys->varcnt ;

# ifndef DYLP_NDEBUG
  dy_opts->print.pivoting = saveprint ;
# endif
# ifdef DYLP_PARANOIA
  if (j < 1 || j > n)
  { errmsg(102,rtnnme,dy_sys->nme,"variable",j,1,n) ;
    return (dyrFATAL) ; }
  if (!(why == dyrSINGULAR || why == dyrMADPIV))
  { errmsg(1,rtnnme,__LINE__) ;
    return (dyrFATAL) ; }
# endif
# ifndef DYLP_NDEBUG
/*
  The default case in this switch is needed to suppress GCC warnings --- it
  doesn't grok the paranoid check.
*/
  if (dy_opts->print.pivreject >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  marking %s (%d) ineligible for pivoting ",
	        consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
    switch (why)
    { case dyrSINGULAR:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"(%s).",dy_prtdyret(why)) ;
	break ; }
      case dyrMADPIV:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"(%s = %g).",dy_prtdyret(why),ratio) ;
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (dyrFATAL) ; } } }
# endif

/*
  Flag the culprit --- the extent of externally visible activity.  Then make
  the entry in the pivot reject list. Check for adequate list length and
  expand if necessary.
*/
  setflg(dy_status[j],vstatNOPIVOT) ;
  ndx = pivrej_ctl.cnt++ ;
  if (ndx >= pivrej_ctl.sze)
  { newsze = minn(2*pivrej_ctl.sze,n+1) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivreject >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: expanding pivot reject list from %d to %d entries.",
		  rtnnme,pivrej_ctl.sze,newsze) ; }
#   endif
    pivrejlst =
      (pivrej_struct *) REALLOC(pivrejlst,newsze*sizeof(pivrej_struct)) ;
    if (pivrejlst == NULL)
    { errmsg(337,rtnnme,dy_sys->nme,pivrej_ctl.sze,newsze) ;
      return (dyrFATAL) ; }
    pivrej_ctl.sze = newsze ; }
  pivrejlst[ndx].ndx = j ;
  pivrejlst[ndx].iter = dy_lp->basis.pivs ;
  pivrejlst[ndx].why = why ;
  switch (why)
  { case dyrSINGULAR:
    { pivrej_ctl.sing++ ;
      break ; }
    case dyrMADPIV:
    { pivrej_ctl.mad++ ;
      ratio = dy_chkpiv(abarij,maxabarij) ;
      pivrejlst[ndx].ratio = ratio*dy_tols->pivot ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (dyrFATAL) ; } }

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL)
  { switch (why)
    { case dyrSINGULAR:
      { dy_stats->pivrej.sing++ ;
	break ; }
      case dyrMADPIV:
      { dy_stats->pivrej.mad++ ;
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (dyrFATAL) ; } }
    if (pivrej_ctl.cnt > dy_stats->pivrej.max)
    { dy_stats->pivrej.max = pivrej_ctl.cnt ; } }
# endif

  return (dyrOK) ; }



dyret_enum dy_dealWithPunt (void)

/*
  This routine decides on the appropriate action(s) when a simplex decides to
  punt. The algorithm is this:
    1) Sort the entries in pivrejlst into two sets: iter == basis.pivs
       (current) and iter != basis.pivs (old). In the current set, count
       the number of mad and singular entries.
    2) If there are any entries in old, remove them from pivrejlst and
       return with an indication to resume pivoting (dyrRESELECT).
    3) If all entries in current are of type singular, return with an
       indication to abort this simplex phase (dyrPUNT) and hope that we can
       alter the constraint system.
    4) For each permissible reduction in pivot tolerance, check for entries
       of type MADPIV that might become acceptable. If there are any, remove
       them from pivrejlst and return dyrRESELECT.
    5) If 4) failed to identify pivots, return dyrPUNT.

  Parameters: none

  Returns: dyrRESELECT if pivoting can resume
	   dyrPUNT to abort this simplex phase
	   dyrFATAL if something goes wrong
*/

{ int j,ndx,last,oldcnt,curcnt,curmad,brk ;
  double maxratio,pivmul ;
  bool clr_retval ;
  dyret_enum retval ;

  int *old,*current ;
  pivrej_struct *pivrej ;

# ifndef DYLP_NDEBUG
  const char *rtnnme = "dy_dealWithPunt" ;
# endif

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->pivrej.puntcall++ ;
# endif

  retval = dyrINV ;
/*
  If there are no rejected pivots, the punt stands.
*/
  if (pivrej_ctl.cnt == 0)
  {
#   ifdef DYLP_STATISTICS
    if (dy_stats != NULL) dy_stats->pivrej.puntret++ ;
#   endif
    return (dyrPUNT) ; }
/*
  Setup and scan pivrejlst as indicated above.
*/
  last = pivrej_ctl.cnt ;
  brk = dy_lp->basis.pivs ;
  old = (int *) MALLOC((last+1)*sizeof(int)) ;
  current = (int *) MALLOC((last+1)*sizeof(int)) ;
  oldcnt = 0 ;
  curcnt = 0 ;
  curmad = 0 ;
  maxratio = 0 ;

  for (ndx = 0 ; ndx < last ; ndx++)
  { pivrej = &pivrejlst[ndx] ;
    if (pivrej->iter != brk)
    { old[++oldcnt] = ndx ; }
    else
    { current[++curcnt] = ndx ;
      if (pivrej->why == dyrMADPIV)
      { curmad++ ;
	if (maxratio < pivrej->ratio) maxratio = pivrej->ratio ; } } }
/*
  If there are old entries, we can always hope the intervening pivots have
  cured the problem. It happens.
*/
  if (oldcnt > 0)
  { old[0] = oldcnt ;
    clr_retval = dy_clrpivrej(old) ;
    if (clr_retval == TRUE)
    { retval = dyrRESELECT ; }
    else
    { retval = dyrFATAL ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivreject >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  restored %d entries queued before iter = %d.",
		  old[0],brk) ; }
#   endif
  }
/*
  Are there any mad pivots that we can press into service by reducing the pivot
  tolerance?
*/
  else
  if (curmad > 0 && maxratio > dy_tols->zero)
  { pivmul = 1/dy_tols->pivot ;
    while (maxratio*pivmul < 1.0) pivmul *= pivrej_ctl.pivmul ;
    if (1/pivmul >= dy_tols->zero*100)
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivreject >= 1)
      { warn(376,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     dy_tols->pivot,1/pivmul) ; }
#     endif
      dy_tols->pivot = 1/pivmul ;
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL)
      { dy_stats->pivrej.pivtol_red++ ;
	if (dy_tols->pivot < dy_stats->pivrej.min_pivtol)
	{ dy_stats->pivrej.min_pivtol = dy_tols->pivot ; } }
#     endif
      j = 0 ;
      for (ndx = 1 ; ndx <= curcnt ; ndx++)
      { pivrej = &pivrejlst[current[ndx]] ;
	if (pivrej->ratio*pivmul > 1.0)
	{ current[++j] = current[ndx] ; } }
      current[0] = j ;
      clr_retval = dy_clrpivrej(current) ;
      if (clr_retval == TRUE)
      { retval = dyrRESELECT ; }
      else
      { retval = dyrFATAL ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivreject >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n  restored %d entries queued at iter = %d at piv. tol = %g",
	      current[0],brk,dy_tols->pivot) ; }
#   endif
    }
    else
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivreject >= 1)
      { warn(383,rtnnme,dy_sys->nme,
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     dy_tols->zero,dy_prtdyret(dyrPUNT)) ; }
#     endif
      retval = dyrPUNT ; } }
  else
  { retval = dyrPUNT ; }
/*
  That's it, we've done our best. Free the old and current arrays and return.
*/
  FREE(old) ;
  FREE(current) ;

# ifndef DYLP_NDEBUG
  if (retval == dyrPUNT && dy_opts->print.pivreject >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  PUNT! mad = %d, singular = %d.",
	        pivrej_ctl.mad,pivrej_ctl.sing) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL && retval == dyrPUNT) dy_stats->pivrej.puntret++ ;
# endif

  return (retval) ; }

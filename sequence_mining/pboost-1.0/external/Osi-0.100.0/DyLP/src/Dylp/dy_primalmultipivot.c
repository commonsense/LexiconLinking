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
  This file contains routines which select the leaving variable under
  somewhat more relaxed primal pivoting rules. The general notion is that we
  can have `soft' and `hard' limits on delta<j> for the incoming variable
  x<j>. To give an example, suppose we're in primal phase I, and we're
  looking at a variable x<i>, BUUB, which will decrease and leave the basis.
  We can pivot x<i> out when it reaches its upper bound (the `soft' limit on
  delta<j>) or allow it to continue to move and pivot it out at the lower
  bound (the `hard' limit).

  The strategy is to scan the candidates to leave, recording each pivot
  opportunity and marking it as soft or hard. Then we sort by nondecreasing
  delta and pick the best looking pivot opportunity from the candidates up to
  the first hard limit.

  The real advantage comes in cases where the candidate with the limiting
  hard delta has a numerically unstable pivot. Often we can promote a
  candidate with a slightly larger delta and a stable pivot, because (tiny
  unstable pivot)*(slightly larger delta) is still too small to cause loss of
  primal feasibility.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_primalmultipivot.c	1.4	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_primalmultipivot.c 269 2009-04-02 05:38:19Z lou $" ;

/*
  The structure to hold pivot candidates

  Field		Description
  -----		-----------
  ndx		The index of the candidate x<k>
  deltakj	Absolute value of delta for this candidate
  abarkj	Absolute value of pivot for this candidate
  ratiokj	Stability ratio for this candidate; 1.0 is the boundary
		between stable and unstable.
  madpiv	TRUE if this pivot is unstable, FALSE otherwise
  dir		Direction of motion of x<k>; -1 to decrease, +1 to increase
  bnd		The bound where x<k> will pivot; -1 for lb, +1 for ub
  hard		TRUE if this is a hard limit, false if it's a soft limit
*/

typedef struct { int ndx ;
		 double deltakj ;
		 double abarkj ;
		 double ratiokj ;
		 bool madpiv ;
		 int dir ;
		 int bnd ;
		 bool hard ; } primcand_struct ;



static int primcand_cmp (const void *p_primcand1, const void *p_primcand2)
/*
  Comparison function for qsort to sort an array of primcand_structs.
  The primary criterion is deltakj, with ratiokj used as the tiebreaker.

  Parameters:
    p_primcand1,2: primcand_structs to be compared

  Returns: -1 if primcand1 < primcand2
	    0 for equality
	    1 if primcand1 > primcand2
*/
{ double delta1,delta2,ratio1,ratio2 ;
  const primcand_struct *primcand1,*primcand2 ;

  primcand1 = (const primcand_struct *) p_primcand1 ;
  primcand2 = (const primcand_struct *) p_primcand2 ;

/*
  The primary criterion is nondecreasing delta. See promoteSanePivot for
  an optimization.
*/
  delta1 = primcand1->deltakj ;
  delta2 = primcand2->deltakj ;
  if (delta1 < delta2)
  { return (-1) ; }
  else
  if (delta1 > delta2)
  { return (1) ; }
/*
  Equal deltas. Order by nonincreasing pivot ratio.
*/
  ratio1 = primcand1->ratiokj ;
  ratio2 = primcand2->ratiokj ;
  if (ratio1 > ratio2)
    return (-1) ;
  else
  if (ratio1 < ratio2)
    return (1) ;
  else
    return (0) ; }


static void promoteSanePivot (primcand_struct *outcands)
/*
  This routine attempts to promote a candidate with a sane pivot to the front
  of the candidate list.

  The observation is this: numerically unstable (mad) pivots are often
  small.  If the delta of a sane pivot is small enough, we won't cause primal
  infeasibility.  Suppose we have two candidates x<s> (sane) and x<m> (mad)
  to leave as x<j> enters. If (delta<s>-delta<m>)*a<m,j> < tol, then it's a
  good bet we can promote x<s> past x<m> in the pivot list.

  We also want to prefer hard pivots over soft pivots, which are often
  degenerate (i.e., x<k> is already at bound and could be pivoted out with no
  motion). But in a pinch, we'll take a degenerate soft pivot over a mad
  pivot.

  So --- our goal is to walk along outcands until we find the first sane hard
  pivot and then try to promote it to the front of the list. For hard limits,
  we have to be under the tolerance. Soft limits can be ignored.

  Parameters:
    outcands:	(i) the candidate list;
		    outcands[0].ndx is the number of candidates
		(o) the candidate list, possibly rearranged as described

  Returns: undefined
*/

{ int ndx,candcnt,firsthardsane,firstsoftsane ;
  double tol ;
  primcand_struct sane,insane ;

/*
  If outcands[1] is sane and hard, all is copacetic and we can return with no
  further effort.
*/
  if (outcands[1].madpiv == FALSE && outcands[1].hard == TRUE) return ;
/*
  It's a mad, mad, mad, mad world. Scan for sane pivots, ending the scan when
  we find the first hard, sane pivot. If there are no sane pivots, we're done.
*/
  candcnt = outcands[0].ndx ;
  firstsoftsane = -1 ;
  firsthardsane = -1 ;
  for (ndx = 1 ; ndx <= candcnt && firsthardsane == -1 ; ndx++)
  { if (outcands[ndx].madpiv == FALSE)
    { if (outcands[ndx].hard == TRUE)
      { firsthardsane = ndx ; }
      else
      if (firstsoftsane == -1)
      { firstsoftsane = ndx ; } } }
  if (firsthardsane < 0 && firstsoftsane < 0) return ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t      ") ;
    if (firstsoftsane > 0)
    { ndx = outcands[firstsoftsane].ndx ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"first soft sane %s (%d) at %d",
		  consys_nme(dy_sys,'v',ndx,FALSE,NULL),ndx,firstsoftsane) ; }
    if (firsthardsane > 0)
    { ndx = outcands[firsthardsane].ndx ;
      if (firstsoftsane > 0) dyio_outfmt(dy_logchn,dy_gtxecho,", ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"first hard sane %s (%d) at %d",
		  consys_nme(dy_sys,'v',ndx,FALSE,NULL),ndx,firsthardsane) ; }
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
# endif
/*
  First try to promote a sane pivot with a hard limit.  We have a sane pivot at
  firsthardsane and mad pivots up to firsthardsane-1. Work back towards the
  front of the list.
*/
  tol = log10(dy_tols->pfeas/dy_tols->zero)/2 ;
  tol = dy_tols->zero*pow(10.0,tol) ;
  if (firsthardsane > 1)
  { for (ndx = firsthardsane ; ndx > 1 ; ndx--)
    { sane = outcands[ndx] ;
      insane = outcands[ndx-1] ;
      if (insane.hard == FALSE ||
	  ((sane.deltakj-insane.deltakj)*insane.abarkj < tol))
      { outcands[ndx-1] = sane ;
	outcands[ndx] = insane ; }
      else
      { break ; } }
#   ifndef DYLP_NDEBUG
    if ((dy_opts->print.pivoting >= 2 && ndx == 1) ||
	dy_opts->print.pivoting >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t      promoted hard sane %s (%d) to %d. ",
		  consys_nme(dy_sys,'v',outcands[1].ndx,FALSE,NULL),
		  outcands[1].ndx,ndx) ; }
#   endif
#   ifdef DYLP_STATISTICS
    if (dy_stats != NULL)
    { if (ndx == 1) dy_stats->pmulti.promote++ ; }
#   endif
    if (ndx == 1) return ; }
/*
  Try again with the soft sane pivot.
*/
  if (firstsoftsane > 1)
  { for (ndx = firstsoftsane ; ndx > 1 ; ndx--)
    { sane = outcands[ndx] ;
      insane = outcands[ndx-1] ;
      if (insane.hard == FALSE ||
	  ((sane.deltakj-insane.deltakj)*insane.abarkj < tol))
      { outcands[ndx-1] = sane ;
	outcands[ndx] = insane ; }
      else
      { break ; } }
#   ifndef DYLP_NDEBUG
    if ((dy_opts->print.pivoting >= 2 && ndx == 1) ||
	dy_opts->print.pivoting >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t      promoted hard sane %s (%d) to %d. ",
		  consys_nme(dy_sys,'v',outcands[1].ndx,FALSE,NULL),
		  outcands[1].ndx,ndx) ; }
#   endif
#   ifdef DYLP_STATISTICS
    if (dy_stats != NULL)
    { if (ndx == 1) dy_stats->pmulti.promote++ ; }
#   endif
  }

  return ; }



static dyret_enum scanForPrimOutCands (primcand_struct *outcands,
				       int j, int indir,
				       double *abarj, double maxabarj)

/*
  This routine scans the basic variables for candidates to become the leaving
  variable; the criteria are outlined below in the scan loop.  Candidates are
  stored in the outcands array. For each candidate, we record the index,
  information about the pivot coefficient (value, stability), and the pivot
  operation (direction of motion, bound, delta, hard/soft limit).

  During the scan, we use the usual primal delta limit to track a `best'
  candidate, with pivot ratio used as the tiebreaker. If this candidate looks
  good (a sane pivot with a hard |delta| > 0) at the end of the search, it's
  returned without further ado. Note that while it's possible to generate two
  candidate entries (a soft and a hard) for a given variable, the hard
  candidate will always be the second one, which is the one we'll test below.

  Otherwise, the entries are sorted. The primary criteria is nonincreasing
  primal delta, with stability ratio used as the tie-breaker. The sorted list
  is passed to promoteSanePivot, which will try its best to make sure that the
  top candidate has a sane pivot coefficient and a hard delta limit.

  Parameters:
    outcands:	(i) pointer to array to hold the candidates for entry
		(o) filled and sorted array; in particular
		    outcands[0].ndx is the number of candidates
		    outcands[1..] are the candidates
    j:		index of the entering variable x<j>
    indir:	direction of motion of the entering variable x<j>
    abarj:	The ftran'd pivot column abar<j>, column j of inv(B)N
    maxabarj:	The maximum value in the pivot column

  Returns: dyrOK if there are candidates and no errors
	   dyrDEGEN if the pivot would be degenerate
	   dyrUNBOUND if there are no candidates
	   dyrFATAL on error (should only occur when we're paranoid)
*/

{ int m,k,kpos,reject,candcnt,lastcandcnt ;
  double abarkj,ratiokj,xk,ubk,lbk ;
  double *vub,*vlb ;
  flags statk ;
  bool hard,sort ;
  primcand_struct *outcand,*curbest,best ;
  dyret_enum retval ;

/*
  Unclear that it's useful to allow soft degenerate pivots, and it's definitely
  prone to cycling. Could be made into an option, but just disallow for now.
*/
  const bool allowsoftdegen = FALSE ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "scanForPrimOutCands" ;
# endif

# if !defined(DYLP_NDEBUG) || defined(DYLP_PARANOIA)
  int print ;

  print = dy_opts->print.pivoting ;	/* suppress print in dy_chkpiv */
  dy_opts->print.pivoting = 0 ;
# endif

# ifdef DYLP_PARANOIA
  if (outcands == NULL)
  { errmsg(2,rtnnme,"outcands array") ;
    return (dyrFATAL) ; }
# endif

  retval = dyrFATAL ;

  memset(&best,0,sizeof(primcand_struct)) ;

  m = dy_sys->concnt ;
  vub = dy_sys->vub ;
  vlb = dy_sys->vlb ;
  outcands[0].ndx = -1 ;

# ifndef DYLP_NDEBUG
  if (print >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    gathering candidates to leave ... ") ;
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tVariable\t  x<k>\t\tabar<k,j>\t  delta\t\tDisp") ; } }
# endif
/*
  Check if the entering variable is eligible to be the leaving variable. The
  normal case here is a bound-to-bound flip, but the code is written to
  handle moving a superbasic variable to bound. ratiokj is set to inf on the
  theory that not pivoting at all is the ultimate in stable pivots (not to
  mention the work we avoid by not needing to update the basis inverse).
*/
  candcnt = 0 ;
  curbest = NULL ;
  if ((vlb[j] > -dy_tols->inf && indir == -1) ||
      (vub[j] < dy_tols->inf && indir == 1))
  { candcnt++ ;
    outcand = &outcands[candcnt] ;
    xk = dy_x[j] ;
    outcand->ndx = j ;
    outcand->abarkj = 1.0 ;
    if (indir == -1)
    { outcand->deltakj = xk-dy_sys->vlb[j] ; }
    else
    { outcand->deltakj = dy_sys->vub[j]-xk ; }
    setcleanzero(outcand->deltakj,dy_tols->zero) ;
    outcand->ratiokj = dy_tols->inf ;
    outcand->madpiv = FALSE ;
    outcand->dir = indir ;
    outcand->bnd = indir ;
    outcand->hard = TRUE ;
    curbest = outcand ;
#   ifndef DYLP_NDEBUG
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%-8s (%d)",
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t%8g",xk,-1.0) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\t%8g\t accepted",outcand->deltakj*indir) ; }
#   endif
  }
/*
  Open a loop and scan for suitable basic variables. First we look for
  obvious reasons to reject x<k>, then determine if x<k> imposes a limit on
  delta<kj>.
*/
  lastcandcnt = candcnt ;
  for (kpos = 1 ; kpos <= m ; kpos++)
  { k = dy_basis[kpos] ;
#   ifdef DYLP_PARANOIA
    if (dy_chkstatus(k) == FALSE)
    { outcands[0].ndx = -1 ;
      dy_opts->print.pivoting = print ;
      return (dyrFATAL) ; }
#   endif
/*
  First check if we can reject x<k> out of hand:
    * We're working a restricted subproblem, and x<k> is not included
    * abar<kj> == 0
    * stat<k> == BFR.
  Then check to see if x<k> actually imposes a limit by looking at x<k>,
  lb<k>, and ub<k>.  The order of tests depends on the direction of motion
  for x<k>. For example, if x<k> is decreasing:

    * (x<k> >= ub<k>) implies ub<k> is finite and means we can decrease x<k>
      to ub<k> and pivot. If lb<k> is finite, then we have a choice to go for
      lb<k> and the delta required to reach ub<k> is a soft limit. If lb<k>
      is -inf, then ub<k> is our only choice and this is a hard limit.

    * (ub<k> > x<k> >= lb<k> and lb<k> > -inf) means we can decrease x<k> to
      lb<k> and pivot. This is a hard limit.

    * (ub<k> > x<k> > lb<k> and lb<k> <= -inf) means x<k> imposes no limit on
      delta<j> and is not a candidate to leave.

    * (!belowbnd(xk,ubk)) We're right at bound, or borderline infeasible, still
      within the feasibility tolerance. In any case, neither phase I or II will
      tolerate loss of feasibility.

    * (otherwise) We're well into infeasibility. If this is phase I, trust that
      cbar<j> is correct and overall this is a good bet. If this is phase II,
      well, this shouldn't happen. An accuracy or sanity check will catch it.
*/
    reject = 0 ;
    statk = dy_status[k] ;
    abarkj = abarj[kpos] ;
    xk = quiet_nan(0) ;
    outcand = NULL ;

    if (dy_lp->degen > 0 && dy_degenset[kpos] != dy_lp->degen)
    { reject = -5 ; }
    else
    if (flgon(statk,vstatBFR))
    { reject = -1 ; }
    else
    if (withintol(abarkj,0.0,dy_tols->zero))
    { reject = -2 ; }
    else
    { ubk = vub[k] ;
      lbk = vlb[k] ;
      xk = dy_xbasic[kpos] ;
      if (ubk >= dy_tols->inf || lbk <= -dy_tols->inf)
      { hard = TRUE ; }
      else
      { hard = FALSE ; }
      if (abarkj*indir > 0)  /* x<k> decreasing */
      { if (xk > ubk || (xk == ubk && hard == FALSE && allowsoftdegen == TRUE))
	{ candcnt++ ;
	  outcand = &outcands[candcnt] ;
	  outcand->ndx = k ;
	  outcand->abarkj = fabs(abarkj) ;
	  outcand->deltakj = fabs((ubk-xk)/outcand->abarkj) ;
	  setcleanzero(outcand->deltakj,dy_tols->zero) ;
	  ratiokj = dy_chkpiv(abarkj,maxabarj) ;
	  outcand->ratiokj = ratiokj ;
	  if (ratiokj < 1.0)
	  { outcand->madpiv = TRUE ; }
	  else
	  { outcand->madpiv = FALSE ; }
	  outcand->dir = -1 ;
	  outcand->bnd = 1 ;
	  outcand->hard = hard ; }
	if (lbk > -dy_tols->inf)
	{ if (xk >= lbk)
	  { candcnt++ ;
	    outcand = &outcands[candcnt] ;
	    outcand->ndx = k ;
	    outcand->abarkj = fabs(abarkj) ;
	    outcand->deltakj = fabs((lbk-xk)/outcand->abarkj) ;
	    setcleanzero(outcand->deltakj,dy_tols->zero) ;
	    ratiokj = dy_chkpiv(abarkj,maxabarj) ;
	    outcand->ratiokj = ratiokj ;
	    if (ratiokj < 1.0)
	    { outcand->madpiv = TRUE ; }
	    else
	    { outcand->madpiv = FALSE ; }
	    outcand->dir = -1 ;
	    outcand->bnd = -1 ;
	    outcand->hard = TRUE ; }
	  else
	  if (!belowbnd(xk,lbk))
	  { candcnt++ ;
	    outcand = &outcands[candcnt] ;
	    outcand->ndx = k ;
	    outcand->deltakj = 0 ;
	    outcand->abarkj = fabs(abarkj) ;
	    ratiokj = dy_chkpiv(abarkj,maxabarj) ;
	    outcand->ratiokj = ratiokj ;
	    if (ratiokj < 1.0)
	    { outcand->madpiv = TRUE ; }
	    else
	    { outcand->madpiv = FALSE ; }
	    outcand->dir = -1 ;
	    outcand->bnd = -1 ;
	    outcand->hard = TRUE ; }
	  else
	  { reject = -3 ; } }
	else
	if (xk < ubk)
	{ reject = -4 ; } }
      else  /* x<k> increasing */
      { if (xk < lbk || (xk == lbk && hard == FALSE && allowsoftdegen == TRUE))
	{ candcnt++ ;
	  outcand = &outcands[candcnt] ;
	  outcand->ndx = k ;
	  outcand->abarkj = fabs(abarkj) ;
	  outcand->deltakj = fabs((lbk-xk)/outcand->abarkj) ;
	  setcleanzero(outcand->deltakj,dy_tols->zero) ;
	  ratiokj = dy_chkpiv(abarkj,maxabarj) ;
	  outcand->ratiokj = ratiokj ;
	  if (ratiokj < 1.0)
	  { outcand->madpiv = TRUE ; }
	  else
	  { outcand->madpiv = FALSE ; }
	  outcand->dir = 1 ;
	  outcand->bnd = -1 ;
	  outcand->hard = hard ; }
	if (ubk < dy_tols->inf)
	{ if (xk <= ubk)
	  { candcnt++ ;
	    outcand = &outcands[candcnt] ;
	    outcand->ndx = k ;
	    outcand->abarkj = fabs(abarkj) ;
	    outcand->deltakj = fabs((ubk-xk)/outcand->abarkj) ;
	    setcleanzero(outcand->deltakj,dy_tols->zero) ;
	    ratiokj = dy_chkpiv(abarkj,maxabarj) ;
	    outcand->ratiokj = ratiokj ;
	    if (ratiokj < 1.0)
	    { outcand->madpiv = TRUE ; }
	    else
	    { outcand->madpiv = FALSE ; }
	    outcand->dir = 1 ;
	    outcand->bnd = 1 ;
	    outcand->hard = TRUE ; }
	  else
	  if (!abovebnd(xk,ubk))
	  { candcnt++ ;
	    outcand = &outcands[candcnt] ;
	    outcand->ndx = k ;
	    outcand->deltakj = 0 ;
	    outcand->abarkj = fabs(abarkj) ;
	    ratiokj = dy_chkpiv(abarkj,maxabarj) ;
	    outcand->ratiokj = ratiokj ;
	    if (ratiokj < 1.0)
	    { outcand->madpiv = TRUE ; }
	    else
	    { outcand->madpiv = FALSE ; }
	    outcand->dir = 1 ;
	    outcand->bnd = 1 ;
	    outcand->hard = TRUE ; }
	  else
	  { reject = -3 ; } }
	else
	if (xk > lbk)
	{ reject = -4 ; } } }
/*
  Fast test for an incumbent candidate.
*/
    if (candcnt != lastcandcnt)
    { if (curbest == NULL)
      { curbest = outcand ; }
      else
      if (curbest->deltakj > outcand->deltakj)
      { curbest = outcand ; }
      else
      if (curbest->deltakj == outcand->deltakj &&
	  curbest->ratiokj < outcand->ratiokj)
      { curbest = outcand ; } }

#   ifndef DYLP_NDEBUG
    if (print >= 4 && candcnt != lastcandcnt)
    { for (lastcandcnt++ ; lastcandcnt <= candcnt ; lastcandcnt++)
      { outcand = &outcands[lastcandcnt] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%-8s (%d)",
		    consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t%8g",xk,abarkj) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\t%8g\t accepted",outcand->deltakj*indir) ;
	if (outcand->madpiv == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (mad)") ; }
	if (outcand->deltakj == 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (degen)") ; }
	if (outcand->hard == FALSE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (soft)") ; } }
      lastcandcnt-- ; }
    else
    if (print >= 5)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%-8s (%d)",
		  consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t%8g",xk,abarkj) ;
      switch (reject)
      { case -1:
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\t\trejected -- status %s",
		      dy_prtvstat(statk)) ;
	  break ; }
	case -2:
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\t\trejected -- zero pivot") ;
	  break ; }
	case -3:
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\t\trejected -- borderline infeasible") ;
	  break ; }
	case -4:
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\t\trejected -- no limiting bound") ;
	  break ; }
	case -5:
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\t\trejected -- not in restricted subproblem") ;
	  break ; } } }
#   endif

    lastcandcnt = candcnt ; } /* End of scan loop */
  outcands[0].ndx = candcnt ;
/*
  If we saved a suitable candidate while scanning, we can skip sorting and
  promotion.  Suitable is defined as a hard delta > 0 and a sane pivot.
*/
  sort = TRUE ;
  if (curbest != NULL)
  { if (curbest->hard == TRUE && curbest->deltakj > 0 &&
	curbest->madpiv == FALSE)
    { best = *curbest ;
      sort = FALSE ; } }
  /* sort = TRUE ; */
/*
  If we have no candidates, we're unbounded. If the candidate we've saved while
  scanning looks good (sort == FALSE), then we simply copy it into
  outcands[1].

  If we didn't save a suitable candidate (sort == TRUE), and we do have a list
  of candidates, sort the list and then try to promote the best pivot to the
  front. If we have only one candidate, well, we know it's bad, but there's
  nothing we can do.

  The assignments to outcands[0].[deltakj,madpiv] are targetted at a quirk of
  qsort: it seems to occasionally sort outcands[0]. An optimization of some
  sort, presumably.

  Degeneracy is declared if the top two candidates both have deltakj == 0 and
  the second candidate is a hard bound. This may miss the odd instance, but
  should catch the troublesome cases.
*/
# ifndef DYLP_NDEBUG
  dy_opts->print.pivoting = print ;
# endif
  if (candcnt > 0)
  { if (sort == TRUE)
    { if (candcnt > 1)
      { 
#	ifndef DYLP_NDEBUG
	if (print >= 2)
	{ if (sort == FALSE) dyio_outfmt(dy_logchn,dy_gtxecho,"!") ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"sorting ... ") ; }
#	endif
	outcands[0].deltakj = -dy_tols->inf ;
	outcands[0].madpiv = FALSE ;
	qsort(&outcands[1],candcnt,sizeof(primcand_struct),primcand_cmp) ;
	promoteSanePivot(outcands) ;
#	ifdef DYLP_STATISTICS
	if (dy_stats != NULL) dy_stats->pmulti.nontrivial++ ;
#	endif
	if (outcands[1].deltakj == 0 &&
	    outcands[2].deltakj == 0 && outcands[2].hard == TRUE)
	{ retval = dyrDEGEN ; }
	else
	{ retval = dyrOK ; } }
      else
      { retval = dyrOK ; } }
    if (sort == FALSE)
    { /*
      if (best.ndx != outcands[1].ndx)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\nMISMATCH (%d)\n",dy_lp->tot.iters) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "    sort %s (%d) delta = %g ratio = %g\n",
		    consys_nme(dy_sys,'v',outcands[1].ndx,0,NULL),
		    outcands[1].ndx,outcands[1].deltakj,outcands[1].ratiokj) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "    best %s (%d) delta = %g ratio = %g\n",
		    consys_nme(dy_sys,'v',best.ndx,0,NULL),best.ndx,
		    best.deltakj,best.ratiokj) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "    diff delta = %g ratio = %g\n",
		    outcands[1].deltakj-best.deltakj,
		    outcands[1].ratiokj-best.ratiokj) ; }
       */
      outcands[1] = best ;
      retval = dyrOK ; } }
  else
  { retval = dyrUNBOUND ; }

# ifndef DYLP_NDEBUG
  if (print >= 1)
  { if (print >= 3 && candcnt > 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n   ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tVariable\tratio<ik>\tdelta<k>") ;
      for (m = 1 ; m <= candcnt ; m++)
      { k = outcands[m].ndx ;
	ratiokj = outcands[m].ratiokj ;
	xk = outcands[m].deltakj ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%-8s (%d)",
		    consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t%8g",ratiokj,xk) ;
	if (outcands[m].madpiv == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (mad)") ; }
	if (xk == 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (degen)") ; }
	if (outcands[m].hard == FALSE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (soft)") ; } }
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,"%d candidates.",candcnt) ;
    if (print >= 2 && (retval == dyrOK || retval == dyrDEGEN))
    { k = outcands[1].ndx ;
      if (j != k)
      { kpos = dy_var2basis[k] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    selected %s (%d) = %g to leave pos'n %d at",
		    consys_nme(dy_sys,'v',k,FALSE,NULL),k,
		    dy_xbasic[kpos],kpos) ;
	if (outcands[1].dir > 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," %s = %g, ",
		      (dy_status[k] != vstatBLLB)?"ub":"lb",
		      (dy_status[k] != vstatBLLB)?vub[k]:vlb[k]) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho," %s = %g, ",
		      (dy_status[k] != vstatBUUB)?"lb":"ub",
		      (dy_status[k] != vstatBUUB)?vlb[k]:vub[k]) ; }
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "abar<%d,%d> = %g, ",j,k,abarj[kpos]) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    selected %s (%d) = %g to change to %s = %g, ",
		    consys_nme(dy_sys,'v',k,FALSE,NULL),k,dy_x[k],
		    (outcands[1].dir > 0)?"ub":"lb",
		    (outcands[1].dir > 0)?vub[k]:vlb[k]) ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,"delta = %g.",outcands[1].deltakj) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  " Returning %s.",dy_prtdyret(retval)) ; } }
  dy_opts->print.pivoting = print ;
# endif

  return (retval) ; }



dyret_enum primmultiout (int j, int indir,
			 double *abarj, double maxabarj,
			 int *p_xindx, int *p_outdir, double *p_deltaj)
/*
  Select a candidate to leave the basis. This routine coordinates the process,
  and really has very little work to do. scanForPrimOutCands does the heavy
  lifting, collecting a list of candidates and sorting them so that the
  candidate returned in position #1 of the candidate array is the correct
  choice.

  This routine doesn't report loss of primal feasibility. scanForPrimOutCands
  understands how to do pivot selection in the presence of primal
  infeasibility. We'll hope that any problems are transitory and leave it to
  dylp's regular feasibility checks to catch anything serious.

  Parameters:
    xjndx:	Index of the entering variable x<j>.
    indir:	Direction of motion of x<j>.
		    -1: decreasing from upper bound
		     1: increasing from lower bound
    abarj:	Ftran'd column inv(B)a<j> associated with x<j>.
    maxabarj:	MAX{j}(abar<ij>)
    p_xindx:	(o) Index of the leaving variable x<i>. Also valid for
		    return code dyrLOSTPFEAS (in which case it is the index
		    of the variable where feasibility loss was discovered)
		    and dyrREQCHK (in which case it is the index of the
		    variable whose pivot a<ij> was declared bogus).
    p_outdir:	(o) Direction of motion of x<i>, coded as:
		    -1: decreasing to lower bound
		     1: increasing to upper bound
    p_deltaj:	(o) Absolute value of the allowable change in x<j>.

  
  Returns: dyret_enum code, as follows:
    dyrOK:	a strictly basic leaving variable was successfully selected
		(this includes dirty degeneracy)
    dyrDEGEN:	a basic at bound leaving variable is selected; the pivot will
		be (cleanly) degenerate
    dyrMADPIV:	the pivot coefficient abar<ij> would be numerically unstable
    dyrREQCHK:	a possibly bogus abar<ij> was selected as the pivot, and
		refactoring seems wise before trying to use it
		(basis.etas > 1 is the criteria)
    dyrUNBOUND:	the problem is unbounded
    dyrFATAL:	fatal confusion
*/

{ int m,candcnt ;
  dyret_enum retval ;

  primcand_struct *outcands,*candk ;

/*
  Setup. Potentially, each basic variable can produce a soft and a hard pivot
  delta, and we need to allow for a bound-to-bound pivot. Hence we need a
  candidate array of length 2*m+1.
*/
  retval = dyrINV ;
  *p_xindx = 0 ;
  *p_outdir = 0 ;
  *p_deltaj = -1 ;
  m = dy_sys->concnt ;
  outcands = (primcand_struct *) MALLOC((2*m+1+1)*sizeof(primcand_struct)) ;
# ifdef DYLP_ZEROFAULT
  /* Alignment padding. */
  memset(outcands,0,(2*m+1+1)*sizeof(primcand_struct)) ;
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->pmulti.cnt++ ;
# endif
/*
  Generate a sorted list of candidates to leave. No candidates means we're
  primal unbounded (dyrUNBOUND). Fatal error (dyrFATAL) is possible only if
  we're paranoid. Otherwise, the candidate in position #1 is the choice.

  Degeneracy is indicated when we have multiple candidates with a hard delta
  of 0.
*/
  retval = scanForPrimOutCands(outcands,j,indir,abarj,maxabarj) ;
  if (retval == dyrOK || retval == dyrDEGEN)
  { candcnt = outcands[0].ndx ;
    candk = &outcands[1] ;
    *p_xindx = candk->ndx ;
    *p_outdir = candk->dir ;
    *p_deltaj = candk->deltakj ;
    if (candk->madpiv == TRUE)
    { retval = dyrMADPIV ; }
#   ifdef DYLP_STATISTICS
    if (dy_stats != NULL) dy_stats->pmulti.cands += outcands[0].ndx ;
#   endif
  }
  else
  { *p_xindx = -1 ; }

  FREE(outcands) ;
  return (retval) ; }


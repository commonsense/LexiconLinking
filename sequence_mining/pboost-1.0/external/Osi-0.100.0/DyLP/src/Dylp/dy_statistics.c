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
  This file contains routines related to dylp's statistics collection
  facilities. To activate statistics collection, it's first necessary to
  define the compile-time symbol DYLP_STATISTICS. The client code should then
  perform this sequence of actions:
    * Before calling dylp, call dy_initstats to set up the statistics data
      structure.
    * The pointer returned from dy_initstats should then be passed as one of
      the parameters to the call to dylp.
    * Once dylp has returned, dy_dumpstats can be used to print the contents
      of the statistics structure to a file.
    * When the structure is no longer required, call dy_freestats to free the
      data structure.
*/

#define DYLP_INTERNAL

#include "dylp.h"
#include <float.h>
#include <limits.h>

static char sccsid[] UNUSED = "@(#)dy_statistics.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_statistics.c 269 2009-04-02 05:38:19Z lou $" ;



/* Handy macros for initialising and freeing an lpstats_struct */

#define FREE_AND_CLEAR(zz_var_zz) \
  if (zz_var_zz != NULL) { FREE(zz_var_zz) ; zz_var_zz = NULL ; }

#define ALLOC_OR_CLEAR(zz_field_zz,zz_sze_zz,zz_type_zz) \
  if (zz_field_zz == NULL) \
  { zz_field_zz = (zz_type_zz *) CALLOC((zz_sze_zz+1),sizeof(zz_type_zz)) ; } \
  else \
  { memset(zz_field_zz,0,(zz_sze_zz+1)*sizeof(zz_type_zz)) ; }

void dy_initstats (lpstats_struct **p_lpstats, consys_struct *orig_sys)

/*
  This routine allocates (as necessary) and initialises the data structure
  that dylp uses to collect detailed statistics.

  Parameters:
    p_lpstats:	(i) pointer to an lpstats structure; if null, one will be
		    allocated
		(o) an initialised lpstats structure
    orig_sys:	the constraint system that will be passed to dylp

  Returns: undefined
*/

{ int k,m,n ;
  lpstats_struct *lpstats ;

# ifdef DYLP_PARANOIA

  const char *rtnnme = "dy_initstats" ;

  if (p_lpstats == NULL)
  { errmsg(2,rtnnme,"&lpstats") ;
    return ; }
# endif

/*
  Allocate the basic structure, if needed.
*/
  if (*p_lpstats != NULL)
  { lpstats = *p_lpstats ; }
  else
  { lpstats = (lpstats_struct *) CALLOC(1,sizeof(lpstats_struct)) ; }
/*
  Check the capacity of the arrays allocated for constraint and variable
  statistics. If allocated arrays are too small, free them.
*/
  m = orig_sys->concnt ;
  n = orig_sys->varcnt ;
  if (lpstats->cons.sze <= m)
  { FREE_AND_CLEAR(lpstats->cons.angle)
    FREE_AND_CLEAR(lpstats->cons.actcnt)
    FREE_AND_CLEAR(lpstats->cons.deactcnt)
    FREE_AND_CLEAR(lpstats->cons.init)
    FREE_AND_CLEAR(lpstats->cons.fin)
    lpstats->cons.sze = 0 ; }
  if (lpstats->vars.sze <= n)
  { FREE_AND_CLEAR(lpstats->vars.actcnt)
    FREE_AND_CLEAR(lpstats->vars.deactcnt)
    lpstats->vars.sze = 0 ; }
/*
  Allocate or clear the constraint and variable arrays.
*/
  lpstats->cons.sze = m ;
  lpstats->vars.sze = n ;
  ALLOC_OR_CLEAR(lpstats->cons.angle,m,double)
  ALLOC_OR_CLEAR(lpstats->cons.actcnt,m,int)
  ALLOC_OR_CLEAR(lpstats->cons.deactcnt,m,int)
  ALLOC_OR_CLEAR(lpstats->cons.init,m,bool)
  ALLOC_OR_CLEAR(lpstats->cons.fin,m,bool)
  ALLOC_OR_CLEAR(lpstats->vars.actcnt,n,int)
  ALLOC_OR_CLEAR(lpstats->vars.deactcnt,n,int)
/*
  Clear state counts.
*/
  memset(lpstats->phasecnts,0,(dyDONE+1)*sizeof(int)) ;
/*
  Clear constraint angle summary statistics
*/
  lpstats->angle.max = -FLT_MAX ;
  lpstats->angle.min = FLT_MAX ;
  memset(lpstats->angle.hist,0,DYSTATS_HISTBINS*sizeof(int)) ;
/*
  Clear refactoring statistics
*/
  lpstats->factor.cnt = 0 ;
  lpstats->factor.prevpiv = 0 ;
  lpstats->factor.avgpivs = 0 ;
  lpstats->factor.maxpivs = 0 ;
/*
  Clear/set pivot rejection statistics
*/
  lpstats->pivrej.max = 0 ;
  lpstats->pivrej.mad = 0 ;
  lpstats->pivrej.sing = 0 ;
  lpstats->pivrej.pivtol_red = 0 ;
  lpstats->pivrej.min_pivtol = DBL_MAX ;
  lpstats->pivrej.puntcall = 0 ;
  lpstats->pivrej.puntret = 0 ;
/*
  Clear dual multipivot statistics
*/
  lpstats->dmulti.flippable = 0 ;
  lpstats->dmulti.cnt = 0 ;
  lpstats->dmulti.cands = 0 ;
  lpstats->dmulti.promote = 0 ;
  lpstats->dmulti.nontrivial = 0 ;
  lpstats->dmulti.evals = 0 ;
  lpstats->dmulti.flips = 0 ;
  lpstats->dmulti.pivrnks = 0 ;
  lpstats->dmulti.maxrnk = 0 ;
/*
  Clear primal multipivot statistics
*/
  lpstats->pmulti.cnt = 0 ;
  lpstats->pmulti.cands = 0 ;
  lpstats->pmulti.nontrivial = 0 ;
  lpstats->pmulti.promote = 0 ;
/*
  Clear infeasibility statistics
*/
  lpstats->infeas.prevpiv = 0 ;
  lpstats->infeas.maxcnt = 0 ;
  lpstats->infeas.totpivs = 0 ;
  lpstats->infeas.maxpivs = 0 ;
  lpstats->infeas.chgcnt1 = 0 ;
  lpstats->infeas.chgcnt2 = 0 ;
/*
  Clear antidegeneracy statistics
*/
  for (k = 0 ; k < DYSTATS_MAXDEGEN ; k++)
  { lpstats->pdegen[k].cnt = 0 ;
    lpstats->pdegen[k].avgsiz = 0 ;
    lpstats->pdegen[k].maxsiz = 0 ;
    lpstats->pdegen[k].totpivs = 0 ;
    lpstats->pdegen[k].avgpivs = 0 ;
    lpstats->pdegen[k].maxpivs = 0 ; }
  for (k = 0 ; k < DYSTATS_MAXDEGEN ; k++)
  { lpstats->ddegen[k].cnt = 0 ;
    lpstats->ddegen[k].avgsiz = 0 ;
    lpstats->ddegen[k].maxsiz = 0 ;
    lpstats->ddegen[k].totpivs = 0 ;
    lpstats->ddegen[k].avgpivs = 0 ;
    lpstats->ddegen[k].maxpivs = 0 ; }

  *p_lpstats = lpstats ;

  return ; }



void dy_freestats (lpstats_struct **p_lpstats)

/*
  This routine frees an lpstats structure.

  Parameters:
    p_lpstats:	(i) lpstats structure to be freed
		(o) NULL

  Returns: undefined
*/

{ lpstats_struct *lpstats ;

# ifdef DYLP_PARANOIA

  const char *rtnnme = "dy_freestats" ;

  if (p_lpstats == NULL)
  { errmsg(2,rtnnme,"&lpstats") ;
    return ; }
# endif

  lpstats = *p_lpstats ;
  *p_lpstats = NULL ;

# ifdef DYLP_PARANOIA
  if (lpstats == NULL)
  { errmsg(2,rtnnme,"lpstats") ;
    return ; }
# endif

  FREE_AND_CLEAR(lpstats->cons.angle)
  FREE_AND_CLEAR(lpstats->cons.actcnt)
  FREE_AND_CLEAR(lpstats->cons.deactcnt)
  FREE_AND_CLEAR(lpstats->cons.init)
  FREE_AND_CLEAR(lpstats->cons.fin)
  FREE_AND_CLEAR(lpstats->vars.actcnt)
  FREE_AND_CLEAR(lpstats->vars.deactcnt)

  FREE(lpstats) ;

  return ; }

#undef ALLOC_OR_CLEAR
#undef FREE_AND_CLEAR



void dy_finalstats (lpstats_struct *lpstats)

/*
  Various finishing activities:
    * Scan the final active constraint system and set the cons.fin array.
    * Copy final pivot and iteration totals from dy_lp.
  Not much at this point, but it's broken out in anticipation that it may
  one day do more.

  Parameters: none

  Returns: undefined
*/

{ int i,k ;

# ifdef DYLP_PARANOIA

  const char *rtnnme = "dy_finalstats" ;

  if (lpstats == NULL)
  { errmsg(2,rtnnme,"lpstats") ;
    return ; }
# endif

/*
  Scan the active constraints and record the ones that are still active.
*/
  for (k = 1 ; k <= dy_sys->concnt ; k++)
  { i = dy_actcons[k] ;
    if (i > 0)
    { lpstats->cons.fin[i] = TRUE ; } }
/*
  Pick up final pivot and iteration totals.
*/
  lpstats->tot.iters = dy_lp->tot.iters ;
  lpstats->tot.pivs = dy_lp->tot.pivs ;

  return ; }



void dy_dumpstats (ioid chn, bool echo,
		   lpstats_struct *lpstats, consys_struct *orig_sys)

/*
  This is a utility routine which will dump an lpstats_struct in a readable
  manner.

  Parameters:
    chn:	file channnel for output
    echo:	TRUE to echo to stdout, FALSE otherwise
    lpstats:	the statistics structure
    orig_sys:	The constraint system
  
  Returns: undefined
*/

{ int i,j,m,n,ndx,tot ;
  int totvact,maxvact,minvact,totvdeact,maxvdeact,minvdeact ;
  int totcact,maxcact,mincact,totcdeact,maxcdeact,mincdeact ;
  int ineqcnt,nearcnt,perpcnt,farcnt ;
  int initallcact,initcact,fincact,right_init,wrong_init,right_fin,wrong_fin ;
  int rinear,winear,rfnear,wfnear,
      riperp,wiperp,rfperp,wfperp,
      rifar,wifar,rffar,wffar ;
  double anglei,degperbin,brklow,brkhi,temp ;
  contyp_enum ctypi ;
  bool initi,fini ;
  const char *rtnnme = "dy_dumpstats" ;

  double outbrks[] = { 30, 60, 85, 90, 95, 120, 150, 180 } ;
  int outbins = sizeof(outbrks)/sizeof(double) ;
  int outhist[sizeof(outbrks)/sizeof(int)+1] ;
  char outbuf[20] ;
  int histfld = 9 ;

  if (lpstats == NULL)
  { dyio_outfmt(chn,echo,"\n\n<< %s: NULL lpstats structure! >>\n",rtnnme) ;
    return ; }
  
  dyio_outfmt(chn,echo,"\n\nLP statistics:") ;

  dyio_outfmt(chn,echo,"\n  Angle of constraints to objective (degrees):") ;
  dyio_outfmt(chn,echo,"\n\tmax: %g\tmin: %g\n",
	      lpstats->angle.max,lpstats->angle.min) ;

  degperbin = 180/(DYSTATS_HISTBINS-1) ;
  j = 0 ;
  brkhi = 0 ;
  for (i = 0 ; i < outbins ; i++)
  { brklow = brkhi ;
    brkhi = outbrks[i] ;
    if (i < outbins/2)
    { (void) dyio_outfxd(outbuf,histfld,'c',
			 "[%d-%d)",((int) brklow),((int) brkhi)) ; }
    else
    { (void) dyio_outfxd(outbuf,histfld,'c',
			 "(%d-%d]",((int) brklow),((int) brkhi)) ; }
    dyio_outfmt(chn,echo," %s",outbuf) ;
    if (brkhi == 90) dyio_outfmt(chn,echo,"   [90]  ") ;
    outhist[i] = 0 ;
    for ( ; j*degperbin < brkhi ; j++)
    { outhist[i] += lpstats->angle.hist[j] ; } }
  dyio_outchr(chn,echo,'\n') ;
  for (i = 0 ; i < outbins/2 ; i++)
  { (void) dyio_outfxd(outbuf,histfld,'c',"%d",outhist[i]) ;
    dyio_outfmt(chn,echo," %s",outbuf) ; }
  dyio_outfxd(outbuf,8,'c',"%d",lpstats->angle.hist[DYSTATS_HISTBINS-1]) ;
  dyio_outfmt(chn,echo," %s",outbuf) ;
  for ( ; i < outbins ; i++)
  { (void) dyio_outfxd(outbuf,histfld,'c',"%d",outhist[i]) ;
    dyio_outfmt(chn,echo," %s",outbuf) ; }

  dyio_outfmt(chn,echo,
	      "\n  Factoring: %d refactorisations",lpstats->factor.cnt) ;
  dyio_outfmt(chn,echo,"\n\trefactor interval (pivots): avg. %.2f\tmax %d",
	      lpstats->factor.avgpivs,lpstats->factor.maxpivs) ;

  dyio_outfmt(chn,echo,
	      "\n  Pivot rejection: %d unstable, %d singular, max len. %d.",
	      lpstats->pivrej.mad,lpstats->pivrej.sing,lpstats->pivrej.max) ;
  if (lpstats->pivrej.max > 0)
  { dyio_outfmt(chn,echo,"\n\tpivot tolerance: %d reductions, min. tol = %g.",
	        lpstats->pivrej.pivtol_red,lpstats->pivrej.min_pivtol) ;
    dyio_outfmt(chn,echo,"\n\tpunts: %d called, %d returned.",
	        lpstats->pivrej.puntcall,lpstats->pivrej.puntret) ; }
  
  dyio_outfmt(chn,echo,"\n  Reduction of infeasibility:") ;
  dyio_outfmt(chn,echo," maximum number of infeasible vars: %d",
	      lpstats->infeas.maxcnt) ;
  if (lpstats->infeas.maxcnt > 0)
  { tot = lpstats->infeas.chgcnt1+lpstats->infeas.chgcnt2 ;
    temp = ((double) lpstats->infeas.totpivs)/tot ;
    dyio_outfmt(chn,echo,"\n\tpivots to reduce: tot. %d, avg. %.2f, max %d",
	        lpstats->infeas.totpivs,temp,lpstats->infeas.maxpivs) ;
    dyio_outfmt(chn,echo,"\n\tsingle change: %d, multiple change: %d",
	        lpstats->infeas.chgcnt1,lpstats->infeas.chgcnt2) ; }

  dyio_outfmt(chn,echo,"\n  Primal Multipivot: %d calls, %d sort, %d promote",
	      lpstats->pmulti.cnt,lpstats->pmulti.nontrivial,
	      lpstats->pmulti.promote) ;
  if (lpstats->pmulti.cnt > 0)
  { temp = ((float) lpstats->pmulti.cands)/lpstats->pmulti.cnt ;
    dyio_outfmt(chn,echo,"\n\tcandidates: tot. %d\tavg. %.2f",
	        lpstats->pmulti.cands,temp) ; }
  
  dyio_outfmt(chn,echo,
	      "\n  Dual Multipivot: %d calls, %d promote, %d multipivot",
	      lpstats->dmulti.cnt,lpstats->dmulti.promote,
	      lpstats->dmulti.nontrivial) ;
  if (lpstats->dmulti.cnt > 0)
  { if (lpstats->dmulti.nontrivial > 0)
    { temp = ((float) lpstats->dmulti.pivrnks)/lpstats->dmulti.nontrivial ; }
    else
    { temp = 0 ; }
    dyio_outfmt(chn,echo,
		", avg. rank %.2f, max %d",temp,lpstats->dmulti.maxrnk) ;
    temp = ((float) lpstats->dmulti.cands)/lpstats->dmulti.cnt ;
    dyio_outfmt(chn,echo,"\n\tcandidates: tot. %d\tavg. %.2f",
	        lpstats->dmulti.cands,temp) ;
    temp = ((float) lpstats->dmulti.evals)/lpstats->dmulti.cnt ;
    dyio_outfmt(chn,echo,"\n\tcolumn evals: tot. %d\tavg. %.2f",
	        lpstats->dmulti.evals,temp) ;
    temp = ((float) lpstats->dmulti.flips)/lpstats->dmulti.cnt ;
    dyio_outfmt(chn,echo,"\n\tbound flips: tot. %d\tavg. %.2f",
	        lpstats->dmulti.flips,temp) ; }
  temp = ((float) lpstats->dmulti.flippable)/orig_sys->varcnt ;
  dyio_outfmt(chn,echo,"\n\tflippable vars: %d (%.2f%%)",
	      lpstats->dmulti.flippable,temp*100) ;
  
  dyio_outfmt(chn,echo,"\n  Degeneracy:") ;
  if (lpstats->pdegen[0].cnt == 0)
  { dyio_outfmt(chn,echo,"\n    Primal: inactive") ; }
  else
  { dyio_outfmt(chn,echo,
		"\n    Primal:\tLevel\tEntries\t     Vars.\t\tPivots") ;
    dyio_outfmt(chn,echo,"\n\t\t\t\tAvg.\tMax\tTot.\tAvg.\tMax") ;
    for (ndx = 1 ; ndx <= lpstats->pdegen[0].cnt ; ndx++)
    { dyio_outfmt(chn,echo,"\n\t\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d",
		  ndx,lpstats->pdegen[ndx].cnt,
		  lpstats->pdegen[ndx].avgsiz,lpstats->pdegen[ndx].maxsiz,
		  lpstats->pdegen[ndx].totpivs,lpstats->pdegen[ndx].avgpivs,
		  lpstats->pdegen[ndx].maxpivs) ; }
    dyio_outchr(chn,echo,'\n') ; }
  if (lpstats->ddegen[0].cnt == 0)
  { dyio_outfmt(chn,echo,"\n    Dual: inactive") ; }
  else
  { dyio_outfmt(chn,echo,"\n    Dual:\tLevel\tEntries\t     Vars.\t\tPivots") ;
    dyio_outfmt(chn,echo,"\n\t\t\t\tAvg.\tMax\tTot.\tAvg.\tMax") ;
    for (ndx = 1 ; ndx <= lpstats->ddegen[0].cnt ; ndx++)
    { dyio_outfmt(chn,echo,"\n\t\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d",
		  ndx,lpstats->ddegen[ndx].cnt,
		  lpstats->ddegen[ndx].avgsiz,lpstats->ddegen[ndx].maxsiz,
		  lpstats->ddegen[ndx].totpivs,lpstats->ddegen[ndx].avgpivs,
		  lpstats->ddegen[ndx].maxpivs) ; } }
  dyio_outchr(chn,echo,'\n') ;
/*
  Run through the constraint management information and calculate totals,
  mins, and maxes, and information about how well we did guessing the correct
  initial constraint system. Note that it's possible, under fanatic constraint
  management, to deactivate/activate equalities. But we don't want to include
  them in information about the initial constraint system, where they're always
  included.
*/
  m = orig_sys->concnt ;
  n = orig_sys->varcnt ;
  totcact = 0 ;
  totcdeact = 0 ;
  maxcact = 0 ;
  mincact = INT_MAX ;
  maxcdeact = 0 ;
  mincdeact = INT_MAX ;
  ineqcnt = 0 ;
  nearcnt = 0 ;
  perpcnt = 0 ;
  farcnt = 0 ;
  initallcact = 0 ;
  initcact = 0 ;
  fincact = 0 ;
  right_init = 0 ;
  wrong_init = 0 ;
  right_fin = 0 ;
  wrong_fin = 0 ;
  rinear = 0 ;
  winear = 0 ;
  rfnear = 0 ;
  wfnear = 0 ;
  riperp = 0 ;
  wiperp = 0 ;
  rfperp = 0 ;
  wfperp = 0 ;
  rifar = 0 ;
  wifar = 0 ;
  rffar = 0 ;
  wffar = 0 ;
  for (i = 1 ; i <= m ; i++)
  { j = lpstats->cons.actcnt[i] ;
    totcact += j ;
    if (maxcact < j) maxcact = j ;
    if (mincact > j) mincact = j ;
    j = lpstats->cons.deactcnt[i] ;
    totcdeact += j ;
    if (maxcdeact < j) maxcdeact = j ;
    if (mincdeact > j) mincdeact = j ;
    initi = lpstats->cons.init[i] ;
    if (initi == TRUE) initallcact++ ;
    ctypi = orig_sys->ctyp[i] ;
    if (ctypi == contypGE || ctypi == contypLE || ctypi == contypRNG)
    { ineqcnt++ ;
      anglei = lpstats->cons.angle[i] ;
      if (anglei < 90)
      { farcnt++ ; }
      else
      if (anglei > 90)
      { nearcnt++ ; }
      else
      { perpcnt++ ; }
      if (initi == TRUE) initcact++ ;
      fini = lpstats->cons.fin[i] ;
      if (fini == TRUE) fincact++ ;
      if (initi == fini)
      { if (initi == TRUE)
	{ right_init++ ;
	  if (anglei < 90)
	  { rifar++ ; }
	  else
	  if (anglei > 90)
	  { rinear++ ; }
	  else
	  { riperp++ ; } }
	else
	{ right_fin++ ;
	  if (anglei < 90)
	  { rffar++ ; }
	  else
	  if (anglei > 90)
	  { rfnear++ ; }
	  else
	  { rfperp++ ; } } }
      else
      if (initi == TRUE)
      { wrong_init++ ;
	if (anglei < 90)
	{ wifar++ ; }
	else
	if (anglei > 90)
	{ winear++ ; }
	else
	{ wiperp++ ; } }
      else
      { wrong_fin++ ;
	if (anglei < 90)
	{ wffar++ ; }
	else
	if (anglei > 90)
	{ wfnear++ ; }
	else
	{ wfperp++ ; } } } }
  totvact = 0 ;
  totvdeact = 0 ;
  maxvact = 0 ;
  minvact = INT_MAX ;
  maxvdeact = 0 ;
  minvdeact = INT_MAX ;
  for (j = 1 ; j <= n ; j++)
  { i = lpstats->vars.actcnt[j] ;
    totvact += i ;
    if (maxvact < i) maxvact = i ;
    if (minvact > i) minvact = i ;
    i = lpstats->vars.deactcnt[j] ;
    totvdeact += i ;
    if (maxvdeact < i) maxvdeact = i ;
    if (minvdeact > i) minvdeact = i ; }
/*
  Dump per phase information.
*/
  dyio_outfmt(chn,echo,"\n  State\tEntry\tActivity\tAverage\n") ;
  for (i = dyPRIMAL1 ; i <= dyFORCEFULL ; i++)
  { j = lpstats->phasecnts[i] ;
    dyio_outfmt(chn,echo,"\n   %2s%2s\t%5d",
	        dy_prtlpphase(((dyphase_enum) i),TRUE),
	        (lpstats->ini_simplex == ((dyphase_enum) i))?"* ":"  ",j) ;
    temp = 0 ;
    switch (i)
    { case dyPRIMAL1:
      { if (j > 0)
	{ temp = ((double) lpstats->p1.pivs)/j ;
	  dyio_outfmt(chn,echo,"\t%8d\t%7.1f",lpstats->p1.pivs,temp) ; }
	break ; }
      case dyPRIMAL2:
      { if (j > 0)
	{ temp = ((double) lpstats->p2.pivs)/j ;
	  dyio_outfmt(chn,echo,"\t%8d\t%7.1f",lpstats->p2.pivs,temp) ; }
	break ; }
      case dyDUAL:
      { if (j > 0)
	{ temp = ((double) lpstats->d2.pivs)/j ;
	  dyio_outfmt(chn,echo,"\t%8d\t%7.1f",lpstats->d2.pivs,temp) ; }
	break ; }
      case dyPURGEVAR:
      { if (j > 0)
	{ dyio_outfmt(chn,echo,
		      "\t%8d\t%7.1f",totvdeact,((float) totvdeact)/j) ; }
	break ; }
      case dyADDVAR:
      { if (j > 0)
	{ dyio_outfmt(chn,echo,"\t%8d\t%7.1f",totvact,((float) totvact)/j) ; }
	break ; }
      case dyPURGECON:
      { if (j > 0)
	{ dyio_outfmt(chn,echo,
		      "\t%8d\t%7.1f",totcdeact,((float) totcdeact)/j) ; }
	break ; }
      case dyADDCON:
      { if (j > 0)
	{ dyio_outfmt(chn,echo,"\t%8d\t%7.1f",totcact-initallcact,
		 ((float) (totcact-initallcact))/j) ; }
	break ; } } }
  dyio_outchr(chn,echo,'\n') ;
/*
  Dump initial and final constraint summary information.
*/
  dyio_outfmt(chn,echo,"\n  Initial/Final Constraint System: ") ;
  dyio_outfmt(chn,echo,"%d inequalities, initial %d, final %d",
	      ineqcnt,initcact,fincact) ;
  dyio_outfmt(chn,echo,"\n\t%d near, %d perp, %d far",nearcnt,perpcnt,farcnt) ;
  dyio_outfmt(chn,echo,"\n\ttotal: %d/%d active, %d/%d inactive.",
	      right_init,wrong_fin,right_fin,wrong_init) ;
  dyio_outfmt(chn,echo,"\n\tnear:  %d/%d active, %d/%d inactive.",
	      rinear,wfnear,rfnear,winear) ;
  dyio_outfmt(chn,echo,"\n\tperp:  %d/%d active, %d/%d inactive.",
	      riperp,wfperp,rfperp,wiperp) ;
  dyio_outfmt(chn,echo,"\n\tfar:   %d/%d active, %d/%d inactive.\n",
	      rifar,wffar,rffar,wifar) ;
/*
  That takes care of the summaries. Now the detailed statistics on constraints.
*/
  dyio_outfmt(chn,echo,
	      "\n    Constraint\tType   Angle\tInit\tAct\tDeact\tFinal\n") ;
  for (i = 1 ; i <= m ; i++)
  { dyio_outfmt(chn,echo,"\n%5d %s\t%2s %9.5f\t",
	        i,consys_nme(orig_sys,'c',i,FALSE,NULL),
	        consys_prtcontyp(orig_sys->ctyp[i]),
	        lpstats->cons.angle[i]) ;
    if (lpstats->cons.init[i] == TRUE) dyio_outfmt(chn,echo,"Y") ;
    dyio_outfmt(chn,echo,"\t%d\t%d\t",
	   lpstats->cons.actcnt[i],lpstats->cons.deactcnt[i]) ;
    if (lpstats->cons.fin[i] == TRUE) dyio_outfmt(chn,echo,"Y") ;
    }
  dyio_outchr(chn,echo,'\n') ;

  return ; }

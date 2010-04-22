/*
  This file is a part of the OsiDylp LP distribution.

        Copyright (C) 2005 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin St., Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
  This file contains la Duenna, which fusses after each pivot. A spin doctor
  might have been better, but I tend toward anachronism.  There are also
  routines that handle refactoring and primal and dual accuracy checks.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_duenna.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_duenna.c 94 2006-06-29 23:06:51Z lou $" ;



static bool groombasis ()

/*
  This routine runs after direct recalculation of the primal variables to
  make sure no basic variables were misclassified due to accumulated numeric
  inaccuracy in the course of iterative updates.  If it finds a misclassified
  variable, it resets the status to the appropriate value.

  The intent of the routine is to roll with minor corruption due to
  accumulated numeric inaccuracy. Deciding whether corruption is minor,
  however, is not entirely simple.  We no longer have the previous variable
  values, so we're just guessing by comparing the old status with what the
  status ought to be. From experience, having the old values wouldn't improve
  things much. The grooming tolerance (i.e., the amount of error we'll
  tolerate before considering an abort) is large:
    1.0e05*(bogus multiplier)*(pfeas tolerance).

  The rules for considering an abort are pretty simple:
    * If we started at a bound, we have to remain within the grooming
      tolerance for the bound value.
    * If we move over a bound, we have to remain within the grooming
      tolerance for the bound value.
    * Anything else is grounds for an error.
  In the end, though, it seems to require informed human judgement. When
  you're debugging, aborts are often useful. When you're doing production
  runs, well, might as well let the accuracy maintenance algorithms earn
  their keep. dylp provides an option, which you can set to silent, warn,
  or abort.

  Status codes BFX and BFR are assigned initially based on the lower and
  upper bounds. BFR should never change, but BFX can change to BLLB or BUUB
  if somewhere the variable goes infeasible.  Groombasis may need to set BFX
  when a variable moves from out-of-bound (BLLB, BUUB) to at-bound (BLB,
  BUB), and the bounds are equal.

  One more complication: If the primal antidegeneracy mechanism is active,
  the true value (dy_x) of variables in the degenerate set (dy_degenset[bpos]
  > 0) must be at bound. For variables not in the degenerate set, the value
  in dy_x should match the value in dy_xbasic.  If either of these assertions
  fails, inaccuracy has accumulated as we pivoted the permuted values.  The
  only way to restore accuracy is to back out the antidegeneracy mechanism,
  begin the basis scan again, fix whatever's necessary, and request a
  refactor.

  The routine always checks the full basis (even if we're going to fail,
  might as well know the extent of the damage).

  Parameters: none

  Returns: TRUE if the basis is correct or only correctable errors were
	   encountered; FALSE otherwise.
*/

{ int i,ipos,staterrs ;
  double xi,*vlb,vlbi,*vub,vubi,tol ;
  flags stati,newstati,quali ;
  char statbuf[32] ;
  bool retval,statok,backout ;
  
  const char *rtnnme = "groombasis" ;

# ifndef DYLP_NDEBUG
  int print ;
# endif

# ifndef DYLP_NDEBUG
  switch (dy_lp->phase)
  { case dyPRIMAL1:
    { print = dy_opts->print.phase1 ;
      break ; }
    case dyPRIMAL2:
    { print = dy_opts->print.phase2 ;
      break ; }
    case dyDUAL:
    { print = dy_opts->print.phase2 ;
      break ; }
    case dyADDCON:
    { print = dy_opts->print.conmgmt ;
      break ; }
    case dyFORCEDUAL:
    case dyFORCEPRIMAL:
    case dyFORCEFULL:
    { print = maxx(dy_opts->print.conmgmt,dy_opts->print.varmgmt) ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# endif

  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
  tol = dy_tols->pfeas*dy_tols->bogus*1.0e05 ;
  backout = FALSE ;
  retval = TRUE ;
/*
  The outer loop is here to allow us to restart the grooming operation
  after backing out the antidegeneracy mechanism. If things go well, we
  won't need it. The inner loop walks the basis, checking the variables.
*/
  while (TRUE)
  { staterrs = 0 ;
    for (ipos = 1 ; ipos <= dy_sys->concnt ; ipos++)
    { i = dy_basis[ipos] ;
      stati = getflg(dy_status[i],vstatSTATUS) ;
      quali = getflg(dy_status[i],vstatQUALS) ;
      vlbi = vlb[i] ;
      vubi = vub[i] ;
      statok = TRUE ;
/*
  If the antidegeneracy mechanism is active in primal simplex, check that the
  true value is as it should be.  Variables in the degenerate set should be
  at bound (but allow for the fact we might have pivoted a free variable into
  the slot); variables not in the degenerate set should be unchanged and thus
  should agree with the values held in xbasic.  If we've lost accuracy, call
  dy_degenout to back out the antidegeneracy perturbations and start grooming
  the basis from the beginning. We can ignore the return value from
  dy_degenout --- it'll be REQCHK, and once we're done grooming the basis
  we'll be doing just that.
*/
      if (dy_lp->degen > 0 &&
	  (dy_lp->phase == dyPRIMAL1 || dy_lp->phase == dyPRIMAL2))
      { xi = dy_x[i] ;
	if (dy_degenset[ipos] > 0)
	{
#	  ifdef PARANOIA
	  if (dy_lp->degen < dy_degenset[ipos])
	  { errmsg(390,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,ipos,dy_degenset[ipos],
		   consys_nme(dy_sys,'v',i,FALSE,NULL),i,dy_prtvstat(stati),
		   dy_lp->degen) ;
	    retval = FALSE ;
	    continue ; }
#	  endif
	  if (atbnd(xi,vlbi))
	  { /* dy_x[i] = vlbi ; */ }
	  else
	  if (atbnd(xi,vubi))
	  { /* dy_x[i] = vubi ; */ }
	  else
	  if (flgoff(stati,vstatBFR))
	  { backout = TRUE ; } }
	else
	{ if (atbnd(dy_xbasic[ipos],xi))
	  { /* dy_xbasic[ipos] = xi ; */ }
	  else
	  { backout = TRUE ; } }
	if (backout == TRUE)
	{
#	  ifndef DYLP_NDEBUG
	  if (print >= 1)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		 "\n\tbacking out antidegeneracy at (%s)%d due to accumulated",
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho," error;\n\ttrue %s (%d) = %g ;",
		        consys_nme(dy_sys,'v',i,FALSE,NULL),i,xi) ;
	    if (dy_degenset[ipos] > 0)
	    { if (fabs(xi-vlbi) < fabs(xi-vubi))
	      { dyio_outfmt(dy_logchn,dy_gtxecho,
			    " lb = %g ; |x-lb| = %g, tol %g.",
			    vlbi,fabs(xi-vlbi),
			    dy_tols->zero*(1+fabs(vlbi))) ; }
	      else
	      { dyio_outfmt(dy_logchn,dy_gtxecho,
			    " ub = %g ; |x-ub| = %g, tol %g.",
			    vubi,fabs(xi-vubi),
			    dy_tols->zero*(1+fabs(vubi))) ; } }
	    else
	    { dyio_outfmt(dy_logchn,dy_gtxecho,
			  "xbasic = %g ; |x-xbasic| = %g, tol %g.",
			  dy_xbasic[ipos],fabs(xi-dy_xbasic[ipos]),
			  dy_tols->zero*(1+fabs(xi))) ; } }
  #	  endif
	  (void) dy_degenout(0) ;
	  break ; } }
/*
  Do the case checking. Nothing complicated here, just a lot of cases.
*/
      xi = dy_xbasic[ipos] ;
      newstati = stati ;
      switch (stati)
      { case vstatBLLB:
	{ if (!belowbnd(xi,vlbi))
	  { if (!withintol(xi,vlbi,tol*(1.0+fabs(vlbi)))) statok = FALSE ;
	    if (atbnd(xi,vlbi))
	    { if (atbnd(vlbi,vubi))
		newstati = vstatBFX ;
	      else
		newstati = vstatBLB ;
	      /* xi = vlbi ; */ }
	    else
	    if (belowbnd(xi,vubi))
	    { newstati = vstatB ; }
	    else
	    if (atbnd(xi,vubi))
	    { newstati = vstatBUB ;
	      /* xi = vubi ; */ }
	    else
	    { newstati = vstatBUUB ; } }
	  break ; }
	case vstatBLB:
	{ if (!atbnd(xi,vlbi))
	  { if (!withintol(xi,vlbi,tol*(1.0+fabs(vlbi)))) statok = FALSE ;
	    if (xi < vlbi)
	    { newstati = vstatBLLB ; }
	    else
	    if (belowbnd(xi,vubi))
	    { newstati = vstatB ; }
	    else
	    if (atbnd(xi,vubi))
	    { newstati = vstatBUB ;
	      /* xi = vubi ; */ }
	    else
	    { newstati = vstatBUUB ; } }
	  else
	  { /* xi = vlbi ; */ }
	  break ; }
	case vstatB:
	{ if (!(abovebnd(xi,vlbi) && belowbnd(xi,vubi)))
	  { if (atbnd(xi,vlbi))
	    { newstati = vstatBLB ;
	      /* xi = vlbi ; */ }
	    else
	    if (atbnd(xi,vubi))
	    { newstati = vstatBUB ;
	      /* xi = vubi ; */ }
	    else
	    if (xi > vubi)
	    { if (!withintol(xi,vubi,tol*(1.0+fabs(vlbi)))) statok = FALSE ;
	      newstati = vstatBUUB ; }
	    else
	    { if (!withintol(xi,vlbi,tol*(1.0+fabs(vlbi)))) statok = FALSE ;
	      newstati = vstatBLLB ; } }
	  break ; }
	case vstatBUB:
	{ if (!atbnd(xi,vubi))
	  { if (!withintol(xi,vubi,tol*(1.0+fabs(vubi)))) statok = FALSE ;
	    if (xi > vubi)
	    { newstati = vstatBUUB ; }
	    else
	    if (abovebnd(xi,vlbi))
	    { newstati = vstatB ; }
	    else
	    if (atbnd(xi,vlbi))
	    { newstati = vstatBLB ;
	      /* xi = vlbi ; */ }
	    else
	    { newstati = vstatBLLB ; } }
	  else
	  { /* xi = vubi ; */ }
	  break ; }
	case vstatBUUB:
	{ if (!abovebnd(xi,vubi))
	  { if (!withintol(xi,vubi,tol*(1.0+fabs(vubi)))) statok = FALSE ;
	    if (atbnd(xi,vubi))
	    { if (atbnd(vlbi,vubi))
		newstati = vstatBFX ;
	      else
		newstati = vstatBUB ;
	      /* xi = vubi ; */ }
	    else
	    if (abovebnd(xi,vlbi))
	    { newstati = vstatB ; }
	    else
	    if (atbnd(xi,vlbi))
	    { newstati = vstatBLB ;
	      /* xi = vlbi ; */ }
	    else
	    { newstati = vstatBLLB ; } }
	  break ; }
	case vstatBFX:
	{ if (!atbnd(xi,vubi))
	  { if (!withintol(xi,vubi,tol*(1.0+fabs(vubi)))) statok = FALSE ;
	    if (xi > vubi)
	      newstati = vstatBUUB ;
	    else
	      newstati = vstatBLLB ; }
	  else
	  { /* xi = vubi ; */ }
	  break ; }
	case vstatBFR:
	{ break ; }
#       ifdef PARANOIA
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (FALSE) ;  }
#       endif
      }
/*
  We've checked and corrected the status. Reset the status vector and (if
  appropriate) dy_x and dy_xbasic.  Then decide how to react, depending on
  the user's choice and whether this was a minor correction (according to the
  rules at the head of the routine).
*/
      setflg(newstati,quali) ;
      dy_status[i] = newstati ;
#     ifndef DYLP_NDEBUG
      if (statok == TRUE && print >= 3 &&
	  stati != getflg(newstati,vstatSTATUS))
      { setflg(stati,quali) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      status of %s (%d) corrected from %s",
		    consys_nme(dy_sys,'v',i,TRUE,NULL),i,dy_prtvstat(stati)) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," to %s.",dy_prtvstat(newstati)) ; }
#     endif

      if (statok == FALSE)
      { staterrs++ ;
	setflg(stati,quali) ;
	strcpy(statbuf,dy_prtvstat(stati)) ;
	switch (dy_opts->groom)
	{ case 0: { break ; }
	  case 1:
	  { if (fabs(vubi-xi) < fabs(vlbi-xi))
	    { warn(372,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters,consys_nme(dy_sys,'v',i,FALSE,NULL),i,
		     statbuf,dy_prtvstat(newstati),
		     vlbi,xi,vubi,"ub",xi-vubi,tol*(1+fabs(vubi))) ; }
	    else
	    { warn(372,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters,consys_nme(dy_sys,'v',i,FALSE,NULL),i,
		     statbuf,dy_prtvstat(newstati),
		     vlbi,xi,vubi,"lb",xi-vlbi,tol*(1+fabs(vlbi))) ; }
	    break ; }
	  default:
	  { if (fabs(vubi-xi) < fabs(vlbi-xi))
	    { errmsg(372,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters,consys_nme(dy_sys,'v',i,FALSE,NULL),i,
		     statbuf,dy_prtvstat(newstati),
		     vlbi,xi,vubi,"ub",xi-vubi,tol*(1+fabs(vubi))) ; }
	    else
	    { errmsg(372,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters,consys_nme(dy_sys,'v',i,FALSE,NULL),i,
		     statbuf,dy_prtvstat(newstati),
		     vlbi,xi,vubi,"lb",xi-vlbi,tol*(1+fabs(vlbi))) ; }
	    break ; } } } }
/*
  If we made it the whole way through the basis, we're done. Otherwise,
  go up and try again.
*/
    if (ipos > dy_sys->concnt) break ; }
/*
  What to do, what to do? If there are no errors, fine. If there are errors,
  force an abort, if that's what the user wants. Otherwise, boost the
  minimum pivot ratio, then use the return value set earlier (it'll be
  TRUE, unless we're in paranoid mode). The current pivot ratio is
  automatically raised to match the minimum. The reason for only raising the
  minimum here is that it's likely we've already failed an accuracy check back
  in dy_accchk, so that the current pivot ratio is already boosted.
*/
# ifndef DYLP_NDEBUG
  if (print >= 1 && staterrs > 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    (%s)%d: %d major status corrections",
	        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,staterrs) ; }
# endif

  if (staterrs > 0)
  { if (dy_opts->groom >= 2)
    { retval = FALSE ; }
    else
    { (void) dy_setpivparms(+1,0) ; } }

  return (retval) ; }



dyret_enum dy_accchk (flags *checks)

/*
  This routine performs refactoring, accuracy, and feasibility checks, as
  indicated by the flags in the checks parameter. The flags are set on return
  to reflect the results of the checks. Refactoring implies that the primal
  and dual variables are recalculated.
  
  If an accuracy check fails, and the basis has not been refactored, the
  routine will call for a refactor and try again. If the check fails again,
  the routine will boost the current pivot selection ratio one step and try
  again. It'll repeat this loop until the pivot ratio is at its maximum.
  Failure of an accuracy check at maximum pivot selection ratio is fatal. If
  we manage to pass the accuracy check(s), but take more than one refactor to
  do it, boost the current pivot selection ratio on the way out.

  If a feasibility check fails, and the checks don't anticipate this, the
  same refactor-and-retry sequence is followed. Anticipation is taken to be
  the presence of the PFQUIET flag for primal feasibility, DFQUIET for dual
  feasibility.

  On the flip side, if we pass all requested checks on the first try, we'll
  relax the current pivot selection tolerance before returning. Since pretty
  much any call except a scheduled accuracy check will force an initial
  refactor, this means we'll only back off the pivot selection tolerances if
  we make it to an accuracy check interval, and pass without refactoring.

  If we refactor, and subsequently pass any requested accuracy checks,
  groombasis is called to make sure the status vector is correct. It's
  purpose is to deal with any changes that need to be made now that we've
  (presumably) corrected numerical inaccuracy in the basis factorisation. It
  also checks that we aren't drifting while the primal antidegeneracy
  mechanism is active. If it is forced to do major status correction, it
  will boost both the minimum and current pivot selection ratio.

  The primal checks are performed on dy_x, which will be equal to dy_xbasic
  except when the primal antidegeneracy mechanism is active. While it's
  active, dy_xbasic is incrementally updated, but dy_x is not (dy_x is
  holding the values which will be used to back out the perturbation). To
  make a meaningful accuracy check under these circumstances, we make a call
  to dy_calcprimals to refresh dy_x before doing the accuracy checks. If
  we're going to refactor first, this isn't necessary because dy_factor will
  do the recalculation.

  Unfortunately, the dual situation is not so clean. We don't have separate
  vectors for all duals and basic duals. In large part this is a consequence
  of faking dual simplex on the primal basis. The portions of dy_accchk that
  perform dual accuracy and feasibility checks compensate on the fly by
  carefully ignoring perturbed values.

  The primal accuracy check is Bx<B> = b - Nx<N>.

  The primal feasibility check is lb<k> <= x<k> <= ub<k> for all k.

  The dual accuracy check is yB = c<B>.

  The dual feasibility check is cbar<j> = c<j>-dot(y,a<j>) > 0 for variables
  at lower bound and cbar<j> < 0 for variables at upper bound (primal
  optimality). We take the opportunity to reset dy_cbar (which has been
  iteratively updated) with a fresh value.

  Parameter:
    checks:	(i) ladFACTOR to request refactoring prior to accuracy or
			      feasibility checks
		    ladEXPAND will force expansion of the basis prior to
		    refactoring
		    ladPRIMALCHK to request a primal accuracy check
		    ladPRIMFEAS to request a primal feasibility check
		    ladPFQUIET to suppress actions to recover primal
		      feasibility
		    ladDUALCHK to request a dual accuracy check
		    ladDUALFEAS to request a dual feasibility check
		    ladDFQUIET to suppress actions to recover dual
		      feasibility
		(o) failure of a test will set the flag for that test.
		    ladFACTOR will be set if the basis was refactored,
		    ladEXPAND if it was expanded, and ladPRIMALS and ladDUALS
		    will be set if they were calculated.
  Returns: dyrOK if the accuracy check calculations are complete, dyrPATCHED if
	   the only bump was a refactor resulting in a patched basis, otherwise
	   the return code from dy_factor or dyrFATAL for problems originating
	   here.
*/

{ int xkpos,xkndx,pkndx,cndx,pfeascnt,refactorcnt ;
  double *primalerrs,normb,normc,dualresid,primalresid,cbarj,pinfeas ;
  pkvec_struct *pkcol ;
  flags results,vstat,factorflags ;
  bool dorefactor,tryagain,dualDegen ;
  dyret_enum factorresult ;
  const char *rtnnme = "dy_accchk" ;

# ifndef DYLP_NDEBUG
  int print,dfeascnt ;
  double *dualerrs,avgerr,*dfeaserrs,dinfeas ;

  dualerrs = NULL ;
  dfeaserrs = NULL ;
  dfeascnt = -1 ;
  dinfeas = quiet_nan(0) ;
# endif

# ifndef DYLP_NDEBUG
  switch (dy_lp->phase)
  { case dyPRIMAL1:
    { print = dy_opts->print.phase1 ;
      break ; }
    case dyPRIMAL2:
    { print = dy_opts->print.phase2 ;
      break ; }
    case dyDUAL:
    { print = dy_opts->print.dual ;
      break ; }
    case dyADDCON:
    { print = dy_opts->print.conmgmt ;
      break ; }
    case dyADDVAR:
    { print = dy_opts->print.varmgmt ;
      break ; }
    case dyINIT:
    { print = dy_opts->print.crash ;
      break ; }
    case dyFORCEDUAL:
    case dyFORCEPRIMAL:
    case dyFORCEFULL:
    { print = maxx(dy_opts->print.conmgmt,dy_opts->print.varmgmt) ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (dyrFATAL) ; } }
  
  if (print >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    dy_accchk: ") ;
    if (flgon(*checks,ladFACTOR))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"factor") ;
      if (flgon(*checks,ladEXPAND)) dyio_outchr(dy_logchn,dy_gtxecho,'*') ;
      dyio_outchr(dy_logchn,dy_gtxecho,' ') ; }
    if (flgon(*checks,ladPRIMALCHK))
      dyio_outfmt(dy_logchn,dy_gtxecho,"pchk ") ;
    if (flgon(*checks,ladDUALCHK))
      dyio_outfmt(dy_logchn,dy_gtxecho,"dchk ") ;
    if (flgon(*checks,ladPRIMFEAS))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"pfeas") ;
      if (flgon(*checks,ladPFQUIET)) dyio_outchr(dy_logchn,dy_gtxecho,'q') ;
      dyio_outchr(dy_logchn,dy_gtxecho,' ') ; }
    if (flgon(*checks,ladDUALFEAS))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"dfeas") ;
      if (flgon(*checks,ladDFQUIET)) dyio_outchr(dy_logchn,dy_gtxecho,'q') ;
      dyio_outchr(dy_logchn,dy_gtxecho,' ') ; } }
  
  if (flgon(*checks,ladDUALCHK))
    dualerrs = (double *) MALLOC((dy_sys->concnt+1)*sizeof(double)) ;
  if (flgon(*checks,ladDUALFEAS))
    dfeaserrs = (double *) CALLOC((dy_sys->varcnt+1),sizeof(double)) ;

  pfeascnt = -1 ;
  pinfeas = quiet_nan(0) ;
# endif

  if (dy_lp->degen > 0 && dy_lp->phase == dyDUAL)
  { dualDegen = TRUE ; }
  else
  { dualDegen = FALSE ; }

/*
  Grab a work vector to do the primal accuracy check, and do a little prep
  work. We need to recalculate the primal variables if the antidegeneracy
  mechanism is active.
*/
  primalerrs = (double *) MALLOC((dy_sys->concnt+1)*sizeof(double)) ;
  tryagain = TRUE ;
  dorefactor = FALSE ;
  results = 0 ;
  pkcol = NULL ;
  refactorcnt = 0 ;

  factorresult = dyrINV ;
  if (flgon(*checks,ladFACTOR))
  { dorefactor = TRUE ; }
  else
  { if (dy_lp->degen > 0 &&
	(dy_lp->phase == dyPRIMAL1 || dy_lp->phase == dyPRIMAL2))
    { if (dy_calcprimals() == FALSE)
      { errmsg(316,rtnnme,dy_sys->nme) ;
	return (dyrFATAL) ; } } }
/*
  Open the check loop and do the tests. Each time we fail a test, we'll
  come back here, tighten the pivot selection tolerance, and try again.
  If we try to tighten, and we're already at the top, then we fail.
*/
  while (tryagain == TRUE)
  { tryagain = FALSE ;
/*
  Refactor? If the user is forcing the initial refactor, use whatever pivot
  selection regime is currently in force. Additional attempts mean that we've
  failed an accuracy or feasibility check further down. Boost the current
  pivot selection ratio one step.  If we're already at the tightest setting,
  dy_setpivparms will return FALSE and we'll escape the loop.
*/
    if (dorefactor == TRUE)
    { if (refactorcnt != 0 || flgoff(*checks,ladFACTOR))
      { if (dy_setpivparms(+1,0) == FALSE) continue ; }
#     ifndef DYLP_NDEBUG
      if (print >= 2 || (print >= 1 && results != 0))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      [%s] refactoring at (%s)%d, %s, ",
		    dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		    dy_lp->tot.iters,dy_prtpivparms(-1)) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"%d pivots since last refactor.",
		    dy_lp->basis.etas) ; }
#     endif
      factorflags = ladPRIMALS|ladDUALS ;
      if (refactorcnt == 0 && flgon(*checks,ladEXPAND))
	setflg(factorflags,ladEXPAND) ;
      factorresult = dy_factor(&factorflags) ;
      refactorcnt++ ;
      if (!(factorresult == dyrOK || factorresult == dyrPATCHED))
      {
#	ifndef DYLP_NDEBUG
	if ((print >= 2) ||
	    (print >= 1 && flgon(results,ladPRIMALCHK|ladDUALCHK)))
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "%sfailed.",(print >= 6)?"\n\t":" ") ;
#	endif
	return (factorresult) ; }
#     ifndef DYLP_NDEBUG
      if ((print >= 2) ||
	  (print >= 1 && flgon(results,ladPRIMALCHK|ladDUALCHK)))
	dyio_outfmt(dy_logchn,dy_gtxecho,"%s%s.",(print >= 6)?"\n\t":" ",
	       (factorresult == dyrOK)?"done":"patched") ;
#     endif
      setflg(results,ladFACTOR|factorflags) ; }
/*
  Do the accuracy checks. Calculate (b - x<N>) - Bx<B> and/or c<B> - yB,
  as requested, as well as the relevant 1-norms and residuals.

  By definition, when dual degeneracy is active the real value of the duals
  involved in the degeneracy is zero.
*/
    normb = 0 ;
    primalresid = 0 ;
    normc = 0 ;
    dualresid = 0 ;
    clrflg(results,ladPRIMALCHK|ladPRIMFEAS|ladDUALCHK|ladDUALFEAS) ;

    if (dy_reducerhs(primalerrs,TRUE) != TRUE)
    { errmsg(340,rtnnme,dy_sys->nme) ;
      return (dyrFATAL) ; }

    if (flgon(*checks,ladPRIMALCHK|ladDUALCHK))
    { for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      { xkndx = dy_basis[xkpos] ;
	if (consys_getcol_pk(dy_sys,xkndx,&pkcol) == FALSE)
	{ errmsg(112,rtnnme,dy_sys->nme,"retrieve","column",
		 consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx) ;
	return (dyrFATAL) ; }

	if (flgon(*checks,ladPRIMALCHK))
	{ normb += fabs(dy_sys->rhs[xkpos]) ;
	  for (pkndx = 0 ; pkndx < pkcol->cnt ; pkndx++)
	  { cndx = pkcol->coeffs[pkndx].ndx ;
	    primalerrs[cndx] -= pkcol->coeffs[pkndx].val*dy_x[xkndx] ; } }

	if (flgon(*checks,ladDUALCHK))
	{ cbarj = dy_sys->obj[xkndx] ;
	  normc += fabs(cbarj) ;
	  for (pkndx = 0 ; pkndx < pkcol->cnt ; pkndx++)
	  { cndx = pkcol->coeffs[pkndx].ndx ;
	    if (dualDegen == TRUE && dy_ddegenset[cndx] > 0)
	    { /* do nothing --- dy_y[cndx] == 0 */ }
	    else
	    { cbarj -= pkcol->coeffs[pkndx].val*dy_y[cndx] ; } }
	  dualresid += fabs(cbarj) ;
#	  ifndef DYLP_NDEBUG
	  dualerrs[xkpos] = fabs(cbarj) ;
#	  endif
	} }
    
      if (flgon(*checks,ladPRIMALCHK))
	for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
	  primalresid += fabs(primalerrs[xkpos]) ;
/*
  And the results?
*/
      normb = maxx(normb,dy_lp->prim.norm1) ;
      if (!withintol(primalresid,0.0,dy_tols->pchk*(1+normb)))
      { warn(341,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"primal",'b',normb,primalresid,
	     dy_tols->pchk*(1+normb)) ;
	setflg(results,ladPRIMALCHK) ; }
      
      normc = maxx(normc,dy_lp->dual.norm1) ;
      normc *= 2 ;
      if (!withintol(dualresid,0.0,dy_tols->dchk*(1+normc)))
      { warn(341,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"dual",'c',normc,dualresid,
	     dy_tols->dchk*(1+normc)) ;
	setflg(results,ladDUALCHK) ; } }

#   ifndef DYLP_NDEBUG
/*
  Information prints on the results of the accuracy checks.
*/
    if ((print >= 3 && flgon(results,ladPRIMALCHK)) ||
	(print >= 3 && flgon(*checks, ladPRIMALCHK)))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s primal accuracy check, ",
		  (flgon(results,ladPRIMALCHK))?"failed":"passed") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "residual/(1+normb) = %g/%g = %g %c %g.",
		  primalresid,(1+normb),primalresid/(1+normb),
		  (flgon(results,ladPRIMALCHK))?'>':'<',dy_tols->pchk) ;
      dyio_outfmt(dy_logchn,dy_gtxecho," (%.2f%%)",
		  primalresid/(dy_tols->pchk*(1+normb))*100) ; }

    if ((print >= 3 && flgon(results,ladDUALCHK)) ||
	(print >= 3 && flgon(*checks, ladDUALCHK)))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s dual accuracy check, ",
		  (flgon(results,ladDUALCHK))?"failed":"passed") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "residual/(1+normc) = %g/%g = %g %c %g.",
		  dualresid,(1+normc),dualresid/(1+normc),
		  (flgon(results,ladDUALCHK))?'>':'<',dy_tols->dchk) ;
      dyio_outfmt(dy_logchn,dy_gtxecho," (%.2f%%)",
		  dualresid/(dy_tols->dchk*(1+normc))*100) ; }

    if (print >= 3)
    { if (flgon(results,ladPRIMALCHK))
      { avgerr = (dy_tols->pchk*(1+normb))/dy_sys->concnt ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    rows exceeding scaled average tolerance %g:",
		    avgerr) ;
	for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
	  if (fabs(primalerrs[xkpos]) > avgerr)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\trow %s (%d), residual %g.",
		        consys_nme(dy_sys,'c',xkpos,FALSE,NULL),xkpos,
		        primalerrs[xkpos]) ; } }

      if (flgon(results,ladDUALCHK))
      { avgerr = (dy_tols->dchk*(1+normc))/dy_sys->concnt ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    columns exceeding scaled average tolerance %g:",
		    avgerr) ;
	for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
	  if (fabs(dualerrs[xkpos]) > avgerr)
	  { xkndx = dy_basis[xkpos] ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,
			"\n\tpos'n %d, column %s (%d), residual %g.",
		        xkpos,consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		        dualerrs[xkpos]) ; } } }
#   endif

/*
  If we didn't pass the accuracy checks, go back and try refactoring to see
  if it will improve the accuracy.
*/
    if (flgon(results,ladPRIMALCHK|ladDUALCHK))
    { tryagain = TRUE ;
      dorefactor = TRUE ;
      continue ; }
  
/*
  We've passed the accuracy check; if we did a refactor, groom the basis.
*/
    if (factorresult != dyrINV)
    { if (groombasis() == FALSE)
      { errmsg(373,rtnnme,dy_sys->nme) ;
	return (dyrFATAL) ; } }
/*
  The primal feasibility check. There are two entirely separate iterations
  here. The nonparanoid one checks the basic variables by iterating over the
  basis. The paranoid one checks all variables. Since we don't need the column
  for this check, separating it out seems a good idea. If we fail the check,
  go up and try to refactor.
*/
    if (flgon(*checks,ladPRIMFEAS))
    { pinfeas = 0 ;
      pfeascnt = 0 ;
#     ifdef PARANOIA
      for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
      { if (!withinbnds(dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx]))
	{ if (belowbnd(dy_x[xkndx],dy_sys->vlb[xkndx]))
	  { pinfeas += dy_sys->vlb[xkndx]-dy_x[xkndx] ;
	    if (flgoff(*checks,ladPFQUIET))
	      warn(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(dy_status[xkndx]),
		   dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx],
		   dy_sys->vlb[xkndx]-dy_x[xkndx],
		   dy_tols->pfeas*(1+fabs(dy_sys->vlb[xkndx]))) ; }
	  else
	  { pinfeas += dy_x[xkndx]-dy_sys->vub[xkndx] ;
	    if (flgoff(*checks,ladPFQUIET))
	      warn(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(dy_status[xkndx]),
		   dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx],
		   dy_x[xkndx]-dy_sys->vub[xkndx],
		   dy_tols->pfeas*(1+fabs(dy_sys->vub[xkndx]))) ; }
	  pfeascnt++ ;
	  setflg(results,ladPRIMFEAS) ; } }
#     else
      for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      { xkndx = dy_basis[xkpos] ;
	if (!withinbnds(dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx]))
	{ if (belowbnd(dy_x[xkndx],dy_sys->vlb[xkndx]))
	  { pinfeas += dy_sys->vlb[xkndx]-dy_x[xkndx] ;
	    if (flgoff(*checks,ladPFQUIET))
	      warn(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(dy_status[xkndx]),
		   dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx],
		   dy_sys->vlb[xkndx]-dy_x[xkndx],
		   dy_tols->pfeas*(1+fabs(dy_sys->vlb[xkndx]))) ; }
	  else
	  { pinfeas += dy_x[xkndx]-dy_sys->vub[xkndx] ;
	    if (flgoff(*checks,ladPFQUIET))
	      warn(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(dy_status[xkndx]),
		   dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx],
		   dy_x[xkndx]-dy_sys->vub[xkndx],
		   dy_tols->pfeas*(1+fabs(dy_sys->vub[xkndx]))) ; }
	  pfeascnt++ ;
	  setflg(results,ladPRIMFEAS) ; } }
#   endif
      dy_lp->infeas = pinfeas ;
      dy_lp->infeascnt = pfeascnt ; }

#   ifndef DYLP_NDEBUG
/*
  Information prints on the results of the primal feasibility check.
*/
    if (print >= 5)
    { if (flgon(*checks,ladPRIMFEAS) && flgoff(results,ladPRIMFEAS))
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tpassed primal feasibility check.") ; }

    if (print >= 3)
    { if (flgon(results,ladPRIMFEAS))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    (%s)%d: %d variables primal infeasible, pinfeas = %g:",
		dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		pfeascnt,pinfeas) ;
	for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
	  if (!withinbnds(dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx]))
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		        "\n\t%s (%d) = %g, status %s, lb = %g, ub = %g,",
		        consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
			dy_x[xkndx],dy_prtvstat(dy_status[xkndx]),
			dy_sys->vlb[xkndx],dy_sys->vub[xkndx]) ;
	    if (dy_x[xkndx] < dy_sys->vlb[xkndx])
	      dyio_outfmt(dy_logchn,dy_gtxecho,"lb-x = %g, tol = %g",
			  dy_sys->vlb[xkndx]-dy_x[xkndx],
			  dy_tols->pfeas*(1+fabs(dy_sys->vlb[xkndx]))) ;
	    else
	      dyio_outfmt(dy_logchn,dy_gtxecho,"x-ub = %g, tol = %g",
			  dy_x[xkndx]-dy_sys->vub[xkndx],
			  dy_tols->pfeas*(1+fabs(dy_sys->vub[xkndx]))) ; } } }
#   endif

    if (flgon(results,ladPRIMFEAS) && flgoff(*checks,ladPFQUIET))
    { tryagain = TRUE ;
      dorefactor = TRUE ;
      continue ; }
/*
  The dual feasibility check. Also separate, as we're interested in checking
  the nonbasic columns.

  If we're running a perturbed, restricted subproblem for dual
  antidegeneracy, the action deserves some explanation. Even though perturbed
  cbar values are propagated to duals during iterative update, the original
  perturbation was applied directly to cbar, so in general a perturbed
  cbar<j> != c<j> - dot(a<j>,y). Still, dual feasibility should be preserved
  for the perturbed values. We just need to be careful not to erase the
  perturbed cbar.
*/
    if (flgon(*checks,ladDUALFEAS))
    {
#     ifndef DYLP_NDEBUG
      dfeascnt = 0 ;
      dinfeas = 0 ;
#     endif
      for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
      { vstat = dy_status[xkndx] ;
#       ifdef PARANOIA
        if (flgon(vstat,vstatBASIC)) continue ;
#       else
        if (flgon(vstat,vstatBASIC|vstatNBFX)) continue ;
#       endif
        cbarj = consys_dotcol(dy_sys,xkndx,dy_y) ;
#       ifdef PARANOIA
	if (isnan(cbarj) == TRUE)
	{ errmsg(320,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,"y",xkndx,"dual feasibility check") ;
	  return (dyrFATAL) ; }
#       endif
	if (dualDegen == TRUE && dy_ddegenset[xkndx] > 0)
	{ cbarj = dy_cbar[xkndx] ; }
	else
	{ cbarj = dy_sys->obj[xkndx]-cbarj ;
	  setcleanzero(cbarj,dy_tols->cost) ;
	  dy_cbar[xkndx] = cbarj ; }
	if ((flgon(vstat,vstatNBLB) && cbarj < -dy_tols->dfeas) ||
	    (flgon(vstat,vstatNBUB) && cbarj > dy_tols->dfeas) ||
	    flgon(vstat,vstatNBFR|vstatSB))
	{ if (flgoff(*checks,ladDFQUIET))
	    warn(347,rtnnme,dy_sys->nme,
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		 consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		 dy_prtvstat(vstat),xkndx,cbarj,dy_tols->dfeas) ;
#	  ifndef DYLP_NDEBUG
/*
  Note that we depend on CALLOC to set dfeaserrs = 0 on allocation, so that
  nonzero entries represent dual feasibility violations.
*/
	  dinfeas += fabs(cbarj) ;
	  dfeascnt++ ;
	  dfeaserrs[xkndx] = cbarj ;
#	  endif
	  setflg(results,ladDUALFEAS) ; } }
	if (flgon(results,ladDUALFEAS))
	  dy_lp->basis.dinf++ ;
	else
	  dy_lp->basis.dinf = 0 ; }

# ifndef DYLP_NDEBUG
/*
  Information prints on the results of the dual feasibility check.
*/
    if (print >= 5)
    { if (flgon(*checks,ladDUALFEAS) && flgoff(results,ladDUALFEAS))
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tpassed dual feasibility check.") ; }

    if (print >= 3)
    { if (flgon(results,ladDUALFEAS))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    (%s)%d: %d variables dual infeasible, dinfeas = %g:",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dfeascnt,dinfeas) ;
	for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
	{ if (dfeaserrs[xkndx] != 0.0)
	  { vstat = dy_status[xkndx] ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t%s (%d) = %g, status %s, cbarj = %g, tol %g.",
		    consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		    (vstat == vstatNBLB)?dy_sys->vlb[xkndx]:dy_sys->vub[xkndx],
		    dy_prtvstat(vstat),dfeaserrs[xkndx],
		    dy_tols->dfeas) ; } } } }
# endif

    if (flgon(results,ladDUALFEAS) && flgoff(*checks,ladDFQUIET))
    { tryagain = TRUE ;
      dorefactor = TRUE ;
      continue ; } }
/*
  End of the accuracy check loop, for good or ill. If we've managed to pass
  all checks on the first try, with no problems, back off the pivot selection
  tolerances. (This implies we actually performed accuracy checks.) Free the
  work vectors, set the return value for checks and we're done.
*/
  if (refactorcnt == 0 && flgon(*checks,ladPRIMALCHK|ladDUALCHK))
    (void) dy_setpivparms(-1,0) ;

  FREE(primalerrs) ;
  if (pkcol != NULL) pkvec_free(pkcol) ;
# ifndef DYLP_NDEBUG
  if (flgon(*checks,ladDUALCHK)) FREE(dualerrs) ;
  if (flgon(*checks,ladDUALFEAS)) FREE(dfeaserrs) ;
# endif

  *checks = results ;

  return (dyrOK) ; }



dyret_enum dy_duenna (dyret_enum pivresult, int xjndx, int xindx,
		      int xjcand, int xicand)

/*
  This routine is the Duenna for dylp's primal and dual simplex algorithms.
  It is called after each pivot, checks that all the right things are being
  done, and deals with major scandal as best it can.

  The `right things' are the mundane precautions of regular accuracy checks
  and refactorisation. These are applied if the pivot result is dyrOK or
  dyrDEGEN. (If the pivot result is dyrOPTIMAL, the simplex routines will
  do their own preoptimality checks.)

  Pivoting scandals understood by the Duenna, and possible actions taken
  here, include:
    * Unbounded problem -- action depends on whether we have primal I, II,
      or dual unboundedness.
      + Primal I unbounded -- flag the variable as NOPIVOT, and force a
	reselect. What can happen here is we get moving in a space of free
	variables that's orthogonal to feasibility.
      + Primal II unbounded -- no remedy here, but in dylp this will often
	occur in the early stages of solving an LP, when only a subset of the
	constraints are active, and is handled by activating more constraints.
      + Dual unbounded -- indicates primal infeasibility, and no action is
	taken here.
    * Loss of feasibility -- really shouldn't occur, and indicates loss of
      accuracy or algorithmic error.
      + primal phase I -- serious internal confusion, converted to a fatal
	error.
      + primal phase II -- refactor is forced, in anticipation that
	it was an accuracy problem. If infeasibility remains, return loss of
	feasibility so we can revert to primal phase I.
      + dual phase II -- refactor is forced, in anticipation that it was
	an accuracy problem. If infeasibility remains, return loss of
	feasibility in hopes dy_dual will do something (currently it punts
	to primal I).
    * Suspected loss of accuracy (dyrREQCHK) -- force a refactor.
    * Mad pivot -- the pivoting routines may reject a pivot if they judge it
      (numerically) unstable. The pivot column (row) is added to the list of
      rejected pivots, and la Duenna returns dyrRESELECT to request that the
      primal (dual) simplex algorithm reselect the entering (leaving) variable.

  Two pivoting problems originate with dy_pivot and the underlying basis
  package:
    * Singular basis -- this occurs when the basis package attempted the
      recommended pivot and discovered that it resulted in a singular basis.
      The basic remedy is to refactor. It could be the basis was already
      singular, but accumulated inaccuracy masked it until now.  dy_factor
      and the basis package will patch it unless the user has disabled that
      option. Or it could be that this choice of pivot is indeed the problem.
      If the refactor completes without problem, blame the pivot and mark the
      column (row) with the NOPIVOT qualifier. In any event, successful
      recovery implies that la Duenna should request the primal (dual) simplex
      algorithm to reselect.
    * Out of space -- this occurs when the basis package ran out of space.
      The remedy is to refactor, recovering the space occupied by eta
      matrices. This will trigger the allocation of additional space, if
      that's what's required.

  When dy_pivot attempts a pivot and fails the basis is corrupted and
  refactoring is mandatory.

  A pivresult of dyrFATAL is fatal, period.

  La Duenna detects a three fatal conditions based on pivot counts:
    * Iteration limit exceeded --  no remedy, dylp has done as much work as 
      its client authorised.
    * Stalling (possibly actual cycling) -- no remedy. This is a heuristic,
      based on the number of pivots that have occurred with no improvement in
      the objective function.
    * Too many rejected pivots -- no remedy, there have been too many
      consecutive rejected pivots.

  Refactoring is handled inside dy_accchk, which see for additional
  comments.  Where refactoring is attempted (either for error recovery or
  simply because basis.etas says it's time) and dy_factor reports that it has
  discovered and patched a singular basis, la Duenna will return dyrRESELECT
  to force the calling simplex algorithm to reselect the entering (leaving)
  variable.  Refactoring can fail; again, see dy_accchk and dy_factor for
  comments. If the basis is refactored, La Duenna will check that the
  preselected entering (leaving) variable is still legal according to primal
  (dual) pivoting rules. If it isn't, dyrRESELECT is returned to request
  selection of a new candidate.

  Pivots where a nonbasic variable swings from one bound to the other don't
  add an eta matrix in the basis representation and are not counted for
  purposes of determining if it's time to do an accuracy check or refactor.

  Parameters:
    pivresult:	The return code from the pivoting routine.
    xjndx:	The entering variable for this pivot.
    xindx:	The leaving variable for this pivot.
    xjcand:	The candidate entering variable for the next pivot (primal)
    xicand:	The candidate entering variable for the next pivot (dual)

  Returns: any of dyrOK, dyrOPTIMAL, dyrRESELECT, dyrPUNT, dyrACCCHK,
	   dyrUNBOUND, dyrSWING, dyrLOSTPFEAS, dyrLOSTDFEAS, dyrSINGULAR,
	   dyrNUMERIC, dyrBSPACE, dyrSTALLED, dyrITERLIM, or dyrFATAL under
	   normal circumstances. dyrINV is possible if a paranoid check fails.
*/

{ dyret_enum retval,accchk ;
  double cbarcand ;
  bool pivok ;
  flags checkflags,outflags,statcand ;
  const char *rtnnme = "La Duenna" ;

# ifndef DYLP_NDEBUG
  int print ;
# endif

# ifndef DYLP_NDEBUG
  switch (dy_lp->phase)
  { case dyPRIMAL1:
    { print = dy_opts->print.phase1 ;
      break ; }
    case dyPRIMAL2:
    { print = dy_opts->print.phase2 ;
      break ; }
    case dyDUAL:
    { print = dy_opts->print.dual ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (dyrFATAL) ; } }
# endif

  retval = dyrINV ;
  outflags = 0 ;
  checkflags = 0 ;
/*
  Bump the various pivot and iteration counts, and see if we're in trouble
  because of the total pivot limit. It's important to get basis.etas and
  basis.pivs correct --- they come into the control of basis factorisation
  and the pivot reject list, and should not change unless we successfully
  pivoted the basis. The others are less critical.
*/
  pivok = dy_lp->pivok ;
  dy_lp->prev_pivok = pivok ;
  dy_lp->pivok = FALSE ;
  if (pivok == TRUE && xindx != xjndx)
  { dy_lp->basis.etas++ ;
    dy_lp->basis.pivs++ ; }
  if (dy_lp->phase == dyPRIMAL1)
  { dy_lp->p1.iters++ ;
    if (pivok == TRUE) dy_lp->p1.pivs++ ;
    if (dy_opts->iterlim > 0 && dy_lp->p1.pivs > dy_opts->iterlim)
      retval = dyrITERLIM ; }
  else
  if (dy_lp->phase == dyPRIMAL2)
  { dy_lp->p2.iters++ ;
    if (pivok == TRUE) dy_lp->p2.pivs++ ;
    if (dy_opts->iterlim > 0 && dy_lp->p2.pivs > dy_opts->iterlim)
      retval = dyrITERLIM ; }
  else
  { dy_lp->d2.iters++ ;
    if (pivok == TRUE) dy_lp->d2.pivs++ ;
    if (dy_opts->iterlim > 0 && dy_lp->d2.pivs > dy_opts->iterlim)
      retval = dyrITERLIM ; }
  if (retval == dyrITERLIM)
  { if (dy_opts->context != cxBANDC)
    { errmsg(328,rtnnme,dy_sys->nme,dy_opts->iterlim) ; }
    return (retval) ; }
  dy_lp->tot.iters++ ;
  if (pivok == TRUE) dy_lp->tot.pivs++ ;
  if (dy_opts->iterlim > 0 && dy_lp->tot.pivs > 3*dy_opts->iterlim)
  { retval = dyrITERLIM ;
    if (dy_opts->context != cxBANDC)
    { errmsg(328,rtnnme,dy_sys->nme,3*dy_opts->iterlim) ; }
    return (retval) ; }
/*
  Deal with the result returned by the pivoting routine. If there's a problem
  we can't fix, we'll return to the caller from whatever case we're in. If
  there's no problem, or the problem can possibly be fixed or papered over,
  we'll fall past the bottom of the switch and do repair and periodic checks.

  Put the cases that are essentially successful pivots (dyrOK, dyrDEGEN,
  and dyrOPTIMAL) first, so we don't waste time.

  If the pivot has been rejected as unstable (dyrMADPIV), it was added to the
  pivot reject list by dy_primalpivot or dy_dualpivot. All we need to do here
  is return dyrRESELECT to trigger selection of a new candidate. Similarly,
  if the pivot routine has directly requested a reselect.

  A punt indicates that the pivot selection routines dualout/dseupdate
  (primalin/pseupdate) couldn't select a leaving (entering) variable, but
  some potentially promising variables were on the pivot rejection list.
  dy_dealWithPunt() will examine the list to see if any variables should be
  released for pivoting and return one of dyrRESELECT, dyrPUNT, or dyrFATAL,
  depending on how it goes.

  There's nothing we can do here about primal unboundedness or pseudo-
  unboundedness, but it may not be a big problem a few levels up (could be
  because only a subset of constraints are active, could be an accuracy
  problem curable by more frequent refactoring). Return quietly. dyrSWING
  typically indicates a successful pivot in which some variable(s) moved a
  little too far.

  If the code has requested an accuracy check (because it saw a bogus number)
  we'll go ahead and do it but won't necessarily reject the current pivot.
  Loss of primal or dual feasibility also indicates accuracy problems. We'll
  refactor and do an accuracy check.

  If the pivot attempt resulted in a singular basis, we definitely need to
  refactor (the present basis representation is corrupt until a successful
  factorization). Assume that the indication of singularity is genuine and
  reject the pivot. If the basis was already singular (and accumulated
  numerical error was masking this) then we've shot the messenger, but the
  code that cleans up after a patch clears the pivot reject list, among other
  actions.

  Out of space entails the same repair actions as singular basis, but without
  the complication of rejecting the pivot element. The notion is that
  refactoring will eliminate all the eta matrices and reclaim space, but if
  we've just refactored, that won't happen, and we need to force expansion.
  In either case, abort the minor iteration, as the pivot attempt has
  failed.

  Finally, fatal is fatal, and so is anything we don't recognise.
*/
  switch (pivresult)
  { case dyrOK:
    case dyrOPTIMAL:
    case dyrDEGEN:
    { retval = dyrOK ;
      break ; }
    case dyrMADPIV:
    case dyrRESELECT:
    { return (dyrRESELECT) ; }
    case dyrPUNT:
    { retval = dy_dealWithPunt() ;
      return (retval) ; }
    case dyrUNBOUND:
    { return (dyrUNBOUND) ; }
    case dyrSWING:
    { return (dyrSWING) ; }
    case dyrREQCHK:
    { setflg(checkflags,ladFACTOR) ;
#     ifndef DYLP_NDEBUG
      if (print >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  [%s] refactor requested at (%s)%d, ",
		    dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		    dy_lp->tot.iters) ;
        dyio_outfmt(dy_logchn,dy_gtxecho,"%d pivots since last refactor.",
		    dy_lp->basis.etas) ; }
#     endif
      retval = dyrREQCHK ;
      break ; } 
    case dyrLOSTPFEAS:
    case dyrLOSTDFEAS:
    { setflg(checkflags,ladFACTOR) ;
#     ifdef PARANOIA
      if (dy_lp->phase == dyPRIMAL1)
      { errmsg(338,rtnnme,dy_sys->nme,
	       (pivresult == dyrLOSTPFEAS)?"primal":"dual") ;
        return (dyrFATAL) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (print >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  [%s] loss of %s feasibility at iteration %d.",
		    dy_sys->nme,(pivresult == dyrLOSTPFEAS)?"primal":"dual",
		    dy_lp->tot.iters) ; }
#     endif
      retval = dyrREQCHK ;
      break ; } 
    case dyrSINGULAR:
    { setflg(checkflags,ladFACTOR) ;
#     ifndef DYLP_NDEBUG
      if (print >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  [%s](%s)%d: pivot attempt produced ",
		    dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		    dy_lp->tot.iters) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "singular basis; attempting recovery.") ; }
#     endif
      if (dy_lp->phase == dyDUAL)
	retval = dy_addtopivrej(xindx,dyrSINGULAR,0.0,0.0) ;
      else
	retval = dy_addtopivrej(xjndx,dyrSINGULAR,0.0,0.0) ;
      if (retval != dyrOK) return (retval) ;
      retval = dyrRESELECT ;
      break ; }
    case dyrBSPACE:
    { setflg(checkflags,ladFACTOR) ;
      if (dy_lp->basis.etas == 0) setflg(checkflags,ladEXPAND) ;
#     ifndef DYLP_NDEBUG
      if (print >= 1)
      { if (flgoff(checkflags,ladEXPAND))
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n  [%s]: (%s)%d: attempting to compress basis.",
		      dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		      dy_lp->tot.iters) ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n  [%s]: (%s)%d: forcing basis expansion.",dy_sys->nme,
		      dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
#     endif
      retval = dyrOK ;
      break ; }
    case dyrFATAL:
    { errmsg(343,rtnnme,dy_sys->nme,dy_lp->tot.iters) ;
      return (dyrFATAL) ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (dyrFATAL) ; } }
# ifdef PARANOIA
/*
  retval should be one of dyrOK, dyrREQCHK, or dyrRESELECT, dyrOK means that
  we may refactor, but have no particular accuracy problems nor do we need to
  reselect if the refactor succeeds. dyrREQCHK should be interpreted to mean
  that we want to refactor and that there is some question of accuracy to
  boot. dyrRESELECT means that we want to refactor and that we'll need to
  reselect the entering/leaving variable once we return to the calling
  simplex algorithm.
*/
  if (!(retval == dyrOK || retval == dyrREQCHK || retval == dyrRESELECT))
  { errmsg(1,rtnnme,__LINE__) ;
    return (dyrFATAL) ; }
# endif
/*
  To get here, either the pivot was uneventful, or there's a hope that
  refactoring will solve whatever problem we ran into.

  Check whether it's time for a regularly scheduled refactorisation or
  accuracy check.  Refactoring always implies an accuracy check (but not vice
  versa), so it's actually buried in dy_accchk, we just have to request it.
  It'll occasionally happen that an accuracy check will repeat because it's
  had the misfortune to fall on a pivot which is followed by pivots which
  don't count toward basis.etas (nonbasic moving from bound to bound is the
  most common cause). Something to think about fixing up --- either a a
  separate count, iterc, for accuracy checks, or remember the basis.etas
  value of the most recent accuracy check.

  While we're here, check on the pivot tolerance. If we've been running with
  reduced tolerance for a while, we should boost it back to normal.
*/
  if (dy_lp->basis.etas >= dy_opts->factor)
  { setflg(checkflags,ladFACTOR) ;
    dy_checkpivtol() ;
#   ifndef DYLP_NDEBUG
    if (print >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    [%s] (%s)%d: scheduled refactor, interval %d, ",
		  dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dy_opts->factor) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%d pivots since last refactor,",
		  dy_lp->basis.etas) ;
      if (dy_lp->phase == dyDUAL)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"yb = %g.",dy_calcdualobj()) ; }
      else
      if (dy_lp->phase == dyPRIMAL1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"infeas = %g.",dy_calcpinfeas()) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"cx = %g.",dy_calcobj()) ; } }
#   endif
  }
  if ((dy_lp->basis.etas%dy_opts->check == 0 && dy_lp->basis.etas != 0) ||
      flgon(checkflags,ladFACTOR))
  { if (dy_lp->phase == dyPRIMAL1)
    { setflg(checkflags,ladPRIMALCHK) ; }
    else
    if (dy_lp->phase == dyPRIMAL2)
    { setflg(checkflags,ladPRIMALCHK|ladPRIMFEAS|ladPFQUIET|ladDUALCHK) ; }
    else
    if (dy_lp->phase == dyDUAL)
    { setflg(checkflags,ladPRIMALCHK|ladDUALCHK|ladDUALFEAS|ladDFQUIET) ; }
#   ifndef DYLP_NDEBUG
    if (dy_lp->basis.etas%dy_opts->check == 0 &&
	dy_lp->basis.etas != 0 && print >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    [%s] (%s)%d: scheduled check, interval %d, ",
		  dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dy_opts->check) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%d pivots since last refactor.",
		  dy_lp->basis.etas) ; }
#   endif
  }
/*
  Do an accuracy check, possibly accompanied by a refactorisation. Any code
  other than dyrOK or dyrPATCHED is a problem; return it to the caller. If
  we've failed any of the accuracy checks, return dyrACCCHK, since dy_accchk
  has tried everthing dylp knows for dealing with the problem.

  Assuming we pass the accuracy check, clear the pivot rejection list.
  Limiting this to the accuracy check frequency gives a little bit of
  persistence to the list. But ... don't clear the list if the reason we're
  here is because the last pivot resulted in a singular basis!
*/
  if (flgon(checkflags,ladFACTOR|ladPRIMALCHK|ladDUALCHK))
  { accchk = dy_accchk(&checkflags) ;
    if (!(accchk == dyrOK || accchk == dyrPATCHED)) return (accchk) ;
    if (flgon(checkflags,ladPRIMALCHK|ladDUALCHK)) return (dyrACCCHK) ;
    if (pivresult != dyrSINGULAR)
      if (dy_clrpivrej(NULL) != TRUE) return (dyrFATAL) ;
/*
  If the return code is dyrPATCHED, we'll want to select a new pivot
  candidate before attempting another simplex iteration. If we refactored,
  check that the candidate entering (leaving) variable still meets primal
  (dual) criteria.
*/
    if (flgon(checkflags,ladFACTOR) && retval != dyrRESELECT)
    { if (accchk == dyrPATCHED)
      { retval = dyrRESELECT ;
#       ifndef DYLP_NDEBUG
	if (print >= 1)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n    [%s] (%s)%d: ",
		      dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		      dy_lp->tot.iters) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "forcing reselect after basis patch.") ; }
#       endif
      }
      else
      { if ((dy_lp->phase == dyPRIMAL1 ||
	     dy_lp->phase == dyPRIMAL2) && xjcand > 0)
	{ statcand = dy_status[xjcand] ;
	  cbarcand = dy_cbar[xjcand] ;
	  if (!((cbarcand+dy_tols->dfeas < 0 &&
		 flgon(statcand,vstatNBLB|vstatNBFR|vstatSB)) ||
		(cbarcand-dy_tols->dfeas > 0 &&
		 flgon(statcand,vstatNBUB|vstatNBFR|vstatSB))))
	  { retval = dyrRESELECT ; } } 
	else
	if (dy_lp->phase == dyDUAL && xicand > 0)
	{ statcand = dy_status[xicand] ;
	  if (flgoff(statcand,vstatBLLB|vstatBUUB))
	  { retval = dyrRESELECT ; } }
#	ifndef DYLP_NDEBUG
	if (print >= 1 && retval == dyrRESELECT)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n    [%s] (%s)%d: ",
		      dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		      dy_lp->tot.iters) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		"candidate %s (%d) no longer suitable; ",
		consys_nme(dy_sys,'v',
			   (dy_lp->phase == dyDUAL)?xicand:xjcand,FALSE,NULL),
			   (dy_lp->phase == dyDUAL)?xicand:xjcand) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "forcing reselect after refactor.") ; }
#       endif
	} }
/*
  What's the story on feasibility? Check the possible outcomes on the
  accuracy test. If we're here because an accuracy check was requested, or
  because we lost feasibility, and all is now ok, set the return value to
  dyrOK. If we've lost primal or dual feasibility, the proper one to report
  depends on the phase (proper, that is, in the sense of which one dominates
  when we've lost both).
*/
    if (flgon(checkflags,ladPRIMFEAS|ladDUALFEAS))
    { if (dy_lp->phase == dyDUAL)
      { if (flgon(checkflags,ladDUALFEAS))
	  return (dyrLOSTDFEAS) ;
	else
	if (flgon(checkflags,ladPRIMFEAS))
	  return (dyrLOSTPFEAS) ; }
      else
      { if (flgon(checkflags,ladPRIMFEAS))
	  return (dyrLOSTPFEAS) ;
	else
	if (flgon(checkflags,ladDUALFEAS))
	  return (dyrLOSTDFEAS) ; } }
    if (retval == dyrREQCHK) retval = dyrOK ; }
# ifdef PARANOIA
/*
  At this point, return values should be one of dyrOK or dyrRESELECT.
*/
  if (!(retval == dyrOK || retval == dyrRESELECT))
  { errmsg(1,rtnnme,__LINE__) ;
    return (dyrFATAL) ; }
# endif
/*
  Not much left. Now that we're as accurate as we're going to get, see if the
  objective's been improving. If it hasn't, we're stalled (and might be
  cycling).

  Note that the cycling count is not 100% precise --- in particular, La Duenna
  does not reach here for the case of a successful pivot flagged dyrPUNT.
*/
  if (pivok == TRUE)
  { if (withintol(dy_lp->z,dy_lp->lastz.piv,dy_tols->dchk))
    { dy_lp->idlecnt++ ;
      if (dy_lp->idlecnt > dy_opts->idlelim)
      { errmsg(339,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,dy_lp->idlecnt) ;
	return (dyrSTALLED) ; } }
    else
    { dy_lp->lastz.piv = dy_lp->z ;
      dy_lp->idlecnt = 0 ; } }
/*
  If, after all the above, we're still at dyrOK, and the
  original pivresult was dyrOPTIMAL, go with that.
*/
  if (pivresult == dyrOPTIMAL && retval == dyrOK) retval = dyrOPTIMAL ;
/*
  That's it.
*/
  return (retval) ; }


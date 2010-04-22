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
  This file contains the routines that handle a warm start, given a basis &
  status passed in by the user.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_warmstart.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_warmstart.c 269 2009-04-02 05:38:19Z lou $" ;



static void correct_for_patch (void)

/*
  This routine scans dy_status looking for architectural variables that are
  recorded as basic but have been booted out of the basis by a patch
  operation. It's a very special-purpose routine, separated out so it doesn't
  clutter up the code in dy_warmstart.
  Parameters: none

  Returns: undefined
*/

{ int j,cnt ;
  flags statj ;
  double *vlb,*vub ;

  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
/*
  Open a loop to scan the status array, checking that variables recorded as
  basic are really basic. dy_patch clears the var2basis entry when it makes
  the patch, so we're looking for basic status with a 0 in var2basis.

  When we find a variable that needs to be corrected, decide an appropriate
  nonbasic status based on the sign of the objective coefficient and the
  presence/absence of finite bounds.
*/

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\tcorrecting status due to basis patch ...") ; }
# endif
  
  cnt = 0 ;
  for (j = dy_sys->concnt+1 ; j <= dy_sys->varcnt ; j++)
  { statj = dy_status[j] ;
    if (flgon(statj,vstatBASIC) && dy_var2basis[j] == 0)
    { if (vlb[j] > -dy_tols->inf && vub[j] < dy_tols->inf)
      { if (vub[j] == vlb[j])
	{ dy_status[j] = vstatNBFX ;
	  dy_x[j] = vub[j] ; }
	else
	if (dy_sys->obj[j] >= 0)
	{ dy_status[j] = vstatNBLB ;
	  dy_x[j] = vlb[j] ; }
	else
	{ dy_status[j] = vstatNBUB ;
	  dy_x[j] = vub[j] ; } }
      else
      if (vlb[j] > -dy_tols->inf)
      { dy_status[j] = vstatNBLB ;
	dy_x[j] = vlb[j] ; }
      else
      if (vub[j] < dy_tols->inf)
      { dy_status[j] = vstatNBUB ;
	dy_x[j] = vub[j] ; }
      else
      { dy_status[j] = vstatNBFR ;
	dy_x[j] = 0 ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.crash >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t  changing status for %s (%d) to %s,",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		    dy_prtvstat(dy_status[j])) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," value %g.",dy_x[j]) ; }
#     endif
      cnt++ ; } }

# ifndef DYLP_NDEBUG
/*
  Given that this routine has been called, there should be corrections to be
  made, but it's possible that the patch involved only logicals. If so,
  dy_warmstart has already dealt with the problem and we simply can't tell.
  (The necessary data structure is not exported from dy_basis.c) The least we
  can do is print a message.
*/
  if (cnt == 0 && dy_opts->print.crash >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  no architecturals corrected.") ; }
# endif

  return ; }



dyret_enum dy_warmstart (lpprob_struct *orig_lp)

/*
  This routine is responsible for recreating the active constraint system,
  basis, and status specified by the user in orig_lp. It will handle even the
  pathological case of 0 active constraints and 0 active variables. If the
  user has supplied an active variable vector, only those variables will be
  activated. Clearly, the supplied basis, status, and active variable vector
  should be consistent, or bad things will happen.

  If we're operating in fullsys mode, we need to check here for additions to
  the constraint system.

  << In the very near future, this routine should also be upgraded to cope
     with the possibility that constraints specified in the warm start basis
     have disappeared. >>

  Parameters:
    orig_lp:	The original lp problem structure

  Returns: dyrOK if the setup completes without error, any of a number of
	   error codes otherwise (dyrFATAL, dyrINV, or a code from dy_factor)
*/

{ int vndx,dyvndx,bpos,cndx,dycndx,dycsze,dyvsze,nbfxcnt ;
  double *vlb,*vub,vlbj,vubj,obj ;
  consys_struct *orig_sys ;
  flags *orig_status,vstat,calcflgs ;
  dyret_enum retval ;
  basisel_struct *orig_basis ;
  bool *orig_actvars,rngseen,noactvarspec ;
  pkvec_struct *pkcol ;
  char nmebuf[50] ;

  flags parts = CONSYS_OBJ|CONSYS_VUB|CONSYS_VLB|CONSYS_RHS|CONSYS_RHSLOW|
		CONSYS_VTYP|CONSYS_CTYP,
	opts = CONSYS_LVARS|CONSYS_WRNATT ;
  
  const char *rtnnme = "dy_warmstart" ;

  extern void dy_setfinalstatus(void) ;		/* dy_hotstart.c */

# if defined(DYLP_PARANOIA) || !defined(DYLP_NDEBUG)
  double xi ;
# endif

  retval = dyrINV ;
  nbfxcnt = -1 ;

/*
  Do a little unpacking.
*/
  orig_sys = orig_lp->consys ;
  orig_status = orig_lp->status ;
  orig_basis = orig_lp->basis->el ;
  if (flgon(orig_lp->ctlopts,lpctlACTVARSIN) && dy_opts->fullsys == FALSE)
  { orig_actvars = orig_lp->actvars ;
    noactvarspec = FALSE ; }
  else
  { orig_actvars = NULL ;
    noactvarspec = TRUE ; }
/*
  Initialise the statistics on loadable/unloadable variables and constraints.
*/
  dy_lp->sys.forcedfull = FALSE ;
  dy_lp->sys.vars.loadable = orig_sys->varcnt ;
  dy_lp->sys.vars.unloadable = 0 ;
  dy_lp->sys.cons.loadable = orig_sys->concnt ;
  dy_lp->sys.cons.unloadable = 0 ;
/*
  Create the dy_sys constraint system to match the user's basis and active
  variables (if specified). We'll create the system with logicals enabled.
  
  For variables, if there is an active variable vector, skim it for a count.
  Otherwise, skim the status array and count the number of nonbasic fixed
  variables (which will never become active).

  For constraints, we need to consider the possibility that the user has
  added cuts and is trusting dylp to deal with it. If we're operating in the
  usual dynamic mode, this will be picked up automatically, and we can size
  the constraint system to the active constraints of the basis. But if we're
  operating in fullsys mode, we need to add them here. In this case, the
  number of constraints is the current size of the constraint system.

  Take this opportunity to clean the bounds arrays, making sure that bounds
  within the feasibility tolerance of one another are set to be exactly
  equal.  (This simplifies handling fixed variables.) For nonbasic variables,
  force the status to NBFX and cancel activation if actvars is present. Basic
  variables which need BFX are picked up later, after the basis is
  established.
*/
  vub = orig_sys->vub ;
  vlb = orig_sys->vlb ;
  dyio_outfxd(nmebuf,-((int) (sizeof(nmebuf)-1)),
	      'l',"%s[actv]",orig_sys->nme) ;
  if (noactvarspec == FALSE)
  { dyvsze = 0 ;
    for (vndx = 1 ; vndx <= orig_sys->varcnt ; vndx++)
    { vlbj = vlb[vndx] ;
      vubj = vub[vndx] ;
      if (atbnd(vlbj,vubj))
      { if (vlbj != vubj)
	{ 
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.setup >= 3)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		        "\n\tForcing equal bound %g for %s (%d)",
		        (vlbj+vubj)/2,consys_nme(orig_sys,'v',vndx,0,0),vndx) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,
		        "\n\t  original lb = %g, ub = %g, diff = %g, tol = %g",
		        vlbj,vubj,vubj-vlbj,dy_tols->pfeas) ; }
#	  endif
	  vlb[vndx] = (vlbj+vubj)/2 ;
	  vub[vndx] = vlb[vndx] ; }
	if (((int) orig_status[vndx]) > 0)
	{ orig_status[vndx] = vstatNBFX ;
	  orig_actvars[vndx] = FALSE ; } }
      if (vlb[vndx] > vub[vndx])
      { dy_lp->lpret = lpINFEAS ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n\tTrivial infeasibility for %s (%d), lb = %g > ub = %g.",
	       consys_nme(orig_sys,'v',vndx,0,0),vndx,vlb[vndx],vub[vndx]) ; }
#     endif
      }
      if (orig_actvars[vndx] == TRUE) dyvsze++ ; } }
  else
  { nbfxcnt = 0 ;
    for (vndx = 1 ; vndx <= orig_sys->varcnt ; vndx++)
    { vlbj = vlb[vndx] ;
      vubj = vub[vndx] ;
      if (atbnd(vlbj,vubj))
      { if (vlbj != vubj)
	{ 
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.setup >= 3)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
			"\n\tForcing equal bound %g for %s (g)",
		        (vlbj+vubj)/2,consys_nme(orig_sys,'v',vndx,0,0),vndx) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,
		   "\n\t  original lb = %g, ub = %g, diff = %g, tol = %g",
		   vlbj,vubj,vubj-vlbj,dy_tols->pfeas) ; }
#	  endif
	  vlb[vndx] = (vlbj+vubj)/2 ;
	  vub[vndx] = vlb[vndx] ; }
	if (((int) orig_status[vndx]) > 0)
	{ orig_status[vndx] = vstatNBFX ; } }
      if (vlb[vndx] > vub[vndx])
      { dy_lp->lpret = lpINFEAS ; }
      if ((((int) orig_status[vndx]) > 0) &&
	  flgon(orig_status[vndx],vstatNBFX))
      { nbfxcnt++ ; } }
    dyvsze = orig_sys->varcnt-nbfxcnt ; }
  if (dy_opts->fullsys == TRUE)
    dycsze = orig_sys->concnt ;
  else
    dycsze = orig_lp->basis->len ;
  dyvsze += dycsze ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  creating constraint system %s (%d x %d+%d)",
		nmebuf,dycsze,dyvsze-dycsze,dycsze) ;
    if (dy_opts->print.setup >= 3)
    { if (flgoff(orig_lp->ctlopts,lpctlACTVARSIN))
        dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      %d nonbasic fixed variables excluded.",
		    nbfxcnt) ; } }
# endif
  dy_sys = consys_create(nmebuf,parts,opts,dycsze,dyvsze,dy_tols->inf) ;
  if (dy_sys == NULL)
  { errmsg(152,rtnnme,nmebuf) ;
    return (dyrFATAL) ; }
/*
  Hang a set of translation vectors onto each system: origcons and origvars
  on orig_sys, and actcons and actvars on dy_sys.
*/
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(int),(void **) &dy_actvars) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"active -> original variable map") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(int),(void **) &dy_actcons) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"active -> original constraint map") ;
    return (dyrFATAL) ; }
  if (consys_attach(orig_sys,CONSYS_ROW,
		    sizeof(int),(void **) &dy_origvars) == FALSE)
  { errmsg(100,rtnnme,orig_sys->nme,"original -> active variable map") ;
    return (dyrFATAL) ; }
  if (consys_attach(orig_sys,CONSYS_COL,
		    sizeof(int),(void **) &dy_origcons) == FALSE)
  { errmsg(100,rtnnme,orig_sys->nme,"original -> active constraint map") ;
    return (dyrFATAL) ; }
/*
  dy_origvars is cleared to 0 as it's attached, indicating that the original
  variables have no predefined status. We need to correct this.

  If the caller's supplied an active variable vector, we can use it to
  activate variables prior to adding constraints. (But in any case don't
  activate nonbasic fixed variables.) It's illegal to declare a formerly
  basic variable to be inactive by the simple expedient of setting
  actvars[vndx] = FALSE, hence the paranoid check.

  Otherwise, we'll need to depend on dy_loadcon to activate the variables
  referenced in the active constraints. We'll still fill in origvars, with
  two purposes:
    * We can avoid activating nonbasic fixed variables.
    * We can use dy_origvars == 0 as a paranoid check from here on out.
  Inactive variables are required to be nonbasic, so in this case the proper
  status for formerly basic variables is SB.
*/
  if (noactvarspec == FALSE)
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  processing active variable list ...") ; }
#   endif
    pkcol = pkvec_new(0) ;
    for (vndx = 1 ; vndx <= orig_sys->varcnt ; vndx++)
    { if (((int) orig_status[vndx]) > 0)
	vstat = orig_status[vndx] ;
      else
	vstat = vstatB ;
      if (orig_actvars[vndx] == TRUE && flgoff(vstat,vstatNBFX))
      { if (consys_getcol_pk(orig_sys,vndx,&pkcol) == FALSE)
	{ errmsg(122,rtnnme,orig_sys->nme,"variable",
		 consys_nme(orig_sys,'v',vndx,TRUE,NULL),vndx) ;
	  retval = dyrFATAL ;
	  break ; }
	if (consys_addcol_pk(dy_sys,vartypCON,pkcol,
			     orig_sys->obj[vndx],vlb[vndx],vub[vndx]) == FALSE)
	{ errmsg(156,rtnnme,"variable",dy_sys->nme,pkcol->nme) ;
	  retval = dyrFATAL ;
	  break ; }
	dyvndx = pkcol->ndx ;
	dy_origvars[vndx] = dyvndx ;
	dy_actvars[dyvndx] = vndx ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.setup >= 3)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tactivating %s variable %s (%d) to index %d.",
		      consys_prtvartyp(orig_sys->vtyp[vndx]),
		      consys_nme(orig_sys,'v',vndx,FALSE,NULL),vndx,dyvndx) ; }
#       endif
      }
      else
      {
#       ifdef DYLP_PARANOIA
	if (flgon(vstat,vstatBASIC))
	{ errmsg(380,rtnnme,orig_sys->nme,
		 consys_nme(orig_sys,'v',vndx,FALSE,NULL),vndx,
		 dy_prtvstat(vstat),"non-basic") ;
	  retval = dyrFATAL ;
	  break ; }
#	endif
	dy_origvars[vndx] = -((int) vstat) ; } }
    pkvec_free(pkcol) ;
    if (retval != dyrINV) return (retval) ; }
  else
  { for (vndx = 1 ; vndx <= orig_sys->varcnt ; vndx++)
    { if (((int) orig_status[vndx]) > 0)
	vstat = orig_status[vndx] ;
      else
	vstat = vstatSB ;
      MARK_INACTIVE_VAR(vndx,-((int) vstat)) ; } }
/*
  Walk the basis and install the constraints in order. When we're finished
  with this, the active system will be up and about. In the case where
  there's no active variable specification, some of the status information
  written into dy_origvars may have been overwritten; only variables with
  vstatNBFX are guaranteed to remain inactive.
*/
  rngseen = FALSE ;
  for (bpos = 1 ; bpos <= orig_lp->basis->len ; bpos++)
  { cndx = orig_basis[bpos].cndx ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 2)
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    activating %s %s (%d) in pos'n %d",
		  consys_prtcontyp(orig_sys->ctyp[cndx]),
		  consys_nme(orig_sys,'c',cndx,FALSE,NULL),cndx,bpos) ;
#   endif
#   ifdef DYLP_STATISTICS
    if (dy_stats != NULL) dy_stats->cons.init[cndx] = TRUE ;
#   endif
    if (dy_loadcon(orig_sys,cndx,noactvarspec,NULL) == FALSE)
    { errmsg(430,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "activate","constraint",
	     consys_nme(orig_sys,'c',cndx,TRUE,NULL),cndx) ;
      return (dyrFATAL) ; }
    if (orig_sys->ctyp[cndx] == contypRNG) rngseen = TRUE ; }
/*
  If we're in fullsys mode, repeat constraint installation actions for any
  cuts added after this basis was assembled.
*/
  if (dy_opts->fullsys == TRUE)
  { for (cndx = orig_lp->basis->len+1 ; cndx <= orig_sys->concnt ; cndx++)
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 2)
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    activating %s %s (%d) in pos'n %d",
		    consys_prtcontyp(orig_sys->ctyp[cndx]),
		    consys_nme(orig_sys,'c',cndx,FALSE,NULL),cndx,cndx) ;
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL) dy_stats->cons.init[cndx] = TRUE ;
#     endif
      if (dy_loadcon(orig_sys,cndx,noactvarspec,NULL) == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "activate","constraint",
	       consys_nme(orig_sys,'c',cndx,TRUE,NULL),cndx) ;
	return (dyrFATAL) ; }
      if (orig_sys->ctyp[cndx] == contypRNG) rngseen = TRUE ; } }
# ifdef DYLP_PARANOIA
/*
  Paranoid checks and informational print statements.
*/
  if (dy_chkdysys(orig_sys) == FALSE) return (dyrINV) ;
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    system %s has %d constraints, %d+%d variables",
	        dy_sys->nme,dy_sys->concnt,dy_sys->archvcnt,dy_sys->logvcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
	  "\n    %d constraints, %d variables remain inactive in system %s.",
	  orig_sys->concnt-dy_sys->concnt,orig_sys->archvcnt-dy_sys->archvcnt,
	  orig_sys->nme) ;
    if (dy_opts->print.setup >= 4)
    { nbfxcnt = 0 ;
      for (vndx = 1 ; vndx <= orig_sys->varcnt ; vndx++)
      { if (INACTIVE_VAR(vndx))
	{ vstat = (flags) (-dy_origvars[vndx]) ;
	  switch (getflg(vstat,vstatSTATUS))
	  { case vstatNBUB:
	    { xi = orig_sys->vub[vndx] ;
	      break ; }
	    case vstatNBLB:
	    case vstatNBFX:
	    { xi = orig_sys->vlb[vndx] ;
	      break ; }
	    case vstatNBFR:
	    { xi = 0 ;
	      break ; }
	    default:
	    { errmsg(433,rtnnme,dy_sys->nme,
		     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		     "inactive",consys_nme(orig_sys,'v',vndx,TRUE,NULL),
		     vndx,dy_prtvstat(vstat)) ;
	      return (dyrINV) ; } }
	  if (xi != 0)
	  { if (nbfxcnt == 0)
	      dyio_outfmt(dy_logchn,dy_gtxecho,
			  "\n\tinactive variables with nonzero values:") ;
	    nbfxcnt++ ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) = %g, status %s",
		        consys_nme(orig_sys,'v',vndx,FALSE,NULL),vndx,xi,
		        dy_prtvstat(vstat)) ; } } }
      if (nbfxcnt == 0)
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tall inactive variables are zero.") ; } }
# endif
/*
  Time to assemble the basis. Attach the basis and inverse basis vectors to
  the constraint system. consys_attach will initialise them to 0.
*/
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(int),(void **) &dy_basis) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"basis vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(int),(void **) &dy_var2basis) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"inverse basis vector") ;
    return (dyrFATAL) ; }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 1)
  { if (dy_opts->print.setup == 0)
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n %s: regenerating the basis ...",rtnnme) ;
    else
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  regenerating the basis.",rtnnme) ; }
# endif
/*
  Load the basis. For variables, we need to translate architecturals using
  dy_origvars, and watch out for logicals (vndx = negative of associated
  constraint index). After all the paranoia, we finally update dy_basis and
  dy_var2basis.

  Because we loaded the constraints in the order they were listed in the
  basis, we should have that dycndx = bpos, hence dy_actcons[bpos] = cndx.

  If we're installing a basic variable, it should be active already.  For
  architectural variables, the check is made in dy_origvars.  For a logical,
  the associated constraint should be active, hence a non-zero entry in
  dy_origcons.  For architecturals, we also check if there are any non-zero
  coefficients remaining in the column (who knows what the user has done to
  the constraint system).  This rates a message if the print level is high
  enough, but the basis pacakge is capable of patching the basis. (Indeed,
  it's hard to do it correctly here.) 
*/
# ifdef DYLP_PARANOIA
  pkcol = pkvec_new(0) ;
  retval = dyrOK ;
# endif
  for (bpos = 1 ; bpos <= orig_lp->basis->len ; bpos++)
  { cndx = orig_basis[bpos].cndx ;
    dycndx = dy_origcons[cndx] ;
    vndx = orig_basis[bpos].vndx ;
    if (vndx < 0)
    { dyvndx = dy_origcons[-vndx] ; }
    else
    { dyvndx = dy_origvars[vndx] ; }

#   ifdef DYLP_PARANOIA
    if (dycndx <= 0)
    { errmsg(369,rtnnme,orig_sys->nme,"constraint",
	     consys_nme(orig_sys,'c',cndx,FALSE,NULL),cndx,
	     "cons",cndx,dycndx) ;
      retval = dyrINV ;
      break ; }
    if (dy_actcons[bpos] != cndx)
    { errmsg(370,rtnnme,dy_sys->nme,
	     consys_nme(orig_sys,'c',cndx,FALSE,NULL),cndx,bpos,
	     consys_nme(orig_sys,'c',dy_actcons[bpos],FALSE,NULL),
	     dy_actcons[bpos]) ;
      if (dycndx != bpos) { errmsg(1,rtnnme,__LINE__) ; }
      retval = dyrINV ;
      break ; }

    if (vndx < 0)
    { if (dyvndx <= 0)
      { errmsg(369,rtnnme,orig_sys->nme,"constraint",
	       consys_nme(orig_sys,'c',-vndx,FALSE,NULL),-vndx,
	       "cons",-vndx,dyvndx) ;
	retval = dyrINV ;
	break ; } }
    else
    { if (dyvndx <= 0)
      { errmsg(369,rtnnme,orig_sys->nme,"variable",
	       consys_nme(orig_sys,'v',vndx,FALSE,NULL),vndx,
	       "vars",vndx,dyvndx) ;
	retval = dyrINV ;
	break ; }
      if (consys_getcol_pk(dy_sys,dyvndx,&pkcol) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,"variable",
	       consys_nme(orig_sys,'v',vndx,TRUE,NULL),vndx) ;
	retval = dyrFATAL ;
	break ; }
      if (pkcol->cnt == 0 && dy_opts->print.crash >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      %s (%d) has no non-zeros in active constraints.",
		    consys_nme(dy_sys,'v',dyvndx,TRUE,NULL),dyvndx) ; } }
#   endif

    dy_basis[dycndx] = dyvndx ;
    dy_var2basis[dyvndx] = dycndx ; }
/*
  If we're in fullsys mode, make the logical basic for any remaining
  constraints.
*/
  if (dy_opts->fullsys == TRUE)
  { for ( ; bpos <= dy_sys->concnt ; bpos++)
    { dy_basis[bpos] = bpos ;
      dy_var2basis[bpos] = bpos ; } }

# ifdef DYLP_PARANOIA
  pkvec_free(pkcol) ;
  if (retval != dyrOK) return (retval) ;
# endif

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\t    Pos'n Variable           Constraint") ;
    for (bpos = 1 ; bpos <= orig_lp->basis->len ; bpos++)
    { vndx = dy_basis[bpos] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t     %3d  (%3d) %-15s",bpos,vndx,
		  consys_nme(dy_sys,'v',vndx,FALSE,NULL)) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%-15s",
		  consys_nme(dy_sys,'c',bpos,FALSE,NULL)) ; } }
# endif

/*
  Factor the basis. We don't want any of the primal or dual variables
  calculated just yet. If this fails we're in deep trouble. Don't do this
  if we're dealing with a constraint system with no constraints!
*/
  if (dy_sys->concnt > 0)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.crash >= 2)
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n    factoring ...") ;
#   endif
    calcflgs = 0 ;
    retval = dy_factor(&calcflgs) ;
    switch (retval)
    { case dyrOK:
      case dyrPATCHED:
      { break ; }
      default:
      { errmsg(309,rtnnme,dy_sys->nme) ;
	return (retval) ; } } }
/*
  Attach and clear the vectors which will hold the status, values of primal and
  dual variables, and reduced costs.
*/
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(flags),(void **) &dy_status) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"status vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(double),(void **) &dy_xbasic) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"basic variable vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(double),(void **) &dy_x) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"primal variable vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(double),(void **) &dy_y) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"dual variable vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(double),(void **) &dy_cbar) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"reduced cost vector") ;
    return (dyrFATAL) ; }
/*
  Calculate dual variables and reduced costs. Might as well make a try for a
  dual feasible start, eh?
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    calculating dual values ...") ;
# endif
  dy_calcduals() ;
  if (dy_calccbar() == FALSE)
  { errmsg(384,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
    return (dyrFATAL) ; }
/*
  Initialise dy_status for logicals, using dy_var2basis and dy_cbar as guides.

  We have to consider the type of constraint so that we can give artificials
  NBFX status (thus avoiding the issue of whether NBLB or NBUB gives dual
  feasibility), and so that we can check the sign of the associated reduced
  cost to determine the proper bound for the logical associated with a range
  constraint.
*/
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    establishing initial status and reference frame ...") ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n      logicals ...") ; }
# endif
  for (dyvndx = 1 ; dyvndx <= dy_sys->concnt ; dyvndx++)
  { if (dy_var2basis[dyvndx] != 0)
    { if (vub[dyvndx] == vlb[dyvndx])
	dy_status[dyvndx] = vstatBFX ;
      else
	dy_status[dyvndx] = vstatB ; }
    else
    { switch (dy_sys->ctyp[dyvndx])
      { case contypLE:
	case contypGE:
	{ dy_status[dyvndx] = vstatNBLB ;
	  dy_x[dyvndx] = 0 ;
	  break ; }
        case contypEQ:
	{ dy_status[dyvndx] = vstatNBFX ;
	  dy_x[dyvndx] = 0 ;
	  break ; }
        case contypRNG:
	{ if (vub[dyvndx] == vlb[dyvndx])
	  { dy_status[dyvndx] = vstatNBFX ;
	    dy_x[dyvndx] = vub[dyvndx] ; }
	  else
	  if (dy_cbar[dyvndx] < 0)
	  { dy_status[dyvndx] = vstatNBUB ;
	    dy_x[dyvndx] = vub[dyvndx] ; }
	  else
	  { dy_status[dyvndx] = vstatNBLB ;
	    dy_x[dyvndx] = vlb[dyvndx] ; }
	  break ; }
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (dyrFATAL) ; } } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.crash >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  %s (%d) %s",
		  consys_nme(dy_sys,'v',dyvndx,FALSE,NULL),dyvndx,
		  dy_prtvstat(dy_status[dyvndx])) ;
      if (flgon(dy_status[dyvndx],vstatNONBASIC|vstatNBFR))
	dyio_outfmt(dy_logchn,dy_gtxecho," with value %g.",dy_x[dyvndx]) ;
      else
	dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
  }
/*
  Scan dy_origvars, with two purposes in mind:
    * For active architectural variables, initialise dy_status from
      orig_status, using the actual status for nonbasic variables, and
      vstatB, vstatBFX, or vstatBFR for basic variables. (We'll tune this
      once we have the values of the basic variables.) Initialise dy_x to the
      proper value for nonbasic variables. We shouldn't see NBFX here, as
      those variables should have been left inactive.
    * For inactive architectural variables, accumulate the objective function
      correction. Nonbasic free variables are assumed to have value 0.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n      architecturals ...") ;
# endif
  dy_lp->inactzcorr = 0 ;
  for (vndx = 1 ; vndx <= orig_sys->varcnt ; vndx++)
  { dyvndx = dy_origvars[vndx] ;
    if (dyvndx < 0)
    { obj = orig_sys->obj[vndx] ;
      switch ((flags) (-dyvndx))
      { case vstatNBFX:
        case vstatNBLB:
	{ dy_lp->inactzcorr += obj*orig_sys->vlb[vndx] ;
	  break ; }
        case vstatNBUB:
	{ dy_lp->inactzcorr += obj*orig_sys->vub[vndx] ;
	  break ; }
#       ifdef DYLP_PARANOIA
	case vstatNBFR:
	{ break ; }
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (dyrINV) ; }
#	endif
      } }
    else
    { if (((int) orig_status[vndx]) < 0)
      { if (vlb[dyvndx] == vub[dyvndx])
	  dy_status[dyvndx] = vstatBFX ;
	else
	if (vlb[dyvndx] <= -dy_tols->inf && vub[dyvndx] >= dy_tols->inf)
	  dy_status[dyvndx] = vstatBFR ;
	else
	  dy_status[dyvndx] = vstatB ; }
      else
      { dy_status[dyvndx] = orig_status[vndx] ;
	switch (dy_status[dyvndx])
	{ case vstatNBLB:
	  { dy_x[dyvndx] = vlb[dyvndx] ;
	    break ; }
	  case vstatNBUB:
	  { dy_x[dyvndx] = vub[dyvndx] ;
	    break ; }
	  case vstatNBFR:
	  { dy_x[dyvndx] = 0 ;
	    break ; }
#	  ifdef DYLP_PARANOIA
	  default:
	  { errmsg(1,rtnnme,__LINE__) ;
	    return (dyrINV) ; }
#	  endif
	} }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.crash >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  %s (%d) %s",
		    consys_nme(dy_sys,'v',dyvndx,FALSE,NULL),dyvndx,
		    dy_prtvstat(dy_status[dyvndx])) ;
	if (flgon(dy_status[dyvndx],vstatNONBASIC|vstatNBFR))
	  dyio_outfmt(dy_logchn,dy_gtxecho," with value %g.",dy_x[dyvndx]) ;
	else
	  dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#     endif
    } }
/*
  Did we patch the basis? If so, we need to scan the status array and correct
  the entries for the architectural variables that were booted out during the
  patch.
*/
  if (retval == dyrPATCHED) correct_for_patch() ;
/*
  Ok, status is set. Now it's time to calculate initial values for the primal
  variables and objective.  Arguably we don't need the true objective for
  phase I, but it's cheap to calculate.  Once we have the primal variables,
  adjust the status for any that are pinned against a bound or out of bounds,
  and see how it looks, in terms of primal infeasibility.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    calculating primal values ...") ;
# endif
  if (dy_calcprimals() == FALSE)
  { errmsg(316,rtnnme,dy_sys->nme) ;
    return (dyrFATAL) ; }
  dy_lp->z = dy_calcobj() ;
  dy_setfinalstatus() ;
/*
  Make the check for primal and/or dual feasibility, and set the initial
  simplex phase accordingly.
*/
  calcflgs = ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
  retval = dy_accchk(&calcflgs) ;
  if (retval != dyrOK)
  { errmsg(304,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
    return (retval) ; }
  if (flgoff(calcflgs,ladPRIMFEAS))
  { dy_lp->simplex.next = dyPRIMAL2 ; }
  else
  if (flgoff(calcflgs,ladDUALFEAS))
  { dy_lp->simplex.next = dyDUAL ; }
  else
  { dy_lp->simplex.next = dyPRIMAL1 ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    phase %s, initial objective %g",
	        dy_prtlpphase(dy_lp->simplex.next,FALSE),dy_lp->z) ;
    if (dy_lp->infeascnt != 0)
      dyio_outfmt(dy_logchn,dy_gtxecho,", %d infeasible vars, infeas = %g",
		  dy_lp->infeascnt,dy_lp->infeas) ;
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
  if (dy_opts->print.crash >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\nPos'n\tConstraint\tDual\t\tPrimal\n") ;
    for (bpos = 1 ; bpos <= dy_sys->concnt; bpos++)
    { cndx = dy_actcons[bpos] ;
      dyvndx = dy_basis[bpos] ;
      if (dyvndx <= dy_sys->concnt)
	vndx = orig_sys->varcnt+dyvndx ;
      else
	vndx = dy_actvars[dyvndx] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%5d\t(%4d) %-8s\t%12.4g\t(%4d) %-8s %12.4g",
		  bpos,cndx,
		  consys_nme(dy_sys,'c',bpos,FALSE,NULL),dy_y[bpos],vndx,
		  consys_nme(dy_sys,'v',dyvndx,FALSE,NULL),dy_x[dyvndx]) ; } }
# endif

  return (dyrOK) ; }

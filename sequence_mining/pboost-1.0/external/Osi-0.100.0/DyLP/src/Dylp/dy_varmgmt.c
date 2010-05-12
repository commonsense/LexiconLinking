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
  This file contains routines for primal architectural variable management.
  The two primitive routines, dy_actNBPrimArch and dy_deactNBPrimArch,
  provide activation and deactivation of a nonbasic primal architectural
  variable. The top-level routines for normal bulk activation and
  deactivation are dy_activateVars and dy_deactivateVars, respectively.
  dy_activateVars pays attention to the target simplex phase and will only
  activate variables which will be feasible for the phase.

  In the normal course of events, dylp will deactivate architecturals with
  unfavourable cbar<j> and activate architecturals with favourable cbar<j>.
  Put another way, variables already at their optimum bound and unlikely to
  ever change are deactivated. Inactive variables not at their optimum bound
  are activated for adjustment. Put another way, loose dual constraints are
  deactivated, violated dual constraints are activated.

  When attempting to recover from pivoting problems in primal simplex, dylp
  will deactivate architecturals with favourable cbar<j> in an attempt to
  force dual feasibility and allow a transition to dual simplex.

  Logical variables cannot be activated or deactivated independently; they are
  manipulated with their associated constraint. See dy_conmgmt.c

  Activating or deactivating a basic primal variable is problematic. Viewed
  from a primal perspective, activation requires we find a basis position,
  while deactivation amounts to forcing the variable to bound. Implementation
  is nontrivial, and both actions can have serious side effects. This package
  does not support activation into the basis. (Note, however, that it's
  standard procedure when activating a constraint to use the associated
  logical as the new basic variable.)

  Don't be taken in by the pseudo-primitive, dy_deactBPrimArch.  It's used
  only in the context of attempting to force primal feasibility, and it's a
  faint hope at best. The specified variable must be primal infeasible, and
  all hell will break loose when this variable is forced to bound as it's
  moved into the nonbasic partition prior to calling dy_deactNBPrimArch.

  Sometimes we want to view the column for x<j> as a dual constraint: In
  terms of the dual problem, we're activating or deactivating a dual
  constraint in order to bound the dual. The routines that control this
  activity (dy_dualaddvars and subordinates) are in dy_bound.c.
*/

/*
  A few words about the activation and deactivation of variables as it relates
  to the PSE and DSE pricing algorithms.

  Variable deactivation targets only nonbasic architectural variables.  For
  DSE pricing, that's not a problem --- the basis inverse row norms
  ||beta<i>|| are not affected.  For PSE pricing, there is again no
  difficulty. The column goes away, but the norms ||abar~<j>|| of other
  columns are not affected.

  Variable activation has no effect on DSE pricing. Again, bringing a
  variable into the nonbasic partition doesn't affect the beta<i>. For PSE
  pricing, the variable must be added to the reference frame and ||abar~<j>||
  must be calculated.

  Activation and deactivation of constraints and associated logicals is a very
  different story. See the discussion in dy_conmgmt.c.
*/


#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_varmgmt.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_varmgmt.c 269 2009-04-02 05:38:19Z lou $" ;



/*
  Reverse integer comparison so we can sort arrays of indices in nonincreasing
  order.

  Returns:	< 0	if i > j
		  0	if i = j
		> 0	if i < j
*/
static int intcompare (const void *p_i, const void *p_j)
{ int i = *((const int *) p_i) ;
  int j = *((const int *) p_j) ;
  return ((j)-(i)) ; }

static bool prepcol_pk (consys_struct *orig_sys, int oxjndx,
			pkvec_struct **p_aj)

/*
  This routine `preps' a column a<j> from the original constraint system for
  use in the active constraint system. It returns a packed vector containing
  only the active coefficients, with row indices translated from the original
  system to the active system.

  Parameters:
    orig_sys:	The original constraint system.
    oxjndx:	The index of x<j> in the original system.
    p_aj:	(i) A pkvec structure, or NULL (in which case a vector will
		    be allocated)
		(o) The prepped column a<j>.

  Returns: TRUE if a<j> is prepped without error, FALSE otherwise.
*/

{ int pkndx,ocndx,cndx ;
  pkvec_struct *aj ;
  pkcoeff_struct *aij ;

  const char *rtnnme = "prepcol_pk" ;


# ifdef DYLP_PARANOIA
/*
  We shouldn't be here if x<j> is already active.
*/
  if (p_aj == NULL)
  { errmsg(2,rtnnme,"&a<j>") ;
    return (FALSE) ; }
  if (ACTIVE_VAR(oxjndx))
  { char onmbuf[128] ;
    (void) consys_nme(orig_sys,'v',oxjndx,TRUE,onmbuf) ;
    errmsg(431,rtnnme,
	   orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable",onmbuf,oxjndx,
	   consys_nme(dy_sys,'v',dy_origvars[oxjndx],TRUE,NULL),
	   dy_origvars[oxjndx]) ;
    return (FALSE) ; }
# endif

/*
  Pull the column from orig_sys. If the user didn't supply a pkvec to hold the
  prepped column, one will be allocated in the process.
*/
  if (consys_getcol_pk(orig_sys,oxjndx,p_aj) == FALSE)
  { errmsg(122,rtnnme,orig_sys->nme,
	   "column",consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
    return (FALSE) ; }
  aj = *p_aj ;
/*
  For coefficients in active constraints, convert the row index.  Delete
  coefficients in inactive constraints.
*/
  for (pkndx = 0,aij = aj->coeffs ; pkndx < aj->cnt ; )
  { ocndx = aij->ndx ;
    if (INACTIVE_CON(ocndx))
    { cndx = dy_origcons[ocndx] ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tdeleting a<%d,%d> = %g; %s constraint %s inactive.",
		    ocndx,oxjndx,aij->val,
		    consys_prtcontyp(orig_sys->ctyp[ocndx]),
		    consys_nme(orig_sys,'c',ocndx,FALSE,NULL)) ; }
#     endif
      aj->cnt-- ;
      if (pkndx < aj->cnt)
      { aij->ndx = aj->coeffs[aj->cnt].ndx ;
	aij->val = aj->coeffs[aj->cnt].val ; } }
    else
    { cndx = dy_origcons[ocndx] ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n\ta<%d,%d> = %g becomes a<%d,%d>; %s constraint %s active.",
	       ocndx,oxjndx,aij->val,cndx,oxjndx,
	       consys_prtcontyp(orig_sys->ctyp[ocndx]),
	       consys_nme(orig_sys,'c',ocndx,FALSE,NULL)) ; }
#     endif
      aij->ndx = cndx ;
      pkndx++ ;
      aij++ ; } }

  return (TRUE) ; }



bool dy_actNBPrimArch (consys_struct *orig_sys, int ovndx)

/*
  This routine activates the nonbasic primal architectural variable x<j> with
  index ovndx in orig_sys, installing it in the active system dy_sys with
  index j and the status it held while inactive. There are a number of
  details to attend to:
    * The rhs and rhslow values of affected constraints are adjusted.
    * The objective correction is adjusted.
    * The variable is added to the PSE reference frame and the projected column
      norm gamma[j] is initialized to 1.

  Note that by convention, NBFR variables have value zero.

  It's left to the client to decide if initializing gamma[j] to 1 is adequate.
  The alternative (calculating the correct projected norm for the current
  reference frame) is best done in bulk rather than one variable at a time.

  Parameters:
    orig_sys:	The original constraint system.
    ovndx:	The index of the variable in the original system.

  Returns: TRUE if the activation completes without error; FALSE otherwise.
*/

{ int pkndx,i,j ;
  double valj,cj,lbj,ubj ;
  flags statj ;
  pkvec_struct *aj ;
  pkcoeff_struct *aij ;
  const char *rtnnme = "dy_actNBPrimArch" ;

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (ovndx <= 0 || ovndx > orig_sys->archvcnt)
  { errmsg(102,rtnnme,"inactive variable",ovndx,1,orig_sys->archvcnt) ;
    return (FALSE) ; }
  j = (orig_sys->varcnt-dy_lp->sys.vars.unloadable) -
	(dy_lp->sys.vars.loadable+dy_sys->archvcnt) ;
  if (j != 0)
  { errmsg(444,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable",orig_sys->varcnt,dy_lp->sys.vars.unloadable,
	   dy_lp->sys.vars.loadable,dy_sys->archvcnt,j) ;
    return (FALSE) ; }
  if (ACTIVE_VAR(ovndx))
  { char onmbuf[128] ;
    j = dy_origvars[ovndx] ;
    (void) consys_nme(orig_sys,'v',ovndx,TRUE,onmbuf) ;
    errmsg(431,rtnnme,
	   orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable",onmbuf,ovndx,
	   consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
    return (FALSE) ; }
# endif
/*
  Get the status. If we're paranoid, confirm that it's normal nonbasic and
  loadable.
*/
  statj = (flags) (-dy_origvars[ovndx]) ;
# ifdef DYLP_PARANOIA
  if (!LOADABLE_VAR(ovndx))
  { errmsg(445,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable",consys_nme(orig_sys,'v',ovndx,TRUE,NULL),ovndx) ;
    return (FALSE) ; }
  if (flgoff(statj,vstatNBLB|vstatNBUB|vstatNBFR))
  { errmsg(433,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "inactive",consys_nme(orig_sys,'v',ovndx,TRUE,NULL),ovndx,
	   dy_prtvstat(statj)) ;
    return (FALSE) ; }
# endif
/*
  Pull the column from orig_sys and prep it for installation in dy_sys.
*/
  aj = NULL ;
  if (prepcol_pk(orig_sys,ovndx,&aj) == FALSE)
  { errmsg(432,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   consys_nme(orig_sys,'v',ovndx,TRUE,NULL),ovndx) ;
    if (aj != NULL) pkvec_free(aj) ;
    return (FALSE) ; }
  cj = orig_sys->obj[ovndx] ;
  ubj = orig_sys->vub[ovndx] ;
  lbj = orig_sys->vlb[ovndx] ;
/*
  Install the column and update the associated structures.
*/
  if (consys_addcol_pk(dy_sys,vartypCON,aj,cj,
		       orig_sys->vlb[ovndx],orig_sys->vub[ovndx]) == FALSE)
  { errmsg(156,rtnnme,"variable",dy_sys->nme,aj->nme) ;
    pkvec_free(aj) ;
    return (FALSE) ; }
  j = aj->ndx ;
  dy_origvars[ovndx] = j ;
  dy_actvars[j] = ovndx ;
  dy_status[j] = statj ;
  dy_var2basis[j] = 0 ;
  if (dy_lp->p1obj.installed == TRUE)
  { dy_lp->p1obj.p2obj[j] = cj ;
    dy_sys->obj[j] = 0.0 ; }
/*
  Get the variable's value. If it's nonzero, scan the column and make any
  necessary adjustments to rhs and rhslow of affected constraints. Also
  adjust the objective correction.
*/
  if (flgon(statj,vstatNBLB))
  { valj = lbj ; }
  else
  if (flgon(statj,vstatNBUB))
  { valj = ubj ; }
  else
  { valj = 0.0 ; }
  dy_x[j] = valj ;

  if (valj != 0)
  { for (pkndx = 0,aij = aj->coeffs ; pkndx < aj->cnt ; pkndx++,aij++)
    { i = aij->ndx ;
      dy_sys->rhs[i] += aij->val*valj ;
      setcleanzero(dy_sys->rhs[i],dy_tols->zero) ;
      if (dy_sys->ctyp[i] == contypRNG)
      { dy_sys->rhslow[i] += aij->val*valj ;
	setcleanzero(dy_sys->rhslow[i],dy_tols->zero) ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tadjusting %s constraint %s (%d), ",
		    consys_prtcontyp(dy_sys->ctyp[i]),
		    consys_nme(dy_sys,'c',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"a<%d,%d> = %g, x<%d> = %g, ",
		    i,j,aij->val,j,valj) ;
	if (dy_sys->ctyp[i] == contypRNG)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"rhslow & ") ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"rhs += %g.", aij->val*valj) ; }
#     endif
    }
  dy_lp->inactzcorr -= cj*valj ; }
  pkvec_free(aj) ;
/*
  Add the variable to the reference frame and initialize gamma[j].
*/
  dy_frame[j] = TRUE ;
  dy_gamma[j] = 1 ;
/*
  And finally, a little bookkeeping.
*/
  dy_lp->sys.vars.loadable-- ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tadjusting objective correction, ") ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"c<%d> = %g, x<%d> = %g, zcorr -= %g.",
	        ovndx,orig_sys->obj[ovndx],ovndx,valj,
		orig_sys->obj[ovndx]*valj) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\t%s %s (%d) = %g copied to index %d, status %s.",
	        consys_prtvartyp(dy_sys->vtyp[j]),
	        consys_nme(orig_sys,'v',ovndx,FALSE,NULL),ovndx,valj,j,
	        dy_prtvstat(statj)) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->vars.actcnt[ovndx]++ ;
# endif


  return (TRUE) ; }




bool dy_actNBPrimArchList (consys_struct *orig_sys, int cnt, int *ovndxs)
/*
  This routine is purely a shell to call dy_actNBPrimArch for each of the
  indices in the vector ovndxs. It performs minimal error checking, relying
  on checking in actNBPrimArch.  One thing it does do is check if the
  variable is already active. This can happen when the list passed in avndxs
  includes duplicate indices, so we don't want it trapped later as an error.

  Parameters:
    orig_sys:	The original constraint system.
    cnt:	The number of indices in ovndxs
    ovndxs:	Vector of variable indices (0-based)

  Returns: TRUE if all variables are successfully activated, FALSE otherwise
*/

{ int j,k ;
  bool retval ;
  const char *rtnnme = "dy_actNBPrimArchList" ;

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (ovndxs == NULL)
  { errmsg(2,rtnnme,"ovndxs") ;
    return (FALSE) ; }
  if (cnt <= 0 || cnt > orig_sys->archvcnt)
  { errmsg(5,rtnnme,"cnt",cnt) ;
    return (FALSE) ; }
# endif

  retval = TRUE ;
  for (k = 0 ; k < cnt && retval == TRUE ; k++)
  { j = ovndxs[k] ;
    if (ACTIVE_VAR(j)) continue ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    activating variable %s (%d)",
		  consys_nme(orig_sys,'v',j,TRUE,NULL),j) ; }
#   endif
    retval = dy_actNBPrimArch(orig_sys,j) ;
    if (retval == FALSE)
    { errmsg(430,rtnnme,
	     orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "activate","variable",
	     consys_nme(orig_sys,'v',j,TRUE,NULL),j) ; } }

# ifdef DYLP_PARANOIA
  if (retval == TRUE)
  { retval = dy_chkdysys(orig_sys) ; }
# endif

  return (retval) ; }



bool dy_deactBPrimArch (consys_struct *orig_sys, int j)

/*
  This routine deactivates the basic primal architectural variable x<j>,
  removing it from dy_sys. The only reason we want to do this is because
  we're trying to force primal feasibility. The chance of achieving primal
  feasibility is pretty well zero, given that we can't really deactivate the
  implicit bound constraint --- as we force x<j> into the nonbasic partition,
  we simultaneously force it to the nearest bound. But in the context of dual
  error recovery, it will serve the purpose of moving the basis.

  We can't deactivate a basic variable, so we need to force x<j> into the
  nonbasic partition. Suppose x<j> occupies basis pos'n i. The obvious
  candidate to fill the slot is the logical for constraint i. If it's
  nonbasic, we swap it in with value 0 and then call dy_deactNBPrimArch to
  deactivate the newly nonbasic x<j>. But it could also be that x<i> is basic
  in some other position, k. In which case we just try again, with the logical
  x<k> as the new target.

  Why, you ask, are we hunting down logicals, instead of just swapping with
  any old variable? Because we'll have enough problems after this without
  accidentally creating a singular basis.

  Parameters:
    orig_sys:	Original constraint system
    j:		The variable to be deactivated

  Returns: TRUE if the deactivation succeeds, FALSE otherwise.
*/

{ int i,k ;
  double valj ;
  flags statj,stati ;

  const char *rtnnme = "dy_deactBPrimArch" ;

/*
  Prep intermixed with paranoia. We need to retrieve the status and check that
  it's what we expect.
*/
# ifdef DYLP_PARANOIA
  if (j <= dy_sys->logvcnt || j > dy_sys->varcnt)
  { errmsg(102,rtnnme,"active variable",j,dy_sys->logvcnt+1,dy_sys->varcnt) ;
    return (FALSE) ; }
# endif

  statj = dy_status[j] ;

# ifdef DYLP_PARANOIA
  if (flgoff(statj,vstatBLLB|vstatBUUB))
  { errmsg(438,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   consys_nme(dy_sys,'v',j,TRUE,NULL),j,dy_prtvstat(statj)) ;
    return (FALSE) ; }
# endif

/*
  Search for a nonbasic logical x<i> which we can swap into the basis to
  replace x<j>. No loop body required, except for paranoia.
*/
  for (i = dy_var2basis[j] ; dy_var2basis[i] != 0 ; i = dy_var2basis[i])
  {
#   ifdef DYLP_PARANOIA
    if (i <= 0 || i > dy_sys->concnt)
    { errmsg(102,rtnnme,"logical variable",i,1,dy_sys->concnt) ;
      return (FALSE) ; }
#   endif
  }
/*
  We have to determine the appropriate new status for both x<i> and x<j>.
  Arguably, we shouldn't see vstatSB here (superbasic and dual feasible just
  don't go together), but it's trivial to handle. Capture x<j>'s new value
  while we're at it.
*/
  statj = getflg(dy_status[j],vstatSTATUS) ;
  switch (statj)
  { case vstatBLLB:
    { statj = vstatNBLB ;
      valj = dy_sys->vlb[j] ;
      break ; }
    case vstatBUUB:
    { statj = vstatNBUB ;
      valj = dy_sys->vub[j] ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
  stati = getflg(dy_status[i],vstatSTATUS) ;
  switch (stati)
  { case vstatNBLB:
    { stati = vstatBLB ;
      break ; }
    case vstatNBUB:
    { stati = vstatBUB ;
      break ; }
    case vstatNBFX:
    { stati = vstatBFX ;
      break ; }
    case vstatNBFR:
    { stati = vstatBFR ;
      break ; }
    case vstatSB:
    { stati = vstatB ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      swapping %s (%d) %s -> ",
	        consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	        dy_prtvstat(dy_status[j])) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"%s ",dy_prtvstat(statj)) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"<=> %s (%d) %s -> ",
	        consys_nme(dy_sys,'v',i,FALSE,NULL),i,
	        dy_prtvstat(dy_status[i])) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"%s.",dy_prtvstat(stati)) ; }
# endif
/*
  Tweak dy_basis, dy_var2basis, dy_status, dy_x, and dy_xbasic to do the
  swap.  The entries for x<j> are going to disappear momentarily, but
  deactNBPrimArch will want to examine them to do the deactivation.
  Guaranteed we change the composition of the PSE and DSE norms with this
  maneuver.
*/
  k = dy_var2basis[j] ;
  dy_basis[k] = i ;
  dy_xbasic[k] = dy_x[i] ;
  dy_var2basis[i] = k ;
  dy_status[i] = stati ;
  dy_var2basis[j] = 0 ;
  dy_status[j] = statj ;
  dy_x[j] = valj ;
  dy_lp->simplex.init_pse = TRUE ;
  dy_lp->simplex.init_dse = TRUE ;
/*
  Reality is altered. Call dy_deactNBPrimArch to do the heavy lifting.
*/
  return (dy_deactNBPrimArch(orig_sys,j)) ; }



bool dy_deactNBPrimArch (consys_struct *orig_sys, int j)

/*
  This routine deactivates the nonbasic primal architectural variable x<j>,
  removing it from dy_sys. There are a number of details to attend to:
    * The rhs and rhslow values for affected constraints are adjusted.
    * The objective correction is adjusted.
  If consys shifts a variable to fill the hole at index j (the common case),
  we also need to
    * Update the dy_origvars entry for the new x<j>.
    * Check if the new x<j> is basic, and, if so, correct dy_basis.

  Note that by convention, NBFR variables have value zero.

  We don't have to remove x<j> from the reference frame --- since dy_frame and
  dy_gamma are attached to dy_sys, removal occurs automatically as a side
  effect of removing the variable from dy_sys.

  Parameters:
    orig_sys:	The original constraint system.
    j:		The index of the variable in the active system.

  orig_sys is used only for printing and paranoid checks, but it gives a nice
  symmetry with dy_actNBPrimArch.

  Returns: TRUE if the installation completes without error; FALSE otherwise.
*/

{ int ovndx,i,pkndx ;
  double valj ;
  flags statj ;
  pkvec_struct *aj ;
  pkcoeff_struct *aij ;
  const char *rtnnme = "dy_deactNBPrimArch" ;

/*
  Prep intermixed with paranoia. We need to retrieve the status and check that
  it's appropriate for deactivation, and that we have a valid index in the old
  system.
*/
# ifdef DYLP_PARANOIA
  if (j <= dy_sys->logvcnt || j > dy_sys->varcnt)
  { errmsg(102,rtnnme,"active variable",j,dy_sys->logvcnt+1,dy_sys->varcnt) ;
    return (FALSE) ; }
# endif

  statj = dy_status[j] ;

# ifdef DYLP_PARANOIA
  if (flgoff(statj,vstatNONBASIC|vstatNBFR))
  { errmsg(433,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "active",consys_nme(dy_sys,'v',j,TRUE,NULL),j,
	   dy_prtvstat(statj)) ;
    return (FALSE) ; }
# endif

  ovndx = dy_actvars[j] ;

# ifdef DYLP_PARANOIA
  if (ovndx <= 0 || ovndx > orig_sys->varcnt)
  { errmsg(102,rtnnme,"original variable",ovndx,1,orig_sys->varcnt) ;
    return (FALSE) ; }
  i = (orig_sys->varcnt-dy_lp->sys.vars.unloadable) -
	(dy_lp->sys.vars.loadable+dy_sys->archvcnt) ;
  if (i != 0)
  { errmsg(444,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable", orig_sys->varcnt,dy_lp->sys.vars.unloadable,
	   dy_lp->sys.vars.loadable,dy_sys->archvcnt,i) ;
    return (FALSE) ; }
# endif

/*
  More prep. Retrieve the column and the value of the variable.
*/
  aj = NULL ;
  if (consys_getcol_pk(dy_sys,j,&aj) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,"variable",
	   consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
    if (aj != NULL) pkvec_free(aj) ;
    return (FALSE) ; }
  valj = dy_x[j] ;
/*
  Ok, things look good --- we can locate the variable in the active and
  original systems, and we have the column in hand. To deactivate, we
  clean the status and set it in dy_origvars to mark the variable inactive.
  If the outgoing variable is NBFX, mark it with the NOLOAD qualifier.
  If x<j> is nonzero, traverse the column and correct the rhs and rhslow
  values, and adjust the objective correction.
*/
  clrflg(statj,vstatQUALS) ;
  if (statj == vstatNBFX) setflg(statj,vstatNOLOAD) ;
  MARK_INACTIVE_VAR(ovndx,-((int) statj)) ;
  if (valj != 0)
  { for (pkndx = 0, aij = aj->coeffs ; pkndx < aj->cnt ; pkndx++, aij++)
    { i = aij->ndx ;
      dy_sys->rhs[i] -= aij->val*valj ;
      setcleanzero(dy_sys->rhs[i],dy_tols->zero) ;
      if (dy_sys->ctyp[i] == contypRNG)
      { dy_sys->rhslow[i] -= aij->val*valj ;
	setcleanzero(dy_sys->rhslow[i],dy_tols->zero) ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tadjusting %s constraint %s (%d), ",
		    consys_prtcontyp(dy_sys->ctyp[i]),
		    consys_nme(dy_sys,'c',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"a<%d,%d> = %g, x<%d> = %g, ",
		    i,j,aij->val,j,valj) ;
	if (dy_sys->ctyp[i] == contypRNG)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"rhslow & ") ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"rhs -= %g.", aij->val*valj) ; }
#     endif
    }
    dy_lp->inactzcorr += dy_sys->obj[j]*valj ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tadjusting objective correction, ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"c<%d> = %g, x<%d> = %g, zcorr += %g.",
		  j,dy_sys->obj[j],j,valj,dy_sys->obj[j]*valj) ; }
#   endif
  }
  pkvec_free(aj) ;
/*
  Delete the column from the active system. Row vectors attached to dy_sys
  will be automatically rearranged by the consys package.
*/
  if (consys_delcol(dy_sys,j) == FALSE)
  { errmsg(112,rtnnme,dy_sys->nme,"delete","variable",
	   consys_nme(dy_sys,'v',j,FALSE,NULL),j) ;
    return (FALSE) ; }

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->vars.deactcnt[ovndx]++ ;
# endif

/*
  Unless x<j> was the last variable in dy_sys, the last variable has been
  moved into slot j to close up the hole. Normally, there will be two
  remaining details to attend to:
    * Consult dy_actvars[j] to find the index ovndx of x<j> in the original
      system, and correct dy_origvars[ovndx].
    * Check dy_var2basis[j], and see if the variable moved into this slot
      is basic. If so, correct dy_basis.
*/
  if (j <= dy_sys->varcnt)
  { ovndx = dy_actvars[j] ;
#   ifdef DYLP_PARANOIA
    if (dy_origvars[ovndx] != dy_sys->varcnt+1)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; }
#   endif
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) shifted from column %d",
		  consys_nme(dy_sys,'v',j,FALSE,NULL),ovndx,
		  dy_origvars[ovndx]) ; }
#   endif
    dy_origvars[ovndx] = j ;

    i = dy_var2basis[j] ;
    if (i != 0)
    { 
#     ifdef DYLP_PARANOIA
      if (dy_basis[i] != dy_sys->varcnt+1)
      { errmsg(1,rtnnme,__LINE__) ;
        return (FALSE) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,", basis entry %d corrected",i) ; }
#     endif
      dy_basis[i] = j ; } }
/*
  And finally ... if the status was anything other than NBFX, we now have a
  variable we can activate. But if it was NBFX, well, the number of unloadable
  variables just increased.
*/
  if (flgon(statj,vstatNBFX))
  { dy_lp->sys.vars.unloadable++ ; }
  else
  { dy_lp->sys.vars.loadable++ ; }

  return (TRUE) ; }



static bool deactNBPrimArchList (consys_struct *orig_sys,
				   int cnt, int *avndxs)
/*
  This routine is purely a shell to call dy_deactNBPrimArch for each of the
  indices in the vector avndxs. It performs minimal error checking, relying
  on checking in dy_deactNBPrimArch.

  Parameters:
    orig_sys:	The original constraint system
    cnt:	The number of indices in avndxs
    avndxs:	Vector of variable indices (0-based)

  orig_sys is used only for printing and paranoid checks, but it gives a nice
  symmetry with actNBPrimArchList.

  Returns: TRUE if all variables are successfully deactivated, FALSE otherwise
*/

{ int k ;
  bool retval ;
  const char *rtnnme = "deactNBPrimArchList" ;

# ifdef DYLP_PARANOIA
  if (avndxs == NULL)
  { errmsg(2,rtnnme,"avndxs") ;
    return (FALSE) ; }
  if (cnt <= 0 || cnt > dy_sys->archvcnt)
  { errmsg(5,rtnnme,"cnt",cnt) ;
    return (FALSE) ; }
# endif

/*
  To make sure consys doesn't shift variables out from under us, we need to
  delete in nonincreasing order.
*/
  qsort(&avndxs[0],cnt,sizeof(int),intcompare) ;
  
  retval = TRUE ;
  for (k = 0 ; k < cnt && retval == TRUE ; k++)
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    deactivating variable %s (%d)",
		  consys_nme(dy_sys,'v',avndxs[k],TRUE,NULL),avndxs[k]) ; }
#   endif
    retval = dy_deactNBPrimArch(orig_sys,avndxs[k]) ;
    if (retval == FALSE)
    { errmsg(430,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "deactivate","variable",
	     consys_nme(dy_sys,'v',avndxs[k],TRUE,NULL),avndxs[k]) ; } }

# ifdef DYLP_PARANOIA
  if (retval == TRUE)
  { retval = dy_chkdysys(orig_sys) ; }
# endif

  return (retval) ; }



static int scanPrimVarStdDeact (int **p_avndxs)
/*
  This routine scans the active constraint system for variables that are
  candidates for normal variable deactivation during primal simplex, i.e.,
  nonbasic architectural variables which are judged unlikely to be useful in
  improving the solution because they have a lousy reduced cost.  The routine
  begins by scanning dy_cbar to determine a deactivation threshold, then
  purges all variables with reduced costs which exceed dy_tols.purgevar times
  the threshold.

  Variables with status NBFX are deactivated if we happen to scan them, but
  they're relatively rare and we won't do a deactivation scan just for them.

  Parameters:
    p_avndxs:	(i) empty vector to hold variable indices; assumed to be
		    sufficiently large; will be allocated if NULL
		(o) indices of variables to be deactivated; may not be
		    allocated if no candidates are identified

  Returns: number of variables to be deactivated, -1 if there's an error during
	   scanning (error is possible only when paranoid)
*/

{ int j,m,n,purgecnt,*avndxs ;
  double maxpcbar,pthresh,maxncbar,nthresh,cbarj ;
  flags statj ;
  bool purge ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "scanPrimVarStdDeact" ;

  if (p_avndxs == NULL)
  { errmsg(2,rtnnme,"avndxs") ;
    return (-1) ; }
# endif
/*
  Scan the nonbasic variables to calculate the deactivation thresholds. We're
  interested in variables with status NBLB and cbar<j> > 0, and variables with
  status NBUB and cbar<j> < 0.
*/
  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
  maxpcbar = 0 ;
  maxncbar = 0 ;
  for (j = m+1 ; j <= n ; j++)
  { if (flgon(dy_status[j],vstatNBLB))
    { cbarj = dy_cbar[j] ;
      if (cbarj > maxpcbar) maxpcbar = cbarj ; }
    else
    if (flgon(dy_status[j],vstatNBUB))
    { cbarj = dy_cbar[j] ;
      if (cbarj < maxncbar) maxncbar = cbarj ; } }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    %g <= cbar<k> <= %g over deactivation candidates.",
		maxncbar,maxpcbar) ; }
# endif

  purge = FALSE ;
  if (maxpcbar > 0)
  { pthresh = dy_tols->purgevar*maxpcbar ;
    if (pthresh > dy_tols->bogus*dy_tols->dfeas)
    { purge = TRUE ; }
    else
    { pthresh = dy_tols->inf ; } }
  else
  { pthresh = dy_tols->inf ; }
  if (maxncbar < 0)
  { nthresh = dy_tols->purgevar*maxncbar ;
    if (-nthresh > dy_tols->bogus*dy_tols->dfeas)
    { purge = TRUE ; }
    else
    { nthresh = -dy_tols->inf ; } }
  else
  { nthresh = -dy_tols->inf ; }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 2)
  { if (purge == FALSE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n    %g <= purge threshold <= %g; tol = %g; scan aborted.",
	     dy_tols->purgevar*maxncbar,dy_tols->purgevar*maxpcbar,
	     dy_tols->bogus*dy_tols->dfeas) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    %g <= purge threshold <= %g; tol = %g; scanning.",
		  nthresh,pthresh,dy_tols->bogus*dy_tols->dfeas) ; } }
# endif
  if (purge == FALSE) return (0) ;
/*
  Scan preparation.  If the user hasn't supplied a vector to hold the
  indices, allocate one now.
*/
  purgecnt = 0 ;
  if (*p_avndxs == NULL)
  { avndxs = (int *) MALLOC(dy_sys->archvcnt*sizeof(int)) ; }
  else
  { avndxs = *p_avndxs ; }
/*
  Scan the architecturals again, this time with an eye to collecting a
  list of variables to be deactivated.
*/
  for (j = m+1 ; j <= n ; j++)
  { purge = FALSE ;
    statj = dy_status[j] ;
    cbarj = dy_cbar[j] ;
/*
  Do the tests for purging: NBFX status, or NBLB or NBUB status and reduced
  cost over the purge threshold.
*/
    if ((flgon(statj,vstatNBLB) && cbarj > pthresh) ||
        (flgon(statj,vstatNBUB) && cbarj < nthresh) || flgon(statj,vstatNBFX))
    { purge = TRUE ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n    queuing %s (%d) for deactivation, %s, cbar<%d> = %g",
	       consys_nme(dy_sys,'v',j,TRUE,NULL),j,
	       dy_prtvstat(statj),j,cbarj) ; }
#     endif
    }
    if (purge == TRUE)
    { avndxs[purgecnt++] = j ; } }
/*
  If we supplied avndxs and found no variables to deactivate, free the vector.
  If we found candidates, set avndxs into p_avndxs.
*/
  if (*p_avndxs == NULL)
  { if (purgecnt == 0)
    { FREE(avndxs) ; }
    else
    { (*p_avndxs) = avndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d variables for deactivation.",purgecnt) ; }
# endif

  return (purgecnt) ; }



static int scanPrimVarStdAct (consys_struct *orig_sys,
			      int **p_ovndxs, int *preset)
/*
  This routine scans the reduced cost of a set of inactive variables,
  selecting suitable variables for activation.

  The set of inactive variables can be specified using the preset parameter.
  If preset is missing, the set of variables to consider is all inactive
  variables.
  
  `Suitable' depends on the simplex phase. In primal simplex, a variable must
  price out with a favourable (nonoptimal) reduced cost. In dual simplex, a
  variable must price out as dual feasible (equivalently, an unfavourable
  (optimal) reduced cost).

  Parameters:
    orig_sys:	The original constraint system.
    p_ovndxs:	(i) empty vector to hold variable indices; assumed to be
		    sufficiently large; will be allocated if NULL
		(o) indices of variables to be activated; may not be
		    allocated if no candidates are identified
    preset	A preset list of candidates to consider. preset[0] should
		contain the number of candidates in the list.

  Returns: number of variables to activate, -1 if there's an error during
  	   scanning.
*/

{ int i,j,k,actcnt,cand_limit ;
  int *ovndxs ;
  int *scanvars,scan_cnt,ndx ;
  double *orig_obj,cbarj,*orig_y ;
  bool fatal,use_all,activate ;
  flags statj ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "scanPrimVarStdAct" ;

  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (-1) ; }
  if (orig_sys->obj == NULL)
  { errmsg(101,rtnnme,orig_sys->nme,consys_assocnme(NULL,CONSYS_OBJ)) ;
    return (-1) ; }
# endif
# ifndef DYLP_NDEBUG
  cbarj = -1 ;		/* Suppress a compiler warning. */
# endif

  fatal = FALSE ;
  orig_obj = orig_sys->obj ;
  actcnt = 0 ;
/*
  Did the client supply a vector for candidates? If not, make one.

  We shouldn't be here if there's no room to activate. If we're paranoid, check
  this. 
*/
  cand_limit = dy_lp->sys.vars.loadable ;
# ifdef DYLP_PARANOIA
  if (cand_limit == 0)
  { errmsg(1,rtnnme,__LINE__) ;
    return (-1) ; }
# endif
  if (dy_opts->addvar > 0)
  { cand_limit = minn(dy_opts->addvar,cand_limit) ; }
  if (*p_ovndxs == NULL)
  { ovndxs = (int *) MALLOC(cand_limit*sizeof(int)) ; }
  else
  { ovndxs = *p_ovndxs ; }
/*
  Make a vector of duals that matches orig_sys, for efficient evaluation of
  columns.
*/
  orig_y = (double *) CALLOC((orig_sys->concnt+1),sizeof(double)) ;
  for (i = 1 ; i <= dy_sys->concnt ; i++)
  { k = dy_actcons[i] ;
    orig_y[k] = dy_y[i] ; }
/*
  Make the candidate list we'll actually evaluate. If not supplied in preset,
  scan origvars to make the list. If we're scanning here, avoid unloadable
  variables, but we can't guarantee the contents of preset and will need to
  screen again in the next loop.

  Even if we're supplied a preset list of variables, we can use only those
  that are dual feasible if we're in dual simplex.

  Note that scanvars is set up with 1-based indexing.
*/
  if (preset != NULL)
  { scanvars = preset ;
    scan_cnt = preset[0] ;
    if (dy_lp->simplex.next == dyDUAL)
      use_all = FALSE ;
    else
      use_all = TRUE ; }
  else
  { use_all = FALSE ;
    scan_cnt = orig_sys->archvcnt-dy_sys->archvcnt ;
    scanvars = (int *) MALLOC((scan_cnt+1)*sizeof(int)) ;
    scanvars[0] = 0 ;
    scan_cnt = 0 ;
    for (j = 1 ; j <= orig_sys->archvcnt ; j++)
    { if (!LOADABLE_VAR(j)) continue ;
#     ifdef DYLP_PARANOIA
      statj = (flags) -dy_origvars[j] ;
      if (flgoff(statj,vstatNONBASIC|vstatNBFR))
      { errmsg(433,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "active",consys_nme(orig_sys,'v',j,TRUE,NULL),j,
	       dy_prtvstat(statj)) ;
	fatal = TRUE ;
	break ; }
#     endif
      scanvars[++scan_cnt] = j ; }
    if (fatal == TRUE)
    { if (scanvars != NULL && scanvars != preset) FREE(scanvars) ;
      if (orig_y != NULL) FREE(orig_y) ;
      return (-1) ; } }
/*
  Open a loop to walk scanvars checking for architectural variables that
  price out favourably.
*/
  for (ndx = 1 ; ndx <= scan_cnt && actcnt < cand_limit ; ndx++)
  { j = scanvars[ndx] ;
/*
  Skip over variables that are ineligible for activation (status NBFX or
  flagged with vstatNOLOAD).
*/
    if (!LOADABLE_VAR(j))
    { activate = FALSE ; }
/*
  If the variable x<j> is inactive, price it as cbar<j> = c<j> - y<i>a<i,j>,
  taking only the active rows into account (by construction of orig_y). If
  we're in primal phase I, c<j> is identically 0 by definition, since
  inactive variables must be at their upper or lower bound.

  If our target is dual simplex, we want variables with optimal reduced costs,
  the negative of the primal case.
*/
    else
    { statj = (flags) -dy_origvars[j] ;
      if (dy_lp->simplex.next == dyPRIMAL1)
	cbarj = 0 ;
      else
	cbarj = orig_obj[j] ;
      cbarj -= consys_dotcol(orig_sys,j,orig_y) ;
      setcleanzero(cbarj,dy_tols->cost) ;
      if (dy_lp->simplex.next == dyDUAL)
      { cbarj = -cbarj ; }
/*
  Now compare the sign of the reduced cost with the status of the variable,
  to decide if we want to bring this variable into the active set. We're
  minimising, eh.

  Arguably, if we want to detect columns with nonzero coefficients only in
  rows that are ineligible for activation, this is the place. But it's not
  clear that it's worth the effort. This particular pathology is rare and we
  would expend a lot of effort prepping columns, even with a filter of
  cbar<j> = 0. The cost for doing nothing is that we redo consys_dotcol for
  the pathological columns with each scan. Given that this sort of column is
  typically rare and sparse, doing nothing special makes sense.
*/
      if (use_all == TRUE)
      { activate = TRUE ; }
      else
      if (cbarj == 0 ||
	  (cbarj > 0 && statj == vstatNBLB) ||
	  (cbarj < 0 && statj == vstatNBUB))
      { activate = FALSE ; }
      else
      { activate = TRUE ; } }

#   ifndef DYLP_NDEBUG
    if (activate == FALSE)
    { if (dy_opts->print.varmgmt >= 3)
      { statj = (flags) -dy_origvars[j] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n    skipping %s %s (%d), status %s",
		    consys_prtvartyp(orig_sys->vtyp[j]),
		    consys_nme(orig_sys,'v',j,FALSE,NULL),j,
		    dy_prtvstat(statj)) ;
	if (flgon(statj,vstatNBFX))
	  dyio_outchr(dy_logchn,dy_gtxecho,'.') ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,", cbar = %g.",cbarj) ; } }
    else
    { if (dy_opts->print.varmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    activating %s %s (%d), status %s, cbar = %g.",
		    consys_prtvartyp(orig_sys->vtyp[j]),
		    consys_nme(orig_sys,'v',j,FALSE,NULL),j,
		    dy_prtvstat(statj),cbarj) ; } }
#   endif
    
    if (activate == TRUE) ovndxs[actcnt++] = j ; }
  if (orig_y != NULL) FREE(orig_y) ;
  if (scanvars != NULL && scanvars != preset) FREE(scanvars) ;
/*
  If we supplied ovndxs and found no candidates to activate, free the vector.
*/
  if (*p_ovndxs == NULL)
  { if (actcnt == 0)
    { FREE(ovndxs) ; }
    else
    { *p_ovndxs = ovndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d variables for activation.",actcnt) ; }
# endif

  return (actcnt) ; }



int dy_deactivateVars (consys_struct *orig_sys)
/*
  This routine is a simple coordination shell for normal variable deactivation
  from primal simplex in phase dyPURGEVAR.

  Parameters:
    orig_sys:	The original constraint system

  orig_sys is only required for debug printing and paranoid checks. Still, it's
  really convenient to have it as a parameter and it provides a nice symmetry
  with activateVars.

  Returns: number of variables deactivated, or -1 if there's an error.
*/

{ int *candidates,candcnt ;
  int retval ;
  const char *rtnnme = "dy_deactivateVars" ;

  retval = -1 ;
/*
  Call scanPrimVarStdDeact to return a list of candidates for deactivation.
  If we find candidates, try to deactivate them.
*/
  candidates = NULL ;
  candcnt = scanPrimVarStdDeact(&candidates) ;
  if (candcnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable","normal deactivation") ; }
  else
  if (candcnt > 0)
  { if (deactNBPrimArchList(orig_sys,candcnt,candidates) == TRUE)
    { retval = candcnt ; } }
  else
  { retval = 0 ; }
/*
  Clean up and return.
*/
  if (candidates != NULL) FREE(candidates) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 1)
  { if (dy_opts->print.varmgmt >= 2)
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ;
    dyio_outfmt(dy_logchn,dy_gtxecho," %d deactivations.",candcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  constraint system %s now %d x %d (%d + %d).",
	        dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
	        dy_sys->logvcnt) ; }
# endif

  return (retval) ; }



int dy_activateVars (consys_struct *orig_sys, int *preset)

/*
  This routine handles normal variable activation into primal simplex in
  phase dyADDVAR. Most of the heavy work is performed by scanPrimVarStdAct
  and actNBPrimArchList.  A little bit of cleanup is then required to
  update PSE norms.

  If a list of candidates is supplied in preset, only those variables are
  considered.

  The principle used in choosing PSE initialisation or update is a level
  playing field.  If we're entering primal simplex from dual simplex, the PSE
  norms have not been maintained, and it's reasonable to initialise all the
  norms to 1. But if we're just taking a break from primal simplex pivoting
  to add variables, it makes sense to calculate the correct norms for the
  variables we're adding to keep them on equal footing with variables already
  active.

  Parameters:
    orig_sys:	The original constraint system.
    preset:	A preset list of candidate variables to consider.

  Returns: number of variables activated, or -1 if there's an error.
*/

{ int candcnt,candndx,j,i,bfcnt,bfndx ;
  int *candidates,*basic_frame ;
  double abarij,gammaj ;
  double *abarj ;
  flags calcflgs ;
  bool actresult,pseresult ;
  dyret_enum factorresult ;
  int retval ;
  const char *rtnnme = "dy_addvars" ;

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (dyINV) ; }
# endif

  retval = -1 ;
  actresult = TRUE ;
  pseresult = TRUE ;
/*
  Make sure we have the correct objective and reduced costs for the target
  simplex.
*/
  if (dy_lp->simplex.next == dyPRIMAL1 && dy_lp->p1obj.installed == FALSE)
  { if (dy_initp1obj() == FALSE)
    { errmsg(318,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "initialise") ;
      return (-1) ; } }
  else
  if (dy_lp->simplex.next == dyPRIMAL2 && dy_lp->p1obj.installed == TRUE)
  { if (dy_swapobjs(dyPRIMAL2) == FALSE)
    { errmsg(318,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "remove") ;
      return (-1) ; }
    dy_calcduals() ;
    if (dy_calccbar() == FALSE)
    { errmsg(384,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters) ;
      return (-1) ; } }
/*
  Call scanPrimVarStdAct to get a list of candidates. If we get candidates
  back, install them. Installing nothing always succeeds.
*/
  candidates = NULL ;
  candcnt = scanPrimVarStdAct(orig_sys,&candidates,preset) ;
  if (candcnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable","normal activation") ;
    actresult = FALSE ; }
  else
  if (candcnt > 0)
  { actresult = dy_actNBPrimArchList(orig_sys,candcnt,candidates) ; }
  else
  { actresult = TRUE ; }

  if (actresult == FALSE)
  { if (candidates != NULL) FREE(candidates) ;
    return (retval) ; }
/*
  If we're just taking a break from primal simplex pivoting to add variables,
  calculate correct projected column norms gamma<j> = ||abar~<j>||^2 + 1 for
  the new variables. (Note that dy_actNBPrimArch has already added x<j> to
  the reference frame and 1nitialised gamma<j> to 1.) To make the update a
  bit more efficient, scan out a vector of reference frame members who are
  currently basic --- these are the variables of interest in calculating the
  projected norms.

  The safest way to make this decision is to look at dy_lp->simplex.init_pse.
  If it's false, we need to to this update.
*/
  if (candcnt > 0 && dy_lp->simplex.init_pse == FALSE)
  { basic_frame = (int *) MALLOC(dy_sys->concnt*sizeof(int)) ;
    bfcnt = 0 ;
    for (i = 1 ; i <= dy_sys->concnt ; i++)
    { j = dy_basis[i] ;
      if (dy_frame[j] == TRUE)
      { basic_frame[bfcnt++] = i ; } }

    abarj = NULL ;
    for (candndx = 0 ; candndx < candcnt ; candndx++)
    { j = dy_origvars[candidates[candndx]] ;
      if (consys_getcol_ex(dy_sys,j,&abarj) == FALSE)
      { errmsg(122,rtnnme,dy_sys->nme,"column",
	       consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
	pseresult = FALSE ;
	break ; }
      dy_ftran(abarj,FALSE) ;
      gammaj = 1 ;
      for (bfndx = 0 ; bfndx < bfcnt ; bfndx++)
      { i = basic_frame[bfndx] ;
	abarij = abarj[i] ;
	if (abarij != 0)
	{ gammaj += abarij*abarij ; } }
      dy_gamma[j] = gammaj ; }
    if (abarj != NULL) FREE(abarj) ;
    if (basic_frame != NULL) FREE(basic_frame) ; }

  if (candidates != NULL) FREE(candidates) ;
  if (pseresult == FALSE) return (retval) ;
/*
  One final point. If we've added variables, we need to get the reduced costs
  correct. The easy way is to call for a dual feasibility check, which we
  should pass if we're headed for dual simplex, and fail if we're headed for
  primal simplex.
*/
  retval = candcnt ;

  if (candcnt > 0)
  { calcflgs = ladDUALFEAS|ladDFQUIET ;
    factorresult = dy_accchk(&calcflgs) ;
#   if defined(DYLP_PARANOIA) || !defined(DYLP_NDEBUG)
    switch (factorresult)
    { case dyrOK:
      { 
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	{ if (factorresult == dyrOK)
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ; }
#       endif
#	ifdef DYLP_PARANOIA
/*
  In spite of the comment above, there's always numerical accuracy, and we
  can run into a situation where we're adding a variable with tol.cost <
  fabs(cbarj) < tol.dfeas. I'd rather let the main simplex deal with this, so
  all we do here is issue a warning.
*/
	if (dy_lp->simplex.next == dyDUAL)
	{ if (flgon(calcflgs,ladDUALFEAS))
	  { warn(439,rtnnme,
		 dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		 "loss","dual") ; } }
	else
	{ if (flgoff(calcflgs,ladDUALFEAS))
	  { warn(439,rtnnme,
		 dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		 "gain","dual") ; } }
#	endif
	break ; }
      default:
      { retval = -1 ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#       endif
	break ; } }
#   endif
  }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 1)
  { if (dy_opts->print.varmgmt >= 2)
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ;
    dyio_outfmt(dy_logchn,dy_gtxecho," %d activations.",candcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    constraint system %s now %d x %d (%d + %d).",
	        dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
	        dy_sys->logvcnt) ; }
# endif

  return (retval) ; }


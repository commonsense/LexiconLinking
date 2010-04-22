/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2005 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

#define DYLP_INTERNAL

#include "dylp.h"
#include <limits.h>

static char sccsid[] UNUSED = "@(#)dy_conmgmt.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_conmgmt.c 269 2009-04-02 05:38:19Z lou $" ;

/*
  This file contains routines for primal constraint management. It provides
  routines to handle the the activation and deactivation of a primal
  constraint with a basic logical variable, and deactivation of a primal
  constraint with a nonbasic logical variable.  There are also routines to
  scan for candidates for activation and deactivation. The normal top-level
  routines for bulk activation and deactivation are dy_activateCons and
  dy_deactivateCons, respectively.

  In terms of the dual problem, we're activating or deactivating a dual
  architectural variable.
  
  In the normal course of events, dylp will deactivate loose primal
  constraints (basic feasible logicals; nonbasic dual architecturals with
  unfavourable reduced costs) and activate violated constraints (basic
  infeasble logicals; nonbasic dual architecturals with favourable reduced
  costs). (Recall the the reduced cost of a dual is the value of the primal
  variable in that basis position.)

  When attempting to recover from pivoting problems in the dual simplex, dylp
  will deactivate dual architecturals on the pivot reject list (corresponding
  to violated constraints; these will have favourable dual reduced cost) in
  an attempt to achieve primal feasibility and allow a transition to primal
  phase II.

  Activating or deactivating a constraint with a nonbasic logical is
  problematic. Viewed from a primal perspective, activation requires we find
  a variable to occupy the basis position, and deactivation requires we deal
  with the variable that's basic for the constraint.

  When attempting to regain dual feasibility, dylp may request deactivation
  of a constraint with a nonbasic logical because the logical is dual
  infeasible.  This case is handled by forcing the current occupant of the
  basis position into the nonbasic partition and replacing it with the
  logical for the constraint.  From there, it's a matter of deactivating a
  constraint with a basic logical.

  Activation of a constraint using a nonbasic logical isn't an issue in dylp.
  Should the need ever arise, the appropriate strategy would be to convert some
  nonbasic variable to basic at bound and add the logical as nonbasic at bound
  or superbasic.

  The bottom routine for constraint activation is dy_loadcon, which handles
  the business of translating a constraint from the original system frame of
  reference and installing it in the active system.  This is used bare during
  initialisation.  Immediately above is dy_actBLogPrimCon, which deals with
  basis and status issues once dylp has finished initialisation and into
  simplex.

  For deactivation, the bottom routine is deactBLogPrimCon.
*/

/*
  A few words about the algebra of constraint addition and deletion, as it
  relates to the PSE and DSE variables.

  First, suppose that no new variables are activated with a constraint, with
  the exception of the associated logical. It becomes basic with a
  coefficient of 1 (by definition; remember the constraints have been
  rewritten to convert >= constraints to <= constraints). It is, however,
  unfortunately true that added constraints may have non-zero coefficients
  for any of the existing basic variables.

  Suppose we're adding constraint k. The new basis has the form
    B' = [[B 0][a<k,B> 1]],
  and the inverse will be
    inv(B') = [[inv(B) 0][-a<k,B>inv(B) 1]].
  The cost vector c<B'> = (c<B> 0).
  
  For reduced costs, we have no changes, given that the objective coefficient
  for the slack variable is 0.  We do, however, recalculate the duals anyway,
  for accuracy (we will be refactoring to get the new basis inverse).

  For PSE,  the gamma<j> = ||abar~<j>||^2 will not change --- the projected
  edge abar~<j> is unchanged because the new slack is basic, hence not part
  of the reference frame.

  When we delete a constraint with a basic logical, PSE norms must be updated
  if the slack was nonbasic (hence added to the reference frame) at some
  point in the past.  Suppose the slack is basic in position k.  Removing it
  will change values of gamma<j> for all columns where abar~<k,j> has a
  nonzero.  The efficient way out is to reset the reference frame.  (It's an
  open question whether adjusting the norms would give better performance.)

  For DSE, the situation is interesting. You'd think that changing the basis
  by adding/deleting both rows (the constraints) and columns (the logicals)
  would change everything. But as the algebra above shows, the special form
  of the partitions means that the existing beta<i> do not change.
  Deactivation requires no changes at all. Activation requires that we
  calculate the norms for the new rows.

  In practice, when dylp adds constraints, it collects a list of inactive
  variables referenced by the constraints and then considers activating them,
  based on the target simplex phase (for primal, we want favourable
  (nonoptimal) reduced cost; for dual, unfavourable (optimal)). But this is
  handled by variable activation routines, which deal with the side effects.
*/



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




bool dy_loadcon (consys_struct *orig_sys, int i,
		 bool genvars, int *inactndxs)

/*
  This routine loads a constraint i from the original system (orig_sys) into
  the active system (dy_sys). It expects dy_origvars to be valid.
  
  If genvars == TRUE, new columns are generated as needed for eligible
  inactive variables x<j>. A variable is eligible if a<ij> != 0 and the
  variable is not fixed.

  If genvars is FALSE, only the coefficients a<ij> corresponding to already
  active variables will be installed when a<i> is loaded. If a vector is
  provided in inactndxs, it is used to return the indices of eligible
  variables x<j> which would have been activated if genvars were TRUE.

  WARNING: Use genvars == TRUE only for loading the initial constraint system.
	   The simple procedure used here assumes that an empty column is
	   generated when a variable x<j> is first encountered. It WILL FAIL
	   if there are non-zero coefficients for x<j> in constraints that are
	   already active.

	   When dylp is given an active variable specification, it will
	   establish the active variable set before loading constraints, and
	   genvars will be FALSE. If there's no active variable specification,
	   the constraints determine the active variables, and genvars will be
	   TRUE.

	   Once dylp is past loading the initial constraint system, there is
	   much more work to do to activate a variable. Use inactndxs to
	   return the candidates and call an appropriate variable activation
	   routine (viz. dy_varmgmt.c).

  The generation of the logical variable associated with the constraint will
  be handled automatically. The routine will correct dy_origvars and
  dy_actvars when an architectural variable is shifted to make room for this
  constraint's logical variable.  Other than dy_origvars and dy_actvars, no
  other data structures are updated.  The assumption is that this routine is
  used during initialization of dy_sys, before other lp data structures are
  established.

  Parameters:
    orig_sys:	The original constraint system.
    i:		Index of the constraint to be activated.
    genvars:	TRUE if inactive variables used in this constraint should
		be activated; FALSE otherwise
    inactndxs:  (i) Vector to hold inactive variable indices, or NULL
		(o) If genvars == FALSE, and inactndxs != NULL, the indices
		of inactive variables with nonzero coefficients in a<i> will
		be loaded into inactndxs. inactndxs[0] is set to the number
		of valid indices stored in inactndxs[1 .. inactndxs[0]]. The
		vector is assumed to be of adequate size.

  Returns: TRUE if the constraint is successfully installed, FALSE otherwise.
*/

{ int ndx,act_ndx,j,act_j,act_i ;
  int inact_ndx ;
  double rhsadj,rhscorr,act_rhs,act_rhslow ;
  flags statj ;
  pkvec_struct *ai,*aj ;
  pkcoeff_struct *aij ;
  bool retval ;
  const char *rtnnme = "dy_loadcon" ;

# ifndef DYLP_NDEBUG
  int print ;

  switch (dy_lp->phase)
  { case dyINIT:
    { print = dy_opts->print.setup ;
      break ; }
    case dyFORCEFULL:
    case dyADDCON:
    { print = dy_opts->print.conmgmt+1 ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# endif

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (genvars == TRUE && dy_lp->phase != dyINIT)
  { errmsg(1,rtnnme,__LINE__) ;
    return (FALSE) ; }
  if (i <= 0 || i > orig_sys->concnt)
  { errmsg(102,rtnnme,orig_sys->nme,"constraint",i,1,orig_sys->concnt) ;
    return (FALSE) ; }
  ndx = (orig_sys->concnt-dy_lp->sys.cons.unloadable) -
	(dy_lp->sys.cons.loadable+dy_sys->concnt) ;
  if (ndx != 0)
  { errmsg(444,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint",orig_sys->concnt,dy_lp->sys.cons.unloadable,
	   dy_lp->sys.cons.loadable,dy_sys->concnt,ndx) ;
    return (FALSE) ; }
  if (ACTIVE_CON(i))
  { char onmbuf[128] ;
    act_j = dy_origcons[i] ;
    (void) consys_nme(orig_sys,'c',i,TRUE,onmbuf) ;
    errmsg(431,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint",onmbuf,i,
	   consys_nme(dy_sys,'c',act_j,TRUE,NULL),act_j) ;
    return (FALSE) ; }
  if (!LOADABLE_CON(i))
  { errmsg(445,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint",consys_nme(orig_sys,'c',i,TRUE,NULL),i) ;
    return (FALSE) ; }
# endif

/*
  Prep work. Make a 0-length vector that we'll use to retrieve and create
  column headers. Make sure the row vector pointer is null so that getrow_pk
  will allocate one for us, and retrieve the constraint.
*/
  if (genvars == TRUE)
  { aj = pkvec_new(0) ; }
  else
  { aj = NULL ; }
  ai = NULL ;
  if (consys_getrow_pk(orig_sys,i,&ai) == FALSE)
  { errmsg(122,rtnnme,orig_sys->nme,
	   "row",consys_nme(orig_sys,'c',i,TRUE,NULL),i) ;
    if (aj != NULL) pkvec_free(aj) ;
    if (ai != NULL) pkvec_free(ai) ;
    return (FALSE) ; }
  retval = TRUE ;
/*
  Walk the constraint and convert the column indices from the original system
  frame of reference to the dylp system frame of reference. If we're allowed
  to activate, create empty columns as necessary for each eligible variable.
  If we're not allowed to activate, inactive variables are removed from a<i>
  by compressing the coefficient vector in place. rhsadj keeps track of the
  effect on the rhs of the constraint.
*/
  rhsadj = 0 ;
  inact_ndx = 0 ;
  act_ndx = 0 ;
  for (ndx = 0 ; ndx < ai->cnt ; ndx++)
  { aij = &ai->coeffs[ndx] ;
    j = aij->ndx ;
#   ifdef DYLP_PARANOIA
    if (j <= 0 || j > orig_sys->varcnt)
    { errmsg(102,rtnnme,orig_sys->nme,"variable",j,1,orig_sys->varcnt) ;
      retval = FALSE ;
      break ; }
    if (dy_origvars[j] == 0)
    { errmsg(1,rtnnme,__LINE__) ;
      retval = FALSE ;
      break ; }
#   endif
/*
  If x<j> is inactive, we need to activate it or step over it.  As explained
  in the opening comments, for activation we simply create an empty column.
  Fixed variables are never activated. All variables are continuous, as far
  as dylp is concerned.

  In spite of what error message 433 says, we need to allow for inactive SB
  variables here. The source is a warm start. During initialisation, all
  variables start out as inactive, hence nonbasic. The variables that the
  warm start thinks should be basic are assigned SB status.
*/
    if (INACTIVE_VAR(j))
    { statj = (flags) (-dy_origvars[j]) ;
#     ifdef DYLP_PARANOIA
      if (dy_lp->phase == dyINIT)
      { retval = flgon(statj,vstatNONBASIC|vstatEXOTIC) ; }
      else
      { retval = flgon(statj,vstatNONBASIC|vstatNBFR) ; }
      if (retval == FALSE)
      { errmsg(433,rtnnme,dy_sys->nme,
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "inactive",consys_nme(orig_sys,'v',j,TRUE,NULL),
	       j,dy_prtvstat(statj)) ;
	break ; }
#     endif
      if (genvars == TRUE && LOADABLE_VAR(j))
      { retval = consys_getcol_pk(orig_sys,j,&aj) ;
	if (retval == FALSE)
	{ errmsg(122,rtnnme,orig_sys->nme,"variable",
		 consys_nme(orig_sys,'v',j,TRUE,NULL),j) ;
	  break ; }
	retval =
	  consys_addcol_pk(dy_sys,vartypCON,aj,orig_sys->obj[j],
			   orig_sys->vlb[j],orig_sys->vub[j]) ;
	if ( retval == FALSE)
	{ errmsg(156,rtnnme,"variable",dy_sys->nme,aj->nme) ;
	  break ; }
	act_j = aj->ndx ;
#       ifndef DYLP_NDEBUG
	if (print >= 6)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n\t    activating %s variable %s (%d) to index %d, status %s.",
	     consys_prtvartyp(orig_sys->vtyp[j]),
	     consys_nme(orig_sys,'v',j,FALSE,NULL),j,act_j,
	     (statj == 0)?"unspecified":dy_prtvstat(statj)) ; }
#       endif
	dy_origvars[j] = act_j ;
	dy_actvars[act_j] = j ;
	dy_lp->sys.vars.loadable-- ; }
/*
  If activation is disallowed, note the contribution to the right-hand-side.
  If the variable is loadable, record the index in inactndxs. Then move on to
  the next variable.
*/
      else
      { if (inactndxs != NULL && LOADABLE_VAR(j)) inactndxs[++inact_ndx] = j ;
	switch (getflg(statj,vstatSTATUS))
	{ case vstatNBLB:
	  { rhscorr = aij->val*orig_sys->vlb[j] ;
	    break ; }
	  case vstatNBUB:
	  case vstatNBFX:
	  { rhscorr = aij->val*orig_sys->vub[j] ;
	    break ; }
	  default:
	  { rhscorr = 0 ;
	    break ; } }
	rhsadj -= rhscorr ;
#       ifndef DYLP_NDEBUG
	if (print >= 6)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		 "\n\t    skipping inactive %s variable %s (%d), status %s.",
		 consys_prtvartyp(orig_sys->vtyp[j]),
		 consys_nme(orig_sys,'v',j,FALSE,NULL),j,dy_prtvstat(statj)) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,", rhs += %g.",-rhscorr) ; }
#       endif
	continue ; } }
/*
  We're going to use this variable, so convert the index from the original
  system frame to the active system frame. We need to copy the value, too,
  just in case we're compressing.
*/
#   ifndef DYLP_NDEBUG
    if (print >= 5)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t  copying %s variable %s (%d) to index %d.",
		  consys_prtvartyp(orig_sys->vtyp[j]),
		  consys_nme(orig_sys,'v',j,FALSE,NULL),j,dy_origvars[j]) ; }
#   endif
    ai->coeffs[act_ndx].ndx = dy_origvars[j] ;
    ai->coeffs[act_ndx].val = aij->val ;
    act_ndx++ ; }

  if (aj != NULL) pkvec_free(aj) ;
  if (inactndxs != NULL) inactndxs[0] = inact_ndx ;
  if (retval == FALSE)
  { if (ai != NULL) pkvec_free(ai) ;
    return (FALSE) ; }
  ai->cnt = act_ndx ;
/*
  The constraint has been converted -- column indices are in the dy_sys frame
  of reference and inactive columns are removed. Add the row to dy_sys.
  Everything is an architectural constraint as far as dy_sys is concerned.

  Note that we should never look at slots in dy_actvars corresponding to
  logicals --- they don't exist in orig_sys. -INT_MAX should guarantee errors
  if we ever use the value.
*/
  act_rhs = orig_sys->rhs[i]+rhsadj ;
  if (orig_sys->ctyp[i] == contypRNG)
  { act_rhslow = orig_sys->rhslow[i]+rhsadj ; }
  else
  { act_rhslow = 0 ; }
  retval = consys_addrow_pk(dy_sys,'a',orig_sys->ctyp[i],
			    ai,act_rhs,act_rhslow,NULL,NULL) ;
  if (retval == FALSE)
  { errmsg(156,rtnnme,"constraint",dy_sys->nme,ai->nme) ;
    if (ai != NULL) pkvec_free(ai) ;
    return (FALSE) ; }
  act_i = ai->ndx ;
  pkvec_free(ai) ;
  dy_origcons[i] = act_i ;
  dy_actcons[act_i] = i ;
  dy_actvars[act_i] = -INT_MAX ;
/*
  Are there architectural variables? The addition of the logical will cause
  the architectural variable that used to occupy position act_i to be moved
  to dy_sys->varcnt, and we need to adjust dy_origvars accordingly.
  (dy_actvars is attached to dy_sys and has been updated automatically.)
*/
  if (dy_sys->archvcnt > 0)
  { act_j = dy_sys->varcnt ;
    j = dy_actvars[act_j] ;
#   ifdef DYLP_PARANOIA
    if (dy_origvars[j] != act_i)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; }
#   endif
    dy_origvars[j] = act_j ;
#   ifndef DYLP_NDEBUG
    if (print >= 6)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t    variable %s (%d) shifted from index %d",
		  consys_nme(dy_sys,'v',act_j,FALSE,NULL),act_j,act_i) ; }
#   endif
  }

# ifndef DYLP_NDEBUG
  if (print >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      %s %s (%d) copied to index %d",
	        consys_prtcontyp(dy_sys->ctyp[act_i]),
	        consys_nme(dy_sys,'c',act_i,FALSE,NULL),i,act_i) ; }
# endif
/*
  Bookkeeping.
*/
  dy_lp->sys.cons.loadable-- ;

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->cons.actcnt[i]++ ;
# endif

  return (TRUE) ; }




bool dy_actBLogPrimCon (consys_struct *orig_sys, int origi, int *inactvars)

/*
  Assume dy_sys has m constraints and n variables (logical plus architectural).
  This routine activates a primal constraint a<origi> to index i, where
  i = m+1. The basis is augmented with a basic logical x<i>. To make room for
  the logical x<i>, the first architectural x<m+1> is moved to index j = n+1.
  The constraint could be loose, tight, or violated, depending on the value
  of x<i>. What's important here is that we can easily generate a new basis by
  augmenting the existing basis with x<i>.

  Once dy_loadcon adds the constraint to dy_sys, we'll need to examine and
  correct the arrays dy_basis and dy_var2basis.

  Parameters:
    orig_sys: The original constraint system.
    origi:    The constraint to be activated.
    inactvars: (i) An array to hold indices of inactive variables referenced
	       by the constraint, or NULL if not desired
	       (o) Indices of inactive referenced variables, if requested.

  Returns: TRUE if activation succeeds, FALSE if there's an error.
*/

{ int i,j ;
  double lhsi,rhsi,rhslowi ;
  contyp_enum ctypi ;

  const char *rtnnme = "dy_actBLogPrimCon" ;

# ifdef DYLP_PARANOIA
/*
  A little paranoia. Check that origi is valid.
*/
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (origi <= 0 || origi > orig_sys->concnt)
  { errmsg(102,rtnnme,"original constraint",origi,1,orig_sys->concnt) ;
    return (FALSE) ; }
# endif

  ctypi = orig_sys->ctyp[origi] ;

# ifndef DYLP_NDEBUG
/*
  A little information, if the user wants it.
*/
  if (dy_opts->print.conmgmt >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    activating ") ;
    if (ctypi == contypRNG)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%g <= ",orig_sys->rhslow[origi]) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,"%s (%d) %s %g",
	        consys_nme(orig_sys,'c',origi,FALSE,NULL),origi,
	        consys_prtcontyp(ctypi),orig_sys->rhs[origi]) ; }
# endif
/*
  Load the constraint into dy_sys. No new architectural variables.
*/
  if (dy_loadcon(orig_sys,origi,FALSE,inactvars) == FALSE)
  { errmsg(430,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "activate","original constraint",
	   consys_nme(orig_sys,'c',origi,TRUE,NULL),origi) ;
    return (FALSE) ; }
/*
  Correct the basis. The new constraint has index i = m+1 and the logical
  x<i> becomes the basic variable. Since it's basic, it's not part of the PSE
  reference frame, and its reduced cost is zero. If there are architectural
  variables, the old x<m> has become x<j>. If it's basic, we need to correct
  dy_basis.
*/
  i = dy_sys->concnt ;
  j = dy_sys->varcnt ;
  dy_basis[i] = i ;
  dy_var2basis[i] = i ;
  dy_frame[i] = FALSE ;
  dy_cbar[i] = 0.0 ;
  if (j > dy_sys->concnt)
  { if (dy_var2basis[j] != 0)
      dy_basis[dy_var2basis[j]] = j ; }
/*
  Finally, we need to set the status of the logical. Evaluate the constraint
  and set it accordingly. Since x<i> is part of the constraint, set it to 0
  before the evaluation. Note that we need the loaded rhs values here to
  capture the correction for nonzero inactive values.
*/
  rhsi = dy_sys->rhs[i] ;
  dy_x[i] = 0 ;
  lhsi = consys_dotrow(dy_sys,i,dy_x) ;
  setcleanzero(lhsi,dy_tols->zero) ;
  if (abovebnd(lhsi,rhsi))
  { dy_status[i] = vstatBLLB ; }
  else
  if (atbnd(lhsi,rhsi))
  { if (ctypi == contypEQ)
    { dy_status[i] = vstatBFX ; }
    else
    { dy_status[i] = vstatBLB ; } }
  else
  if (ctypi != contypRNG)
  { dy_status[i] = vstatB ; }
  else
  { rhslowi = dy_sys->rhslow[i] ;
    if (belowbnd(lhsi,rhslowi))
    { dy_status[i] = vstatBUUB ; }
    else
    if (atbnd(lhsi,rhslowi))
    { dy_status[i] = vstatBUB ; }
    else
    { dy_status[i] = vstatB ; } }
/*
  And finally, a little paranoia. Generally, we should only be activating
  violated or tight constraints, so the logical should be at or outside its
  bounds. There are two exceptions: If we're trying to bound an unbounded
  primal, or we're forcing a full system.
*/
# ifdef DYLP_PARANOIA
  if (flgon(dy_status[i],vstatB) ||
      (dy_opts->con.actlvl = 0 && flgon(dy_status[i],vstatBLLB|vstatBUUB)))
  { if (dy_lp->phase == dyFORCEFULL ||
	(dy_lp->phase == dyADDCON &&
	 (dy_lp->lpret == lpSWING || dy_lp->lpret == lpUNBOUNDED)))
    { /* ok */ }
    else
    { warn(442,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   consys_nme(dy_sys,'c',i,FALSE,NULL),i,dy_prtvstat(dy_status[i]),
	   (dy_opts->con.actlvl == 0)?"violated":"tight or violated") ; } }
# endif

  return (TRUE) ; }


bool dy_actBLogPrimConList (consys_struct *orig_sys,
			    int cnt, int *ocndxs, int **p_inactvars)

/*
  This routine is a shell to call actBLogPrimCon for each of the indices in
  the vector ocndxs. It performs minimal error checking, relying on checking in
  actBLogPrimCon.

  If p_inactvars is non-NULL, the routine accumulates a list of indices of
  inactive variables referenced by the activated constraints. There is no
  attempt to remove duplicates, hence there's no good limit on the size of this
  vector. It may be realloc'd over the course of activation. We can bound the
  number of indices collected for a given constraint, however.

  Parameters:
    orig_sys:	The original constraint system
    cnt:	The number of indices in ocndxs
    ocndxs:	Vector of constraint indices (0-based)
    p_inactvars: (i) If p_inactvars == NULL, collection of indices of
		 referenced variables is suppressed. If *p_inactvars == NULL,
		 a vector will be allocated.  If *p_inactvars != NULL,
		 (*p_inactvars)[0] must contain the allocated size.
		 (o) If collection of indices is taking place, the vector
		 of indices is returned, with (*p_inactvars)[0] set to the
		 number of indices collected. If no referenced variables were
		 encountered, *p_inactvars will remain NULL if no vector was
		 supplied.

  Returns: TRUE if all constraints are successfully activated, FALSE otherwise.
*/

{ int j,k,ndx,act_n,inact_n ;
  int *onecon,onecon_cnt,*collection,coll_cnt,coll_sze ;
  bool *seen ;
  bool with_vars,retval ;
  const char *rtnnme = "dy_actBLogPrimConList" ;

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (ocndxs == NULL)
  { errmsg(2,rtnnme,"ocndxs") ;
    return (FALSE) ; }
  if (cnt <= 0 || cnt > orig_sys->concnt)
  { errmsg(5,rtnnme,"cnt",cnt) ;
    return (FALSE) ; }
# endif

  retval = TRUE ;
/*
  Are we collecting indices of referenced variables? If so, set up to do it.
  We can guarantee that any one constraint will reference no more than
  (inact_n - act_n) variables. To avoid duplicates, we need to keep track
  of the indices already collected. Then we can guarantee that the collection
  itself will be no more than (inact_n - act_n) variables.
*/
  if (p_inactvars != NULL)
  { with_vars = TRUE ;
    act_n = dy_sys->archvcnt ;
    inact_n = orig_sys->archvcnt ;
    coll_sze = inact_n-act_n+1 ;
    seen = (bool *) CALLOC((inact_n+1),sizeof(bool)) ;
    if (*p_inactvars == NULL)
    { collection = (int *) MALLOC(coll_sze*sizeof(int)) ; }
    else
    { collection = *p_inactvars ;
      coll_sze = collection[0] ; }
    collection[0] = 0 ;
    coll_cnt = 0 ;
    onecon = (int *) MALLOC(coll_sze*sizeof(int)) ; }
  else
  { with_vars = FALSE ;
    onecon = NULL ;
    collection = NULL ;
    coll_cnt = -1 ;
    seen = NULL ; }
/*
  Open a loop to request activation of each constraint in ocndxs. On return,
  concatenate the indices from onecon to collection.
*/
  for (k = 0 ; k < cnt && retval == TRUE ; k++)
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.conmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    activating constraint %s (%d)",
		  consys_nme(orig_sys,'c',ocndxs[k],TRUE,NULL),ocndxs[k]) ;
      if (with_vars == FALSE || dy_opts->print.conmgmt < 4)
      { dyio_outchr(dy_logchn,dy_gtxecho,'.') ; } }
#   endif
    retval = dy_actBLogPrimCon(orig_sys,ocndxs[k],onecon) ;
    if (retval == FALSE)
    { errmsg(430,rtnnme,
	     orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "activate","constraint",
	     consys_nme(orig_sys,'c',ocndxs[k],TRUE,NULL),ocndxs[k]) ; }
    if (with_vars == TRUE)
    { onecon_cnt = onecon[0] ;
      for (ndx = 1 ; ndx <= onecon_cnt ; ndx++)
      { j = onecon[ndx] ;
	if (seen[j] == FALSE)
	{ collection[++coll_cnt] = j ;
	  seen[j] = TRUE ; } }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.conmgmt >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,", %d referenced variables queued.",
		    coll_cnt-collection[0]) ;
	collection[0] = coll_cnt ; }
#     endif
    } }
/*
  Clean up from collection of indices of referenced variables.
*/
  if (with_vars == TRUE)
  { collection[0] = coll_cnt ;
    if (*p_inactvars == NULL)
    { if (coll_cnt == 0)
      { FREE(collection) ; }
      else
      { *p_inactvars = collection ; } }
    if (onecon != NULL) FREE(onecon) ;
    if (seen != NULL) FREE(seen) ; }

# ifdef DYLP_PARANOIA
  if (retval == TRUE)
  { retval = dy_chkdysys(orig_sys) ; }
# endif

  return (retval) ; }




bool dy_deactNBLogPrimCon (consys_struct *orig_sys, int i)

/*
  This routine deactivates a tight primal constraint a<i> with a nonbasic
  logical x<i>. The only reason we want to do this is because we're trying to
  force dual feasibility and we need to deactivate the logical.
  Realistically, though, we're deleting the (dual) column associated with a
  nonzero dual, so dual feasibility will be a stroke of luck. But
  deactivating a primal constraint can't harm primal feasibility, eh? The
  routine takes the attitude that if we have primal feasibility, it'll try to
  retain it and force the basic variable out with superbasic status.

  What's important here is that the job is complicated considerably by the
  presence of some other variable x<j> basic in pos'n i. The easiest way to
  overcome this problem is to alter reality.

  We'll take the variable x<j> and mark it nonbasic at bound, unless we have
  primal feasibility, in which case we'll opt for SB status in an attempt to
  retain feasibility.  Then we'll take x<i>, the nonbasic logical for a<i>,
  and mark it basic in pos'n i.  Then we'll call deactBLogPrimCon to do the
  heavy lifting.

  Parameters:
    orig_sys: The original constraint system.
    i:	      The constraint to be deactivated.

  orig_sys is used only for printing, statistics, and paranoid checks, but
  it gives a nice symmetry with actBLogPrimCon.

  Returns: TRUE if deactivation succeeds, FALSE if there's an error.
*/

{ int j,m,n ;
  double lbj,ubj,valj ;
  flags stati,statj ;

  const char *rtnnme = "dy_deactNBLogPrimCon" ;

/*
  A little paranoia, mixed with prep. Check that i is valid and nonbasic.
*/
  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
# ifdef DYLP_PARANOIA
  if (i <= 0 || i > m)
  { errmsg(102,rtnnme,"constraint",i,1,m) ;
    return (FALSE) ; }
# endif
  stati = getflg(dy_status[i],vstatSTATUS) ;
# ifdef DYLP_PARANOIA
  if (flgoff(stati,vstatNBLB|vstatNBUB|vstatNBFX))
  { errmsg(437,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   consys_nme(dy_sys,'c',i,TRUE,NULL),i,dy_prtvstat(stati)) ;
    return (FALSE) ; }
# endif
/*
  Grab hold of x<j>, the variable that's basic in pos'n i, and push it out into
  the nonbasic partition. If we're unlucky, it'll go out with SB status. If
  we're doubly unlucky, it's part of the reference frame and we'll need to
  reinitialize the PSE reference frame.

  SB status only applies in primal phase II. If we don't have feasibility, no
  need to use exotic status to preserve it. Just pick a finite bound.
*/
  j = dy_basis[i] ;
  statj = getflg(dy_status[j],vstatSTATUS) ;
  lbj = dy_sys->vlb[j] ;
  ubj = dy_sys->vub[j] ;
  switch (statj)
  { case vstatB:
    { if (dy_lp->simplex.active == dyPRIMAL2)
      { statj = vstatSB ;
	valj = dy_x[j] ; }
      else
      { if (lbj > -dy_tols->inf)
	{  statj = vstatNBLB ;
	   valj = lbj ; }
	else
	{ statj = vstatNBUB ;
	  valj = ubj ; } }
      break ; }
    case vstatBUB:
    case vstatBUUB:
    { statj = vstatNBUB ;
      valj = ubj ;
      break ; }
    case vstatBLB:
    case vstatBLLB:
    { statj = vstatNBLB ;
      valj = lbj ;
      break ; }
    case vstatBFX:
    { statj = vstatNBFX ;
      valj = lbj ;
      break ; }
    case vstatBFR:
    { statj = vstatNBFR ;
      valj = dy_x[j] ; 
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
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
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      swapping %s (%d) %s -> ",
	        consys_nme(dy_sys,'v',i,FALSE,NULL),i,
		dy_prtvstat(dy_status[i])) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"%s ",dy_prtvstat(stati)) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"<=> %s (%d) %s -> ",
	        consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		dy_prtvstat(dy_status[j])) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"%s.",dy_prtvstat(statj)) ; }
# endif
/*
  Tweak dy_basis, dy_var2basis, and dy_status to do the swap. The entries for
  x<i> are going to disappear momentarily, but deactBLogPrimCon will want to
  examine them. If x<j> is part of the reference frame, moving it to the
  nonbasic partition changes the composition of the PSE norms. The change in
  basis changes the DSE norms.
*/
  dy_var2basis[j] = 0 ;
  dy_status[j] = statj ;
  if (dy_frame[j] == TRUE) dy_lp->simplex.init_pse = TRUE ;
  dy_x[j] = valj ;
  dy_lp->simplex.init_dse = TRUE ;
  dy_basis[i] = i ;
  dy_var2basis[i] = i ;
  dy_status[i] = stati ;
/*
  Reality is altered. Call deactBLogPrimCon to do the heavy lifting.
*/
  return (dy_deactBLogPrimCon(orig_sys,i)) ; }



bool dy_deactBLogPrimCon (consys_struct *orig_sys, int i)

/*
  This routine deactivates a primal constraint a<i> with a basic logical x<i>.
  The constraint could be loose, tight, or violated, depending on the value
  of the logical. What's important here is that we can easily fix up the
  basis because we'll be deleting one basis position and one basic variable.

  To make cleanup easier, if x<i> does not already occupy basis position i,
  do a swap to make it so. Assuming m constraints and n variables (including
  logicals) before deletion, the pattern of motion becomes

    a<i> --> deleted
    a<m> --> a<i>

    x<i> --> deleted
    x<m> --> x<i>
    x<n> --> x<m>

  Once the constraint is deleted, we'll need to examine and correct the
  arrays dy_basis and dy_var2basis, as well as dy_origcons and dy_origvars.

  It's entirely possible to find that there are no constraints left after
  removing loose constraints. It follows that there will be no variables, and
  in fact we can end up with a completely empty constraint system. This is
  not uncommon deep in a branch-and-bound search tree, when dylp can be handed
  a system with many fixed variables.

  If we delete a logical that's part of the current PSE index, we need to
  correct the projected column norms. But ... in the context of dylp, there's
  a good chance we'll proceed from constraint deletion into dual simplex. So
  we'll just indicate that the column norms need to be initialized when we
  reenter primal simplex.

  Parameters:
    orig_sys: The original constraint system.
    i:	      The constraint to be deactivated.

  orig_sys is used only for printing, statistics, and paranoid checks, but
  it gives a nice symmetry with actBLogPrimCon.

  Returns: TRUE if deactivation succeeds, FALSE if there's an error.
*/

{ int j,k,m,n,bposi,origi ;
  flags stati ;

  const char *rtnnme = "dy_deactBLogPrimCon" ;

/*
  A little paranoia, mixed with prep. Check that i is valid, that we can locate
  a<i> in the original system, and that the logical is basic. Also check the
  constraint invariant.
*/
  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;

# ifdef DYLP_PARANOIA
  if (i <= 0 || i > m)
  { errmsg(102,rtnnme,"constraint",i,1,m) ;
    return (FALSE) ; }
# endif

  stati = dy_status[i] ;
  bposi = dy_var2basis[i] ;

# ifdef DYLP_PARANOIA
  if (flgoff(stati,vstatBASIC))
  { errmsg(436,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   consys_nme(dy_sys,'c',i,TRUE,NULL),i,dy_prtvstat(stati)) ;
    return (FALSE) ; }
  if (bposi <= 0 || dy_basis[bposi] != i)
  { errmsg(330,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   consys_nme(dy_sys,'v',i,FALSE,NULL),i,dy_prtvstat(stati),
	   bposi,dy_basis[bposi]) ;
    return (FALSE) ; }
# endif

  origi = dy_actcons[i] ;

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (origi <= 0 || origi > orig_sys->concnt)
  { errmsg(102,rtnnme,"original constraint",origi,1,orig_sys->concnt) ;
    return (FALSE) ; }
  k = (orig_sys->concnt-dy_lp->sys.cons.unloadable) -
	(dy_lp->sys.cons.loadable+dy_sys->concnt) ;
  if (k != 0)
  { errmsg(444,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint",orig_sys->concnt,dy_lp->sys.cons.unloadable,
	   dy_lp->sys.cons.loadable,dy_sys->concnt,k) ;
    return (FALSE) ; }
# endif

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->cons.deactcnt[origi]++ ;
# endif
/*
  Before this entry is deleted, check if x<i> is part of the PSE reference
  frame. If it is, indicate that we'll need to initialize the projected column
  norms the next time we enter primal simplex.
*/
  if (dy_frame[i] == TRUE) dy_lp->simplex.init_pse = TRUE ;
/*
  We've verified that the constraint is suitable for deletion. To make the
  basis touchup easier, ensure that x<i> occupies basis position i by
  swapping with the current occupant. We don't need to set the new values in
  dy_basis and dy_var2basis (they're about to disappear, along with the ith
  entry of all other vectors attached to dy_sys).
*/
  if (bposi != i)
  { k = dy_basis[i] ;
    dy_basis[bposi] = k ;
    dy_var2basis[k] = bposi ; }
/*
  In a similar vein, if x<m> is basic, and we will shift a<m> to fill the
  hole left by deleting a<i>, do a swap to ensure x<m> is basic in pos'n m.
  Then the automatic compression of basis and var2basis will move them in
  tandem. One less complication to deal with after the fact.
*/
  if (i < m)
  { k = dy_var2basis[m] ;
    if (k > 0 && k != m)
    { j = dy_basis[m] ;
      dy_var2basis[j] = k ;
      dy_basis[k] = j ;
      dy_var2basis[m] = m ;
      dy_basis[m] = m ; } }
/*
  At this point, we can say the following:
    + Logical x<i> occupies basis pos'n i and will be deleted along with
      constraint a<i>
    + If logical x<m> is basic, it occupies pos'n m.
      If logical x<m> is nonbasic, then some variable x<q> occupies basis
      pos'n m. It's possible that q = n.
    + If architectural x<n> is basic, it occupies pos'n p. It's possible
      p = m.

  Now delete the constraint a<i>. Mark it as inactive in origcons and then
  perform the actual deletion from dy_sys. Deletion of the associated logical
  x<i> will occur automatically, as will shifts to compact the constraint
  system and the various attached arrays.

  In particular:
    + dy_basis and dy_actcons will be compacted as a<i> is deleted.  The
      information in pos'n m will be moved to pos'n i.
    + dy_var2basis and dy_actvars will be compacted as x<i> is deleted. The
      information in pos'n m will be moved to pos'n i, and the information in
      pos'n n will be moved to pos'n m.
*/
  k = dy_actcons[i] ;
  MARK_INACTIVE_CON(k) ;
  if (consys_delrow(dy_sys,i) == FALSE)
  { errmsg(112,rtnnme,dy_sys->nme,"delete","constraint",
	   consys_nme(dy_sys,'c',i,FALSE,NULL),i) ;
    return (FALSE) ; }

# ifdef DYLP_PARANOIA
/*
  A few checks to make sure that things are in the expected places. If there
  are constraints remaining and a<i> was not the last constraint, then the
  old constraint a<m> should occupy position i. Similarly, if there are
  architecturals remaining, the old variable x<n> should now be x<m>.
*/
  if (i <= dy_sys->concnt)
  { k = dy_actcons[i] ;
    if (dy_origcons[k] != m)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
  if (m <= dy_sys->varcnt)
  { k = dy_actvars[m] ;
    if (dy_origvars[k] != n)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# endif

/*
  Now to repair the collateral damage. The trick is to avoid fixing things
  more than once.  Let's start with the direct effects of constraint motion.
  If a<i> was not the last constraint, a<m> has been shifted to position i,
  and the associated slack x<m> is now x<i>. We need to correct origcons,
  var2basis, and basis. Avoid corrections that involve x<n>; we'll get to them
  in the next block of code.
*/
  if (i <= dy_sys->concnt)
  { k = dy_actcons[i] ;
    dy_origcons[k] = i ;
    k = dy_basis[i] ;
    if (k == m)
    { dy_basis[i] = i ;
      dy_var2basis[i] = i ; }
    else
    if (k != n)
    { dy_var2basis[k] = i ; } }
/*
  If we shifted an architectural (which will happen unless there were none to
  shift), then we need to correct origvars, var2basis, and basis. Index shift
  is n -> m, and var2basis is compressed. So we check var2basis[m] to see if
  our architectural is basic. If it was basic in pos'n m, we need to change
  that to i, where the information now lives. Then go and correct the basis
  entry to contain our variable's new index.
*/
  if (m <= dy_sys->varcnt)
  { k = dy_actvars[m] ;
    dy_origvars[k] = m ;
    bposi = dy_var2basis[m] ;
    if (bposi != 0)
    { if (bposi == m)
      { dy_basis[i] = m ;
	dy_var2basis[m] = i ; }
      else
      { dy_basis[bposi] = m ; } } }
/*
  And finally, a little bookkeeping.
*/
  dy_lp->sys.cons.loadable++ ;

/*
  We're done. Do some printing, if requested, then return.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tactive now %d x %d (%d+%d).",
	        dy_sys->concnt,dy_sys->varcnt,
		dy_sys->archvcnt,dy_sys->concnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\tconstraint %s (%d) and logical deleted from pos'n %d.",
	        consys_nme(orig_sys,'c',origi,FALSE,NULL),origi,i) ;
    if (i <= dy_sys->concnt)
    { k = dy_actcons[i] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tconstraint %s (%d) shifted from pos'n %d, ",
		  consys_nme(orig_sys,'c',k,FALSE,NULL),k,m) ;
      k = dy_basis[i] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"basis[%d] = %s (%d)",
		  i,consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
      k = dy_var2basis[i] ;
      if (k != 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tbasis pos'n %d updated to %s (%d).",
		    k,consys_nme(dy_sys,'v',dy_basis[k],FALSE,NULL),
		    dy_basis[k]) ; } }
    if (m <= dy_sys->varcnt)
    { k = dy_actvars[m] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tarchitectural %s (%d) shifted from pos'n %d.",
		  consys_nme(orig_sys,'v',k,FALSE,NULL),k,n) ;
      k = dy_var2basis[m] ;
      if (k != 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tbasis pos'n %d updated to %s (%d).",
		    k,consys_nme(dy_sys,'v',dy_basis[k],FALSE,NULL),
		    dy_basis[k]) ; } } }
# endif

  return (TRUE) ; }


static bool deactBLogPrimConList (consys_struct *orig_sys,
				   int cnt, int *acndxs)
/*
  This routine is purely a shell to call deactBLogPrimCon for each of the
  indices in the vector acndxs. It performs minimal error checking, relying
  on checking in deactBLogPrimCon.

  Parameters:
    orig_sys:	The original constraint system
    cnt:	The number of indices in acndxs
    acndxs:	Vector of constraint indices (0-based)

  orig_sys is used only for printing and paranoid checks, but it gives a nice
  symmetry with dy_actBLogPrimConList.

  Returns: TRUE if all constraints are successfully activated, FALSE otherwise
*/

{ int k ;
  bool retval ;
  const char *rtnnme = "deactBLogPrimConList" ;

# ifdef DYLP_PARANOIA
  if (acndxs == NULL)
  { errmsg(2,rtnnme,"acndxs") ;
    return (FALSE) ; }
  if (cnt <= 0 || cnt > dy_sys->concnt)
  { errmsg(5,rtnnme,"cnt",cnt) ;
    return (FALSE) ; }
# endif

/*
  To make sure consys doesn't shift constraints out from under us, we need to
  delete in nonincreasing order.
*/
  qsort(&acndxs[0],cnt,sizeof(int),intcompare) ;
  
  retval = TRUE ;
  for (k = 0 ; k < cnt && retval == TRUE ; k++)
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.conmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    deactivating constraint %s (%d)",
		  consys_nme(dy_sys,'c',acndxs[k],TRUE,NULL),acndxs[k]) ; }
#   endif
    retval = dy_deactBLogPrimCon(orig_sys,acndxs[k]) ;
    if (retval == FALSE)
    { errmsg(430,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "deactivate","constraint",
	     consys_nme(dy_sys,'c',acndxs[k],TRUE,NULL),acndxs[k]) ; } }

# ifdef DYLP_PARANOIA
  if (retval == TRUE)
  { retval = dy_chkdysys(orig_sys) ; }
# endif

  return (retval) ; }



static int scanPrimConStdAct (consys_struct *orig_sys, int **p_ocndxs)

/*
  This routine scans the original constraint system looking for
  inactive constraints to add to the active system. There are two settings
  for opts.con.actlvl:
    0: activate violated constraints
    1: activate tight and violated constraints

  Parameters:
    orig_sys:	The original constraint system
    p_ocndxs:	(i) empty vector to hold constraint indices; assumed
		    sufficiently large if non-NULL; if NULL, allocated if
		    necessary
		(o) indices of constraints to be activated; may not be
		    allocated if no constraints are identified

  Returns: number of candidates for activation, -1 if error.
*/

{ int i,j,k,m,n,actcnt,cand_limit ;
  int *ocndxs ;
  double *orig_x,*orig_rhs,*orig_rhslow,*orig_vub,*orig_vlb ;
  double lhsi,rhsi,rhslowi ;
  contyp_enum *orig_ctyp,ctypi ;
  flags statj ;
  bool activate ;

  const char *rtnnme = "scanPrimConStdAct" ;

# ifdef DYLP_PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (-1) ;  }
  if (p_ocndxs == NULL)
  { errmsg(2,rtnnme,"&ocndxs") ;
    return (-1) ; }
# endif

  m = orig_sys->concnt ;
  n = orig_sys->varcnt ;

/*
  Did the client supply a vector for candidate indices? If not, make one.

  We shouldn't be here if there's no room to activate. Check this if we're
  paranoid.
*/
  cand_limit = dy_lp->sys.cons.loadable ;
# ifdef DYLP_PARANOIA
  if (cand_limit == 0)
  { errmsg(1,rtnnme,__LINE__) ;
    return (-1) ; }
# endif
  if (dy_opts->con.actlim > 0)
  { cand_limit = minn(dy_opts->con.actlim,cand_limit) ; }
  if (*p_ocndxs == NULL)
  { ocndxs = (int *) MALLOC(cand_limit*sizeof(int)) ; }
  else
  { ocndxs = *p_ocndxs ; }
/*
  Create a solution vector that matches orig_sys, to make the scan a bit more
  efficient.
*/
  orig_vub = orig_sys->vub ;
  orig_vlb = orig_sys->vlb ;
  orig_x = (double *) CALLOC((n+1),sizeof(double)) ;
  for (j = 1 ; j <= n ; j++)
  { k = dy_origvars[j] ;
    if (k > 0)
    { orig_x[j] = dy_x[k] ; }
    else
    { statj = (flags) -k ;
#     ifdef DYLP_PARANOIA
      if (flgoff(statj,vstatNONBASIC|vstatNBFR))
      { errmsg(433,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "inactive",consys_nme(orig_sys,'v',j,TRUE,NULL),j,
	       dy_prtvstat(statj)) ;
	if (orig_x != NULL) FREE(orig_x) ;
	if (*p_ocndxs == NULL) FREE(ocndxs) ;
	return (-1) ; }
#     endif
      if (flgon(statj,vstatNBUB))
      { orig_x[j] = orig_vub[j] ; }
      else
      if (flgon(statj,vstatNBLB|vstatNBFX))
      { orig_x[j] = orig_vlb[j] ; } } }
/*
  Now we can step through the constraints. Evaluate each loadable inactive
  constraint and check to see if it's violated.
*/
  orig_ctyp = orig_sys->ctyp ;
  orig_rhs = orig_sys->rhs ;
  orig_rhslow = orig_sys->rhslow ;
  actcnt = 0 ;
  for (i = 1 ; i <= m && actcnt < cand_limit ; i++)
  { if (!LOADABLE_CON(i)) continue ;
    ctypi = orig_ctyp[i] ;
    lhsi = consys_dotrow(orig_sys,i,orig_x) ;
    setcleanzero(lhsi,dy_tols->zero) ;
/*
  Check the lhs against the rhs. There are two levels of activation, specified
  by opts.con.actlvl:
    * strict (0) activates only when lhs < rhslow or lhs > rhs
    * tight (1) activates when lhs <= rhslow or lhs >= rhs
*/
    rhsi = orig_rhs[i] ;
    if (ctypi == contypRNG)
    { rhslowi = orig_rhslow[i] ; }
    else
    if (ctypi == contypEQ)
    { rhslowi = rhsi ; }
    else
    { rhslowi = -dy_tols->inf ; }
    switch (dy_opts->con.actlvl)
    { case 0:
      { if (abovebnd(lhsi,rhsi) || belowbnd(lhsi,rhslowi))
	{ activate = TRUE ; }
	else
	{ activate = FALSE ; }
	break ; }
      case 1:
      { if (!(belowbnd(lhsi,rhsi) && abovebnd(lhsi,rhslowi)))
	{ activate = TRUE ; }
	else
	{ activate = FALSE ; }
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	if (orig_x != NULL) FREE(orig_x) ;
	if (*p_ocndxs == NULL) FREE(ocndxs) ;
	return (-1) ; } }
#   ifndef DYLP_NDEBUG
    if (activate == FALSE)
    { if (dy_opts->print.conmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    skipping %s constraint %s (%d), %g <= %g <= %g.",
		    consys_prtcontyp(orig_ctyp[i]),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i,
		    rhslowi,lhsi,rhsi) ; } }
    else
    { if (dy_opts->print.conmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    queued %s constraint %s (%d),",
		    consys_prtcontyp(orig_ctyp[i]),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i) ;
	if (abovebnd(lhsi,rhsi))
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " lhs - rhs = %g - %g = %g, tol %g.",
		      lhsi,rhsi,lhsi-rhsi,dy_tols->zero*(1+fabs(rhsi))) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " rhslow - lhs = %g - %g = %g, tol %g.",
		      rhslowi,lhsi,rhslowi-lhsi,
		      dy_tols->zero*(1+fabs(rhslowi))) ; } } }
#   endif

    if (activate == TRUE) ocndxs[actcnt++] = i ; }
  if (orig_x != NULL) FREE(orig_x) ;
/*
  If we supplied ocndxs and found no candidates to activate, free it.
*/
  if (*p_ocndxs == NULL)
  { if (actcnt == 0)
    { FREE(ocndxs) ; }
    else
    { *p_ocndxs = ocndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d %s constraints for activation.",actcnt,
	        dy_opts->con.actlvl == 0?"violated":"tight or violated") ; }
# endif

  return (actcnt) ; }



static int scanPrimConStdDeact (int **p_acndxs)

/*
  This routine scans the active constraint system for constraints that are
  candidates for deactivation. The primary criterion is that the logical for
  the constraint must be basic. Within this group of candidates, selection
  is based on the value of opts.con.deactlvl:

    * 0: (normal) Queues for deactivation inequalities a<i> which are
	 strictly loose.  Logicals x<i> must be strictly within bound (status
	 vstatB). The notion is that we've left these constraints behind and
	 will never return.
    * 1: (aggressive) Adds inequalities a<i> which are tight but the value of
	 the associated dual variable y<i> = 0.  (Logicals x<i> are at bound
	 with status vstatBUB or vstatBLB.) The notion is that these
	 inequalities aren't contributing (in terms of satisfying the dual,
	 yA >= c) and are more than likely just loitering about being
	 degenerate.
    * 2: (fanatic) Adds equalities a<i> with y<i> = 0. (Logicals
	 (artificials) x<i> have status vstatBFX.) Aimed at purging
	 constraints in large set covering problems but has relatively little
	 effect because it's hard to catch x<i> with status vstatBFX. Likely
	 should consider x<i> with status NBFX. Needs more thought.

  Logicals always have at least one finite bound and  should never have
  status vstatBFR. Logicals for a range constraint can have two finite
  bounds, but they should never be equal, so they can never have status
  vstatBFX. The logical (artificial) for a satisfied equality can only have
  status vstatBFX, but we'd be lucky to see it before dylp's pivoting
  priorities move it to the nonbasic partition, never to return. All of this
  means that we don't actually have to check constraint types, just look for
  the appropriate status codes.

  Parameters:
    p_acndxs:	(i) empty vector to hold constraint indices; assumed
		    sufficiently large if non-NULL; if NULL, allocated as
		    necessary
		(o) indices of constraints to be deactivated; may not be
		    allocated if no constraints are identified

  Returns: number of candidates for deactivation; -1 if error. Errors are
	   possible only if we're paranoid.
*/

{ int j,m,purgecnt ;
  int *acndxs ;
  flags statj ;
  bool purge ;

  const char *rtnnme = "scanPrimConStdDeact" ;

# ifdef DYLP_PARANOIA
  if (p_acndxs == NULL)
  { errmsg(2,rtnnme,"&acndxs") ;
    return (-1) ; }
# endif

  m = dy_sys->concnt ;
  if (*p_acndxs == NULL)
  { acndxs = (int *) MALLOC(m*sizeof(int)) ; }
  else
  { acndxs = *p_acndxs ; }
/*
  Open a loop to search for candidates for deactivation. It's pretty
  straightforward, as all we need to do is examine the status of the logical
  and (perhaps) the dual.
*/
  purgecnt = 0 ;
  for (j = 1 ; j <= m ; j++)
  { statj = dy_status[j] ;
    purge = FALSE ;
    if (flgon(statj,vstatB|vstatBLB|vstatBUB|vstatBFX))
    { switch (dy_opts->con.deactlvl)
      { case 0: /* normal */
	{ if (flgon(statj,vstatB)) purge = TRUE ;
	  break ; }
	case 1: /* aggressive */
	{ if (flgon(statj,vstatB) ||
	      (flgon(statj,vstatBLB|vstatBUB) && dy_y[j] == 0))
	    purge = TRUE ;
	  break ; }
	case 2: /* fanatic */
	{ if (flgon(statj,vstatB) ||
	      (flgon(statj,vstatBLB|vstatBUB|vstatBFX) && dy_y[j] == 0))
	    purge = TRUE ;
	  break ; }
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  if (*p_acndxs == NULL && acndxs != NULL) FREE(acndxs) ;
	  return (-1) ; } } }

#   ifndef DYLP_NDEBUG
    if (purge == FALSE)
    { if (dy_opts->print.conmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    skipped %s %s (%d), ",
		    consys_prtcontyp(dy_sys->ctyp[j]),
		    consys_nme(dy_sys,'c',j,TRUE,NULL),j) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "%s (%d) = %g, status %s, basis pos'n %d.",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j,
		    dy_x[j],dy_prtvstat(statj),dy_basis[j]) ; } }
    else       
    { if (dy_opts->print.conmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    queued %s %s (%d), ",
		    consys_prtcontyp(dy_sys->ctyp[j]),
		    consys_nme(dy_sys,'c',j,TRUE,NULL),j) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "%s (%d) = %g, status %s, basis pos'n %d.",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j,
		    dy_x[j],dy_prtvstat(statj),dy_basis[j]) ; } }
#   endif

    if (purge == TRUE)
    { acndxs[purgecnt++] = j ; } }

  if (*p_acndxs == NULL)
  { if (purgecnt <= 0)
    { FREE(acndxs) ; }
    else
    { *p_acndxs = acndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 1)
  { const char *strat ;
    switch (dy_opts->con.deactlvl)
    { case 0:
      { strat = "normal" ;
	break ; }
      case 1:
      { strat = "aggressive" ;
	break ; }
      case 2:
      { strat = "fanatic" ;
	break ; }
      default:
      { strat = "invalid" ;
	break ; } }
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  %s scan queued %d constraints for deactivation.",
	        strat,purgecnt) ; }
# endif

  return (purgecnt) ; }



int dy_deactivateCons (consys_struct *orig_sys)

/*
  This routine coordinates normal constraint deactivation in phase dyPURGECON.
  In addition to the actual scan and deactivation, it sees to rebuilding the
  basis and solution. The heavy lifting is performed in scanPrimConStdDeact
  and deactBLogPrimCon.

  See the comments at the head of the file for the effects on PSE and DSE
  norms.

  Parameters:
    orig_sys:	The original constraint system

  orig_sys is used only for debug printing and paranoid checks. Still, it's
  really convenient to have it as a parameter and it provides a nice symmetry
  with activateCons.

  Returns: number of constraints deactivated; -1 if there's an error.
*/

{ int *candidates,cand_cnt ;
  dyret_enum factorresult ;
  flags factorflags ;
  int retval ;

  const char *rtnnme = "dy_deactivateCons" ;

  retval = -1 ;
/*
  Call scanPrimConStdDeact to return a list of candidates for deactivation.
  If we find candidates, try to deactivate them.
*/
  candidates = NULL ;
  cand_cnt = scanPrimConStdDeact(&candidates) ;
  if (cand_cnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint","normal deactivation") ; }
  else
  if (cand_cnt > 0)
  { if (deactBLogPrimConList(orig_sys,cand_cnt,candidates) == TRUE)
/*
  If we purged constraints and there are still some left, we need to refactor
  and recalculate the primal and dual variables. It has happened that the only
  tight constraints are variables at bound, so that no explicit constraints
  remain.

  Arguably this is overkill --- in most cases, we'll move to dyADDCON, add
  constraints, and refactor again.
*/
    { if (dy_sys->concnt > 0)
      { factorflags = ladDUALS|ladPRIMALS ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n      refactoring ...") ;
#       endif
	factorresult = dy_factor(&factorflags) ;
	switch (factorresult)
	{ case dyrOK:
	  case dyrPATCHED:
	  { retval = cand_cnt ;
#           ifndef DYLP_NDEBUG
	    if (dy_opts->print.conmgmt >= 3)
	    { if (factorresult == dyrOK)
		dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ;
	      else
		dyio_outfmt(dy_logchn,dy_gtxecho,"\n    patched.") ; }
#           endif
	    break ; }
	  default:
	  { retval = -1 ;
#           ifndef DYLP_NDEBUG
	    if (dy_opts->print.conmgmt >= 3)
	      dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#           endif
	    break ; } } }
      else
      { retval = cand_cnt ; } }
    else
    { retval = -1 ; } }
  else
  { retval = cand_cnt ; }

  if (candidates != NULL) FREE(candidates) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  %d constraints deactivated.",cand_cnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  constraint system %s now %d x %d (%d + %d).",
	        dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
	        dy_sys->logvcnt) ; }
# endif

  return (retval) ; }



int dy_activateCons (consys_struct *orig_sys, bool with_vars)

/*
  This routine coordinates normal constraint activation in phase dyADDCON. In
  addition to the actual scan and activation, it sees to rebuilding the basis
  and solution. The heavy lifting is performed in scanPrimConStdAct and
  actBLogPrimCon.

  See the comments at the head of the file for the algebra and (lack of) side
  effects associated with constraint activation without variable activation.
  Fortunately, if we're activating variables, too, the details are handled by
  the variable activation routines.

  Parameter:
    orig_sys:	The original constraint system
    with_vars:	If true, consider activating inactive variables associated
		with newly activated constraints.

  Returns: number of constraints activated; -1 if there's an error.
*/

{ int *candidates,cand_cnt,*inactvars,inact_cnt ;
  int retval,var_retval ;
  bool actresult ;
  flags calcflgs ;
  dyret_enum factorresult ;

  const char *rtnnme = "dy_activateCons" ;

  retval = -1 ;

/*
  Call scanPrimConStdAct to return a list of candidates for activation, then
  call actBLogPrimConList to install them. Installing nothing always succeeds.
  If with_vars == TRUE, pass in a pointer that will return loaded with a vector
  of indices of inactive variables referenced by the candidate constraints.
*/
  candidates = NULL ;
  inactvars = NULL ;
  inact_cnt = 0 ;
  cand_cnt = scanPrimConStdAct(orig_sys,&candidates) ;
  if (cand_cnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint","normal activation") ;
    actresult = FALSE ; }
  else
  if (cand_cnt > 0)
  { if (with_vars == TRUE)
    { actresult = dy_actBLogPrimConList(orig_sys,cand_cnt,
					candidates,&inactvars) ;
      if (inactvars != NULL)
      { inact_cnt = inactvars[0] ; } }
    else
    { actresult = dy_actBLogPrimConList(orig_sys,cand_cnt,
					candidates,NULL) ; } }
  else
  { actresult = TRUE ; }
  if (candidates != NULL) FREE(candidates) ;
  if (actresult == FALSE)
  { if (inactvars != NULL) FREE(inactvars) ;
    return (retval) ; }
/*
  If we added constraints, we need to refactor and recalculate the primal and
  dual variables. It's unlikely we'll have primal feasibility, but very likely
  we'll have dual feasibility.
*/
  retval = cand_cnt ;
  if (cand_cnt > 0)
  { dy_lp->simplex.init_dse = TRUE ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.conmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      factoring, calculating variables, ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"and checking feasibility ...") ; }
#   endif
    calcflgs = ladFACTOR|ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
    factorresult = dy_accchk(&calcflgs) ;
    switch (factorresult)
    { case dyrOK:
      case dyrPATCHED:
      { retval = cand_cnt ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	{ if (factorresult == dyrOK)
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ;
	  else
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    patched.") ; }
#       endif
	if (flgoff(calcflgs,ladDUALFEAS))
	{ dy_lp->simplex.next = dyDUAL ; }
	else
	if (flgoff(calcflgs,ladPRIMFEAS))
	{ dy_lp->simplex.next = dyPRIMAL2 ; }
	else
	{ dy_lp->simplex.next = dyPRIMAL1 ; }
	break ; }
      default:
      { retval = -1 ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#       endif
	break ; } }
/*
  Are we activating referenced variables? If so, invoke the variable management
  routines.
*/
    if (with_vars == TRUE && inact_cnt > 0)
    { var_retval = dy_activateVars(orig_sys,inactvars) ;
      if (var_retval < 0)
      { errmsg(440,rtnnme,dy_sys->nme,
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
	retval = var_retval ; } } }
  if (inactvars != NULL) FREE(inactvars) ;


# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  activated %d constraints",cand_cnt) ;
    if (with_vars == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  " with %d referenced variables",inact_cnt) ; }
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  constraint system %s now %d x %d (%d + %d).",
	        dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
	        dy_sys->logvcnt) ; }
# endif

  return (retval) ; }

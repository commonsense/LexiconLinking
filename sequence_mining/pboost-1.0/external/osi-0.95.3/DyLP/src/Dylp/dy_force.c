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

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_force.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_force.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  This file contains a pair of routines that attempt to force a transition
  from primal to dual simplex or from dual to primal simplex, and a final
  routine which loads the full constraint system. These routines are used as
  part of dylp's error recovery strategy. When one simplex encounters
  pivoting difficulties which cannot be easily resolved within that simplex,
  a possible solution is to force a transition to the other simplex. If this
  can't be achieved in a way that guarantees some progress toward the optimal
  solution, the last ditch fall-back position is to load the entire constraint
  system.

  A specific example will help. Suppose that primal phase II is running, and
  it detects that all desireable entering variables have been rejected
  because the pivot results in a singular basis. In this case, it will return
  a punt.

  One way to resolve the problem would be to find some other entering
  variable. This is easy: it's the standard activate variable step. We can
  resume pivoting without changing simplex algorithms.

  But suppose we can't find any variables to activate? A second way to
  resolve the problem would be to change the basis matrix by adding
  constraints.  There's no sense adding loose constraints, however, and
  adding violated constraints will cost us primal feasibility. But there's a
  way out of this bind: force dual feasibility. Then we can add as many
  violated constraints as we like and resume pivoting in dual simplex.  Note
  that we really do need to find violated constraints. If we're primal
  feasible, and force dual feasibility, we're optimal and done!

  But because we're working with a partial system, it can happen that there
  are no violated constraints. We need to add both constraints and variables
  in order to make progress. (Trust me. It's happened.) One can think of lots
  of clever things to try, but the computational cost to determine an
  intelligent action is generally high. Not so clever, but guaranteed to work
  if anything will, is to load the full constraint system.

  (All the above works in primal phase I, with the observation that infeasible
  means we're `optimal' under the phase I objective but haven't attained
  feasibility.)

  How do we force dual feasibility? By removing nonbasic variables with
  favourable (hence dual infeasible) reduced costs. For architecturals, this
  is trivial.
  
  For logicals, it's a nightmare. Deleting a nonbasic logical implies
  deleting a (tight) constraint which has an associated nonzero (basic) dual.
  In general, this changes the values of all dual variables, so our current
  notion of dual feasibility/infeasibility goes out the window. Further,
  there's a dangling basic variable to be dealt with. If it's feasible, we
  can set it to superbasic status. But if it's infeasible, we can only set it
  to the nearest bound, and that'll change our primal feasibility status.

  Forcing primal feasibility has the same feel. Deactivating a violated
  explicit constraint is trivial. Deactivating a violated bound constraint is
  impossible, because bounds on variables are enforced as part of the primal
  simplex algorithm, and deactivated variables are assumed to be at bound.

  Experience during code development says that if the dual <-> primal
  transition cannot be accomplished by trivial deactivations, it's not worth
  the attempt. If you want to allow dylp to try, set the heroics option to
  TRUE. In the absence of heroics, dylp loads the full constraint system.

  Failure of a transition attempt can occur in two places:
    * During the first half of the transition, in forcePrimal2Dual
      (forceDual2Primal), because heroics are forbidden or because they were
      tried and didn't work.
    * During the second half of the transition, in gen/add constraints
      (gen/add variables), when it is determined that there are no violated
      constraints (variables with favourable reduced costs) to add to the
      active system. I.e., we're optimal in this subproblem and need to add
      both constraints and variables to make progress.
*/

/*
  Structure to hold candidates for deactivation or flipping when attempting to
  force dual feasibility.

  Field		Description
  -----		-----------
  ndx		variable index
  flippable	TRUE if the variable can be flipped to the opposite bound
  delta		delta for a bound to bound flip
*/

typedef struct { int ndx ;
		 bool flippable ;
		 double delta ; } fdcand_struct ;



static int fdcandcompare (const void *p_fdcand1, const void *p_fdcand2)
/*
  Comparison function to sort an array of fdcand_struct into nonincreasing
  order of indices.

  Returns:	< 0	if i > j
		  0	if i = j
		> 0	if i < j
*/
{ int i,j ;
  const fdcand_struct *fdcand1,*fdcand2 ;

  fdcand1 = (const fdcand_struct *) p_fdcand1 ;
  fdcand2 = (const fdcand_struct *) p_fdcand2 ;

  i = fdcand1->ndx ;
  j = fdcand2->ndx ;

  return ((j)-(i)) ; }

static int intcompare (const void *p_i, const void *p_j)
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



static int scanPrimVarDualInfeas (fdcand_struct **p_fdcands)
/*
  This routine scans the active constraint system for variables that are dual
  infeasible and thus must be deactivated or flipped in order to force dual
  feasibility.  Put another way, it looks for nonbasic variables with
  favourable reduced costs: status NBLB and cbar<j> < 0 and status NBUB and
  cbar<j> > 0.

  Parameters:
    p_fdcands:	(i) empty vector to hold candidates; assumed to be
		    sufficiently large; will be allocated if NULL
		(o) information on variables to be flipped or deactivated;
		    fdcands[0].flippable is set TRUE/FALSE to indicate if
		    any variables can be flipped; may be NULL if no candidates
		    are found.

  Returns: number of variables to be deactivated/flipped, -1 if there's an
	   error during scanning (error is possible only when paranoid)
*/

{ int j,n,purgecnt,flipcnt ;
  fdcand_struct *fdcands ;
  double cbarj,lbj,ubj ;
  double *vlb,*vub ;
  flags statj ;
  bool purge ;

# ifdef PARANOIA

  const char *rtnnme = "scanPrimVarDualInfeas" ;

  if (p_fdcands == NULL)
  { errmsg(2,rtnnme,"fdcands") ;
    return (-1) ; }
# endif

/*
  Scan preparation.  If the user hasn't supplied a vector to hold the
  indices, allocate one now.
*/
  n = dy_sys->varcnt ;
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
  purgecnt = 0 ;
  if (*p_fdcands == NULL)
  { fdcands =
      (fdcand_struct *) MALLOC((dy_sys->archvcnt+1)*sizeof(fdcand_struct)) ; }
  else
  { fdcands = *p_fdcands ; }
/*
  Scan the variables to collect the list of variables to be deactivated or
  flipped: NBLB status with cbarj < 0 or NBUB status with cbarj > 0.
*/
  flipcnt = 0 ;
  for (j = 1 ; j <= n ; j++)
  { purge = FALSE ;
    statj = dy_status[j] ;
    cbarj = dy_cbar[j] ;
    if ((flgon(statj,vstatNBLB) && cbarj < -dy_tols->dfeas) ||
        (flgon(statj,vstatNBUB) && cbarj > dy_tols->dfeas))
    { fdcands[++purgecnt].ndx = j ;
      fdcands[purgecnt].flippable = FALSE ;
      lbj = vlb[j] ;
      ubj = vub[j] ;
      if (flgon(statj,vstatNBLB))
      { if (ubj < dy_tols->inf)
	{ fdcands[purgecnt].flippable = TRUE ;
	  flipcnt++ ;
	  fdcands[purgecnt].delta = ubj-lbj ; } }
      else
      { if (lbj > -dy_tols->inf)
	{ fdcands[purgecnt].flippable = TRUE ;
	  flipcnt++ ;
	  fdcands[purgecnt].delta = lbj-ubj ; } }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      queuing %s (%d) for %s, %s, cbar<%d> = %g",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j,
		    (fdcands[purgecnt].flippable == TRUE)?"flip":"deactivation",
		    dy_prtvstat(statj),j,cbarj) ; }
#     endif
    } }
/*
  Prepare for return. If we found candidates, set fdcands[0] and return the
  vector. If we found none, and we allocated the candidate vector, then free
  it.
*/
  if (purgecnt == 0)
  { if (*p_fdcands == NULL)
    { FREE(fdcands) ; } }
  else
  { if (flipcnt > 0)
      fdcands[0].flippable = TRUE ;
    else
      fdcands[0].flippable = FALSE ;
    if (*p_fdcands == NULL)
      *p_fdcands = fdcands ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    queued %d variables, %d deactivate, %d flip.",
	        purgecnt,purgecnt-flipcnt,flipcnt) ; }
# endif

  return (purgecnt) ; }



dyphase_enum dy_forcePrimal2Dual (consys_struct *orig_sys)

/*
  This routine attempts to force dual feasibility.  The approach is to
  deactivate or repair violated dual constraints (i.e., nonbasic primal
  variables with favourable reduced costs).  By `repair', we mean flipping
  the variable to its opposite bound, thus restoring dual feasibility.
  Variables can be flipped only if both bounds are finite.

  Deactivating or flipping nonbasic architectural variables with favourable
  reduced costs is trivial. Flipping nonbasic primal logicals associated with
  range constraints (hence lower and upper bounded) is no different than
  flipping nonbasic architecturals.

  Deactivating nonbasic primal logicals, which are associated with tight
  constraints, which are associated with basic dual architecturals, is a can
  of worms. In general, the associated dual will be nonzero. The primal
  constraint and logical must be deactivated together, which amounts to
  forcing the dual to zero, which has the potential to change the feasibility
  of every dual constraint. Plus, there's a dangling basic variable to be
  dealt with. If we're in primal II, we can make it superbasic (guaranteeing
  we won't achieve dual feasibility, but also not losing primal feasibility).
  If we're in primal I, well, who cares, just force the variable to the
  nearest bound.  The best we can do is deactivate the constraint and, at the
  end, refactor, recalculate primals and duals, check for primal and dual
  feasibility, and go with the flow.

  So, what to do? Experience says that forcing deactivation of nonbasic
  logicals is at best a crap shoot, and often wasted effort. It's
  controllable. If the heroics.p2d option is set to TRUE, tight constraints
  will be deactivated. Experience also says that flipping nonbasic variables
  tends to result in cycling. Flipping is also controllable (heroics.flip);
  if disallowed, a flippable variable will be deactivated instead.

  What's the result? If we're converting from primal to dual and deactivated
  only nonbasic architecturals, then we'll have dual feasibility, hence
  primal optimality, at the end of the routine. In this case, we can reset
  the lp result to dyOPTIMAL and head for dyGENCON, the normal route to dual
  simplex.

  If we suppressed deactivation of tight constraints, then we won't have dual
  feasibility.  Moreover, we haven't dealt with all the problem pivot
  candidates, so there's no sense simply returning to primal simplex --- we
  need to do something to change the situation. No sense wasting time on half
  measures. We're in trouble, so force the full system.

  If we deactivated tight constraints, we've certainly changed the basis, and
  anything is possible. We'll need to refactor and do feasibility checks, and
  go with the flow.

  This routine is also used to attempt to restore dual feasibility while
  running dual simplex. In this case, if we're successful we'll head back to
  dual simplex, detouring via dyADDVARS in an attempt to add some dual
  constraints to balance the ones we've just deactivated. If we've failed to
  restore dual feasibility, we can't return to dual simplex. We already have
  primal infeasibility, so we head for dyPRIMAL1, again detouring via ADDVARS,
  this time with the objective of adding some variables to help us gain
  feasibility.

  The result of processing is communicated in dy_lp->simplex.next, which will
  be set to dyDUAL or dyPRIMAL[1|2] as appropriate, and in the returned phase
  code.

  Parameters:
    orig_sys:	The original constraint system

  Returns: appropriate next phase (see end of routine), or dyINV if an
	   error occurs.
*/

{ int j,m ;
  int cand_cnt,ndx,varcnt,concnt,suppressed ;
  fdcand_struct *candidates ;
  double *vlb,*vub ;
  dyret_enum factorresult ;
  flags calcflgs,statj ;
  bool retval ;
  dyret_enum upd_retval ;
  dyphase_enum next_phase ;

  const char *rtnnme = "dy_forcePrimal2Dual" ;

  next_phase = dyINV ;

# ifdef PARANOIA
  retval = FALSE ;
  if (dy_lp->simplex.active == dyDUAL && dy_lp->lpret == lpLOSTFEAS)
    retval = TRUE ;
  if ((dy_lp->simplex.active == dyPRIMAL1 ||
       dy_lp->simplex.active == dyPRIMAL2) &&
      (dy_lp->lpret == lpPUNT || dy_lp->lpret == lpSTALLED))
    retval = TRUE ;
  if (retval == FALSE)
  { errmsg(441,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   dy_prtlpphase(dy_lp->simplex.active,TRUE),
	   dy_prtlpret(dy_lp->lpret)) ;
    return (dyINV) ; }
# endif

/*
 If the primal phase I objective is installed, boot it out and reinstall the
 phase II objective.
*/
  if (dy_lp->p1obj.installed == TRUE)
  { if (dy_swapobjs(dyPRIMAL2) == FALSE)
    { errmsg(318,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "remove") ;
      return (dyINV) ; }
    dy_calcduals() ;
    if (dy_calccbar() == FALSE)
    { errmsg(384,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters) ;
      return (dyINV) ; } }
/*
  Call scanPrimVarDualInfeas to return a list of candidates for deactivation.
*/
  candidates = NULL ;
  cand_cnt = scanPrimVarDualInfeas(&candidates) ;
  if (cand_cnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable","forced primal -> dual transition") ;
    if (candidates != NULL) FREE(candidates) ;
    return (dyINV) ; }
# ifdef PARANOIA
/*
  We should always find candidates, otherwise why are we here?
*/
  if (cand_cnt == 0)
  { errmsg(1,rtnnme,__LINE__) ;
    if (candidates != NULL) FREE(candidates) ;
    return (dyINV) ; }
# endif
/*
  We have candidates. Sort the list and then open a loop to do the
  deactivations. In the case of logicals, that entails deactivating the
  constraint.
*/
  qsort(&candidates[1],cand_cnt,sizeof(fdcand_struct),fdcandcompare) ;
  m = dy_sys->concnt ;
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
  varcnt = 0 ;
  concnt = 0 ;
  suppressed = 0 ;
  retval = TRUE ;
  for (ndx = 1 ; ndx <= cand_cnt && retval == TRUE ; ndx++)
  { j = candidates[ndx].ndx ;
#   ifdef PARANOIA
    if (j < 1 || j > dy_sys->varcnt)
    { errmsg(102,rtnnme,dy_sys->nme,"variable",j,1,dy_sys->varcnt) ;
      retval = FALSE ;
      break ; }
#   endif
/*
  The easiest case: we can just flip the variable. dy_updateprimals can return
  dyrOK, dyrREQCHK, dyrSWING, or dyrFATAL. Only dyrFATAL is of concern, as
  we'll be refactoring, etc., on our way back from forcing dual feasibility.
*/
    if (candidates[ndx].flippable == TRUE && dy_opts->heroics.flips == TRUE)
    { varcnt++ ; 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    flipping variable %s (%d)",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j) ; }
#     endif
      upd_retval = dy_updateprimals(j,candidates[ndx].delta,NULL) ;
      if (upd_retval == dyrFATAL)
      { retval = FALSE ;
	break ; }
      statj = dy_status[j] ;
      if (flgon(statj,vstatNBLB))
      { dy_x[j] = vub[j] ; }
      else
      { dy_x[j] = vlb[j] ; }
      comflg(dy_status[j],vstatNBLB|vstatNBUB) ; }
/*
  The easy case: we have a nonbasic architectural. Just call
  dy_deactNBPrimArch to deal with it.
*/
    else
    if (j > m)
    { varcnt++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    deactivating variable %s (%d)",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j) ; }
#     endif
      retval = dy_deactNBPrimArch(orig_sys,j) ;
      if (retval == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "deactivate","variable",
	       consys_nme(dy_sys,'v',j,TRUE,NULL),j) ; } }
/*
  The hard case: we have a nonbasic logical. Not that it looks any more
  difficult from here. If we're primal feasible, this call will create
  superbasic variables when it forces the dangling basic variable into the
  nonbasic partition.
*/
    else
    if (dy_opts->heroics.p2d == TRUE)
    { concnt++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    deactivating constraint %s (%d)",
		    consys_nme(dy_sys,'c',j,TRUE,NULL),j) ; }
#     endif
      retval = dy_deactNBLogPrimCon(orig_sys,j) ;
      if (retval == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "deactivate","constraint",
	       consys_nme(dy_sys,'c',j,TRUE,NULL),j) ; } }
    else
    { suppressed++ ; } }

  FREE(candidates) ;
  if (retval == FALSE) return (dyINV) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 1)
  { if (dy_opts->print.force >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho," %d+%d deletions.",concnt,varcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  constraint system %s now %d x %d (%d + %d).",
		dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
		dy_sys->logvcnt) ; }
# endif

# ifdef PARANOIA
  if (dy_chkdysys(orig_sys) == FALSE) return (dyINV) ;
# endif

/*
  Time to clean up a bit. Keeping in mind that we got here because we were in
  trouble (lpPUNT or lpLOSTFEAS), it isn't going to hurt to refactor and do
  accuracy and feasibility checks whether we need it or not.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n      factoring, checking accuracy and feasibility ...") ; }
# endif
  calcflgs = ladFACTOR|ladPRIMALCHK|ladDUALCHK|
	     ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
  factorresult = dy_accchk(&calcflgs) ;
  switch (factorresult)
  { case dyrOK:
    case dyrPATCHED:
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { if (factorresult == dyrOK)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    patched.") ;
	dyio_outfmt(dy_logchn,dy_gtxecho," Feasibility:") ;
	if (flgoff(calcflgs,ladPRIMFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," primal") ; }
	if (flgoff(calcflgs,ladDUALFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," dual") ; }
	if (flgall(calcflgs,ladPRIMFEAS|ladDUALFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," none") ; } }
#     endif
      break ; }
    default:
    { next_phase = dyINV ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#     endif
      break ; } }

# ifdef PARANOIA
/*
  The crucial question is whether we have dual feasibility. See the long
  explanation at the head of the routine for the detailed reasons behind the
  phase transitions and associated actions. Here on the front line, the truth
  of the matter is that we've refactored and numerical wierdness could put us
  anywhere. If we're feeling paranoid, check that feasibility is what we
  expect.
*/
  if (suppressed == 0)
  { if (concnt == 0 && flgon(calcflgs,ladDUALFEAS))
    { warn(439,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "absence","dual") ; } }
  else
  { if (flgoff(calcflgs,ladDUALFEAS))
    { warn(439,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "presence","dual") ; } }
  if (dy_lp->simplex.active == dyDUAL && flgoff(calcflgs,ladPRIMFEAS))
  { warn(439,rtnnme,
	 dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	 "presence","primal") ; }
# endif

/*
  We're here in an attempt to restore dual feasibility. If we've changed the
  basis, we need to reinitialise the DSE norms. I can't conceive of why we'd
  suddenly acquire primal feasibility, but it won't hurt to code for it.
  Leave the lp return code as lpLOSTFEAS. The expectation is that dyADDVAR will
  determine if the correct primal objective is installed and take action if
  needed.
*/
  if (dy_lp->simplex.active == dyDUAL)
  { if (flgoff(calcflgs,ladDUALFEAS))
    { dy_lp->simplex.next = dyDUAL ;
      if (concnt != 0) dy_lp->simplex.init_dse = TRUE ;
      next_phase = dyADDVAR ; }
    else
    { dy_lp->simplex.init_pse = TRUE ;
      if (flgoff(calcflgs,ladPRIMFEAS))
      { dy_lp->simplex.next = dyPRIMAL2 ; }
      else
      { dy_lp->simplex.next = dyPRIMAL1 ; }
      next_phase = dyADDVAR ; } }
/*
  We're attempting a primal -> dual transition. Dual feasibility says we were
  successful in forcing dual feasibility, the first step of the transition.
  Otherwise, we're in deep trouble; might as well force the full system.
*/
  else
  { dy_lp->lpret = lpFORCEDUAL ;
    if (flgoff(calcflgs,ladDUALFEAS))
    { dy_lp->simplex.next = dyDUAL ;
      dy_lp->simplex.init_dse = TRUE ;
      next_phase = dyGENCON ; }
    else
    { next_phase = dyFORCEFULL ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n   next phase %s, next simplex %s.",
	        dy_prtlpphase(next_phase,FALSE),
	        dy_prtlpphase(dy_lp->simplex.next,FALSE)) ; }
# endif

  return (next_phase) ; }




static int scanPrimConForceDeact (int **p_acndxs)
/*
  This routine scans the active constraint system for constraints (including
  implicit bound constraints) that are violated and thus must be deactivated
  in order to force primal feasibility.  Put another way, it looks for basic
  variables that are outside of bounds.

  Parameters:
    p_acndxs:	(i) empty vector to hold constraint indices; assumed to be
		    sufficiently large; will be allocated if NULL
		(o) indices of constraint to be deactivated; may not be
		    allocated if no candidates are identified

  Returns: number of constraint to be deactivated, -1 if there's an error
	   during scanning (error is possible only when paranoid)
*/

{ int bpos,j,m,purgecnt ;
  int *acndxs ;
  flags statj ;


# ifdef PARANOIA

  const char *rtnnme = "scanPrimConForceDeact" ;

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
  straightforward, as all we need to do is examine the status of the
  basic variable.
*/
  purgecnt = 0 ;
  for (bpos = 1 ; bpos <= m ; bpos++)
  { j = dy_basis[bpos] ;
    statj = dy_status[j] ;
    if (flgon(statj,vstatBLLB|vstatBUUB))
    { acndxs[purgecnt++] = j ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 3)
      { if (j <= m)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n    queued %s %s (%d) for deactivation, ",
		      consys_prtcontyp(dy_sys->ctyp[j]),
		      consys_nme(dy_sys,'c',j,TRUE,NULL),j) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "%s (%d) = %g, status %s, basis pos'n %d.",
		      consys_nme(dy_sys,'v',j,TRUE,NULL),j,
		      dy_x[j],dy_prtvstat(statj),bpos) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n    queued %s (%d) = %g for deactivation, ",
		      consys_nme(dy_sys,'v',j,TRUE,NULL),j,dy_x[j]) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"status %s, basis pos'n %d.",
		      dy_prtvstat(statj),bpos) ; } }
#     endif

    } }

  if (*p_acndxs == NULL)
  { if (purgecnt <= 0)
    { FREE(acndxs) ; }
    else
    { *p_acndxs = acndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    queued %d constraints for deactivation.",purgecnt) ; }
# endif

  return (purgecnt) ; }



dyphase_enum dy_forceDual2Primal (consys_struct *orig_sys)

/*
  This routine attempts to force primal feasibility. The approach is to
  deactivate violated primal constraints with basic logicals (i.e., nonbasic
  dual architectural variables with favourable reduced costs).

  There's one major fly in the ointment: basic primal architecturals which
  violate one of their (implicit) bound constraints. If this bound were an
  explicit constraint, we'd have a basic logical and life would be simple.
  But it isn't, and we don't, and life is complicated. We don't want to shoot
  the messenger (the explicit constraint in this basis pos'n). The real
  problem is that the implicit bound constraint is wired into the simplex
  algorithm and it's not possible to deactivate it. Pushing the variable
  into the nonbasic partition, or going further and deactivating it, simply
  enforces the bound constraint. Still, there's little else we can do if we
  want to deactivate them.

  So, what to do? Experience says that forcing deactivation of basic
  architecturals is wasted effort, but it's controllable. If the heroics.d2p
  option is set to TRUE, basic architecturals are forced out.

  What's the result? If the only violated constraints are explicit
  constraints, we'll have primal feasibility at the end of the routine. If
  there are violated implicit bound constraints, and heroics is FALSE, we
  won't have primal feasibility, and there's no point in returning to dual
  simplex because we haven't dealt with all the problem pivot candidates. If
  there are violated implicit bound constraints and heroics is TRUE, anything
  is possible and we just have to check.

  The result of processing is communicated in dy_lp->simplex.next, which will
  be set to dyDUAL or dyPRIMAL[1|2] as appropriate, and the returned phase
  code.

  Parameters:
    orig_sys:	The original constraint system

  Returns: appropriate next phase, or dyINV if an error occurs
*/

{ int j,m ;
  int *candidates,cand_cnt,ndx,varcnt,concnt,suppressed ;
  dyret_enum factorresult ;
  flags calcflgs ;
  bool retval ;
  dyphase_enum next_phase ;

  const char *rtnnme = "dy_forceDual2Primal" ;

# ifdef PARANOIA
  retval = FALSE ;
  if (dy_lp->simplex.active == dyDUAL &&
      (dy_lp->lpret == lpPUNT || dy_lp->lpret == lpSTALLED))
    retval = TRUE ;
  if (retval == FALSE)
  { errmsg(441,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   dy_prtlpphase(dy_lp->simplex.active,TRUE),
	   dy_prtlpret(dy_lp->lpret)) ; }
# endif

  next_phase = dyINV ;

/*
  Have we made any progress since the last time we were here? If not, we're
  quite possibly looping through constraint deactivation/activation. Just
  head for primal phase I.
*/
  if (dy_lp->z-dy_lp->lastz.fp < dy_tols->purge*(1.0+fabs(dy_lp->z)))
  { dy_lp->simplex.next = dyPRIMAL1 ;
    dy_lp->simplex.init_pse = TRUE ;
    return (dyPRIMAL1) ; }
  dy_lp->lastz.fp = dy_lp->z ;
/*
  Call scanPrimConForceDeact to return a list of candidates for deactivation.
*/
  candidates = NULL ;
  cand_cnt = scanPrimConForceDeact(&candidates) ;
  if (cand_cnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint","forced dual -> primal transition") ;
    if (candidates != NULL) FREE(candidates) ;
    return (dyINV) ; }
# ifdef PARANOIA
/*
  We should always find candidates, otherwise why are we here?
*/
  if (cand_cnt == 0)
  { errmsg(1,rtnnme,__LINE__) ;
    if (candidates != NULL) FREE(candidates) ;
    return (dyINV) ; }
# endif
/*
  We have candidates. Sort the list and then open a loop to do the
  deactivations. In the case of explicit constraints with basic logicals,
  it's easy: we just deactivate the constraint. If the violated constraint is
  an implicit bound, we proceed as discussed in the comments at the head of
  the routine.
*/
  qsort(&candidates[0],cand_cnt,sizeof(int),intcompare) ;
  m = dy_sys->concnt ;
  varcnt = 0 ;
  concnt = 0 ;
  suppressed = 0 ;
  retval = TRUE ;
  for (ndx = 0 ; ndx < cand_cnt && retval == TRUE ; ndx++)
  { j = candidates[ndx] ;
#   ifdef PARANOIA
    if (j < 1 || j > dy_sys->varcnt)
    { errmsg(102,rtnnme,dy_sys->nme,"variable",j,1,dy_sys->varcnt) ;
      retval = FALSE ;
      break ; }
#   endif
/*
  The easy case: we have a basic logical. dy_deactBLogPrimCon will deal with
  this neatly.
*/
    if (j <= m)
    { concnt++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    deactivating constraint %s (%d)",
		    consys_nme(dy_sys,'c',j,TRUE,NULL),j) ; }
#     endif
      retval = dy_deactBLogPrimCon(orig_sys,j) ;
      if (retval == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "deactivate","constraint",
	       consys_nme(dy_sys,'c',j,TRUE,NULL),j) ; } }
/*
  The hard case: we have a basic architectural. The ugliness is hidden in
  dy_deactBPrimArch.
*/
    else
    if (dy_opts->heroics.d2p == TRUE)
    { varcnt++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    deactivating variable %s (%d)",
		    consys_nme(dy_sys,'v',j,TRUE,NULL),j) ; }
#     endif
      retval = dy_deactBPrimArch(orig_sys,j) ;
      if (retval == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "deactivate","variable",
	       consys_nme(dy_sys,'v',j,TRUE,NULL),j) ; } }
    else
    { suppressed++ ; } }

  FREE(candidates) ;
  if (retval == FALSE) return (dyINV) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 1)
  { if (dy_opts->print.conmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho," %d+%d deletions.",concnt,varcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  constraint system %s now %d x %d (%d + %d).",
		dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
		dy_sys->logvcnt) ; }
# endif

# ifdef PARANOIA
  if (dy_chkdysys(orig_sys) == FALSE) return (dyINV) ;
# endif

/*
  Time to clean up a bit. Whatever we've done, it's changed the basis, so we
  might as well refactor now. While we're there, might as well do accuracy
  and feasibility checks.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n      factoring, checking accuracy and feasibility.") ; }
# endif
  calcflgs = ladFACTOR|ladPRIMALCHK|ladDUALCHK|
	     ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
  factorresult = dy_accchk(&calcflgs) ;
  switch (factorresult)
  { case dyrOK:
    case dyrPATCHED:
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { if (factorresult == dyrOK)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    patched.") ;
	dyio_outfmt(dy_logchn,dy_gtxecho," Feasibility:") ;
	if (flgoff(calcflgs,ladPRIMFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," primal") ; }
	if (flgoff(calcflgs,ladDUALFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," dual") ; }
	if (flgall(calcflgs,ladPRIMFEAS|ladDUALFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," none") ; } }
#     endif
      break ; }
    default:
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#     endif
      return (dyINV) ; } }
/*
  We've successfully refactored & calculated primal and dual feasibility.
  The crucial question is whether we have primal feasibility.
  
  If we're primal feasible, we've completed the first half of the dual ->
  primal transition.  We can set simplex.next = dyPRIMAL2 and head for
  dyADDVAR to look for variables with favourable reduced cost.

  If we failed to gain primal feasibility, well, specify primal phase I
  instead and see how it goes.

  If we're paranoid, we check that the feasibility is what we expect.
*/
# ifdef PARANOIA
  if (suppressed == 0)
  { if (varcnt == 0 && flgon(calcflgs,ladPRIMFEAS))
    { warn(439,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "absence","primal") ; } }
  else
  { if (flgoff(calcflgs,ladPRIMFEAS))
    { warn(439,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "presence","primal") ; } }
# endif
  
  dy_lp->lpret = lpFORCEPRIMAL ;
  if (flgoff(calcflgs,ladPRIMFEAS))
  { dy_lp->simplex.next = dyPRIMAL2 ;
    next_phase = dyADDVAR ; }
  else
  { dy_lp->simplex.next = dyPRIMAL1 ;
    next_phase = dyADDVAR ; } 
  dy_lp->simplex.init_pse = TRUE ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n   next phase %s, next simplex %s.",
	        dy_prtlpphase(next_phase,FALSE),
	        dy_prtlpphase(dy_lp->simplex.next,FALSE)) ; }
# endif

  return (next_phase) ; }




static int scanPrimVarForceAct (consys_struct *orig_sys, int **p_ovndxs)

/*
  This routine scans the original constraint system looking for inactive
  variables to add to the active system.

  Parameters:
    orig_sys:	The original constraint system
    p_ovndxs:	(i) empty vector to hold constraint indices; assumed
		    sufficiently large if non-NULL; if NULL, allocated if
		    necessary
		(o) indices of constraints to be activated; may not be
		    allocated if no constraints are identified

  Returns: number of candidates for activation, -1 if error.
*/

{ int j,n,actcnt ;
  flags statj ;
  int *ovndxs ;


# ifdef PARANOIA

  const char *rtnnme = "scanPrimVarForceAct" ;

  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (-1) ;  }
  if (p_ovndxs == NULL)
  { errmsg(2,rtnnme,"&ovndxs") ;
    return (-1) ; }
# endif

  n = orig_sys->varcnt ;

/*
  Did the client supply a vector for candidate indices? If not, make one.
*/
  actcnt = n-dy_sys->archvcnt ;
  if (*p_ovndxs == NULL)
  { ovndxs = (int *) MALLOC(actcnt*sizeof(int)) ; }
  else
  { ovndxs = *p_ovndxs ; }
/*
  Now step through the variables, remembering the inactive ones. As always, we
  never activate fixed variables.
*/
  actcnt = 0 ;
  for (j = 1 ; j <= n ; j++)
  { if (dy_origvars[j] > 0) continue ;
    statj = (flags) (-dy_origvars[j]) ;
    if (flgon(statj,vstatNBFX)) continue ;
    ovndxs[actcnt++] = j ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.force >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    queued %s variable %s (%d),",
		  consys_prtvartyp(orig_sys->vtyp[j]),
		  consys_nme(orig_sys,'v',j,FALSE,NULL),j) ; }
#   endif
  }
/*
  If we supplied ovndxs and found no candidates to activate, free it.
*/
  if (*p_ovndxs == NULL)
  { if (actcnt == 0)
    { FREE(ovndxs) ; }
    else
    { *p_ovndxs = ovndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d variables for activation.",actcnt) ; }
# endif

  return (actcnt) ; }


static int scanPrimConForceAct (consys_struct *orig_sys, int **p_ocndxs)

/*
  This routine scans the original constraint system looking for inactive
  constraints to add to the active system.

  Parameters:
    orig_sys:	The original constraint system
    p_ocndxs:	(i) empty vector to hold constraint indices; assumed
		    sufficiently large if non-NULL; if NULL, allocated if
		    necessary
		(o) indices of constraints to be activated; may not be
		    allocated if no constraints are identified

  Returns: number of candidates for activation, -1 if error.
*/

{ int i,m,actcnt ;
  int *ocndxs ;


# ifdef PARANOIA

  const char *rtnnme = "scanPrimConForceAct" ;

  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (-1) ;  }
  if (p_ocndxs == NULL)
  { errmsg(2,rtnnme,"&ocndxs") ;
    return (-1) ; }
# endif

  m = orig_sys->concnt ;

/*
  Did the client supply a vector for candidate indices? If not, make one.
*/
  actcnt = m-dy_sys->concnt ;
  if (*p_ocndxs == NULL)
  { ocndxs = (int *) MALLOC(actcnt*sizeof(int)) ; }
  else
  { ocndxs = *p_ocndxs ; }
/*
  Now step through the constraints, remembering the inactive ones.
*/
  actcnt = 0 ;
  for (i = 1 ; i <= m ; i++)
  { if (dy_origcons[i] > 0) continue ;
    ocndxs[actcnt++] = i ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.force >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    queued %s constraint %s (%d),",
		  consys_prtcontyp(orig_sys->ctyp[i]),
		  consys_nme(orig_sys,'c',i,FALSE,NULL),i) ; }
#   endif
     }
/*
  If we supplied ocndxs and found no candidates to activate, free it.
*/
  if (*p_ocndxs == NULL)
  { if (actcnt == 0)
    { FREE(ocndxs) ; }
    else
    { *p_ocndxs = ocndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d constraints for activation.",actcnt) ; }
# endif

  return (actcnt) ; }

dyphase_enum dy_forceFull (consys_struct *orig_sys)

/*
  This routine activates all inactive variables and constraints. It's the
  last resort for error recovery by constraint sytem modification.

  Parameters:
    orig_sys:	The original constraint system

  Returns: appropriate next phase (see end of routine), or dyINV if an
	   error occurs.
*/

{ int *candidates,varcnt,concnt ;
  dyret_enum factorresult ;
  flags calcflgs ;
  bool retval ;
  dyphase_enum next_phase ;

  const char *rtnnme = "dy_forceFull" ;

  next_phase = dyINV ;

# ifdef PARANOIA
  if (!(dy_lp->lpret == lpFORCEDUAL || dy_lp->lpret == lpFORCEPRIMAL ||
	dy_lp->lpret == lpPUNT || dy_lp->lpret == lpACCCHK))
  { errmsg(4,rtnnme,"simplex return code",dy_prtlpret(dy_lp->lpret)) ;
    return (dyINV) ; }
# endif

/*
  Call scanPrimConForceAct to return a list of candidates for activation,
  then call dy_actBLogPrimConList to activate them.
*/
  candidates = NULL ;
  concnt = scanPrimConForceAct(orig_sys,&candidates) ;
  if (concnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint","forced full activation") ;
    retval = FALSE ; }
  else
  if (concnt > 0)
  { retval =  dy_actBLogPrimConList(orig_sys,concnt,candidates,NULL) ; }
  else
  { retval = TRUE ; }
  if (candidates != NULL) FREE(candidates) ;
  if (retval == FALSE)
  { return (dyINV) ; }
/*
  And repeat for inactive variables.
*/
  candidates = NULL ;
  varcnt = scanPrimVarForceAct(orig_sys,&candidates) ;
  if (varcnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "variable","forced full activation") ; }
  else
  if (varcnt > 0)
  { retval =  dy_actNBPrimArchList(orig_sys,varcnt,candidates) ; }
  else
  { retval = TRUE ; }
  if (candidates != NULL) FREE(candidates) ;
  if (concnt < 0 || retval == FALSE)
  { return (dyINV) ; } 

# ifdef PARANOIA
/*
  We should have activated at least one constraint or variable. Check the
  constraint system while we're here.
*/
  if (concnt+varcnt == 0)
  { errmsg(1,rtnnme,__LINE__) ;
    return (dyINV) ; }
  if (dy_chkdysys(orig_sys) == FALSE) return (dyINV) ;
# endif

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  %d+%d activations.",concnt,varcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  constraint system %s now %d x %d (%d + %d).",
	        dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
	        dy_sys->logvcnt) ; }
# endif

/*
  Time to clean up a bit.  Refactor and do accuracy and feasibility checks
  whether we need it or not.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n      factoring, checking accuracy and feasibility ...") ; }
# endif
  calcflgs = ladFACTOR|ladPRIMALCHK|ladDUALCHK|
	     ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
  factorresult = dy_accchk(&calcflgs) ;
  switch (factorresult)
  { case dyrOK:
    case dyrPATCHED:
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
      { if (factorresult == dyrOK)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    patched.") ;
	dyio_outfmt(dy_logchn,dy_gtxecho," Feasibility:") ;
	if (flgoff(calcflgs,ladPRIMFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," primal") ; }
	if (flgoff(calcflgs,ladDUALFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," dual") ; }
	if (flgall(calcflgs,ladPRIMFEAS|ladDUALFEAS))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," none") ; } }
#     endif
      break ; }
    default:
    { next_phase = dyINV ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.force >= 2)
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#     endif
      break ; } }
/*
  Where we go next should depend on our feasibility status. Basically, we
  want to head for the appropriate simplex phase.
*/
  dy_lp->lpret = lpFORCEFULL ;
  if (flgoff(calcflgs,ladPRIMFEAS))
  { dy_lp->simplex.next = dyPRIMAL2 ;
    dy_lp->simplex.init_pse = TRUE ;
    next_phase = dyPRIMAL2 ; }
  else
  if (flgoff(calcflgs,ladDUALFEAS) && dy_opts->usedual == TRUE)
  { dy_lp->simplex.next = dyDUAL ;
    dy_lp->simplex.init_dse = TRUE ;
    next_phase = dyDUAL ; }
  else
  { dy_lp->simplex.next = dyPRIMAL1 ;
    dy_lp->simplex.init_pse = TRUE ;
    next_phase = dyPRIMAL1 ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.force >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n   next phase %s, next simplex %s.",
	        dy_prtlpphase(next_phase,FALSE),
	        dy_prtlpphase(dy_lp->simplex.next,FALSE)) ; }
# endif

  return (next_phase) ; }


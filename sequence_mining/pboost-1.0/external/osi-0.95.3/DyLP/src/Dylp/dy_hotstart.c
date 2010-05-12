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
  This file contains the routines that handle a hot start. The assumption is
  that the dylp data structures already exist, and the previous run
  terminated with one of optimal, infeasible, or unbounded status.

  The guiding principle, when thinking about this bit of code, was ``thou
  shalt not change the basis.'' For convenience with this initial
  implementation, I've actually held to ``thou shalt not change the set of
  constraints or variables in the constraint system''.

  This seemed like a good boundary because the user doesn't have access to
  the active system. Making it visible would violate all kinds of modularity
  principles, and writing a bunch of utility routines to make changes didn't
  strike me as fun. It's not possible to allow the user to change inactive
  constraints or variables because this could change the indices of the
  remaining constraints or variables (as they're moved to make room or fill
  holes) and hence indices in dy_actcons and dy_actvars could become
  invalid.

  The limitation that this imposes is that if you want to add or delete a
  constraint or variable, you have to drop back to warm start.

  What is allowed is changes to the variable upper and lower bounds, and the
  right hand side and objective function coefficients. The user must set
  flags indicating which arrays have been modified.

  Changes to the rhs and bounds interact, because the rhs values in the
  active system (dy_sys) must be corrected to reflect the values of inactive
  variables. A change to either will trigger a reinstallation of the rhs
  array followed by recalculation of the correction required for inactive
  variables.  (This seems like overkill until you consider that it's no more
  work to do the calculation from scratch than it is to accumulate
  corrections to the existing values.) A change to the rhs or bounds implies
  that we need to reset the status array and recalculate the primal
  variables.

  Changes to the objective and bounds interact, through the contribution of
  inactive variables to the objective.  A change to the objective also
  implies that we need to recalculate the dual variables.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_hotstart.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_hotstart.c 94 2006-06-29 23:06:51Z lou $" ;



/*
  This routine is exported only to dy_coldstart and dy_warmstart, as well
  as being used in dy_hotstart.
*/

void dy_setfinalstatus (void)

/*
  This code is common to all three start routines (coldstart, warmstart,
  hotstart). It scans the newly calculated basic variables and assigns
  them their final status. In the process, it calculates the number of
  infeasible variables, and the total infeasibility.

  Parameters: none

  Returns: undefined
*/

{ int aindx, xkndx ;
  double xk,lbk,ubk ;

# ifdef PARANOIA
  const char *rtnnme = "dy_setfinalstatus" ;
# endif

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\testablishing final status ...") ;
# endif

  dy_lp->infeas = 0.0 ;
  dy_lp->infeascnt = 0 ;

/*
  Step through the constraints, and have a look at the basic variable in each
  position.

  The paranoid check will complain if the basis is corrupt, but since nothing
  can go wrong if we're not paranoid, it just complains and moves to the next
  entry.
*/
  for (aindx = 1 ; aindx <= dy_sys->concnt ; aindx++)
  { xkndx = dy_basis[aindx] ;
    xk = dy_xbasic[aindx] ;
    lbk = dy_sys->vlb[xkndx] ;
    ubk = dy_sys->vub[xkndx] ;
#   ifdef PARANOIA
    if (xkndx <= 0 || xkndx > dy_sys->varcnt)
    { errmsg(303,rtnnme,dy_sys->nme,aindx,1,xkndx,dy_sys->varcnt) ;
      continue ; }
#   endif
    switch (dy_status[xkndx])
    { case vstatB:
      { if (atbnd(xk,lbk))
	{ dy_status[xkndx] = vstatBLB ; }
	else
	if (belowbnd(xk,lbk))
	{ dy_lp->infeascnt++ ;
	  dy_lp->infeas += lbk-xk ;
	  dy_status[xkndx] = vstatBLLB ; }
	else
	if (atbnd(xk,ubk))
	{ dy_status[xkndx] = vstatBUB ; }
	else
	if (abovebnd(xk,ubk))
	{ dy_lp->infeascnt++ ;
	  dy_lp->infeas += xk-ubk ;
	  dy_status[xkndx] = vstatBUUB ; }
	break ; }
      case vstatBFX:
      { if (!atbnd(xk,lbk))
	{ if (belowbnd(xk,lbk))
	  { dy_lp->infeascnt++ ; 
	    dy_lp->infeas += lbk-xk ;
	    dy_status[xkndx] = vstatBLLB ; }
	  else
	  { dy_lp->infeascnt++ ;
	    dy_lp->infeas += xk-ubk ;
	    dy_status[xkndx] = vstatBUUB ; } }
	break ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.crash >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  %s (%d) %s",
		  consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		  dy_prtvstat(dy_status[xkndx])) ;
      if (lbk > -dy_tols->inf)
	dyio_outfmt(dy_logchn,dy_gtxecho,", lb = %g",lbk) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,", val = %g",xk) ;
      if (ubk < dy_tols->inf)
	dyio_outfmt(dy_logchn,dy_gtxecho,", ub = %g",ubk) ;
      if (flgon(dy_status[xkndx],vstatBLLB|vstatBUUB))
      { dyio_outfmt(dy_logchn,dy_gtxecho,", infeasibility = ") ;
	if (flgon(dy_status[xkndx],vstatBLLB))
	  dyio_outfmt(dy_logchn,dy_gtxecho,"%g",lbk-xk) ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,"%g",xk-ubk) ; }
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
  }
  setcleanzero(dy_lp->infeas,dy_tols->zero) ;

  return ; }



static void hot_updateMiscState (lpret_enum lpret)
/*
  Reset various bits of miscellaneous state before we resume simplex pivoting.

  Pulled out of dy_hotstart so that I can call this in two places.

  Parameters: none

  Returns: undefined
*/

{ dy_lp->lpret = lpret ;
  dy_lp->tot.iters = 0 ;
  dy_lp->tot.pivs = 0 ;
  dy_lp->prev_pivok = dy_lp->pivok ;
  dy_lp->pivok = FALSE ;
  dy_lp->degenpivcnt = 0 ;
  dy_lp->idlecnt = 0 ;

  return ; }


static bool process_inactive (lpprob_struct *orig_lp, int oxkndx)

/*
  This routine handles the data structure updates for an inactive variable
  x<k>.  We need to have a look at the bounds l<k> and u<k>, and perhaps
  update the status kept in dy_origvars. We need to add the contribution
  c<k>l<k> or c<k>u<k> to the objective function. Finally, if we've reloaded
  b & blow due to a bound or rhs change, we need to walk the column a<k>
  and adjust b<i> (and perhaps blow<i>) for each nonzero a<ik> in the active
  system.

  Parameters:
    orig_lp:	the original lp problem
    oxkndx:	index of x<k> in orig_sys
  
  Returns: TRUE if the update is made without incident, FALSE otherwise.
*/

{ int oaindx,aindx,ndx ;
  double xk,lk,uk ;
  pkvec_struct *ak ;
  pkcoeff_struct *aik ;
  consys_struct *orig_sys ;
  flags xkstatus ;
  const char *rtnnme = "process_inactive" ;

  orig_sys = orig_lp->consys ;

# ifdef PARANOIA
/*
  Any inactive variable should be nonbasic, and the paranoid check is looking
  to make sure of this.
*/
  xkstatus = orig_lp->status[oxkndx] ;
  if (!VALID_STATUS(xkstatus))
  { errmsg(300,rtnnme,(int) xkstatus,
	   consys_nme(orig_sys,'v',oxkndx,FALSE,NULL),oxkndx) ;
    return (FALSE) ; }
  if (flgoff(xkstatus,vstatNONBASIC|vstatNBFR))
  { errmsg(433,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "inactive",consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx,
	   dy_prtvstat(xkstatus)) ;
    return (FALSE) ; }
# endif
/*
  The bounds can change arbitrarily, and the client may not be maintaining
  the status vector, but we're limited in what we can do --- our only source
  of a value for inactive variables is the bounds. (Contrast with the
  equivalent section in process_active.) For the typical circumstance where
  both bounds are finite, and a bound has changed, but the variable is not
  fixed or free, the best we can do is use orig_lp->status. If the variable
  was previously fixed, we can only guess; arbitrarily choose NBLB.

  The default case in the switch below should never execute, but it serves
  for paranoia and lets gcc conclude xk will always have a value.
*/
  lk = orig_sys->vlb[oxkndx] ;
  uk = orig_sys->vub[oxkndx] ;
  if (lk > -dy_tols->inf && uk < dy_tols->inf)
  { if (lk == uk)
    { xkstatus = vstatNBFX ; }
    else
    if (orig_lp->status[oxkndx] == vstatNBFX)
    { xkstatus = vstatNBLB ; }
    else
    { xkstatus = orig_lp->status[oxkndx] ; } }
  else
  if (lk > -dy_tols->inf)
  { xkstatus = vstatNBLB ; }
  else
  if (uk < dy_tols->inf)
  { xkstatus = vstatNBUB ; }
  else
  { xkstatus = vstatNBFR ; }
  switch (xkstatus)
  { case vstatNBLB:
    case vstatNBFX:
    { xk = lk ;
      break ; }
    case vstatNBUB:
    { xk = uk ;
      break ; }
    case vstatNBFR:
    { xk = 0 ;
      break ; }
    default:
    { xk = 0 ;
      errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
  orig_lp->status[oxkndx] = xkstatus ;
  dy_origvars[oxkndx] = -((int) xkstatus) ;
  dy_lp->inactzcorr += xk*orig_sys->obj[oxkndx] ;
/*
  Now it's time to walk the column and make any necessary corrections to the
  rhs (and perhaps rhslow).
*/
  if (flgon(orig_lp->ctlopts,lpctlRHSCHG|lpctlLBNDCHG|lpctlUBNDCHG))
  { ak = NULL ;
    if (consys_getcol_pk(orig_sys,oxkndx,&ak) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,"variable",
	     consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
      if (ak != NULL) pkvec_free(ak) ;
      return (FALSE) ; }
    for (ndx = 0, aik = &ak->coeffs[0] ; ndx < ak->cnt ; ndx++, aik++)
    { oaindx = aik->ndx ;
      aindx = dy_origcons[oaindx] ;
      if (aindx > 0)
      { dy_sys->rhs[aindx] -= aik->val*xk ;
	if (dy_sys->ctyp[aindx] == contypRNG)
	  dy_sys->rhslow[aindx] -= aik->val*xk ; } }
    pkvec_free(ak) ; }
/*
  And we're done. Print some information and return.
*/

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  %s (%d) %s inactive with value ",
	        consys_nme(orig_sys,'v',oxkndx,FALSE,NULL),oxkndx,
	        dy_prtvstat(xkstatus)) ;
    switch (xkstatus)
    { case vstatNBFX:
      case vstatNBLB:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"%g.",orig_sys->vlb[oxkndx]) ;
	break ; }
      case vstatNBUB:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"%g.",orig_sys->vub[oxkndx]) ;
	break ; }
      default:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"??.") ;
	break ; } } }
# endif

return (TRUE) ; }



static void process_active (lpprob_struct *orig_lp, int oxkndx)

/*
  This routine handles the data structure updates for an active variable
  x<k>.  We need to copy the new values for l<k>, u<k>, and c<k> into the
  active system. For nonbasic variables, we need to choose a status based on
  the bounds.  For basic variables, the status vector encodes the basis
  index, so we need to decide on an initial status --- either B, BFX, or BFR.

  The routine expects that bounds have been groomed (i.e., if the difference
  between l<k> and u<k> is less than the feasibility tolerance, they have been
  forced to exact equality).

  Parameters:
    orig_lp:	the original lp problem
    oxkndx:	index of x<k> in orig_sys
  
  Returns: undefined (the only possible error is a paranoid check)
*/

{ int xkndx ;
  double lk,uk,xk ;
  flags xkstatus ;
  consys_struct *orig_sys ;

# ifdef PARANOIA
  const char *rtnnme = "process_active" ;
# endif

  orig_sys = orig_lp->consys ;
/*
  Get the index of the variable in the active system, and the status.
  The paranoid check is that we're not attempting to convert between
  basic and nonbasic status.
*/
  xkndx = dy_origvars[oxkndx] ;
  xkstatus = dy_status[xkndx] ;

# ifdef PARANOIA
  if ((flgon(xkstatus,vstatBASIC) &&
       ((int) orig_lp->status[oxkndx]) > 0) ||
      (flgon(xkstatus,vstatNONBASIC|vstatNBFR) &&
       ((int) orig_lp->status[oxkndx]) < 0))
  { char buf[30] ;
    if (((int) orig_lp->status[oxkndx]) > 0)
      strcpy(buf,dy_prtvstat(orig_lp->status[oxkndx])) ;
    else
      strcpy(buf,"unspecified basic") ;
    errmsg(398,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',xkndx,FALSE,NULL),
	   xkndx,dy_prtvstat(xkstatus),buf) ;
    return ; }
# endif
/*
  Update the bounds and objective coefficient.
*/
  lk = orig_sys->vlb[oxkndx] ;
  dy_sys->vlb[xkndx] = lk ;
  uk = orig_sys->vub[oxkndx] ;
  dy_sys->vub[xkndx] = uk ;
  dy_sys->obj[xkndx] = orig_sys->obj[oxkndx] ;
/*
  For nonbasic variables, set the proper status based on the bounds and put
  the proper value in dy_x. Because the bounds can change arbitrarily and the
  client may not be maintaining the status vector, it's easiest to start from
  scratch, using the value from dy_x to decide the best new status.

  For basic variables, just decide between strictly basic (B), basic fixed
  (BFX), and basic free (BFR).  This will be correct, in the absence of bound
  changes, and the values held in dy_x and dy_xbasic are unchanged. If bounds
  have changed, we'll recalculate the primal variables and then decide on the
  final status of basic variables (which could be BLLB or BUUB).
*/
  if (flgon(dy_status[xkndx],vstatNONBASIC|vstatNBFR))
  { if (lk > -dy_tols->inf && uk < dy_tols->inf)
    { if (lk == uk)
      { xkstatus = vstatNBFX ;
	xk = lk ; }
      else
      if ((dy_x[xkndx] - lk) < (uk-dy_x[xkndx]))
      { xkstatus = vstatNBLB ;
	xk = lk ; }
      else
      { xkstatus = vstatNBUB ;
	xk = uk ; } }
    else
    if (lk > -dy_tols->inf)
    { xkstatus = vstatNBLB ;
      xk = lk ; }
    else
    if (uk < dy_tols->inf)
    { xkstatus = vstatNBUB ;
      xk = uk ; }
    else
    { xkstatus = vstatNBFR ;
      xk = 0 ; }
    dy_x[xkndx] = xk ; }
  else
  { if (lk == uk)
      xkstatus = vstatBFX ;
    else
    if (lk <= -dy_tols->inf && uk >= dy_tols->inf)
      xkstatus = vstatBFR ;
    else
      xkstatus = vstatB ; }
  dy_status[xkndx] = xkstatus ;
/*
  We're done. Print some information and return.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  %s (%d) %s active",
		consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		dy_prtvstat(dy_status[xkndx])) ;
    if (flgon(xkstatus,vstatNONBASIC|vstatNBFR))
      dyio_outfmt(dy_logchn,dy_gtxecho," with value %g.",dy_x[xkndx]) ;
    else
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
# endif

  return ; }




dyret_enum dy_hotstart (lpprob_struct *orig_lp)

/*
  This routine is responsible for handling a hot start. The assumption is
  that all data structures are in place, and that the user is allowed to change
  the bounds on variables and any of the rhs and objective coefficients. See
  the comments at the head of the file.

  Changes to the rhs and bounds are handled first. We reinstall the rhs
  array, then scan the variables, updating bounds and status and making the
  rhs corrections required for inactive variables.  If the bounds or rhs
  change, we need new primals. After we calculate new primals, we'll need to
  scan the basic variables and make sure their final status is correct.

  If the objective or bounds change, we need to recalculate the contribution
  to the objective from inactive variables. If the objective changes, we need
  new duals.  (It's also true that if the objective changes, we need new
  reduced costs, but that's handled in commonstart.)

  The most likely situation is that we haven't pivoted since refactoring as
  part of the preoptimality sequence, so we shouldn't need to refactor here.
  Instead, we leave it to dy_duenna to pick this up with the next pivot, as
  well as any possible accuracy check.

  Once all the changes have been incorporated, calculate primals and duals to
  determine primal and dual feasibility, and select the appropriate simplex
  phase in dy_lp->simplex.next.

  Parameters:
    orig_lp:	The original lp problem structure

  Returns: dyrOK if the setup completes without error, dyrINV or dyrFATAL
	   otherwise.
*/

{ int oxkndx,xkndx,oaindx,aindx ;
  double *ogvlb,*dyvlb,*ogvub,*dyvub,*ogobj,*dyobj,*dyrhs,*ogrhs ;
  consys_struct *orig_sys ;
  flags *ogstatus,calcflgs ;
  dyret_enum retval ;
  lpret_enum lpret ;
  const char *rtnnme = "dy_hotstart" ;

  /* dy_scaling.c */
  extern void dy_refreshlclsystem(flags what) ;

/*
  It could happen that there are no changes, in which case there's no point
  in going through the motions.
*/
  if (flgoff(orig_lp->ctlopts,
	     lpctlLBNDCHG|lpctlUBNDCHG|lpctlOBJCHG|lpctlRHSCHG))
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.crash >= 1)
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  no data structure changes at hot start.") ;
#   endif
    hot_updateMiscState(lpINV) ;
    return (dyrOK) ; }
/*
  But it's far more likely there are changes, and we need to get on with them.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  updating data structures at hot start ...") ;
    if (dy_opts->print.crash >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    scanning changes to") ;
      if (flgon(orig_lp->ctlopts,lpctlRHSCHG))
	dyio_outfmt(dy_logchn,dy_gtxecho," rhs") ;
      if (flgon(orig_lp->ctlopts,lpctlLBNDCHG))
	dyio_outfmt(dy_logchn,dy_gtxecho," vlb") ;
      if (flgon(orig_lp->ctlopts,lpctlUBNDCHG))
	dyio_outfmt(dy_logchn,dy_gtxecho," vub") ;
      if (flgon(orig_lp->ctlopts,lpctlOBJCHG))
	dyio_outfmt(dy_logchn,dy_gtxecho," obj") ;
      dyio_outfmt(dy_logchn,dy_gtxecho," ...") ; } }
#   endif

/*
  Transfer any changes from the client's system to the scaled local copy, if
  it exists. Then set up convenient handles for the various vectors.
*/
  dy_refreshlclsystem(orig_lp->ctlopts) ;

  orig_sys = orig_lp->consys ;
  dyrhs = dy_sys->rhs ;
  ogrhs = orig_sys->rhs ;
  ogvlb = orig_sys->vlb ;
  dyvlb = dy_sys->vlb ;
  ogvub = orig_sys->vub ;
  dyvub = dy_sys->vub ;
  ogobj = orig_sys->obj ;
  dyobj = dy_sys->obj ;
  ogstatus = orig_lp->status ;
/*
  If any of the rhs or bounds have been changed, we need to reinstall the rhs
  and bounds.  Begin by scanning the orig_sys rhs array, updating the dy_sys
  entries for the active constraints. If a range constraint comes by, we also
  need to set the upper bound of the associated logical.
*/
  if (flgon(orig_lp->ctlopts,lpctlLBNDCHG|lpctlUBNDCHG|lpctlRHSCHG))
  { 
    for (aindx = 1 ; aindx <= dy_sys->concnt ; aindx++)
    { oaindx = dy_actcons[aindx] ;
      if (oaindx > 0)
      { dyrhs[aindx] = ogrhs[oaindx] ;
	if (dy_sys->ctyp[aindx] == contypRNG)
	{ dy_sys->rhslow[aindx] = orig_sys->rhslow[oaindx] ;
	  dyvub[aindx] = dyrhs[aindx]-dy_sys->rhslow[aindx] ; } } } }
/*
  We need to scan the columns no matter what changed.  Objective coefficient
  changes are just copied into the active system as needed. The real action
  is updating bounds and dealing with the side effects of bounded variables.
    * Recalculate the contribution to inactzcorr for each inactive variable.
    * Update dy_sys->vlb, dy_sys->vub, and dy_sys->obj for each active
      variable.
    * Update dy_status for each active variable.
    * Update dy_x for each nonbasic active variable.

  Take the opportunity to groom the bounds, making sure that if the bounds
  are within the feasibility tolerance of one another, they are set exactly
  equal. This simplifies the handling of fixed variables.
*/
  dy_lp->inactzcorr = 0 ;
  lpret = lpINV ;
  for (oxkndx = 1 ; oxkndx <= orig_sys->varcnt ; oxkndx++)
  { xkndx = dy_origvars[oxkndx] ;
/*
  Force equality for bounds within the feasibility tolerance of one another.
  Then check for infeasibility (crossed bounds).
*/
    if (atbnd(ogvlb[oxkndx],ogvub[oxkndx]) && ogvlb[oxkndx] != ogvub[oxkndx])
    { 
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tForcing equal bound %g for %s (%d)",
		    (ogvlb[oxkndx]+ogvub[oxkndx])/2,
		    consys_nme(orig_sys,'v',oxkndx,0,0),oxkndx) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t  original lb = %g, ub = %g, diff = %g, tol = %g",
		    ogvlb[oxkndx],ogvub[oxkndx],ogvub[oxkndx]-ogvlb[oxkndx],
		    dy_tols->pfeas) ; }
#     endif
      ogvlb[oxkndx] = (ogvlb[oxkndx]+ogvub[oxkndx])/2 ;
      ogvub[oxkndx] = ogvlb[oxkndx] ; }
    if (ogvlb[oxkndx] > ogvub[oxkndx])
    { lpret = lpINFEAS ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		 "\n\tTrivial infeasibility for %s (%d), lb = %g > ub = %g.",
		 consys_nme(orig_sys,'v',oxkndx,0,0),oxkndx,
		 ogvlb[oxkndx],ogvub[oxkndx]) ; }
#     endif
    }
/*
  Inactive variables: update the status in dy_origvars and calculate the
  contribution to inactzcorr. If we've reloaded rhs and rhslow, correct
  them to account for the value of the variable.

  Active variables: touch up bounds for fixed variables, update vlb, vub,
  and obj arrays for dy_sys, update dy_status, and update dy_x for nonbasic
  variables.
*/
    if (xkndx < 0)
    { if (process_inactive(orig_lp,oxkndx) == FALSE) return (dyrFATAL) ; }
    else
    { process_active(orig_lp,oxkndx) ; } }
/*
  Now, what do we need? Calculate primal values first.  If we calculate new
  primal variables, we need to reset the status of the basic variables, which
  means we need to do a quick scan of the logicals to reset their status.
  Arguably this is not necessary if only the objective changed, but overall
  it's a good investment of our time.
*/
  if (dy_calcprimals() == FALSE)
  { errmsg(316,rtnnme,dy_sys->nme) ;
    return (dyrFATAL) ; }
  for (xkndx = 1 ; xkndx <= dy_sys->concnt ; xkndx++)
  { if (dy_var2basis[xkndx] != 0)
    { if (dyvub[xkndx] == dyvlb[xkndx])
	dy_status[xkndx] = vstatBFX ;
      else
	dy_status[xkndx] = vstatB ; } }
  dy_setfinalstatus() ;
/*
  Is the phase I objective installed? If so, remove it. This hurts a bit,
  particularly if we ultimately end up targetting primal phase I as the
  starting simplex, but it's the only way to test for a dual feasible start.
  And if we have dual feasibility, it's a big win.
*/
  if (dy_lp->p1obj.installed == TRUE)
  { if (dy_swapobjs(dyPRIMAL2) == FALSE)
    { errmsg(318,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"remove") ;
      return (dyrFATAL) ; } }
/*
  Calculate duals and reduced costs and see if we're primal or dual feasible.
  Calculate the objective just for kicks.
*/
  dy_calcduals() ;
  if (dy_calccbar() == FALSE)
  { errmsg(384,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
    return (dyrFATAL) ; }
  dy_lp->z = dy_calcobj() ;

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
/*
  Reset a few control variables and counts in dy_lp.
*/
  hot_updateMiscState(lpret) ;
/*
  And that should do it. Let's make a paranoid check or two, then we're
  off and running.
*/
# ifdef PARANOIA
  if (dy_chkdysys(orig_sys) == FALSE) return (dyrFATAL) ;
# endif

  return (dyrOK) ; }

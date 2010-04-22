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
  This file contains routines for pricing columns and for calculating up and
  down penalties associated with forcing basic variables to new bounds.
  They're useful in the context of using dylp in an lp-based branch-and-bound
  MILP code. These routines shield the client from the the problems of
  translating between the original system and the active system (selection of
  active constraints and variables, and scaling).

  To state the obvious, these routines are useable only when dylp's data
  structures have been left intact after solving the lp.
*/

/*
  In penalty calculations, the core calculation is the cost of the first dual
  pivot back toward feasibility after imposing new bounds on the variable.
  Typically we're trying to estimate the result of branching on x<i>, and
  we're evaluating the mutually exclusive cases of x<i> >= ub<i> and x<i> <=
  lb<i>. The nonbasic variables x<k> are scanned to find the variable x<j>
  which can drive x<i> to bound with minimum degradation of the objective.

  For the common case where x<i> is an integer variable and we want to know
  the standard up/down penalties, ub<i> = ceil(x<i> and lb<i> = floor(x<i>).

  If x<j> is an integer variable, it will most likely end up with a
  non-integral value.  One can make the argument that x<j> really needs to
  move all the way to its next integral value in order to be a feasible
  solution to the ILP, hence we should calculate the penalty as if we moved
  to the next integral value. The problem, of course, is that abar<j> might
  not really be the best direction to move in order to achieve this.
  Moreover, the penalty calculated in this way may well be more than the
  actual cost of regaining optimality for the lp relaxation. More often,
  however, the first dual pivot is an underestimate of the real degradation
  and the strengthened penalty only partially corrects this.

  dylp can go either way. If the compile-time constant STRONG_PENALTIES is
  defined, you get penalties strengthened as described for x<j> integer.
  Otherwise, the basic penalty is calculated.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_penalty.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_penalty.c 269 2009-04-02 05:38:19Z lou $" ;



/*
  This is a dummy stub routine for the benefit of optimised builds. It'll do
  until I rewrite the original.
*/

static bool dy_unscale_betai(consys_struct *orig_sys, int oxindx,
			       double **betai, double **ai)

{ return (FALSE) ; }



/*
  This routine has been rewritten to provide the the full unscaled vector of
  reduced costs. I've moved the original code over here and renamed it
  for future reference.
*/

static void dy_orig_cbarLocal (int nbcnt, double *cbar, int *vndx)

/*
  This is a special purpose routine which unscales the vector of selected
  reduced costs produced by dy_pricenbvars. All we do here is walk the vectors
  and apply the unscaling.

  sc_cbar<j> = sc_c<j> - sc_c<B>sc_inv(B)sc_a<j>
	     = c<j>S<j> - c<B>S<B>inv(S<B>)inv(B)inv(R)Ra<j>S<j>
	     = c<j>S<j> - c<B>inv(B)a<j>S<j>
	     = cbar<j>S<j>

  To unscale sc_cbar<j>, we simply multiply by 1/S<j>, keeping in mind that
  if x<j> is a logical for row i, the appropriate factor is R<i>.

  Parameters:
    nbcnt:	number of entries in cbar, nbvars
    cbar:	vector of reduced costs
    vndx:	corresponding variable indices

  Note that cbar and vndx are indexed from 0, and the indices specified in
  vndx are in the frame of the original constraint system, which is what we
  need for accesses to the scaling vectors.

  Returns: undefined
*/

{ int j,k ;
  double cbarj ;
  const double *rscale,*cscale ;

/*
  Is unscaling required? If so, acquire the vectors and go to it.
*/
  if (dy_isscaled() == FALSE) return ;
  dy_scaling_vectors(&rscale,&cscale) ;
/*
  Get on with the calculation.  Recall that vndx encodes the index of a
  logical as -i.
*/
  for (k = 0 ; k < nbcnt ; k++)
  { j = vndx[k] ;
    cbarj = cbar[k] ;
    if (j > 0)
    { cbarj /= cscale[j] ; }
    else
    { cbarj *= rscale[-j] ; }
    setcleanzero(cbarj,dy_tols->dfeas) ;
    cbar[k] = cbarj ; }

  return ; }



bool dy_pricenbvars (lpprob_struct *orig_lp, flags priceme,
		     double **p_ocbar, int *p_nbcnt, int **p_nbvars)

/*
  This routine will calculate the reduced costs of nonbasic variables in the
  original constraint system, using the information available in the dylp
  data structures. For active variables, the value from dy_cbar is used.  For
  inactive variables, the column is priced. The set of variables priced is
  all variables whose status is specified with priceme. On return, nbcnt has
  the number of columns actually priced, nbvars holds the indices, and ocbar
  holds the corresponding reduced costs. E.g., if x<j> is the ith variable
  priced, nbvars[i] = j and ocbar[i] = cbar<j>.

  Variable indices returned in nbvars are in the orig_sys frame, which does
  not include logicals. Priced logicals are indicated by the negative of the
  index of their associated constraint.

  The overall flow of the routine is to calculate the scaled cbar's, then
  unscale them as a final step.

  Parameters:
    orig_lp:	the lp problem structure
    priceme:	variables with a status included in priceme will be priced
    p_ocbar:	(i) a vector to hold the reduced costs; if NULL, one will
		be allocated
		(o) the reduced costs for priced nonbasic variables; only
		entries for the indices in nbvars are valid
    p_nbcnt:	(o) the number of variables priced
    p_nbvars	(i) a vector to hold the indices of the variables which are
		priced; if NULL, one will be allocated
		(o) the indices of the priced variables

  Returns: TRUE if there are no problems with pricing, FALSE otherwise.
*/

{ int oxjndx,xjndx,pkndx,cndx,nbcnt,*nbvars ;
  double cbarj,*ocbar ;
  flags statj ;
  bool retval ;
  consys_struct *orig_sys ;
  pkvec_struct *aj ;
  pkcoeff_struct *aij ;

  const char *rtnnme = "dy_pricenbvars" ;

  /* dy_unscaling.c */
  extern void dy_orig_cbarLocal (int nbcnt, double *cbar, int *vndx) ;

# ifdef DYLP_PARANOIA
  if (p_ocbar == NULL)
  { errmsg(2,rtnnme,"&cbar") ;
    return (FALSE) ; }
  if (p_nbcnt == NULL)
  { errmsg(2,rtnnme,"&nbcnt") ;
    return (FALSE) ; }
  if (p_nbvars == NULL)
  { errmsg(2,rtnnme,"&nbvars") ;
    return (FALSE) ; }
  if (orig_lp == NULL)
  { errmsg(2,rtnnme,"orig_lp") ;
    return (FALSE) ; }
  if (orig_lp->consys == NULL)
  { errmsg(2,rtnnme,"orig_lp->consys") ;
    return (FALSE) ; }
# endif

/*
  Check for valid data structures, then pull out the constraint system. The
  call to initlclsystem will replace the client's copy of the original
  constraint system with the local scaled copy, if it exists.
*/
  if (flgoff(orig_lp->ctlopts,lpctlDYVALID))
  { errmsg(396,rtnnme,orig_lp->consys->nme,"price nonbasic columns") ;
    return (FALSE) ; }
  (void) dy_initlclsystem(orig_lp,TRUE) ;
  orig_sys = orig_lp->consys ;
/*
  Did the client give us vectors, or do we need to allocate them? There can
  be at most varcnt nonbasic variables.
*/
  if (*p_ocbar == NULL)
    *p_ocbar = (double *) CALLOC(orig_sys->varcnt,sizeof(double)) ;
  ocbar = *p_ocbar ;
  if (*p_nbvars == NULL)
    *p_nbvars = (int *) CALLOC(orig_sys->varcnt,sizeof(int)) ;
  nbvars = *p_nbvars ;
/*
  Open a loop to walk the columns of orig_sys.

  Active variables have a valid dy_sys index (> 0) in dy_origvars. We simply
  check dy_status and grab the entry in dy_cbar. To simplify later use, set
  cbar<j> to 0 if it's less than the dual feasibility tolerance.
*/
  retval = TRUE ;
  aj = NULL ;
  nbcnt = 0 ;
  for (oxjndx = 1 ; oxjndx <= orig_sys->varcnt ; oxjndx++)
  { if (ACTIVE_VAR(oxjndx))
    { xjndx = dy_origvars[oxjndx] ;
      statj = dy_status[xjndx] ;
      if (flgon(statj,priceme))
      { cbarj = dy_cbar[xjndx] ;
	setcleanzero(cbarj,dy_tols->dfeas) ;
#       ifdef DYLP_PARANOIA
	if ((flgon(statj,vstatNBUB) && cbarj > 0) ||
	    (flgon(statj,vstatNBLB) && cbarj < 0))
	{ errmsg(739,rtnnme,dy_sys->nme,"active",
		 consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx,
		 dy_prtvstat(statj),cbarj) ;
	  retval = FALSE ;
	  break ; }
#       endif
	ocbar[nbcnt] = cbarj ;
	nbvars[nbcnt++] = oxjndx ; } }
/*
  Inactive variables have the negative of their status (nonbasic at upper or
  lower bound, or fixed) in dy_origvars. We have to price the column in this
  case. Calculate cbar<j> = c<j> - y<i>a<i,j>, using only the active rows of
  the column.
*/
    else
    { statj = (flags) -dy_origvars[oxjndx] ;
#     ifdef DYLP_PARANOIA
      if (flgoff(statj,vstatNBFX|vstatNBUB|vstatNBLB|vstatNBFR))
      { errmsg(433,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "inactive",consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),
	       oxjndx,dy_prtvstat(statj)) ;
	retval = FALSE ;
	break ; }
#     endif
      if (flgon(statj,priceme))
      { cbarj = orig_sys->obj[oxjndx] ;
	if (consys_getcol_pk(orig_sys,oxjndx,&aj) == FALSE)
	{ errmsg(122,rtnnme,orig_sys->nme,"column",
		 consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
	  retval = FALSE ;
	  break ; }
	for (pkndx = 0,aij = aj->coeffs ; pkndx < aj->cnt ; pkndx++,aij++)
	{ if (ACTIVE_CON(aij->ndx))
	  { cndx = dy_origcons[aij->ndx] ;
	    cbarj -= dy_y[cndx]*aij->val ; } }
	setcleanzero(cbarj,dy_tols->dfeas) ;
#       ifdef DYLP_PARANOIA
	if ((flgon(statj,vstatNBUB) && cbarj > 0) ||
	    (flgon(statj,vstatNBLB) && cbarj < 0))
	{ errmsg(739,rtnnme,dy_sys->nme,"inactive",
		 consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx,
		 dy_prtvstat(statj),cbarj) ;
	  retval = FALSE ;
	  break ; }
#       endif
	ocbar[nbcnt] = cbarj ;
	nbvars[nbcnt++] = oxjndx ; } } }
/*
  Last but not least, step through the logicals (in the dy_sys frame) and add
  them, using the negative of their associated constraint (in the orig_sys
  frame) as the index.
*/
  for (xjndx = 1 ; xjndx <= dy_sys->concnt ; xjndx++)
  { statj = dy_status[xjndx] ;
    if (flgon(statj,priceme))
    { cbarj = dy_cbar[xjndx] ;
      setcleanzero(cbarj,dy_tols->dfeas) ;
#     ifdef DYLP_PARANOIA
      if ((flgon(statj,vstatNBUB) && cbarj > 0) ||
	  (flgon(statj,vstatNBLB) && cbarj < 0))
      { errmsg(739,rtnnme,dy_sys->nme,"logical",
	       consys_nme(dy_sys,'v',xjndx,TRUE,NULL),xjndx,
	       dy_prtvstat(statj),cbarj) ;
	retval = FALSE ;
	break ; }
#     endif
      ocbar[nbcnt] = cbarj ;
      nbvars[nbcnt++] = -dy_actcons[xjndx] ; } }
  *p_nbcnt = nbcnt ;
/*
  Unscale the reduced costs.
*/
  (void) dy_orig_cbarLocal(nbcnt,ocbar,nbvars) ;
/*
  And we're done. Clean up and return.
*/
  dy_freelclsystem(orig_lp,FALSE) ;
  if (aj != NULL) pkvec_free(aj) ;

  return (retval) ; }



static bool pricedualpiv (consys_struct *orig_sys, double *betai, double *ai,
			  int oxindx, double nubi, double xi, double nlbi,
			  int nbcnt, int *nbvars, double *nbcbar,
			  double *p_upeni, double *p_dpeni)

/*
  This routine performs the core pricing calculation for the up and down
  penalties associated with the dual pivot that will force x<i> to decrease
  and leave the basis at nub<i> or increase and leave the basis at nlb<i>.
  The assumption is that we're pricing a prospective disjunction.

  If x<i> is a logical, and the associated constraint is inactive (i.e., not
  part of the current basis), a<i> must contain the coefficients of the
  constraint. (Note that this will also work in the case that we want to
  price a new constraint, not currently in orig_sys.)

  It is expected that any necessary setup for the calculation has been taken
  care of by a front-end routine (at present, only dy_pricedualpiv).  The
  most likely source of a scaled copy of the original constraint system is
  the local copy maintained by dy_scaling.c The most likely source of an
  unscaled copy is the client's copy of the original constraint system.
  Neither system is likely to have explicit logicals, hence this routine
  handles the columns for logicals implicitly, to avoid any problems.

  Here's the computation. For convenience, assume that x<i> is basic in basis
  position i, and this row is constraint i.

  dpen<i> = MIN{j}(-(nlb(i) - x<i>)cbar<j>/abar<ij>)
  
    for j in nbvars s.t.
	abar<ij> > 0 and x<j> < u<j> or abar<ij> < 0 and x<j> > l<j>
  
  upen<i> = MIN{j}(-(nub(i) - x<i>)cbar<j>/abar<ij>)

    for j in nbvars s.t.
	abar<ij> < 0 and x<j> < u<j> or abar<ij> > 0 and x<j> > l<j>

  See the written documentation for development of the formulas. Briefly, the
  penalty calculation gives the deterioration in the objective for the first
  dual pivot back toward feasibility.

  The routine calculates abar<ij> = dot(beta<i>,a<j>) for each j in nbvars.
  It will keep going until dpen<i> = upen<i> = 0, or all variables in nbvars
  have been scanned.

  nbvars is assumed to contain the indices of the nonbasic variables `of
  interest', meaning that they're not fixed. nbcbar is in correspondence, so
  that if nbvars[i] = j, nbcbar[i] = cbar<j>. Logicals (which don't exist in
  the orig_sys frame of reference) are expected to be coded as the negative
  of the index of their associated constraint. (See dy_pricenbvars.) nbvars
  and nbcbar are indexed from 0.

  Parameters:
    orig_sys:	the unscaled original constraint system
    betai:	unscaled row i of the basis inverse
    ai:		unscaled coefficients of constraint i, if it is not active
    oxindx:	index of x<i> in orig_sys; negative values are assumed to
		specify logicals (as the negative of the index of the
		associated constraint); 0 means a<i> holds a new row.
    nubi:	unscaled new upper bound for x<i>
    xi:		current value of x<i> in optimal solution
    nlbi:	unscaled new lower bound for x<i>
    nbcnt:	the number of variable indices in nbvars
    nbvars:	indices of the nonbasic variables of interest (indices in the
		orig_sys frame)
    nbcbar:	the reduced costs, in correspondence with nbvars
    p_upeni:	(o) the up penalty when x<i> is forced up to nlb<i>
    p_dpeni:	(o) the down penalty when x<i> is forced down to nub<i>

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise
*/

{ 

#if 0

  int oxjndx,nbndx,pkndx,cndx,vndx ;
  double abarij,upenij,dpenij,upeni,dpeni,cbarj ;
  flags statj ;
  vartyp_enum *vtyp ;
  pkvec_struct *aj ;
  pkcoeff_struct *aij ;
  bool activexi,retval ;

  const char *rtnnme = "pricedualpiv" ;


# ifdef DYLP_PARANOIA
  int i,ipos ;
  double *abarj ;

  /* dy_tableau.c */
  bool dy_abarj(consys_struct *orig_sys, int j_orig, double **p_abarj) ;

  abarj = NULL ;
# endif

  vtyp = orig_sys->vtyp ;

/*
  Do we have an active variable or an inactive variable? If oxindx < 0 it's a
  logical, so we need to ask whether or not the constraint is active. We
  don't deal with inactive architectural variables.
*/
  activexi = TRUE ;
  if (oxindx == 0)
  { activexi = FALSE ; }
  else
  if (oxindx < 0)
  { if (INACTIVE_CON(-oxindx)) activexi = FALSE ; }
  else
  { if (INACTIVE_VAR(oxindx))
    { errmsg(737,rtnnme,orig_sys->nme,
	     consys_nme(orig_sys,'v',oxindx,FALSE,NULL),oxindx) ;
      return (FALSE) ; } }
# ifdef DYLP_PARANOIA
  if (activexi == FALSE && ai == NULL)
  { errmsg(2,rtnnme,"a<i> (inactive)") ;
    return (FALSE) ; }
# endif
/*
  Now open a loop to walk nbvars.
  
  Architectural variables are specified with their index in orig_sys, and may
  or may not be active.  For each architectural variable, pull the unscaled
  column, calculate abar<ij> = dot(beta<i>,a<j>), and acquire the status.

  Logicals are specified with the index of their associated constraint in
  orig_sys, negated, and will certainly be active, but we deal with them
  implicitly. (This is important, because orig_sys might well be the client's
  copy and logicals may not be explicitly represented.)

  If we're dealing with an inactive row, remember that it may also contain a
  coefficient for the variable in this column.

  If abar<ij> is 0, that's all the further we need to go with this x<j>,
  because we can't use it to move x<i>. The escape comes after the paranoid
  check.
*/
  dpeni = dy_tols->inf ;
  upeni = dy_tols->inf ;
  aj = NULL ;
  retval = TRUE ;
  for (nbndx = 0 ; nbndx < nbcnt ; nbndx++)
  { oxjndx = nbvars[nbndx] ;
    if (oxjndx < 0)
    { cndx = dy_origcons[-oxjndx] ;
      if (dy_sys->ctyp[cndx] == contypGE)
      { abarij = -betai[cndx] ; }
      else
      { abarij = betai[cndx] ; }
      statj = dy_status[cndx] ; }
    else
    { if (consys_getcol_pk(orig_sys,oxjndx,&aj) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,"column",
	       consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
	retval = FALSE ;
	break ; }
      abarij = 0 ;
      for (pkndx = 0,aij = aj->coeffs ; pkndx < aj->cnt ; pkndx++,aij++)
      { if (ACTIVE_CON(aij->ndx))
	{ cndx = dy_origcons[aij->ndx] ;
	  abarij += betai[cndx]*aij->val ; } }
      vndx = dy_origvars[oxjndx] ;
      if (vndx > 0)
      { statj = dy_status[vndx] ; }
      else
      { statj = (flags) -vndx ; }
      if (activexi == FALSE && oxjndx > 0)
	abarij += betai[dy_sys->concnt+1]*ai[oxjndx] ; }
    setcleanzero(abarij,dy_tols->zero) ;

#   ifdef DYLP_PARANOIA
/*
  We can do a check if the row is active, by ftran'ing the column and checking
  that we get the same value for abar<ij> with both calculations. The tolerance
  on the check --- 1000*dy_tols.zero ---  is pretty loose, but remember that
  the original system is unscaled and could be pretty ugly.

  This code is broken at the moment, because I'm changing the semantics of
  dy_abarj to return the full abar<j> in the context of the original system.
  See dy_tableau.c.	-- lh, 080515 --
*/
    if (activexi)
    { if (dy_abarj(orig_sys,oxjndx,&abarj) == FALSE)
      { if (oxjndx < 0)
	{ vndx = orig_sys->varcnt-oxjndx ; }
	else
	{ vndx = oxjndx ; }
	errmsg(740,rtnnme,orig_sys->nme,"ftran'd column",
	       consys_nme(orig_sys,'v',vndx,TRUE,NULL),oxjndx) ;
	retval = FALSE ;
	break ; }
      if (oxindx < 0)
      { i = dy_origcons[-oxindx] ; }
      else
      { i = dy_origvars[oxindx] ; }
      ipos = dy_var2basis[i] ;
    if (!withintol(abarj[ipos],abarij,1000*dy_tols->zero))
    { errmsg(741,rtnnme,orig_sys->nme,oxindx,oxjndx,abarij,ipos,abarj[ipos],
	     abarj[ipos]-abarij,1000*dy_tols->zero) ;
      retval = FALSE ;
      break ; } }
#   endif
      
    if (abarij == 0) continue ;
/*
  Calculate the penalty dpen<ij> for using x<j> to force x<i> to nub<i>.  We
  can use x<j> only if abar<ij> has the right sign to match the status of
  x<j> (dual pivoting rules, essentially).  dpen<i> becomes
  min(dpen<i>,dpen<ij>).  If dpen<i> is already 0, it can't get worse, so
  there's no need to do the calculation.

  The drill is the same for upen<i>, using nlb<i>.
*/
    cbarj = nbcbar[nbndx] ;
#   ifdef DYLP_PARANOIA
    if (cbarj != 0 && fabs(cbarj) < dy_tols->dfeas)
    { int tmpndx ;
      if (oxjndx < 0)
	tmpndx = -oxjndx ;
      else
	tmpndx = oxjndx ;
      errmsg(738,rtnnme,orig_sys->nme,tmpndx,cbarj,dy_prtvstat(statj),
	     consys_nme(orig_sys,'v',tmpndx+orig_sys->varcnt,FALSE,NULL),
	     tmpndx,dy_tols->dfeas,dy_tols->dfeas-fabs(cbarj)) ;
      retval = FALSE ;
      break ; }
#   endif
    if (dpeni > 0)
    { if ((abarij > 0 && flgoff(statj,vstatNBUB)) ||
	  (abarij < 0 && flgoff(statj,vstatNBLB)))
      { dpenij = (nubi-xi)*cbarj/(-abarij) ;
	setcleanzero(dpenij,dy_tols->zero) ;
#	ifdef STRONG_PENALTIES
	if (oxjndx > 0)
	{ if (INT_VARTYPE(vtyp[oxjndx]) && dpenij < fabs(cbarj))
	    dpenij = fabs(cbarj) ; }
#	endif
#       ifdef DYLP_PARANOIA
	if (dpenij < 0)
	{ errmsg(736,rtnnme,orig_sys->nme,"dpen",oxindx,oxjndx,dpenij) ;
	  retval = FALSE ;
	  break ; }
#	endif
	if (dpenij < dpeni) dpeni = dpenij ; } }
    if (upeni > 0)
    { if ((abarij < 0 && flgoff(statj,vstatNBUB)) ||
	  (abarij > 0 && flgoff(statj,vstatNBLB)))
      { upenij = (nlbi-xi)*cbarj/(-abarij) ;
	setcleanzero(upenij,dy_tols->zero) ;
#	ifdef STRONG_PENALTIES
	if (oxjndx > 0)
	{ if (INT_VARTYPE(vtyp[oxjndx]) && upenij < fabs(cbarj))
	    upenij = fabs(cbarj) ; }
#	endif
#       ifdef DYLP_PARANOIA
	if (upenij < 0)
	{ errmsg(736,rtnnme,orig_sys->nme,"upen",oxindx,oxjndx,upenij) ;
	  retval = FALSE ;
	  break ; }
#	endif
	if (upenij < upeni) upeni = upenij ; } }
/*
  If both upeni and dpeni have been reduced to 0, we can quit.
*/
    if (upeni == 0 && dpeni == 0) break ; }
/*
  We've finished the penalty calculation loop. Clean up, set the return
  values, and return.
*/
  if (aj != NULL) pkvec_free(aj) ;

# ifdef DYLP_PARANOIA
  if (abarj != NULL) FREE(abarj) ;
# endif

  *p_upeni = upeni ;
  *p_dpeni = dpeni ;

  return (retval) ;

#endif

  return (FALSE) ; }




bool dy_pricedualpiv (lpprob_struct *orig_lp, int oxindx,
		      double nubi, double xi, double nlbi,
		      int nbcnt, int *nbvars,
		      double *nbcbar, double *p_upeni, double *p_dpeni)

/*
  This routine is a generalised pricing routine that calculates the up and
  down penalties associated with the dual pivot that will force x<i> (basic,
  infeasible) to either rise to nlb<i> and leave the basis, or fall to nub<i>
  and leave the basis. It's generalised (e.g., with respect to a standard
  up/down penalty calculation) in the sense that it doesn't assume floor and
  ceiling for nub<i> and nlb<i>, respectively, and it can handle the
  calculation when x<i> is the slack associated with an inactive constraint.
  (This is handy for evaluating candidates for branch-on-hyperplane.)

  dpen<i> = MIN{j}(-(nlb(i) - x<i>)cbar<j>/abar<ij>)
  
    for j in nbvars s.t.
	abar<ij> > 0 and x<j> < u<j> or abar<ij> < 0 and x<j> > l<j>
  
  upen<i> = MIN{j}(-(nub(i) - x<i>)cbar<j>/abar<ij>)

    for j in nbvars s.t.
	abar<ij> < 0 and x<j> < u<j> or abar<ij> > 0 and x<j> > l<j>

  Each abar<ij> is calculated as dot(beta<i>,a<j>).

  See the written documentation for development of the formulas. Briefly, 
  the penalty calculation gives the deterioration in the objective for the
  first dual pivot back toward feasibility.

  The calculation of an unscaled basis inverse row beta<i> is handled by
  dy_unscale_betai. The remainder of the heavy lifting is handled by
  pricedualpiv.

  nbvars is assumed to contain the nonbasic variables `of interest', meaning
  that they're not fixed. nbcbar is in correspondence, so that if nbvars[i] =
  j, nbcbar[i] = cbar<j>. Logicals (which don't exist in the orig_sys frame
  of reference) are expected to be coded as the negative of the index of
  their associated constraint. (See dy_pricenbvars.)

  Parameters:
    orig_lp:	the lp problem structure; contains a pointer to the unscaled
		original constraint system
    oxindx:	index of x<i> in orig_sys; negative values are assumed to
		specify logicals (as the negative of the index of the
		associated constraint)
    nubi:	new upper bound for x<i>
    xi:		current value for x<i>
    nlbi:	new lower bound for x<i>
    nbcnt:	the number of variable indices in nbvars
    nbvars:	indices of the nonbasic variables of interest (indices in the
		orig_sys frame)
    nbcbar:	the reduced costs, in correspondence with nbvars
    p_upeni:	(o) the up penalty for x<i>
    p_dpeni:	(o) the down penalty for x<i>

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise
*/

{ double *betai,*ai ;
  consys_struct *orig_sys ;
  bool retval ;

  const char *rtnnme = "dy_pricedualpiv" ;

  /* dy_scaling.c */
  extern bool dy_unscale_betai(consys_struct *orig_sys, int oxindx,
			       double **betai, double **ai) ;
# ifdef DYLP_PARANOIA
  if (p_upeni == NULL)
  { errmsg(2,rtnnme,"&upen<i>") ;
    return (FALSE) ; }
  if (p_dpeni == NULL)
  { errmsg(2,rtnnme,"&dpen<i>") ;
    return (FALSE) ; }
  if (orig_lp == NULL)
  { errmsg(2,rtnnme,"orig_lp") ;
    return (FALSE) ; }
  if (orig_lp->consys == NULL)
  { errmsg(2,rtnnme,"orig_lp->consys") ;
    return (FALSE) ; }
# endif

/*
  Check for valid data structures, then pull out the constraint system. We want
  the client's (unscaled) copy, so we don't swap in the scaled local copy.
*/
  if (flgoff(orig_lp->ctlopts,lpctlDYVALID))
  { errmsg(396,rtnnme,orig_lp->consys->nme,"calculate penalty") ;
    return (FALSE) ; }
  orig_sys = orig_lp->consys ;
/*
  Get the unscaled row of the basis inverse.
*/
  betai = NULL ;
  ai = NULL ;
  if (dy_unscale_betai(orig_sys,oxindx,&betai,&ai) == FALSE)
  { return (FALSE) ;  }
/*
  Call pricedualpiv to do the rest of the heavy lifting.
*/
  retval = pricedualpiv(orig_sys,betai,ai,oxindx,nubi,xi,nlbi,
			nbcnt,nbvars,nbcbar,p_upeni,p_dpeni) ;
/*
  Clean up and return.
*/
  if (ai != NULL) FREE(ai) ;
  if (betai != NULL) FREE(betai) ;

  return (retval) ; }

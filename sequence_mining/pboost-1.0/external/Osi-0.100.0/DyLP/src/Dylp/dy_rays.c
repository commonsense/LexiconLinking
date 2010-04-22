/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2008 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This file contains routines to return primal and dual rays, expressed in
  terms of the original system.

  If you think back to LP duality, (primal, dual) unbounded implies (dual,
  primal) infeasibility.  The reverse does not hold: it's possible to be
  primal and dual infeasible.  Still, it takes some pretty special pathology
  to manage simultaneous infeasibility (primal infeasibility defined by
  parallel constraints, so that we'd be unbounded if only we could achieve
  feasibility).

  If dylp returns primal or dual unbounded, then we know there's a ray in
  there somewhere. In fact, you'll never see dual unbounded, as return codes
  are always expressed in terms of the primal problem, and dual unbounded is
  converted to primal infeasible.
  
  Dylp's dual simplex doesn't have a phase I, so it's impossible to get to
  simultaneous dual and primal infeasibility by that route.  Primal simplex
  does have a phase I and may well return primal infeasible without ever
  running dual simplex. Still, given the rarity of simultaneous dual and
  primal infeasibility, I'm not inclined to try and juggle return codes to
  explicitly distinguish {primal, dual} {unbounded, infeasible}.

  These routines assume that dylp has run to completion, and returned
    * unbounded, in which case there will be primal rays, or
    * infeasible, in which case there will (with high probability) be dual
      rays.
  The infeasible case can be refined a bit by checking dy_lp->active to see
  if the last simplex phase was dual or primal. If it was dual, we have a ray
  with probability 1.

  We're only going to return rays that emanate from the current basic
  solution. A search for all rays is a nontrivial and expensive task.

  With these assumptions, we can ignore any interactions with the inactive
  parts of the constraint system when testing for ray-ness.  That doesn't
  mean that ray components associated with inactive constraints or variables
  are all zero, only that the inactive components will not block the ray. If
  they could, dylp would have activated them and we wouldn't be here.

  There's another bit of pathology to consider: It's possible to define a
  problem where we can move infinitely far in either direction. For example,

    min x1  s.t.  x1 + x2 = 0   x1, x2 free

  However, only the direction that improves the objective qualifies as a
  ray.  This follows by extension from the very common case where the
  constraints of an lp form a polyhedral cone with a finite optimum for, say,
  minimisation, but changing to maximisation results in an unbounded problem.
  Clearly, moving in a direction that degrades the objective does not
  normally qualify as a ray.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char svnid[] UNUSED = "$Id: dy_rays.c 269 2009-04-02 05:38:19Z lou $" ;

/* dy_tableau.c */

#if DYLP_PARANOIA > 0
extern bool dy_std_paranoia(const lpprob_struct *orig_lp,const char *rtnnme) ;
#endif


static bool testForPrimalRay (int j, int *p_dir, double **p_abarj)

/*
  This routine evaluates an active column abar<j> = inv(B)a<j> to determine
  if it constitutes a primal ray.  It's a much-simplified version of the
  algorithm that selects the leaving variable for primal simplex.

    * Instead of searching for the smallest limit on the change in x<j>, we
      just want to know if there's no limit. Hence we can abandon the
      evaluation as soon as any basic primal variable limits the change
      in x<j>.

    * There's no need to worry about degeneracy, numerical stability and all
      the other little considerations that go into selecting a good pivot.

  NOTE: This routine works in the active system reference frame. The vector
	returned is abar<j>. This is not quite a ray: it needs a coefficient
	(-1.0) for x<j> and it must be negated to become the corresponding
	ray. It makes more sense to do this in the client, once we know how
	the ray is to be presented. See, for example, dy_primalRays, which
	transforms the ray for use by the outside world.

  Parameters:
    j:		Index of the column to be evaluated
    p_dir:	(o) Direction of motion,
		  1: ray in the direction abar<j> (x<j> increasing)
		  0: not a ray
		 =1: ray in the direction -abar<j> (x<j> decreasing)
    p_abarj:	(o) if non-NULL, used to return abar<j> as an expanded vector
		    in basis order iff abar<j> is a primal ray, otherwise
		    *p_abarj is set to NULL
  
  Returns: TRUE if the evaluation completes without error, FALSE otherwise.
*/

{ int k,m,kpos,dir ;
  flags statj,statk ;
  double *abarj ;
  double abarkj,cbarj ;
  bool rayUp,rayDown ;
  
  const char *rtnnme = "testForPrimalRay" ;

# ifndef DYLP_NDEBUG
  int v ;

  if (dy_opts->print.rays >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n      Testing if column %s (%d) is a primal ray",
	        consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
# endif

  if (p_abarj != NULL) *p_abarj = NULL ;
  *p_dir = 0 ;
  rayUp = TRUE ;
  rayDown = TRUE ;
/*
  Start by checking the reduced cost.  If cbar<j> = 0, we can change x<j> all
  we like, but we won't go unbounded.  Use the reduced cost / objective
  tolerance here, otherwise we could be in disagreement with dylp. If the
  reduced cost is not zero, then any ray must improve the objective. Dylp is
  always minimising: min  z = c<B>inv(B)b + cbar<N>x<N>, where
  cbar<N> = (c<N> - c<B>inv(B)N).
*/
  cbarj = dy_cbar[j] ;
  if (withintol(cbarj,0,dy_tols->cost))
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  ".\n\tcbar<%d> = %g; no ray.",j,cbarj) ; }
#   endif
    return (TRUE) ; }
  else
  if (cbarj < 0)
  { rayDown = FALSE ;
    dir = 1 ; }
  else
  { rayUp = FALSE ;
    dir = -1 ; }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.rays >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,".\n\tcbar<%d> = %g allows %s ray",
		j,cbarj,((dir < 0)?"down":"up")) ; }
# endif
/*
  Reduce the possibilities based on the status of x<j>. We shouldn't be called
  for basic variables, but bail out now if that's the case.
*/
  statj = getflg(dy_status[j],vstatSTATUS) ;
  if (flgon(statj,vstatBASIC))
  { rayUp = FALSE ;
    rayDown = FALSE ; }
  else
  { switch (statj)
    { case vstatNBFX:
      { rayUp = FALSE ;
	rayDown = FALSE ;
	break ; }
      case vstatNBLB:
      { rayDown = FALSE ;
	if (dy_sys->vub[j] < dy_tols->inf)
	{ rayUp = FALSE ; }
	break ; }
      case vstatNBUB:
      { rayUp = FALSE ;
	if (dy_sys->vlb[j] > -dy_tols->inf)
	{ rayDown = FALSE ; }
	break ; }
      case vstatNBFR:
      { break ; }
      case vstatSB:
      { if (dy_sys->vlb[j] > -dy_tols->inf)
	{ rayDown = FALSE ; }
	if (dy_sys->vub[j] < dy_tols->inf)
	{ rayUp = FALSE ; }
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (FALSE) ; } } }
  if (rayUp == FALSE && rayDown == FALSE)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 4)
    { if (flgon(statj,vstatBASIC|vstatNBFX))
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; status %s; no ray.",
		    j,cbarj,dy_prtvstat(statj)) ; }
      else
      if (dir == 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; status %s, ub = %g; no ray.",
		    dy_prtvstat(statj),dy_sys->vub[j]) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; status %s, lb = %g; no ray.",
		    dy_prtvstat(statj),dy_sys->vlb[j]) ; } }
#   endif
    return (TRUE) ; }
/*
  We'll have to work for it. Retrieve and ftran column a<j>.
*/
  abarj = NULL ;
  if (consys_getcol_ex(dy_sys,j,&abarj) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,
	   "column",consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
    if (abarj != NULL) FREE(abarj) ;
    return (FALSE) ; }
  dy_ftran(abarj,FALSE) ;
  m = dy_sys->concnt ;
/*
  Separate the loops for an up ray and a down ray for efficiency and clarity.
  Only one case will apply.

  Open a loop to step through the basic variables. As soon as we find a limit
  on delta, we're done. Clearly, if abar<kj> = 0, or x<k> is free, no limit
  is implied. Otherwise, if increasing x<j> moves x<k> to a finite bound,
  we're bounded.  For x<j> increasing, if abar<kj> > 0, we're moving x<k>
  toward lb<k>, and if abar<kj> < 0, we're moving x<k> toward ub<k>.
*/
  if (rayUp == TRUE)
  { for (kpos = 1 ; kpos <= m && rayUp == TRUE ; kpos++)
    { abarkj = abarj[kpos] ;
      if (withintol(abarkj,0,dy_tols->zero)) continue ;
      k = dy_basis[kpos] ;
      statk = dy_status[k] ;
      if (flgon(statk,vstatBFR)) continue ;
      if ((abarkj > 0.0 && dy_sys->vlb[k] > -dy_tols->inf) ||
	  (abarkj < 0.0 && dy_sys->vub[k] < dy_tols->inf))
      { rayUp = FALSE ;
	break ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 4)
    { if (rayUp == FALSE)
      { kpos-- ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; basis pos'n %d: abar<%d,%d> = %g; %s (%d) ",
		    kpos,k,j,abarkj,consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
	if (abarkj < 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"ub = %g ; no ray up.",
		      dy_sys->vub[k]) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"lb = %g ; no ray up.",
		      dy_sys->vlb[k]) ; } }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; confirmed.") ; } }
#   endif
  }
/*
  The same, for a down ray. Here, abar<kj> > 0 will move x<k> toward ub<k> and
  abar<kj> < 0 will move x<k> toward lb<k>.
*/
  if (rayDown == TRUE)
  { for (kpos = 1 ; kpos <= m && rayDown == TRUE ; kpos++)
    { abarkj = abarj[kpos] ;
      if (withintol(abarkj,0,dy_tols->zero)) continue ;
      k = dy_basis[kpos] ;
      statk = dy_status[k] ;
      if (flgon(statk,vstatBFR)) continue ;
      if ((abarkj > 0.0 && dy_sys->vub[k] < dy_tols->inf) ||
	  (abarkj < 0.0 && dy_sys->vlb[k] > -dy_tols->inf))
      { rayDown = FALSE ;
	break ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 4)
    { if (rayDown == FALSE)
      { kpos-- ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; basis pos'n %d: abar<%d,%d> = %g; %s (%d) ",
		    kpos,k,j,abarkj,consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
	if (abarkj < 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"lb = %g ; no ray down.",
		      dy_sys->vlb[k]) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"ub = %g ; no ray down.",
		      dy_sys->vub[k]) ; } }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; confirmed.") ; } }
#   endif
  }

# ifndef DYLP_NDEBUG
  if ((rayUp == TRUE || rayDown == TRUE) && dy_opts->print.rays >= 6)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    active ray %s (%d)\n      non-zeros:",
		consys_nme(dy_sys,'v',j,FALSE,NULL),j) ;
    v = 0 ;
    for (kpos = 1 ; kpos <= m ; kpos++)
    { abarkj = abarj[kpos] ;
      if (withintol(abarkj,0,dy_tols->zero)) continue ;
      k = dy_basis[kpos] ;
      dyio_outfmt(dy_logchn,dy_gtxecho," (%s (%d) %g)",
		  consys_nme(dy_sys,'v',k,FALSE,NULL),k,abarkj) ;
      v++ ;
      if (v%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t") ; } }
# endif
  
/*
  That's it. If this is a ray, set the direction and return abar<j> if the
  client's requested it, otherwise free it.
*/
  if (rayUp == TRUE || rayDown == TRUE)
  { *p_dir = dir ;
    if (p_abarj != NULL)
    { *p_abarj = abarj ; }
    else
    { if (abarj != NULL) FREE(abarj) ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays == 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,": yes.") ; }
#   endif
  }
  else
  { if (abarj != NULL) FREE(abarj) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays == 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,": no.") ; }
#   endif
  }
  
  return (TRUE) ; }



bool dy_primalRays (lpprob_struct *orig_lp, int *p_numRays, double ***p_rays)

/*
  This routine returns the primal rays emanating from the current basic
  solution. A call to this routine is productive only when the previous call
  to dylp returned a result of unbounded and dylp's internal data structures
  are still valid. A call when the previous simplex ended in anything other
  than optimal, infeasible, or unbounded is considered an error (in judgment,
  at the least). A call when the result of optimisation was anything other
  than unbounded will return zero rays.

  Parameters:
    orig_lp:	the lp problem structure
    p_numRays:	(i) the maximum number of rays to return
		(o) the actual number of rays returned
    p_rays:	(i) vector of (double *) or NULL; 0-based indexing
		    If supplied by client, must be capable of holding at least
		    p_numRays rays.
		    If not supplied by client, allocated if necessary; in
		    particular, not allocated unless at least one ray is
		    returned
		(o) p_numRays entries will point to rays; each ray is an
		    n-vector in original system column order.

  Returns: TRUE if no errors occurred while searching for rays; FALSE
	   otherwise.
*/

{ int m,n,i,j,m_orig,n_orig,j_orig ;
  int retval ;
  bool error ;
  double *sc_abarj ;

  consys_struct *orig_sys ;
  bool scaled ;
  const double *rscale,*cscale ;
  double Sj,dir ;

  int numCols, numRays,maxRays,rayDir,j_ray,j_orig_ray,i_orig_ray ;
  flags statj_ray ;
  double **rayCollection ;
  bool ourCollection,logical ;
  double *ray ;

  char *rtnnme = "dy_primalRays" ;

# if DYLP_PARANOIA > 0
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return (FALSE) ; }
  if (p_numRays == NULL)
  { errmsg(2,rtnnme,"&numRays") ;
    return (FALSE) ; }
  if (p_rays == NULL)
  { errmsg(2,rtnnme,"&rays") ;
    return (FALSE) ; }
# endif

/*
  Do enough setup for a valid return with no rays.
*/
  maxRays = *p_numRays ;
  if (maxRays == 0)
  { return (TRUE) ; }
  *p_numRays = 0 ;
  rayCollection = *p_rays ;
  if (rayCollection != NULL)
  { ourCollection = FALSE ; }
  else
  { ourCollection = TRUE ; }
/*
  What was the result of the last lp? If it was unbounded, we have some rays.
  If it was optimal or infeasible, by definition we have no rays and can
  return TRUE. Any other code indicates an error; return FALSE.
*/
  orig_sys = orig_lp->consys ;
  switch (orig_lp->lpret)
  { case lpUNBOUNDED:
    { break ; }
    case lpOPTIMAL:
    case lpINFEAS:
    { warn(954,rtnnme,orig_sys->nme,"primal",dy_prtlpret(orig_lp->lpret)) ;
      return (TRUE) ; }
    default:
    { errmsg(954,rtnnme,orig_sys->nme,"primal",dy_prtlpret(orig_lp->lpret)) ;
      return (FALSE) ; } }
/*
  The lp was unbounded, so we have one sure ray specified in orig_lp.  Since
  we've coopted sign to indicate ray direction, the logical for constraint i
  is coded as n+i.
*/
  n_orig = orig_sys->varcnt ;
  m_orig = orig_sys->concnt ;
  i_orig_ray = -1 ;
  j_orig_ray = (int) orig_lp->obj ;
  if (j_orig_ray < 0)
  { j_orig_ray = -j_orig_ray ;
    rayDir = -1 ; }
  else
  { rayDir = 1 ; }
  if (j_orig_ray > n_orig)
  { i_orig_ray = j_orig_ray-n_orig ;
    logical = TRUE ; }
  else
  { logical = FALSE ; }
/*
  Set up a header to hold the collection, if the client did not supply it. Also
  determine if we'll need to do unscaling, and if so acquire the scaling
  matrices.
*/
  if (ourCollection == TRUE)
  { rayCollection = (double **) MALLOC(maxRays*sizeof(double *)) ; }
  sc_abarj = NULL ;
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }
/*
  Set up to walk the columns of the active system, starting from the known
  unbounded column j and wrapping around. For each column that tests out as a
  ray, testForPrimalRay will return abar<j> in basis order in the active
  reference frame. We'll unscale it and simultaneously drop the coefficients
  into a vector in variable order in the original reference frame.
*/
  error = FALSE ;
  numRays = 0 ;
  if (logical ==  TRUE)
  { j_ray = dy_origcons[i_orig_ray] ; }
  else
  { j_ray = dy_origvars[j_orig_ray] ; }
  n = dy_sys->varcnt ;
  m = dy_sys->concnt ;
/*
  Start to walk the columns of the active system. Basic variables can't be
  rays. Neither can NBFX variables.
*/
  for (numCols = 1 ; numCols <= n ; numCols++, j_ray = (j_ray%n)+1)
  { statj_ray = dy_status[j_ray] ;
    if (flgon(statj_ray,vstatBASIC|vstatNBFX)) continue ;
/*
  The column is not obviously unqualified, so call for a thorough check.
*/
    retval = testForPrimalRay(j_ray,&rayDir,&sc_abarj) ;
    if (retval == FALSE)
    { errmsg(447,rtnnme,dy_sys->nme,
	     consys_nme(dy_sys,'v',j_ray,FALSE,NULL),j_ray,"primal") ;
      error = TRUE ;
      break ; }
    if (rayDir == 0) continue ;
/*
  We have a ray. We need to unscale it and translate it from active system
  basis order to original system variable order. Begin by allocating a vector
  to hold the ray.

  In terms of scaling, we have sc_abar<j> = inv(S)inv(B)a<j>S<j> and we need
  to remove the leading and trailing column scaling, if present.

  Getting the ray pointed in the right direction takes some work.  We're
  testing x<B> = inv(B)b - abar<j>x<j>, so ray<j> = -abar<j>. Then, if x<j>
  is decreasing, that's another factor of -1 (encoded in rayDir). Finally, if
  the constraint in question is a >= constraint in the original system,
  there's another factor of -1 folded into the row scaling (but we can't use
  that directly; we want only the sign).

  In terms of change of reference frame, we're moving from active system
  basis order to original system variable order. The translation is basis
  position to basic variable (active) to variable (original), i -> j ->
  j_orig. Drop logicals, as they're not present in the original system.

  Last, but not least, x<j> is itself moving and must be part of the ray. Add
  a coefficient of -1.0, adjusted for direction, directly in the original
  reference frame.
*/
    ray = CALLOC((n_orig+1),sizeof(double)) ;
    rayCollection[numRays] = ray ;
    numRays++ ;
    dir = -1.0*rayDir ;
    if (logical == TRUE)
    { i_orig_ray = dy_actcons[j_ray] ;
      if (orig_sys->ctyp[i_orig_ray] == contypGE)
      { dir = -dir ; } }
    else
    { j_orig_ray = dy_actvars[j_ray] ; }

    if (scaled == TRUE)
    { if (logical == TRUE)
      { Sj = rscale[i_orig_ray]*dir ; }
      else
      { Sj = 1/cscale[j_orig_ray]*dir ; }
      for (i = 1 ; i <= m ; i++)
      { if (sc_abarj[i] == 0) continue ;
	j = dy_basis[i] ;
	if (j <= dy_sys->concnt) continue ;
        j_orig = dy_actvars[j] ;
	ray[j_orig] = Sj*sc_abarj[i]*cscale[j_orig] ;
	setcleanzero(ray[j_orig],dy_tols->zero) ; } }
    else
    { for (i = 1 ; i <= m ; i++)
      { if (sc_abarj[i] == 0) continue ;
	j = dy_basis[i] ;
	if (j <= dy_sys->concnt) continue ;
	j_orig = dy_actvars[j] ;
	ray[j_orig] = sc_abarj[i]*dir ;
	setcleanzero(ray[j_orig],dy_tols->zero) ; } }
    if (sc_abarj != NULL)
    { FREE(sc_abarj) ;
      sc_abarj = NULL ; }
    if (logical == FALSE)
    { ray[j_orig_ray] = -1.0*dir ; }

#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 5)
    { if (logical == TRUE)
      { j_orig = orig_sys->varcnt+i_orig_ray ; }
      else
      { j_orig = j_orig_ray ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n    ray<%d>: %s (%d)\n      non-zeros:",
	      numRays,consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig) ;
      i = 0 ;
      for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
      { if (ray[j_orig] != 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (%s (%d) %g)",
		      consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig,
		      ray[j_orig]) ; }
	i++ ;
	if (i%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t") ; } }
#   endif

    if (numRays >= maxRays) break ; }
/*
  End of scanning loop. We've either found the requested number of rays,
  scanned all columns, or encountered an error. Time for cleanup. Free the
  vector used for abar<j>. If we have an error, free the ray collection.
*/
  if (sc_abarj != NULL) FREE(sc_abarj) ;
  if (error == TRUE)
  { if (rayCollection != NULL)
    { for (i = 0 ; i < numRays ; i++)
      { if (rayCollection[i] != NULL)
	{ FREE(rayCollection[i]) ;
	  rayCollection[i] = NULL ; } }
      if (ourCollection == TRUE) FREE(rayCollection) ; }
    return (FALSE) ; }
/*
  That's it --- finish up and return.
*/
  *p_rays = rayCollection ;
  *p_numRays = numRays ;

  return (TRUE) ; }



static void testForDualRay (int i, int *p_dir, double **p_abari)

/*
  This routine evaluates an active row abar<i> = e<i>(inv(B)N) to determine
  if it constitutes a dual ray.  As explained at the head of the file and in
  dy_dualRays, there's a nonzero possibility this is a pathological problem
  that's primal and dual infeasible. Still, we can assume that dylp has
  activated all variables that might help.

  The algorithm is a much-simplified version of the algorithm that selects the
  leaving dual variable (entering primal variable) for dual simplex.

    * Instead of searching for the smallest limit on the change in y<i>, we
      just want to know if there's no limit. Hence we can abandon the
      evaluation as soon as any basic dual variable limits the change in y<i>.

    * There's no need to worry about degeneracy, numerical stability and all
      the other little considerations that go into a selecting a good pivot.

  NOTE: This routine works in the active system reference frame. The vector
	returned is abar<i>. This is not quite a ray: it needs a coefficient
	(-1.0) for y<i> and it must be negated to become the corresponding
	ray. It makes more sense to do this in the client, once we know how
	the ray is to be presented. See, for example, dy_dualRays, which
	transforms the ray for use by the outside world.

  Parameters:
    i:		Index of the row to be evaluated
    p_dir:	(o) Direction of motion,
		  1: ray in the direction abar<i> (y<i> increasing)
		  0: not a ray
		 -1: ray in the direction -abar<i> (y<i> decreasing)
    p_abari:	(o) if non-NULL, used to return abar<i> as an expanded vector
		    in column order iff abar<i> is a dual ray, otherwise
		    *p_abari is set to NULL.

  Returns: undefined
*/

{ int bvi,j,k,m,n,dir ;
  flags stati,statk ;
  double abarik ;
  bool rayUp,rayDown ;

  double *betai,*abari ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.rays >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n      Testing if row %s (%d) is a dual ray",
	        consys_nme(dy_sys,'c',i,FALSE,NULL),i) ; }
# endif

/*
  Initial setup assumes no ray.
*/
  betai = NULL ;
  abari = NULL ;
  if (p_abari != NULL) *p_abari = NULL ;
  *p_dir = 0 ;
  rayUp = FALSE ;
  rayDown = FALSE ;
  dir = 0 ;
/*
  Start by checking the dual reduced cost, bbar<i> = e<i>inv(B)b. (a.k.a. the
  value of the basic variable in pos'n i).  If bbar<i> = 0, we can change
  y<i> all we like, but we won't go unbounded.  If bbar<i> is not zero, then
  any ray must improve the objective. If you think of the dual as minimising,
  then for bbar<i> > 0, y<i> must enter decreasing. For bbar<i> < 0, y<i>
  must enter increasing. To avoid the awkward questions that arise once you
  start to think hard about this while running dual simplex on the primal
  structures, we'll just go with the status of the basic variable:
    * If it's BUUB, then the dual will enter decreasing.
    * If it's BLLB, then the dual will enter increasing.
    * Otherwise, the dual reduced cost is 0 or unfavourable, and we don't have
      a ray.
*/
  bvi = dy_basis[i] ;
  stati = getflg(dy_status[bvi],vstatSTATUS) ;
  switch (stati)
  { case vstatBUUB:
    { rayDown = TRUE ;
      dir = -1 ;
      break ; }
    case vstatBLLB:
    { rayUp = TRUE ;
      dir = 1 ;
      break ; }
    default:
    { break ; } }
/*
  The action for no ray is simply to return TRUE. This is obscured by all the
  printing in the debug version.
*/
# ifdef DYLP_NDEBUG
  if (dir == 0)
  { return (TRUE) ; }
# else
  if (dir == 0)
  { if (dy_opts->print.rays >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,".\n\tbasic var %s (%d) %s; no ray.",
		  consys_nme(dy_sys,'v',bvi,FALSE,NULL),bvi,
		  dy_prtvstat(stati)) ; }
    return ; }
  else
  if (dy_opts->print.rays >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,".\n\tbasic var %s (%d) %s allows %s ray",
		consys_nme(dy_sys,'v',bvi,FALSE,NULL),bvi,dy_prtvstat(stati),
		((dir < 0)?"down":"up")) ; }
# endif
/*
  We'll have to work for it. We can use dy_btran to retrieve beta<i> as
  e<i>inv(B). It's annoying, but allocate abari too, even though we may throw
  it away, 'cause otherwise we'll have no place to put coefficients while
  we're working.
*/
  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
  betai = CALLOC((m+1),sizeof(double)) ;
  abari = CALLOC((n+1),sizeof(double)) ;
  betai[i] = 1.0 ;
  dy_btran(betai) ;
/*
  Separate the loops for an up ray and a down ray for efficiency and clarity.
  Only one case will apply.

  Walk the columns (skipping basic variables) looking for coefficients that
  will block motion.  We're looking at y<k> = cbar<k> - abar<ik>y<i> and
  asking whether we can drive y<k> to 0.  If we make it through all the
  columns without turning up a blocking coefficient, we have a ray.

  For an up ray (y<i> increasing), the blocking conditions are:
    * x<k> NBLB (hence cbar<k> >= 0) and abarik < 0
    * x<k> NBUB (hence cbar<k> <= 0) and abarik > 0

  If we happen to run across an NBFR variable with abarik != 0, we're
  immediately done. Dual feasibility is possible only with cbar<k> = 0.
  NBFX variables, on the other hand, can never block --- we can regard them as
  NBLB or NBUB as needed.
*/
  if (rayUp == TRUE)
  { for (k = 1 ; k <= n ; k++)
    { statk = getflg(dy_status[k],vstatSTATUS) ;
      if (flgon(statk,vstatBASIC|vstatNBFX)) continue ;
      abarik = consys_dotcol(dy_sys,k,betai) ;
      if (withintol(abarik,0,dy_tols->zero)) continue ;
      if ((flgon(statk,vstatNBLB|vstatNBFR) && abarik < 0) ||
	  (flgon(statk,vstatNBUB|vstatNBFR) && abarik > 0))
      { rayUp = FALSE ;
	break ; }
      abari[k] = abarik ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 4)
    { if (rayUp == FALSE)
      { k-- ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; %s %s (%d): abar<%d,%d> = %g; %s (%d) ",
		    dy_prtvstat(statk),consys_nme(dy_sys,'v',k,FALSE,NULL),k,
		    i,k,abarik) ;
	if (flgon(statk,vstatNBFR))
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"; no ray.") ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"; no ray up.") ; } }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; confirmed.") ; } }
#   endif
  }
/*
  The same, for a down ray.
*/
  if (rayDown == TRUE)
  { for (k = 1 ; k <= n ; k++)
    { statk = getflg(dy_status[k],vstatSTATUS) ;
      if (flgon(statk,vstatBASIC|vstatNBFX)) continue ;
      abarik = consys_dotcol(dy_sys,k,betai) ;
      if (withintol(abarik,0,dy_tols->zero)) continue ;
      if ((flgon(statk,vstatNBLB|vstatNBFR) && abarik > 0) ||
	  (flgon(statk,vstatNBUB|vstatNBFR) && abarik < 0))
      { rayDown = FALSE ;
	break ; }
      abari[k] = abarik ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 4)
    { if (rayDown == FALSE)
      { k-- ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; %s %s (%d): abar<%d,%d> = %g; %s (%d) ",
		    dy_prtvstat(statk),consys_nme(dy_sys,'v',k,FALSE,NULL),k,
		    i,k,abarik) ;
	if (flgon(statk,vstatNBFR))
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"; no ray.") ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"; no ray down.") ; } }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"; confirmed.") ; } }
#   endif
  }
  if (betai != NULL) FREE(betai) ;

# ifndef DYLP_NDEBUG
  if ((rayUp == TRUE || rayDown == TRUE) && dy_opts->print.rays >= 6)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    active ray %s (%d)\n      non-zeros:",
		consys_nme(dy_sys,'c',i,FALSE,NULL),i) ;
    k = 0 ;
    for (j = 1 ; j <= n ; j++)
    { abarik = abari[j] ;
      if (withintol(abarik,0,dy_tols->zero)) continue ;
      dyio_outfmt(dy_logchn,dy_gtxecho," (%s (%d) %g)",
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j,abarik) ;
      k++ ;
      if (k%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t") ; } }
# endif
  
/*
  That's it. If this is a ray, set the direction and return abar<i> if the
  client's requested it, otherwise free it.
*/
  if (rayUp == TRUE || rayDown == TRUE)
  { *p_dir = dir ;
    if (p_abari != NULL)
    { *p_abari = abari ; }
    else
    { if (abari != NULL) FREE(abari) ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays == 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,": yes.") ; }
#   endif
  }
  else
  { if (abari != NULL) FREE(abari) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays == 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,": no.") ; }
#   endif
  }
  
  return ; }


bool dy_dualRays (lpprob_struct *orig_lp, int *p_numRays, double ***p_rays)

/*
  This routine returns the dual rays emanating from the current basic
  solution. A call to this routine can be productive only when the previous
  call to dylp returned a result of (primal) infeasible(*) and dylp's
  internal data structures are still valid. A call when the previous simplex
  ended in anything other than optimal, infeasible, or unbounded is
  considered an error (in judgment, at the least). A call when the result of
  optimisation was anything other than infeasible will return zero rays.

  The ray returned is an m-vector in the original system frame of reference,
  with the coefficients listed in row order. We're keeping our head firmly in
  the sand with respect to the duals that would be associated with the upper
  and lower bounds on variables, were those constraints to be made explicit.
  Nor are we interested in basic dual logicals.

  (*) Not to belabour the point, but it's possible for a problem to be dual
      and primal infeasible. In this case, you'll get no rays.

  Parameters:
    orig_lp:	the lp problem structure
    p_numRays:	(i) the maximum number of rays to return
		(o) the actual number of rays returned
    p_rays:	(i) vector of (double *) or NULL; 0-based indexing
		    If supplied by client, must be capable of holding at least
		    p_numRays rays.
		    If not supplied by client, allocated if necessary; in
		    particular, not allocated unless at least one ray is
		    returned
		(o) p_numRays entries will point to rays; each ray is an
		    m-vector in original system row order.

  Returns: TRUE if no errors occurred while searching for rays; FALSE
	   otherwise.
*/

{ int m,n,i,m_orig,n_orig,i_orig,j_orig ;
  bool error ;
  double *sc_abari ;

  consys_struct *orig_sys ;
  bool scaled ;
  const double *rscale,*cscale ;
  double Si,dir ;

  int numRows, numRays,maxRays,rayDir,i_ray,i_orig_ray,bv_ray,bv_orig_ray ;
  flags statbv_ray ;
  double **rayCollection ;
  bool ourCollection,logical ;
  double *ray ;

  char *rtnnme = "dy_dualRays" ;

# if DYLP_PARANOIA > 0
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return (FALSE) ; }
  if (p_numRays == NULL)
  { errmsg(2,rtnnme,"&numRays") ;
    return (FALSE) ; }
  if (p_rays == NULL)
  { errmsg(2,rtnnme,"&rays") ;
    return (FALSE) ; }
# endif

/*
  Do enough setup for a valid return with no rays.
*/
  maxRays = *p_numRays ;
  if (maxRays == 0)
  { return (TRUE) ; }
  *p_numRays = 0 ;
  rayCollection = *p_rays ;
  if (rayCollection != NULL)
  { ourCollection = FALSE ; }
  else
  { ourCollection = TRUE ; }
/*
  What was the result of the last lp? If it was infeasible, we probably have
  some rays. If it was optimal or unbounded, by definition we have no rays
  and can return TRUE. Any other code indicates an error; return FALSE.
*/
  orig_sys = orig_lp->consys ;
  switch (orig_lp->lpret)
  { case lpINFEAS:
    { break ; }
    case lpOPTIMAL:
    case lpUNBOUNDED:
    { warn(954,rtnnme,orig_sys->nme,"dual",dy_prtlpret(orig_lp->lpret)) ;
      return (TRUE) ; }
    default:
    { errmsg(954,rtnnme,orig_sys->nme,"dual",dy_prtlpret(orig_lp->lpret)) ;
      return (FALSE) ; } }
/*
  The lp was infeasible, so with high probability we'll have a ray.

  Set up a header to hold the collection, if the client did not supply it.
  Also determine if we'll need to do unscaling, and if so acquire the scaling
  matrices.
*/
  if (ourCollection == TRUE)
  { rayCollection = (double **) MALLOC(maxRays*sizeof(double *)) ; }
  sc_abari = NULL ;
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }
/*
  Set up to walk the rows of the active system.  For each row that tests out
  as a ray, testForDualRay will return abar<i> in column order in the active
  reference frame. We'll unscale it and simultaneously drop the coefficients
  into a vector in row order in the original reference frame.
*/
  n_orig = orig_sys->varcnt ;
  m_orig = orig_sys->concnt ;
  n = dy_sys->varcnt ;
  m = dy_sys->concnt ;
  error = FALSE ;
  numRays = 0 ;
  rayDir = 0 ;
/*
  Start to walk the rows of the active system. Unless the row would be a
  candidate in a dual pivot, there's clearly no ray. Test for BLLB or BUUB
  status for the associated basic variable.
*/
  for (numRows = 1, i_ray = 1 ; numRows <= m ; numRows++, i_ray = (i_ray%m)+1)
  { bv_ray = dy_basis[i_ray] ;
    statbv_ray = dy_status[bv_ray] ;
    if (!flgon(statbv_ray,vstatBLLB|vstatBUUB)) continue ;
/*
  The row is not obviously unqualified, so call for a thorough check.
*/
    testForDualRay(i_ray,&rayDir,&sc_abari) ;
    if (rayDir == 0) continue ;
/*
  We have a ray. We need to unscale it and translate it from active system
  column order to original system row order. Begin by allocating a vector
  to hold the ray.

  In terms of scaling, we have sc_abar<i> = e<i>(inv(S)inv(B)NS<j>) and we need
  to remove the leading and trailing column scaling, if present.

  Getting the ray pointed in the right direction takes some work. Playing
  fast and loose with notation, we're testing y = cbar - abar<i>y<i>, so
  ray<i> = -abar<i>. Then, if y<i> is decreasing, that's another factor of -1
  (encoded in rayDir). Finally, if the constraint in question is a >=
  constraint in the original system, there's another factor of -1 folded into
  the row scaling (but we can't use that directly; we want only the sign).

  In terms of change of reference frame, we're moving from active system
  column order to original system row order. But since we're only interested
  in the values associated with nonbasic logicals, in practice we scan only
  the first m positions of abar<i> and the translation is simply active row
  to original row, i -> i_orig. As mentioned above, we're keeping our head
  firmly in the sand when it comes to bounded architecturals.

  Last, but not least, y<i> is itself moving and, if it's associated with an
  architectural constraint, must be part of the ray. Add a coefficient of
  -1.0, adjusted for direction, directly in the original reference frame.
*/
    ray = CALLOC((m_orig+1),sizeof(double)) ;
    rayCollection[numRays] = ray ;
    numRays++ ;
    dir = -1.0*rayDir ;
    i_orig_ray = dy_actcons[i_ray] ;
    if (orig_sys->ctyp[i_orig_ray] == contypGE)
    { dir = -dir ; }

    if (bv_ray <= m)
    { logical = TRUE ;
      bv_orig_ray = dy_actcons[bv_ray] ; }
    else
    { logical = FALSE ;
      bv_orig_ray = dy_actvars[bv_ray] ; }

    if (scaled == TRUE)
    { if (logical == TRUE)
      { Si = 1/rscale[bv_orig_ray]*dir ; }
      else
      { Si = cscale[bv_orig_ray]*dir ; }
      for (i = 1 ; i <= m ; i++)
      { if (sc_abari[i] == 0) continue ;
        i_orig = dy_actcons[i] ;
	ray[i_orig] = rscale[i_orig]*sc_abari[i]*Si ;
	setcleanzero(ray[i_orig],dy_tols->zero) ; } }
    else
    { for (i = 1 ; i <= m ; i++)
      { if (sc_abari[i] == 0) continue ;
	i_orig = dy_actcons[i] ;
	ray[i_orig] = sc_abari[i]*dir ;
	setcleanzero(ray[i_orig],dy_tols->zero) ; } }
    if (sc_abari != NULL)
    { FREE(sc_abari) ;
      sc_abari = NULL ; }
    if (logical == TRUE)
    { ray[bv_orig_ray] = -1.0*dir ; }

#   ifndef DYLP_NDEBUG
    if (dy_opts->print.rays >= 5)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n    ray<%d>: %s (%d)\n      non-zeros:",numRays,
	      consys_nme(orig_sys,'c',i_orig_ray,FALSE,NULL),i_orig_ray) ;
      j_orig = 0 ;
      for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
      { if (ray[i_orig] != 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (%s (%d) %g)",
		      consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig,
		      ray[i_orig]) ; }
	j_orig++ ;
	if (j_orig%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t") ; } }
#   endif

    if (numRays >= maxRays) break ; }
/*
  End of scanning loop. We've either found the requested number of rays,
  scanned all columns, or encountered an error. Time for cleanup. Free the
  vector used for abar<j>. If we have an error, free the ray collection.
*/
  if (sc_abari != NULL) FREE(sc_abari) ;
  if (error == TRUE)
  { if (rayCollection != NULL)
    { for (i = 0 ; i < numRays ; i++)
      { if (rayCollection[i] != NULL)
	{ FREE(rayCollection[i]) ;
	  rayCollection[i] = NULL ; } }
      if (ourCollection == TRUE) FREE(rayCollection) ; }
    return (FALSE) ; }
/*
  That's it --- finish up and return.
*/
  *p_rays = rayCollection ;
  *p_numRays = numRays ;

  return (TRUE) ; }

/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2005 -- 2008 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This file contains routines which return primal and dual solutions, unscaled,
  in the original system frame of reference. For primal variables, there are
  three routines:

    * dy_rowPrimals:	the primal basic variables x<B> = inv(B)b, in basis
			(row) order
    * dy_colPrimals:	the primal architectural variables, in column order
    * dy_logPrimals:	the primal logical variables, in row order

  Because we're working with the primal system it's easy to do a clean
  separation between logical and architectural variables. It's also useful to
  have the basic variables x<B> in row order; this will in general be a mix of
  primal architectural and logical variables. dy_rowPrimals also returns a
  vector of variable indices matching x<B>.

  As a handy adjunct, there are two routines to return the status of primal
  variables:

    * dy_colStatus:	the status of the primal architectural variables, in
			column order
    * dy_logStatus:	the status of the primal logical variables, in row
			order

  For dual variables, there are two routines:

    * dy_rowDuals:	the dual variables y = c<B>inv(B) associated with the
			architectural constraints, in basis (row) order
    * dy_colDuals:	the dual variables cbar<N> = c<N> - yN associated with
			implicit bound constraints, in column order
  
  Because we're running dual simplex on the primal constraint system, we
  don't have the same clean separation into architectural and logical duals,
  nor can we easily separate out the dual basic variables. The values
  returned by rowDuals are a mix of architectural and logical duals.  The
  values returned by colDuals are commonly called the reduced costs of the
  nonbasic primal architectural variables and are a mixture of architectural
  and logical duals.

  There is another utility, dy_orig_soln, which is (at present a fairly
  specialised routine used by dylp_utils:buildsoln. It should be rewritten to
  conform to the pattern of the previous four.

  There's a routine, dy_expandxopt, which takes the primal solution generated
  by dy_orig_soln (which contains basic variables only) and expands it to a
  full solution (all architectural variables). This, too, should be pulled from
  dylp_utils and moved here (or perhaps made obsolete).
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char svnid[] UNUSED = "$Id: dy_solutions.c 269 2009-04-02 05:38:19Z lou $" ;

#ifdef DYLP_PARANOIA
extern bool dy_std_paranoia (const lpprob_struct *orig_lp,
			     const char *rtnnme) ;
#endif



void dy_colDuals (lpprob_struct *orig_lp, double **p_cbar)

/*
  Returns the unscaled vector of duals associated with architectural columns
  (aka reduced costs), in the original system frame of reference.

  These are the duals associated with implicit bound constraints. See
  dy_rowDuals for the duals associated with explicit (architectural)
  constraints. (These latter are the usual notion of dual variables.)

  The algorithm is to walk the columns of orig_sys, copying over the reduced
  cost from dy_cbar when the variable is active, otherwise calculting cbar<j>
  on the spot.

  For active variables, we have

  sc_cbar<j> = sc_c<j> - sc_c<B>sc_inv(B)sc_a<j>
	     = c<j>S<j> - c<B>S<B>inv(S<B>)inv(B)inv(R)Ra<j>S<j>
	     = c<j>S<j> - c<B>inv(B)a<j>S<j>
	     = cbar<j>S<j>

  To unscale sc_cbar<j>, we simply multiply by 1/S<j>, keeping in mind that
  if x<j> is a logical for row i, the appropriate factor is R<i>.

  For inactive variables, we calculate dot(y,a<j>) using the scaled version
  of the original system, which leaves us with the same sc_abar<j>.

  Why not use the client's original system and the vector of unscaled duals
  returned by dy_rowDuals?  That would certainly be an option. One argument
  against it is the additional work involved to get the unscaled duals. The
  other argument is that maximising the independence of the two calculations
  means that the test routine (which confirms cbar<j> = c<j> - dot(y,a<j>)
  in the external frame) is marginally more convincing.

  Parameters:
    orig_lp:	the original lp problem
    p_cbar:	(i) pointer to vector; if NULL, a vector of the appropriate
		    size will be allocated
		(o) vector of reduced costs

  Returns: undefined
*/

{ int i,j,m,n,i_orig,j_orig,m_orig,n_orig ;
  flags statj ;
  consys_struct *orig_sys ;

  double *orig_y ;
  consys_struct *scaled_orig_sys ;
  bool scaled ;

  double cbarj ;
  double *cbar ;
  const double *rscale,*cscale ;

# ifdef DYLP_PARANOIA
  char *rtnnme = "dy_colDuals" ;

  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_cbar == NULL)
  { errmsg(2,rtnnme,"cbar") ;
    return ; }
# endif
/*
  Is unscaling required? Acquire the scaling vectors and set up scaled_orig_sys
  accordingly.
*/
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ;
    scaled_orig_sys = dy_scaled_origsys() ; }
  else
  { scaled_orig_sys = NULL ; }

  orig_sys = orig_lp->consys ;
  n_orig = orig_sys->varcnt ;
  m_orig = orig_sys->concnt ;
  n = dy_sys->varcnt ;
  m = dy_sys->concnt ;
/*
  Do we need a vector?
*/
  if (*p_cbar != NULL)
  { cbar = *p_cbar ;
    memset(cbar,0,(n_orig+1)*sizeof(double)) ; }
  else
  { cbar = (double *) CALLOC((n+1),sizeof(double)) ; }
/*
  Make a vector of duals that matches orig_sys, for efficient pricing of
  inactive columns.
*/
  orig_y = (double *) CALLOC((m_orig+1),sizeof(double)) ;
  for (i = 1 ; i <= m ; i++)
  { i_orig = dy_actcons[i] ;
    orig_y[i_orig] = dy_y[i] ; }
/*
  Get on with the calculation. For an active variable, we can pull the value
  from dy_cbar. For an inactive variable, we need to calculate dot(y,a<j>).
  Then we unscale and drop the result into the proper place in the result
  vector.  Since we're starting from orig_sys, we'll never reference a column
  for a logical variable.
*/
  for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
  { if (ACTIVE_VAR(j_orig))
    { j = dy_origvars[j_orig] ;
      statj = getflg(dy_status[j],vstatSTATUS) ;
      if (flgon(statj,vstatBASIC))
      { cbarj = 0.0 ; }
      else
      { cbarj = dy_cbar[j]/cscale[j_orig] ; } }
    else
    { if (scaled == TRUE)
      { cbarj = scaled_orig_sys->obj[j_orig] ; 
	cbarj -= consys_dotcol(scaled_orig_sys,j_orig,orig_y) ;
	cbarj /= cscale[j_orig] ; }
      else
      { cbarj = orig_sys->obj[j_orig] ;
	cbarj -= consys_dotcol(orig_sys,j_orig,orig_y) ; } }
    setcleanzero(cbarj,dy_tols->cost) ;
    cbar[j_orig] = cbarj ; }
/*
  Clean up a bit and we're done.
*/
  if (orig_y != NULL) FREE(orig_y) ;
  *p_cbar = cbar ;

  return ; }



void dy_rowDuals (lpprob_struct *orig_lp, double **p_y)

/*
  This routine returns the unscaled vector of row duals, commonly referred to
  as the dual variables, c<B>inv(B). The values are unscaled and returned in a
  vector matching the original system frame of reference. Duals associated with
  inactive rows are always zero. The relevant bit of unscaling is:

  sc_y<i> = - sc_c<B>sc_inv(B)
	  = - c<B>S<B>inv(S<B>)inv(B)inv(R)
	  = - c<B>inv(B)inv(R)

  So, to recover y<i> we need to postmultiply by inv(R). The appropriate row
  factor is the one associated with the original row.

  Parameters:
    orig_lp:	the original lp problem
    p_y:	(i) vector to hold the dual variables; if NULL, a vector of
		    appropriate size will be allocated
		(o) values of the dual variables, unscaled, in the original
		    system frame of reference

  Returns: undefined
*/

{ int i,m,n,i_orig,m_orig,n_orig ;
  double yi ;
  double *y ;

  consys_struct *orig_sys ;

  bool scaled ;
  const double *rscale,*cscale ;

# ifndef DYLP_NDEBUG
  int j,v ;
# endif
# ifdef DYLP_PARANOIA
  char *rtnnme = "dy_rowDuals" ;

  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_y == NULL)
  { errmsg(2,rtnnme,"y") ;
    return ; }
# endif

/*
  Is unscaling required? Acquire the scaling vectors.
  accordingly.
*/
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }

  orig_sys = orig_lp->consys ;
  n_orig = orig_sys->varcnt ;
  m_orig = orig_sys->concnt ;
  n = dy_sys->varcnt ;
  m = dy_sys->concnt ;
/*
  Do we need a vector?
*/
  if (*p_y != NULL)
  { y = *p_y ;
    memset(y,0,(m_orig+1)*sizeof(double)) ; }
  else
  { y = (double *) CALLOC((m_orig+1),sizeof(double)) ; }
/*
  Step through the constraints of the original system. For active constraints,
  acquire and unscale the dual value.
*/
  for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
  { if (ACTIVE_CON(i_orig))
    { i = dy_origcons[i_orig] ;
      yi = dy_y[i] ;
      if (scaled == TRUE)
      { yi *= rscale[i_orig] ; }
      setcleanzero(yi,dy_tols->cost) ; }
    else
    { yi = 0.0 ; }
    y[i_orig] = yi ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.soln >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\ty =") ;
    v = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (y[i_orig] != 0)
      { if ((++v)%4 == 0)
	{ v = 0 ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; }
	i = dy_origcons[i_orig] ;
	j = dy_basis[i] ;
	dyio_outfmt(dy_logchn,dy_gtxecho," (%d %g %s %d)",
		    i_orig,y[i_orig],
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; } } }
# endif

/*
  That's it. Return the vector.
*/
  *p_y = y ;

  return ; }



void dy_colPrimals (lpprob_struct *orig_lp, double **p_x)

/*
  This routine returns the values of the primal architectural variables
  (basic and nonbasic), unscaled, in the frame of reference of the original
  system. Unscaling is straightforward. For basic variables, we have

  sc_x<B> = sc_inv(B)sc_b
	  = inv(S<B>)inv(B)inv(R)Rb
	  = inv(S<B>)(inv(B)b)

  so all that's needed to recover x<B> = inv(B)b is to multiply by S<B>.
  Upper and lower bounds on variables have the same scaling (inv(S)).

  Parameters:
    orig_lp:	the original lp problem
    p_x:	(i) vector to hold the primal architectural variables;
		    if NULL, a vector of appropriate size will be allocated
		(o) values of the primal architectural variables, unscaled,
		    in the original system frame of reference

  Returns: undefined
*/

{ int j,j_orig,n_orig ;
  double xj ;
  flags statj ;

  consys_struct *orig_sys ;
  double *x ;

  bool scaled ;
  const double *rscale,*cscale ;

  char *rtnnme = "dy_colPrimals" ;

# ifndef DYLP_NDEBUG
  int v ;
# endif

# ifdef DYLP_PARANOIA
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_x == NULL)
  { errmsg(2,rtnnme,"x") ;
    return ; }
# endif

/*
  Is unscaling required? Acquire the scaling vectors.
  accordingly.
*/
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }

  orig_sys = orig_lp->consys ;
  n_orig = orig_sys->varcnt ;
/*
  Do we need a vector?
*/
  if (*p_x != NULL)
  { x = *p_x ;
    memset(x,0,(n_orig+1)*sizeof(double)) ; }
  else
  { x = (double *) CALLOC((n_orig+1),sizeof(double)) ; }
/*
  Walk the columns of the original system. For each variable that's active
  (basic or nonbasic), we can obtain the value from dy_x and unscale. For
  each variable that's inactive, we have to do a bit of work to decode the
  status and look up the appropriate bound value.
*/
  for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
  { if (ACTIVE_VAR(j_orig))
    { j = dy_origvars[j_orig] ;
      if (scaled == TRUE)
      { xj = cscale[j_orig]*dy_x[j] ; }
      else
      { xj = dy_x[j] ; } }
    else
    { statj = (flags)(-dy_origvars[j_orig]) ;
      switch (statj)
      { case vstatNBFX:
	case vstatNBLB:
	{ xj = orig_sys->vlb[j_orig] ;
	  break ; }
	case vstatNBUB:
	{ xj = orig_sys->vub[j_orig] ;
	  break ; }
	case vstatNBFR:
	{ xj = 0 ;
	  break ; }
	default:
	{ warn(359,rtnnme,orig_sys->nme,
	       consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig,
		 dy_prtvstat(statj)) ;
	  xj = 0.0 ;
	  break ; } } }

    setcleanzero(xj,dy_tols->zero) ;
    x[j_orig] = xj ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.soln >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\tx =") ;
    v = 0 ;
    for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
    { if (x[j_orig] != 0)
      { if ((++v)%4 == 0)
	{ v = 0 ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; }
	dyio_outfmt(dy_logchn,dy_gtxecho," (%d %g %s)",
		    j_orig,x[j_orig],
		    consys_nme(orig_sys,'v',j_orig,FALSE,NULL)) ; } } }
# endif

/*
  That's it. Return the vector.
*/
  *p_x = x ;

  return ; }



void dy_rowPrimals (lpprob_struct *orig_lp, double **p_xB, int **p_indB)

/*
  This routine returns the values of the primal basic variables, unscaled, in
  row (basis) order in the frame of reference of the original system.

  Unscaling is straightforward:

  sc_x<B> = sc_inv(B)sc_b
	  = inv(S<B>)inv(B)inv(R)Rb
	  = inv(S<B>)(inv(B)b)

  so all that's needed to recover x<B> = inv(B)b is to multiply by S<B>. For
  logicals, recall that S<i> = 1/R<i>.

  By construction, the basic variable for inactive constraints is the logical
  for the constraint. Obtaining beta<i> for an inactive row and calculating
  dot(beta<i>,b) is a lot of work. Use b<i> - dot(a<i>,x) instead.

  Parameters:
    orig_lp:	the original lp problem
    p_xB:	(i) vector to hold the values of the primal basic variables;
		    if NULL, a vector of appropriate size will be allocated
		(o) values of the primal basic variables, unscaled, in the
		    original system frame of reference
    p_indB:	(i) vector to hold the indices of the primal basic variables;
		    if NULL, a vector of appropriate size will be allocated
		(o) indices of the primal basic variables, unscaled, in the
		    original system frame of reference

  Returns: undefined
*/

{ int i,j,m,i_orig,j_orig,m_orig,n_orig ;
  double xj,lhs ;

  consys_struct *orig_sys ;
  double *x,*xB ;
  int *indB ;

  bool scaled ;
  const double *rscale,*cscale ;

# ifndef DYLP_NDEBUG
  int v ;
# endif
# ifdef DYLP_PARANOIA
  char *rtnnme = "dy_rowPrimals" ;

  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_xB == NULL)
  { errmsg(2,rtnnme,"x") ;
    return ; }
  if (p_indB == NULL)
  { errmsg(2,rtnnme,"x") ;
    return ; }
# endif

/*
  Is unscaling required? Acquire the scaling vectors.

  If there are inactive constraints, we'll need the primal architecturals in
  order to calculate the value of the associated (basic) logical.
*/
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }

  orig_sys = orig_lp->consys ;
  n_orig = orig_sys->varcnt ;
  m_orig = orig_sys->concnt ;
  m = dy_sys->concnt ;

  x = NULL ;
  if (m < m_orig)
  { dy_colPrimals(orig_lp,&x) ; }

/*
  Do we need vectors? Do the necessary setup.
*/
  if (*p_xB != NULL)
  { xB = *p_xB ;
    memset(xB,0,(m_orig+1)*sizeof(double)) ; }
  else
  { xB = (double *) CALLOC((m_orig+1),sizeof(double)) ; }
  if (*p_indB != NULL)
  { indB = *p_indB ;
    memset(indB,0,(m_orig+1)*sizeof(int)) ; }
  else
  { indB = (int *) CALLOC((m_orig+1),sizeof(int)) ; }
/*
  Walk the columns of the original system. For each constraint that's active,
  we can obtain the value from dy_xbasic. For each inactive constraint, we
  need to calculate the value of the logical.

  Indices of logicals are recorded in indB as the negative of the constraint
  index.
*/
  for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
  { if (ACTIVE_CON(i_orig))
    { i = dy_origcons[i_orig] ;
      j = dy_basis[i] ;
      if (j <= m)
      { j_orig = dy_actcons[j] ; }
      else
      { j_orig = dy_actvars[j] ; }
      if (scaled == TRUE)
      { if (j <= m)
	{ xj = (1/rscale[j_orig])*dy_xbasic[i] ; }
	else
	{ xj = cscale[j_orig]*dy_xbasic[i] ; } }
      else
      { xj = dy_xbasic[i] ; }
      if (j <= m)
      { indB[i_orig] = -j_orig ; }
      else
      { indB[i_orig] = j_orig ; } }
    else
    { lhs = consys_dotrow(orig_sys,i_orig,x) ;
      xj = orig_sys->rhs[i_orig]-lhs ;
      indB[i_orig] = -i_orig ; }

    setcleanzero(xj,dy_tols->zero) ;
    xB[i_orig] = xj ; }

  if (x != NULL) FREE(x) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.soln >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\txB =") ;
    v = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if ((++v)%4 == 0)
      { v = 0 ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; }
      j_orig = indB[i_orig] ;
      if (j_orig < 0)
      { j = n_orig-j_orig ; }
      else
      { j = j_orig ; }
      dyio_outfmt(dy_logchn,dy_gtxecho," (%d %g %s %d)",
		  i_orig,xB[i_orig],
		  consys_nme(orig_sys,'v',j,FALSE,NULL),j_orig) ; } }
# endif

/*
  That's it. Return the vectors.
*/
  *p_xB = xB ;
  *p_indB = indB ;

  return ; }




void dy_logPrimals (lpprob_struct *orig_lp, double **p_logx)

/*
  This routine returns the values of the primal logical variables, unscaled,
  in the frame of reference of the original system (i.e., the value of the
  logical for constraint i is in position i of the vector). Unscaling is
  straightforward:

  sc_x<B> = sc_inv(B)sc_b
	  = inv(S<B>)inv(B)inv(R)Rb
	  = inv(S<B>)(inv(B)b)

  so all that's needed to recover x<B> = inv(B)b is to multiply by S<B>. We
  just have to remember that for a logical, S<i> = 1/R<i>. It's more work to
  get the value of the logical for an inactive constraint --- we have to
  actually calculate b - dot(a<i>,x).

  Parameters:
    orig_lp:	the original lp problem
    p_logx:	(i) vector to hold the primal logical variables;
		    if NULL, a vector of appropriate size will be allocated
		(o) values of the primal logical variables, unscaled,
		    in the original system frame of reference

  Returns: undefined
*/

{ int j,m,i_orig,m_orig ;
  double xj,lhs ;

  consys_struct *orig_sys ;
  double *logx,*x ;

  bool scaled ;
  const double *rscale,*cscale ;

# ifndef DYLP_NDEBUG
  int v,n_orig ;
# endif
# ifdef DYLP_PARANOIA
  char *rtnnme = "dy_logPrimals" ;

  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_logx == NULL)
  { errmsg(2,rtnnme,"logx") ;
    return ; }
# endif

/*
  Is unscaling required? Acquire the scaling vectors. If we have inactive
  constraints, we'll need the values of the architecturals in order to
  calculate the value of the associated logical.
*/
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }

  orig_sys = orig_lp->consys ;
  m_orig = orig_sys->concnt ;
  m = dy_sys->concnt ;

  x = NULL ;
  if (m < m_orig)
  { dy_colPrimals(orig_lp,&x) ; }
/*
  Do we need a vector?
*/
  if (*p_logx != NULL)
  { logx = *p_logx ;
    memset(logx,0,(m_orig+1)*sizeof(double)) ; }
  else
  { logx = (double *) CALLOC((m_orig+1),sizeof(double)) ; }
/*
  Walk the rows of the original system. For each constraint that's active, we
  can obtain the value of the associated logical from dy_x. For each
  constraint that's inactive, we have to actually calculate the row activity
  dot(x,a<i>) and do the arithmetic.
*/
  for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
  { if (ACTIVE_CON(i_orig))
    { j = dy_origcons[i_orig] ;
      if (scaled == TRUE)
      { xj = (1/rscale[i_orig])*dy_x[j] ; }
      else
      { xj = dy_x[j] ; } }
    else
    { lhs = consys_dotrow(orig_sys,i_orig,x) ;
      xj = orig_sys->rhs[i_orig]-lhs ; }

    setcleanzero(xj,dy_tols->zero) ;
    logx[i_orig] = xj ; }

  if (x != NULL) FREE(x) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.soln >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\tlogx =") ;
    n_orig = orig_sys->varcnt ;
    v = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (logx[i_orig] != 0)
      { if ((++v)%4 == 0)
	{ v = 0 ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; }
	dyio_outfmt(dy_logchn,dy_gtxecho," (%d %g %s)",
		    i_orig,logx[i_orig],
		    consys_nme(orig_sys,'v',n_orig+i_orig,FALSE,NULL)) ; } } }
# endif

/*
  That's it. Return the vector.
*/
  *p_logx = logx ;

  return ; }



void dy_colStatus (lpprob_struct *orig_lp, flags **p_colstat)

/*
  This routine returns the status of the primal architectural variables, in
  column order for the original system. The routine reports out the full set of
  dylp status codes.

  Parameters:
    orig_lp:	the original lp problem
    p_colstat:	(i) vector to hold the status of the primal architectural
		    variables; if NULL, a vector of appropriate size will
		    be allocated
		(o) status of the primal architectural variables, in the
		    original system frame of reference

  Returns: undefined
*/

{ int j,j_orig,n_orig ;
  flags statj ;

  consys_struct *orig_sys ;
  flags *colstat ;

# ifndef DYLP_NDEBUG
  int v ;
# endif
# ifdef DYLP_PARANOIA
  char *rtnnme = "dy_colStatus" ;

  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_colstat == NULL)
  { errmsg(2,rtnnme,"colstat") ;
    return ; }
# endif

  orig_sys = orig_lp->consys ;
  n_orig = orig_sys->varcnt ;
/*
  Do we need a vector?
*/
  if (*p_colstat != NULL)
  { colstat = *p_colstat ;
    memset(colstat,0,(n_orig+1)*sizeof(flags)) ; }
  else
  { colstat = (flags *) CALLOC((n_orig+1),sizeof(flags)) ; }
/*
  Walk the columns of the original system. For active variables, copy the
  status from dy_status. For inactive variables, we acquire it from
  dy_origvars.
*/
  for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
  { if (ACTIVE_VAR(j_orig))
    { j = dy_origvars[j_orig] ;
      statj = dy_status[j] ; }
    else
    { statj = (flags)(-dy_origvars[j_orig]) ; }
    colstat[j_orig] = statj ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.soln >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\tcolstat =") ;
    v = 0 ;
    for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
    { if ((++v)%4 == 0)
      { v = 0 ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %s)",
		  consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig,
		  dy_prtvstat(colstat[j_orig])) ; } }
# endif

/*
  That's it. Return the vector.
*/
  *p_colstat = colstat ;

  return ; }




void dy_logStatus (lpprob_struct *orig_lp, flags **p_logstat)

/*
  This routine returns the status of the primal logical variables, in row
  order for the original system. The routine reports out the full set of dylp
  status codes.

  It's actually a fair bit of work to get the status right for inactive
  constraints. Because we're reporting the full set of dylp status codes, and
  the client might be calling in a situation where the outcome was infeasible
  or unbounded, wee need to calculate the value and assign the appropriate
  status code. 

  Parameters:
    orig_lp:	the original lp problem
    p_logstat:	(i) vector to hold the status of the primal logical variables;
		    if NULL, a vector of appropriate size will be allocated
		(o) status of the primal logical variables, in the
		    original system frame of reference

  Returns: undefined
*/

{ int i,m,i_orig,m_orig ;
  flags stati ;
  double rhsi,rhslowi,lhsi,xi,lbi,ubi ;

  consys_struct *orig_sys ;
  flags *logstat ;
  double *x ;

  char *rtnnme = "dy_logStatus" ;

# ifndef DYLP_NDEBUG
  int v,n_orig ;
# endif

# ifdef DYLP_PARANOIA
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return ; }
  if (p_logstat == NULL)
  { errmsg(2,rtnnme,"logstat") ;
    return ; }
# endif

  orig_sys = orig_lp->consys ;
  m_orig = orig_sys->concnt ;
  m = dy_sys->concnt ;
/*
  If we're not playing with a full deck, we'll need the values of the
  architecturals to determine the appropriate status for the logical.
*/
  x = NULL ;
  if (m < m_orig)
  { dy_colPrimals(orig_lp,&x) ; }
/*
  Do we need a vector?
*/
  if (*p_logstat != NULL)
  { logstat = *p_logstat ;
    memset(logstat,0,(m_orig+1)*sizeof(flags)) ; }
  else
  { logstat = (flags *) CALLOC((m_orig+1),sizeof(flags)) ; }
/*
  Walk the rows of the original system. For active constraints, copy the
  status of the logical from dy_status. For inactive constraints, we need to
  actually calculate the value of the logical and assign the appropriate
  status. This is more work than you'd think, because we need to determine the
  appropriate bounds for the logical based on the constraint type, and we need
  to allow for the possibility that the problem was infeasible or unbounded and
  the logical is not within bounds. We also need to allow for the possibility
  that dylp deactivated a tight constraint with y<i> = 0. The convention for
  logicals in the original system is that all have a coefficient of 1.0. Thus
  we have bounds of (0,infty) for a slack (contypLE), (0,0) for an artificial
  (contypEQ), (-infty,0) for a surplus (contypGE), and (0,rhs-rhslow) for a
  bounded slack (contypRNG).
*/
  for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
  { if (ACTIVE_CON(i_orig))
    { i = dy_origcons[i_orig] ;
      stati = dy_status[i] ; }
    else
    { lhsi = consys_dotrow(orig_sys,i_orig,x) ;
      rhsi = orig_sys->rhs[i_orig] ;
      xi = rhsi-lhsi ;
      setcleanzero(xi,dy_tols->zero) ;
      lbi = -dy_tols->inf ;
      ubi = dy_tols->inf ;
      switch (orig_sys->ctyp[i_orig])
      { case contypLE:
	{ lbi = 0.0 ;
	  break ; }
	case contypEQ:
	{ lbi = 0.0 ;
	  ubi = 0.0 ;
	  break ; }
        case contypGE:
	{ ubi = 0.0 ;
	  break ; }
	case contypRNG:
	{ rhslowi = orig_sys->rhslow[i_orig] ;
	  lbi = 0 ;
	  ubi = rhsi-rhslowi ;
	  break ; }
	case contypNB:
	{ continue ; }
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  break ; } }
      if (belowbnd(xi,lbi))
      { stati = vstatBLLB ; }
      else
      if (atbnd(xi,lbi))
      { stati = vstatBLB ; }
      else
      if (atbnd(xi,ubi))
      { stati = vstatBUB ; }
      else
      if (abovebnd(xi,ubi))
      { stati = vstatBUUB ; }
      else
      { stati = vstatB ; } }
    logstat[i_orig] = stati ; }

  if (x != NULL) FREE(x) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.soln >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\trowstat =") ;
    n_orig = orig_sys->varcnt ;
    v = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if ((++v)%4 == 0)
      { v = 0 ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %s)",
		  consys_nme(orig_sys,'v',i_orig+n_orig,FALSE,NULL),i_orig,
		  dy_prtvstat(logstat[i_orig])) ; } }
# endif

/*
  That's it. Return the vector.
*/
  *p_logstat = logstat ;

  return ; }




void dy_orig_soln (double *x, double *y)

/*
  This routine unscales the primal and dual variable values associated with the
  rows of the active system before returning them to the user. The necessary
  unscaling is as follows:

    primal architectural:	x<j>S<j>
    primal logical:		x<i>/R<i>

    dual:			y<i>R<i>

  The vectors are indexed by basis position.

  This routine isn't really general purpose --- it's called only from
  dylp_utils:buildsoln and assumes that x and y are already populated with
  scaled values. It should get a makeover to match the interface conventions of
  the other routines in the package.

  Parameters:
    x:	basic primal variables
    y:	dual variables

  Returns: undefined.
*/

{ int i,j,i_orig,j_orig ;
  double xi,yi ;
  const double *rscale,*cscale ;

/*
  Did we scale? If not, return right off. Otherwise, acquire the scaling
  vectors.
*/
  if (dy_isscaled() == FALSE) return ;
  dy_scaling_vectors(&rscale,&cscale) ;
/*
  Since we're only dealing with duals and basic primal variables, it suffices
  to step through the constraints (equivalently, basis positions).
*/
  for (i = 1 ; i <= dy_sys->concnt ; i++)
  { i_orig = dy_actcons[i] ;
    j = dy_basis[i] ;
    xi = x[i] ;
    if (j <= dy_sys->concnt)
      xi /= rscale[i_orig] ;
    else
    { j_orig = dy_actvars[j] ;
      xi *= cscale[j_orig] ; }
    setcleanzero(xi,dy_tols->zero) ;
    x[i] = xi ;
    
    yi = y[i] ;
    yi *= rscale[i_orig] ;
    setcleanzero(yi,dy_tols->cost) ;
    y[i] = yi ; }

  return ; }



bool dy_expandxopt (lpprob_struct *lp, double **p_xopt)

/*
  This is a utility routine to load an expanded vector with the optimal
  solution to an lp relaxation. If the client supplies the vector, it's
  assumed it's large enough to hold the result.

  Note that unscaling is not required here. lp->x should have been unscaled
  when it was generated, and the client's constraint system (lp->consys) is
  not touched when dylp scales.

  Parameters:
    lp:		lpprob_struct with optimal solution attached
    p_xopt:	(i) vector to be filled in (created if null)
		(o) vector filled with optimal solution from lp

  Returns: TRUE if there's no problem translating the solution, FALSE
	   otherwise.
*/

{ int j,jpos ;
  consys_struct *consys ;
  flags *status,jstat ;
  double *xopt ;

  const char *rtnnme = "dy_expandxopt" ;

# ifdef DYLP_PARANOIA
  if (p_xopt == NULL)
  { errmsg(2,rtnnme,"&x<opt>") ;
    return (FALSE) ; }
  if (lp == NULL)
  { errmsg(2,rtnnme,"lp problem") ;
    return (FALSE) ; }
  if (lp->lpret != lpOPTIMAL)
  { errmsg(4,rtnnme,"lp return code",dy_prtlpret(lp->lpret)) ;
    return (FALSE) ; }
  if (lp->consys == NULL)
  { errmsg(2,rtnnme,"lp constraint system") ;
    return (FALSE) ; }
  if (lp->basis == NULL)
  { errmsg(2,rtnnme,"lp basis") ;
    return (FALSE) ; }
  if (lp->basis->el == NULL)
  { errmsg(2,rtnnme,"lp basis vector") ;
    return (FALSE) ; }
  if (lp->status == NULL)
  { errmsg(2,rtnnme,"lp status") ;
    return (FALSE) ; }
# endif

  consys = lp->consys ;
  status = lp->status ;
/*
  If the user didn't supply a solution vector, allocate one now.
*/
  if (*p_xopt == NULL)
  { xopt = (double *) MALLOC((consys->varcnt+1)*sizeof(double)) ; }
  else
  { xopt = *p_xopt ; }

  for (j = 1 ; j <= consys->varcnt ; j++)
  { if (((int ) status[j]) < 0)
    { jstat = vstatB ;
      jpos = -((int) status[j]) ;
      xopt[j] = lp->x[jpos] ; }
    else
    { jstat = status[j] ;
      switch (jstat)
      { case vstatNBFX:
	case vstatNBLB:
	{ xopt[j] = consys->vlb[j] ;
	  break ; }
	case vstatNBUB:
	{ xopt[j] = consys->vub[j] ;
	  break ; }
	case vstatNBFR:
	{ xopt[j] = 0 ;
	  break ; }
	default:
	{ errmsg(359,rtnnme,consys->nme,
		 consys_nme(consys,'v',j,FALSE,NULL),j,dy_prtvstat(jstat)) ;
	  if (*p_xopt == NULL) FREE(xopt) ;
	  return (FALSE) ; } } } }

  *p_xopt = xopt ;

  return (TRUE) ; }

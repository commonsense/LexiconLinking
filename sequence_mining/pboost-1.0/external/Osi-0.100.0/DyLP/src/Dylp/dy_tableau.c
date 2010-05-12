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
  This file contains routines which implement `tableau' functions. They will
  calculate:
    * a column of the basis inverse, beta<j> = inv(B)e<j>
    * a row of the basis inverse, beta<i> = e<i>inv(B)
    * a column of inv(B)A, abar<j> = inv(B)Ae<j> = inv(B)a<j>
    * a row of inv(B)A, abar<i> = e<i>inv(B)A

  The requested column or row should be given in the context of the original
  system, and the result will be returned in this context.

  Since dylp's active system is not always the full original system, we need to
  take some care. Let B be the basic partition of the active system, and N the
  inactive partition. Let G be the matrix composed of coefficients of inactive
  rows and basic columns, and let H be the matrix composed of coefficients of
  inactive rows and nonbasic columns.

  Now, if we activate the inactive rows G and declare the logical variable
  for each row to be basic, the basic component of the full system will be
  the matrix [[ B  0 ] [ G  I ]] The full basis inverse can be calculated as
  [[ inv(B)  0 ] [ -Ginv(B)  I ]]. See the typeset documentation for a decent
  presentation of all this.

  Dylp deletes its pointer to the original system when it returns --- this is
  the only safe course, because we have no control over it. Following that
  logic, the routines here require the client to pass an lpprob_struct as a
  parameter and a pointer to the unscaled original system is taken from it.
  Clearly, though, things will go badly wrong if there have been changes from
  the original system used on the last call to dylp.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char svnid[] UNUSED = "$Id: dy_tableau.c 269 2009-04-02 05:38:19Z lou $" ;



#if DYLP_PARANOIA > 0

bool dy_std_paranoia (const lpprob_struct *orig_lp, const char *rtnnme)
/*
  Some standard paranoid checks for a call from outside, collected into one
  place.

  Parameters:
    orig_lp:	lp problem structure
    rtnnme:	the client routine

  Returns: FALSE if the paranoid checks fail, TRUE otherwise.
*/

{ consys_struct *orig_sys ;

/*
  Check for presence.
*/
  if (orig_lp == NULL)
  { errmsg(2,rtnnme,"orig_lp") ;
    return (FALSE) ; }
  orig_sys = orig_lp->consys ;
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
/*
  Check for a corrupt constraint system.
*/
  if (flgon(orig_sys->opts,CONSYS_CORRUPT))
  { errmsg(115,rtnnme,orig_lp->consys->nme) ;
    return (FALSE) ; }
/*
  Check that dylp and the lpprob_struct agree on whether dylp retains valid
  data structures.
*/
  if ((flgoff(orig_lp->ctlopts,lpctlDYVALID) && dy_retained == TRUE) ||
      (flgon(orig_lp->ctlopts,lpctlDYVALID) && dy_retained == FALSE))
  { errmsg(1,rtnnme,__LINE__) ;
    return (FALSE) ; }

  return (TRUE) ; }

#endif    /* DYLP_PARANOIA */





bool dy_betaj (lpprob_struct *orig_lp, int tgt_j, double **p_betaj)
/*
  Given a basic variable x<j>, this routine returns the corresponding
  unscaled column of the basis inverse, beta<j>.

  Of course, it's not quite that simple. The client only knows about the
  original system, so tgt_j is the index of x<j> in the original system.
  We need to:
    1) Find the index j in the active system, determine the basis pos'n k,
       and retrieve the portion of beta<j> (really, column k of the basis)
       that corresponds to the scaled active system.
    2) Unscale beta<j> and translate it to the original system frame of
       reference.
    3) Calculate the remaining coefficients of beta<j> due to inactive rows.

  It is assumed that orig_sys is unscaled. It's an error if tgt_j is not a
  basic variable.

  As pointed out at the head of the file, there are two components to be
  calculated: The part of beta<j> that's drawn from the active basis B, and
  the part that's drawn from the inactive matrix G. Clearly, things will go
  wrong if the constraint system passed in through orig_lp has been modified
  and no longer matches dylp's idea of the original system.

  In particular, note that dy_origvars and dy_origcons may well be attached to
  the scaled local copy of the original system. The WILL NOT be updated by
  changes to the client's unscaled copy.

  Suppose that x<j> is basic in pos'n k, which correspondes to row i_orig.
  The approach is to calculate the partially unscaled basis column sc_beta<k>
  as inv(B)R<i_orig>e<k>, then finish unscaling with S<B> as we drop the
  coefficients into their proper positions in a vector indexed in the
  orig_sys frame of reference. Then we add the coefficients due to -G inv(B)

  Parameters:
    orig_lp:	lp problem structure
    tgt_j:	column index in original system; negative values are assumed
		to specify logicals as the negative of the index of the
		associated constraint.
    p_betaj:	(i) vector to hold beta<j>; if NULL, one will be allocated;
		if non-NULL, will be cleared to zero.
		(o) inv(B)e<j>, unscaled

  Returns: TRUE if the calculation is successful, FALSE otherwise.
*/

{ int m_orig,n_orig,i_orig,j_orig,k_orig ;
  int m,n,i,j,k,j_bpos,v ;

  bool scaled,active,logical,natural ;
  double *sc_betaj,*betaj ;
  const double *rscale, *cscale ;
  double betaij ;

  pkvec_struct *ai ;

  consys_struct *orig_sys ;

  char *rtnnme = "dy_betaj" ;

# if DYLP_PARANOIA > 0
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return (FALSE) ; }
  if (p_betaj == NULL)
  { errmsg(2,rtnnme,"betaj") ;
    return (FALSE) ; }
# endif
/*
  Always check for valid data structures.
*/
  if (flgoff(orig_lp->ctlopts,lpctlDYVALID))
  { errmsg(396,rtnnme,orig_lp->consys->nme,
	   "calculate column of basis inverse") ;
    return (FALSE) ; }

  orig_sys = orig_lp->consys ;
  m_orig = orig_sys->concnt ;
  n_orig = orig_sys->varcnt ;

  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  generating column beta<%d>",tgt_j) ; }
# endif
/*
  Determine what sort of variable we're looking at, and do some validity checks
  appropriate for the type.
*/
  j = 0 ;
  if (tgt_j < 0)
  { i_orig = -tgt_j ;
    if (i_orig > m_orig)
    { errmsg(102,rtnnme,orig_sys->nme,"row",i_orig,1,m_orig) ;
      return (FALSE) ; }
    logical = TRUE ;
    if (ACTIVE_CON(i_orig))
    { active = TRUE ;
      j = dy_origcons[i_orig] ; }
    else
    { active = FALSE ; } }
  else
  if (tgt_j > 0)
  { j_orig = tgt_j ;
    if (j_orig > n_orig)
    { errmsg(102,rtnnme,orig_sys->nme,"column",j_orig,1,n_orig) ;
      return (FALSE) ; }
    logical = FALSE ;
    if (ACTIVE_VAR(j_orig))
    { active = TRUE ;
      j = dy_origvars[j_orig] ; }
    else
    { active = FALSE ; } }
  else
  { errmsg(102,rtnnme,orig_sys->nme,"column",tgt_j,1,n_orig) ;
    return (FALSE) ; }
/*
  If the variable is active, it better be basic. For a logical, check if it's
  in the natural basis position. Note that an architectural will
  automatically fail the `natural position' test.

  Inactive architecturals are by definition nonbasic, hence an error here.
  Inactive logicals are by definition basic in the natural positon and we can
  synthesize the column.
*/
  natural = FALSE ;
  if (active == TRUE)
  { j_bpos = dy_var2basis[j] ;
    if (j_bpos == 0)
    { errmsg(951,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	     "calculate column of basis inverse") ;
    return (FALSE) ; }
    i_orig = dy_actcons[j_bpos] ;
    if (j == j_bpos)
    { natural = TRUE ; } }
  else
  { if (logical == TRUE)
    { j_bpos = i_orig ;
      natural = TRUE ; }
    else
    { errmsg(950,rtnnme,"architectural variable",
	     consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig,
	     "calculate column of basis inverse") ;
      return (FALSE) ; } }
/*
  Special case: The basis inverse column for a logical is simply a unit
  vector with 1.0 in the appropriate position, *if* the logical is in it's
  `natural' position as the basic variable for the associated constraint.
  This holds whether the logical is active or inactive. After the analysis
  above, i_orig holds the correct position.
*/
  if (logical == TRUE && natural == TRUE)
  { if (*p_betaj == NULL)
    { *p_betaj = (double *) CALLOC((m_orig+1),sizeof(double)) ; }
    else
    { memset(*p_betaj,0.0,((size_t) (m_orig+1)*sizeof(double))) ; }
    (*p_betaj)[i_orig] = 1.0 ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.tableau >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", logical for ") ;
      if (active == FALSE)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"inactive ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,"constraint %s (%d)",
		  consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig) ;
      if (active == TRUE)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", basis pos'n %d, constraint %s (%d)",j_bpos,
		    consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig) ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,".") ;
      if (dy_opts->print.tableau >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  non-zeros: (%d, %g)",
		    i_orig,(*p_betaj)[i_orig]) ; } }
#   endif
    return (TRUE) ; }
/*
  We have an architectural variable or an unnatural logical. In either case,
  the variable is active with index j and basic in basis position j_bpos.
  i_orig is the index of the corresponding row in orig_sys.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { if (logical == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  ", logical for active constraint %s (%d)",
		  consys_nme(dy_sys,'c',j,FALSE,NULL),j) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,", architectural %s (%d)",
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,
		", basis pos'n %d, constraint %s (%d).",j_bpos,
		consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig) ; }
# endif
/*
  Set up and retrieve the portion of beta<j> corresponding to the active
  matrix. The actual unscaling looks like inv(B) = S<B> sc_inv(B) R, and then
  we're extracting the column at pos'n j_bpos (which corresponds to some row
  i_orig). It's convenient to fold R<i_orig> into the unit vector and get
  half of the unscaling done as we extract beta<j>.
*/
  sc_betaj = (double *) CALLOC((m+1),sizeof(double)) ;
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ;
    sc_betaj[j_bpos] = rscale[i_orig] ; }
  else
  { sc_betaj[j_bpos] = 1.0 ; }
  dy_ftran(sc_betaj,FALSE) ;

# ifndef DYLP_NDEBUG
  /* Still in dy_sys reference frame. */
  if (dy_opts->print.tableau >= 6)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  dy_sys nonzeros:") ;
    k = 0 ;
    for (i = 1 ; i <= m ; i++)
    { if (sc_betaj[i] != 0)
      { j = dy_basis[i] ;
	dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d",
		    consys_nme(dy_sys,'c',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," %s %d %g)",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,sc_betaj[i]) ;
	k++ ;
	if (k%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t  ") ; } } }
# endif

/*
  Change reference frame, and complete the unscaling, if necessary. Recall that
  for the logical for row i, the column scaling factor is 1/R[i]. So that
  we're not testing for scaling in the loop body, replicate the loop.
*/
  if (*p_betaj == NULL)
  { betaj = (double *) CALLOC((m_orig+1),sizeof(double)) ;
    *p_betaj = betaj ; }
  else
  { memset(*p_betaj,0.0,((size_t) (m_orig+1)*sizeof(double))) ;
    betaj = *p_betaj ; }
  if (scaled == TRUE)
  { for (i = 1 ; i <= m ; i++)
    { i_orig = dy_actcons[i] ;
      k = dy_basis[i] ;
      if (k <= m)
      { k_orig = dy_actcons[k] ;
	betaj[i_orig] = sc_betaj[i]/rscale[k_orig] ;
	setcleanzero(betaj[i_orig],dy_tols->zero) ; }
      else
      { j_orig = dy_actvars[k] ;
	betaj[i_orig] = sc_betaj[i]*cscale[j_orig] ; } } }
  else
  { for (i = 1 ; i <= m ; i++)
    { i_orig = dy_actcons[i] ;
      betaj[i_orig] = sc_betaj[i] ;  } }

# ifndef DYLP_NDEBUG
  /* Now in orig_sys frame of reference. */
  if (dy_opts->print.tableau >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  active nonzeros:") ;
    k = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (betaj[i_orig] != 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %g)",
		    consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig,
		    betaj[i_orig]) ;
	k++ ;
	if (k%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t  ") ; } } }
# endif

  if (sc_betaj != NULL) FREE(sc_betaj) ;
/*
  Ok, that was the easy part. Now we need to fill in the portions of beta<j>
  contributed by inactive constraints --- the G submatrix. In some respects
  this is straightforward. We have the unscaled portion of beta<j>
  contributed by the active system, and we need to calculate -G beta<j>. Note
  that we are now working completely in the original frame of reference,
  except for a quick excursion to determine if a variable is basic in the
  active system.  Of course, if there are no loadable constraints, we can
  skip all this.

  The trek between reference frames is arduous. Given a<g,j_orig>, to find the
  appropriate row of beta<j>, we do j_orig -> j -> j_bpos -> k_orig. In
  words, column in original system to column in active system to basis
  position (row) in active system to row in original system, which is the
  element we want in beta<j>.
*/
  if (dy_lp->sys.cons.loadable > 0)
  { ai = pkvec_new(orig_sys->maxrowlen) ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (ACTIVE_CON(i_orig)) continue ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.tableau >= 5)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    Processing inactive row %s (%d)",
		    consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig) ; }
#     endif
      if (consys_getrow_pk(orig_sys,i_orig,&ai) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,"row",
	       consys_nme(orig_sys,'c',i_orig,TRUE,NULL),i_orig) ;
	if (ai != NULL) pkvec_free(ai) ;
	if (betaj != NULL) FREE(betaj) ;
	return (FALSE) ; }
      betaij = 0 ;
      for (v = 0 ; v < ai->cnt ; v++)
      { j_orig = ai->coeffs[v].ndx ;
	if (INACTIVE_VAR(j_orig)) continue ;
	j = dy_origvars[j_orig] ;
	j_bpos = dy_var2basis[j] ;
	if (j_bpos > 0)
	{ k_orig = dy_actcons[j_bpos] ;
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.tableau >= 5)
	  { dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %d %g)",
			consys_nme(orig_sys,'v',j_orig,FALSE,NULL),
			j_orig,k_orig,ai->coeffs[v].val) ; }
#	  endif
	  betaij += ai->coeffs[v].val*betaj[k_orig] ; } }
      setcleanzero(betaij,dy_tols->zero) ;
      betaj[i_orig] = -betaij ; }
    if (ai != NULL) pkvec_free(ai) ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  nonzeros:") ;
    k = 0 ;
    for (i = 1 ; i <= m_orig ; i++)
    { if (betaj[i] != 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (%d, %g)",i,betaj[i]) ;
	k++ ;
	if (k%5 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; } } }
# endif

/*
  That should do it.
*/
  return (TRUE) ; }




bool dy_abarj (lpprob_struct *orig_lp, int tgt_j, double **p_abarj)

/*
  This routine returns the unscaled ftran'd column inv(B)a<j>.
  
  Of course, it's not quite that simple. The client only knows about the
  original system, so j_orig is the index of x<j> in the original system.
  We need to:
    1) Find the index j in the active system, retrieve a<j>, and calculate
       the portion of abar<j> that corresponds to the scaled active system.
       It's entirely possible that x<j_orig> is not active, in which case
       we need to do a fair bit more work to cobble up something that looks
       like a scaled active column.
    2) Unscale abar<j> and translate it to the original system frame of
       reference.
    3) Calculate the remaining coefficients of abar<j> due to inactive rows.

  It is assumed that orig_sys is unscaled.

  Check the written documentation to get a good handle on the math. The
  relevant outer unscaling is:

  x<k> architectural basic in pos'n i: abar<ij> = S<k> sc_abar<ij> (1/S<j>)

  s<k> logical basic in pos'n i:       abar<ij> = (1/R<i>) sc_abar<ij> (1/S<j>)

  To cancel a factor of inv(R) attached to the scaled basis inverse, it'll
  be convenient to apply row scaling to an unscaled column prior to doing
  the ftran.

  Finally, we'll need to calculate the coefficients of abar<j> that belong to
  inactive constraints. Recall that we could extend the basis as
  [[B 0] [G I]] with an inverse of [[inv(B) 0] [-Ginv(B) I]]. Then for a full
  column a<j> = [a<B,j> a<G,j>], the value of abar<j> will be
  [ inv(B)a<B,j>  a<G,j> - G inv(B)a<B,j> ] = [ abarj<B> a<G,j> - G abarj<B> ].

  Parameters:
    orig_lp:	lp problem structure
    tgt_j:	column index in original system; negative values are assumed
		to specify logicals as the negative of the index of the
		associated constraint.
    p_abarj:	(i) vector to hold abar<j>; if NULL, one will be allocated
		(o) inv(B)a<j>, unscaled

  Returns: TRUE if the calculation is successful, FALSE otherwise.
*/

{ int n,m,i,j,k,j_bpos,n_orig,m_orig,i_orig,j_orig,k_orig,v ;
  double *sc_abarj,*abarj ;
  const double *rscale,*cscale ;
  double Sj,abarij,agj ;
  pkvec_struct *aj_pk ;
  bool scaled,active,logical ;

  pkvec_struct *ai ;
  consys_struct *orig_sys ;

  const char *rtnnme = "dy_abarj" ;

# if DYLP_PARANOIA > 0
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return (FALSE) ; }
  if (p_abarj == NULL)
  { errmsg(2,rtnnme,"abarj") ;
    return (FALSE) ; }
# endif
/*
  Always check for valid data structures.
*/
  if (flgoff(orig_lp->ctlopts,lpctlDYVALID))
  { errmsg(396,rtnnme,orig_lp->consys->nme,
	   "calculate column of basis inverse") ;
    return (FALSE) ; }

  orig_sys = orig_lp->consys ;
  m_orig = orig_sys->concnt ;
  n_orig = orig_sys->varcnt ;

  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n  generating column abar<%d>, ",tgt_j) ; }
# endif
/*
  If we're scaled, grab the scaling vectors.
*/
  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }
/*
  Is the column active? Remember, the client can point us at logicals with a
  negative index, so we have to consider the question of active rows. An active
  logical has the same index as its associated row.
*/
  j = 0 ;
  if (tgt_j < 0)
  { i_orig = -tgt_j ;
    if (i_orig > m_orig)
    { errmsg(102,rtnnme,orig_sys->nme,"row",i_orig,1,m_orig) ;
      return (FALSE) ; }
    logical = TRUE ;
    if (ACTIVE_CON(i_orig))
    { active = TRUE ;
      j = dy_origcons[i_orig] ; }
    else
    { active = FALSE ; } }
  else
  if (tgt_j > 0)
  { j_orig = tgt_j ;
    if (j_orig > n_orig)
    { errmsg(102,rtnnme,orig_sys->nme,"column",j_orig,1,n_orig) ;
      return (FALSE) ; }
    logical = FALSE ;
    if (ACTIVE_VAR(j_orig))
    { active = TRUE ;
      j = dy_origvars[j_orig] ; }
    else
    { active = FALSE ; } }
  else
  { errmsg(102,rtnnme,orig_sys->nme,"column",tgt_j,1,n_orig) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { if (logical == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"logical %s (%d) for ",
		  consys_nme(orig_sys,'v',n_orig+i_orig,FALSE,NULL),i_orig) ;
      if (active == FALSE)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"inactive ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,"constraint %s (%d)",
		  consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig) ; }
    else
    { if (active == FALSE)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"inactive ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,"variable %s (%d)",
		  consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig) ; } }
# endif
/*
  Special case: If the specified column represents the logical for an inactive
  constraint, the algebra says the answer is a unit vector.
*/
  if (active == FALSE && logical == TRUE)
  { if (*p_abarj == NULL)
    { abarj = (double *) CALLOC((m_orig+1),sizeof(double)) ;
      *p_abarj = abarj ; }
    else
    { abarj = *p_abarj ;
      memset(abarj,0,((size_t) (m_orig+1)*sizeof(double))) ; }
    abarj[i_orig] = 1.0 ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.tableau >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  nonzeros: (%d, %g)",
		  i_orig,abarj[i_orig]) ; }
#   endif
    return (TRUE) ; }
/*
  Set up a working vector for abar<j>.
*/
  sc_abarj = (double *) CALLOC((m+1),sizeof(double)) ;
/*
  If the column is active, we can get sc_abarj with little effort. Do half
  the outer unscaling (the final 1/S<j>) to get to the same point as the else
  case, where we start with an inactive and unscaled column.
*/
  if (active == TRUE)
  { if (logical == TRUE)
    { sc_abarj[j] = 1.0 ; }
    else
    if (consys_getcol_ex(dy_sys,j,&sc_abarj) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,"column",
	     consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
      if (sc_abarj != NULL) FREE(sc_abarj) ;
      return (FALSE) ; }
    dy_ftran(sc_abarj,FALSE) ;
    if (scaled == TRUE)
    { if (logical == TRUE)
      { Sj = rscale[i_orig] ; }
      else
      { Sj = 1/cscale[j_orig] ; }
      for (k = 1 ; k <= m ; k++)
      { sc_abarj[k] *= Sj ; } } }
/*
  An inactive column. This is an architectural (we disposed of inactive
  logicals above). We need to acquire the unscaled column from orig_sys and
  filter for the active coefficients.  Because orig_sys is unscaled, we need
  to premultiply with R (i.e., make the column look scaled) to cancel the 1/R
  attached to the basis inverse. On the other side, we don't need to multiply
  by 1/S<j> to nullify column scaling because the column is unscaled to start
  with. Once we've prepped the column, do the ftran.
*/
  else
  { aj_pk = NULL ;
    if (consys_getcol_pk(orig_sys,j_orig,&aj_pk) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,"column",
	     consys_nme(orig_sys,'v',j_orig,TRUE,NULL),j_orig) ;
      if (aj_pk != NULL) pkvec_free(aj_pk) ;
      return (FALSE) ; }
    if (scaled == TRUE)
    { for (k = 0 ; k < aj_pk->cnt ; k++)
      { i_orig = aj_pk->coeffs[k].ndx ;
	if (ACTIVE_CON(i_orig))
	{ i = dy_origcons[i_orig] ;
	  sc_abarj[i] = rscale[i_orig]*aj_pk->coeffs[k].val ; } } }
    else
    { for (k = 0 ; k < aj_pk->cnt ; k++)
      { i_orig = aj_pk->coeffs[k].ndx ;
	if (ACTIVE_CON(i_orig))
	{ i = dy_origcons[i_orig] ;
	  sc_abarj[i] = aj_pk->coeffs[k].val ; } } }
    pkvec_free(aj_pk) ;
    dy_ftran(sc_abarj,FALSE) ; }
/*
  We've reached a common point for active and inactive columns.  abar<j> is
  mostly unscaled and still in the active system frame of reference. Allocate
  a vector to hold the final values in the original system's frame of
  reference.
*/
  if (*p_abarj == NULL)
  { abarj = (double *) CALLOC((m_orig+1),sizeof(double)) ;
    *p_abarj = abarj ; }
  else
  { abarj = *p_abarj ;
    memset(abarj,0,((size_t) (m_orig+1)*sizeof(double))) ; }
/*
  Copy over the values, doing the final unscaling if needed.  This cancels a
  scaling factor (1/S) attached to the scaled basis inverse. The only trick
  here is that we need to account for logicals out of natural position and
  acquire the correct scaling factor for the logical actually occupying
  position i of the basis.
*/
  if (scaled == TRUE)
  { for (i = 1 ; i <= m ; i++)
    { if (sc_abarj[i] == 0) continue ;
      j = dy_basis[i] ;
      i_orig = dy_actcons[i] ;
      if (j <= dy_sys->concnt)
      { j_orig = dy_actcons[j] ;
	abarj[i_orig] = sc_abarj[i]/rscale[j_orig] ; }
      else
      { j_orig = dy_actvars[j] ;
	abarj[i_orig] = sc_abarj[i]*cscale[j_orig] ; }
      setcleanzero(abarj[i_orig],dy_tols->zero) ; } }
  else
  { for (i = 1 ; i <= m ; i++)
    { if (sc_abarj[i] == 0) continue ;
      i_orig = dy_actcons[i] ;
      abarj[i_orig] = sc_abarj[i] ;
      setcleanzero(abarj[i_orig],dy_tols->zero) ; } }

# ifndef DYLP_NDEBUG
  /* Now in orig_sys frame of reference. */
  if (dy_opts->print.tableau >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  active nonzeros:") ;
    k = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (abarj[i_orig] != 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %g)",
		    consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig,
		    abarj[i_orig]) ;
	k++ ;
	if (k%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t  ") ; } } }
# endif

  if (sc_abarj != NULL) FREE(sc_abarj) ;

/*
  Ok, that was the easy part. Now we need to fill in the portion of abar<j>
  contributed by inactive constraints --- the G matrix. In some respects this
  is straightforward. We have the unscaled portion of abarj<j> contributed by
  the active system, and we need to calculate a<G,j> - G abarj<B>. Note that
  we are now working completely in the original frame of reference, except
  for a quick excursion to determine if a variable is basic in the active
  system.  Of course, if there are no loadable constraints, we can skip all
  this.

  The translation here is pretty ugly. Given a<g,j_orig>, to find the
  appropriate row of abar<j>, we do j_orig -> j -> j_bpos -> k_orig. In words,
  column in original system to column in active system to basis position in
  active position (row) to row in original system, which is the element we want
  in abarj<j>. Nonbasic variables are not of interest except that we need the
  coefficient for our target column, tgt_j.
*/
  if (dy_lp->sys.cons.loadable > 0)
  { ai = pkvec_new(orig_sys->maxrowlen) ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (ACTIVE_CON(i_orig)) continue ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.tableau >= 5)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    Processing inactive row %s (%d)",
		    consys_nme(orig_sys,'c',i_orig,FALSE,NULL),i_orig) ; }
#     endif
      if (consys_getrow_pk(orig_sys,i_orig,&ai) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,"row",
	       consys_nme(orig_sys,'c',i_orig,TRUE,NULL),i_orig) ;
	if (ai != NULL) pkvec_free(ai) ;
	if (abarj != NULL) FREE(abarj) ;
	return (FALSE) ; }
      abarij = 0 ;
      agj = 0 ;
      for (v = 0 ; v < ai->cnt ; v++)
      { j_orig = ai->coeffs[v].ndx ;
	if (j_orig == tgt_j)
	{ agj = ai->coeffs[v].val ; }
	if (INACTIVE_VAR(j_orig))
	{ continue ; }
	j = dy_origvars[j_orig] ;
	j_bpos = dy_var2basis[j] ;
	if (j_bpos > 0)
	{ k_orig = dy_actcons[j_bpos] ;
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.tableau >= 5)
	  { dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %d %g)",
			consys_nme(orig_sys,'v',j_orig,FALSE,NULL),
			j_orig,k_orig,ai->coeffs[v].val) ; }
#	  endif
	  abarij += ai->coeffs[v].val*abarj[k_orig] ; } }
      abarj[i_orig] = agj-abarij ; }
    if (ai != NULL) pkvec_free(ai) ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  nonzeros:") ;
    k = 0 ;
    for (i = 1 ; i <= m_orig ; i++)
    { if (abarj[i] != 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (%d, %g)",i,abarj[i]) ;
	k++ ;
	if (k%5 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t   ") ; } } }
# endif

/*
  That should do it.
*/

  return (TRUE) ; }




bool dy_betai (lpprob_struct *orig_lp, int tgt_i, double **p_betai)
/*
  Given a row i, this routine returns the corresponding unscaled row of the
  basis inverse, beta<i>.

  Of course, it's not quite that simple. The client only knows about the
  original system, so tgt_i is the index of i in the original system.

  There are two cases:

    1) If tgt_i is active, we need to determine its position in the active
       system and extract the corresponding row of the basis inverse,
       e<i> inv(B). Then we need to translate this row into the original
       system frame of reference, padding it out with zeros.

    2) If tgt_i is inactive, we need to synthesize the row that would result
       if the constraint were activated, g<i> inv(B).  The logical for the
       constraint is used as the basic variable.  This is accomplished by
       translating g<i> into the active frame of reference, executing the
       btran, and then translating back to the original system frame of
       reference, adding padding and a coefficient for the slack.

  It is assumed that orig_sys is unscaled.

  Clearly, things will go wrong if the constraint system passed in through
  orig_lp has been modified and no longer matches dylp's idea of the original
  system.

  In particular, note that dy_origvars and dy_origcons may well be attached to
  the scaled local copy of the original system. The WILL NOT be updated by
  changes to the client's unscaled copy.

  Parameters:
    orig_lp:	lp problem structure
    tgt_i:	constraint (row) index in original system
    p_betai:	(i) vector to hold beta<i>; if NULL, one will be allocated;
		if non-NULL, will be cleared to zero.
		(o) e<i> inv(B), unscaled

  Returns: TRUE if the calculation is successful, FALSE otherwise.
*/

{ int m_orig,n_orig,i_orig,j_orig ;
  int m,n,i,j,j_bpos,v ;

  bool scaled,active ;
  double *sc_betai,*betai ;
  const double *rscale, *cscale ;
  double Sj,gij ;

  pkvec_struct *ai ;

  consys_struct *orig_sys ;

  char *rtnnme = "dy_betai" ;

# ifndef DYLP_NDEBUG
  int k_orig ;
# endif

# if DYLP_PARANOIA > 0
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return (FALSE) ; }
  if (p_betai == NULL)
  { errmsg(2,rtnnme,"betai") ;
    return (FALSE) ; }
# endif
/*
  Always check for valid data structures.
*/
  if (flgoff(orig_lp->ctlopts,lpctlDYVALID))
  { errmsg(396,rtnnme,orig_lp->consys->nme,"calculate row of basis inverse") ;
    return (FALSE) ; }
/*
  Do a bit of setup. Pull constraint system sizes for convenient use. Grab the
  scaling vectors if we're scaled.
*/
  orig_sys = orig_lp->consys ;
  m_orig = orig_sys->concnt ;
  n_orig = orig_sys->varcnt ;

  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;

  scaled = dy_isscaled() ;
  if (scaled == TRUE)
  { dy_scaling_vectors(&rscale,&cscale) ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  generating row beta<%d>,",tgt_i) ; }
# endif
/*
  What sort of constraint do we have?
*/
  if (ACTIVE_CON(tgt_i))
  { active = TRUE ;
    i = dy_origcons[tgt_i] ; }
  else
  { active = FALSE ;
    i = -1 ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { if (active == FALSE)
    { dyio_outfmt(dy_logchn,dy_gtxecho," inactive") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho," constraint %s (%d)",
		consys_nme(orig_sys,'c',tgt_i,FALSE,NULL),tgt_i) ;
    if (active == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", basis pos'n %d",i) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,".") ; }
# endif

/*
  For an active constraint, we can retrieve the row as e<i> inv(B). But, we
  really have the scaled basis inverse inv(S) inv(B) inv(R).  Hence it's
  convenient to use a vector with S<i> in place of a unit coefficient to
  cancel the leading scale factor. We have to be careful to get the right
  scale factor --- the column scale factor for the logical for constraint i
  is 1/R<i>, and logicals need not be in natural position.
*/
  sc_betai = (double *) CALLOC((m+1),sizeof(double)) ;
  if (active == TRUE)
  { if (scaled == TRUE)
    { j = dy_basis[i] ;
      if (j > m)
      { j_orig = dy_actvars[j] ;
	Sj = cscale[j_orig] ; }
      else
      { i_orig = dy_actcons[j] ;
	Sj = 1/rscale[i_orig] ; }
      sc_betai[i] = Sj ; }
    else
    { sc_betai[i] = 1.0 ; }
    dy_btran(sc_betai) ; }
/*
  For an inactive constraint, we have more work to do. We need to pull the
  row from orig_sys, apply column scaling, and drop the coefficients into the
  vector in basis order so that we can use btran.  But since this is an
  inactive constraint, we don't have to worry about logicals.
*/
  else
  { ai = NULL ;
    if (consys_getrow_pk(orig_sys,tgt_i,&ai) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,"row",
	     consys_nme(orig_sys,'c',tgt_i,FALSE,NULL),tgt_i) ;
      if (ai != NULL) pkvec_free(ai) ;
      if (sc_betai != NULL) FREE(sc_betai) ;
      return (FALSE) ; }
    if (scaled == TRUE)
    { for (v = 0 ; v < ai->cnt ; v++)
      { j_orig = ai->coeffs[v].ndx ;
	if (ACTIVE_VAR(j_orig))
	{ j = dy_origvars[j_orig] ;
	  j_bpos = dy_var2basis[j] ;
	  if (j_bpos > 0)
	  { gij = cscale[j_orig]*ai->coeffs[v].val ;
	    sc_betai[j_bpos] = -gij ; } } } }
    else
    { for (v = 0 ; v < ai->cnt ; v++)
      { j_orig = ai->coeffs[v].ndx ;
	if (ACTIVE_VAR(j_orig))
	{ j = dy_origvars[j_orig] ;
	  j_bpos = dy_var2basis[j] ;
	  if (j_bpos > 0)
	  { sc_betai[j_bpos] = -ai->coeffs[v].val ; } } } }
    if (ai != NULL)
    { pkvec_free(ai) ; }
    dy_btran(sc_betai) ; }
/*
  At this point, we have a row beta<i> which is partially unscaled and in
  basis order.  First order of business is to allocate a working array for
  the final product.
*/
  if (*p_betai == NULL)
  { betai = (double *) CALLOC((m_orig+1),sizeof(double)) ;
    *p_betai = betai ; }
  else
  { betai = *p_betai ;
    memset(betai,0,((size_t) (m_orig+1)*sizeof(double))) ; }
/*
  To complete the unscaling, we need to postmultiply by R. The array is in
  basis order, which is correct, but we need to reposition so that the row
  order matches the original system.
*/
  if (scaled == TRUE)
  { for (i = 0 ; i <= m ; i++)
    { i_orig = dy_actcons[i] ;
      betai[i_orig] = sc_betai[i]*rscale[i_orig] ; } }
  else
  { for (i = 0 ; i <= m ; i++)
    { i_orig = dy_actcons[i] ;
      betai[i_orig] = sc_betai[i] ; } }
  if (active == FALSE)
  { betai[tgt_i] = 1.0 ; }

  if (sc_betai != NULL) FREE(sc_betai) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  nonzeros:") ;
    v = 0 ;
    for (i_orig = 1 ; i_orig <= m_orig ; i_orig++)
    { if (betai[i_orig] != 0)
      { if (ACTIVE_CON(i_orig))
	{ i = dy_origcons[i_orig] ;
	  j = dy_basis[i] ;
	  if (j <= m)
	  { k_orig = dy_actcons[j] ;
	    dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %g)",
			consys_nme(orig_sys,'v',n_orig+k_orig,FALSE,NULL),
			k_orig,betai[i_orig]) ; }
	  else
	  { j_orig = dy_actvars[j] ;
	    dyio_outfmt(dy_logchn,dy_gtxecho," (%s %d %g)",
			consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig,
			betai[i_orig]) ; } }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho, " (%s %d %g)",
		      consys_nme(orig_sys,'v',n_orig+i_orig,FALSE,NULL),
		      i_orig,betai[i_orig]) ; }
	v++ ;
	if (v%3 == 0) dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t\t  ") ; } } }
# endif

/*
  That's it, we're done.
*/
  return (TRUE) ; }



bool dy_abari (lpprob_struct *orig_lp, int tgt_i, double **p_abari,
	       double **p_betai)
/*
  This routine returns the value of row i of inv(B)A = inv(B) [ B N ] in
  p_abari.

  If p_betai is non-NULL, the routine returns row i of inv(B) [ A I ] where
  I is the identity matrix of coefficients of logicals. Row i of inv(B)A is
  still returned in p_abari, and e<i> inv(B) I = beta<i> is returned in
  p_betai.

  Given the primitives we have available (ftran, btran), the best we can do
  here is extract the relevant row of the basis inverse and calculate
  dot(beta<i>,a<j>) for j in N. Fortunately, we have a handy routine to
  calculate beta<i>.

  Parameters:
    orig_lp:	lp problem structure
    tgt_i:	constraint (row) index in original system
    p_abari:	(i) vector to hold abar<i>; if NULL, one will be allocated;
		if non-NULL, will be cleared to zero.
		(o) e<i> inv(B) A,  unscaled
    p_betai:	(i) vector to hold beta<i>; if NULL, one will be allocated;
		if non-NULL, will be cleared to zero.
		(o) e<i> inv(B) I,  unscaled

  Returns: TRUE if the calculation is successful, FALSE otherwise.
*/

{ int m_orig,n_orig,j_orig ;
  int i,j,j_bpos ;

  bool active,dologicals ;

  double *betai,*abari ;

  consys_struct *orig_sys ;

  char *rtnnme = "dy_betai" ;

# if DYLP_PARANOIA > 0
  if (dy_std_paranoia(orig_lp,rtnnme) == FALSE)
  { return (FALSE) ; }
  if (p_abari == NULL)
  { errmsg(2,rtnnme,"abari") ;
    return (FALSE) ; }
# endif

  if (p_betai != NULL)
  { dologicals = TRUE ; }
/*
  Always check for valid data structures.
*/
  if (flgoff(orig_lp->ctlopts,lpctlDYVALID))
  { errmsg(396,rtnnme,orig_lp->consys->nme,"calculate row of basis inverse") ;
    return (FALSE) ; }
/*
  Do a bit of setup. Pull constraint system sizes for convenient use.
*/
  orig_sys = orig_lp->consys ;
  m_orig = orig_sys->concnt ;
  n_orig = orig_sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  generating row abar<%d>,",tgt_i) ; }
# endif

/*
  What sort of constraint do we have?
*/
  if (ACTIVE_CON(tgt_i))
  { active = TRUE ;
    i = dy_origcons[tgt_i] ; }
  else
  { active = FALSE ;
    i = -1 ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.tableau >= 1)
  { if (active == FALSE)
    { dyio_outfmt(dy_logchn,dy_gtxecho," inactive") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho," constraint %s (%d)",
		consys_nme(orig_sys,'c',tgt_i,FALSE,NULL),tgt_i) ;
    if (active == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", basis pos'n %d",i) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,".") ; }
# endif

/*
  Call dy_betai to get row beta<i> of the basis inverse.
*/
  betai = *p_betai ;
  if (dy_betai(orig_lp,tgt_i,&betai) == FALSE)
  { errmsg(952,rtnnme,orig_sys->nme,"row",tgt_i,"constraint",
	   consys_nme(orig_sys,'c',tgt_i,FALSE,NULL),tgt_i) ;
    if (betai != NULL) FREE(betai) ;
    return (FALSE) ; }
/*
  Get a vector to return abar<i>.
*/
  if (*p_abari == NULL)
  { abari = (double *) CALLOC((n_orig+1),sizeof(double)) ;
    *p_abari = abari ; }
  else
  { abari = *p_abari ;
    memset(abari,0,((size_t) (n_orig+1)*sizeof(double))) ; }
/*
  Now walk the columns of orig_sys calculating abar<ij> = dot(beta<i>,a<j>).

  We can help ourselves a bit here by recognising basic columns, which will
  resolve to 0 or 1, depending on whether the variable is basic for this row.
  Other than that, however, active or inactive is irrelevant.
*/
  for (j_orig = 1 ; j_orig <= n_orig ; j_orig++)
  { if (ACTIVE_VAR(j_orig))
    { j = dy_origvars[j_orig] ;
      j_bpos = dy_var2basis[j] ;
      if (j_bpos > 0)
      { if (j_bpos == i)
	{ abari[j_orig] = 1.0 ; }
	else
	{ abari[j_orig] = 0.0 ; }
	continue ; } }
    abari[j_orig] = consys_dotcol(orig_sys,j_orig,betai) ; }
/*
  Did the client ask for the columns corresponding to logicals? If so, hand
  back beta<i>. Otherwise, we're done with it.
*/
  if (dologicals == TRUE)
  { *p_betai = betai ; }
  else
  { if (betai != NULL) FREE(betai) ; }
/*
  That's it, we're done.
*/
  return (TRUE) ; }

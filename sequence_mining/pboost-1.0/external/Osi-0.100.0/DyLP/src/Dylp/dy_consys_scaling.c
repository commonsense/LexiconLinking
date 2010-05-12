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
  This file contains utility routines for scaling a constraint matrix A.
  There are routines to calculate scaling matrices using geometric mean and
  equilibration scaling.

  Let a<imax> = max{j} |a<ij>|, a<imin> = min{j s.t. a<ij> != 0} |a<ij>|.

  Geometric mean scaling scales a row (column) by dividing all coefficients
  by sqrt(a<imax>*a<imin>). The process iterates until the change in
  sqrt(a<max>/a<min>) for the matrix as a whole is less than the tolerance
  or exceeds the maximum allowable number of iterations.

  << Interestingly, in glpk the scaling is sqrt(a<imax>/a<imin>)! Clp uses a
     fixed iteration limit of 3, no tolerance on change from iteration to
     iteration. >>

  Equilibration scaling scales the row (column) by dividing all coefficients
  by a<imax>, so that the largest element in the row (column) is 1.

  The overall scaling algorithm is iterated geometric scaling followed by a
  final equilibration scaling. The result is a pair of scaling vectors R and
  S which can be treated as diagonal matrices to produce the scaled matrix
  scaled(A) = R*A*S.
*/



#include "dylib_errs.h"
#include "dylib_std.h"
#include "dy_consys.h"

#define CONSYS_SCALING_DEBUG 0

#if (CONSYS_SCALING_DEBUG > 0)
  extern ioid dy_logchn ;
  extern bool dy_gtxecho ;
#endif

static char sccsid[] UNUSED = "@(#)dy_consys_scaling.c	4.8	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_consys_scaling.c 269 2009-04-02 05:38:19Z lou $" ;



bool consys_evalsys (consys_struct *consys, double *p_scm, int *p_gecnt)

/*
  This routine evaluates the constraint system given as a parameter,
  determining the minimum and maximum coefficients and calculating an
  initial value for the geometric mean figure of merit.

  The maximum coefficient is defined as amax = max{i,j} |a<ij>|.
  The minimum coefficient is defined as amin = min{i,j, a<ij> != 0} |a<ij>|.

  The figure of merit is sqrt(amax/amin).

  Parameters:
    consys:	constraint system to be evaluated
    scm:	(o) sqrt(amax/amin)
    gecnt:	(o) the number of >= inequalities in the constraint system

  Returns: TRUE if the evaluation concludes without error, FALSE otherwise.
	   A FALSE return is possible only when we're paranoid.
*/

{ int i,gecnt ;
  double amax,amin,aij ;
  double *rsc,*csc ;

  rowhdr_struct *rowi ;
  coeff_struct *coeffij ;

# ifdef DYLP_PARANOIA

  char *rtnnme = "consys_evalsys" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(2,rtnnme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(2,rtnnme,"column header") ;
    return (FALSE) ; }
  if (consys->ctyp == NULL)
  { errmsg(101,rtnnme,consys->nme,"constraint type vector") ;
    return (FALSE) ; }
  if (p_scm == NULL)
  { errmsg(2,rtnnme,"scm") ;
    return (FALSE) ; }
  if (p_gecnt == NULL)
  { errmsg(2,rtnnme,"gecnt") ;
    return (FALSE) ; }
# endif

  *p_scm = quiet_nan(0) ;
  *p_gecnt = -1 ;

  rsc = consys->rowscale ;
  csc = consys->colscale ;
  amax = 0.0 ;
  amin = consys->inf ;
  gecnt = 0 ;

/*
  Open a loop and scan the rows of the constraint matrix.
*/
  for (i = 1 ; i <= consys->concnt ; i++)
  { rowi = consys->mtx.rows[i] ;
#   ifdef DYLP_PARANOIA
    if (rowi == NULL)
    { errmsg(103,rtnnme,consys->nme,"row",i) ;
      return (FALSE) ; }
    if (rowi->ndx != i)
    { errmsg(126,rtnnme,consys->nme,"row",rowi,rowi->ndx,i,rowi) ;
      return (FALSE) ; }
    if ((rowi->coeffs == NULL && rowi->len != 0) ||
	(rowi->coeffs != NULL && rowi->len == 0))
    { errmsg(134,rtnnme,consys->nme,"row",rowi->nme,i,rowi->len,
	     (rowi->coeffs == NULL)?"null":"non-null") ;
      return (FALSE) ; }
#   endif
    if (consys->ctyp[i] == contypGE)
    { gecnt++ ; }

    for (coeffij = rowi->coeffs ; coeffij != NULL ; coeffij = coeffij->rownxt)
    {
#     ifdef DYLP_PARANOIA
      if (coeffij->rowhdr != rowi)
      { errmsg(125,rtnnme,"rowhdr",coeffij,"row",rowi->nme,i) ;
	return (FALSE) ; }
      if (coeffij->colhdr == NULL)
      { errmsg(125,rtnnme,"colhdr",coeffij,"row",rowi->nme,i) ;
	return (FALSE) ; }
      if (coeffij->colhdr->ndx <= 0 || coeffij->colhdr->ndx > consys->varcnt)
      { errmsg(102,rtnnme,consys->nme,"column",coeffij->colhdr->ndx,1,
	       consys->varcnt) ;
	return (FALSE) ; }
#     endif

      aij = coeffij->val ;
      if (aij == 0.0) continue ;

      aij = fabs(aij) ;
      if (rsc != NULL) aij *= rsc[i] ;
      if (csc != NULL) aij *= csc[coeffij->colhdr->ndx] ;
      if (aij > amax) amax = aij ;
      if (aij < amin) amin = aij ; } }
/*
  Record the results and return. Allow for 0x0 systems; they happen for (more
  or less) legitimate reasons.
*/
  if (consys->concnt == 0)
  { *p_gecnt = 0 ;
    *p_scm = 1.0 ;
    consys->maxaij = 0 ;
    consys->minaij = 0 ; }
  else
  { *p_gecnt = gecnt ;
    *p_scm = sqrt(amax/amin) ;
    consys->maxaij = amax ;
    consys->minaij = amin ; }

  return (TRUE) ; }



bool consys_geomscale (consys_struct *consys,
		       double **p_rowscale, double **p_colscale)
/*
  Given the constraint matrix A (consys), this routine will calculate diagonal
  scaling matrices R (rowscale) and S (colscale) using geometric mean scaling.

  The routine assumes it is scaling the matrix R*A*S, and the parameters
  rowscale and colscale are updated with new scaling coefficients. The
  constraint system is >> not << scaled.

  This routine is light on paranoia, on the assumption that the matrix has
  just been scanned with consys_evalsys, which is heavy on paranoia.

  Parameters:
    consys:	constraint system
    p_rowscale:	(i) initial row scaling coefficients; created if null
		(o) revised row scaling coefficients
    p_colscale:	(i) initial column scaling coefficients; created if null
		(o) revised column scaling coefficients

  Returns: TRUE if no errors occurred while calculating the scaling
	   coefficients, FALSE otherwise.
*/

{ int i,j,iter ;
  double sqm,sqm_old,eps,aij,rcmax,rcmin,maxaij,minaij ;
  double *rowscale,*colscale ;
  coeff_struct *coeffij ;

# if defined(DYLP_PARANOIA) || CONSYS_SCALING_DEBUG >= 1

  char *rtnnme = "consys_geomscale" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(2,rtnnme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(2,rtnnme,"column header") ;
    return (FALSE) ; }
  if (p_rowscale == NULL)
  { errmsg(2,rtnnme,"&rowscale") ;
    return (FALSE) ; }
  if (p_colscale == NULL)
  { errmsg(2,rtnnme,"&colscale") ;
    return (FALSE) ; }
# endif

/*
  If the client didn't supply initial scaling matrices, create them now and
  initialise them to 1.
*/
  if (*p_rowscale == NULL)
  { rowscale = (double *) MALLOC((consys->concnt+1)*sizeof(double)) ;
    rowscale[0] = 0 ;
    for (i = 1 ; i <= consys->concnt ; i++) rowscale[i] = 1.0 ; }
  else
  { rowscale = *p_rowscale ; }
  if (*p_colscale == NULL)
  { colscale = (double *) MALLOC((consys->varcnt+1)*sizeof(double)) ;
    colscale[0] = 0 ;
    for (j = 1 ; j <= consys->varcnt ; j++) colscale[j] = 1.0 ; }
  else
  { colscale = *p_colscale ; }

  sqm_old = sqrt(consys->maxaij/consys->minaij) ;
  eps = 1.0 ;

/*
  Open up the outer loop to control the number of scaling iterations. For
  each scaling iteration, scale by row, then by column.
*/
  for (iter = 1 ; iter <= 20 && eps > .05  ; iter++)
  { maxaij = 0.0 ;
    minaij = consys->inf ;
    for (i = 1 ; i <= consys->concnt ; i++)
    { coeffij = consys->mtx.rows[i]->coeffs ;
      if (coeffij == NULL) continue ;
      rcmax = 0.0 ;
      rcmin = consys->inf ;
      for ( ; coeffij != NULL ; coeffij = coeffij->rownxt)
      { aij = fabs(coeffij->val) ;
	if (aij == 0) continue ;
	j = coeffij->colhdr->ndx ;
	aij *= colscale[j] ;
	if (aij > rcmax) rcmax = aij ;
	if (aij < rcmin) rcmin = aij ; }
      rowscale[i] = 1/sqrt(rcmax*rcmin) ;
      if (rowscale[i]*rcmax > maxaij) maxaij = rowscale[i]*rcmax ;
      if (rowscale[i]*rcmin < minaij) minaij = rowscale[i]*rcmin ; }
#   if (CONSYS_SCALING_DEBUG >= 1)
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  %s: iter %d: %g <= a<ij> <= %g, geom = %g",
		rtnnme,iter,minaij,maxaij,sqrt(maxaij/minaij)) ;
#   endif

    maxaij = 0.0 ;
    minaij = consys->inf ;
    for (j = 1 ; j <= consys->varcnt ; j++)
    { coeffij = consys->mtx.cols[j]->coeffs ;
      if (coeffij == NULL) continue ;
      rcmax = 0.0 ;
      rcmin = consys->inf ;
      for ( ; coeffij != NULL ; coeffij = coeffij->colnxt)
      { aij = fabs(coeffij->val) ;
	if (aij == 0) continue ;
	i = coeffij->rowhdr->ndx ;
	aij *= rowscale[i] ;
	if (aij > rcmax) rcmax = aij ;
	if (aij < rcmin) rcmin = aij ; }
      colscale[j] = 1/sqrt(rcmax*rcmin) ;
      if (colscale[j]*rcmax > maxaij) maxaij = colscale[j]*rcmax ;
      if (colscale[j]*rcmin < minaij) minaij = colscale[j]*rcmin ; }
    sqm = sqrt(maxaij/minaij) ;
    eps = (sqm_old-sqm)/sqm_old ;
#   if (CONSYS_SCALING_DEBUG >= 1)
    dyio_outfmt(dy_logchn,dy_gtxecho,
	   "\n  %s: iter %d: %g <= a<ij> <= %g, geom = %g, eps = %g",
	   rtnnme,iter,minaij,maxaij,sqm,eps) ;
#   endif
    sqm_old = sqm ;
  }

  consys->maxaij = maxaij ;
  consys->minaij = minaij ;

  *p_rowscale = rowscale ;
  *p_colscale = colscale ;

  return (TRUE) ; }



bool consys_equiscale (consys_struct *consys,
		       double **p_rowscale, double **p_colscale)
/*
  Given the constraint matrix A (consys), this routine will calculate diagonal
  scaling matrices R (rowscale) and S (colscale) using equilibration scaling.

  The routine assumes it is scaling the matrix R*A*S, and the parameters
  rowscale and colscale are updated with new scaling coefficients. The
  constraint system is >> not << scaled.

  This routine is light on paranoia, on the assumption that the matrix has
  just been scanned with consys_evalsys, which is heavy on paranoia.

  Parameters:
    consys:	constraint system
    p_rowscale:	(i) initial row scaling coefficients; created if null
		(o) revised row scaling coefficients
    p_colscale:	(i) initial column scaling coefficients; created if null
		(o) revised column scaling coefficients

  Returns: TRUE if no errors occurred while calculating the scaling
	   coefficients, FALSE otherwise.
*/

{ int i,j ;
  double sqm,sqm_old,eps,aij,rcmax,rcmin,maxaij,minaij ;
  double *rowscale,*colscale ;
  coeff_struct *coeffij ;

# if defined(DYLP_PARANOIA) || CONSYS_SCALING_DEBUG >= 1

  char *rtnnme = "consys_equiscale" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(2,rtnnme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(2,rtnnme,"column header") ;
    return (FALSE) ; }
  if (p_rowscale == NULL)
  { errmsg(2,rtnnme,"&rowscale") ;
    return (FALSE) ; }
  if (p_colscale == NULL)
  { errmsg(2,rtnnme,"&colscale") ;
    return (FALSE) ; }
# endif

/*
  If the client didn't supply initial scaling matrices, create them now and
  initialise them to 1.
*/
  if (*p_rowscale == NULL)
  { rowscale = (double *) MALLOC((consys->concnt+1)*sizeof(double)) ;
    rowscale[0] = 0 ;
    for (i = 1 ; i <= consys->concnt ; i++) rowscale[i] = 1.0 ; }
  else
  { rowscale = *p_rowscale ; }
  if (*p_colscale == NULL)
  { colscale = (double *) MALLOC((consys->varcnt+1)*sizeof(double)) ;
    colscale[0] = 0 ;
    for (j = 1 ; j <= consys->varcnt ; j++) colscale[j] = 1.0 ; }
  else
  { colscale = *p_colscale ; }

  sqm_old = sqrt(consys->maxaij/consys->minaij) ;
  eps = 1.0 ;

/*
  Update the column scaling vector.
*/
  maxaij = 0.0 ;
  minaij = consys->inf ;
  for (j = 1 ; j <= consys->varcnt ; j++)
  { coeffij = consys->mtx.cols[j]->coeffs ;
    if (coeffij == NULL) continue ;
    rcmax = 0.0 ;
    rcmin = consys->inf ;
    for ( ; coeffij != NULL ; coeffij = coeffij->colnxt)
    { aij = fabs(coeffij->val) ;
      if (aij == 0) continue ;
      i = coeffij->rowhdr->ndx ;
      aij *= rowscale[i] ;
      if (aij > rcmax) rcmax = aij ;
      if (aij < rcmin) rcmin = aij ; }
    colscale[j] = 1/rcmax ;
    if (colscale[j]*rcmax > maxaij) maxaij = colscale[j]*rcmax ;
    if (colscale[j]*rcmin < minaij) minaij = colscale[j]*rcmin ; }
  sqm = sqrt(maxaij/minaij) ;
  eps = (sqm_old-sqm)/sqm_old ;
# if (CONSYS_SCALING_DEBUG >= 1)
  dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n  %s: %g <= a<ij> <= %g, geom = %g, eps = %g",
	      rtnnme,minaij,maxaij,sqm,eps) ;
# endif

  consys->maxaij = maxaij ;
  consys->minaij = minaij ;

  *p_rowscale = rowscale ;
  *p_colscale = colscale ;

  return (TRUE) ; }



bool consys_applyscale (consys_struct *consys, bool convctyp,
			double *rowscale, double *colscale)

/*
  This routine applies the scaling matrices rowscale and colscale to the
  constraint system in consys. In addition to the coefficient matrix, the
  scaling is applied to the objective, right-hand-side, and bounds, if they
  are attached.

  Note that constraint upper and lower bounds are NOT scaled. They should
  be recalculated.

  Note that dylp uses negative rowscale values to convert >= constraints to
  <= constraints. If this is occuring, convctyp should be TRUE for exactly
  one call.

  Arguably, this routine could be made more efficient in the case where we're
  scaling for the sole purpose of flipping >= constraints to <= constraints.
  Just scan rowscale and change constraint system values only for values that
  are -1.0. It's not clear that this case occurs often enough to be worth the
  trouble. Binary coefficient matrices, maybe.

  Parameters:
    consys:	constraint system to be scaled
    convctyp:	TRUE if the routine should scan rowscale for negative values
		and convert >= constraints to <= constraints; FALSE otherwise
    rowscale:	row scaling matrix
    colscale:	column scaling matrix (may be NULL under the special condition
		explained above)

  Returns: TRUE if the system is successfully scaled, FALSE otherwise.
*/

{ int i,j ;
  double aij,maxaij,minaij ;
  coeff_struct *coeffij ;

# ifdef DYLP_PARANOIA

  char *rtnnme = "consys_applyscale" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(2,rtnnme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(2,rtnnme,"column header") ;
    return (FALSE) ; }
  if (rowscale == NULL)
  { errmsg(2,rtnnme,"rowscale") ;
    return (FALSE) ; }
  if (colscale == NULL)
  { errmsg(2,rtnnme,"colscale") ;
    return (FALSE) ; }
# endif

/*
  Perform row scaling on the coefficient matrix and right-hand-side vectors.
*/
  for (i = 1 ; i <= consys->concnt ; i++)
  { coeffij = consys->mtx.rows[i]->coeffs ;
    for ( ; coeffij != NULL ; coeffij = coeffij->rownxt)
      coeffij->val *= rowscale[i] ; }
  if (consys->rhs != NULL)
  { for (i = 1 ; i <= consys->concnt ; i++)
      consys->rhs[i] *= rowscale[i] ; }
  if (consys->rhslow != NULL)
  { for (i = 1 ; i <= consys->concnt ; i++)
      consys->rhslow[i] *= rowscale[i] ; }
  if (convctyp == TRUE && consys->ctyp != NULL)
  { for (i = 1 ; i <= consys->concnt ; i++)
    { if (rowscale[i] < 0)
      { 
#       ifdef DYLP_PARANOIA
	if (consys->ctyp[i] != contypGE)
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (FALSE) ; }
#       endif
	consys->ctyp[i] = contypLE ; } } }
/*
  Perform column scaling on the coefficient matrix, objective coefficients,
  and bounds. The bounds are scaled by 1/S so that the implicit coefficients
  remain 1.0.

  Finite infinity is a pain here: We have to test for it in bounds vectors
  and avoid scaling, because after scaling, it won't be infinity any more.
*/
  maxaij = 0.0 ;
  minaij = consys->inf ;
  for (j = 1 ; j <= consys->varcnt ; j++)
  { coeffij = consys->mtx.cols[j]->coeffs ;
    for ( ; coeffij != NULL ; coeffij = coeffij->colnxt)
    { coeffij->val *= colscale[j] ;
      if (coeffij->val != 0.0)
      { aij = fabs(coeffij->val) ;
	if (aij < minaij) minaij = aij ;
	if (aij > maxaij) maxaij = aij ; } } }
  if (consys->obj != NULL)
  { for (j = 1 ; j <= consys->varcnt ; j++)
      consys->obj[j] *= colscale[j] ; }
  if (flgoff(consys->opts,CONSYS_FININF))
  { if (consys->vlb != NULL)
    { for (j = 1 ; j <= consys->varcnt ; j++)
	consys->vlb[j] /= colscale[j] ; }
    if (consys->vub != NULL)
    { for (j = 1 ; j <= consys->varcnt ; j++)
	consys->vub[j] /= colscale[j] ; } }
  else
  { if (consys->vlb != NULL)
    { for (j = 1 ; j <= consys->varcnt ; j++)
      { if (consys->vlb[j] > -consys->inf)
	  consys->vlb[j] /= colscale[j] ; } }
    if (consys->vub != NULL)
    { for (j = 1 ; j <= consys->varcnt ; j++)
      { if (consys->vub[j] < consys->inf)
	  consys->vub[j] /= colscale[j] ; } } }

  consys->maxaij = maxaij ;
  consys->minaij = minaij ;

  return (TRUE) ; }

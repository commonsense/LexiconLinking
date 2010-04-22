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
  This file contains utility routines for the packed vector structure, as
  well as a few utilities for general expanded vectors.
*/



#include "dylib_errs.h"
#include "dylib_std.h"
#include "dylib_strrtns.h"
#include "dy_vector.h"

static char sccsid[] UNUSED = "@(#)vector_utils.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_vector_utils.c 269 2009-04-02 05:38:19Z lou $" ;

static const char *noname = "<<n/a>>" ;



pkvec_struct *pkvec_new (int sze)

/*
  This routine allocates a new packed vector structure.

  Parameters:
    sze:	the allocated size of the coeffs array

  Returns: pointer to the packed vector, or null if allocation fails
*/

{ pkvec_struct *pkvec ;

  const char *rtnnme = "pkvec_new" ;

  if (sze < 0) sze = 0 ;
  pkvec = (pkvec_struct *) CALLOC(1,sizeof(pkvec_struct)) ;
  pkvec->sze = sze ;
  pkvec->nme = noname ;
  if (sze == 0)
  { pkvec->coeffs = NULL ; }
  else
  { pkvec->coeffs = (pkcoeff_struct *) MALLOC(sizeof(pkcoeff_struct)*sze) ;
    if (pkvec->coeffs == NULL)
    { errmsg(8,rtnnme,__LINE__,sizeof(pkcoeff_struct)*sze) ;
      return (NULL) ; } }
  
  return (pkvec) ; }



bool pkvec_resize (pkvec_struct *pkvec, int sze)

/*
  This routine resizes a packed vector as specified by sze. If sze is 0, the
  default is to expand by 10% of the current size, with a minimum expansion
  of 10. It's an error if the new size is less than the current number of
  coefficients in the vector.

  Parameters:
    sze:	the new allocated size of the vector
    pkvec:	the packed vector

  Returns: TRUE if the expansion succeeds, FALSE otherwise.
*/

{ pkcoeff_struct *coeffs ;

  const char *rtnnme = "pkvec_resize" ;

# ifdef DYLP_PARANOIA
  if (pkvec == NULL)
  { errmsg(2,rtnnme,"pkvec") ;
    return (FALSE) ; }
  if ((pkvec->coeffs == NULL && (pkvec->sze > 0 || pkvec->cnt > 0)) ||
      (pkvec->cnt > pkvec->sze))
  { errmsg(90,rtnnme,(pkvec->nme == NULL)?"<<null>>":pkvec->nme,
	   pkvec->ndx,pkvec->sze,pkvec->cnt,(pkvec->coeffs == NULL)?"un":"") ;
    return (FALSE) ; }
#endif

  if (sze == 0) sze = minn((int)(pkvec->sze*1.1),pkvec->sze+10) ;
  if (sze < pkvec->cnt)
  { errmsg(91,rtnnme,(pkvec->nme == NULL)?"<<null>>":pkvec->nme,
	   pkvec->ndx,pkvec->cnt,sze) ;
    return (FALSE) ; }

  coeffs = pkvec->coeffs ;
  pkvec->coeffs = 
      (pkcoeff_struct *) REALLOC(pkvec->coeffs,sizeof(pkcoeff_struct)*sze) ;
  if (pkvec->coeffs == NULL)
  { errmsg(8,rtnnme,__LINE__,sizeof(pkcoeff_struct)*sze) ;
    pkvec->coeffs = coeffs ;
    return (FALSE) ; }
  pkvec->sze = sze ;

  return (TRUE) ; }



void pkvec_free (pkvec_struct *pkvec)

/*
  This routine frees a packed vector.

  Parameters:
    pkvec:	packed vector

  Returns: undefined.
*/

{ 

# ifdef DYLP_PARANOIA
  const char *rtnnme = "pkvec_free" ;

  if (pkvec == NULL)
  { errmsg(2,rtnnme,"pkvec") ;
    return ; }
# endif

  if (pkvec->coeffs != NULL) FREE(pkvec->coeffs) ;
  FREE(pkvec) ;

  return ; }



bool pkvec_check (pkvec_struct *pkvec, const char *caller)

/*
  This routine goes over pkvec and decides if it's consistent, according to
  the following rules:
  * nme != NULL (this is a sort of sanity check -- failure here likely means
    that the vector has been corrupted somewhere along the way, or was
    sloppily crafted without calling pkvec_new)
  * sze >= 0
  * sze == 0 iff coeffs == NULL
  * ndx >= 0 and isnan(dflt) == FALSE

  A vector can have sze == 0 and cnt > 0 in the case where it's used to acquire
  header information without reading the coefficients of the row or column. The
  allowance for ndx == 0 is made with reluctance, but it's necessary when
  loading a new matrix -- the index is filled in once the row/column is
  installed in the matrix.

  If sze > 0, the coeffs array is assumed to be valid, and the following two
  checks are also applied:
  * 0 <= cnt <= sze
  * coeffs[*].ndx > 0 and isnan(coeffs[*].val) == FALSE

  If any problems are found, the pkvec_check issues an error message on behalf
  of its calling routine.

  Parameters:
    pkvec:	packed vector
    caller:	name of calling routine
  
  Returns: TRUE if the vector checks out, FALSE otherwise.
*/

{ int ndx ;
  const char *rtnnme = "pkvec_check" ; 

# ifdef DYLP_PARANOIA
  if (pkvec == NULL)
  { errmsg(2,rtnnme,"pkvec") ;
    return (FALSE) ; }
# endif
  if (caller == NULL) caller = rtnnme ;

  if (pkvec->nme == NULL)
  { errmsg(95,caller,pkvec) ;
    return (FALSE) ; }

  if (pkvec->sze < 0 ||
      (pkvec->sze == 0 && pkvec->coeffs != NULL) ||
      (pkvec->sze != 0 && pkvec->coeffs == NULL))
  { errmsg(90,caller,pkvec->nme,pkvec->ndx,
	   pkvec->sze,pkvec->cnt,(pkvec->coeffs == NULL)?"un":"") ;
    return (FALSE) ; }
  
  if (pkvec->ndx < 0 ||
      isnan(pkvec->dflt))
  { errmsg(93,caller,pkvec->nme,pkvec->ndx,pkvec->dflt) ;
    return (FALSE) ; }

  if (pkvec->sze == 0) return (TRUE) ;

  if (pkvec->cnt < 0 || pkvec->cnt > pkvec->sze)
  { errmsg(90,caller,pkvec->nme,pkvec->ndx,
	   pkvec->sze,pkvec->cnt,(pkvec->coeffs == NULL)?"un":"") ;
    return (FALSE) ; }
  
  for (ndx = 0 ; ndx < pkvec->cnt ; ndx++)
  { if (pkvec->coeffs[ndx].ndx < 0 ||
	isnan(pkvec->coeffs[ndx].val))
    { errmsg(94,caller,pkvec->nme,pkvec->ndx,ndx,pkvec->coeffs[ndx].ndx,
	     ndx,pkvec->coeffs[ndx].val) ;
      return (FALSE) ; } }

  return (TRUE) ; }



double exvec_1norm (double *vec, int len)

/*
  Simple utility routine to calculate the 1-norm SUM{j} |vec<j>| of an
  expanded vector. The vector is assumed to be indexed from 1 to len.

  Parameters:
    vec:	expanded vector
    len:	length of the vector

  Returns: the 1-norm of the vector, or NaN if there's a problem.
*/

{ int ndx ;
  double norm ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "exvec_1norm" ;

  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (ndx = 1 ; ndx <= len ; ndx++) norm += fabs(vec[ndx]) ;

  return (norm) ; }



double exvec_ssq (double *vec, int len)

/*
  Simple utility routine to calculate the sum of squares SUM{j} vec<j>**2 of
  an expanded vector. It's often more useful to have this than the actual
  2-norm. The vector is assumed to be indexed from 1 to len.

  Parameters:
    vec:	expanded vector
    len:	length of the vector

  Returns: the sum of squares of the vector, or NaN if there's a problem.
*/

{ int ndx ;
  double norm ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "exvec_ssq" ;

  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (ndx = 1 ; ndx <= len ; ndx++) norm += vec[ndx]*vec[ndx] ;

  return (norm) ; }


double exvec_2norm (double *vec, int len)

/*
  Simple utility routine to calculate the 2-norm sqrt(SUM{j} vec<j>**2) of
  an expanded vector. The vector is assumed to be indexed from 1 to len.

  Parameters:
    vec:	expanded vector
    len:	length of the vector

  Returns: the 2-norm of the vector, or NaN if there's a problem.
*/

{ int ndx ;
  double norm ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "exvec_2norm" ;

  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (ndx = 1 ; ndx <= len ; ndx++) norm += vec[ndx]*vec[ndx] ;

  return (sqrt(norm)) ; }


double pkvec_2norm (pkvec_struct *vec)

/*
  Simple utility routine to calculate the 2-norm sqrt(SUM{j} vec<j>**2) of
  a packed vector.

  Parameters:
    vec:	packed vector

  Returns: the 2-norm of the vector, or NaN if there's a problem.
*/

{ int ndx ;
  pkcoeff_struct *coeffs ;
  double norm ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "pkvec_2norm" ;

  if (vec == NULL)
  { errmsg(2,rtnnme,"pkvec") ;
    return (quiet_nan(0)) ; }
  if (vec->coeffs == NULL)
  { errmsg(2,rtnnme,"pkvec coeffs") ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  coeffs = vec->coeffs ;
  for (ndx = 0 ; ndx < vec->cnt ; ndx++)
    norm += coeffs[ndx].val*coeffs[ndx].val ;

  return (sqrt(norm)) ; }



double exvec_infnorm (double *vec, int len, int *p_jmax)

/*
  Simple utility routine to calculate the infinity-norm MAX{j} |vec<j>| of
  an expanded vector. The vector is assumed to be indexed from 1 to len.

  Parameters:
    vec:	expanded vector
    len:	length of the vector
    p_jmax:	(o) if non-null, will be set to index of max value

  Returns: the inf-norm (max) of the vector, or NaN if there's a problem.
*/

{ int j,jmax ;
  double norm ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "exvec_infnorm" ;

  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    if (p_jmax != NULL) *p_jmax = -1 ;
    return (quiet_nan(0)) ; }
# endif

/*
  Initialising (norm, jmax) to (0.0, len) will give the proper result for
  a vector of length 0, and also suppresses a compiler `possible rui' warning
  for jmax.
*/
  norm = 0.0 ;
  if (p_jmax != NULL)
  { jmax = len ;
    for (j = 1 ; j <= len ; j++)
    { if (fabs(vec[j]) > norm)
      { norm = fabs(vec[j]) ;
	jmax = j ; } }
    *p_jmax = jmax ; }
  else
  { for (j = 1 ; j <= len ; j++) norm = maxx(fabs(vec[j]),norm) ; }

  return (norm) ; }



double pkvec_dotexvec (pkvec_struct *pkvec, double *exvec)

/*
  Utility routine to calculate the dot product SUM{j} pkvec<j>*exvec<j>
  of a packed vector and an expanded vector.

  Parameters:
    pkvec:	packed vector
    exvec:	expanded vector

  Returns: the dot product of the two vectors, or NaN if there's a problem
*/

{ int pkndx ;
  double dot ;
  pkcoeff_struct *coeffs ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "pkvec_dotexvec" ;

  if (pkvec == NULL)
  { errmsg(2,rtnnme,"pkvec") ;
    return (quiet_nan(0)) ; }
  if (pkvec->coeffs == NULL)
  { errmsg(2,rtnnme,"pkvec coeffs") ;
    return (quiet_nan(0)) ; }
  if (exvec == NULL)
  { errmsg(2,rtnnme,"exvec") ;
    return (quiet_nan(0)) ; }
# endif

  dot = 0 ;
  coeffs = pkvec->coeffs ;
  for (pkndx = 0 ; pkndx < pkvec->cnt ; pkndx++)
    dot += coeffs[pkndx].val*exvec[coeffs[pkndx].ndx] ;

  return (dot) ; }

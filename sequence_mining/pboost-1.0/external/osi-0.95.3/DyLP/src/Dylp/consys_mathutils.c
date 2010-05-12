/*
  This file is a portion of the OsiDylp LP distribution.

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
  This file contains math utility routines for the constraint system data
  structure. These handle things like dot products, norms, etc.
*/



#include "dylib_errs.h"
#include "dylib_std.h"
#include "consys.h"

static char sccsid[] UNUSED = "@(#)consys_mathutils.c	4.5	11/11/04" ;
static char svnid[] UNUSED = "$Id: consys_mathutils.c 71 2006-06-09 04:21:15Z andreasw $" ;



double consys_dotrow (consys_struct *consys, int rowndx, double *vec)

/*
  This routine computes the dot product of the specified row with the expanded
  vector passed in vec.

  Parameters:
    consys:	constraint system
    rowndx:	row
    vec:	vector

  Returns: dot product, or NaN if the calculation goes awry.
*/

{ double dotprod ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_dotrow" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (quiet_nan(0)) ; }
# endif
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (quiet_nan(0)) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (quiet_nan(0)) ; }
  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (quiet_nan(0)) ; }
# endif

  dotprod = 0 ;
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    dotprod += coeff->val*vec[coeff->colhdr->ndx] ; }

  return (dotprod) ; }



double consys_dotcol (consys_struct *consys, int colndx, double *vec)

/*
  This routine computes the dot product of the specified column with the
  expanded vector passed in vec.

  Parameters:
    consys:	constraint system
    colndx:	column
    vec:	vector

  Returns: dot product, or NaN if the calculation goes awry.
*/

{ double dotprod ;
  colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_dotcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (quiet_nan(0)) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (quiet_nan(0)) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (quiet_nan(0)) ; }
  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (quiet_nan(0)) ; }
# endif

  dotprod = 0 ;
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->concnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->concnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    dotprod += coeff->val*vec[coeff->rowhdr->ndx] ; }

  return (dotprod) ; }



double consys_1normrow (consys_struct *consys, int rowndx)

/*
  This routine computes the 1-norm of a row: SUM{j} |a<i,j>|

  Parameters:
    consys:	constraint system
    rowndx:	row

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_1normrow" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (quiet_nan(0)) ; }
# endif
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (quiet_nan(0)) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#endif
    norm += fabs(coeff->val) ; }

  return (norm) ; }


double consys_ssqrow (consys_struct *consys, int rowndx)

/*
  This routine computes the sum of squares of a row: SUM{j} a<i,j>**2. It's
  sometimes more useful to have this than the actual 2-norm.

  Parameters:
    consys:	constraint system
    rowndx:	row

  Returns: value of the sum of squares, or NaN if the calculation goes awry
*/

{ double norm ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_ssqrow" ;
# endif

/*
  The usual paranoia, plus an honest index check.
*/
# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (quiet_nan(0)) ; }
# endif
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (quiet_nan(0)) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm += coeff->val*coeff->val ; }

  return (norm) ; }


double consys_2normrow (consys_struct *consys, int rowndx)

/*
  This routine computes the 2-norm of a row: sqrt(SUM{j} a<i,j>**2)

  Parameters:
    consys:	constraint system
    rowndx:	row

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_2normrow" ;
# endif

/*
  The usual paranoia, plus an honest index check.
*/
# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (quiet_nan(0)) ; }
# endif
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (quiet_nan(0)) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm += coeff->val*coeff->val ; }

  return (sqrt(norm)) ; }


double consys_infnormrow (consys_struct *consys, int rowndx)

/*
  This routine computes the infinity-norm of a row: MAX{j} |a<i,j>|

  Parameters:
    consys:	constraint system
    rowndx:	row

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_infnormrow" ;
# endif

/*
  The usual paranoia, plus an honest index check.
*/
# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (quiet_nan(0)) ; }
# endif
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (quiet_nan(0)) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm = maxx(fabs(coeff->val),norm) ; }

  return (norm) ; }



double consys_1normcol (consys_struct *consys, int colndx)

/*
  This routine computes the 1-norm of a column: SUM{i} |a<i,j>|.

  Parameters:
    consys:	constraint system
    colndx:	column

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_1normcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (quiet_nan(0)) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (quiet_nan(0)) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm += fabs(coeff->val) ; }

  return (norm) ; }


double consys_ssqcol (consys_struct *consys, int colndx)

/*
  This routine computes the sum of squares of a column: SUM{i} a<i,j>**2.
  It's sometimes more useful to have this than the actual 2-norm.

  Parameters:
    consys:	constraint system
    colndx:	column

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_ssqcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (quiet_nan(0)) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (quiet_nan(0)) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm += coeff->val*coeff->val ; }

  return (norm) ; }


double consys_2normcol (consys_struct *consys, int colndx)

/*
  This routine computes the 2-norm of a column: sqrt(SUM{i} a<i,j>**2).

  Parameters:
    consys:	constraint system
    colndx:	column

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_2normcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (quiet_nan(0)) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (quiet_nan(0)) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm += coeff->val*coeff->val ; }

  return (sqrt(norm)) ; }


double consys_infnormcol (consys_struct *consys, int colndx)

/*
  This routine computes the infinity-norm of a column: MAX{i} |a<i,j>|.

  Parameters:
    consys:	constraint system
    colndx:	column

  Returns: value of the norm, or NaN if the calculation goes awry
*/

{ double norm ;
  colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_infnormcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (quiet_nan(0)) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (quiet_nan(0)) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (quiet_nan(0)) ; }
# endif

  norm = 0 ;
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->varcnt) ;
      return (quiet_nan(0)) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    norm = maxx(fabs(coeff->val),norm) ; }

  return (norm) ; }



bool consys_mulrow (consys_struct *consys, int rowndx, double scalar)

/*
  This routine multiplies a row i by a scalar q. It deals with the coefficients
  a<i>, and also with b<i>, blow<i>, cub<i>, and clb<i>, if they exist. If
  q < 0, the type of constraint is changed accordingly (>= swapped with <=)
  and clb<i> is swapped with cub<i>.

  Note that range constraints always take the form blow <= ax <= b, so if we
  multiply a range constraint by q < 0, the resulting constraint is
  qblow >= (qa)x >= qb => qb <= (qa)x <= qblow.

  Attempting to multiply a constraint by 0 gets you a warning if the
  CONSYS_WRNZERO flag is set in consys->opts.

  The routine will work with clb<i> and cub<i> only if both are present. It's
  difficult to define consistent changes otherwise.

  Parameters:
    consys:	constraint system
    rowndx:	row to be modified
    scalar:	the multiplicative scalar

  Returns: TRUE if no problems are encountered, FALSE otherwise.
*/

{ double tmprhs ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  conbnd_struct tmpbnd ;
  bool do_conbnds ;

  const char *rtnnme = "consys_mulrow" ;

/*
  The usual paranoia, plus an honest index check.
*/
# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (FALSE) ; }
# endif

  rowhdr = consys->mtx.rows[rowndx] ;

# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (FALSE) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (scalar == 0 && flgon(consys->opts,CONSYS_WRNZERO))
  { warn(132,rtnnme,consys->nme,"row",rowhdr->nme,rowndx) ; }
# endif
  if (consys->cub != NULL && consys->clb != NULL)
    do_conbnds = TRUE ;
  else
    do_conbnds = FALSE ;
/*
  The straightforward part. Multiply the coefficients by the scalar.
*/
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (FALSE) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (FALSE) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (FALSE) ; }
#   endif
    coeff->val *= scalar ; }
/*
  If we did get a 0 for the scalar, we can be done in no time.
*/
  if (scalar == 0)
  { if (consys->rhs != NULL) consys->rhs[rowndx] = 0 ;
    if (consys->rhslow != NULL) consys->rhslow[rowndx] = 0 ;
    if (do_conbnds == TRUE)
    { tmpbnd.revs = 0 ;
      tmpbnd.inf = 0 ;
      tmpbnd.bnd = 0 ;
      consys->cub[rowndx] = tmpbnd ;
      consys->clb[rowndx] = tmpbnd ; }
    return (TRUE) ; }
/*
  For q != 0, it's a little more work. Correct b<i>, blow<i>, cub<i>, and
  clb<i>, if they exist.
*/
  if (consys->rhs != NULL) consys->rhs[rowndx] *= scalar ;
  if (consys->rhslow != NULL) consys->rhslow[rowndx] *= scalar ;
  if (do_conbnds == TRUE)
  { consys->cub[rowndx].bnd *= scalar ;
    consys->clb[rowndx].bnd *= scalar ; }
/*
  And now the complicated bit. If q < 0, swap the constraint bounds, then take
  additional action as needed, depending on the constraint type.
*/
  if (scalar < 0)
  { if (do_conbnds == TRUE)
    { tmpbnd = consys->cub[rowndx] ;
      consys->cub[rowndx] = consys->clb[rowndx] ;
      consys->clb[rowndx] = tmpbnd ; }
    switch (consys->ctyp[rowndx])
    { case contypLE:
      { consys->ctyp[rowndx] = contypGE ;
	break ; }
      case contypGE:
      { consys->ctyp[rowndx] = contypLE ;
	break ; }
      case contypRNG:
      { tmprhs = consys->rhs[rowndx] ;
	consys->rhs[rowndx] = consys->rhslow[rowndx] ;
	consys->rhslow[rowndx] = tmprhs ;
	break ; }
      case contypEQ:
      case contypNB:
      { break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (FALSE) ; } } }

  return (TRUE) ; }



bool consys_divrow (consys_struct *consys, int rowndx, double scalar)

/*
  This routine divides a row i by a scalar q. It deals with the coefficients
  a<i>, and also with b<i>, blow<i>, cub<i>, and clb<i>, if they exist. If
  q < 0, the type of constraint is changed accordingly (>= swapped with <=)
  and clb<i> is swapped with cub<i>. It's a separate routine (rather than
  using consys_mulrow to multiply by 1/scalar) to try and retain accuracy.

  Note that range constraints always take the form blow <= ax <= b, so if we
  divide a range constraint by q < 0, the resulting constraint is
  qblow >= (qa)x >= qb => qb <= (qa)x <= qblow.

  Attempting to divide a constraint by 0 is an error.

  The routine will work with clb<i> and cub<i> only if both are present. It's
  difficult to define consistent changes otherwise.

  Parameters:
    consys:	constraint system
    rowndx:	row to be divided
    scalar:	the dividing scalar

  Returns: TRUE if no problems are encountered, FALSE otherwise.
*/

{ double tmprhs ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  conbnd_struct tmpbnd ;
  bool do_conbnds ;

  const char *rtnnme = "consys_divrow" ;

/*
  The usual paranoia, plus an honest index check.
*/
# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (FALSE) ; }
# endif

  rowhdr = consys->mtx.rows[rowndx] ;

# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (FALSE) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (FALSE) ; }
  if (scalar == 0)
  { errmsg(5,rtnnme,"scalar",(int) scalar) ;
    return (FALSE) ; }
# endif
  if (consys->cub != NULL && consys->clb != NULL)
    do_conbnds = TRUE ;
  else
    do_conbnds = FALSE ;
/*
  The straightforward part. Divide the coefficients by the scalar.
*/
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (FALSE) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (FALSE) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (FALSE) ; }
#   endif
    coeff->val /= scalar ; }
/*
  Correct b<i>, blow<i>, cub<i>, and clb<i>, if they exist.
*/
  if (consys->rhs != NULL) consys->rhs[rowndx] /= scalar ;
  if (consys->rhslow != NULL) consys->rhslow[rowndx] /= scalar ;
  if (do_conbnds == TRUE)
  { consys->cub[rowndx].bnd /= scalar ;
    consys->clb[rowndx].bnd /= scalar ; }
/*
  And now the complicated bit. If q < 0, swap the constraint bounds, then take
  additional action as needed, depending on the constraint type.
*/
  if (scalar < 0)
  { if (do_conbnds == TRUE)
    { tmpbnd = consys->cub[rowndx] ;
      consys->cub[rowndx] = consys->clb[rowndx] ;
      consys->clb[rowndx] = tmpbnd ; }
    switch (consys->ctyp[rowndx])
    { case contypLE:
      { consys->ctyp[rowndx] = contypGE ;
	break ; }
      case contypGE:
      { consys->ctyp[rowndx] = contypLE ;
	break ; }
      case contypRNG:
      { tmprhs = consys->rhs[rowndx] ;
	consys->rhs[rowndx] = consys->rhslow[rowndx] ;
	consys->rhslow[rowndx] = tmprhs ;
	break ; }
      case contypEQ:
      case contypNB:
      { break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (FALSE) ; } } }

  return (TRUE) ; }



int consys_gcdrow (consys_struct *consys, int rowndx)

/*
  This routine calculates the gcd of the coefficients of the specified row
  using the euclidean algorithm. Note that explicit zeros should not appear
  in the coefficient matrix.
  
  Obviously, the coefficients should be integer. If they're not, the routine
  returns 0.

  If the row is empty, the routine returns 0, on the theory that whatever you
  were trying to do with this row, it's probably not suitable.

  The code uses the following statement of the euclidean algorithm, courtesy
  of Martin, Large Scale Linear and Integer Optimization, p. 106. Assume
  a<1> and a<2> positive integer, a<1> > a<2>.

    while (a<1> > 0 && a<2> > 0)
    { q = floor(a<1>/a<2>) ;
      r = a<1> - q*a<2> ;
      a<1> = a<2> ;
      a<2> = r ; }

  Parameters:
    consys:	constraint system
    rowndx:	row to be evaluated
  
  Returns: gcd(a<1>, ..., a<n>), 0 if the coefficients aren't integer, -1 if
	   anything else goes wrong.
*/

{ double gcd,a1,a2,q,r ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_gcdrow" ;
# endif

/*
  The usual paranoia, plus an honest index check.
*/
# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (-1) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (-1) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (-1) ; }
# endif
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (-1) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (-1) ; }
# endif
/*
  Trivial cases: 0 or 1 coefficients.
*/
  if (rowhdr->len == 0) return (1) ;
  coeff = rowhdr->coeffs ;
# ifdef PARANOIA
  if (coeff == NULL)
  { errmsg(116,rtnnme,consys->nme,rowhdr->nme,rowhdr->ndx,rowhdr->len,0) ;
    return (-1) ; }
  if (coeff->colhdr == NULL)
  { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	   consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
    return (-1) ; }
  if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	   1,consys->varcnt) ;
    return (-1) ; }
  if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
  { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	   coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
    return (-1) ; }
# endif
  a1 = coeff->val ;
  if (a1 < 0) a1 = -a1 ;
  if (floor(a1) != a1) return (0) ;
  if (rowhdr->len == 1) return ((int) a1) ;
/*
  Two or more coefficients. We work through them, calculating gcd(gcd,a<i>).
  We first do a quick test for a<i>/gcd integer (in which case we can keep
  gcd and move on to the next coefficient). When the gcd drops to 1, we bail
  out.
*/
  gcd = a1 ;
  for (coeff = coeff->rownxt ; gcd > 1 && coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (-1) ; }
    if (coeff->colhdr->ndx <= 0 || coeff->colhdr->ndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",coeff->colhdr->ndx,
	     1,consys->varcnt) ;
      return (-1) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (-1) ; }
#   endif
    a1 = coeff->val ;
    if (a1 < 0) a1 = -a1 ;
    if (floor(a1) != a1) return (0) ;
    if (a1 > gcd)
    { if (floor(a1/gcd) == a1/gcd) continue ;
      a2 = gcd ; }
    else
    { a2 = a1 ;
      a1 = gcd ; }
/*
  We need to do a gcd calculation.
*/
    while (a1 > 0 && a2 > 0)
    { q = floor(a1/a2) ;
      r = a1 - q*a2 ;
      a1 = a2 ;
      a2 = r ; }
    gcd = a1 ; }

  return ((int) gcd) ; }



bool consys_accumcol (consys_struct *consys, int colndx, double *vec)

/*
  This routine adds the column specified by colndx to the expanded vector
  passed in vec.

  Parameters:
    consys:	constraint system
    colndx:	column
    vec:	vector

  Returns: TRUE if there are no problems, FALSE otherwise.
*/

{ colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_accumcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (FALSE) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (FALSE) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (FALSE) ; }
  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (FALSE) ; }
# endif

  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (FALSE) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->concnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->concnt) ;
      return (FALSE) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (FALSE) ; }
#   endif
    vec[coeff->rowhdr->ndx] += coeff->val ; }

  return (TRUE) ; }



bool consys_mulaccumcol (consys_struct *consys, int colndx,
			 double scalar, double *vec)

/*
  This routine multiplies the column specified by colndx by scalar and then
  adds it to the expanded vector passed in vec. Identical to consys_accumcol,
  except for the multiplication.

  Parameters:
    consys:	constraint system
    colndx:	column
    scalar:	scalar multiplier for column
    vec:	vector

  Returns: TRUE if there are no problems, FALSE otherwise.
*/

{ colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_accumcol" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (FALSE) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (FALSE) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (FALSE) ; }
  if (vec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (FALSE) ; }
# endif

  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (FALSE) ; }
    if (coeff->rowhdr->ndx <= 0 || coeff->rowhdr->ndx > consys->concnt)
    { errmsg(102,rtnnme,consys->nme,"row",coeff->rowhdr->ndx,
	     1,consys->concnt) ;
      return (FALSE) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (FALSE) ; }
#   endif
    vec[coeff->rowhdr->ndx] += scalar*coeff->val ; }

  return (TRUE) ; }


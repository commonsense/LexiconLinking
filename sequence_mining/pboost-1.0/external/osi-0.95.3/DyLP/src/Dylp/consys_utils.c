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
  This file contains utility routines for maintaining a constraint system 
  data structure, tailored for LP-based branch-and-cut MILP algorithms.

  'Architectural constraints' are the constraints included in the original
  problem statement; by analogy, 'architectural variables' are the variables
  in those constraints. 'Logical variables' are the various slack, surplus,
  and artificial variables that are added to the constraint system for ease
  of solving the LP. 'Cuts', or 'cutting planes', are additional constraints
  used in the course of (M)ILP. Column generation algorithms may add new
  architectural variables. Cutting planes may also add new architectural
  variables, and will come with their own associated logical variables.

  Consider now a system with n architectural variables, m constraints
  (architectural and/or cut), and m logical variables associated with the
  constraints.

  The design goals for supporting a LP algorithm are:
    * Efficient access by row (constraint) or by column (variable),
      in either sparse or expanded vector format. 
    * Efficient 'compact' access -- constraints should be indexed 1..m with
      no gaps, and architectural variables should be indexed 1..n with no
      gaps.
    * A standard convention for deriving the index of the associated logical
      variable from the index of a constraint.
    * Efficient substitution of alternate objective, rhs, and bound vectors.
    * A compact handle for passing around a constraint system.

  The above are sufficient for simple LP, which sees a constant constraint
  system and solves it. Dynamic LP is more challenging, requiring the ability
  to work with a dynamically changing subsystem of the original constraint
  system. The implementation must provide support for addition and deletion
  of constraints and variables in such a way that maintenance of a basis and
  tracking the relationship between the original and partial systems can be
  done in a reasonably efficient manner.

  The branch-and-cut level imposes some additional goals:
    * Efficient addition and deletion of cutting planes (both the constraint
      and its associated logical variable).
    * Efficient access to a specified constraint or variable.

  The consys package maintains two classes of variables, architectural and
  logical. For reasons of efficiency (which will become clear below) we
  reverse the usual convention, numbering logicals from 1..m and
  architecturals from m+1..m+n. This permits efficient addition/deletion of
  constraints and architectural variables by reordering.

  The consys package also maintains two classes of constraints, architectural
  and cut.  The constraint indices are a single compact set, with the indices
  of architectural constraints preceding the indices of cut constraints.
  Unlike the architectural/logical variable distinction, the architectural/
  cut constraint distinction can be ignored by a client, if it so chooses.

  The addition/deletion of an architectural variable is handled as follows:
   * To add a variable at index k, the current occupant is moved to a newly
     created slot at index n+1 and the new variable is inserted at index k.
   * To delete a variable at index k, the variable is removed and the variable
     at index n is moved into the empty slot.
  In short, addition/deletion is handled by swapping between the add/delete
  location and the last architectural location. This gives a constant cost,
  independent of the size of the system.

  Addition and deletion of constraints can be handled in the same way, by
  swapping between the add/delete location and the last location. The logical
  variables are reordered to match and the appropriate move is made in the
  set of architectural variables. The extension to maintaining architectural
  and cut constraints as separate classes is straightforward.

  It's important to note that this would fail if consys actually stored
  architectural variables at indices 1..n and logicals at indices n+1..n+m.
  Adding an architectural variable would require freeing index n+1, which
  would imply moving the logical in slot n+1 to slot n+m+1.  But because
  logical indices are derived from constraint indices, we can't reorder
  logicals without reordering the constraints. Now consider the effect on the
  constraint indices of moving the logical at position n+1 to position n+m+1
  so that an architectural variable can take index n+1. Or moving from n+m to
  n-1 on deletion of an architectural. The only way to get it right is to
  physically move the whole set of logicals to open/close a hole at the end
  of the set of architecturals. This is just too expensive.

  In practice, an external client cannot specify the column for adding a
  variable or constraint. Additions always take place at the last position
  of a set; the client can only specify if the addition is an architectural
  variable, architectural constraint, or cut constraint.

  As a final comment, specifying access to a particular constraint or
  variable from the branch-and-cut level is trickier than it looks.
  Scenario: Cut (a) is added to the constraint array; suppose this gives a
  total of m constraints.  Then cut (b) is deleted from pos'n k < m, and cut
  (a) is moved from pos'n m to pos'n k to maintain compactness. When it comes
  time to delete cut (a), we need to be able to specify it uniquely and find
  it, now that it's in pos'n k.) There's no solution in place yet; guess I'll
  cross that bridge when I come to it.


  The constraint system structure can be grouped into four components:

    * A header structure, which is the 'handle' used to pass the data
      structure;
    * a constraint coefficient matrix;
    * a set of associated vectors
    * a list of attached vectors

  The set of associated vectors is comprised of:
    * constraint and variable type vectors;
    * an objective function vector;
    * a rhs vector, and (when required) a 'range' rhs vector,
      rhslow, for range constraints;
    * variable upper and lower bounds;
    * constraint upper and lower bounds (for the lhs of constraints; this
      is used by the arc consistency algorithms).
  
  The constraint matrix is the only integral part. The header contains
  pointers to the other components, which can be swapped with relative
  ease. The only real difference between associated and attached vectors is
  that there are dedicated pointers in the constraint system header for the
  associated vectors.
  
  The list of attached vectors is intended to address the problem of dynamic
  resizing. On one hand, it's desireable to have no a priori limit on, say,
  the number of cuts that can be added. But if the client has a number of
  alternate rhs vectors being used with this constraint system, they must
  also be resized when the constraint system is resized. The list of attached
  vectors addresses this -- vectors on the list are resized with the
  constraint system. Since the vector might be moved by realloc, each vector
  has a list of pointers which should be rewritten if the vector is moved.

  Extreme coefficients are always a problem, at both ends. The infinitesimal
  is easier, in some sense. By fiat, |a<ij>| < 1.0e-20 is dropped. The infinite
  is harder. There's far too much code out there that thinks that infinity
  should be a large finite number (typically, the maximum double precision
  floating point value, DBL_MAX) rather than the IEEE FP standard value for
  infinity. Finite and infinite infinity simply do not play well together.
  So that we can implement scaling and avoid scaling finite infinity, you need
  to tell consys_create what you plan to use for infinity. At some point, it
  may become obvious that the same should apply at the infinitesimal end.

  Don't get me started on NaN. Or numeric_limits in C++.
*/



#include "dylib_errs.h"
#include "dylib_io.h"
#include "dylib_std.h"
#include "dylib_strrtns.h"
#include "consys.h"

static char sccsid[] UNUSED = "@(#)consys_utils.c	4.7	10/15/05" ;
static char svnid[] UNUSED = "$Id: consys_utils.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  Default constraint system constants
*/

#define ARCHC_DFLT 2500
#define CUT_DFLT    500
#define ARCHV_DFLT 1000

#define EXPAND_MIN_DFLT 10
#define EXPAND_PCT_DFLT .1



static bool empty_col (consys_struct *consys, int colndx, bool *rescan)

/*
  This routine removes the coefficients of a column. If one of the deleted
  coefficients belongs to the current maximum row, rescan will be set to
  TRUE.

  If DYLP_NDEBUG is not defined, a check for zero length columns can be activated
  with the CONSYS_WRNZERO flag. Since this is an internal routine, index
  range checks are active only when we're paranoid.

  Parameters:
    consys:	The constraint system.
    colndx:	Column to be emptied.
    rescan:	(o) Set to TRUE if the column deletion affected the current
		    maximum length row, FALSE otherwise.
  
  Returns: TRUE if the column is emptied without error, FALSE otherwise
	   (paranoia only)
*/

{ int rowndx ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *ccoeff,*rcoeff ;


# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "empty_column" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_COLHDR)) ;
    return (FALSE) ; }
  if (rescan == NULL)
  { errmsg(2,rtnnme,"rescan") ;
    return (FALSE) ; }
  if (colndx <= 0 || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
    return (FALSE) ; }
# endif
/*
  Empty the column by repeatedly delinking and freeing the coefficient next to
  the header. Flag for a rescan if we're unlucky enough to delete an element
  from the maximum length row.
*/
  *rescan = FALSE ;
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (FALSE) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
    return (FALSE) ; }
# endif
  for (ccoeff = colhdr->coeffs ; ccoeff != NULL ; ccoeff = colhdr->coeffs)
  { colhdr->coeffs = ccoeff->colnxt ;
    rowhdr = ccoeff->rowhdr ;
#   ifdef PARANOIA
    if (rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",ccoeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (FALSE) ; }
#   endif
    rowndx = rowhdr->ndx ;
#   ifdef PARANOIA
    if (rowndx <= 0 || rowndx > consys->concnt)
    { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
      return (FALSE) ; }
    if (consys->mtx.rows[rowndx] != rowhdr)
    { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowndx,rowndx,
	     consys->mtx.rows[rowndx]) ;
      return (FALSE) ; }
#   endif
    if (rowndx == consys->maxrowndx) *rescan = TRUE ;
    if (rowhdr->coeffs == ccoeff)
      rowhdr->coeffs = ccoeff->rownxt ;
    else
    { for (rcoeff = rowhdr->coeffs ; rcoeff != NULL ; rcoeff = rcoeff->rownxt)
	if (rcoeff->rownxt == ccoeff) break ;
#     ifdef PARANOIA
      if (rcoeff == NULL)
      { errmsg(119,rtnnme,consys->nme,rowndx,colndx,ccoeff->val,
	       "column",colndx,"row",rowndx) ;
	return (FALSE) ; }
#     endif
      rcoeff->rownxt = ccoeff->rownxt ; }
    rowhdr->len-- ;
#   ifndef DYLP_NDEBUG
    if (rowhdr->len == 0 && flgon(consys->opts,CONSYS_WRNZERO))
      warn(118,rtnnme,consys->nme,"row",rowhdr->nme,rowndx) ;
#   endif
    FREE(ccoeff) ; }
  consys->mtx.coeffcnt -= colhdr->len ;
  colhdr->len = 0 ;

  return (TRUE) ; }



static bool empty_row (consys_struct *consys, int rowndx, bool *rescan)

/*
  This routine removes the coefficients of a row. The structure and flow
  parallels empty_col, but there's really no way to reduce the amount of code
  by merging the two routines.

  If DYLP_NDEBUG is not defined, a check for zero length columns can be activated
  with the CONSYS_WRNZERO flag. Since this is an internal routine, index
  range checks are active only when we're paranoid.

  Parameters:
    consys:	The constraint system.
    rowndx:	Row to be emptied.
    rescan:	(o) Set to TRUE if the row deletion affected the current
		    maximum length column, FALSE otherwise.
  
  Returns: TRUE if the row is emptied without error, FALSE otherwise
	   (paranoia only).
*/

{ int colndx ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *ccoeff,*rcoeff ;


# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "empty_row" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_COLHDR)) ;
    return (FALSE) ; }
  if (rescan == NULL)
  { errmsg(2,rtnnme,"rescan") ;
    return (FALSE) ; }
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,consys->concnt) ;
    return (FALSE) ; }
# endif
/*
  Empty the row by repeatedly delinking and freeing the coefficient next to
  the header. Flag for a rescan if we're unlucky enough to delete an element
  from the maximum length column.
*/
  *rescan = FALSE ;
  rowhdr = consys->mtx.rows[rowndx] ;
# ifdef PARANOIA
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (FALSE) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
    return (FALSE) ; }
# endif
  for (rcoeff = rowhdr->coeffs ; rcoeff != NULL ; rcoeff = rowhdr->coeffs)
  { rowhdr->coeffs = rcoeff->rownxt ;
    colhdr = rcoeff->colhdr ;
#   ifdef PARANOIA
    if (colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",rcoeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (FALSE) ; }
#   endif
    colndx = colhdr->ndx ;
#   ifdef PARANOIA
    if (colndx <= 0 || colndx > consys->varcnt)
    { errmsg(102,rtnnme,consys->nme,"column",colndx,1,consys->varcnt) ;
      return (FALSE) ; }
    if (consys->mtx.cols[colndx] != colhdr)
    { errmsg(126,rtnnme,consys->nme,"column",colhdr,colndx,colndx,
	     consys->mtx.cols[colndx]) ;
      return (FALSE) ; }
#   endif
    if (colndx == consys->maxcolndx) *rescan = TRUE ;
    if (colhdr->coeffs == rcoeff)
      colhdr->coeffs = rcoeff->colnxt ;
    else
    { for (ccoeff = colhdr->coeffs ; ccoeff != NULL ; ccoeff = ccoeff->colnxt)
	if (ccoeff->colnxt == rcoeff) break ;
#     ifdef PARANOIA
      if (ccoeff == NULL)
      { errmsg(119,rtnnme,consys->nme,rowndx,colndx,rcoeff->val,
	       "row",rowndx,"column",colndx) ;
	return (FALSE) ; }
#     endif
      ccoeff->colnxt = rcoeff->colnxt ; }
    colhdr->len-- ;
#   ifndef DYLP_NDEBUG
    if (colhdr->len == 0 && flgon(consys->opts,CONSYS_WRNZERO))
      warn(118,rtnnme,consys->nme,"column",colhdr->nme,colndx) ;
#   endif
    FREE(rcoeff) ; }
  consys->mtx.coeffcnt -= rowhdr->len ;
  rowhdr->len = 0 ;

  return (TRUE) ; }



static bool move_col (consys_struct *consys, int fndx, int tndx)

/*
  This routine relocates a column from position fndx to position tndx. It
  is assumed that there is space at tndx. The routine adjusts the column
  header and any relevant attached arrays.

  Since this is an internal routine, index range checks are active only if
  we're paranoid.

  Parameters:
    consys:	The constraint matrix
    fndx:	The 'from' index
    tndx:	The 'to' index

  Returns: TRUE if the transfer completes without error, FALSE otherwise.
	   The routine will not fail unless paranoid or debugging checks
	   are enabled.
*/

{ attvhdr_struct *attvhdr ;
  colhdr_struct *colhdr ;

# ifdef PARANOIA

  const char *rtnnme = "move_col" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_COLHDR)) ;
    return (FALSE) ; }
  if (fndx <= 0 || fndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",fndx,1,consys->varcnt) ;
    return (FALSE) ; }
  if (tndx <= 0 || tndx > consys->varcnt+1)
  { errmsg(102,rtnnme,consys->nme,"column",tndx,1,consys->varcnt+1) ;
    return (FALSE) ; }
# endif
/*
  Get the the column header and relocate the column.
*/
  colhdr = consys->mtx.cols[fndx] ;
  consys->mtx.cols[tndx] = colhdr ;
  colhdr->ndx = tndx ;
  if (consys->maxcolndx == fndx) consys->maxcolndx = tndx ;
  if (consys->xzndx == fndx) consys->xzndx = tndx ;
/*
  Walk the attached vectors list and make the corresponding change in any of
  the row vectors.
*/
  for (attvhdr = consys->attvecs ; attvhdr != NULL ; attvhdr = attvhdr->nxt)
  { 
#   ifdef PARANOIA
    if (!(VALID_ATTVTYPE(attvhdr->what)))
    { errmsg(114,rtnnme,consys->nme,"attached",attvhdr->what) ;
      return (FALSE) ; }
    if (attvhdr->vec == NULL)
    { errmsg(127,rtnnme,consys->nme,attvhdr,
	     consys_assocnme(NULL,attvhdr->what)) ;
      return (FALSE) ; }
#   endif

    if (flgon(attvhdr->what,CONSYS_ROWVEC))
      memcpy(((char *) attvhdr->vec)+tndx*attvhdr->elsze,
	     ((char *) attvhdr->vec)+fndx*attvhdr->elsze,attvhdr->elsze) ; }

  return (TRUE) ; }



static bool move_row (consys_struct *consys, int fndx, int tndx)

/*
  This routine relocates a row from position fndx to position tndx. It
  is assumed that there is space at tndx. The routine adjusts the row
  header and any relevant attached arrays.

  Since this is an internal routine, index range checks are active only if
  we're paranoid.

  Parameters:
    consys:	The constraint matrix
    fndx:	The 'from' index
    tndx:	The 'to' index

  Returns: TRUE if the transfer completes without error, FALSE otherwise.
	   The routine will not fail unless paranoid or debugging checks
	   are enabled.
*/

{ attvhdr_struct *attvhdr ;
  rowhdr_struct *rowhdr ;

# ifdef PARANOIA

  const char *rtnnme = "move_row" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
  if (fndx <= 0 || fndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",fndx,1,consys->concnt) ;
    return (FALSE) ; }
  if (tndx <= 0 || tndx > consys->concnt+1)
  { errmsg(102,rtnnme,consys->nme,"row",tndx,1,consys->concnt+1) ;
    return (FALSE) ; }
# endif
/*
  Get the the row header and relocate the row.
*/
  rowhdr = consys->mtx.rows[fndx] ;
  consys->mtx.rows[tndx] = rowhdr ;
  rowhdr->ndx = tndx ;
  if (consys->maxrowndx == fndx) consys->maxrowndx = tndx ;
  if (consys->objndx == fndx) consys->objndx = tndx ;
/*
  Walk the attached vectors and make the corresponding change in any of the
  column vectors.
*/
  for (attvhdr = consys->attvecs ; attvhdr != NULL ; attvhdr = attvhdr->nxt)
  { 
#   ifdef PARANOIA
    if (!(VALID_ATTVTYPE(attvhdr->what)))
    { errmsg(114,rtnnme,consys->nme,"attached",attvhdr->what) ;
      return (FALSE) ; }
    if (attvhdr->vec == NULL)
    { errmsg(127,rtnnme,consys->nme,attvhdr,
	     consys_assocnme(NULL,attvhdr->what)) ;
      return (FALSE) ; }
#   endif

    if (flgon(attvhdr->what,CONSYS_COLVEC))
      memcpy(((char *) attvhdr->vec)+tndx*attvhdr->elsze,
	     ((char *) attvhdr->vec)+fndx*attvhdr->elsze,attvhdr->elsze) ; } 


  return (TRUE) ; }



static bool find_maxes (consys_struct *consys, bool scan_cols, bool scan_rows)

/*
  This routine scans the constraint matrix for the maximum length column
  and/or row, as indicated.

  Parameters:
    consys:	The constraint matrix.
    scan_cols:	TRUE if the maximum length column needs to be located.
    scan_rows:	TRUE if the maximum length row needs to be located.
  
  Returns: TRUE if the scan completes without error, FALSE otherwise. The
	   routine will not fail unless paranoid checks are enabled.
*/

{ int colndx,rowndx ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "find_maxes" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_COLHDR)) ;
    return (FALSE) ; }
# endif
/*
  Do the scans, as requested.
*/
  if (scan_cols == TRUE)
  { consys->maxcollen = 0 ;
    for (colndx = 1 ; colndx <= consys->varcnt ; colndx++)
    {
#     ifdef PARANOIA
      if (consys->mtx.cols[colndx] == NULL)
      { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
	return (FALSE) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (consys->mtx.cols[colndx]->len == 0 &&
	  flgon(consys->opts,CONSYS_WRNZERO))
	warn(118,rtnnme,
	     consys->nme,"column",consys->mtx.cols[colndx]->nme,colndx) ;
#     endif
      if (consys->mtx.cols[colndx]->len > consys->maxcollen)
      { consys->maxcollen = consys->mtx.cols[colndx]->len ;
	consys->maxcolndx = colndx ; } } }
  if (scan_rows == TRUE)
  { consys->maxrowlen = 0 ;
    for (rowndx = 1 ; rowndx <= consys->concnt ; rowndx++)
    {
#     ifdef PARANOIA
      if (consys->mtx.rows[rowndx] == NULL)
      { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
	return (FALSE) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (consys->mtx.rows[rowndx]->len == 0 &&
	  flgon(consys->opts,CONSYS_WRNZERO))
	warn(118,rtnnme,
	     consys->nme,"row",consys->mtx.rows[rowndx]->nme,rowndx) ;
#     endif
      if (consys->mtx.rows[rowndx]->len > consys->maxrowlen)
      { consys->maxrowlen = consys->mtx.rows[rowndx]->len ;
	consys->maxrowndx = rowndx ; } } }
  
  return (TRUE) ; }



static bool add_logical (consys_struct *consys, int rowndx)

/*
  This routine adds a logical variable appropriate for the constraint
  specified by rowndx. The routine assumes that its caller has arranged for
  the column at rowndx to be empty. It will create a new column header and
  coefficient, and fill in the various associated arrays.

  The conventions are as follows:

    Constraint	Logical		Coefficient	Bounds
      LE	 slack		  1.0		0 <= slack <= +inf
      EQ	 artificial	  1.0		0 <= artificial <= 0
      GE	 surplus	 -1.0		0 <= surplus <= +inf
      RNG	 slack		  1.0		0 <= slack <= (rhs-rhslow)

  To properly set up an artificial variable, the vub array must be present.
  To properly set up a slack for a range constraint, the rhs, rhslow, and vub
  arrays must be present.

  Parameters:
    consys:	constraint system
    rowndx:	constraint receiving the logical variable.

  Returns: TRUE if the logical is successfully created, FALSE if there are
	   problems.
*/

{ int colndx ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  contyp_enum contyp ;
  double val, lb, ub ;
  char *nme ;
  const char *rtnnme = "add_logical" ;

  /* consys_io.c */

  extern char *consys_lognme(consys_struct *consys,
			     int rowndx, char *clientbuf) ;

# ifdef PARANOIA
/*
  The usual, and checks that some necessary associated arrays are present.
  This is an internal routine, so index bounds checks are dropped unless
  we're being paranoid.
*/
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
  if (consys->ctyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CTYP)) ;
    return (FALSE) ; }
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (FALSE) ; }
  if (consys->mtx.rows[rowndx] == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_COLHDR)) ;
    return (FALSE) ; }
  if (consys->vtyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VTYP)) ;
    return (FALSE) ; }
# endif
/*
  Construct a name, coefficient and upper and lower bounds, based on the
  type of constraint.
*/
  nme = consys_lognme(consys,rowndx,NULL) ;
  rowhdr = consys->mtx.rows[rowndx] ;
  val = 1.0 ;
  lb = 0.0 ;
  ub = consys->inf ;
  contyp = consys->ctyp[rowndx] ;
  switch (contyp)
  { case contypLE:
    { break ; }
    case contypEQ:
    { if (consys->vub != NULL)
	ub = 0.0 ;
      else
      { errmsg(120,rtnnme,consys->nme,nme,rowhdr->nme,rowndx) ;
	return (FALSE) ; }
      break ; }
    case contypGE:
    { val = -1.0 ;
      break ; }
    case contypRNG:
    { if (consys->rhs != NULL &&
	  consys->rhslow != NULL && consys->vub != NULL)
      { ub = consys->rhs[rowndx]-consys->rhslow[rowndx] ; }
      else
      { errmsg(120,rtnnme,consys->nme,nme,rowhdr->nme,rowndx) ;
	return (FALSE) ; }
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# ifdef PARANOIA
  if (consys->varcnt+1 > consys->colsze)
  { errmsg(129,rtnnme,consys->nme,nme,rowndx,consys->varcnt+1,
	   "variables",consys->colsze) ;
    return (FALSE) ; }
# endif
/*
  Initialise the column header, install the coefficient in the constraint,
  and set values in the associated arrays.  The previous checks guarantee
  that the arrays are present if they're really needed. You wouldn't think
  we'd need to check maxcollen here, but creating constraint systems in the
  presence of lots of fixed variables can leads to seriously degenerate
  situations.
*/
  colndx = rowndx ;
  colhdr = (colhdr_struct *) CALLOC(1,sizeof(colhdr_struct)) ;
  consys->mtx.cols[colndx] = colhdr ;
  colhdr->ndx = colndx ;
  colhdr->nme = STRALLOC(nme) ;

  coeff = (coeff_struct *) MALLOC(sizeof(coeff_struct)) ;
  coeff->rowhdr = rowhdr ;
  coeff->colhdr = colhdr ;
  coeff->val = val ;

  coeff->colnxt = colhdr->coeffs ;
  colhdr->coeffs = coeff ;
  colhdr->len = 1 ;
  if (consys->maxcollen == 0)
  { consys->maxcollen = 1 ;
    consys->maxcolndx = colndx ; }

  coeff->rownxt = rowhdr->coeffs ;
  rowhdr->coeffs = coeff ;
  rowhdr->len++ ;
  if (rowhdr->len > consys->maxrowlen)
  { consys->maxrowlen = rowhdr->len ;
    consys->maxrowndx = rowndx ; }

  consys->mtx.coeffcnt++ ;
  consys->logvcnt++ ;
  consys->varcnt++ ;

  consys->vtyp[colndx] = vartypCON ;
  if (consys->obj != NULL) consys->obj[colndx] = 0.0 ;
  if (consys->vlb != NULL) consys->vlb[colndx] = 0.0 ;
  if (consys->vub != NULL) consys->vub[colndx] = ub ;

  return (TRUE) ; }



consys_struct *consys_create (const char *nme, flags parts, flags opts,
			      int concnt, int varcnt, double infinity)

/*
  This routine initialises a constraint system data structure, returning a
  pointer to the header. The parameters varcnt and concnt can be used to
  set the initial size of the structure; they are interpreted as follows:
    <= 0: use default sizes
     > 0: use parameter value
  
  Parameters:
    nme:	the name of the constraint system
    parts:	flags indicating which of the associated vectors should be
		created
    opts:	flags indicating which options should be set
    varcnt:	the expected number of variables
    concnt:	the expected number of constraints
    infinity:	the value to be used as infinity

  Returns: pointer to the header.
*/

{ int colsze,rowsze ;
  consys_struct *consys ;
  const char *rtnnme = "consys_create" ;

/*
  Set the row and column sizes.
*/
  if (concnt > 0)
    rowsze = concnt ;
  else
    rowsze = ARCHC_DFLT+CUT_DFLT ;
  if (varcnt > 0)
    colsze = varcnt ;
  else
    colsze = ARCHV_DFLT ;
  if (flgon(opts,CONSYS_LVARS)) colsze += rowsze ;
/*
  Create the header and initialise the non-pointer fields as needed. We're
  counting on calloc to clear the allocated space, so that the parts field
  is clear and the associated vector pointers are null.
*/
  consys = (consys_struct *) CALLOC(1,sizeof(consys_struct)) ; 
  if (nme == NULL)
    consys->nme = STRALLOC("<<unnamed>>") ;
  else
    consys->nme = STRALLOC(nme) ;
  consys->opts = opts ;
  consys->inf = infinity ;
  consys->tiny = 1.0e-20 ;
  if (finite(infinity)) setflg(consys->opts,CONSYS_FININF) ;
  consys->colsze = colsze ;
  consys->rowsze = rowsze ;
/*
  Do up the constraint matrix headers next.
*/
  consys->mtx.cols =
      (colhdr_struct **) CALLOC(colsze+1,sizeof(colhdr_struct *)) ;
  consys->mtx.rows =
      (rowhdr_struct **) CALLOC(rowsze+1,sizeof(rowhdr_struct *)) ;
/*
  Now add any requested associated vectors.
*/

#define CONSYS_MKASSOC(zz_flg_zz,zz_type_zz,zz_ptr_zz) \
  if (flgon(parts,zz_flg_zz)) \
  { if (consys_attach(consys,zz_flg_zz, \
		      sizeof(zz_type_zz),(void **) &zz_ptr_zz) == FALSE) \
    { errmsg(100,rtnnme,consys->nme,consys_assocnme(NULL,zz_flg_zz)) ; \
      return (NULL) ; } \
    setflg(consys->parts,zz_flg_zz) ; \
    clrflg(parts,zz_flg_zz) ; }

  CONSYS_MKASSOC(CONSYS_OBJ,double,consys->obj)
  CONSYS_MKASSOC(CONSYS_VTYP,vartyp_enum,consys->vtyp)
  CONSYS_MKASSOC(CONSYS_VUB,double,consys->vub)
  CONSYS_MKASSOC(CONSYS_VLB,double,consys->vlb)
  CONSYS_MKASSOC(CONSYS_RHS,double,consys->rhs)
  CONSYS_MKASSOC(CONSYS_RHSLOW,double,consys->rhslow)
  CONSYS_MKASSOC(CONSYS_CTYP,contyp_enum,consys->ctyp)
  CONSYS_MKASSOC(CONSYS_CUB,conbnd_struct,consys->cub)
  CONSYS_MKASSOC(CONSYS_CLB,conbnd_struct,consys->clb)
  CONSYS_MKASSOC(CONSYS_RSCALE,double,consys->rowscale)
  CONSYS_MKASSOC(CONSYS_CSCALE,double,consys->colscale)

#undef CONSYS_MKASSOC

# ifndef DYLP_NDEBUG
  if (parts != 0) { errmsg(114,rtnnme,consys->nme,"associated",parts) ; }
# endif

  return (consys) ; }



bool consys_dupsys (consys_struct *src, consys_struct **p_dst, flags dstvecs)

/*
  This routine duplicates a constraint system, including the coefficient
  matrix and those associated vectors specified by dstvecs which exist in
  the source system.

  The algorithm is to first create an empty constraint system with allocated
  size equal to the original. The new header is touched up with information
  from the source. Empty rows are created in the new system by walking the
  row header of the source and calling addrow_pk for each row. The columns
  are filled in using getcol_pk and addcol_pk to copy from the source to the
  destination.  Finally, the specified associated vectors are block copied.

  Parameters:
    src:	the constraint system to be duplicated
    p_dst:	(o) the duplicate constraint system
    dstvecs:	flags indicating the associated vectors that are to be
		duplicated.

  Returns: TRUE if the duplication completes without error, FALSE otherwise.
*/


{ int i,j ;
  char ac ;
  pkvec_struct *pkvec ;

  consys_struct *dst ;

  const char *rtnnme = "consys_dupsys" ;

# ifdef PARANOIA
  if (src == NULL)
  { errmsg(2,rtnnme,"src") ;
    return (FALSE) ; }
  if (p_dst == NULL)
  { errmsg(2,rtnnme,"&dst") ;
    return (FALSE) ; }
# endif

  *p_dst = NULL ;

/*
  For a start, create an empty system with the proper allocated size and
  associated vectors. We're only interested in the vectors that actually
  exist in the source.
*/
  dstvecs = dstvecs & src->parts ;
  dst = consys_create(src->nme,dstvecs,0,src->rowsze,src->colsze,src->inf) ;
  if (dst == NULL)
  { errmsg(152,rtnnme,src->nme) ;
    return (FALSE) ; }
/*
  Copy various header fields from the source system. The reason that src->opts
  is not used to create dst is that it may contain CONSYS_LVARS, which would
  complicate the calculation of allocated size.
*/
  dst->opts = src->opts ;

  dst->maxaij = src->maxaij ;
  dst->minaij = src->minaij ;

  if (src->objnme != NULL) dst->objnme = STRALLOC(src->objnme) ;
  dst->objndx = src->objndx ;
  dst->xzndx = src->xzndx ;
/*
  Walk the row header of src and add empty rows to dst.
*/
  pkvec = pkvec_new(0) ;
  for (i = 1 ; i <= src->concnt ; i++)
  { pkvec->nme = src->mtx.rows[i]->nme ;
    if (i <= src->archccnt)
      ac = 'a' ;
    else
      ac = 'c' ;
    if (consys_addrow_pk(dst,ac,src->ctyp[i],pkvec,0.0,0.0,NULL,NULL) == FALSE)
    { errmsg(156,rtnnme,"row",dst->nme,pkvec->nme) ;
      if (pkvec != NULL) pkvec_free(pkvec) ;
      consys_free(dst) ;
      return (FALSE) ; }
#  ifdef PARANOIA
   if (i != pkvec->ndx)
   { errmsg(136,rtnnme,src->nme,"loss of synch","rows",i,pkvec->ndx) ;
      if (pkvec != NULL) pkvec_free(pkvec) ;
      consys_free(dst) ;
      return (FALSE) ; }
#  endif
  }
  if (pkvec != NULL) pkvec_free(pkvec) ;
/*
  Add the coefficients, column by column.
*/
  pkvec = pkvec_new(src->maxcollen) ;
  for (j = 1 ; j <= src->varcnt ; j++)
  { if (consys_getcol_pk(src,j,&pkvec) == FALSE)
    { errmsg(122,rtnnme,src->nme,"column",consys_nme(src,'v',j,TRUE,NULL),j) ;
      if (pkvec != NULL) pkvec_free(pkvec) ;
      consys_free(dst) ;
      return (FALSE) ; }
    if (consys_addcol_pk(dst,src->vtyp[j],pkvec,0.0,0.0,0.0) == FALSE)
    { errmsg(156,rtnnme,"column",dst->nme,pkvec->nme) ;
      if (pkvec != NULL) pkvec_free(pkvec) ;
      consys_free(dst) ;
      return (FALSE) ; }
#   ifdef PARANOIA
    if (j != pkvec->ndx)
    { errmsg(136,rtnnme,src->nme,"loss of synch","columns",i,pkvec->ndx) ;
      if (pkvec != NULL) pkvec_free(pkvec) ;
      consys_free(dst) ;
      return (FALSE) ; }
#   endif
  }
  if (pkvec != NULL) pkvec_free(pkvec) ;
# ifdef PARANOIA
  if (src->mtx.coeffcnt != dst->mtx.coeffcnt)
  { errmsg(136,rtnnme,src->nme,"coefficient count mismatch","columns",
	   src->mtx.coeffcnt,dst->mtx.coeffcnt) ;
    consys_free(dst) ;
    return (FALSE) ; }
# endif
/*
  Force copy to use the same maximum row and column indices as the original.
  (The difference arises from insertion order of columns.)
*/
  dst->maxrowndx = src->maxrowndx ;
  dst->maxcolndx = src->maxcolndx ;
/*
  Copy the associated vectors.
*/
  if (flgon(dstvecs,CONSYS_OBJ))
    memcpy(dst->obj,src->obj,(src->varcnt+1)*sizeof(double)) ;
  if (flgon(dstvecs,CONSYS_VUB))
    memcpy(dst->vub,src->vub,(src->varcnt+1)*sizeof(double)) ;
  if (flgon(dstvecs,CONSYS_VLB))
    memcpy(dst->vlb,src->vlb,(src->varcnt+1)*sizeof(double)) ;
  if (flgon(dstvecs,CONSYS_VTYP))
    memcpy(dst->vtyp,src->vtyp,(src->varcnt+1)*sizeof(vartyp_enum)) ;
  if (flgon(dstvecs,CONSYS_CSCALE))
    memcpy(dst->colscale,src->colscale,(src->varcnt+1)*sizeof(double)) ;

  if (flgon(dstvecs,CONSYS_RHS))
    memcpy(dst->rhs,src->rhs,(src->concnt+1)*sizeof(double)) ;
  if (flgon(dstvecs,CONSYS_RHSLOW))
    memcpy(dst->rhslow,src->rhslow,(src->concnt+1)*sizeof(double)) ;
  if (flgon(dstvecs,CONSYS_CUB))
    memcpy(dst->cub,src->cub,(src->concnt+1)*sizeof(conbnd_struct)) ;
  if (flgon(dstvecs,CONSYS_CLB))
    memcpy(dst->clb,src->clb,(src->concnt+1)*sizeof(conbnd_struct)) ;
  if (flgon(dstvecs,CONSYS_CTYP))
    memcpy(dst->ctyp,src->ctyp,(src->concnt+1)*sizeof(contyp_enum)) ;
  if (flgon(dstvecs,CONSYS_RSCALE))
    memcpy(dst->rowscale,src->rowscale,(src->concnt+1)*sizeof(double)) ;

  *p_dst = dst ;

  return (TRUE) ; }



void consys_free (consys_struct *consys)

/*
  This routine frees the space allocated for a constraint system. Sigh. The
  trouble with complex data structures is they're a pain to disassemble.

  Parameter:
    consys:	constraint system
  
  Returns: undefined
*/

{ int ndx ;
  attvhdr_struct *attvhdr ;
  lnk_struct *lnk ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# ifdef PARANOIA

  const char *rtnnme = "consys_free" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return ; }
# endif
/*
  Start by dismantling the attached vector list.
*/
  for (attvhdr = consys->attvecs ; attvhdr != NULL ; attvhdr = consys->attvecs)
  { consys->attvecs = attvhdr->nxt ;
    for (lnk = attvhdr->pveclst ; lnk != NULL ; lnk = attvhdr->pveclst)
    { attvhdr->pveclst = lnk->llnxt ;
      FREE(lnk) ; }
    FREE(attvhdr) ; }
/*
  Take care of the associated vectors.
*/
  if (consys->obj != NULL) FREE(consys->obj) ;
  if (consys->vtyp != NULL) FREE(consys->vtyp) ;
  if (consys->vub != NULL) FREE(consys->vub) ;
  if (consys->vlb != NULL) FREE(consys->vlb) ;
  if (consys->rhs != NULL) FREE(consys->rhs) ;
  if (consys->rhslow != NULL) FREE(consys->rhslow) ;
  if (consys->ctyp != NULL) FREE(consys->ctyp) ;
  if (consys->cub != NULL) FREE(consys->cub) ;
  if (consys->clb != NULL) FREE(consys->clb) ;
  if (consys->rowscale != NULL) FREE(consys->rowscale) ;
  if (consys->colscale != NULL) FREE(consys->colscale) ;
/*
  Now dismantle the constraint matrix. Walk the column header array, freeing
  the coefficients and header for each column. Then walk the row header
  array, freeing the header for each row. Finally, free the column and row
  header arrays and the matrix header.
*/
  for (ndx = 1 ; ndx <= consys->varcnt ; ndx++)
  { colhdr = consys->mtx.cols[ndx] ;
    for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = colhdr->coeffs)
    { colhdr->coeffs = coeff->colnxt ;
      FREE(coeff) ; }
    if (colhdr->nme != NULL) STRFREE(colhdr->nme) ;
    FREE(colhdr) ; }
  for (ndx = 1 ; ndx <= consys->concnt ; ndx++)
  { rowhdr = consys->mtx.rows[ndx] ;
    if (rowhdr->nme != NULL) STRFREE(rowhdr->nme) ;
    FREE(rowhdr) ; }
  FREE(consys->mtx.cols) ;
  FREE(consys->mtx.rows) ;
/*
  Finish up with the constraint system header.
*/
  if (consys->objnme != NULL) STRFREE(consys->objnme) ;
  if (consys->nme != NULL) STRFREE(consys->nme) ;
  FREE(consys) ;

  return ; }



bool consys_attach (consys_struct *consys, flags what, int elsze, void **pvec)

/*
  This routine 'attaches' a vector (more specifically, a reference to it) to
  the constraint system so that the vector and references will be
  automatically resized when the constraint system is resized. If *pvec
  is null, a vector of appropriate size is created. If the vector
  pointed to by pvec isn't already on the attvecs list, a new entry is
  created. Then pvec is added to the list of pointers referring to the
  vector.

  The expectation is that this list will never become very big, so there's no
  particular attempt to keep the entries in order.

  If DYLP_NDEBUG is not set and the WRNATT flag is set, a warning will be issued if
  the attach request specifies a (pointer,vector) pair already attached to the
  system.

  Parameters:
    consys:	constraint system
    what:	type of vector being attached
    elsze:	the size of an element in the vector
    pvec:	(i) address of pointer to the vector; if *pvec == NULL, a
		    vector will be created
		(o) *pvec holds address of vector

  Returns: TRUE if the attachment is successful, FALSE otherwise (paranoia
	   only).
*/

{ int ndx ;
  attvhdr_struct *attvhdr ;
  lnk_struct *lnk ;
  double *dvec ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_attach" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (pvec == NULL)
  { errmsg(2,rtnnme,"&vec") ;
    return (FALSE) ; }
  if (!(VALID_ATTVTYPE(what)))
  { errmsg(6,rtnnme,"vector type",what) ;
    return (FALSE) ; }
  if (elsze <= 0)
  { errmsg(5,rtnnme,"element size",elsze) ;
    return (FALSE) ; }
# endif
/*
  Have a look at vec, and create one if we need to. Otherwise, look through
  the list of attached vectors to see if there's already a header for this
  vector (and we're just adding another pointer to it). Remember, colsze is
  the allocated column capacity, hence the size of a row; similarly for
  rowsze.
*/
  if (*pvec == NULL)
  { if (flgon(what,CONSYS_ROWVEC))
    { *pvec = (void *) CALLOC(consys->colsze+1,elsze) ; }
    else
    { *pvec = (void *) CALLOC(consys->rowsze+1,elsze) ; }
    if (flgon(what,CONSYS_VUB))
    { dvec = (double *) *pvec ;
      dvec[0] = 0.0 ;
      for (ndx = 1 ; ndx <= consys->colsze ; ndx++)
      { dvec[ndx] = consys->inf ; } }
    else
    if (flgon(what,CONSYS_CSCALE))
    { dvec = (double *) *pvec ;
      dvec[0] = 0.0 ;
      for (ndx = 1 ; ndx <= consys->rowsze ; ndx++) dvec[ndx] = 1.0 ; }
    else
    if (flgon(what,CONSYS_RSCALE))
    { dvec = (double *) *pvec ;
      dvec[0] = 0.0 ;
      for (ndx = 1 ; ndx <= consys->colsze ; ndx++) dvec[ndx] = 1.0 ; }
    attvhdr = NULL ; }
  else
  { for (attvhdr = consys->attvecs ; attvhdr != NULL ; attvhdr = attvhdr->nxt)
      if (attvhdr->vec == *pvec) break ; }
/*
  If we didn't find a header, create one now.
*/
  if (attvhdr == NULL)
  { attvhdr = (attvhdr_struct *) MALLOC(sizeof(attvhdr_struct)) ;
    attvhdr->what = what ;
    attvhdr->elsze = elsze ;
    attvhdr->vec = *pvec ;
    attvhdr->pveclst = NULL ;
    attvhdr->nxt = consys->attvecs ;
    consys->attvecs = attvhdr ; }
# ifdef PARANOIA
  else
  { if (attvhdr->what != what)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; } }
# endif
/*
  Make an entry for pvec. If there's already an entry, consider warning the
  user.
*/
  for (lnk = attvhdr->pveclst ; lnk != NULL ; lnk = lnk->llnxt)
    if (pvec == (void **) lnk->llval) break ;
  if (lnk == NULL)
  { lnk = (lnk_struct *) MALLOC(sizeof(lnk_struct)) ;
    lnk->llval = (void *) pvec ;
    lnk->llnxt = attvhdr->pveclst ;
    attvhdr->pveclst = lnk ; }
# ifndef DYLP_NDEBUG
  else
  { if (flgon(consys->opts,CONSYS_WRNATT))
      warn(107,rtnnme,consys_assocnme(consys,what),*pvec,pvec) ; }
# endif
  
  return (TRUE) ; }



bool consys_update (consys_struct *consys, void *old, void *new)

/*
  This routine looks up the entry for the vector specified by old and updates
  all the references with the vector in new. In effect, the old vector is
  detached from the constraint system and the new one attached in its place.

  It's the user's responsibility to make sure that some reference to old is
  retained, should that be necessary.

  Parameters:
    consys:	constraint system
    old:	old vector
    new:	new vector

  Returns: TRUE if the update is successful, FALSE otherwise.
*/

{ attvhdr_struct *attvhdr ;
  lnk_struct *lnk ;
  const char *rtnnme = "consys_update" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (old == NULL)
  { errmsg(2,rtnnme,"old") ;
    return (FALSE) ; }
  if (new == NULL)
  { errmsg(2,rtnnme,"new") ;
    return (FALSE) ; }
# endif
/*
  Look for an entry for old. It's an error if we can't find it.
*/
  for (attvhdr = consys->attvecs ; attvhdr != NULL ; attvhdr = attvhdr->nxt)
    if (attvhdr->vec == old) break ;
  if (attvhdr == NULL)
  { setflg(consys->opts,CONSYS_CORRUPT) ;
    errmsg(104,rtnnme,consys->nme,old) ;
    return (FALSE) ; }
/*
  Do the update.
*/
# ifdef PARANOIA
  if (attvhdr->pveclst == NULL)
  { errmsg(108,rtnnme,consys_assocnme(NULL,attvhdr->what),attvhdr->vec) ;
    return (FALSE) ; }
# endif
  attvhdr->vec = new ;
  for (lnk = attvhdr->pveclst ; lnk != NULL ; lnk = lnk->llnxt)
  {
# ifdef PARANOIA
    if (lnk->llval == NULL)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; }
# endif
    *((void **) lnk->llval) = new ; }

  return (TRUE) ; }



bool consys_detach (consys_struct *consys, void **pvec, bool all)

/*
  This routine detaches the reference specified by pvec by removing it from
  the attached vector list. If all == TRUE, all references to the vector
  *pvec are removed and the vector is removed from the attached vector list.

  It's an error if pvec doesn't show up on the reference list for the vector.

  Parameters:
    consys:	constraint system
    pvec:	address of pointer to vector; *pvec must also be valid
    all:	TRUE to detach all references to this vector except the
		associated vector entry, FALSE if only pvec is to be
		detached.

  Returns: TRUE if the entry is found and detached, FALSE otherwise
*/

{ attvhdr_struct *attvhdr,**pattvhdr ;
  lnk_struct *lnk,**plnk ;
  bool pvec_seen ;
  void *vec ;
  const char *rtnnme = "consys_detach" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->attvecs == NULL)
  { errmsg(101,rtnnme,consys->nme,"attached vector list") ;
    return (FALSE) ; }
  if (pvec == NULL)
  { errmsg(2,rtnnme,"&vec") ;
    return (FALSE) ; }
  if (*pvec == NULL)
  { errmsg(2,rtnnme,"vec") ;
    return (FALSE) ; }
# endif
/*
  Find the entry for the vector. It's an error to try and detach a vector that
  isn't attached.
*/
  vec = *pvec ;
  for (pattvhdr = &consys->attvecs, attvhdr = *pattvhdr ;
       attvhdr != NULL ;
       pattvhdr = &attvhdr->nxt, attvhdr = *pattvhdr)
    if (attvhdr->vec == vec) break ;
  if (attvhdr == NULL)
  { setflg(consys->opts,CONSYS_CORRUPT) ;
    errmsg(104,rtnnme,consys->nme,vec) ;
    return (FALSE) ; }
# ifdef PARANOIA
  if (attvhdr->pveclst == NULL)
  { errmsg(108,rtnnme,consys_assocnme(NULL,attvhdr->what),attvhdr->vec) ;
    return (FALSE) ; }
# endif
/*
  Work over the pvec list, detaching entries as appropriate. If there are no
  references left, delete the vector's entry too.
*/
  if (all == TRUE)
  { for (lnk = attvhdr->pveclst ; lnk != NULL ; lnk = attvhdr->pveclst)
    { attvhdr->pveclst = lnk->llnxt ;
      FREE(lnk) ; }
    *pattvhdr = attvhdr->nxt ;
    FREE(attvhdr) ; }
  else
  { pvec_seen = FALSE ;
    plnk = &attvhdr->pveclst ;
    for (lnk = *plnk ; lnk != NULL ; lnk = *plnk)
    { if (pvec == (void **) lnk->llval)
      { pvec_seen = TRUE ;
	*plnk = lnk->llnxt ;
	FREE(lnk) ; }
      else
      { plnk = &lnk->llnxt ; } }
    if (pvec_seen == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(109,rtnnme,consys->nme,pvec,consys_assocnme(NULL,attvhdr->what),
	     attvhdr->vec) ;
      return (FALSE) ; }
    if (attvhdr->pveclst == NULL)
    { *pattvhdr = attvhdr->nxt ;
      FREE(attvhdr) ; } }

  return (TRUE) ; }



bool consys_realloc (consys_struct *consys, char rowcol, int incr)

/*
  This routine is responsible for expanding the constraint matrix and the
  attached arrays when additional room is needed to add rows or columns.
  It does nothing except create open space at the end of arrays (increasing
  colsze or rowsze).
  
  Expanding the column capacity entails expanding the column header of the
  constraint matrix, and any attached row vectors.
  
  Expanding the row capacity entails expanding the row header of the
  constraint matrix, and any attached column vectors. If coupling is
  enabled, a column expansion is forced to make room for the logical
  variables that will go with the new rows.

  In practice, it's all pretty straightforward: The row and/or column headers
  are expanded as needed, then we walk the attached vector list, expanding as
  needed and updating pointers when the expanded vector moves in memory.

  If incr is <= 0, a default expansion of 10% (minimum 10 rows/columns) is
  used.

  Parameters:
    consys:	Constraint system.
    rowcol:	'c' for column expansion, 'r' for row expansion
    incr:	The number of additional rows/columns.

  Returns: TRUE if the expansion is successful, FALSE otherwise.
*/

{ int rowsze,oldrowsze,colsze,oldcolsze,ndx,elsze ;
  bool expcols,exprows ;
  flags what ;
  char *vec ;
  double *dvec ;
  attvhdr_struct *attvhdr ;
  lnk_struct *lnk ;
  const char *rtnnme = "consys_realloc" ;

/*
  Suppress `might be used uninitialized' from gcc.
*/
  rowsze = -1 ;
  oldrowsze = -1 ;
  colsze = -1 ;
  oldcolsze = -1 ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
# endif
/*
  Figure out the size of the increase, if the user didn't specify. If we're
  asked to do a row increase and coupling is enabled, check that the capacity
  remaining on the column side is sufficient to accommodate the additional
  logicals. If not, calculate the needed increase and force a column
  expansion.
*/
  expcols = FALSE ;
  exprows = FALSE ;
  switch (rowcol)
  { case 'c':
    { if (incr <= 0) incr = maxx(10,(int)(consys->colsze*EXPAND_PCT_DFLT)) ;
      colsze = consys->colsze+incr ;
      expcols = TRUE ;
      break ; }
    case 'r':
    { if (incr <= 0) incr = maxx(10,(int)(consys->rowsze*EXPAND_PCT_DFLT)) ;
      rowsze = consys->rowsze+incr ;
      exprows = TRUE ;
      if (flgon(consys->opts,CONSYS_LVARS))
	if (consys->colsze < consys->varcnt+incr)
	{ expcols = TRUE ;
	  colsze = consys->varcnt+incr ; }
      break ; }
    default:
    { errmsg(3,rtnnme,"rowcol",rowcol) ;
      return (FALSE) ; } }
/*
  Expand the row and column header as necessary.
*/
  if (exprows == TRUE)
  { oldrowsze = consys->rowsze ;
    consys->mtx.rows = (rowhdr_struct **)
		REALLOC(consys->mtx.rows,(rowsze+1)*sizeof(rowhdr_struct **)) ;
    memset(&consys->mtx.rows[oldrowsze+1],0,
	   (rowsze-oldrowsze)*sizeof(rowhdr_struct **)) ;
    consys->rowsze = rowsze ; }
  if (expcols == TRUE)
  { oldcolsze = consys->colsze ;
    consys->mtx.cols = (colhdr_struct **)
		REALLOC(consys->mtx.cols,(colsze+1)*sizeof(colhdr_struct **)) ;
    memset(&consys->mtx.cols[oldcolsze+1],0,
	   (colsze-oldcolsze)*sizeof(colhdr_struct **)) ;
    consys->colsze = colsze ; }
/*
  Now, walk the attached vector list and do any necessary expansions. Where
  the expansion results in the vector moving to a different block of memory,
  walk the reference list and update the pointers.  Note that we don't care
  why vec didn't change -- it could be that the vector wasn't expanded at
  all, or that it was expanded but didn't move.
*/
  for (attvhdr = consys->attvecs ; attvhdr != NULL ; attvhdr = attvhdr->nxt)
  { vec = (char *) attvhdr->vec ;
    what = attvhdr->what ;
    elsze = attvhdr->elsze ;
#   ifdef PARANOIA
    if (elsze <= 0)
    { errmsg(106,rtnnme,consys->nme,vec,elsze) ;
      return (FALSE) ; }
    if (attvhdr->pveclst == NULL)
    { errmsg(108,rtnnme,consys_assocnme(NULL,attvhdr->what),attvhdr->vec) ;
      return (FALSE) ; }
    if (!(VALID_ATTVTYPE(what)))
    { errmsg(114,rtnnme,consys->nme,"attached",attvhdr->what) ;
      return (FALSE) ; }
#   endif
    if (flgon(what,CONSYS_ROWVEC))
    { if (expcols == TRUE)
      { vec = (char *) REALLOC(attvhdr->vec,(colsze+1)*elsze) ;
	if (flgon(what,CONSYS_VUB))
	{ dvec = (double *) vec ;
	  for (ndx = oldcolsze+1 ; ndx <= colsze ; ndx++)
	  { dvec[ndx] = consys->inf ; } }
	else
	if (flgon(what,CONSYS_CSCALE))
	{ dvec = (double *) vec ;
	  dvec[0] = 0.0 ;
	  for (ndx = oldcolsze+1 ; ndx <= colsze ; ndx++) dvec[ndx] = 1.0 ; }
	else
	{ memset(vec+(oldcolsze+1)*elsze,0,(colsze-oldcolsze)*elsze) ; } } }
    else
    { if (exprows == TRUE)
      { vec = (char *) REALLOC(attvhdr->vec,(rowsze+1)*elsze) ;
	if (flgon(what,CONSYS_RSCALE))
	{ dvec = (double *) vec ;
	  dvec[0] = 0.0 ;
	  for (ndx = oldrowsze+1 ; ndx <= rowsze ; ndx++) dvec[ndx] = 1.0 ; }
	else
	{ memset(vec+(oldrowsze+1)*elsze,0,(rowsze-oldrowsze)*elsze) ; } } }
    if (vec != attvhdr->vec)
    { for (lnk = attvhdr->pveclst ; lnk != NULL ; lnk = lnk->llnxt)
      {
#       ifdef PARANOIA
	if (lnk->llval == NULL)
	{ errmsg(110,rtnnme,lnk,consys_assocnme(NULL,what),attvhdr->vec) ;
	  return (FALSE) ; }
#       endif
	*((void **) lnk->llval) = (void *) vec ; }
      attvhdr->vec = (void *) vec ; } }

# ifdef PARANOIA
/*
  Just checking to see if we dropped anything.
*/
  if (flgon(consys->parts,CONSYS_RHS) && consys->rhs == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_RHS)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_RHSLOW) && consys->rhslow == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_RHSLOW)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_CTYP) && consys->ctyp == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CTYP)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_CUB) && consys->cub == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CUB)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_CLB) && consys->clb == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CLB)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_OBJ) && consys->obj == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_OBJ)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_VTYP) && consys->vtyp == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VTYP)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_VUB) && consys->vub == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VUB)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_VLB) && consys->vlb == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VLB)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_RSCALE) && consys->rowscale == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_RSCALE)) ;
    return (FALSE) ; }
  if (flgon(consys->parts,CONSYS_CSCALE) && consys->rowscale == NULL)
  { errmsg(113,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CSCALE)) ;
    return (FALSE) ; }
# endif

  return (TRUE) ; }



bool consys_addcol_pk (consys_struct *consys,
		       vartyp_enum vartyp, pkvec_struct *pkcol,
		       double obj, double vlb, double vub)

/*
  This routine adds an architectural variable to the constraint matrix
  (expanding the allocated size if necessary). If the appropriate associated
  vectors exist, the objective coefficient and lower and upper bounds will
  be set.

  If DYLP_NDEBUG is not defined and the WRNZERO flag is set, the routine will
  warn about columns with no non-zero coefficients. NO WARNING is issued if
  an associated vector doesn't exist.

  If pkcol is a vector of capacity 0 (pkcol->sze == 0) this is taken
  to be a call for the sole purpose of creating the header, and no warning
  is issued.

  Parameters:
    consys:	Constraint system.
    vartyp:	The type of variable.
    pkcol:	(i) name (optional) and coefficients for the column.
		(o) column index; name, if not supplied.
    obj:	Coefficient for the objective function.
    vlb, vub	Lower and upper bounds.

  Returns: TRUE if the insertion is successful, FALSE otherwise.
*/

{ int colndx,vecndx,avail,nzcnt ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  pkcoeff_struct *pkcoeff ;
  char nmebuf[20] ;
  const char *rtnnme = "consys_addcol_pk" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (consys->vtyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VTYP)) ;
    return (FALSE) ; }
  if (!VALID_VARTYPE(vartyp))
  { errmsg(5,rtnnme,"vartyp",(int) vartyp) ;
    return (FALSE) ; }
  if (pkcol == NULL)
  { errmsg(2,rtnnme,"pkcol") ;
    return (FALSE) ; }
  pkcol->ndx = 0 ;
# endif
/*
  Calculate the column index (and a name, if necessary). Then make sure we
  have space; acquire some if needed.
*/
  colndx = consys->varcnt+1 ;
  pkcol->ndx = colndx ;
  if (pkcol->nme == NULL)
  { dyio_outfxd(nmebuf,-((int) (sizeof(nmebuf)-1)),'l',"var<%d>",colndx) ;
    pkcol->nme = nmebuf ; }
# ifdef PARANOIA
/*
  We need to postpone this check to here because pkvec_check requires that
  nme != NULL.
*/
  if (pkvec_check(pkcol,rtnnme) == FALSE) return (FALSE) ;
# endif
  if (flgon(consys->opts,CONSYS_LVARS))
    avail = consys->colsze-consys->rowsze ;
  else
    avail = consys->colsze ;
  if (avail < consys->varcnt+1)
    if (consys_realloc(consys,'c',0) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,
	     "capacity expansion","column",pkcol->nme,pkcol->ndx) ;
      return (FALSE) ; }
/*
  Initialise the column header and set the type. Adjust the variable counts.
*/
  colhdr = (colhdr_struct *) CALLOC(1,sizeof(colhdr_struct)) ;
  consys->mtx.cols[colndx] = colhdr ;
  colhdr->ndx = colndx ;
  colhdr->nme = STRALLOC(pkcol->nme) ;
  consys->vtyp[colndx] = vartyp ;
  consys->archvcnt++ ;
  consys->varcnt++ ;
  if (vartyp == vartypINT)
    consys->intvcnt++ ;
  else
  if (vartyp == vartypBIN)
    consys->binvcnt++ ;
/*
  If we're using our private buffer, hide it again.
*/
  if (pkcol->nme == nmebuf) pkcol->nme = colhdr->nme ;
/*
  Add the column, keeping an eye out for any changes in row or column maxima.
  It's an error to try and insert a coefficient into a constraint that's not
  already in existence. If pkcol->sze is 0, then there are no coefficients
  (the purpose of the call was to create the header).
*/
  if (pkcol->sze > 0)
  {
#   ifndef DYLP_NDEBUG
    if (pkcol->cnt == 0 && flgon(consys->opts,CONSYS_WRNZERO))
      warn(118,rtnnme,consys->nme,"column",colhdr->nme,colndx) ;
#   endif
    nzcnt = 0 ;
    pkcoeff = pkcol->coeffs ;
    for (vecndx = 0 ; vecndx < pkcol->cnt ; vecndx++)
    { if (pkcoeff->ndx <= 0 || pkcoeff->ndx > consys->concnt)
      { errmsg(102,rtnnme,consys->nme,"row",pkcoeff->ndx,1,consys->concnt) ;
	return (FALSE) ; }
      if (fabs(pkcoeff->val) >= consys->inf)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(128,rtnnme,consys->nme,pkcoeff->ndx,colndx,pkcoeff->val,
	       "column",colhdr->nme) ;
	return (FALSE) ; }
      if (fabs(pkcoeff->val) > consys->tiny)
      { rowhdr = consys->mtx.rows[pkcoeff->ndx] ;
#       ifdef PARANOIA
	if (rowhdr == NULL)
	{ errmsg(103,rtnnme,consys->nme,"row",pkcoeff->ndx) ;
	  return (FALSE) ; }
	if (pkcoeff->ndx != rowhdr->ndx)
	{ errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,
		 pkcoeff->ndx,rowhdr) ;
	  return (FALSE) ; }
#       endif
	coeff = (coeff_struct *) MALLOC(sizeof(coeff_struct)) ;
	coeff->rowhdr = rowhdr ;
	coeff->colhdr = colhdr ;
	coeff->val = pkcoeff->val ;
	coeff->rownxt = rowhdr->coeffs ;
	rowhdr->coeffs = coeff ;
	coeff->colnxt = colhdr->coeffs ;
	colhdr->coeffs = coeff ;
	rowhdr->len++ ;
	nzcnt++ ;
	if (rowhdr->len > consys->maxrowlen)
	{ consys->maxrowlen = rowhdr->len ;
	  consys->maxrowndx = pkcoeff->ndx ; } }
#     ifndef DYLP_NDEBUG
      else
      { warn(130,rtnnme,consys->nme,pkcoeff->ndx,colndx,pkcoeff->val,
	     consys->tiny,"column",colhdr->nme) ; }
#     endif
      pkcoeff++ ; }

    colhdr->len = nzcnt ;
    consys->mtx.coeffcnt += nzcnt ;
    if (colhdr->len > consys->maxcollen)
    { consys->maxcollen = colhdr->len ;
      consys->maxcolndx = colndx ; } }
/*
  Check for the relevant associated arrays -- obj, vlb, vub -- and set values
  if they're present.
*/
  if (consys->obj != NULL) consys->obj[colndx] = obj ;
  if (consys->vlb != NULL) consys->vlb[colndx] = vlb ;
  if (consys->vub != NULL) consys->vub[colndx] = vub ;

  return (TRUE) ; }



bool consys_addcol_ex (consys_struct *consys,
		       vartyp_enum vartyp, const char **nme, double *excol,
		       double obj, double vlb, double vub)

/*
  This routine adds an architectural variable to the constraint matrix
  (expanding the allocated size if necessary). If the appropriate associated
  vectors exist, the objective coefficient and lower and upper bounds will
  be set.

  It is assumed that excol is of size consys->concnt. Only non-zero entries
  are installed in the constraint system; zeros are quietly suppressed.

  Parameters:
    consys:	Constraint system.
    vartyp:	The type of variable.
    nme:	(i) The name for the variable (will be generated if null)
		(o) The name for the variable
    excol:	(i) coefficients for the column.
		(o) excol[0] is set to the column index.
    obj:	Coefficient for the objective function.
    vlb, vub	Lower and upper bounds.

  Returns: TRUE if the insertion is successful, FALSE otherwise.
*/

{ int colndx,rowndx,avail ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  char nmebuf[20] ;
  const char *rtnnme = "consys_addcol_ex" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (consys->vtyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VTYP)) ;
    return (FALSE) ; }
  if (!VALID_VARTYPE(vartyp))
  { errmsg(5,rtnnme,"vartyp",(int) vartyp) ;
    return (FALSE) ; }
  if (excol == NULL)
  { errmsg(2,rtnnme,"excol") ;
    return (FALSE) ; }
  if (nme == NULL)
  { errmsg(2,rtnnme,"nme") ;
    return (FALSE) ; }
# endif
/*
  Calculate the column index (and a name, if necessary). Then make sure we
  have space; acquire some if needed.
*/
  colndx = consys->varcnt+1 ;
  excol[0] = colndx ;
  if (nme == NULL)
  { dyio_outfxd(nmebuf,-((int) (sizeof(nmebuf)-1)),'l',"var<%d>",colndx) ;
    *nme = nmebuf ; }
  if (flgon(consys->opts,CONSYS_LVARS))
    avail = consys->colsze-consys->rowsze ;
  else
    avail = consys->colsze ;
  if (avail < consys->varcnt+1)
    if (consys_realloc(consys,'c',0) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,
	     "capacity expansion","column",*nme,excol[0]) ;
      return (FALSE) ; }
/*
  Initialise the column header and set the type. Adjust the variable counts.
*/
  colhdr = (colhdr_struct *) CALLOC(1,sizeof(colhdr_struct)) ;
  consys->mtx.cols[colndx] = colhdr ;
  colhdr->ndx = colndx ;
  colhdr->nme = STRALLOC(*nme) ;
  consys->vtyp[colndx] = vartyp ;
  consys->archvcnt++ ;
  consys->varcnt++ ;
  if (vartyp == vartypINT)
    consys->intvcnt++ ;
  else
  if (vartyp == vartypBIN)
    consys->binvcnt++ ;
/*
  If we're using our private buffer, hide it again.
*/
  if (*nme == nmebuf) *nme = colhdr->nme ;
/*
  Add the column, keeping an eye out for any changes in row or column maxima.
  It's an error to try and insert a coefficient into a constraint that's not
  already in existence.
*/
  for (rowndx = 1 ; rowndx <= consys->concnt ; rowndx++)
  { if (fabs(excol[rowndx]) >= consys->inf)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(128,rtnnme,consys->nme,rowndx,colndx,excol[rowndx],
	     "column",colhdr->nme) ;
      return (FALSE) ; }
    if (fabs(excol[rowndx]) >= consys->tiny)
    { colhdr->len++ ;
      rowhdr = consys->mtx.rows[rowndx] ;
#     ifdef PARANOIA
      if (rowhdr == NULL)
      { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
	return (FALSE) ; }
      if (rowndx != rowhdr->ndx)
      { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,rowhdr) ;
	return (FALSE) ; }
#     endif
      coeff = (coeff_struct *) MALLOC(sizeof(coeff_struct)) ;
      coeff->rowhdr = rowhdr ;
      coeff->colhdr = colhdr ;
      coeff->val = excol[rowndx] ;
      coeff->rownxt = rowhdr->coeffs ;
      rowhdr->coeffs = coeff ;
      coeff->colnxt = colhdr->coeffs ;
      colhdr->coeffs = coeff ;
      rowhdr->len++ ;
      if (rowhdr->len > consys->maxrowlen)
      { consys->maxrowlen = rowhdr->len ;
	consys->maxrowndx = rowndx ; } }
#   ifndef DYLP_NDEBUG
    else
    if (excol[rowndx] != 0.0)
    { warn(130,rtnnme,consys->nme,rowndx,colndx,excol[rowndx],
	   consys->tiny,"row",rowhdr->nme) ; }
#   endif
  }

  consys->mtx.coeffcnt += colhdr->len ;
  if (colhdr->len > consys->maxcollen)
  { consys->maxcollen = colhdr->len ;
    consys->maxcolndx = colndx ; }
/*
  Check for the relevant associated arrays -- obj, vlb, vub -- and set values
  if they're present.
*/
  if (consys->obj != NULL) consys->obj[colndx] = obj ;
  if (consys->vlb != NULL) consys->vlb[colndx] = vlb ;
  if (consys->vub != NULL) consys->vub[colndx] = vub ;

  return (TRUE) ; }



bool consys_addrow_pk (consys_struct *consys, char class,
		       contyp_enum contyp, pkvec_struct *pkrow,
		       double rhs, double rhslow,
		       conbnd_struct *cub, conbnd_struct *clb)

/*
  This routine adds a constraint to the constraint matrix, expanding the
  allocated size if necessary and reordering constraints as appropriate.  If
  the relevant associated vectors exist, the right-hand-side, range (rhslow),
  and upper and lower bounds will be set.

  If coupling is on (i.e., logical variables are enabled), a logical variable
  will be created according to the usual rules and any reordering of
  constraints will be reflected in the logicals. IT'S IMPORTANT that the
  logical not be added until after pkrow is unloaded -- after all, pkrow
  could have a coefficient in the column occupied by the first architectural
  variable, and we don't want to screw up the column index stored with the
  coefficient.

  If DYLP_NDEBUG is not defined and the WRNZERO flag is set, the routine will
  warn about columns with no non-zero coefficients. NO WARNING is issued if
  an associated vector doesn't exist.

  Parameters:
    consys:	Constraint system.
    class:	Class of the constraint -- 'a' for architectural, 'c' for cut
    contyp:	The type of constraint.
    pkrow:	(i) name (optional) and coefficients for the row
		(o) row index; name, if not supplied
    rhs:	right-hand-side coefficient
    rhslow:	alternate right-hand-side coefficient, for range constraints
    cub,clb:	upper and lower bounds for the value of the left-hand-side

  Returns: TRUE if the insertion is successful, FALSE otherwise.
*/

{ int rowndx,vecndx,nzcnt ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  pkcoeff_struct *pkcoeff ;
  char nmebuf[20] ;

  const char *rtnnme = "consys_addrow_pk" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (!(class == 'a' || class == 'c'))
  { errmsg(3,rtnnme,"constraint class",class) ;
    return (FALSE) ; }
  if (consys->ctyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CTYP)) ;
    return (FALSE) ; }
  if (!(VALID_CONTYPE(contyp) || contyp == contypNB))
  { errmsg(5,rtnnme,"contyp",(int) contyp) ;
    return (FALSE) ; }
  if (pkrow == NULL)
  { errmsg(2,rtnnme,"pkrow") ;
    return (FALSE) ; }
  pkrow->ndx = 0 ;
  if (flgon(consys->opts,CONSYS_LVARS))
    if (consys->logvcnt != consys->concnt)
    { errmsg(131,rtnnme,consys->nme,consys->logvcnt,consys->concnt) ;
      return (FALSE) ; }
# endif
/*
  Calculate the row index and name.  A new architectural constraint will be
  placed at the end of the block of architectural constraints, and a new cut
  constraint will be placed at the end of the rows.
*/
  if (class == 'a')
    rowndx = consys->archccnt+1 ;
  else
    rowndx = consys->concnt+1 ;
  pkrow->ndx = rowndx ;
  if (pkrow->nme == NULL)
  { dyio_outfxd(nmebuf,-((int) (sizeof(nmebuf)-1)),'l',"%s<%d>",
	        (class == 'a')?"con":"cut",rowndx) ;
    pkrow->nme = nmebuf ; }
# ifdef PARANOIA
/*
  We need to postpone this check to here because pkvec_check requires that
  nme != NULL.
*/
  if (pkvec_check(pkrow,rtnnme) == FALSE) return (FALSE) ;
# endif
/*
  If necessary, boost the allocated capacity to open up space for the new
  row. (Space for the logical will also be allocated if coupling is enabled.)
  Then move the first cut out of the way if that's necessary to make a space
  for a new architectural constraints.
*/
  if (consys->rowsze < consys->concnt+1)
    if (consys_realloc(consys,'r',0) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,
	     "capacity expansion","row",pkrow->nme,pkrow->ndx) ;
      return (FALSE) ; }
  if (rowndx < consys->concnt+1)
    if (move_row(consys,rowndx,consys->concnt+1) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (FALSE) ; }
/*
  Initialise the row header and type.  Bump the appropriate constraint
  counts.
*/
  rowhdr = (rowhdr_struct *) CALLOC(1,sizeof(rowhdr_struct)) ;
  consys->mtx.rows[rowndx] = rowhdr ;
  rowhdr->ndx = rowndx ;
  rowhdr->nme = STRALLOC(pkrow->nme) ;
  consys->ctyp[rowndx] = contyp ;
  if (class == 'a')
    consys->archccnt++ ;
  else
    consys->cutccnt++ ;
  consys->concnt++ ;
/*
  If we're using our private buffer, hide it again.
*/
  if (pkrow->nme == nmebuf) pkrow->nme = rowhdr->nme ;
/*
  Add the row, keeping an eye out for any change in the maximum length
  column.  It's an error to try and insert a coefficient into a column that's
  not already in existence. If pkrow->sze == 0, we're just here to create the
  header.
*/
  if (pkrow->sze > 0)
  {
#   ifndef DYLP_NDEBUG
    if (pkrow->cnt == 0 && flgon(consys->opts,CONSYS_WRNZERO))
      warn(118,rtnnme,consys->nme,"row",rowhdr->nme,rowndx) ;
#   endif
    nzcnt = 0 ;
    pkcoeff = pkrow->coeffs ;
    for (vecndx = 0 ; vecndx < pkrow->cnt ; vecndx++)
    { if (pkcoeff->ndx <= 0 || pkcoeff->ndx > consys->varcnt)
      { errmsg(102,rtnnme,consys->nme,"column",pkcoeff->ndx,1,consys->varcnt) ;
	return (FALSE) ; }
      if (fabs(pkcoeff->val) >= consys->inf)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(128,rtnnme,consys->nme,rowndx,pkcoeff->ndx,pkcoeff->val,
	       "row",rowhdr->nme) ;
	return (FALSE) ; }
      if (fabs(pkcoeff->val) > consys->tiny)
      { colhdr = consys->mtx.cols[pkcoeff->ndx] ;
#       ifdef PARANOIA
	if (colhdr == NULL)
	{ errmsg(103,rtnnme,consys->nme,"column",pkcoeff->ndx) ;
	  return (FALSE) ; }
	if (pkcoeff->ndx != colhdr->ndx)
	{ errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,
		 pkcoeff->ndx,colhdr) ;
	  return (FALSE) ; }
#       endif
	coeff = (coeff_struct *) MALLOC(sizeof(coeff_struct)) ;
	coeff->colhdr = colhdr ;
	coeff->rowhdr = rowhdr ;
	coeff->val = pkcoeff->val ;
	coeff->rownxt = rowhdr->coeffs ;
	rowhdr->coeffs = coeff ;
	coeff->colnxt = colhdr->coeffs ;
	colhdr->coeffs = coeff ;
	colhdr->len++ ;
	nzcnt++ ;
	if (colhdr->len > consys->maxcollen)
	{ consys->maxcollen = colhdr->len ;
	  consys->maxcolndx = pkcoeff->ndx ; } }
#     ifndef DYLP_NDEBUG
      else
      { warn(130,rtnnme,consys->nme,rowndx,pkcoeff->ndx,pkcoeff->val,
	     consys->tiny,"row",rowhdr->nme) ; }
#     endif
      pkcoeff++ ; }
    rowhdr->len = nzcnt ;
    consys->mtx.coeffcnt += nzcnt ;
/*
  Check if this row should become the maximum length row.
*/
    if (rowhdr->len > consys->maxrowlen)
    { consys->maxrowlen = rowhdr->len ;
      consys->maxrowndx = rowndx ; } }
/*
  Check for the relevant associated arrays -- rhs, rhslow, clb, cub -- and
  set values if they're present.
*/
  if (consys->rhs != NULL) consys->rhs[rowndx] = rhs ;
  if (consys->rhslow != NULL) consys->rhslow[rowndx] = rhslow ;
  if (consys->clb != NULL && clb != NULL) consys->clb[rowndx] = *clb ;
  if (consys->cub != NULL && cub != NULL) consys->cub[rowndx] = *cub ;
/*
  Add a logical? If so, we'll move the first architectural out of the way.
  We might also have to move the logical associated with the first cut, if
  we're adding an architectural constraint when cuts are already present.
*/
  if (flgon(consys->opts,CONSYS_LVARS))
  { if (consys->archvcnt > 0)
      if (move_col(consys,consys->logvcnt+1,consys->varcnt+1) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->logvcnt+1,FALSE,NULL),
	       consys->logvcnt+1) ;
	return (FALSE) ; }
    if (rowndx < consys->concnt)
      if (move_col(consys,rowndx,consys->logvcnt+1) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
	return (FALSE) ; }
    if (add_logical(consys,rowndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(121,rtnnme,consys->nme,rowhdr->nme,rowndx) ;
      return (FALSE) ; } }
  
  return (TRUE) ; }



bool consys_getcol_pk (consys_struct *consys, int colndx, pkvec_struct **pkvec)

/*
  This routine copies a column out of the constraint matrix into the vector
  pkvec (creating a vector if needed). The routine will complain if a vector
  is supplied and it's too small for the column. As a special case, a vector
  of size 0 will retrieve only the column header information.

  Parameters:
    consys:	constraint system
    colndx:	column to be fetched
    pkvec:	(i) packed vector (if NULL, one will be created)
		(o) packed vector, loaded with column

  Returns: TRUE if the column is fetched without incident, FALSE otherwise
	   (paranoia or debug only).
*/

{ colhdr_struct *colhdr ;
  coeff_struct *coeff ;
  pkcoeff_struct *pkcoeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_getcol_pk" ;
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
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,
	   colhdr) ;
    return (FALSE) ; }
  if (pkvec == NULL)
  { errmsg(2,rtnnme,"&pkvec") ;
    return (FALSE) ; }
# endif
/*
  Get a fresh vector, if need be. Step along the column and copy the
  coefficients and row indices to the packed vector. Finally, fill in the
  header information for the packed vector.
*/
  if (*pkvec == NULL)
    *pkvec = pkvec_new(consys->maxcollen) ;
# ifdef PARANOIA
  else
    if (pkvec_check(*pkvec,rtnnme) == FALSE) return FALSE ;
# endif
  if ((*pkvec)->sze != 0)
  {
#   ifdef PARANOIA
    if ((*pkvec)->sze < colhdr->len)
    { errmsg(92,rtnnme,((*pkvec)->nme == NULL)?"<<null>>":(*pkvec)->nme,
	     (*pkvec)->ndx,(*pkvec)->sze,"column",
	     consys_nme(consys,'v',colndx,TRUE,NULL),colhdr->len) ;
      return (FALSE) ; }
#   endif
    for (coeff = colhdr->coeffs,pkcoeff = (*pkvec)->coeffs ;
	 coeff != NULL ;
	 coeff = coeff->colnxt,pkcoeff++)
    {
#     ifdef PARANOIA
      if (coeff->rowhdr == NULL)
      { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	       consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
	return (FALSE) ; }
      if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
      { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	       coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
	return (FALSE) ; }
#     endif
      pkcoeff->ndx = coeff->rowhdr->ndx ;
      pkcoeff->val = coeff->val ; } }
  (*pkvec)->ndx = colndx ;
  (*pkvec)->nme = colhdr->nme ;
  (*pkvec)->dim = consys->concnt ;
  (*pkvec)->dflt = 0 ;
  (*pkvec)->cnt = colhdr->len ;

  return (TRUE) ; }



bool consys_getcol_ex (consys_struct *consys, int colndx, double **vec)

/*
  This routine copies a column out of the constraint matrix into the vector
  vec (creating a vector if needed).

  Parameters:
    consys:	constraint system
    colndx:	column to be fetched
    vec:	(i) vector (if NULL, one will be created)
		(o) vector loaded with column

  Returns: TRUE if the column is fetched without incident, FALSE otherwise
	   (paranoia or debug only).
*/

{ colhdr_struct *colhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_getcol_ex" ;
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
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,
	   colhdr) ;
    return (FALSE) ; }
  if (vec == NULL)
  { errmsg(2,rtnnme,"&vec") ;
    return (FALSE) ; }
# endif
/*
  Make a vector, if need be, or clear the one we're handed. Then step along
  the column and set the necessary coefficients.
*/
  if (*vec == NULL)
    *vec = (double *) CALLOC(consys->concnt+1,sizeof(double)) ;
  else
    memset(*vec,0,(consys->concnt+1)*sizeof(double)) ;
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  {
#   ifdef PARANOIA
    if (coeff->rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (FALSE) ; }
    if (coeff->rowhdr != consys->mtx.rows[coeff->rowhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"row",coeff->rowhdr,coeff->rowhdr->ndx,
	     coeff->rowhdr->ndx,consys->mtx.rows[coeff->rowhdr->ndx]) ;
      return (FALSE) ; }
 #  endif   
    (*vec)[coeff->rowhdr->ndx] = coeff->val ; }

  return (TRUE) ; }




bool consys_getrow_pk (consys_struct *consys, int rowndx, pkvec_struct **pkvec)

/*
  This routine copies a row out of the constraint matrix into the vector
  pkvec (creating a vector if needed). The routine will complain if a vector
  is supplied and it's not big enough to hold the row. As a special case,
  a vector of size 0 will retrieve just the row header information.

  Parameters:
    consys:	constraint system
    rowndx:	row to be fetched
    pkvec:	(i) packed vector (if NULL, one will be created)
		(o) packed vector, loaded with row

  Returns: TRUE if the row is fetched without incident, FALSE otherwise
	   (paranoia or debug only).
*/

{ rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;
  pkcoeff_struct *pkcoeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_getrow_pk" ;
# endif

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
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,
	   rowhdr) ;
    return (FALSE) ; }
  if (pkvec == NULL)
  { errmsg(2,rtnnme,"&pkvec") ;
    return (FALSE) ; }
# endif
/*
  Get a fresh vector, if need be. Step along the row and copy the
  coefficients and column indices to the packed vector. Finally, fill in the
  header information for the packed vector.
*/
  if (*pkvec == NULL)
    *pkvec = pkvec_new(consys->maxrowlen) ;
# ifdef PARANOIA
  else
    if (pkvec_check(*pkvec,rtnnme) == FALSE) return FALSE ;
# endif
  if ((*pkvec)->sze != 0)
  {
#   ifdef PARANOIA
    if ((*pkvec)->sze < rowhdr->len)
    { errmsg(92,rtnnme,((*pkvec)->nme == NULL)?"<<null>>":(*pkvec)->nme,
	     (*pkvec)->ndx,(*pkvec)->sze,"row",
	     consys_nme(consys,'c',rowndx,TRUE,NULL),rowhdr->len) ;
      return (FALSE) ; }
#   endif
    for (coeff = rowhdr->coeffs,pkcoeff = (*pkvec)->coeffs ;
	 coeff != NULL ;
	 coeff = coeff->rownxt,pkcoeff++)
    {
#     ifdef PARANOIA
      if (coeff->colhdr == NULL)
      { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	       consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
	return (FALSE) ; }
      if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
      { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	       coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
	return (FALSE) ; }
 #    endif   
      pkcoeff->ndx = coeff->colhdr->ndx ;
      pkcoeff->val = coeff->val ; } }
  (*pkvec)->ndx = rowndx ;
  (*pkvec)->nme = rowhdr->nme ;
  (*pkvec)->dim = consys->varcnt ;
  (*pkvec)->dflt = 0 ;
  (*pkvec)->cnt = rowhdr->len ;

  return (TRUE) ; }



bool consys_getrow_ex (consys_struct *consys, int rowndx, double **vec)

/*
  This routine copies a row out of the constraint matrix into the vector
  vec (creating a vector if needed).

  Parameters:
    consys:	constraint system
    rowndx:	row to be fetched
    vec:	(i) vector (if NULL, one will be created)
		(o) vector loaded with row

  Returns: TRUE if the row is fetched without incident, FALSE otherwise
	   (paranoia or debug only).
*/

{ rowhdr_struct *rowhdr ;
  coeff_struct *coeff ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_getrow_ex" ;
# endif

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
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,
	   rowhdr) ;
    return (FALSE) ; }
  if (vec == NULL)
  { errmsg(2,rtnnme,"&vec") ;
    return (FALSE) ; }
# endif
/*
  Make a vector, if need be, or clear the one we're handed. Then step along
  the row and set the necessary coefficients.
*/
  if (*vec == NULL)
    *vec = (double *) CALLOC(consys->varcnt+1,sizeof(double)) ;
  else
    memset(*vec,0,(consys->varcnt+1)*sizeof(double)) ;
  for (coeff = rowhdr->coeffs ; coeff != NULL ; coeff = coeff->rownxt)
  {
#   ifdef PARANOIA
    if (coeff->colhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"colhdr",coeff,"row",
	     consys_nme(consys,'c',rowndx,FALSE,NULL),rowndx) ;
      return (FALSE) ; }
    if (coeff->colhdr != consys->mtx.cols[coeff->colhdr->ndx])
    { errmsg(126,rtnnme,consys->nme,"column",coeff->colhdr,coeff->colhdr->ndx,
	     coeff->colhdr->ndx,consys->mtx.cols[coeff->colhdr->ndx]) ;
      return (FALSE) ; }
 #  endif   
    (*vec)[coeff->colhdr->ndx] = coeff->val ; }

  return (TRUE) ; }



bool consys_delcol (consys_struct *consys, int colndx)

/*
  This routine removes a column from the constraint matrix. If needed, the
  last architectural variable is swapped into the vacated slot.

  Parameters:
    consys:	Constraint system
    colndx:	Index of the column to be deleted.

  Returns: TRUE if the deletion proceeded without error, FALSE otherwise.
*/

{ colhdr_struct *colhdr ;
  vartyp_enum vartyp ;
  bool rescan_rows,rescan_cols ;
  const char *rtnnme = "consys_delcol" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (consys->vtyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VTYP)) ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (colndx <= consys->logvcnt || colndx > consys->varcnt)
  { errmsg(102,rtnnme,consys->nme,"column",
	   colndx,consys->logvcnt+1,consys->varcnt) ;
    return (FALSE) ; }
# endif
  colhdr = consys->mtx.cols[colndx] ;
# ifdef PARANOIA
  if (colhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
    return (FALSE) ; }
  if (colndx != colhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,
	   colhdr) ;
    return (FALSE) ; }
# endif
/*
  Empty the column, forcing the column rescan flag if we're deleting the
  current maximum length column.
*/
  if (empty_col(consys,colndx,&rescan_rows) == FALSE)
  { setflg(consys->opts,CONSYS_CORRUPT) ;
    errmsg(112,rtnnme,consys->nme,"empty","column",colhdr->nme,colndx) ;
    return (FALSE) ; }
  if (colndx == consys->maxcolndx)
    rescan_cols = TRUE ;
  else
    rescan_cols = FALSE ;
/*
  Before we fill in the hole, correct the integer variable counts and check
  if we've removed the variable associated with the objective constraint.
  This is best done before columns are shifted.
*/
  if (colndx == consys->xzndx) consys->xzndx = -1 ;
  if (consys->vtyp != NULL)
  { vartyp = consys->vtyp[colndx] ;
    if (vartyp == vartypINT)
      consys->intvcnt-- ;
    else
    if (vartyp == vartypBIN)
      consys->binvcnt-- ; }
/*
  If this wasn't the last architectural variable, shift the last variable down
  to fill the hole, then correct the main variable counts. (Changing the main
  variable counts before moving the column will trigger a range error.)
*/
  if (colndx < consys->varcnt)
    if (move_col(consys,consys->varcnt,colndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","column",
	     consys_nme(consys,'v',consys->varcnt,FALSE,NULL),consys->varcnt) ;
      return (FALSE) ; }
  consys->archvcnt-- ;
  consys->varcnt-- ;
/*
  Rescan for maxima if necessary, and free the column header. That pretty well
  takes care of it.
*/
  if (rescan_rows == TRUE || rescan_cols == TRUE)
    if (find_maxes(consys,rescan_cols,rescan_rows) == FALSE)
    { errmsg(112,rtnnme,consys->nme,"maxima update","column",colhdr->nme,
	     colhdr->ndx) ;
      return (FALSE) ; }
  if (colhdr->nme != NULL) STRFREE(colhdr->nme) ;
  FREE(colhdr) ;

  return (TRUE) ; }



bool consys_delrow (consys_struct *consys, int rowndx)

/*
  This routine removes a row from the constraint matrix. If coupling is
  enabled, the corresponding logical variable is also removed. Any required
  row and column swaps are made to keep the constraint system compact.

  Parameters:
    consys:	Constraint system
    rowndx:	Index of the row to be deleted.

  Returns: TRUE if the deletion proceeded without error, FALSE otherwise.
*/

{ int colndx ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  bool rescan_rows,rescan_cols ;
  const char *rtnnme = "consys_delrow" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (flgon(consys->opts,CONSYS_LVARS))
    if (consys->logvcnt != consys->concnt)
    { errmsg(131,rtnnme,consys->nme,consys->logvcnt,consys->concnt) ;
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
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,
	   rowhdr) ;
    return (FALSE) ; }
# endif
/*
  If coupling is on, delete the logical variable (this avoids a pointless
  warning about a 0-length column in the usual case where the logical
  variable appears only in its associated constraint). We'll need to swap the
  last logical into the space vacated by the deleted logical, and swap the
  last architectural variable into the space vacated by the last logical. We
  may need to do an additional swap within the logicals, to maintain the
  separation between architectural and cut constraints.
*/
  rescan_rows = FALSE ;
  if (flgon(consys->opts,CONSYS_LVARS))
  { colndx = rowndx ;
    colhdr = consys->mtx.cols[colndx] ;
#   ifdef PARANOIA
    if (colhdr == NULL)
    { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
      return (FALSE) ; }
    if (colhdr->ndx != colndx)
    { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
      return (FALSE) ; }
#   endif
    if (empty_col(consys,colndx,&rescan_rows) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"empty","column",colhdr->nme,colndx) ;
      return (FALSE) ; }
    if (colhdr->nme != NULL) STRFREE(colhdr->nme) ;
    FREE(colhdr) ;
    if (colndx < consys->archccnt && consys->cutccnt > 0)
    { if (move_col(consys,consys->archccnt,colndx) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->archccnt,FALSE,NULL),
	       consys->archccnt) ;
	return (FALSE) ; }
      colndx = consys->archccnt ; }
    if (colndx < consys->logvcnt)
      if (move_col(consys,consys->logvcnt,colndx) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->logvcnt,FALSE,NULL),
	       consys->logvcnt) ;
	return (FALSE) ; }
    if (consys->archvcnt > 0)
      if (move_col(consys,consys->varcnt,consys->logvcnt) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->varcnt,FALSE,NULL),
	       consys->varcnt) ;
	return (FALSE) ; }
    consys->logvcnt-- ;
    consys->varcnt-- ; }
/*
  Now empty the row, forcing the row rescan flag if we're deleting the
  current maximum length row, and invalidating objndx if we've deleted the
  objective constraint. We'll need to do the same reordering of the
  constraints that we applied to the logicals.
*/
  if (empty_row(consys,rowndx,&rescan_cols) == FALSE)
  { setflg(consys->opts,CONSYS_CORRUPT) ;
    errmsg(112,rtnnme,consys->nme,"empty","row",rowhdr->nme,rowndx) ;
    return (FALSE) ; }
  if (rowndx == consys->maxrowndx) rescan_rows = TRUE ;
  if (rowndx == consys->objndx) consys->objndx = -1 ;
  if (rowndx < consys->archccnt && consys->cutccnt > 0)
  { if (move_row(consys,consys->archccnt,rowndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","row",
	     consys_nme(consys,'c',consys->archccnt,FALSE,NULL),
	     consys->archccnt) ;
      return (FALSE) ; }
    rowndx = consys->archccnt ; }
  if (rowndx < consys->concnt)
    if (move_row(consys,consys->concnt,rowndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","row",
	     consys_nme(consys,'c',consys->concnt,FALSE,NULL),consys->concnt) ;
      return (FALSE) ; }
  if (rowhdr->ndx <= consys->archccnt)
    consys->archccnt-- ;
  else
    consys->cutccnt-- ;
  consys->concnt-- ;
/*
  Rescan for maxima if necessary, and free the row header. That pretty well
  takes care of it.
*/
  if (rescan_rows == TRUE || rescan_cols == TRUE)
    if (find_maxes(consys,rescan_cols,rescan_rows) == FALSE)
    { errmsg(112,rtnnme,consys->nme,"maxima update","row",rowhdr->nme,
	     rowhdr->ndx) ;
      return (FALSE) ; }
  if (rowhdr->nme != NULL) STRFREE(rowhdr->nme) ;
  FREE(rowhdr) ;

  return (TRUE) ; }



bool consys_delrow_stable (consys_struct *consys, int rowndx)

/*
  [Jun 01, 01] Same as consys_delrow, except it maintains the order of
  constraints (thus the name 'stable').... TODO

  This routine removes a row from the constraint matrix. If coupling is
  enabled, the corresponding logical variable is also removed. Any required
  row and column swaps are made to keep the constraint system compact.

  Parameters:
    consys:	Constraint system
    rowndx:	Index of the row to be deleted.

  Returns: TRUE if the deletion proceeded without error, FALSE otherwise.
*/

{ int colndx,i ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;
  bool rescan_rows,rescan_cols ;
  const char *rtnnme = "consys_delrow" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (flgon(consys->opts,CONSYS_LVARS))
    if (consys->logvcnt != consys->concnt)
    { errmsg(131,rtnnme,consys->nme,consys->logvcnt,consys->concnt) ;
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
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,
	   rowhdr) ;
    return (FALSE) ; }
# endif
/*
  If coupling is on, delete the logical variable (this avoids a pointless
  warning about a 0-length column in the usual case where the logical
  variable appears only in its associated constraint). We'll need to swap the
  last logical into the space vacated by the deleted logical, and swap the
  last architectural variable into the space vacated by the last logical. We
  may need to do an additional swap within the logicals, to maintain the
  separation between architectural and cut constraints.
*/
  rescan_rows = FALSE ;
  if (flgon(consys->opts,CONSYS_LVARS))
  { colndx = rowndx ;
    colhdr = consys->mtx.cols[colndx] ;
#   ifdef PARANOIA
    if (colhdr == NULL)
    { errmsg(103,rtnnme,consys->nme,"column",colndx) ;
      return (FALSE) ; }
    if (colhdr->ndx != colndx)
    { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,colhdr) ;
      return (FALSE) ; }
#   endif
    if (empty_col(consys,colndx,&rescan_rows) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"empty","column",colhdr->nme,colndx) ;
      return (FALSE) ; }
    if (colhdr->nme != NULL) STRFREE(colhdr->nme) ;
    FREE(colhdr) ;
    if (colndx < consys->archccnt && consys->cutccnt > 0)
    { if (move_col(consys,consys->archccnt,colndx) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->archccnt,FALSE,NULL),
	       consys->archccnt) ;
	return (FALSE) ; }
      colndx = consys->archccnt ; }
    if (colndx < consys->logvcnt)
      if (move_col(consys,consys->logvcnt,colndx) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->logvcnt,FALSE,NULL),
	       consys->logvcnt) ;
	return (FALSE) ; }
    if (consys->archvcnt > 0)
      if (move_col(consys,consys->varcnt,consys->logvcnt) == FALSE)
      { setflg(consys->opts,CONSYS_CORRUPT) ;
        errmsg(112,rtnnme,consys->nme,"swap","column",
	       consys_nme(consys,'v',consys->varcnt,FALSE,NULL),
	       consys->varcnt) ;
	return (FALSE) ; }
    consys->logvcnt-- ;
    consys->varcnt-- ; }
/*
  Now empty the row, forcing the row rescan flag if we're deleting the
  current maximum length row, and invalidating objndx if we've deleted the
  objective constraint. We'll need to do the same reordering of the
  constraints that we applied to the logicals.
*/
  if (empty_row(consys,rowndx,&rescan_cols) == FALSE)
  { setflg(consys->opts,CONSYS_CORRUPT) ;
    errmsg(112,rtnnme,consys->nme,"empty","row",rowhdr->nme,rowndx) ;
    return (FALSE) ; }
  if (rowndx == consys->maxrowndx) rescan_rows = TRUE ;
  if (rowndx == consys->objndx) consys->objndx = -1 ;

  if (rowndx < consys->archccnt && consys->cutccnt > 0)
  { if (move_row(consys,consys->archccnt,rowndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","row",
	     consys_nme(consys,'c',consys->archccnt,FALSE,NULL),
	     consys->archccnt) ;
      return (FALSE) ; }
    rowndx = consys->archccnt ; }
  if (rowndx < consys->concnt)
    if (move_row(consys,consys->concnt,rowndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","row",
	     consys_nme(consys,'c',consys->concnt,FALSE,NULL),consys->concnt) ;
      return (FALSE) ; }
  for (i=rowndx; i<consys->concnt; i++)
  { if (move_row(consys,i+1,i) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(112,rtnnme,consys->nme,"swap","row",
	     consys_nme(consys,'c',i+1,FALSE,NULL),i+1) ;
      return (FALSE) ; } }
  if (rowhdr->ndx <= consys->archccnt)
    consys->archccnt-- ;
  else
    consys->cutccnt-- ;
  consys->concnt-- ;
/*
  Rescan for maxima if necessary, and free the row header. That pretty well
  takes care of it.
*/
  if (rescan_rows == TRUE || rescan_cols == TRUE)
    if (find_maxes(consys,rescan_cols,rescan_rows) == FALSE)
    { errmsg(112,rtnnme,consys->nme,"maxima update","row",rowhdr->nme,
	     rowhdr->ndx) ;
      return (FALSE) ; }
  if (rowhdr->nme != NULL) STRFREE(rowhdr->nme) ;
  FREE(rowhdr) ;

  return (TRUE) ; }



double consys_getcoeff (consys_struct *consys, int rowndx, int colndx)

/*
  This routine retrieves a single coefficient a<rowndx,colndx> from the
  constraint matrix. The specified row and column must exist. If the
  coefficient doesn't exist, the routine assumes it's zero.

  consys_getcoeff searches for the coefficient by scanning the column, on the
  observation that in lp constraint systems the number of variables is
  usually much greater than the number of constraints.

  Parameters:
    consys:	the constraint system
    rowndx:	row index
    colndx:	column index

  Returns: a<rowndx,colndx> if the coefficient exists,
	   0 if rowndx and colndx are valid, but a<rowndx,colndx> is absent,
	   NaN otherwise (paranoia only)
*/

{ int lclndx ;
  coeff_struct *coeff ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "consys_getcoeff" ;
# endif

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (quiet_nan(0)) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (quiet_nan(0)) ; }
  if (flgon(consys->opts,CONSYS_LVARS))
    if (consys->logvcnt != consys->concnt)
    { errmsg(131,rtnnme,consys->nme,consys->logvcnt,consys->concnt) ;
      return (quiet_nan(0)) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (quiet_nan(0)) ; }
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
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,
	   colhdr) ;
    return (quiet_nan(0)) ; }
  rowhdr = consys->mtx.rows[rowndx] ;
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (quiet_nan(0)) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,
	   rowhdr) ;
    return (quiet_nan(0)) ; }
# endif
/*
  Scan the column for the requested coefficient.
*/
  for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
  { rowhdr = coeff->rowhdr ;
#   ifdef PARANOIA
    if (rowhdr == NULL)
    { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	     consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
      return (quiet_nan(0)) ; }
#   endif
    lclndx = rowhdr->ndx ;
#   ifdef PARANOIA
    if (lclndx <= 0 || lclndx > consys->concnt)
    { errmsg(102,rtnnme,consys->nme,"row",lclndx,1,consys->concnt) ;
      return (quiet_nan(0)) ; }
    if (consys->mtx.rows[lclndx] != rowhdr)
    { errmsg(126,rtnnme,consys->nme,"row",rowhdr,lclndx,lclndx,
	     consys->mtx.rows[lclndx]) ;
      return (quiet_nan(0)) ; }
#   endif
    if (lclndx == rowndx) break ; }
/*
  Return the value, or zero.
*/
  if (coeff == NULL)
    return (0.0) ;
  else
    return (coeff->val) ; }



bool consys_setcoeff (consys_struct *consys,
		      int rowndx, int colndx, double val)

/*
  This routine modifies a single coefficient a<rowndx,colndx> in the constraint
  matrix. If rowndx and colndx are valid but the coefficient doesn't exist,
  it will be created. If the coefficient exists, but val is zero, it will be
  removed.

  consys_setcoeff searches for the coefficient by scanning the column, on the
  observation that in lp constraint systems the number of variables is
  usually much greater than the number of constraints.

  Parameters:
    consys:	the constraint system
    rowndx:	row index
    colndx:	column index
    val:	new value of the coefficient

  Returns: TRUE if the indices are valid and the coefficient is properly set,
	   FALSE otherwise.
*/

{ int lclndx ;
  bool scanrows,scancols ;
  coeff_struct *coeff,**pcoeff ;
  colhdr_struct *colhdr ;
  rowhdr_struct *rowhdr ;

  const char *rtnnme = "consys_chgcoeff" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (flgon(consys->opts,CONSYS_LVARS))
    if (consys->logvcnt != consys->concnt)
    { errmsg(131,rtnnme,consys->nme,consys->logvcnt,consys->concnt) ;
      return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (FALSE) ; }
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
  { errmsg(126,rtnnme,consys->nme,"column",colhdr,colhdr->ndx,colndx,
	   colhdr) ;
    return (FALSE) ; }
  rowhdr = consys->mtx.rows[rowndx] ;
  if (rowhdr == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (FALSE) ; }
  if (rowndx != rowhdr->ndx)
  { errmsg(126,rtnnme,consys->nme,"row",rowhdr,rowhdr->ndx,rowndx,
	   rowhdr) ;
    return (FALSE) ; }
# endif
  if (fabs(val) >= consys->inf)
  { setflg(consys->opts,CONSYS_CORRUPT) ;
    errmsg(128,rtnnme,consys->nme,rowndx,colndx,val,
	   "coefficient","<no name>") ;
    return (FALSE) ; }
/*
  Scan the column for the requested coefficient. If we find it, we can change
  the value and we're finished.
*/
  if (val != 0)
  { for (coeff = colhdr->coeffs ; coeff != NULL ; coeff = coeff->colnxt)
    { rowhdr = coeff->rowhdr ;
#     ifdef PARANOIA
      if (rowhdr == NULL)
      { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	       consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
	return (FALSE) ; }
#     endif
      lclndx = rowhdr->ndx ;
#     ifdef PARANOIA
      if (lclndx <= 0 || lclndx > consys->concnt)
      { errmsg(102,rtnnme,consys->nme,"row",lclndx,1,consys->concnt) ;
	return (FALSE) ; }
      if (consys->mtx.rows[lclndx] != rowhdr)
      { errmsg(126,rtnnme,consys->nme,"row",rowhdr,lclndx,lclndx,
	       consys->mtx.rows[lclndx]) ;
	return (FALSE) ; }
#     endif
      if (lclndx == rowndx) break ; }
/*
  Found it. Set the value and we're done.
*/
    if (coeff != NULL)
    { coeff->val = val ;
      return (TRUE) ; }
/*
  Sigh. No coefficient, we'll have to create one.
*/
    rowhdr = consys->mtx.rows[rowndx] ;
    coeff = (coeff_struct *) MALLOC(sizeof(coeff_struct)) ;
    coeff->rowhdr = rowhdr ;
    coeff->colhdr = colhdr ;
    coeff->val = val ;
    coeff->rownxt = rowhdr->coeffs ;
    rowhdr->coeffs = coeff ;
    rowhdr->len++ ;
    if (rowhdr->len > consys->maxrowlen)
    { consys->maxrowlen = rowhdr->len ;
      consys->maxrowndx = rowndx ; }
    colhdr->len++ ;
    if (colhdr->len > consys->maxcollen)
    { consys->maxcollen = colhdr->len ;
      consys->maxcolndx = colndx ; } }
/*
  If val is zero, we have to run a slightly more complicated search loop which
  will allow us to delete the coefficient.
*/
  else
  { for (pcoeff = &colhdr->coeffs, coeff = *pcoeff ;
	 coeff != NULL ;
	 pcoeff = &coeff->colnxt, coeff = *pcoeff)
    { rowhdr = coeff->rowhdr ;
#     ifdef PARANOIA
      if (rowhdr == NULL)
      { errmsg(125,rtnnme,consys->nme,"rowhdr",coeff,"column",
	       consys_nme(consys,'v',colndx,FALSE,NULL),colndx) ;
	return (FALSE) ; }
#     endif
      lclndx = rowhdr->ndx ;
#     ifdef PARANOIA
      if (lclndx <= 0 || lclndx > consys->concnt)
      { errmsg(102,rtnnme,consys->nme,"row",lclndx,1,consys->concnt) ;
	return (FALSE) ; }
      if (consys->mtx.rows[lclndx] != rowhdr)
      { errmsg(126,rtnnme,consys->nme,"row",rowhdr,lclndx,lclndx,
	       consys->mtx.rows[lclndx]) ;
	return (FALSE) ; }
#     endif
      if (lclndx == rowndx) break ; }
/*
  Best if we don't find the coefficient, actually.
*/
    if (coeff == NULL) return (TRUE) ;
/*
  But if we do, delink it from the column, then find it row-wise and delink it
  from the row, then free the space.
*/
    *pcoeff = coeff->colnxt ;
    for (pcoeff = &rowhdr->coeffs ;
	 *pcoeff != NULL ;
	 pcoeff = &(*pcoeff)->rownxt)
    { if (*pcoeff == coeff) break ; }
#   ifdef PARANOIA
    if (*pcoeff == NULL)
    { errmsg(119,rtnnme,consys->nme,rowndx,colndx,coeff->val,
	     "column",colhdr->ndx,"row",rowhdr->ndx) ;
      return (FALSE) ; }
#   endif
    *pcoeff = coeff->rownxt ;
    FREE(coeff) ;
/*
  Correct the column and row lengths, and the total coefficient count, and
  rescan for maximum column and row if necessary.
*/
    consys->mtx.coeffcnt-- ;
    colhdr->len-- ;
    if (colndx == consys->maxcolndx)
      scancols = TRUE ;
    else
      scancols = FALSE ;
    rowhdr->len-- ;
    if (rowndx == consys->maxrowndx)
      scanrows = TRUE ;
    else
      scanrows = FALSE ;
    if (scancols == TRUE || scanrows == TRUE)
    { if (find_maxes(consys,scancols,scanrows) == FALSE)
      { errmsg(112,rtnnme,consys->nme,"maxima update","column",colhdr->nme,
	       colhdr->ndx) ;
	return (FALSE) ; } } }
    
  return (TRUE) ; }



bool consys_logicals (consys_struct *consys)

/*
  This routine establishes logical variables for a constraint system.

  Parameters:
    consys:	the constraint system

  Returns: TRUE if logicals are set up successfully, FALSE otherwise.
*/

{ int ndx ;
  const char *rtnnme = "consys_logicals" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,"row header") ;
    return (FALSE) ; }
  if (consys->mtx.cols == NULL)
  { errmsg(101,rtnnme,consys->nme,"column header") ;
    return (FALSE) ; }
  if (flgon(consys->opts,CONSYS_LVARS) || consys->logvcnt > 0)
  { errmsg(123,rtnnme,consys->nme) ;
    return (FALSE) ; }
# endif

/*
  Make sure the constraint system has enough room for the logicals.
*/
  ndx = consys->archvcnt+consys->concnt-consys->colsze ;
  if (ndx > 0)
    if (consys_realloc(consys,'c',ndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(124,rtnnme,consys->nme) ;
      return (FALSE) ; }
/*
  Now step through the rows and create logical variables.
*/
  for (ndx = 1 ; ndx <= consys->concnt ; ndx++)
    if (add_logical(consys,ndx) == FALSE)
    { setflg(consys->opts,CONSYS_CORRUPT) ;
      errmsg(121,rtnnme,consys->nme,
	     consys_nme(consys,'c',ndx,FALSE,NULL),ndx) ;
      return (FALSE) ; }
/*
  Turn on coupling and we're done.
*/
  setflg(consys->opts,CONSYS_LVARS) ;

  return (TRUE) ; }


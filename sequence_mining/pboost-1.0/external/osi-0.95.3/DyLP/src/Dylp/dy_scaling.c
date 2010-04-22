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
  This file contains routines which handle scaling of the original constraint
  system and unscaling when producing output. Intimately tied to this is the
  question of whether dylp will make a local copy of the original constraint
  system or simply refer to the system supplied by the client.

  When dylp is first called, it will call dy_initsystem to evaluate the
  situation and establish the `original system' seen by the rest of dylp.
  The rules are this:
    * If the client forces a local copy (copyorigsys == TRUE), we make a local
      copy of the original system.
    * If the client supplied scaling matrices, we make a scaled local copy of
      the original system.
    * If scaling is allowed, and an evaluation of the constraint system
      establishes that it's necessary, we make a scaled local copy of the
      original system.
  Otherwise, no local copy is required, and dylp can reference the original
  system supplied by the client.
*/


#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_scaling.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_scaling.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  Define this symbol to enable checks on the calculation of rows of the
  unscaled basis inverse.

  #define CHECK_UNSCALED_BETAI
*/

/*
  Constraint system pointers maintained by this module.

  local_sys: Pointer to the local copy of the original constraint system.
	     Initialized by dy_initlclsystem, destroyed by dy_freelclsystem.
	     NULL if there is no local copy.

  client_sys: Pointer to the system originally provided by the client. Used
	      to remember the client's system while we're using a local copy.

*/

static consys_struct *local_sys, *client_sys ;

/*
  If we own scaling vectors, they are stored in these variables.
*/

static double *lcl_rowscale,*lcl_colscale ;

bool dy_isscaled (void)
/*
  A trivial routine to allow other parts of dylp to learn if the system is
  scaled.
*/

{ return (local_sys != NULL) ; }



bool dy_initlclsystem (lpprob_struct *orig_lp, bool hotstart)

/*
  This routine looks at the constraint system and options provided by the
  client and decides whether dylp can reference the supplied constraint
  system directly or whether a local copy should be made. A local copy is
  required if the constraint system is to be scaled. Options specified by the
  client can force the creation of a local copy, and force or forbid the use
  of scaling. If scaling is allowed but not forced, this routine will evaluate
  the constraint matrix and apply scaling if necessary.

  If we're doing a hot start, this call is strictly for information hiding.
  All the data structures should exist, and it's just a question of swapping
  orig_lp->consys, should we need to do it. This form of the call is also used
  when we're freeing data structures from previous runs.

  Parameters:
    orig_lp:	(i) the original lp problem, as supplied by the client
		(o) orig_lp->consys may be replaced with a local copy of the
		    constraint system
    hotstart:	TRUE if we're doing a hot start, FALSE otherwise

  Returns: TRUE if all goes well, FALSE if there is a failure.
*/

{ int scalefactor ;
  bool localcopy,scale ;
  flags scaled_vecs ;
  double orig_scm,scaled_scm ;

  consys_struct *orig_sys ;

  const char *rtnnme = "dy_initlclsystem" ;

# ifdef PARANOIA
  if (orig_lp == NULL)
  { errmsg(2,rtnnme,"orig_lp") ;
    return (FALSE) ; }
  if (orig_lp->consys == NULL)
  { errmsg(2,rtnnme,"orig_lp->consys") ;
    return (FALSE) ; }
# endif

/*
  If this is a hot start, all we need to do here is (possibly) set orig_lp to
  point to the existing local copy. More is likely needed, but that'll be
  handled via dy_refreshlclsystem, called from dy_hotstart.
*/
  if (hotstart == TRUE)
  { if (local_sys != NULL)
    { client_sys = orig_lp->consys ;
      orig_lp->consys = local_sys ; }
    return (TRUE) ; }
/*
  It's not a hot start. Calculate the geometric mean for the original system.
*/
  orig_sys = orig_lp->consys ;
  orig_scm = consys_evalsys(orig_sys) ;
/*
  Decide if we need a local copy, and if we need to scale the constraint
  matrix. To decide if we need a scaled copy, we need to scan the constraint
  matrix.
*/
  local_sys = NULL ;
  client_sys = NULL ;
  lcl_rowscale = NULL ;
  lcl_colscale = NULL ;

  localcopy = dy_opts->copyorigsys ;

  switch (dy_opts->scaling)
  { case 0: /* scaling prohibited */
    { scale = FALSE ;
      break ; }
    case 1: /* scale with user-supplied matrices */
    { localcopy = TRUE ;
      scale = TRUE ;
#     ifdef PARANOIA
      if (orig_sys->rowscale == NULL)
      { errmsg(101,rtnnme,"row scaling vector") ;
	return (FALSE) ; }
      if (orig_sys->colscale == NULL)
      { errmsg(101,rtnnme,"column scaling vector") ;
	return (FALSE) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.scaling >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      [%s]: scaling with client vectors",
		    orig_sys->nme) ; }
#     endif
      break ; }
    case 2: /* scale if necessary */
    { if (orig_sys->minaij >= .5 && orig_sys->maxaij <= 2)
      { scale = FALSE ; }
      else
      { localcopy = TRUE ;
	scale = TRUE ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.scaling >= 2)
      { if (scale == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n      [%s]: system will be scaled;",
		      orig_sys->nme) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n      [%s]: scaling not required;",
		      orig_sys->nme) ; }
	dyio_outfmt(dy_logchn,dy_gtxecho," %g <= |a<ij>| <= %g, metric = %g.",
		    orig_sys->minaij,orig_sys->maxaij,orig_scm) ; }
#     endif
      break ; }
    default:
    { errmsg(7,rtnnme,__LINE__,"scaling option code",dy_opts->scaling) ;
      return (FALSE) ; } }
# ifdef PARANOIA
  if (scale == TRUE && localcopy == FALSE)
  { errmsg(1,rtnnme,__LINE__) ;
    return (FALSE) ; }
# endif
/*
  If we need a local copy, make one.  We don't need to duplicate all the
  attached vectors of the original system, only the ones relevant for LP.
*/
  if (localcopy == TRUE)
  { scaled_vecs = CONSYS_OBJ|
		  CONSYS_CTYP|CONSYS_RHS|CONSYS_RHSLOW|CONSYS_RSCALE|
		  CONSYS_VTYP|CONSYS_VLB|CONSYS_VUB|CONSYS_CSCALE ;
    if (consys_dupsys(orig_sys,&local_sys,scaled_vecs) == FALSE)
    { errmsg(137,rtnnme,orig_sys->nme) ;
      return (FALSE) ; }
    client_sys = orig_lp->consys ;
    orig_lp->consys = local_sys ; }

/*
  Do we need to scale? If the client provided us with scaling vectors, store
  them with the local copy. A full attach isn't required, as the local copy
  will be static through the dylp run and any hot starts.
*/
  if (scale == TRUE)
  { if (dy_opts->scaling == 1)
    { local_sys->rowscale = client_sys->rowscale ;
      local_sys->colscale = client_sys->colscale ; }
/*
  If no scaling vectors were supplied, call consys_geomscale and
  consys_equiscale to calculate the scaling vectors.
*/
    else
    { if (consys_geomscale(local_sys,&local_sys->rowscale,
					  &local_sys->colscale) == FALSE)
      { errmsg(135,rtnnme,local_sys->nme) ;
	return (FALSE) ; }
      lcl_rowscale = local_sys->rowscale ;
      lcl_colscale = local_sys->colscale ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.scaling >= 2)
      { scaled_scm = sqrt(local_sys->maxaij/local_sys->minaij) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      [%s]: after geometric scaling",
		    local_sys->nme) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," %g <= |a<ij>| <= %g, geom = %g.",
		    local_sys->minaij,local_sys->maxaij,scaled_scm) ; }
#     endif
      if (consys_equiscale(local_sys,&local_sys->rowscale,
					  &local_sys->colscale) == FALSE)
      { errmsg(135,rtnnme,local_sys->nme) ;
	return (FALSE) ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.scaling >= 2)
      { scaled_scm = sqrt(local_sys->maxaij/local_sys->minaij) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      [%s]: after equilibration scaling",
		    local_sys->nme) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," %g <= |a<ij>| <= %g, geom = %g.",
		    local_sys->minaij,local_sys->maxaij,scaled_scm) ; }
#       endif
    }
/*
  Apply the scaling vectors and report the result.
*/
    if (consys_applyscale(local_sys,local_sys->rowscale,
					  local_sys->colscale) == FALSE)
    { errmsg(135,rtnnme,local_sys->nme) ;
      return (FALSE) ; }
    scaled_scm = sqrt(local_sys->maxaij/local_sys->minaij) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.scaling >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      [%s]: after scaling",
		    local_sys->nme) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," %g <= |a<ij>| <= %g, geom = %g.",
		    local_sys->minaij,local_sys->maxaij,scaled_scm) ; }
#   endif
  }
  else
  { scaled_scm = orig_scm ; }
/*
  How'd we do? We may still want to adjust the zero tolerance and feasibility
  scaling factors. Use the ratio of the original geometric mean to the scaled
  geometric mean.
*/
  scalefactor = ((int)(log10(orig_scm/scaled_scm)+.5)-1) ;
  if (scalefactor > 1)
  { dy_tols->pfeas_scale *= pow(10.0,(double) scalefactor) ;
    dy_tols->dfeas_scale *= pow(10.0,(double) scalefactor) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.scaling >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    [%s]: adjusting pfeas scale by 1.0e+%d to %g.",
		  local_sys->nme,scalefactor,dy_tols->pfeas_scale) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    [%s]: adjusting dfeas scale by 1.0e+%d to %g.",
		  local_sys->nme,scalefactor,dy_tols->dfeas_scale) ; } 
#   endif
  }
  if (scalefactor > 2)
  { scalefactor -= 2 ;
    dy_tols->zero /= pow(10.0,(double) scalefactor) ;
    dy_tols->cost /= pow(10.0,(double) scalefactor) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.scaling >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    [%s]: adjusting primal zero by 1.0e-%d to %g.",
		  local_sys->nme,scalefactor,dy_tols->zero) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    [%s]: adjusting dual zero by 1.0e-%d to %g.",
		  local_sys->nme,scalefactor,dy_tols->cost) ; }
#   endif
  }

  return (TRUE) ; }



void dy_refreshlclsystem (flags what)

/*
  This routine is called by dy_hotstart to transfer the client's changes in
  rhs, bounds, or objective from client_sys to local_sys. It's moderately
  brutal. For each vector that's changed, we do a wholesale copy.

  Parameters:
    what:	the vectors that are to be refreshed

  Returns: undefined
*/

{ int i,j,m,n ;
  double infinity ;

# ifdef PARANOIA
  const char *rtnnme = "dy_refreshlclsystem" ;
# endif

/*
  If there's no scaling, nothing needs to be done.
*/
  if (local_sys == NULL) return ;

# ifdef PARANOIA
  if (getflg(client_sys->opts,CONSYS_FININF) !=
      getflg(local_sys->opts,CONSYS_FININF))
  { errmsg(1,rtnnme,__LINE__) ;
    return ; }
# endif

/*
  For each of the vectors that can change, do a refresh if requested. When
  we're doing the bounds, watch out for finite infinity.
*/
  m = client_sys->concnt ;
  n = client_sys->varcnt ;
  infinity = client_sys->inf ;

  if (flgon(what,lpctlOBJCHG))
  { for (j = 1 ; j <= n ; j++)
      local_sys->obj[j] = client_sys->obj[j]*lcl_colscale[j] ; }

  if (flgon(what,lpctlRHSCHG))
  { for (i = 1 ; i <= m ; i++)
    { local_sys->rhs[i] = client_sys->rhs[i]*lcl_rowscale[i] ;
      local_sys->rhslow[i] = client_sys->rhslow[i]*lcl_rowscale[i] ; } }

  if (flgon(what,lpctlLBNDCHG))
  { if (flgon(client_sys->opts,CONSYS_FININF))
    { for (j = 1 ; j <= n ; j++)
      { if (client_sys->vlb[j] > -infinity)
	{ local_sys->vlb[j] = client_sys->vlb[j]/lcl_colscale[j] ; }
	else
	{ local_sys->vlb[j] = -infinity ; } } }
    else
    { for (j = 1 ; j <= n ; j++)
	local_sys->vlb[j] = client_sys->vlb[j]/lcl_colscale[j] ; } }

  if (flgon(what,lpctlUBNDCHG))
  { if (flgon(client_sys->opts,CONSYS_FININF))
    { for (j = 1 ; j <= n ; j++)
      { if (client_sys->vub[j] < infinity)
	{ local_sys->vub[j] = client_sys->vub[j]/lcl_colscale[j] ; }
	else
	{ local_sys->vub[j] = infinity ; } } }
    else
    { for (j = 1 ; j <= n ; j++)
	local_sys->vub[j] = client_sys->vub[j]/lcl_colscale[j] ; } }

  return ; }



void dy_freelclsystem (lpprob_struct *orig_lp, bool freesys)

/*
  This routine cleans up the local copy of the constraint system. If there's no
  local copy, the routine is a noop. The minimal action for a local copy is to
  correct the pointer in orig_lp. If free is TRUE, then the local copy is
  released.

  Parameters:
    orig_lp:	(i) the original lp problem, as supplied by the client
		(o) orig_lp->consys will be restored to the client copy of the
		    constraint system
    freesys:	TRUE to free the local system

  Returns: undefined
*/

{

/*
  Do we even have a local copy? If not, we're done already.
*/
  if (local_sys == NULL) return ;
/*
  Replace the consys pointer in orig_lp with a pointer to the client's
  constraint system.
*/
  orig_lp->consys = client_sys ;
  if (freesys == FALSE) return ;
/*
  We have work to do. Clear the row and column scaling vectors, if necessary,
  then free the constraint system. It's important to set these pointers back to
  NULL, so that we don't try redundant free's in an environment (like COIN)
  where multiple clients may try to free data structures. Likewise, the scaling
  vectors are associated with some constraint system, so don't free them here.
*/
  if (lcl_rowscale != NULL)
  { lcl_rowscale = NULL ; }
  if (lcl_colscale != NULL)
  { lcl_colscale = NULL ; }
  if (local_sys != NULL)
  { consys_free(local_sys) ;
    local_sys = NULL ; }

  return  ; }



void dy_unscale_cbar (int nbcnt, double *cbar, int *vndx)

/*
  This is a special purpose routine which unscales the vector of selected
  reduced costs produced by dy_pricenbvars. See the head of the file for a
  description of the math. All we do here is walk the vectors and apply the
  unscaling.

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

/*
  Is unscaling required?
*/
  if (lcl_colscale == NULL) return ;
/*
  Get on with the calculation. Architectural columns j must be scaled by
  inv(S)<j>. Columns for logical variables must be scaled by R<i>, where i is
  the index of the associated constraint. Recall that vndx encodes the index
  of a logical as -i.
*/
  for (k = 0 ; k < nbcnt ; k++)
  { j = vndx[k] ;
    cbarj = cbar[k] ;
    if (j > 0)
    { cbarj /= lcl_colscale[j] ; }
    else
    { cbarj *= lcl_rowscale[-j] ; }
    setcleanzero(cbarj,dy_tols->dfeas) ;
    cbar[k] = cbarj ; }

  return ; }



#ifdef CHECK_UNSCALED_BETAI

static bool check_unscalebetai (consys_struct *orig_sys,
				int ipos, double *betai)

/*
  This routine performs two checks on betai. First, it calculates
  dot(beta<i>,a<j>) for each basic column a<j>, to make sure that we have
  e<i>. Then it calculates the entries of beta<i> by extracting the relevant
  entry of the column vector beta<k> = inv(B)e<k>, k = 1 .. dy_sys->concnt.
  Valid only for active rows of dy_sys.

  The tolerances on the checks are pretty wide --- 1000*dy_tols.zero --- but
  then we wouldn't be scaling if the unscaled system was well-behaved.

  Parameters:
    orig_sys:	the original (unscaled) constraint system.
    ipos:	index (position) of basis inverse row to be checked
    betai:	the basis inverse row

  Returns: TRUE if the values match, FALSE otherwise.
*/

{ int i_orig,j,j_orig,kpos,k ;
  double *betak,betaik,*aj,abarij ;
  bool retval ;
  pkvec_struct *aj_pk ;

  const char *rtnnme = "check_unscalebetai" ;

  retval = TRUE ;

/*
  Allocate beta<k>. No sense using calloc, we have to clear it for each
  basis inverse column extraction. Similarly we need to clear aj for each
  original system column extraction.
*/
  betak = (double *) MALLOC((dy_sys->concnt+1)*sizeof(double)) ;
  aj = (double *) MALLOC((orig_sys->concnt+1)*sizeof(double)) ;
/*
  Open a loop to run through the basis columns. The first action is to acquire
  the unscaled column for the variable in basis position kpos. For logicals,
  we simply clear a<j> and drop a 1.0 in the right place.
*/
  for (kpos = 1 ; kpos <= dy_sys->concnt && retval == TRUE ; kpos++)
  { j = dy_basis[kpos] ;
    if (j > dy_sys->concnt)
    { j_orig = dy_actvars[j] ;
      if (consys_getcol_ex(orig_sys,j_orig,&aj) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,"column",
	       consys_nme(orig_sys,'v',j_orig,TRUE,NULL),j_orig) ;
	retval = FALSE ;
	break ; } }
    else
    { memset(aj,0,(orig_sys->concnt+1)*sizeof(double)) ;
      i_orig = dy_actcons[j] ;
      if (orig_sys->ctyp[i_orig] == contypGE)
      { aj[i_orig] = -1.0 ; }
      else
      { aj[i_orig] = 1.0 ; } }
/*
  At this point, we have an unscaled a<j>, with rows indexed in the orig_sys
  frame, and an unscaled beta<i>, with rows indexed in the dy_sys frame. To
  do the dot product, we walk beta<i>. For each nonzero entry, we look up
  the original constraint and do the multiplication.
*/
    abarij = 0 ;
    for (k = 1 ; k <= dy_sys->concnt ; k++)
    { betaik = betai[k] ;
      if (betaik == 0) continue ;
      i_orig = dy_actcons[k] ;
      abarij += betaik*aj[i_orig] ; }
    setcleanzero(abarij,dy_tols->zero) ;
    if (kpos == ipos)
    { if (!withintol(abarij,1.0,1000*dy_tols->zero))
      { errmsg(743,rtnnme,dy_sys->nme,ipos,kpos,1.0,abarij,
	       1.0-abarij,1000*dy_tols->zero) ;
	retval = FALSE ; } }
    else
    { if (!withintol(abarij,0.0,1000*dy_tols->zero))
      { errmsg(743,rtnnme,dy_sys->nme,ipos,kpos,0.0,abarij,
	       abarij,1000*dy_tols->zero) ;
	retval = FALSE ; } }
    if (retval == FALSE) break ;
/*
  Now the second test. Extract the column beta<k> of the basis inverse and then
  compare the two values of beta<ik>. Begin by clearing betak.
*/
    memset(betak,0,(dy_sys->concnt+1)*sizeof(double)) ;
/*
  If not for scaling, we'd just load 1.0 in position kpos of betak. But with
  scaling, we need to determine the original row and load the row scaling
  coefficient.
*/
    if (lcl_rowscale == NULL)
    { betak[kpos] = 1.0 ; }
    else
    { i_orig = dy_actcons[kpos] ;
      betak[kpos] = lcl_rowscale[i_orig] ; }
/*
  Ftran the vector, producing column beta<k> of the basis inverse.
*/
    dy_ftran(betak,FALSE) ;
/*
  And finish unscaling the entry in position ipos.
*/
    if (lcl_rowscale != NULL)
    { j = dy_basis[ipos] ;
	if (j <= dy_sys->concnt)
	{ i_orig = dy_actcons[j] ;
	  betak[ipos] /= lcl_rowscale[i_orig] ; }
	else
	{ j_orig = dy_actvars[j] ;
	  betak[ipos] *= lcl_colscale[j_orig] ; } }
/*
  And compare. 
*/
    if (!withintol(betak[ipos],betai[kpos],1000*dy_tols->zero))
    { errmsg(742,rtnnme,dy_sys->nme,ipos,kpos,betai[kpos],betak[ipos],
	     betai[kpos]-betak[ipos],1000*dy_tols->zero) ;
      retval = FALSE ; } }

  if (betak != NULL) FREE(betak) ;
  if (aj != NULL) FREE(aj) ;

  return (retval) ; }



#endif /* CHECK_UNSCALED_BETAI */

bool dy_unscale_betai (consys_struct *orig_sys, int j_orig,
		       double **p_betai, double **p_ai)

/*
  Generally speaking, the routine returns an unscaled row of the basis inverse.
  The exact behaviour varies with the parameters.

  j_orig != 0 : if 1 <= j_orig <= orig_sys->varcnt, j_orig specifies an
		architectural variable in orig_sys
		if 1 <= -j_orig <= orig_sys->concnt, j_orig specifies a
		logical variable and the associated constraint in orig_sys.

  For x<j_orig> (architectural or logical) active and basic in position i,
  the unscaled basis row beta<i> is returned.
  
  If x<j_orig> specifies a logical, hence a constraint i_orig = j_orig, that
  is currently inactive, the routine will calculate the row of the basis
  inverse which would result if constraint i_orig were added to the basis in
  position dy_sys->concnt+1 with x<i_orig> as the basic variable. A little
  math shows that the proper row is
     beta<i> = [ -a<i_orig,B>inv(B)  1].
  (If the constraint is a >= constraint, the row is negated.)

  j_orig == 0 : *p_ai specifies a new constraint a<i_new> not in orig_sys;
		the constraint is assumed to be <= and the coefficients are
		assumed to be unscaled

  The routine calculates beta<i> = [ -a<i_new,B>inv(B)  1]. For this case,
  a<i_new> should NOT include the new logical, and the architecturals should
  be numbered to match orig_sys. (If we were to pretend we'd already added
  the constraint and logical, the indexing of architecturals would be off by
  1 due to insertion of the logical in orig_sys. I'm not sure this is the
  right choice, but only further work with cuts will tell that tale.)

  The math associated with unscaling a basis inverse row isn't really any
  uglier than the math associated with unscaling variables or reduced costs,
  but it's much harder to explain. See the accompanying documentation for
  the linear algebra in all it's gory detail. For here, the following will
  suffice. First, we only need to unscale the elements beta_scaled<ik> where
  x<B(k)> is an architectural variable. Let k_orig be the index of x<B(k)>
  in the original system.

    * If x<j_orig> is an architectural and basic in position i,
	 beta<ik> = colscale<k_orig>*beta_scaled<ik>*rowscale<k_orig> 
    
    * If x<j_orig> is a logical and basic in position i,
	 beta<ik> = (1/rowscale<k_orig>)*beta_scaled<ik>*rowscale<k_orig> 
 
  If we're working with a new constraint a<i_new>, then the row and column
  scaling factors are both 1.0.


  Parameters:
    orig_sys:	unscaled copy of the original constraint system
    j_orig:	index of x<j_orig> in the original system
    p_betai:	(i) pointer to a vector; if NULL, one will be allocated
		(o) unscaled row of the basis inverse
    p_ai:	(i) pointer to a vector; can be NULL if j_orig != 0, in
		    which case one will be allocated if necessary; unused
		    if x<j_orig> is active and basic;
		    If j_orig == 0, then a<i> must contain a new constraint.
		(o) unscaled row a<i>

  Returns: TRUE if the calculation succeeds, FALSE otherwise.
*/

{ int i,i_orig,j,p,k,k_orig ;
  double *betai,*ai,sign ;
  bool inbasis,retval ;

  const char *rtnnme = "dy_unscale_betai" ;

# ifdef PARANOIA
  if (j_orig == 0)
  { if (p_ai == NULL)
    { errmsg(2,rtnnme,"&a<i_new>") ;
      return (FALSE) ; }
    if (*p_ai == NULL)
    { errmsg(2,rtnnme,"a<i_new>") ;
      return (FALSE) ; } }
# endif

  retval = TRUE ;
  betai = NULL ;
  ai = NULL ;

/*
  Do we have a variable that's active and basic, or will we need to expand the
  basis? If j_orig < 0 it's a logical, so we need to ask whether or not the
  associated constraint and logical are active. In no event do we deal with an
  active, nonbasic variable.
*/
  inbasis = TRUE ;
  if (j_orig > 0)
  { j = dy_origvars[j_orig] ;
    if (j <= 0)
    { errmsg(737,rtnnme,orig_sys->nme,
	     consys_nme(orig_sys,'v',j_orig,FALSE,NULL),j_orig) ;
      return (FALSE) ; }
    i = dy_var2basis[j] ;
    if (i <= 0)
    { errmsg(380,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	     dy_prtvstat(dy_status[j]),"basic") ;
      return (FALSE) ; }
    i_orig = dy_actcons[i] ; }
  else
  if (j_orig < 0)
  { i_orig = -j_orig ;
    i = dy_origcons[i_orig] ;
    if (i <= 0)
    { inbasis = FALSE ;
      j = 0 ; }
    else
    { j = i ;
      i = dy_var2basis[j] ;
      if (i <= 0)
      { errmsg(380,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	       dy_prtvstat(dy_status[j]),"basic") ;
	return (FALSE) ; } } }
  else
  { inbasis = FALSE ;
    j = 0 ;
    i = 0 ;
    i_orig = 0 ; }
/*
  If we need to allocate betai, always make it large enough to allow
  expansion of the basis by one row/column.
*/
  if (*p_betai == NULL)
  { betai = (double *) CALLOC(dy_sys->concnt+2,sizeof(double)) ; }
  else
  { betai = *p_betai ; }
/*
  If we need to expand the basis, then either we'll need to recover the
  unscaled constraint from orig_sys (i_orig != 0) or the client must have
  supplied one as a<i_new>. sign will handle the negation of the basis inverse
  row in the case of a >= constraint.
*/
  sign = quiet_nan(0) ;
  if (inbasis == FALSE)
  { if (i_orig == 0)
    { ai = *p_ai ;
      sign = 1.0 ; }
    else
    { ai = NULL ;
      if (consys_getrow_ex(orig_sys,i_orig,&ai) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,"row",
	       consys_nme(orig_sys,'c',i_orig,TRUE,NULL),i_orig) ;
	if (betai != NULL) FREE(betai) ;
	if (ai != NULL) FREE(ai) ;
	return (FALSE) ; }
      if (orig_sys->ctyp[i_orig] == contypGE)
      { sign = -1.0 ; }
      else
      { sign = 1.0 ; } } }
/*
  The easy case: x<j> is active in basis pos'n i. We'll be able to extract it
  by btran'ing a unit vector. Apply the relevant scaling factor for the row
  (see the explanation at the head of the routine).
*/
  if (inbasis == TRUE)
  { if (lcl_colscale == NULL)
    { betai[i] = 1.0 ; }
    else
    if (j <= dy_sys->concnt)
    { if (i_orig == 0)
      { betai[i] = 1.0 ; }
      else
      { betai[i] = lcl_rowscale[i_orig] ; } }
    else
    { betai[i] = lcl_colscale[j_orig] ; } }
/*
  The harder case --- this row is not in the active basis. We already have
  the unscaled constraint in a<i> (either a<i_orig> or a<i_new>). As
  mentioned above, we want to btran a vector with column scaling in place in
  order to cancel the inv(S) scaling present in the basis inverse. Do this as
  we copy the basic coefficients to beta<i>. Again, scaling for logicals is
  always 1.0. Duplicate the loop so that the test for scaling isn't in the
  inner loop. The final detail is to add +/-1 in the expansion position.
*/
  else
  { if (lcl_colscale != NULL)
    { for (p = 1 ; p <= dy_sys->concnt ; p++)
      { k = dy_basis[p] ;
	if (k > dy_sys->concnt)
	{ k_orig = dy_actvars[k] ;
	  betai[p] = -ai[k_orig]*lcl_colscale[k_orig]*sign ; } } }
    else
    { for (p = 1 ; p <= dy_sys->concnt ; p++)
      { k = dy_basis[p] ;
	if (k > dy_sys->concnt)
	{ k_orig = dy_actvars[k] ;
	  betai[p] = -ai[k_orig]*sign ; } } }
    betai[p] = sign ; }
/*
  btran our vector.
*/
  dy_btran(betai) ;
/*
  If scaling is active, beta<i> needs one more unscaling step to remove the
  row scaling present in the basis inverse. Walk the basis and unscale each
  coefficient when the basis position corresponds to an architectural
  variable.
*/
  if (lcl_rowscale != NULL)
  { for (p = 1 ; p <= dy_sys->concnt ; p++)
    { if (betai[p] == 0) continue ;
      j = dy_basis[p] ;
      k_orig = dy_actcons[p] ;
      betai[p] *= lcl_rowscale[k_orig] ; } }
/*
  That should do it. Return the results.
*/
  *p_betai = betai ;
  if (inbasis == FALSE) *p_ai = ai ;

# ifdef CHECK_UNSCALED_BETAI
  if (inbasis == TRUE)
  { retval = check_unscalebetai(orig_sys,i,betai) ; }
# endif

  return (retval) ; }



extern void dy_unscale_soln (double *x, double *y)

/*
  This routine unscales the primal and dual variable values before returning
  them to the user. The necessary unscaling is as follows:

    primal architectural:	x<j>S<j>
    primal logical:		x<i>/R<i>

    dual:			y<i>R<i>

  Parameters:
    x:	basic primal variables
    y:	dual variables

  Returns: undefined.
*/

{ int i,j,i_orig,j_orig ;
  double xi,yi ;

/*
  Did we scale? If not, return right off.
*/
  if (lcl_colscale == NULL) return ;
/*
  Since we're only dealing with basic primal variables, it suffices to step
  through the constraints (equivalently, basis positions).
*/
  for (i = 1 ; i <= dy_sys->concnt ; i++)
  { i_orig = dy_actcons[i] ;
    j = dy_basis[i] ;
    xi = x[i] ;
    if (j <= dy_sys->concnt)
      xi /= lcl_rowscale[i_orig] ;
    else
    { j_orig = dy_actvars[j] ;
      xi *= lcl_colscale[j_orig] ; }
    setcleanzero(xi,dy_tols->zero) ;
    x[i] = xi ;
    
    yi = y[i] ;
    yi *= lcl_rowscale[i_orig] ;
    setcleanzero(yi,dy_tols->cost) ;
    y[i] = yi ; }

  return ; }



bool dy_unscaleabarj (consys_struct *orig_sys, int j_orig, double **p_abarj)

/*
  This routine returns the unscaled ftran'd column inv(B)a<j>, given an
  unscaled a<j> and a scaled inv(B).  Check the written documentation to get
  a good handle on the math. The relevant scaling is:

  x<j> architectural basic in pos'n i: abar<ij> = S<k> sc_abar<ij> (1/S<j>)

  s<i> logical basic in pos'n i:       abar<ij> = (1/R<i>) sc_abar<ij> (1/S<j>)

  To cancel a factor of inv(R) attached to the scaled basis inverse, it'll
  be convenient to apply row scaling to an unscaled column prior to doing
  the ftran.

  This routine might be generally useful, but for the nonce it's just a
  debugging routine, used to check the values abar<ij> calculated in
  dy_penalty.c:pricedualpiv using an unscaled beta<i> and an unscaled a<j>.
  (In effect, a check on the unscaling of beta<i>.)

  Parameters:
    orig_sys:	unscaled copy of the original constraint system
    j_orig:	column index in original system; negative values are assumed
		to specify logicals as the negative of the index of the
		associated constraint.
    p_abarj:	(i) vector to hold abar<j>; if NULL, one will be allocated
		(o) inv(B)a<j>, unscaled

  Returns: TRUE if the calculation is successful, FALSE otherwise.
*/

{ int i,i_orig,j,k ;
  double *abarj ;
  pkvec_struct *aj_pk ;

  const char *rtnnme = "dy_unscaleabarj" ;

/*
  Acquire the unscaled column.
*/
  if (j_orig > 0)
  { aj_pk = NULL ;
    if (consys_getcol_pk(orig_sys,j_orig,&aj_pk) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,"column",
	     consys_nme(orig_sys,'v',j_orig,TRUE,NULL),j_orig) ;
      if (aj_pk != NULL) pkvec_free(aj_pk) ;
      return (FALSE) ; } }
  else
  { aj_pk = pkvec_new(1) ;
    aj_pk->coeffs[0].ndx = -j_orig ;
    if (orig_sys->ctyp[-j_orig] == contypGE)
    { aj_pk->coeffs[0].val = -1.0 ; }
    else
    { aj_pk->coeffs[0].val = 1.0 ; }
    aj_pk->cnt = 1 ; }
/*
  Acquire space for abarj, if necessary, else clear the client's vector.
*/
  if (*p_abarj == NULL)
  { abarj = (double *) CALLOC((dy_sys->concnt+1),sizeof(double)) ; }
  else
  { abarj = *p_abarj ;
    memset(abarj,0,(dy_sys->concnt+1)*sizeof(double)) ; }
/*
  Walk the column, placing scaled coefficients for the active rows into
  abarj.
*/
  if (lcl_rowscale != NULL)
  { for (k = 0 ; k < aj_pk->cnt ; k++)
    { i_orig = aj_pk->coeffs[k].ndx ;
      i = dy_origcons[i_orig] ;
      if (i > 0)
      { abarj[i] = lcl_rowscale[i_orig]*aj_pk->coeffs[k].val ; } } }
  else
  { for (k = 0 ; k < aj_pk->cnt ; k++)
    { i_orig = aj_pk->coeffs[k].ndx ;
      i = dy_origcons[i_orig] ;
      if (i > 0)
      { abarj[i] = aj_pk->coeffs[k].val ; } } }
  pkvec_free(aj_pk) ;
/*
  Do the ftran.
*/
  dy_ftran(abarj,FALSE) ;
/*
  Now walk the basis and apply the proper unscaling factor for each basis
  position.
*/
  if (lcl_rowscale != NULL)
  { for (k = 1 ; k <= dy_sys->concnt ; k++)
    { j = dy_basis[k] ;
      if (j <= dy_sys->concnt)
      { i_orig = dy_actcons[j] ;
	abarj[k] /= lcl_rowscale[i_orig] ; }
      else
      { j_orig = dy_actvars[j] ;
	abarj[k] *= lcl_colscale[j_orig] ; } } }
/*
  Set the return vector and we're done.
*/
  *p_abarj = abarj ;

  return (TRUE) ; }


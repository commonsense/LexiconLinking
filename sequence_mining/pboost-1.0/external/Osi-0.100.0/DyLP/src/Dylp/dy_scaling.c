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
  This file contains routines which handle scaling of the original constraint
  system. Intimately tied to this is the question of whether dylp will make a
  local copy of the original constraint system or simply refer to the system
  supplied by the client.

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
    * Even if scaling is forbidden, if the constraint system contains >=
      constraints dylp will make use of scaling to convert them to <=
      constraints. (The rest of the code assumes no >= constraints.)

  Otherwise, no local copy is required, and dylp can reference the original
  system supplied by the client. THIS IS A HAZARD! It's also a feature. From an
  efficiency standpoint, it's a clear win. But dylp attaches vectors to the
  original system. If the client does something radical with it (deletes it,
  for instance), all bets are off.

  The net effect is that the rest of dylp is completely unaware that scaling
  has ever happened. In particular, logical variables are inserted with the
  coefficient +/- 1.0, as appropriate. This isn't an issue for the rest of
  dylp, but it *is* an issue if we ever have to supply unscaled results to
  the outside world (e.g., the routines that generate a solution, or the
  tableau routines). To make this work, if column i represents the logical
  for constraint i, colscale[i] must be 1/rowscale[i]. The colscale array is
  not physically lengthened to make this happen.  When processing columns,
  it's necessary to check whether the column requested corresponds to a
  logical and scale accordingly.

*/


#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_scaling.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_scaling.c 269 2009-04-02 05:38:19Z lou $" ;

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
  The scaling vectors. These may have been supplied by the client, in which
  case they are owned by client_sys, or we may have created them here, in which
  case they are owned by local_sys.
*/

static double *lcl_rowscale,*lcl_colscale ;



/*
  A few utility routines for use in other files that do scaling and unscaling.
*/

bool dy_isscaled (void)
/*
  Returns: TRUE if scaling is active, FALSE otherwise.
*/

{ return (!(lcl_rowscale == NULL && lcl_colscale == NULL)) ; }


void dy_scaling_vectors (const double **rscale, const double **cscale)
/*
  Exports the scaling vectors.

  Parameters:
    rscale:	(o) the row scaling vector
    cscale:	(o) the column scaling vector

  Returns: undefined
*/

{ *rscale = lcl_rowscale ;
  *cscale = lcl_colscale ;

  return ; }


consys_struct *dy_scaled_origsys ()
/*
  This routine exposes the scaled original system.

  Parameters: none

  Returns: the scaled original system, or NULL if no scaled local copy exists.
*/

{ return (local_sys) ; }
    




bool dy_initlclsystem (lpprob_struct *orig_lp, bool hotstart)

/*
  This routine looks at the constraint system and options provided by the
  client and decides whether dylp can reference the supplied constraint
  system directly or whether a local copy should be made. A local copy is
  required if the constraint system is to be scaled or if the constraint
  system contains >= constraints. Options specified by the client can force
  the creation of a local copy, and force or (mostly) forbid the use of
  scaling. If scaling is allowed but not forced, this routine will evaluate
  the constraint matrix and apply scaling if necessary.

  Scaling vectors are attached to the local copy of the original system. This
  is not really necessary, but it's cheap insurance against the day in the
  future when I decide to do something that actually changes the local copy of
  the original system.

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

{ int i,scalefactor,orig_gecnt ;
  bool localcopy,scale,pmone ;
  flags scaled_vecs ;
  double orig_scm,scaled_scm ;

  consys_struct *orig_sys ;

  const char *rtnnme = "dy_initlclsystem" ;

# ifdef DYLP_PARANOIA
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
  An error return is possible only if we're paranoid.
*/
  orig_sys = orig_lp->consys ;
  if (consys_evalsys(orig_sys,&orig_scm,&orig_gecnt) == FALSE)
  { errmsg(138,rtnnme,orig_sys->nme) ;
    return (FALSE) ; }
  if (orig_gecnt > 0)
  { pmone = TRUE ; }
  else
  { pmone = FALSE ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.scaling >= 2 && pmone == TRUE)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	    "\n      [%s]: found %d '>=' inequalities; at least +/-1 scaling.",
	    orig_sys->nme,orig_gecnt) ; }
# endif
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
  { case 0: /* numeric scaling prohibited */
    { scale = FALSE ;
      break ; }
    case 1: /* scale with user-supplied matrices */
    { localcopy = TRUE ;
      scale = TRUE ;
#     ifdef DYLP_PARANOIA
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
		  "\n      [%s]: scaling for numeric accuracy not required;",
		  orig_sys->nme) ; }
	dyio_outfmt(dy_logchn,dy_gtxecho," %g <= |a<ij>| <= %g, metric = %g.",
		    orig_sys->minaij,orig_sys->maxaij,orig_scm) ; }
#     endif
      break ; }
    default:
    { errmsg(7,rtnnme,__LINE__,"scaling option code",dy_opts->scaling) ;
      return (FALSE) ; } }
# ifdef DYLP_PARANOIA
  if (scale == TRUE && localcopy == FALSE)
  { errmsg(1,rtnnme,__LINE__) ;
    return (FALSE) ; }
# endif
  if (pmone == TRUE)
  { localcopy = TRUE ; }
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
  Do we need to scale for numeric accuracy? If the client provided us with
  scaling vectors, use them. We can use them directly if there's no local
  copy. Given a local copy, transfer the vectors to the local system.
*/
  if (scale == TRUE)
  { 
    if (dy_opts->scaling == 1)
    { if (localcopy == TRUE)
      { memcpy(local_sys->rowscale,client_sys->rowscale,
	       ((size_t) (local_sys->concnt+1)*sizeof(double))) ;
        memcpy(local_sys->colscale,client_sys->colscale,
	       ((size_t) (local_sys->varcnt+1)*sizeof(double))) ;
	lcl_rowscale = local_sys->rowscale ;
	lcl_colscale = local_sys->colscale ; }
      else
      { lcl_rowscale = client_sys->rowscale ;
	lcl_colscale = client_sys->colscale ; } }
/*
  If no scaling vectors were supplied, call consys_geomscale and
  consys_equiscale to calculate the scaling vectors, and store them as the
  active vectors.
*/
    else
    { if (consys_geomscale(local_sys,&local_sys->rowscale,
					  &local_sys->colscale) == FALSE)
      { errmsg(135,rtnnme,local_sys->nme) ;
	return (FALSE) ; }
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
#     endif
      lcl_rowscale = local_sys->rowscale ;
      lcl_colscale = local_sys->colscale ; }
/*
  If we have a local copy, attach it.
*/
    if (localcopy == TRUE)
    { if (consys_attach(local_sys,CONSYS_RSCALE,sizeof(double),
			(void **) &lcl_rowscale) == FALSE)
      { errmsg(100,rtnnme,local_sys->nme,
	       consys_assocnme(NULL,CONSYS_RSCALE)) ;
	return (FALSE) ; }
      if (consys_attach(local_sys,CONSYS_CSCALE,sizeof(double),
			(void **) &lcl_colscale) == FALSE)
      { errmsg(100,rtnnme,local_sys->nme,
	       consys_assocnme(NULL,CONSYS_CSCALE)) ;
	return (FALSE) ; } } }
/*
  Do we need to scale by -1 to convert >= constraints to <= constraints? If
  we need to do this, we're guaranteed to have a local copy of the constraint
  system, but we may not have scaling vectors yet. The attach will initialise
  them to 1.0.
*/
  if (pmone == TRUE)
  { if (scale == FALSE)
    { if (lcl_rowscale == NULL)
      { if (consys_attach(local_sys,CONSYS_RSCALE,sizeof(double),
			  (void **) &lcl_rowscale) == FALSE)
	{ errmsg(100,rtnnme,local_sys->nme,
		 consys_assocnme(NULL,CONSYS_RSCALE)) ;
	  return (FALSE) ; }
	local_sys->rowscale = lcl_rowscale ;}
      if (lcl_colscale == NULL)
      { if (consys_attach(local_sys,CONSYS_CSCALE,sizeof(double),
			  (void **) &lcl_colscale) == FALSE)
	{ errmsg(100,rtnnme,local_sys->nme,
		 consys_assocnme(NULL,CONSYS_CSCALE)) ;
	  return (FALSE) ; }
	local_sys->colscale = lcl_colscale ; } }
    for (i = 1 ; i <= orig_sys->concnt ; i++)
    { if (orig_sys->ctyp[i] == contypGE)
      { lcl_rowscale[i] *= -1.0 ; } } }
/*
  Apply the scaling vectors and report the result. This call will actually
  convert the constraint type, along with applying the scaling factors.
*/
  if (scale == TRUE || pmone == TRUE)
  { if (consys_applyscale(local_sys,pmone,lcl_rowscale,lcl_colscale) == FALSE)
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

# ifdef DYLP_PARANOIA
  const char *rtnnme = "dy_refreshlclsystem" ;
# endif

/*
  If there's no scaling, nothing needs to be done.
*/
  if (local_sys == NULL) return ;

# ifdef DYLP_PARANOIA
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
  client_sys = NULL ;
  if (freesys == FALSE) return ;
/*
  Free the constraint system. The scaling vectors are associated with
  some constraint system, so don't need to free them here, but we do need to
  clear the local pointers to be safe.
*/
  lcl_rowscale = NULL ;
  lcl_colscale = NULL ;
  if (local_sys != NULL)
  { consys_free(local_sys) ;
    local_sys = NULL ; }

  return  ; }


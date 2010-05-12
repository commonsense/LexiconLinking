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
  This file contains routines which attempt to deal with an unbounded simplex
  problem. For both primal and dual simplex, this entails adding constraints.
  The goal, in each case, is to find constraints that will bound the problem
  without losing feasibility.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_bound.c	4.6	11/11/04" ;
static char svnid[] UNUSED = "$Id: dy_bound.c 94 2006-06-29 23:06:51Z lou $" ;



# if 0  /* debug routine */

int dbg_scanPrimConBndAct (consys_struct *orig_sys, int act_j,
			   int **p_ocndxs)

/*
  An modified version of scanPrimConBndAct that simply scans all the
  constraints in orig_sys.

  This routine scans the original constraint system looking for constraints
  that can bound motion in the direction -eta<j>delta<j>, where delta<j> is
  the change in the variable x<j> and eta<j> is the jth column of the matrix
  trans(inv(B)N -I). (trans(*) is matrix transpose.)

  This derives from the relation
    trans(x<B> x<N>) = trans(inv(B)b l/u) - trans(inv(B)N -I)deltax<N>
  where l/u is the upper or lower bound, as appropriate, for the nonbasic
  variables.
  
  When we choose the entering variable x<j>, deltax<N> becomes a vector of
  zeros, with delta<j> in the proper position. This selects column j of
  trans(inv(B)N -I) = eta<j>.

  Parameters:
    orig_sys:	The original constraint system
    act_j:	index (in dy_sys) of the offending column; negated if the
		variable is decreasing
    p_ocndxs:	(i) empty vector to hold constraint indices; assumed
		    sufficiently large if non-NULL; if *p_ocndxs = NULL,
		    allocated if necessary; if p_ocndxs = NULL, nothing
		    is returned.
		(o) indices of constraints to be activated; may not be
		    allocated if no constraints are identified

  Returns: number of candidates for activation, -1 if error.
*/

{ int i,j,k,bpos,m,n,act_m,actcnt,cand_limit,dir,retval,save_print ;
  int *ocndxs ;
  double *abarj ;
  double *orig_x,*orig_rhs,*orig_rhslow,*orig_vub,*orig_vlb,*orig_etaj ;
  double idotj,lhsi,rhsi,rhslowi ;
  contyp_enum *orig_ctyp,ctypi ;
  flags statj ;
  bool activate ;

  const char *rtnnme = "dgb_scanPrimConBndAct" ;

# ifdef PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (-1) ;  }
# endif

  save_print = dy_opts->print.conmgmt ;
  dy_opts->print.conmgmt = 4 ;

/*
  Set the multiplier for eta<j>. If x<j> is increasing, we'd normally multiply
  by -1. To compensate for decreasing (negative) motion, change dir to +1.
*/
  if (act_j < 0)
  { dir = 1 ;
    act_j = -act_j ; }
  else
  { dir = -1 ; }

# ifdef PARANOIA
  if (act_j < 1 || act_j > dy_sys->varcnt)
  { errmsg(102,rtnnme,"active variable",act_j,1,dy_sys->varcnt) ;
    return (-1) ; }
# endif

/*
  The first thing to do is to calculate abar<j> = inv(B)a<j>.
*/
  abarj = NULL ;
  if (consys_getcol_ex(dy_sys,act_j,&abarj) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,
	   "column",consys_nme(dy_sys,'v',act_j,TRUE,NULL),act_j) ;
    if (abarj != NULL) FREE(abarj) ;
    return (-1) ; }
  dy_ftran(abarj,FALSE) ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    eta<j> for %s (%d) %s:",
		consys_nme(dy_sys,'v',act_j,FALSE,NULL),act_j,
		(dir < 0)?"increasing":"decreasing") ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%8s%20s%16s","pos'n","var (ndx)","eta<ij>") ;
    for (bpos = 1 ; bpos <= dy_sys->concnt ; bpos++)
    { if (abarj[bpos] != 0)
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n%8d%14s (%3d)%16.8g",bpos,
		    consys_nme(dy_sys,'v',dy_basis[bpos],FALSE,NULL),
		    dy_basis[bpos],-abarj[bpos]) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n%8s%14s (%3d)%16.8g","n/a",
	        consys_nme(dy_sys,'v',act_j,FALSE,NULL),act_j,1.0) ; }
# endif
/*
  Now load abar<j> into a vector orig_eta<j> that we can use directly to form
  dot(orig_a<i>,orig_eta<j>) in the original system. If x<j> is decreasing,
  we need to negate eta<j>, which is handled by multiplying by dir. Remember
  that logicals do not exist in the original system.
*/
  retval = 0 ;
  m = orig_sys->concnt ;
  n = orig_sys->varcnt ;
  act_m = dy_sys->concnt ;
  orig_etaj = (double *) CALLOC((n+1),sizeof(double)) ;
  for (bpos = 1 ; bpos <= act_m ; bpos++)
  { k = dy_basis[bpos] ;
    if (k > act_m)
    { j = dy_actvars[k] ;
#     ifdef PARANOIA
      if (j < 1 || j > n)
      { errmsg(102,rtnnme,"original variable",j,1,n) ;
	retval = -1 ;
	break ; }
#     endif
      orig_etaj[j] = dir*abarj[bpos] ; } }
  if (act_j > act_m)
  { j = dy_actvars[act_j] ;
    orig_etaj[j] = -1.0*dir ; }
  if (abarj != NULL) FREE(abarj) ;
# ifdef PARANOIA
  if (j < 1 || j > n)
  { errmsg(102,rtnnme,"original variable",j,1,n) ;
    retval = -1 ; }
  if (retval < 0)
  { if (orig_etaj != NULL) FREE(orig_etaj) ;
    return (-1) ; }
# endif
/*
  Similarly, form the solution vector in terms of the original system.
*/
  orig_vub = orig_sys->vub ;
  orig_vlb = orig_sys->vlb ;
  orig_x = (double *) CALLOC((n+1),sizeof(double)) ;
  for (j = 1 ; j <= n ; j++)
  { k = dy_origvars[j] ;
    if (k > 0)
    { 
#     ifdef PARANOIA
      if (k <= act_m || k > dy_sys->varcnt)
      { errmsg(102,rtnnme,"original variable",j,act_m+1,dy_sys->varcnt) ;
	retval = -1 ;
	break ; }
#     endif
      orig_x[j] = dy_x[k] ; }
    else
    { statj = (flags) -k ;
#     ifdef PARANOIA
      if (flgoff(statj,vstatNONBASIC|vstatNBFR))
      { errmsg(433,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "inactive",consys_nme(orig_sys,'v',j,TRUE,NULL),j,
	       dy_prtvstat(statj)) ;
	retval = -1 ;
	break ; }
#     endif
      if (flgon(statj,vstatNBUB))
      { orig_x[j] = orig_vub[j] ; }
      else
      if (flgon(statj,vstatNBLB|vstatNBFX))
      { orig_x[j] = orig_vlb[j] ; } } }
# ifdef PARANOIA
  if (retval < 0)
  { if (orig_etaj != NULL) FREE(orig_etaj) ;
    if (orig_x != NULL) FREE(orig_x) ;
    return (-1) ; }
# endif
/*
  Does the client want indices returned? Did the client supply a vector for
  candidate indices? If not, make one.
*/
  cand_limit = orig_sys->concnt ;
  if (p_ocndxs == NULL || *p_ocndxs == NULL)
  { ocndxs = (int *) MALLOC(cand_limit*sizeof(int)) ; }
  else
  { ocndxs = *p_ocndxs ; }
/*
  Now we can step through the constraints in the original system. For each
  inactive constraint, we first check idotj = dot(orig_a<i>,orig_eta<j>) to
  see if we have a bounding candidate.  For a <= constraint, we need idotj >
  0; for an equality or range constraint, idotj != 0 is sufficient.
*/
  orig_ctyp = orig_sys->ctyp ;
  orig_rhs = orig_sys->rhs ;
  orig_rhslow = orig_sys->rhslow ;
  actcnt = 0 ;
  for (i = 1 ; i <= m && actcnt <= cand_limit ; i++)
  { 
    idotj = consys_dotrow(orig_sys,i,orig_etaj) ;
    ctypi = orig_ctyp[i] ;
    setcleanzero(idotj,dy_tols->zero) ;
    if (idotj == 0 || (ctypi == contypLE && idotj < 0))
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.conmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    skipping %s %s (%d), dot(a<i>,eta<j>) = %g, ",
		    consys_prtcontyp(ctypi),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i,idotj) ; }
#     endif
      continue ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.conmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    considering %s %s (%d), dot(a<i>,eta<j>) = %g, ",
		  consys_prtcontyp(ctypi),
		  consys_nme(orig_sys,'c',i,FALSE,NULL),i,idotj) ; }
#   endif
/*
  We have a bounding candidate. Is it violated at the current solution?
*/
    lhsi = consys_dotrow(orig_sys,i,orig_x) ;
    setcleanzero(lhsi,dy_tols->zero) ;
    rhsi = orig_rhs[i] ;
    if (ctypi == contypRNG)
    { rhslowi = orig_rhslow[i] ; }
    else
    if (ctypi == contypEQ)
    { rhslowi = rhsi ; }
    else
    { rhslowi = -dy_tols->inf ; }
    if (abovebnd(lhsi,rhsi) || belowbnd(lhsi,rhslowi))
    { activate = FALSE ; }
    else
    { activate = TRUE ; }

#   ifndef DYLP_NDEBUG
    if (activate == TRUE)
    { if (dy_opts->print.conmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    queued %s %s (%d), %g <= %g <= %g.",
		    consys_prtcontyp(orig_ctyp[i]),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i,
		    rhslowi,lhsi,rhsi) ; } }
    else
    { if (dy_opts->print.conmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    skipping %s constraint %s (%d),",
		    consys_prtcontyp(orig_ctyp[i]),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i) ;
	if (abovebnd(lhsi,rhsi))
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " lhs - rhs = %g - %g = %g, tol %g.",
		      lhsi,rhsi,lhsi-rhsi,dy_tols->zero*(1+fabs(rhsi))) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " rhslow - lhs = %g - %g = %g, tol %g.",
		      rhslowi,lhsi,rhslowi-lhsi,
		      dy_tols->zero*(1+fabs(rhslowi))) ; } } }
#   endif

    if (activate == TRUE) ocndxs[actcnt++] = i ; }
  if (orig_etaj != NULL) FREE(orig_etaj) ;
  if (orig_x != NULL) FREE(orig_x) ;
/*
  If we supplied ocndxs and found no candidates to activate, free it.
*/
  if (p_ocndxs != NULL)
  { if (*p_ocndxs == NULL)
    { if (actcnt == 0)
      { FREE(ocndxs) ; }
      else
      { *p_ocndxs = ocndxs ; } } }
  else
  { FREE(ocndxs) ; }


# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d constraints for activation.",actcnt) ; }
# endif

  dy_opts->print.conmgmt = save_print ;

  return (actcnt) ; }

#endif /* debug routine */



/* Routines to bound an unbounded primal problem. */


static int scanPrimConBndAct (consys_struct *orig_sys,
			      int act_j, int **p_ocndxs)

/*
  This routine scans the original constraint system looking for constraints
  that can bound motion in the direction -eta<j>delta<j>, where delta<j> is
  the change in the variable x<j> and eta<j> is the jth column of the matrix
  trans(inv(B)N -I). (trans(*) is matrix transpose.)

  This derives from the relation
    trans(x<B> x<N>) = trans(inv(B)b l/u) - trans(inv(B)N -I)deltax<N>
  where l/u is the upper or lower bound, as appropriate, for the nonbasic
  variables.
  
  When we choose the entering variable x<j>, deltax<N> becomes a vector of
  zeros, with delta<j> in the proper position. This selects column j of
  trans(inv(B)N -I) = eta<j>.

  Parameters:
    orig_sys:	The original constraint system
    act_j:	index (in dy_sys) of the offending column; negated if the
		variable is decreasing
    p_ocndxs:	(i) empty vector to hold constraint indices; assumed
		    sufficiently large if non-NULL; if NULL, allocated if
		    necessary
		(o) indices of constraints to be activated; may not be
		    allocated if no constraints are identified

  Returns: number of candidates for activation, -1 if error.
*/

{ int i,j,k,bpos,m,n,act_m,actcnt,cand_limit,dir,retval ;
  int *ocndxs ;
  double *abarj ;
  double *orig_x,*orig_rhs,*orig_rhslow,*orig_vub,*orig_vlb,*orig_etaj ;
  double idotj,lhsi,rhsi,rhslowi ;
  contyp_enum *orig_ctyp,ctypi ;
  flags statj ;
  bool activate ;

  const char *rtnnme = "scanPrimConBndAct" ;

# ifdef PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (-1) ;  }
  if (p_ocndxs == NULL)
  { errmsg(2,rtnnme,"&ocndxs") ;
    return (-1) ; }
# endif

/*
  Set the multiplier for eta<j>. If x<j> is increasing, we'd normally multiply
  by -1. To compensate for decreasing (negative) motion, change dir to +1.
*/
  if (act_j < 0)
  { dir = 1 ;
    act_j = -act_j ; }
  else
  { dir = -1 ; }

# ifdef PARANOIA
  if (act_j < 1 || act_j > dy_sys->varcnt)
  { errmsg(102,rtnnme,"active variable",act_j,1,dy_sys->varcnt) ;
    return (-1) ; }
# endif

/*
  The first thing to do is to calculate abar<j> = inv(B)a<j>.
*/
  abarj = NULL ;
  if (consys_getcol_ex(dy_sys,act_j,&abarj) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,
	   "column",consys_nme(dy_sys,'v',act_j,TRUE,NULL),act_j) ;
    if (abarj != NULL) FREE(abarj) ;
    return (-1) ; }
  dy_ftran(abarj,FALSE) ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    eta<j> for %s (%d) %s:",
	        consys_nme(dy_sys,'v',act_j,FALSE,NULL),act_j,
	        (dir < 0)?"increasing":"decreasing") ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%8s%20s%16s","pos'n","var (ndx)","eta<ij>") ;
    for (bpos = 1 ; bpos <= dy_sys->concnt ; bpos++)
    { if (abarj[bpos] != 0)
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n%8d%14s (%3d)%16.8g",bpos,
		    consys_nme(dy_sys,'v',dy_basis[bpos],FALSE,NULL),
		    dy_basis[bpos],-abarj[bpos]) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n%8s%14s (%3d)%16.8g","n/a",
	        consys_nme(dy_sys,'v',act_j,FALSE,NULL),act_j,1.0) ; }
# endif
/*
  Now load abar<j> into a vector orig_eta<j> that we can use directly to form
  dot(orig_a<i>,orig_eta<j>) in the original system. If x<j> is decreasing,
  we need to negate eta<j>, which is handled by multiplying by dir. Remember
  that logicals do not exist in the original system.
*/
  retval = 0 ;
  m = orig_sys->concnt ;
  n = orig_sys->varcnt ;
  act_m = dy_sys->concnt ;
  orig_etaj = (double *) CALLOC((n+1),sizeof(double)) ;
  for (bpos = 1 ; bpos <= act_m ; bpos++)
  { k = dy_basis[bpos] ;
    if (k > act_m)
    { j = dy_actvars[k] ;
#     ifdef PARANOIA
      if (j < 1 || j > n)
      { errmsg(102,rtnnme,"original variable",j,1,n) ;
	retval = -1 ;
	break ; }
#     endif
      orig_etaj[j] = dir*abarj[bpos] ; } }
  if (act_j > act_m)
  { j = dy_actvars[act_j] ;
    orig_etaj[j] = -1.0*dir ; }
  if (abarj != NULL) FREE(abarj) ;
# ifdef PARANOIA
  if (j < 1 || j > n)
  { errmsg(102,rtnnme,"original variable",j,1,n) ;
    retval = -1 ; }
  if (retval < 0)
  { if (orig_etaj != NULL) FREE(orig_etaj) ;
    return (-1) ; }
# endif
/*
  Similarly, form the solution vector in terms of the original system.
*/
  orig_vub = orig_sys->vub ;
  orig_vlb = orig_sys->vlb ;
  orig_x = (double *) CALLOC((n+1),sizeof(double)) ;
  for (j = 1 ; j <= n ; j++)
  { k = dy_origvars[j] ;
    if (k > 0)
    { 
#     ifdef PARANOIA
      if (k <= act_m || k > dy_sys->varcnt)
      { errmsg(102,rtnnme,"original variable",j,act_m+1,dy_sys->varcnt) ;
	retval = -1 ;
	break ; }
#     endif
      orig_x[j] = dy_x[k] ; }
    else
    { statj = (flags) -k ;
#     ifdef PARANOIA
      if (flgoff(statj,vstatNONBASIC|vstatNBFR))
      { errmsg(433,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "inactive",consys_nme(orig_sys,'v',j,TRUE,NULL),j,
	       dy_prtvstat(statj)) ;
	retval = -1 ;
	break ; }
#     endif
      if (flgon(statj,vstatNBUB))
      { orig_x[j] = orig_vub[j] ; }
      else
      if (flgon(statj,vstatNBLB|vstatNBFX))
      { orig_x[j] = orig_vlb[j] ; } } }
# ifdef PARANOIA
  if (retval < 0)
  { if (orig_etaj != NULL) FREE(orig_etaj) ;
    if (orig_x != NULL) FREE(orig_x) ;
    return (-1) ; }
# endif
/*
  Did the client supply a vector for candidate indices? If not, make one.
*/
  cand_limit = m-act_m ;
  if (dy_opts->con.actlim > 0)
  { cand_limit = minn(dy_opts->con.actlim,cand_limit) ; }
  if (*p_ocndxs == NULL)
  { ocndxs = (int *) MALLOC(cand_limit*sizeof(int)) ; }
  else
  { ocndxs = *p_ocndxs ; }
/*
  Now we can step through the constraints in the original system. For each
  inactive constraint, we first check idotj = dot(orig_a<i>,orig_eta<j>) to
  see if we have a bounding candidate.  For a <= constraint, we need idotj >
  0; for an equality or range constraint, idotj != 0 is sufficient.
*/
  orig_ctyp = orig_sys->ctyp ;
  orig_rhs = orig_sys->rhs ;
  orig_rhslow = orig_sys->rhslow ;
  actcnt = 0 ;
  for (i = 1 ; i <= m && actcnt <= cand_limit ; i++)
  { if (dy_origcons[i] > 0) continue ;
    idotj = consys_dotrow(orig_sys,i,orig_etaj) ;
    ctypi = orig_ctyp[i] ;
    setcleanzero(idotj,dy_tols->zero) ;
    if (idotj == 0 || (ctypi == contypLE && idotj < 0))
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.conmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    skipping %s %s (%d), dot(a<i>,eta<j>) = %g, ",
		    consys_prtcontyp(ctypi),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i,idotj) ; }
#     endif
      continue ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.conmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    considering %s %s (%d), dot(a<i>,eta<j>) = %g, ",
		  consys_prtcontyp(ctypi),
		  consys_nme(orig_sys,'c',i,FALSE,NULL),i,idotj) ; }
#   endif
/*
  We have a bounding candidate. Is it violated at the current solution?
*/
    lhsi = consys_dotrow(orig_sys,i,orig_x) ;
    setcleanzero(lhsi,dy_tols->zero) ;
    rhsi = orig_rhs[i] ;
    if (ctypi == contypRNG)
    { rhslowi = orig_rhslow[i] ; }
    else
    if (ctypi == contypEQ)
    { rhslowi = rhsi ; }
    else
    { rhslowi = -dy_tols->inf ; }
    if (abovebnd(lhsi,rhsi) || belowbnd(lhsi,rhslowi))
    { activate = FALSE ; }
    else
    { activate = TRUE ; }

#   ifndef DYLP_NDEBUG
    if (activate == TRUE)
    { if (dy_opts->print.conmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    queued %s %s (%d), %g <= %g <= %g.",
		    consys_prtcontyp(orig_ctyp[i]),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i,
		    rhslowi,lhsi,rhsi) ; } }
    else
    { if (dy_opts->print.conmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    skipping %s constraint %s (%d),",
		    consys_prtcontyp(orig_ctyp[i]),
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i) ;
	if (abovebnd(lhsi,rhsi))
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " lhs - rhs = %g - %g = %g, tol %g.",
		      lhsi,rhsi,lhsi-rhsi,dy_tols->zero*(1+fabs(rhsi))) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " rhslow - lhs = %g - %g = %g, tol %g.",
		      rhslowi,lhsi,rhslowi-lhsi,
		      dy_tols->zero*(1+fabs(rhslowi))) ; } } }
#   endif

    if (activate == TRUE) ocndxs[actcnt++] = i ; }
  if (orig_etaj != NULL) FREE(orig_etaj) ;
  if (orig_x != NULL) FREE(orig_x) ;
/*
  If we supplied ocndxs and found no candidates to activate, free it.
*/
  if (*p_ocndxs == NULL)
  { if (actcnt == 0)
    { FREE(ocndxs) ; }
    else
    { *p_ocndxs = ocndxs ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  queued %d constraints for activation.",actcnt) ; }
# endif

  return (actcnt) ; }



int dy_activateBndCons (consys_struct *orig_sys)

/*
  This routine coordinates bounding constraint activation in phase dyADDCON.
  In addition to the actual scan and activation, it sees to rebuilding the
  basis and solution. The heavy lifting is performed in scanPrimConBndAct and
  actBLogPrimCon.

  See the comments in dy_conmgmt.c for the effects on PSE and DSE norms.
  Notwithstanding, the approach taken here is to simply set the init_dse
  flag. The reason is that dylp does not have a phase transition sequence
  which adds constraints and returns to dual simplex. It's always primal ->
  act/deact constraints -> dual -> act/deact variables -> primal The DSE
  norms are not maintained through primal simplex, so there's no motivation
  to do an update here.

  Parameter:
    orig_sys:	The original constraint system

  Returns: number of constraints activated; -1 if there's an error.
*/

{ int *candidates,cand_cnt,act_j ;
  int retval ;
  bool actresult ;
  flags calcflgs ;
  dyret_enum factorresult ;

  const char *rtnnme = "dy_activateBndCons" ;

  retval = -1 ;

/*
  Call scanPrimConBndAct to return a list of candidates for activation, then
  call actBLogPrimConList to install them. Installing nothing always succeeds.
*/
  candidates = NULL ;
  act_j = dy_lp->ubnd.ndx ;
  cand_cnt = scanPrimConBndAct(orig_sys,act_j,&candidates) ;
  if (cand_cnt < 0)
  { errmsg(434,rtnnme,
	   dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "constraint","bounding activation") ;
    actresult = FALSE ; }
  else
  if (cand_cnt > 0)
  { actresult = dy_actBLogPrimConList(orig_sys,cand_cnt,candidates,NULL) ; }
  else
  { actresult = TRUE ; }
  if (candidates != NULL) FREE(candidates) ;
  if (actresult == FALSE) return (retval) ;
/*
  If we added constraints, we need to refactor and recalculate the primal and
  dual variables. Then decide on the simplex phase. If we came in with primal
  feasibility, we should still have it.
*/
  if (cand_cnt > 0)
  { dy_lp->simplex.init_dse = TRUE ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.conmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      factoring, calculating variables, ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"and checking feasibility ...") ; }
#   endif
    calcflgs = ladFACTOR|ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
    factorresult = dy_accchk(&calcflgs) ;
    switch (factorresult)
    { case dyrOK:
      case dyrPATCHED:
      { retval = cand_cnt ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	{ if (factorresult == dyrOK)
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    done.") ;
	  else
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n    patched.") ; }
#       endif
#	ifdef PARANOIA
	if (dy_lp->simplex.active == dyPRIMAL2 && factorresult == dyrOK &&
	    (flgon(calcflgs,ladPRIMFEAS)))
	{ errmsg(1,rtnnme,__LINE__) ;
	  retval = -1 ;
	  break ; }
#	endif
	if (flgoff(calcflgs,ladPRIMFEAS))
	{ dy_lp->simplex.next = dyPRIMAL2 ; }
	else
	if (flgoff(calcflgs,ladDUALFEAS))
	{ dy_lp->simplex.next = dyDUAL ; }
	else
	{ dy_lp->simplex.next = dyPRIMAL1 ; }
	break ; }
      default:
      { retval = -1 ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.conmgmt >= 3)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n    failed.") ;
#       endif
	break ; } } }
  else
  { retval = cand_cnt ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.conmgmt >= 1)
  { if (dy_opts->print.conmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho," %d activations.",cand_cnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    constraint system %s now %d x %d (%d + %d).",
	        dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
	        dy_sys->logvcnt) ; }
# endif

  return (retval) ; }




/* Routines to bound an unbounded dual problem. */


static int type1var(consys_struct *orig_sys,
		    int xindx, int diri, int oxkndx, flags xkstatus,
		    double abarik, double cbark)
/*
  This routine evaluates x<k> to see if it qualifies as a type 1 variable.
  A type 1 variable will bound the dual problem (its associated dual can be
  driven to 0 as y<i> enters) and it can be activated into the nonbasic
  partition while retaining dual feasibility.

  Parameters:
    orig_sys:	the original constraint system
    xindx:	index of the entering dual y<i> (leaving primal x<i>)
    diri:	direction of motion for y<i> (x<i>)
		  +1 to increase from 0 (to lb<i>)
		  -1 to decrease from 0 (to ub<i>)
    oxkndx:	index of candidate for activation y<k> (x<k>)
    xkstatus:	status for x<k>
    abarik:	pivot element abar<ik>
    cbark:	value of y<k> (reduced cost of x<k>)
  
  Returns: 1 if x<k> is activated,
	   0 if x<i> is not activated
	  -1 if something goes fatally wrong
*/

{ int xkndx ;
  const char *rtnnme = "type1var" ;

/*
  Will this variable be dual feasible? If not, it can't be a type 1. If it
  is dual feasible, it'll work, because we tested that cbar<k> and abar<ik>
  had the proper signs before calling type1var.
*/
  if ((flgon(xkstatus,vstatNBLB) && cbark < 0) ||
      (flgon(xkstatus,vstatNBUB) && cbark > 0))
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tt1eval: %s %s (%d) not dual feasible; cbar<k> = %g.",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx,cbark) ; }
#   endif
    return (0) ; }
/*
  x<k> satisfies the type 1 criteria. Activate it, and insert the reduced cost
  in dy_cbar.
*/
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    type 1 activation %s %s (%d), ",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"cbar<k> = %g, abar<ik> = %g.",
		  cbark,abarik) ; }
#   endif

    if (dy_actNBPrimArch(orig_sys,oxkndx) == FALSE)
    { errmsg(430,rtnnme,
	     orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "activate","variable",
	     consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
      return (-1) ; }
    xkndx = dy_origvars[oxkndx] ;
    dy_cbar[xkndx] = cbark ;
  
    return (1) ; }



static void type2eval (consys_struct *orig_sys,int xindx, int diri,
		       int oxkndx, flags xkstatus, double abarik, double cbark,
		       int *oxjndx, double *deltaj, double *cbarj, int *dirj)

/*
  This routine evaluates x<k> as a possible type 2 variable, and replaces the
  incumbent type 2 candidate if that's appropriate. We're looking for variables
  which bound the dual problem but would not be dual feasible if they were
  activated as nonbasic.

  A type 2 variable has the signs of cbar<k> and abar<ik> reversed from the
  usual conventions for a dual pivot. Because of this, the primal move
  is in the `wrong' direction. E.g., if x<i> is rising to its lower bound
  and leaving, and x<k> is NBLB with abar<ik> > 0, x<k> has to >decrease<
  in order to bring x<i> to its lower bound.

  Parameters:
    orig_sys:	the original constraint system
    xindx:	index of the entering dual y<i> (leaving primal x<i>)
    diri:	direction of motion for y<i> (x<i>)
		  +1 to increase from 0 (to lb<i>)
		  -1 to decrease from 0 (to ub<i>)
    oxkndx:	index of candidate for activation y<k> (x<k>)
    xkstatus:	status for x<k>
    abarik:	pivot element abar<ik>
    cbark:	value of y<k> (reduced cost of x<k>)
    oxjndx:	(i) if nonzero, index of the type 2 incumbent x<j>
		(o) replaced with new incumbent, if appropriate
    deltaj:	(i) dual delta for y<i> given x<j>
		(o) replaced with new delta, if x<j> is replaced
    cbarj:	(i) reduced cost of x<j>
    		(o) replaced with new reduced cost, if x<j> is replaced
    dirj:	(i) direction of motion of x<j>
		(o) replaced with new direction of motion if x<j> is replaced
  
  Returns: undefined
*/

{ double deltak ;

/*
  Will this variable be dual feasible? If so, it can't be a type 2. If it
  isn't dual feasible, it'll work, because we tested that cbar<k> and abar<ik>
  had the proper signs before calling type2var.
*/
  if ((flgon(xkstatus,vstatNBLB) && cbark > 0) ||
      (flgon(xkstatus,vstatNBUB) && cbark < 0) || flgon(xkstatus,vstatNBFR))
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tt2eval: %s %s (%d) dual feasible; cbar<k> = %g.",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx,cbark) ; }
#   endif
    return ; }
/*
  Is x<k> a better choice than the current incumbent x<j>? If so, replace it.
  It takes a little work to figure out the proper direction of motion,
  particularly if x<k> is NBFR. Remember, the primal direction of change is
  backwards from the normal motion.
*/
  deltak = fabs(cbark/abarik) ;
  if (deltak < *deltaj)
  { *oxjndx = oxkndx ;
    *deltaj = deltak ;
    *cbarj = cbark ;
    if (flgon(xkstatus,vstatNBLB))
    { *dirj = -1 ; }
    else
    if (flgon(xkstatus,vstatNBUB))
    { *dirj = 1 ; }
    else
    { if (diri > 0)
      { if (abarik > 0)
	{ *dirj = 1 ; }
	else
	{ *dirj = -1 ; } }
      else
      { if (abarik > 0)
	{ *dirj = -1 ; }
	else
	{ *dirj = 1 ; } } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      t2eval: choosing %s %s (%d), delta %g,",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',*oxjndx,TRUE,NULL),*oxjndx,*deltaj) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"cbar<j> = %g, abar<ij> = %g.",
		  *cbarj,abarik) ; }
#   endif
  }
# ifndef DYLP_NDEBUG
  else
  { if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      t2eval: skipping %s %s (%d), delta %g,",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx,deltak) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"cbar<k> = %g, abar<ik> = %g.",
		  cbark,abarik) ; } }
# endif

  return ; }



static void type3eval (consys_struct *orig_sys,
		       int xindx, int diri, int oxkndx, flags xkstatus,
		       double abarik, double cbark, int *oxjndx, double *distj,
		       double *cbarj, int *dirj, double *abarij)

/*
  This routine evaluates x<k> as a possible type 3 variable, and replaces the
  incumbent type 3 candidate if that's appropriate. We're looking for variables
  which can be activated with a bound-to-bound pivot and drive the reduced
  cost of y<i> to 0 (drive x<i> toward feasibility). We'll pick the variable
  that puts us on the correct side of the bound, with a preference to minimise
  |bnd<i> - delta<k>*abar<ik>|.

  Parameters:
    orig_sys:	the original constraint system
    xindx:	index of the entering dual y<i> (leaving primal x<i>)
    diri:	direction of motion for y<i> (x<i>)
		  +1 to increase from 0 (to lb<i>)
		  -1 to decrease from 0 (to ub<i>)
    oxkndx:	index of candidate for activation y<k> (x<k>)
    xkstatus:	status for x<k>
    abarik:	pivot element abar<ik>
    cbark:	value of y<k> (reduced cost of x<k>)
    oxjndx:	(i) if nonzero, index of the type 3 incumbent x<j>
		(o) replaced with new incumbent, if appropriate
    distj:	(i) |bnd<i> - (x<i>+delta<j>)|, with sign assigned to be
		    positive if we're within bounds, and negative if we're
		    out of bound.
		(o) replaced with new delta, if x<j> is replaced
    cbarj:	(i) reduced cost of x<j>
    		(o) replaced with new reduced cost, if x<j> is replaced
    dirj:	(i) direction of motion of x<j>
		(o) replaced with new direction of motion if x<j> is replaced
    abarij:	(i) pivot for x<j>
		(o) replaced with new pivot, if x<j> is replaced
  
  Returns: undefined
*/

{ int dirk ;
  double lbk,ubk,deltak,distk,bndi,xival ;
  bool newxj ;

/*
  Get the bounds on x<k> and calculate the possible delta. If x<k> doesn't
  have both bounds, it can't be a type 3 variable.
*/
  lbk = orig_sys->vlb[oxkndx] ;
  ubk = orig_sys->vub[oxkndx] ;
  if (lbk <= -dy_tols->inf || ubk >= dy_tols->inf)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      t3eval: skipping %s %s (%d)",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
      if (lbk <= -dy_tols->inf)
	dyio_outfmt(dy_logchn,dy_gtxecho,", lb = %g",lbk) ;
      if (ubk >= dy_tols->inf)
	dyio_outfmt(dy_logchn,dy_gtxecho,", ub = %g",ubk) ;
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
    return ; }
/*
  Now look at x<k>'s status and reduced cost, and from that decide if we can
  flip and activate it.
*/
  if (flgon(xkstatus,vstatNBLB))
  { if (cbark > 0)
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n      t3eval: skipping %s %s (%d), cbar = %g, can't flip.",
		dy_prtvstat(xkstatus),
		consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx,cbark) ; }
#     endif
      return ; }
    dirk = 1 ;
    deltak = -abarik*(ubk-lbk) ; }
  else
  { if (cbark < 0)
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n      t3eval: skipping %s %s (%d), cbar = %g, can't flip.",
	       dy_prtvstat(xkstatus),
	       consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx,cbark) ; }
#     endif
      return ; }
    dirk = -1 ;
    deltak = -abarik*(lbk-ubk) ; }
  setcleanzero(deltak,dy_tols->zero) ;
/*
  Does delta<k> move x<i> toward feasibility?
*/
  if (deltak == 0)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      t3eval: skipping %s %s (%d), delta = 0.",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ; }
#   endif
    return ; }
  if ((diri > 0 && deltak < 0) || (diri < 0 && deltak > 0))
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n      t3eval: skipping %s %s (%d), direction; delta = %g.",
	     dy_prtvstat(xkstatus),consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),
	     oxkndx,deltak) ; }
#   endif
    return ; }
/*
  The only remaining question is whether x<k> brings x<i> closer to
  feasibility than the incumbent. The preference is for feasibility. If both
  x<k> and x<j> are feasible or infeasible (sign(dist<k>) = sign(dist<j>)),
  we want to get as close to the bound as possible. If the distance is equal,
  keep the incumbent.  The sign of dist<k> is set to be negative if the new
  value of x<i> is out-of-bound in either direction.
*/
  xival = dy_x[xindx]+deltak ;
  setcleanzero(xival,dy_tols->zero) ;

  newxj = FALSE ;
  if (diri < 0)
  { bndi = dy_sys->vub[xindx] ;
    distk = bndi-xival ;
    if (belowbnd(xival,dy_sys->vlb[xindx])) distk = -distk ; }
  else
  { bndi = dy_sys->vlb[xindx] ;
    distk = xival-bndi ;
    if (abovebnd(xival,dy_sys->vub[xindx])) distk = -distk ; }
  setcleanzero(distk,dy_tols->zero) ;
  if (distk > 0 && *distj < 0)
  { newxj = TRUE ; }
  else
  if (fabs(distk) < fabs(*distj))
  { newxj = TRUE ; }

  if (newxj == TRUE)
  { *oxjndx = oxkndx ;
    *distj = distk ;
    *cbarj = cbark ;
    *dirj = dirk ;
    *abarij = abarik ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      t3eval: choosing %s %s (%d), ",
		  dy_prtvstat(xkstatus),
		  consys_nme(orig_sys,'v',*oxjndx,TRUE,NULL),*oxjndx) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"cbar = %g, delta = %g, dist = %g.",
		  *cbarj,deltak,*distj) ; }
#   endif
  }

  return ; }



static int type2activate (consys_struct *orig_sys,
			  int xindx, int diri,
			  int oxjndx, int dirj, double cbarj)

/*
  This routine performs the activate and pivot step required for a type 2
  variable. It has some error recovery capability, but makes no attempt at
  real sophistication. If things get too rough, we can always fall back on
  primal simplex.

  Because we've flipped the signs on cbar<j> and abar<ij>, the primal move
  is in the `wrong' direction. E.g., if x<i> is rising to its lower bound
  and leaving, x<j> NBLB, and abar<ij> > 0, x<j> has to >decrease< in order
  to bring x<i> to its lower bound.

  Parameters:
    orig_sys:	the original constraint system
    xindx:	index of the leaving variable x<i>
    diri:	direction of motion of x<i>
    oxjndx:	index (in orig_sys) of the entering variable x<j>
    dirj:	direction of motion of x<j>
    cbarj:	reduced cost of x<j>

  Returns: 1 if a variable is activated and pivoted into the basis without
	     error
	   0 if there's a nonfatal problem (this will cause a reversion
	     to primal simplex).
	  -1 if there's a fatal problem
*/

{ int xjndx,xkndx ;
  double abarij,deltaj ;
  bool pivoted ;
  dyret_enum pivresult,duennaresult ;
  int retval ;
  const char *rtnnme = "type2activate" ;

# ifndef DYLP_NDEBUG
  flags xjstatus ;

  if (dy_opts->print.varmgmt >= 1)
  { xjstatus = (flags) -dy_origvars[oxjndx] ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    activating and pivoting %s %s (%d), ",
	        dy_prtvstat(xjstatus),
	        consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
    dyio_outfmt(dy_logchn,dy_gtxecho," cbar = %g, %s.",
	        cbarj,(dirj < 0)?"falling":"rising") ; }
# endif
/*
  Try to activate the variable.
*/
  if (dy_actNBPrimArch(orig_sys,oxjndx) == FALSE)
  { errmsg(430,rtnnme,
	   orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	   "activate","variable",
	   consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
    return (-1) ; }
  xjndx = dy_origvars[oxjndx] ;
  dy_cbar[xjndx] = cbarj ;
/*
  Now attempt the pivot.
*/
  pivresult = dy_dualpivot(xindx,diri,&xjndx,&dirj,&cbarj,&abarij,
			   &deltaj,&xkndx) ;
  switch (pivresult)
  { case dyrOK:
    case dyrDEGEN:
    case dyrOPTIMAL:
    case dyrPUNT:
    case dyrREQCHK:
    { pivoted = TRUE ;
      break ; }
    default:
    { pivoted = FALSE ;
      break ; } }
# ifndef DYLP_NDEBUG
  if ((dy_opts->print.varmgmt >= 3) ||
      (dy_opts->print.varmgmt >= 2 && pivresult != dyrOK))
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n      pivot attempt %s, pivot return code = %s.",
	        (pivoted == TRUE)?"succeeded":"failed",
	        dy_prtdyret(pivresult)) ; }
  if (dy_opts->print.dual >= 4)
    dy_logpivot(pivresult,xjndx,dirj,cbarj,xindx,diri,abarij,deltaj) ;
# endif
/*
  Call La Duenna after the pivot as usual. It boils down to these cases:
    * The pivot went through, and La Duenna had no complaints. In this
      case, we can go back to dual simplex and try for normal pivots.
    * The pivot went through, but La Duenna had non-fatal problems. We'll
      revert to primal simplex and see if it goes better.
    * The pivot didn't go through, but La Duenna managed to salvage the
      situation. In this case we haven't managed to change the situation
      and returning to dual simplex would be futile. Again, try primal simplex.
      If `salvage the situation' meant dealing with a singular pivot, we need
      to remove the leaving variable from the pivot reject list (there will be
      no other entries, given we're in dyADDVAR).
    * There was a fatal error, either in the pivot attempt or from La Duenna.
      Return dyINV.
  If the pivot resulted in a fatal error, that'll be reflected in the return
  code from dy_duenna. We have to lie about the phase for a moment here, to
  avoid running afoul of all kinds of checks.
*/
  dy_lp->phase = dyDUAL ;
  duennaresult = dy_duenna(pivresult,xjndx,xindx,-1,-1) ;
  if (pivresult == dyrSINGULAR)
  { if (dy_clrpivrej(NULL) != TRUE) return (-1) ; }
  dy_lp->phase = dyADDVAR ;
  switch (duennaresult)
  { case dyrOK:
    case dyrRESELECT:
    case dyrOPTIMAL:
    { if (pivoted == TRUE)
      { retval = 1 ; }
      else
      { retval = 0 ; }
      break ; }
    case dyrPUNT:
    case dyrLOSTDFEAS:
    { retval = 0 ;
      break ; }
    default:
    { retval = -1 ; } }

# ifndef DYLP_NDEBUG
  if ((dy_opts->print.varmgmt >= 3) ||
      (dy_opts->print.varmgmt >= 2 && duennaresult != dyrOK))
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      La Duenna return code %s.",
	        dy_prtdyret(duennaresult)) ; }
  if (dy_opts->print.varmgmt >= 1)
  { if (dy_opts->print.varmgmt >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  ") ; }
    if (retval == 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho," 1 activation and pivot.") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    constraint system %s now %d x %d (%d + %d).",
		  dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
		  dy_sys->logvcnt) ; }
    else
    if (retval == 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  " activate & pivot failed; reverting to primal.") ; } }
# endif

  return (retval) ; }



static int type3activate (consys_struct *orig_sys, double *betai,
			  int xindx, int diri,
			  int oxjndx, int dirj,
			  double cbarj, double abarij)
/*
  This routine is responsible for handling type 3 variables. It expects to
  be passed a type 3 variable x<j> which can be activated with a
  bound-to-bound pivot.

  Once the initial variable is dealt with, we recalculate the primal variables
  and then run dualout to see if the same leaving variable x<i> is selected.
  If not, we're done, and we can return to dual simplex. Otherwise, we scan
  for another type 3 variable and repeat. If we run out of type 3 variables
  before we select a different leaving variable, we're in trouble. In theory
  that's it, but for now revert to the primal to see if it can figure it out.

  The overall flow of each pivot is much as type2activate --- we do the pivot
  (easy, it's just a matter of updating the primal variables), then call
  dy_duenna to check things over and do some bookkeeping. If anything goes
  mildly wrong, we bail out back to the primal simplex.

  In the case where a sequence of type 3 activations hasn't helped, the return
  value is a lie --- it's forced to 0. There's probably a better fix, but
  this'll do for now. I'm not entirely convinced this routine is correct.

  Parameters:
    orig_sys:	the original constraint system
    betai:	row i of the basis inverse
    xindx:	index of the leaving variable x<i>
    diri:	direction of motion of x<i>
    oxjndx:	index (in orig_sys) of the swinging variable x<j>
    dirj:	direction of motion of x<j>
    cbarj:	reduced cost of x<j>
    abarij:	pivot for x<j>
  
  Returns: > 0: one or more bound-to-bound pivots has succeeded in changing
	        the variable selected by dy_dualout; a normal dual pivot is
		possible
	     0: if we run out of type 3 variables and still get the same
		leaving variable, or if we encounter recoverable error on the
		pivot
	    -1: if things go badly wrong
*/

{ int candxi,oxkndx,xjndx,xqndx,pkndx ;
  double cbark,abarik,deltak,distj ;
  flags xjstatus,xkstatus ;
  bool fatal ;
  int actcnt,retval ;
  dyret_enum duennaresult,outresult ;
  pkvec_struct *ak ;
  pkcoeff_struct *aqk ;
  const char *rtnnme = "type3activate" ;

# ifndef DYLP_NDEBUG
  double deltaj ;
# endif

  retval = -1 ;
  actcnt = 0 ;
  ak = NULL ;

/*
  Dive right into the loop that handles activation, pivot, and selection
  of a new type 3 variable. We do this as long as we're still selecting x<i>
  to leave and have a variable x<j> to activate and pivot.
*/
  candxi = xindx ;
  while  (oxjndx != 0)
  {

#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 1)
    { xjstatus = (flags) -dy_origvars[oxjndx] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    activating and flipping %s %s (%d), ",
		  dy_prtvstat(xjstatus),
		  consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
      dyio_outfmt(dy_logchn,dy_gtxecho," cbar = %g, %s.",
		  cbarj,(dirj < 0)?"falling":"rising") ; }
#   endif
/*
  Try to activate the variable.
*/
    if (dy_actNBPrimArch(orig_sys,oxjndx) == FALSE)
    { errmsg(430,rtnnme,
	     orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "activate","variable",
	     consys_nme(orig_sys,'v',oxjndx,TRUE,NULL),oxjndx) ;
      return (-1) ; }
    actcnt++ ;
/*
  Do the bound-to-bound swing. This is inefficient, but if the algorithm
  works I can always come back and change this to an incremental update.
*/
    xjndx = dy_origvars[oxjndx] ;
    xjstatus = dy_status[xjndx] ;
    if (flgon(xjstatus,vstatNBLB))
    { dy_status[xjndx] = vstatNBUB ;
      dy_x[xjndx] = dy_sys->vub[xjndx] ; }
    else
    { dy_status[xjndx] = vstatNBLB ;
      dy_x[xjndx] = dy_sys->vlb[xjndx] ; }
    dy_cbar[xjndx] = cbarj ;
    if (dy_calcprimals() == FALSE)
    { errmsg(316,rtnnme,dy_sys->nme) ;
      return (-1) ; }
    dy_setbasicstatus() ;
/*
  Log the pivot and run it past La Duenna.
*/
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.dual >= 4)
    { deltaj = dirj*(dy_sys->vub[xjndx]-dy_sys->vlb[xjndx]) ;
      dy_logpivot(dyrOK,xjndx,dirj,cbarj,xjndx,dirj,abarij,deltaj) ; }
#   endif

    dy_lp->phase = dyDUAL ;
    duennaresult = dy_duenna(dyrOK,xjndx,xindx,-1,-1) ;
    dy_lp->phase = dyADDVAR ;
    switch (duennaresult)
    { case dyrOK:
      case dyrRESELECT:
      case dyrOPTIMAL:
      { retval = 1 ;
	break ; }
      case dyrPUNT:
      case dyrLOSTDFEAS:
      { retval = 0 ;
	break ; }
      default:
      { retval = -1 ; } }
/*
  If the pivot went through without problem, keep going. The next thing to do
  is call dy_dualout to see what we select for the leaving variable.
  dy_dualout can return one of dyrOK (found a candidate), dyrOPTIMAL or
  dyrPUNT (no candidates, or what there were are are flagged with the NOPIVOT
  qualifier). If candxi isn't the same as xindx, we're done with these wierd
  pivots and can go back to normal dual simplex.

  It's obvious we want to return to dual simplex (retval = 1) when we've
  successfully pivoted and will select a new leaving variable. The reason
  for returning when we get dyrOPTIMAL or dyrPUNT from dy_dualout is that it's
  easier to allow dy_dual to take clear of cleanly winding up the dual simplex
  run.
*/
    if (retval != 1) break ;

    outresult = dy_dualout(&candxi) ;
    if (outresult != dyrOK)
    { retval = 1 ;
      break ; }
    if (xindx != candxi)
    { retval = actcnt ;
      break ; }
/*
  Sigh. We need another type 3 candidate. Open a loop to scan the inactive
  variables and see if we can find another one. By the time we get here, we've
  been through the inactive variables once already, so drop the paranoia.
*/

    oxjndx = 0 ;
    distj = dy_tols->inf ;
    cbarj = 0.0 ;
    abarij = 0.0 ;
    fatal = FALSE ;
    for (oxkndx = 1 ; oxkndx <= orig_sys->varcnt ; oxkndx++)
    { if (dy_origvars[oxkndx] > 0) continue ;
/*
  Status check. We never activate fixed variables.
*/
      xkstatus = (flags) -dy_origvars[oxkndx] ;
      if (flgon(xkstatus,vstatNBFX))
      {
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.varmgmt >= 3)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n      skipping %s %s (%d).",dy_prtvstat(xkstatus),
		      consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ; }
#       endif
	continue ;
      }
/*
  Fetch the column for x<k> and calculate abar<ik> = beta<i>a<k> and
  cbar<k> = c<k> - ya<k>. We use only the active elements of a<k>. If
  abar<ik> = 0, there's no point in activating the variable.
*/
      if (consys_getcol_pk(orig_sys,oxkndx,&ak) == FALSE)
      { errmsg(122,rtnnme,orig_sys->nme,
	       "column",consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
	fatal = TRUE ;
	break ; }
      abarik = 0.0 ;
      cbark = orig_sys->obj[oxkndx] ;
      for (pkndx = 0,aqk = ak->coeffs ; pkndx < ak->cnt ; pkndx++,aqk++)
      { if (dy_origcons[aqk->ndx] > 0)
	{ xqndx = dy_origcons[aqk->ndx] ;
	  abarik += betai[xqndx]*aqk->val ;
	  cbark -= dy_y[xqndx]*aqk->val ; } }
      setcleanzero(abarik,dy_tols->zero) ;
      if (abarik == 0)
      {
  #    ifndef DYLP_NDEBUG
	if (dy_opts->print.varmgmt >= 3)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n      skipping %s %s (%d), abarik = 0.",
		      dy_prtvstat(xkstatus),
		      consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ; }
  #     endif
	continue ; }
      setcleanzero(cbark,dy_tols->zero) ;
      deltak = cbark/abarik ;
      setcleanzero(deltak,dy_tols->zero) ;
/*
  Is this suitable for a type 3 variable? If we get a distance from bound of
  0.0, it won't get any better. Break out of the search loop and do the pivot.
*/
    if (diri*(-deltak) <= 0)
    { type3eval(orig_sys,xindx,diri,oxkndx,xkstatus,abarik,cbark,
	        &oxjndx,&distj,&cbarj,&dirj,&abarij) ;
      if (distj == 0) break ; } } }
/*
  We've dropped out of the main loop. There are lots of reasons, but the only
  one we need to deal with here is oxjndx == 0, which means that we didn't
  find another type 3 variable to pivot in. In that case, we'll revert to the
  primal just to make sure we havn't missed anything.
*/
  if (ak != NULL) pkvec_free(ak) ;
  if (oxjndx == 0) retval = 0 ;

  return (retval) ; }



int dy_dualaddvars (consys_struct *orig_sys)

/*
  dy_dualaddvars is called when the dual simplex reports unbounded (hence
  primal infeasible), and we need to add dual constraints (primal variables)
  to bound the dual. The problem, in dual terms, is that y<i> can increase
  without bound. The offending y<i> is actually identified by placing the
  index of the leaving primal associated with basis position i into
  dy_lp.ubnd.ndx.
  
  dualaddvars can also be called as part of dylp's initialisation sequence,
  adding variables prior to beginning dual simplex iterations. In this case,
  the routine considers only type 1 activations --- it's too early to pivot.

  There are three possibilities:
    1) There are inactive variables which could bound the dual problem and
       can be activated into the nonbasic partition while retaining dual
       feasibility.

       In this case, we'll activate all such variables and return to the
       dual simplex. If there are no such variables, perhaps ...

    2) There are inactive variables which could bound the dual problem but
       would be dual infeasible if activated into the nonbasic partition.

       In this case, we choose one such variable, according to a permuted
       version of the usual dual pivoting rules, activate it, pivot it into
       the basis, and then return to dual simplex. There's a real chance
       that we'll be right back here on the next pivot, but c'est la vie.
       If there are no such variables, then
    
    3) There are inactive variables with a<ij> != 0, but they don't satisfy
       the conditions for 1) or 2). We're at a point that's primal infeasible
       and dual unbounded for the subset we've been working with, and adding
       any inactive variables will make us dual infeasible. What's called for
       here is bound-to-bound (b2b) pivots as we activate variables, chosen
       in such a way that we eventually want to chose some other dual
       variable to enter. I.e., we change y<i>'s reduced cost until some
       other dual looks like a better candidate to enter. In primal terms, we
       want to drive x<i> toward feasibility.

  After some preliminaries, the inactive variables are scanned, looking for
  type 1, 2, and 3 variables. Once we've found the first type 1, we quit
  looking for the other two types. Any type 1 variables are activated as they
  are found.

  If we find no type 1 variables, but have a suitable type 2 variable, we do
  the pivot here and then return to dual simplex. It's a pretty costly way to
  do a dual pivot (in terms of the overhead of leaving and reentering dual
  simplex) but it shouldn't happen too often.

  If we find no type 1 or 2 variables, we do the b2b pivot on the type 3
  variable and then run dualout to see if we get a different entering dual.
  If so, we return to dual simplex. If not, we enter a loop that repeats
  the process. (We need a pivot or a different y<i> before we can hope to
  find any type 1 or type 2 variables.)

  To see how the type 2 pivot works, recall that in dual simplex one looks
  for an entering variable x<j> s.t. j = arg min{k} |cbar<k>/abar<ik>| and
  the signs of cbar<k> and abar<ik> are appropriate. So, for x<i> leaving at
  lb<i> and x<j> entering from lb<i>, we need cbar<j> >= 0 and abar<ij> < 0
  so that cbar<i> = -cbar<j>/abar<ij> >= 0 after the pivot. But if we had
  cbar<j> <= 0 and abar<ij> > 0, cbar<i> would still be ok. But we couldn't
  activate x<j> into the nonbasic partition, because for x<j> NBLB, we need
  cbar<j> >= 0 for dual feasibility. See dy_dualpivot.c:dy_dualin for a more
  complete explanation of the dual pivot rules. We apply them here with the
  bound reversed for the entering variable.

  Parameters:
    orig_sys:	The original constraint system
  
  Returns: number of variables activated; -1 if there's an error
*/

{ int acttype,newcnt,xindx,xipos,diri,evalcode,
      oxkndx,pkndx,xqndx,ox2ndx,dir2,ox3ndx,dir3 ;
  double *betai,abarik,cbark,deltak,delta2,cbar2,dist3,cbar3,abari3 ;
  flags xkstatus ;
  bool fatal ;
  pkvec_struct *ak ;
  pkcoeff_struct *aqk ;
  int retval ;

  const char *rtnnme = "dy_dualaddvars" ;

  retval = -1 ;

# ifdef PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (retval) ; }
  if (dy_lp->simplex.next != dyDUAL)
  { errmsg(4,rtnnme,"phase",dy_prtlpphase(dy_lp->simplex.next,FALSE)) ;
    return (retval) ; }
  if (!dy_lp->lpret == lpINFEAS)
  { errmsg(4,rtnnme,"lp return code",dy_prtlpret(dy_lp->lpret)) ;
    return (retval) ; }
  if (dy_lp->ubnd.ndx == 0 || abs(dy_lp->ubnd.ndx) > dy_sys->varcnt )
  { errmsg(102,rtnnme,dy_sys->nme,"variable",abs(dy_lp->ubnd.ndx),
	   1,dy_sys->varcnt) ;
    return (retval) ; }
# endif
/*
  Preliminaries. First figure out how the entering variable is leaving. Then
  calculate beta<i>.
*/
  xindx = dy_lp->ubnd.ndx ;
  if (xindx < 0)
  { xindx = -xindx ;
    diri = -1 ; }
  else
  { diri = 1 ; }
  xipos = dy_var2basis[xindx] ;
# ifdef PARANOIA
  if (xipos <= 0 || xipos > dy_sys->concnt)
  { errmsg(102,rtnnme,dy_sys->nme,"constraint",xipos,1,dy_sys->concnt) ;
    return (retval) ; }
# endif
  betai = (double *) CALLOC(dy_sys->concnt+1,sizeof(double)) ;
  betai[xipos] = 1.0 ;
  dy_btran(betai) ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.varmgmt >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    leaving variable %s (%d) ",
	        consys_nme(dy_sys,'v',xindx,TRUE,NULL),xindx) ;
    if (diri > 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "rising to lb = %g.",dy_sys->vlb[xindx]) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "falling to ub = %g.",dy_sys->vub[xindx]) ; } }
# endif
/*
  Now open up a loop to walk the variables in orig_sys and check the inactive
  ones for activation. If we're here as part of initialisation, we're only
  interested in type 1 activations.

  (This routine not currently used in initialisation. --lh, 051203 --)
*/
  if (dy_lp->phase == dyINIT)
  { acttype = 1 ; }
  else
  { acttype = dy_opts->dualadd ; }
  ak = NULL ;
  newcnt = 0 ;
  ox2ndx = 0 ;
  delta2 = dy_tols->inf ;
  cbar2 = 0.0 ;
  ox3ndx = 0 ;
  dist3 = -dy_tols->inf ;
  cbar3 = 0.0 ;
  abari3 = 0 ;
  fatal = FALSE ;
  for (oxkndx = 1 ; oxkndx <= orig_sys->varcnt ; oxkndx++)
  { if (dy_origvars[oxkndx] > 0) continue ;
/*
  Status check. We never activate fixed variables.
*/
    xkstatus = (flags) -dy_origvars[oxkndx] ;
#   ifdef PARANOIA
    if (flgoff(xkstatus,vstatNBFX|vstatNBUB|vstatNBLB|vstatNBFR))
    { errmsg(433,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "inactive",consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),
	     oxkndx,dy_prtvstat(xkstatus)) ;
      fatal = TRUE ;
      break ; }
#   endif
    if (flgon(xkstatus,vstatNBFX))
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      skipping %s %s (%d).",dy_prtvstat(xkstatus),
		    consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ; }
#     endif
      continue ;
    }
/*
  Fetch the column for x<k> and calculate abar<ik> = beta<i>a<k> and
  cbar<k> = c<k> - ya<k>. We use only the active elements of a<k>. If
  abar<ik> = 0, there's no point in activating the variable.
*/
    if (consys_getcol_pk(orig_sys,oxkndx,&ak) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,
	     "column",consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
      fatal = TRUE ;
      break ; }
    abarik = 0.0 ;
    cbark = orig_sys->obj[oxkndx] ;
    for (pkndx = 0,aqk = ak->coeffs ; pkndx < ak->cnt ; pkndx++,aqk++)
    { if (dy_origcons[aqk->ndx] > 0)
      { xqndx = dy_origcons[aqk->ndx] ;
	abarik += betai[xqndx]*aqk->val ;
	cbark -= dy_y[xqndx]*aqk->val ; } }
    setcleanzero(abarik,dy_tols->zero) ;
    if (abarik == 0)
    {
#    ifndef DYLP_NDEBUG
      if (dy_opts->print.varmgmt >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n      skipping %s %s (%d), abarik = 0.",
		    dy_prtvstat(xkstatus),
		    consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ; }
#     endif
      continue ; }
    setcleanzero(cbark,dy_tols->zero) ;
    deltak = cbark/abarik ;
    setcleanzero(deltak,dy_tols->zero) ;
/*
  What do we have? If dir<i>*(-cbar<k>/abar<ik>) >= 0, x<k> is a possible
  type 1 or type 2.  Check first for a type 1 variable. type1var will make
  the necessary checks and activate x<j> if it qualifies. If this isn't a
  type 1, and we have none to date, try for a type 2. type2var will make the
  necessary checks, and replace the incumbent type 2 variable if appropriate.
*/
    if (diri*(-deltak) >= 0)
    { evalcode = type1var(orig_sys,xindx,diri,oxkndx,xkstatus,abarik,cbark) ;
      if (evalcode < 0)
      { errmsg(400,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,1,
	       consys_nme(orig_sys,'v',oxkndx,TRUE,NULL),oxkndx) ;
	fatal = TRUE ;
	break ; }
      else
      if (evalcode > 0)
      { newcnt++ ; }
      else
      if (newcnt == 0 && acttype >= 2)
      { type2eval(orig_sys,xindx,diri,oxkndx,xkstatus,abarik,cbark,
		 &ox2ndx,&delta2,&cbar2,&dir2) ; } }
/*
  Maybe this variable is a possible type 3. If we don't have anything better
  going, call type3eval to do the evaluation.
*/
    else
    if (newcnt == 0 && ox2ndx == 0 && acttype >= 3)
    { type3eval(orig_sys,xindx,diri,oxkndx,xkstatus,abarik,cbark,
	        &ox3ndx,&dist3,&cbar3,&dir3,&abari3) ; }
/*
  Are we over our limit? If so, abort the loop.
*/
    if (dy_opts->addvar > 0 && newcnt >= dy_opts->addvar) break ;
  }
/*
  Free the storage we've acquired. Bail out if we've had a fatal error. If
  we've activated variables, return to dual simplex.
*/
  if (ak != NULL) pkvec_free(ak) ;
  if (fatal == TRUE) return (retval) ;
  if (newcnt > 0)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.varmgmt >= 1)
    { if (dy_opts->print.varmgmt >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho," %d activations.",newcnt) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    constraint system %s now %d x %d (%d + %d).",
		  dy_sys->nme,dy_sys->concnt,dy_sys->varcnt,dy_sys->archvcnt,
		  dy_sys->logvcnt) ; }
#   endif
    retval = newcnt ; }
/*
  We didn't activate any type 1 variables. Do we have a type 2 variable which
  we can pivot in?
*/
  else
  if (ox2ndx > 0 && acttype >= 2)
  { retval = type2activate(orig_sys,xindx,diri,ox2ndx,dir2,cbar2) ; }
/*
  Well, do we have a type 3 variable? If so, call type3activate to do the
  bound-to-bound pivot. type3activate will keep on with the bound-to-bound
  pivots until x<i> is no longer selected as the leaving variable or until
  there are no type 3 variables remaining.
*/
  else
  if (ox3ndx > 0 && acttype >= 3)
  { retval = type3activate(orig_sys,
			   betai,xindx,diri,ox3ndx,dir3,cbar3,abari3) ; }
/*
  Nothing! Guess we're done, eh?
*/
  else
  { retval = 0 ; }

  if (betai != NULL) FREE(betai) ;

/*
  A little paranoia and we're out of here.
*/

# ifdef PARANOIA
  if (dy_chkdysys(orig_sys) == FALSE) retval = -1 ;
# endif

  return (retval) ; }

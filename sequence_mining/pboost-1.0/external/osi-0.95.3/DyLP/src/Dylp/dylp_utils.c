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
  This file contains utility routines for the dylp dynamic simplex package.
*/

#define DYLP_INTERNAL

#include "dylib_strrtns.h"
#include "dylp.h"
#include <limits.h>

static char sccsid[] UNUSED = "@(#)dylp_utils.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dylp_utils.c 94 2006-06-29 23:06:51Z lou $" ;



lpret_enum dyret2lpret (dyret_enum dyret)

/*
  Simple-minded utility routine to translate dyret_enum codes into lpret_enum
  codes. A consequence of trying to contain ripples from a rewrite of the 
  simplex routines. dyret codes without direct correspondence translate into
  lpINV, but in context this isn't a problem.

  Parameter:
    dyret:	dyret_enum code

  Returns: corresponding lpret_enum code, as best it can figure it out.
*/

{ switch (dyret)
  { case dyrOK:
    { return (lpINV) ; }
    case dyrOPTIMAL:
    { return (lpOPTIMAL) ; }
    case dyrUNBOUND:
    { return (lpUNBOUNDED) ; }
    case dyrSWING:
    { return (lpSWING) ; }
    case dyrINFEAS:
    { return (lpINFEAS) ; }
    case dyrREQCHK:
    { return (lpINV) ; }
    case dyrACCCHK:
    { return (lpACCCHK) ; }
    case dyrLOSTPFEAS:
    { return (lpLOSTFEAS) ; }
    case dyrLOSTDFEAS:
    { return (lpLOSTFEAS) ; }
    case dyrDEGEN:
    { return (lpINV) ; }
    case dyrRESELECT:
    { return (lpINV) ; }
    case dyrMADPIV:
    { return (lpINV) ; }
    case dyrPUNT:
    { return (lpPUNT) ; }
    case dyrPATCHED:
    { return (lpINV) ; }
    case dyrNUMERIC:
    { return (lpFATAL) ; }
    case dyrBSPACE:
    { return (lpNOSPACE) ; }
    case dyrSTALLED:
    { return (lpSTALLED) ; }
    case dyrITERLIM:
    { return (lpITERLIM) ; }
    case dyrFATAL:
    { return (lpFATAL) ; }
    case dyINV:
    { return (lpINV) ; }
    default:
    { return (lpINV) ; } } }




bool dy_reducerhs (double *rhs, bool init)

/*
  This routine calculates the reduced right-hand-side vector b - Nx<N>, which
  accounts for nonbasic variables which have nonzero values (either because
  they are at a nonzero bound, or are superbasic).  With that in mind, the
  routine needs a valid basis, status, and values for the nonbasic
  variables when it's called.

  Parameters:
    rhs:        the rhs vector to be reduced
    init:	if TRUE, rhs will be initialised to dy_sys->rhs

  Returns: true, basically; a false return indicates serious confusion
*/

{ int vndx,pkndx,cndx ;
  pkvec_struct *pkcol ;
  const char *rtnnme = "dy_reducerhs" ;

#ifdef PARANOIA
  if (rhs == NULL)
  { errmsg(2,rtnnme,"rhs") ;
    return (FALSE) ; }
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (FALSE) ; }
  if (init == TRUE && dy_sys->rhs == NULL)
  { errmsg(2,rtnnme,"dy_sys->rhs") ;
    return (FALSE) ; }
  if (dy_x == NULL)
  { errmsg(2,rtnnme,"dy_x") ;
    return (FALSE) ; }
  if (dy_status == NULL)
  { errmsg(2,rtnnme,"dy_status") ;
    return (FALSE) ; }
#endif

  if (init == TRUE)
    memcpy(rhs,dy_sys->rhs,(dy_sys->concnt+1)*sizeof(double)) ;
/*
  Walk the status array. For each nonbasic variable with a non-zero value,
  retrieve the column and reduce the rhs accordingly. Note that we have to
  consider the logical variables, as a range constraint produces an
  upper-bounded slack. NBFR variables are 0 by definition.
*/
  pkcol = NULL ;
  for (vndx = 1 ; vndx <= dy_sys->varcnt ; vndx++)
    if (flgon(dy_status[vndx],vstatNONBASIC) && dy_x[vndx] != 0.0)
    {
#     ifdef PARANOIA
      if (fabs(dy_x[vndx]) >= dy_tols->inf)
      { errmsg(315,rtnnme,consys_nme(dy_sys,'v',vndx,TRUE,NULL),vndx,
	       dy_prtvstat(dy_status[vndx])) ;
	return (FALSE) ; }
#     endif
      if (consys_getcol_pk(dy_sys,vndx,&pkcol) == FALSE)
      { errmsg(122,rtnnme,dy_sys->nme,"column",
	       consys_nme(dy_sys,'v',vndx,TRUE,NULL),vndx) ;
	return (FALSE) ; }
      for (pkndx = 0 ; pkndx < pkcol->cnt ; pkndx++)
      { cndx = pkcol->coeffs[pkndx].ndx ;
	rhs[cndx] -= pkcol->coeffs[pkndx].val*dy_x[vndx] ; } }
/*
  And a little quick grooming.
*/
  for (cndx = 1 ; cndx <= dy_sys->concnt ; cndx++)
    setcleanzero(rhs[cndx],dy_tols->zero) ;
  
  if (pkcol != NULL) pkvec_free(pkcol) ;

  return (TRUE) ; }



bool dy_calcprimals (void)

/*
  This routine calculates the values of the basic variables, leaving the
  result in dy_x. dy_xbasic is updated if the antidegeneracy machinery is
  idle. There's little to do here - dy_reducerhs and dy_ftran do all the
  heavy lifting.

  Returns: TRUE unless there's serious confusion.
*/

{ int vndx,bndx ;
  double *xvec ;
  bool degenActive ;
  const char *rtnnme = "dy_calcprimals" ;

# ifndef DYLP_NDEBUG
  int print ;

  switch (dy_lp->phase)
  { case dyPRIMAL1:
    { print = dy_opts->print.phase1 ;
      break ; }
    case dyPRIMAL2:
    { print = dy_opts->print.phase2 ;
      break ; }
    case dyDUAL:
    { print = dy_opts->print.dual ;
      break ; }
    case dyADDCON:
    { if (dy_opts->print.conmgmt >= 3)
	print = 5 ;
      else
	print = 0 ;
      break ; }
    case dyADDVAR:
    { if (dy_opts->print.varmgmt >= 3)
	print = 5 ;
      else
	print = 0 ;
      break ; }
    case dyINIT:
    { if (dy_opts->print.crash >= 4)
	print = 7 ;
      else
	print = 0 ;
      break ; }
    default:
    { print = 0 ;
      break ; } }
# endif

#ifdef PARANOIA
  if (dy_xbasic == NULL)
  { errmsg(2,rtnnme,"xbasic") ;
    return (FALSE) ; }
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (FALSE) ; }
  if (dy_x == NULL)
  { errmsg(2,rtnnme,"dy_x") ;
    return (FALSE) ; }
  if (dy_basis == NULL)
  { errmsg(2,rtnnme,"dy_basis") ;
    return (FALSE) ; }
#endif

  xvec = NULL ;
  if (dy_lp->degen > 0 &&
      (dy_lp->phase == dyPRIMAL1 || dy_lp->phase == dyPRIMAL2))
  { degenActive = TRUE ; }
  else
  { degenActive = FALSE ; }

/*
  For the normal case of a non-empty basis ...
*/
  if (dy_sys->concnt > 0)
  {
/*
  Is the antidegeneracy machinery active? If so, we need to allocate a
  vector for the calculation so as not to erase the perturbation in dy_xbasic.
*/
    if (degenActive == TRUE)
    { xvec = (double *) MALLOC((dy_sys->concnt+1)*sizeof(double)) ; }
    else
    { xvec = dy_xbasic ; }
/*
  Calculate x<B> = inv(B)(rhs - A<N>x<N>).  dy_reducerhs accounts for
  nonbasic variables with nonzero values. dy_ftran takes care of
  premultiplying by the basis inverse. Then update the primal feasibility
  tolerance.
*/
    if (dy_reducerhs(xvec,TRUE) == FALSE)
    { errmsg(340,rtnnme,dy_sys->nme) ;
      if (degenActive == TRUE) FREE(xvec) ;
      return (FALSE) ; }
    dy_ftran(xvec,FALSE) ;

    dy_lp->prim.norm1 = exvec_1norm(xvec,dy_sys->concnt) ;
    dy_lp->prim.norm2 = exvec_2norm(xvec,dy_sys->concnt) ;
    dy_lp->prim.max = exvec_infnorm(xvec,dy_sys->concnt,&dy_lp->prim.maxndx) ;
    dy_tols->pfeas = dy_lp->prim.max ;
    if (dy_tols->pfeas < 10.0)
      dy_tols->pfeas = dy_tols->zero ;
    else
      dy_tols->pfeas = log10(dy_tols->pfeas)*dy_tols->zero ;
    dy_tols->pfeas = dy_tols->pfeas_scale*dy_tols->pfeas ;

/*
  Update dy_x.
*/
    if (dy_lp->phase == dyPRIMAL1 || dy_lp->phase == dyPRIMAL2)
    { for (bndx = 1 ; bndx <= dy_sys->concnt ; bndx++)
      { vndx = dy_basis[bndx] ;
	if (flgon(dy_status[vndx],vstatBFX))
	{ 
#         ifdef PARANOIA
	  if (!withintol(xvec[bndx],dy_sys->vub[vndx],dy_tols->zero))
	  { if (!withintol(xvec[bndx],dy_sys->vub[vndx],dy_tols->pfeas*100))
	    { errmsg(333,rtnnme,dy_sys->nme,
		     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		     consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		     dy_prtvstat(dy_status[vndx]),xvec[bndx],
		     dy_x[vndx],"fix",fabs(dy_sys->vub[vndx]-xvec[bndx]),
		     dy_tols->pfeas*100) ; }
	    else
	    { warn(333,rtnnme,dy_sys->nme,
		   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		   consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		   dy_prtvstat(dy_status[vndx]),xvec[bndx],
		   dy_x[vndx],"fix",fabs(dy_sys->vub[vndx]-xvec[bndx]),
		   dy_tols->zero) ; } }
#	  endif
	  dy_x[vndx] = xvec[bndx] ; }
	else
	{ dy_x[vndx] = xvec[bndx] ; } } }
    else
    { for (bndx = 1 ; bndx <= dy_sys->concnt ; bndx++)
      { vndx = dy_basis[bndx] ;
	dy_x[vndx] = xvec[bndx] ; } } }
/*
  And for the pathological case ...  Lie about the norms, given their typical
  usage.
*/
  else
  { dy_tols->pfeas = dy_tols->pfeas_scale*dy_tols->zero ;
    dy_lp->prim.norm1 = 1 ;
    dy_lp->prim.norm2 = 1 ;
    dy_lp->prim.max = 1 ;
    dy_lp->prim.maxndx = -1 ; }

# ifndef DYLP_NDEBUG
  if (print >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\tprim.max = %g, scale = %g, pzero = %g, pfeas = %g.",
	        dy_lp->prim.max,dy_tols->pfeas_scale,
	        dy_tols->zero,dy_tols->pfeas) ; }
# endif

# ifdef PARANOIA
/*
  Check the nonbasic variables to see that their value agrees with their
  status.
*/
  for (vndx = 1 ; vndx <= dy_sys->varcnt ; vndx++)
  { switch (dy_status[vndx])
    { case vstatNBUB:
      case vstatNBFX:
      { if (!atbnd(dy_x[vndx],dy_sys->vub[vndx]))
	{ errmsg(333,rtnnme,dy_sys->nme,
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		 consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		 dy_prtvstat(dy_status[vndx]),dy_x[vndx],
		 dy_sys->vub[vndx],"ub",fabs(dy_sys->vub[vndx]-dy_x[vndx]),
		 dy_tols->pfeas*(1+fabs(dy_sys->vub[vndx]))) ;
	  if (degenActive == TRUE) FREE(xvec) ;
	  return (FALSE) ; }
	break ; }
      case vstatNBLB:
      { if (!atbnd(dy_x[vndx],dy_sys->vlb[vndx]))
	{ errmsg(333,rtnnme,dy_sys->nme,
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		 consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		 dy_prtvstat(dy_status[vndx]),dy_x[vndx],
		 dy_sys->vlb[vndx],"lb",fabs(dy_x[vndx]-dy_sys->vlb[vndx]),
		 dy_tols->pfeas*(1+fabs(dy_sys->vlb[vndx]))) ;
	  if (degenActive == TRUE) FREE(xvec) ;
	  return (FALSE) ; }
	break ; } } }
# endif
# ifndef DYLP_NDEBUG
/*
  Do any requested printing. Level 5 gets the basic variables, level 6 adds
  the nonbasic variables.
*/
  if (print >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: recalculated primal variables:",rtnnme) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%8s%20s%16s%16s%16s%8s","pos'n","var (ndx)",
	        "lb","val","ub","status") ;
    if (degenActive) dyio_outfmt(dy_logchn,dy_gtxecho,"%16s","perturbation") ;
    for (bndx = 1 ; bndx <= dy_sys->concnt ; bndx++)
    { vndx = dy_basis[bndx] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%8d%14s (%3d)%16.8g%16.8g%16.8g%8s",bndx,
		  consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		  dy_sys->vlb[vndx],dy_x[vndx],dy_sys->vub[vndx],
		  dy_prtvstat(dy_status[vndx])) ;
      if (degenActive == TRUE && dy_degenset[bndx] > 0)
      { if (dy_brkout[bndx] > 0)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"%16.8g",
		      dy_xbasic[bndx]-dy_sys->vlb[vndx]) ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,"%16.8g",
		      dy_sys->vub[vndx]-dy_xbasic[bndx]) ; } }
    if (print >= 6)
    { for (vndx = 1 ; vndx <= dy_sys->varcnt ; vndx++)
      { if (dy_var2basis[vndx] != 0) continue ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n%8s%14s (%3d)%16.8g%16.8g%16.8g%8s"," ",
		    consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		    dy_sys->vlb[vndx],dy_x[vndx],dy_sys->vub[vndx],
		    dy_prtvstat(dy_status[vndx])) ; } } }
# endif

/*
  Clean up and return.
*/
  if (degenActive == TRUE) FREE(xvec) ;

  return (TRUE) ; }



dyret_enum dy_updateprimals (int j, double deltaj, double *p_abarj)

/*
  This routine updates the value and status of the primal basic variables
  to reflect the change in x<j>. The calculation is straightforward:
    x<B> = inv(B)b - (inv(B)a<j>)*delta
  It's assumed that deltaj and abarj have the correct signs. If the client
  doesn't supply abarj, it will be calculated here, but not returned to the
  client. This routine is a bit more specialized than dy_calcprimals. It's
  intended to be used in the primal update portion of a pivot, and updates
  the status vector along with the value.

  The routine makes an effort to snap the updated value of x<k> to lb<k>, 0,
  or ub<k> whenever possible. It calculates separate tolerances for each
  value based on the magnitude of the value and the distance travelled
  (deltak).

  If deltaj == 0, the routine will quietly return.

  We have to be careful to preserve the pivot reject flag w.r.t. variables on
  the pivot reject list.

  Parameters:
    j:		index of entering variable x<j>
    deltaj:	amount of change in x<j>
    p_abarj:	inv(B)a<j>; if NULL, will be calculated here

  Returns: dyrOK if the update is successful
	   dyrREQCHK if a bogus value is created
	   dyrSWING if we detected a massive move into infeasibility
	   dyrFATAL is a positive indication of failure
*/

{ int kpos,k,m ;
  double *abarj,abarkj,deltak,xk,ubk,lbk,newxk ;
  flags statk,newstatk,qualsk ;

  double eps0,epsu,epsl ;

  int swingndx ;
  double swingratio,maxswing ;
  bool swing ;

  dyret_enum retval ;

  const char *rtnnme = "dy_updateprimals" ;

/*
  Setup
*/
  if (withintol(deltaj,0,dy_tols->zero)) return (dyrOK) ;
  retval = dyrINV ;
  m = dy_sys->concnt ;
  swing = FALSE ;
  maxswing = 0 ;
  swingndx = -1 ;
/*
  Did the user provide the ftran'd column? If not, retrieve a<j> and
  calculate abar<j>.
*/
  if (p_abarj != NULL)
  { abarj = p_abarj ; }
  else
  { abarj = NULL ;
    if (consys_getcol_ex(dy_sys,j,&abarj) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,
	     "column",consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
      if (abarj != NULL) FREE(abarj) ;
      return (dyrFATAL) ; }
    dy_ftran(abarj,FALSE) ; }
/*
  Now walk the basis. For each abar<kj> != 0, update x<k>.
*/
  for (kpos = 1 ; kpos <= m ; kpos++)
  { abarkj = abarj[kpos] ;
    if (abarkj == 0) continue ;
    deltak = abarkj*deltaj ;
    if (withintol(deltak,0,dy_tols->zero)) continue ;
    k = dy_basis[kpos] ;
    xk = dy_xbasic[kpos] ;
    statk = getflg(dy_status[k],vstatSTATUS) ;
    qualsk = getflg(dy_status[k],vstatQUALS) ;
    eps0 = dy_tols->zero ;
    ubk = dy_sys->vub[k] ;
    if (ubk < dy_tols->inf)
    { epsu = dy_tols->pfeas*(1.0+fabs(ubk)) ; }
    else
    { epsu = 0 ; }
    lbk = dy_sys->vlb[k] ;
    if (-dy_tols->inf < lbk)
    { epsl = dy_tols->pfeas*(1.0+fabs(lbk)) ; }
    else
    { epsl = 0 ; }
    newxk = xk-deltak ;
    setcleanzero(newxk,eps0) ;
/*
  Choose a new status for x<k>. The tests for abovebnd and belowbnd will use
  the feasibility tolerance.
*/
    if (statk == vstatBFR)
    { newstatk = vstatBFR ; }
    else
    if (belowbnd(newxk,ubk))
    { if (abovebnd(newxk,lbk))
      { newstatk = vstatB ; }
      else
      if (belowbnd(newxk,lbk))
      { newstatk = vstatBLLB ; }
      else
      { newstatk = vstatBLB ; } }
    else
    { if (abovebnd(newxk,ubk))
      { newstatk = vstatBUUB ; }
      else
      { newstatk = vstatBUB ; } }

    if (flgon(newstatk,vstatBLB|vstatBUB) && lbk == ubk)
    { newstatk = vstatBFX ; }
/*
  Check for bogus values?

  Clean zeros and values judged to be at bound cannot be bogus, so we can
  skip the test.  If bogus === 1.0, there's no point in doing these tests,
  because the bogus region is identical to the tolerance region.  The
  corrective action for a bogus value is to refactor. If basis.etas == 1,
  we've just refactored, so again, no sense in checking.
*/
    if (dy_lp->basis.etas > 1 && dy_tols->bogus > 1.0 &&
	newxk != 0 && flgoff(newstatk,vstatBLB|vstatBUB|vstatBFX))
    { if (fabs(newxk) < eps0*dy_tols->bogus)
      { retval = dyrREQCHK ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.pivoting >= 1)
	  warn(374,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,"x",k,fabs(newxk),eps0*dy_tols->bogus,
	       eps0*dy_tols->bogus-newxk) ;
#       endif
      }
      if (flgon(newstatk,vstatB|vstatBLLB) &&
	  fabs(lbk-newxk) < epsl*dy_tols->bogus)
      { retval = dyrREQCHK ;
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.pivoting >= 1)
	  warn(375,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,"lb",k,"x",k,
	       lbk,newxk,lbk-newxk,epsl*dy_tols->bogus,
	       epsl*dy_tols->bogus-(lbk-newxk)) ;
#	endif
      }
      if (flgon(newstatk,vstatB|vstatBUB) &&
	  fabs(newxk-ubk) < epsu*dy_tols->bogus)
      { retval = dyrREQCHK ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.pivoting >= 1)
	{ warn(375,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,"ub",k,"x",k,
	       ubk,newxk,ubk-newxk,epsu*dy_tols->bogus,
	       epsu*dy_tols->bogus-(ubk-newxk)) ; }
#       endif
      } }
/*
  Are we headed off into the ozone? Make a note.
*/
    if (flgon(newstatk,vstatBLLB|vstatBUUB))
    { swingratio = (fabs(newxk)+1)/(fabs(dy_xbasic[kpos])+1) ;
      if (swingratio > dy_tols->swing)
      { swing = TRUE ;
	if (swingratio > maxswing)
	{ maxswing = swingratio ;
	  swingndx = k ; } } }
/*
  Set the new status and value into the status and value arrays.
*/
    dy_xbasic[kpos] = newxk ;
    dy_x[k] = newxk ;
    dy_status[k] = newstatk|qualsk ; }
/*
  All basic variables are updated. Did we detect a swing into the unknown?
*/
  if (swing == TRUE)
  { retval = dyrSWING ;
    dy_lp->ubnd.ndx = swingndx ;
    dy_lp->ubnd.ratio = maxswing ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.dual >= 2 || dy_opts->print.phase1 >= 2 ||
	dy_opts->print.phase2 >= 2)
    { k = dy_lp->ubnd.ndx ;
      statk = dy_status[k] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    Pseudo-unbounded: growth %e for %s (%d) %s = %g",
		  dy_lp->ubnd.ratio,consys_nme(dy_sys,'v',k,FALSE,NULL),k,
		  dy_prtvstat(statk),dy_x[k]) ;
      if (flgon(statk,vstatBUUB))
      { dyio_outfmt(dy_logchn,dy_gtxecho," > %g.",dy_sys->vub[k]) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho," < %g.",dy_sys->vlb[k]) ; } }
#   endif
  }
# ifndef DYLP_NDEBUG
/*
  Do we need to print the updates? If so, do it here.
*/
  if (dy_opts->print.dual >= 6)
  { bool first = TRUE ;
    for (kpos = 1 ; kpos <= dy_sys->concnt ; kpos++)
      if (abarj[kpos] != 0)
      { if (first == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n\trevised primal variables:") ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n%8s%20s%16s%16s%16s %s","pos'n","var (ndx)",
		      "lb","val","ub","status") ;
	  first = FALSE ; }
	k = dy_basis[kpos] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n%8d%14s (%3d)%16.8g%16.8g%16.8g %s",kpos,
		    consys_nme(dy_sys,'v',k,FALSE,NULL),k,
		    dy_sys->vlb[k],dy_xbasic[kpos],dy_sys->vub[k],
		    dy_prtvstat(dy_status[k])) ; }
    if (first == TRUE)
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tno change to primal variables.") ; }
# endif
/*
  That's it. If retval is still dyrINV, nothing out of the ordinary happened,
  so set it to dyrOK. Clean up and we can return.
*/
  if (retval == dyrINV) retval = dyrOK ;
  if (p_abarj == NULL) FREE(abarj) ;

  return (retval) ; }



void dy_calcduals (void)

/*
  This routine calculates the values of the dual variables as y = c<B>inv(B).
  If antidegeneracy is active, updates are restricted to duals that correspond
  to logicals included in the restricted subproblem.

  The routine also updates dfeas, the scaled zero tolerance for reduced costs.
*/

{ int xkpos,xkndx ;
  bool degenActive ;
  double *yvec ;

# ifndef DYLP_NDEBUG
  int print ;

  switch (dy_lp->phase)
  { case dyPRIMAL1:
    { print = dy_opts->print.phase1 ;
      break ; }
    case dyPRIMAL2:
    { print = dy_opts->print.phase2 ;
      break ; }
    case dyDUAL:
    { print = dy_opts->print.dual ;
      break ; }
    case dyADDCON:
    { if (dy_opts->print.conmgmt >= 3)
	print = 7 ;
      else
	print = 0 ;
      break ; }
    case dyINIT:
    { if (dy_opts->print.crash >= 4)
	print = 7 ;
      else
	print = 0 ;
      break ; }
    default:
    { print = 0 ;
      break ; } }
# endif

  if (dy_lp->phase == dyDUAL && dy_lp->degen > 0)
  { degenActive = TRUE ; }
  else
  { degenActive = FALSE ; }
/*
  For a normal, non-empty basis:
*/
  if (dy_sys->concnt > 0)
  { if (degenActive == TRUE)
    { yvec = (double *) CALLOC((dy_sys->concnt+1),sizeof(double)) ; }
    else
    { yvec = dy_y ; }

    for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
    { xkndx = dy_basis[xkpos] ;
      yvec[xkpos] = dy_sys->obj[xkndx] ; }
    dy_btran(yvec) ;

    dy_lp->dual.norm1 = exvec_1norm(yvec,dy_sys->concnt) ;
    dy_lp->dual.norm2 = exvec_2norm(yvec,dy_sys->concnt) ;
    dy_lp->dual.max = exvec_infnorm(yvec,dy_sys->concnt,&dy_lp->dual.maxndx) ;

    dy_tols->dfeas = dy_lp->dual.max ;
    if (dy_tols->dfeas < 10.0)
      dy_tols->dfeas = dy_tols->cost ;
    else
      dy_tols->dfeas = log10(dy_tols->dfeas)*dy_tols->cost ;
    dy_tols->dfeas = dy_tols->dfeas_scale*dy_tols->dfeas ;
    
    if (degenActive == TRUE)
    { for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      { if (dy_ddegenset[xkpos] == 0)
	  dy_y[xkpos] = yvec[xkpos] ; }
      FREE(yvec) ; } }
/*
  And the pathological case of an empty basis. We lie about the norms, given
  their typical use.
*/
  else
  { dy_tols->dfeas = dy_tols->dfeas_scale*dy_tols->cost ;
    dy_lp->dual.norm1 = 1 ;
    dy_lp->dual.norm2 = 1 ;
    dy_lp->dual.max = 1 ; }

# ifndef DYLP_NDEBUG
  if (print >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\tdual.max = %g, scale = %g, dzero = %g, dfeas = %g.",
	        dy_lp->dual.max,dy_tols->dfeas_scale,
		dy_tols->cost,dy_tols->dfeas) ;
    if (print >= 7)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\n%8s%20s%16s","pos'n","constraint","val") ;
      if (degenActive)
	dyio_outfmt(dy_logchn,dy_gtxecho,"%16s","perturbation") ;
      for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%8d%20s%16.8g",xkpos,
		    consys_nme(dy_sys,'c',xkpos,FALSE,NULL),dy_y[xkpos]) ;
	if (degenActive && dy_ddegenset[xkpos] > 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"%16.8g",dy_y[xkpos]) ; } }
      dyio_outchr(dy_logchn,dy_gtxecho,'\n') ; } }
# endif

  return ; }



double dy_calcobj (void)

/*
  This routine calculates the value of the objective function using the
  primal variables. During primal phase I, suppress inactzcorr, as c<j> = 0
  for nonbasic variables by definition.

  Returns: objective value, or NaN in the event of a paranoid error.
*/

{ int vndx ;
  double z ;

#ifdef PARANOIA
  const char *rtnnme = "dy_calcobj" ;

  if (dy_x == NULL)
  { errmsg(2,rtnnme,"dy_x") ;
    return (quiet_nan(0)) ; }
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (quiet_nan(0)) ; }
  if (dy_sys->obj == NULL)
  { errmsg(2,rtnnme,"dy_sys->obj") ;
    return (quiet_nan(0)) ; }
#endif

  z = 0.0 ;
  for (vndx = 1 ; vndx <= dy_sys->varcnt ; vndx++)
    if (dy_x[vndx] != 0.0) z += dy_sys->obj[vndx]*dy_x[vndx] ;
  if (dy_lp->p1obj.installed == FALSE) z += dy_lp->inactzcorr ;
  setcleanzero(z,dy_tols->zero) ;

  return (z) ; }


double dy_calcdualobj (void)

/*
  This routine calculates the dual objective z = yb, taking into account
  nonbasic variables at nonzero bounds. NOTE that the value returned will not
  reflect the contribution due to inactive variables. Also, unlike the primal
  case, we need to correct dy_y on the fly to compensate for duals involved in
  a perturbed subproblem when antidegeneracy is active.

  Parameters: none

  Returns: dual objective value, or NaN if there's an error during the
	   calculation.
*/

{ int i,j,m,n ;
  flags statj ;

  double z,bndj,cbarj ;
  double *rhs,*vlb,*vub ;

  const char *rtnnme = "dy_calcdualobj" ;

#ifdef PARANOIA
  if (dy_x == NULL)
  { errmsg(2,rtnnme,"dy_x") ;
    return (quiet_nan(0)) ; }
  if (dy_y == NULL)
  { errmsg(2,rtnnme,"dy_y") ;
    return (quiet_nan(0)) ; }
  if (dy_var2basis == NULL)
  { errmsg(2,rtnnme,"dy_var2basis") ;
    return (quiet_nan(0)) ; }
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (quiet_nan(0)) ; }
  if (dy_sys->obj == NULL)
  { errmsg(2,rtnnme,"dy_sys->obj") ;
    return (quiet_nan(0)) ; }
  if (dy_sys->rhs == NULL)
  { errmsg(2,rtnnme,"dy_sys->rhs") ;
    return (quiet_nan(0)) ; }
#endif

  z = 0.0 ;

  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
  rhs = dy_sys->rhs ;
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;

/*
  Start with the easy part: dot(y,b) for active explicit constraints. We need
  to exclude duals incorporated in a degenerate subproblem, which are by
  definition zero.
*/
  for (i = 1 ; i <= m ; i++)
  { j = dy_basis[i] ;
    if (dy_ddegenset[j] == 0 && dy_y[i] != 0)
    { z += dy_y[i]*rhs[i] ; } }
/*
  Now the harder part: the contribution due to variables nonbasic at nonzero
  bounds. For the standard max primal/min dual setup, the math tells us that
  the reduced cost cbar<j> of a variable x<j> that is NBUB is in fact the
  negative of the dual that would be basic if the bound constraint on x<j>
  were explicit in the constraint system. For x<j> NBLB, we're really dealing
  with a >= constraint, and a sign flip is needed to get the value of the
  dual that would be associated with the explicit constraint of -x<j> <=
  -l<j>.  But we'll acheive that by using l<j> instead of -l<j>. Finally,
  because dylp actually uses a min primal/min dual setup, we need one more
  sign flip.

  We need to scan the logicals because range constraints can result in a
  slack with a finite nonzero upper bound.
*/
  for (j = 1 ; j <= n ; j++)
  { if (dy_ddegenset[j] > 0) continue ;
    statj = getflg(dy_status[j],vstatSTATUS) ;
    if (flgon(statj,vstatBASIC)) continue ;
    cbarj = dy_cbar[j] ;
    if (cbarj == 0) continue ;
    switch (statj)
    { case vstatNBLB:
      case vstatNBFX:
      { bndj = vlb[j] ;
	break ; }
      case vstatNBUB:
      { bndj = vub[j] ;
	break ; }
      case vstatNBFR:
      case vstatSB:
      { bndj = 0 ;
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (quiet_nan(0)) ; } }
    if (bndj != 0) z += cbarj*bndj ; }
  if (dy_lp->p1obj.installed == FALSE) z += dy_lp->inactzcorr ;
/*
  Final cleanup and return.
*/
  setcleanzero(z,dy_tols->zero) ;

  return (z) ; }


double dy_calcpinfeas (void)

/*
  This routine calculates primal infeasibility from scratch. It scans the
  status vector and totals infeasibility for all variables flagged as BLLB or
  BUUB.

  Parameters: none

  Returns: total primal infeasibility.
*/

{ int k,n ;
  double infeas,xk,ubk,lbk ;
  flags statk ;
  double *vub,*vlb ;

  n = dy_sys->varcnt ;
  vub = dy_sys->vub ;
  vlb = dy_sys->vlb ;

  infeas = 0 ;
  for (k = 1 ; k <= n ; k++)
  { statk = dy_status[k] ;
    if (flgon(statk,vstatBLLB))
    { lbk = vlb[k] ;
      xk = dy_x[k] ;
      infeas += lbk-xk ; }
    else
    if (flgon(statk,vstatBUUB))
    { ubk = vub[k] ;
      xk = dy_x[k] ;
      infeas += xk-ubk ; } }

  return (infeas) ; }





bool dy_calccbar (void)

/*
  This routine calculates cbar from scratch for all nonbasic variables, as
    cbar<j> = c<j> - dot(y,a<j>).

  The basic restriction for calculating reduced costs is that we must have
  dual variables, and we need to be able to distinguish basic variables from
  nonbasic variables. dy_var2basis is the earliest indicator we can use for
  this --- it's valid early in dyINIT, so that we can use it to try and
  establish a dual feasible start. It remains valid in later phases.

  On the other hand, if we're feeling paranoid, this is a good opportunity to
  check the status of all variables. Hence the two separate code blocks below.

  Parameters: none

  Returns: TRUE if the values are calculated without error, FALSE otherwise.
*/

{ int xjndx ;
  double cbarj ;

# ifdef PARANOIA
  flags xjstatus ;
  const char *rtnnme = "dy_calccbar" ;
# endif

/*
  Open a loop and walk the columns, looking for nonbasic variables.
*/
  for (xjndx = 1 ; xjndx <= dy_sys->varcnt ; xjndx++)
  {
#   ifndef PARANOIA
/* The expedient decision: basic or not? */
    if (dy_var2basis[xjndx] > 0)
    { dy_cbar[xjndx] = 0 ;
      continue ; }
#   else
/* The cautious decision (when available). */
    if (dy_lp->phase == dyINIT)
    { if (dy_var2basis[xjndx] > 0)
      { dy_cbar[xjndx] = 0 ;
	continue ; } }
    else
    { xjstatus = dy_status[xjndx] ;
      if (dy_chkstatus(xjndx) == FALSE) return (FALSE) ;
      if (flgon(xjstatus,vstatBASIC))
      { dy_cbar[xjndx] = 0.0 ;
	continue ; } }
#   endif
/*
  Calculate the reduced cost and store it in cbar<j>.
*/
    cbarj = consys_dotcol(dy_sys,xjndx,dy_y) ;
#   ifdef PARANOIA
    if (isnan(cbarj) == TRUE)
    { errmsg(320,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"y",xjndx,"reduced cost") ;
      return (FALSE) ; }
#   endif
    dy_cbar[xjndx] = dy_sys->obj[xjndx]-cbarj ;
    setcleanzero(dy_cbar[xjndx],dy_tols->cost) ; }
  
  return (TRUE) ; }



void dy_setbasicstatus (void)

/*
  This routine runs through the basis and resets the status of basic variables.
  It's handy after arbitrary changes to the nonbasic variables (as when we
  patch the basis, or flip variables for activation during dual simplex).

  Parameters: none

  Returns: undefined
*/

{ int xindx,bpos ;
  double xi,ubi,lbi ;
  flags xistatus ;

/*
  Just scan the basis, checking the values against the bounds.
*/
  for (bpos = 1 ; bpos <= dy_sys->concnt ; bpos++)
  { xindx = dy_basis[bpos] ;
    xi = dy_x[xindx] ;
    xistatus = dy_status[xindx] ;
    lbi = dy_sys->vlb[xindx] ;
    ubi = dy_sys->vub[xindx] ;
    if (flgoff(xistatus,vstatBFR))
    { if (belowbnd(xi,lbi))
      { dy_status[xindx] = vstatBLLB ; }
      else
      if (atbnd(xi,lbi))
      { if (lbi == ubi)
	  dy_status[xindx] = vstatBFX ;
	else
	  dy_status[xindx] = vstatBLB ;
	dy_x[xindx] = lbi ;
	dy_xbasic[bpos] = lbi ; }
      else
      if (belowbnd(xi,ubi))
      { dy_status[xindx] = vstatB ; }
      else
      if (atbnd(xi,ubi))
      { dy_status[xindx] = vstatBUB ;
	dy_x[xindx] = ubi ;
	dy_xbasic[bpos] = ubi ; }
      else
      { dy_status[xindx] = vstatBUUB ; } }
#   ifndef DYLP_NDEBUG
    if (xistatus != dy_status[xindx])
    { if ((dy_lp->phase != dyADDVAR && dy_opts->print.basis >= 3) ||
	  (dy_lp->phase == dyADDVAR && dy_opts->print.varmgmt >= 3))
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) = %g, status %s ",
		    consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,xi,
		    dy_prtvstat(xistatus)) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"corrected to %s.",
		    dy_prtvstat(dy_status[xindx])) ; } }
#   endif
  }

  return ; }



void dy_dseinit (void)

/*
  This routine calculates from scratch the basis inverse row norms used in
  DSE pricing. The algorithm is straightforward --- use dy_btran to
  calculate e<k>inv(B) = beta<k>, then calculate ||beta<k>||^2.

  Parameters: none

  Returns: undefined
*/

{ int xkpos ;
  double *betak ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.dual >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n   (%s)%d: initialising ||beta<k>||^2 for DSE.",
	        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
# endif

/*
  Grab space for a working vector.
*/
  betak = (double *) MALLOC((dy_sys->concnt+1)*sizeof(double)) ;
/*
  And now the loop to do the calculation. The drill is
    * clear the vector
    * insert a 1 in the proper place
    * btran it
    * calculate the norm.
*/
  for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
  { memset(betak,0,(dy_sys->concnt+1)*sizeof(double)) ;
    betak[xkpos] = 1.0 ;
    dy_btran(betak) ;
    dy_rho[xkpos] = exvec_ssq(betak,dy_sys->concnt) ; }
/*
  Clean up and return.
*/
  FREE(betak) ;

  return ; }



void dy_pseinit (void)

/*
  This routine resets the reference frame for PSE pricing. It's a trivial
  routine, but it seemed like a good idea to encapsulate this bit of code.
  It gets used in several places.

  Parameters: none

  Returns: undefined
*/

{ int xindx ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.phase1 >= 2 || dy_opts->print.phase2 >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n   (%s)%d: initialising ||abar~<k>||^2 for PSE.",
	        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
# endif

  memset(dy_frame,0,(dy_sys->varcnt+1)*sizeof(bool)) ;
  memset(dy_gamma,0,(dy_sys->varcnt+1)*sizeof(double)) ;

  for (xindx = 1 ; xindx <= dy_sys->varcnt ; xindx++)
  { if (flgon(dy_status[xindx],vstatNONBASIC|vstatNBFR))
    { dy_frame[xindx] = TRUE ;
      dy_gamma[xindx] = 1.0 ; } }
  
  return ; }



#ifdef PARANOIA

bool dy_chkstatus (int vndx)

/*
  This routine is a paranoid check routine -- it does a number of consistency
  checks based on the status of the variable. Major basis errors (errors of
  logic and big-time numerical inaccuracy) rate a return value of FALSE.

  There is a relatively common case where accumulated numerical inaccuracy
  results in errors on the order of 1.0e-10 or so; these get a return of
  TRUE; the normal accuracy control machinery should catch these eventually.

  If vndx is 0, all variables are checked.

  Parameters:  vndx:       variable index; if 0, all variables are checked

  Returns: TRUE if status is acceptable, FALSE if major errors turn up.
*/

{ int first,last,xkndx,iter ;
  flags xkstatus ;
  double xk,xbk,ubk,lbk,normxb,tolu,bogusu,toll,bogusl ;
  bool retval,primDegen,dualDegen ;
  char *statk ;
  const char *nmek,*phase ;
  const char *rtnnme = "dy_chkstatus" ;

  if (vndx < 0 || vndx > dy_sys->varcnt)
  { errmsg(102,rtnnme,dy_sys->nme,"variable",vndx,0,dy_sys->varcnt) ;
    return (FALSE) ; }

  if (vndx == 0)
  { first = 1 ;
    last = dy_sys->varcnt ; }
  else
  { first = vndx ;
    last = vndx ; }

  retval = TRUE ;

  if (dy_lp->degen > 0  && dy_lp->phase == dyDUAL)
  { dualDegen = TRUE ; }
  else
  { dualDegen = FALSE ; }
  if (dy_lp->degen > 0 &&
      (dy_lp->phase == dyPRIMAL1 || dy_lp->phase == dyPRIMAL2))
  { primDegen = TRUE ; }
  else
  { primDegen = FALSE ; }

  phase = dy_prtlpphase(dy_lp->phase,TRUE) ;
  iter = dy_lp->tot.iters ;
  normxb = exvec_2norm(dy_xbasic,dy_sys->concnt) ;

  for (xkndx = first ; xkndx <= last ; xkndx++)
  { xkstatus = dy_status[xkndx] ;
    nmek = consys_nme(dy_sys,'v',xkndx,TRUE,NULL) ;
/*
  Clear any qualifiers before we go further.
*/
    clrflg(xkstatus,vstatQUALS) ;
/*
  Check that the status is a single status, and not some mix. If we fail this,
  skip to the next iteration of the loop.
*/
    if (!(xkstatus == vstatBFX || xkstatus == vstatBUB ||
	  xkstatus == vstatB   || xkstatus == vstatBLB ||
	  xkstatus == vstatBLLB || xkstatus == vstatBUUB ||
	  xkstatus == vstatBFR ||
	  xkstatus == vstatNBFX || xkstatus == vstatNBUB ||
	  xkstatus == vstatNBLB || xkstatus == vstatNBFR ||
	  xkstatus == vstatSB))
    { errmsg(300,rtnnme,(int) xkstatus,nmek,xkndx) ;
      retval = FALSE ;
      continue ; }
    statk = dy_prtvstat(xkstatus) ;
/*
  For basic variables, check that basis and var2basis agree. For nonbasic
  and superbasic variables, check that var2basis is 0. Make sure that the
  value in dy_var2basis is valid first.

  While we're here, might as well check that dy_xbasic and dy_x agree for
  basic variables, unless the anti-degeneracy algorithm is working. (Note
  that this > should be < a logical question, not a tolerance question. But
  somewhere, small discrepancies creep in, and I haven't been able to track
  down the cause just yet.)
*/
    xk = dy_x[xkndx] ;
    if (dy_var2basis[xkndx] < 0 || dy_var2basis[xkndx] > dy_sys->varcnt)
    { errmsg(330,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
	     dy_var2basis[xkndx],-1) ;
      retval = FALSE ; }
    else
    if (flgon(xkstatus,vstatBASIC))
    { if (dy_basis[dy_var2basis[xkndx]] != xkndx)
      { errmsg(330,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
	       dy_var2basis[xkndx],dy_basis[dy_var2basis[xkndx]]) ;
	retval = FALSE ; }
      xk = dy_x[xkndx] ;
      xbk = dy_xbasic[dy_var2basis[xkndx]] ;
      if (primDegen == FALSE)
      { if (!withintol(xbk,xk,dy_tols->zero*(1.0+fabs(xbk))))
	{ if (withintol(xbk,xk,10*dy_tols->zero*(1.0+fabs(xbk))))
	  { warn(332,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,xk,xbk,
		   fabs(xk-xbk),dy_tols->zero*(1.0+fabs(xbk))) ; }
	  else
	  { errmsg(332,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,xk,xbk,
		   fabs(xk-xbk),dy_tols->zero*(1.0+fabs(xbk))) ;
	    retval = FALSE ; } }
#	ifdef UNDEFINED
/*
  This bit of code will, on a fairly regular basis, spit out warnings about
  |xbk-xk| around 1.0e-15 -- 1.0e-17. Typically in primal phase II. Logically,
  dy_x and dy_xbasic should be equivalent, and a quick check seems to say that
  they're updated simultaneously. 10^-17 is roughly 2^-56, so this could just
  be noise in the least significant bits.
*/
	else
	if (xbk != xk)
	{ warn(332,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,xk,xbk,
	       fabs(xk-xbk),0.0) ; }
#	endif
      } }
    else
    { if (dy_var2basis[xkndx] != 0)
      { errmsg(330,rtnnme,dy_sys->nme,phase, iter,nmek,xkndx,statk,
	       dy_var2basis[xkndx],dy_basis[dy_var2basis[xkndx]]) ;
	retval = FALSE ; }
      xk = dy_x[xkndx] ; }
/*
  Check that the bounds are consistent with the status.
    * Fixed variables (BFX, NBFX) should have equal bounds. Failure of this
      condition is a logical error. Conversely, variables with equal bounds
      should have status BFX, NBFX, BLLB, or BUUB. The latter two aren't
      legal in P2, but the normal error recovery machinery should cope with
      loss of feasibility as long as the status and value agree.
    * Nonbasic variables should be exactly equal to their bounds --- they
      are supposed to be set to the bound as they're pivoted out. Failure of
      this condition is a logical error.
    * Variables which are basic at bound (BFX, BLB, BUB) should have a
      value which is equal (within tolerance) of the appropriate bound.
    * Strictly basic (B) and superbasic (SB) variables should be strictly
      within bound. For basic variables, if anti-degeneracy is active this
      test must be made against dy_xbasic. Since we've tested above for dy_x
      = dy_xbasic when anti-degeneracy isn't active, we can test here using
      xbasic. Superbasics (by virtue of not being basic) must be tested
      using xk.
    * Free variables should have no bounds.
*/
    ubk = dy_sys->vub[xkndx] ;
    lbk = dy_sys->vlb[xkndx] ;
    if (ubk == lbk &&
	!(flgon(xkstatus,vstatBFX|vstatNBFX) ||
	  (flgon(xkstatus,vstatBLLB|vstatBUUB) && dy_lp->phase != dyPRIMAL2)))
    { errmsg(379,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,ubk,statk) ;
      retval = FALSE ; }
    tolu = dy_tols->pfeas*(1.0+maxx(normxb,fabs(ubk))) ;
    bogusu = 1000*tolu ;
    toll = dy_tols->pfeas*(1.0+maxx(normxb,fabs(lbk))) ;
    bogusl = 1000*toll ;
    switch (xkstatus)
    { case vstatBFX:
      { if (lbk != ubk)
	{ errmsg(331,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,xbk,statk,
		 lbk,ubk,ubk-lbk,0.0) ;
	  retval = FALSE ;
	  break ; }
	if (!atbnd(xbk,ubk))
	{ if (!withintol(xbk,ubk,bogusu))
	  { errmsg(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   xbk,ubk,"lb = ub",fabs(xbk-ubk),bogusu) ;
	    retval = FALSE ; }
	  else
	  { warn(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 xbk,ubk,"lb = ub",fabs(xbk-ubk),tolu) ; } }
	break ; }
      case vstatBUUB:
      { if (!abovebnd(xbk,ubk))
	{ if (xbk < ubk-bogusu)
	  { errmsg(335,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   lbk,xbk,ubk,ubk-xbk,bogusu+tolu) ;
	    retval = FALSE ; }
	  else
	  { warn(335,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 lbk,xbk,ubk,ubk-xbk,tolu) ; } }
	else
	if (dy_lp->phase == dyPRIMAL2)
	  warn(323,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
	       lbk,xbk,ubk,xbk-ubk,tolu) ;
	break ; }
      case vstatBUB:
      { if (!atbnd(xbk,ubk))
	{ if (!withintol(xbk,ubk,bogusu))
	  { errmsg(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   xbk,ubk,"ub",ubk-xbk,bogusu) ;
	    retval = FALSE ; }
	  else
	  { warn(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 xbk,ubk,"ub",ubk-xbk,tolu) ; } }
	break ; }
      case vstatB:
      { if (!(abovebnd(xbk,lbk) && belowbnd(xbk,ubk)))
	{ if (xbk < lbk-bogusl)
	  { errmsg(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   lbk,xbk,ubk,lbk-xbk,bogusl+toll) ;
	    retval = FALSE ; }
	  else
	  if (xbk > ubk+bogusu)
	  { errmsg(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   lbk,xbk,ubk,xbk-ubk,bogusu+tolu) ;
	    retval = FALSE ; }
	  else
	  if (xbk < lbk+toll)
	  { warn(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 lbk,xbk,ubk,lbk-xbk,toll) ; }
	  else
	  if (xbk > ubk-tolu)
	  { warn(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 lbk,xbk,ubk,xbk-ubk,tolu) ; } }
	break ; }
      case vstatSB:
      { if (!(abovebnd(xk,lbk) && belowbnd(xk,ubk)))
	{ if (xk < lbk-bogusl)
	  { errmsg(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   lbk,xk,ubk,lbk-xk,bogusl+toll) ;
	    retval = FALSE ; }
	  else
	  if (xk > ubk+bogusu)
	  { errmsg(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   lbk,xk,ubk,xk-ubk,bogusu+tolu) ;
	    retval = FALSE ; }
	  else
	  if (xk < lbk+toll)
	  { warn(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 lbk,xk,ubk,lbk-xk,toll) ; }
	  else
	  if (xk > ubk-tolu)
	  { warn(322,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 lbk,xk,ubk,xbk-ubk,tolu) ; } }
	break ; }
      case vstatBLB:
      { if (!atbnd(xbk,lbk))
	{ if (!withintol(xbk,lbk,bogusl))
	  { errmsg(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   xbk,lbk,"lb",xbk-lbk,bogusl) ;
	    retval = FALSE ; }
	  else
	  { warn(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 xbk,lbk,"lb",xbk-lbk,toll) ; } }
	break ; }
      case vstatBLLB:
      { if (!belowbnd(xbk,lbk))
	{ if (xbk > lbk+bogusl)
	  { errmsg(335,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		   lbk,xbk,ubk,xbk-lbk,bogusl+toll) ;
	    retval = FALSE ; }
	  else
	  { warn(335,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 lbk,xbk,ubk,xbk-lbk,toll) ; } }
	else
	if (dy_lp->phase == dyPRIMAL2)
	  warn(323,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
	       lbk,xbk,ubk,lbk-xbk,toll) ;
	break ; }
      case vstatBFR:
      case vstatNBFR:
      { if (ubk < dy_tols->inf || lbk > -dy_tols->inf)
	{ errmsg(334,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,lbk,ubk,statk) ;
	  retval = FALSE ; }
	break ; }
      case vstatNBFX:
      { if (lbk != ubk)
	{ errmsg(331,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,xbk,statk,
		 lbk,ubk,ubk-lbk,0.0) ;
	  retval = FALSE ;
	  break ; }
	if (xk != ubk)
	{ errmsg(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 xk,ubk,"lb = ub",ubk-xk,0.0) ;
	  retval = FALSE ; }
	break ; }
      case vstatNBUB:
      { if (xk != ubk)
	{ errmsg(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 xk,ubk,"ub",ubk-xk,0.0) ;
	  retval = FALSE ; }
	break ; }
      case vstatNBLB:
      { if (xk != lbk)
	{ errmsg(333,rtnnme,dy_sys->nme,phase,iter,nmek,xkndx,statk,
		 xk,lbk,"lb",xk-lbk,0.0) ;
	  retval = FALSE ; }
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (FALSE) ; } }
/*
  Final checks:
  * Superbasic variables should only show up in primal phase II, where we're
    concerned about maintaining feasibility. During primal phase I or dual
    simplex, we force nonbasic variables to bound and live with any loss of
    feasibility.
  * Nonbasic free variables should not occur during dual simplex -- we
    shouldn't have been able to make a feasible start if they're present.
*/
    switch (xkstatus)
    { case vstatSB:
      { if (dy_lp->phase != dyPRIMAL2)
	{ errmsg(346,rtnnme,dy_sys->nme,phase,iter,statk,nmek,xkndx) ;
	  retval = FALSE ; }
	break ; }
      case vstatNBFR:
      { if (dy_lp->phase == dyDUAL)
	{ errmsg(346,rtnnme,dy_sys->nme,phase,iter,statk,nmek,xkndx) ;
	  retval = FALSE ; }
	break ; } } }

  return (retval) ; }




void dy_chkdual (int lvl)

/*
  This routine checks the dual variables and reduced costs for correctness,
  and also checks for dual feasibility (agreement of reduced costs and status).

  If dual antidegeneracy is inactive, the routine is checking for
    * numerical agreement between the iteratively updated duals and reduced
      costs (dy_y, dy_cbar) and the fresh values calculated here;
    * numerical agreement between the reduced costs of logicals and the
      associated duals; and
    * logical agreement between duals, reduced costs, and primal status.

  When dual antidegeneracy is active, the routine allows for differences
  between the iterative and calculated values when the column is part of a
  restricted subproblem. It checks that the freshly calculated values are
  zero (otherwise the column should not be involved in the restricted
  subproblem). Unfortunately, given the way that the perturbation is applied
  (directly to dy_cbar), there's no way I know of to check for correct values
  of perturbed reduced costs and duals.

  Parameters:
    lvl:	0: suppresses all printing (allows or convenient suppression
		   of messages without recompiling)
		1: prints a summary only if problems are detected
		2: prints a warning for each problem detected
*/

{ int i,j,m,n ;
  flags statj ;
  double cbarj,yi,diff,tol ;
  bool degenActive ;
  double *y,*cbar ;

  int yerrcnt,cbarerrcnt,ycbarerrcnt,cbarstatuserrcnt ;
  double yerrtot,cbarerrtot ;

/*
  Scaling can tighten the base value of tols.cost, and we don't want that
  here. Hardwire the default value.
*/
  const double base_tol = 1.0e-11 ;

  const char *rtnnme = "dy_chkdual" ;

  if (dy_lp->degen > 0 && dy_lp->phase == dyDUAL)
  { degenActive = TRUE ; }
  else
  { degenActive = FALSE ; }

  m = dy_sys->concnt ;
  n = dy_sys->varcnt ;
/*
  Start out by calculating what we think the duals and reduced costs should be,
  independent of the current state of the dual degeneracy algorithm.
*/
  y = (double *) MALLOC((m+1)*sizeof(double)) ;
  cbar = (double *) MALLOC((n+1)*sizeof(double)) ;

  for (i = 1 ; i <= m ; i++)
  { j = dy_basis[i] ;
    y[i] = dy_sys->obj[j] ; }
  dy_btran(y) ;

  for (j = 1 ; j <= n ; j++)
  { cbarj = consys_dotcol(dy_sys,j,y) ;
    cbar[j] = dy_sys->obj[j]-cbarj ; }
/*
  Now see if we agree on the values for the duals. When dual antidegeneracy
  is inactive, or for columns not involved in the restricted subproblem,
  we're simply checking our fresh values against the iteratively updated
  values.

  When dual antidegeneracy is active, in a position that's involved in a
  restricted subproblem, the real dual should be 0.  Since reduced
  subproblems are nested, it's sufficient to check dy_ddegenset == 0.
*/
  tol = 1000*base_tol ;
  yerrcnt = 0 ;
  yerrtot = 0.0 ;
  for (i = 1 ; i <= m ; i++)
  { j = dy_basis[i] ;
    if (dy_ddegenset[j] == 0)
    { diff = fabs(dy_y[i]-y[i]) ;
      if (diff > tol)
      { if (lvl >= 2)
	{ warn(321,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       (dy_lp->degen == 0)?"real":"uninvolved",
	       "dy_y",i,dy_y[i],y[i],diff,tol) ; }
	yerrcnt++ ;
	yerrtot += diff ; } }
    else
    { if (fabs(y[i]) > tol)
      { if (lvl >= 2)
	{ warn(321,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "real","y",i,y[i],0.0,y[i],tol) ; }
	yerrcnt++ ;
	yerrtot += fabs(y[i]) ; } } }
/*
  Now try for the reduced costs. As with the duals, if the column is involved
  in the reduced subproblem, the real reduced cost should be 0. If not, they
  should be identical.
*/
  tol = 1000*base_tol ;
  cbarerrcnt = 0 ;
  cbarerrtot = 0.0 ;
  for (j = 1 ; j <= n ; j++)
  { if (dy_ddegenset[j] == 0)
    { diff = fabs(dy_cbar[j]-cbar[j]) ;
      if (diff > tol)
      { if (lvl >= 2)
	{ warn(321,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       (dy_lp->degen == 0)?"real":"uninvolved",
	       "dy_cbar",j,dy_cbar[j],cbar[j],diff,tol) ; }
      cbarerrcnt++ ;
      cbarerrtot += diff ; } }
    else
    { if (fabs(cbar[j]) > tol)
      { if (lvl >= 2)
	{ warn(321,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "real","cbar",j,cbar[j],0.0,cbar[j],tol) ; }
      cbarerrcnt++ ;
      cbarerrtot += fabs(cbar[j]) ; } } }
/*
  Check that, for 1 < i < m, that the reduced cost of the logical
  equals the negative of the dual variable for the constraint for variables
  not involved in a degenerate set. Since we perturbed cbar<i> as a surrogate
  for perturbing c<i>, it's no longer true that cbar<i> = -dot(y,a<i>); we now
  have cbar<i> = c<i> - dot(y,a<i>). Unfortunately, we don't know c<i> (though
  we could calculate it).
*/
  ycbarerrcnt = 0 ;
  for (i = 1 ; i <= m ; i++)
  { statj = dy_status[i] ;
    if (dy_ddegenset[i] == 0)
    { cbarj = cbar[i] ;
      yi = y[i] ;
      if (fabs(cbarj+yi) > base_tol)
      { if (lvl >= 2)
	{ warn(336,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       consys_nme(dy_sys,'v',i,FALSE,NULL),i,dy_prtvstat(statj),
	       "real",-yi,cbarj,fabs(cbarj+yi),base_tol) ; }
	ycbarerrcnt++ ; } } }
/*
  Check that primal status agrees with reduced costs. If antidegeneracy is
  active, the real values should still be in agreement (zero is neutral, hence
  should not cause an error).
*/
  cbarstatuserrcnt = 0 ;
  for (j = 1 ; j <= n ; j++)
  { statj = dy_status[j] ;
    cbarj = dy_cbar[j] ;
    if ((flgon(statj,vstatNBLB) && cbarj < -dy_tols->dfeas) ||
	(flgon(statj,vstatNBUB) && cbarj > dy_tols->dfeas) ||
	flgon(statj,vstatNBFR|vstatSB))
    { if (lvl >= 2)
      { warn(347,rtnnme,dy_sys->nme,
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	       dy_prtvstat(statj),j,cbarj,dy_tols->dfeas) ; }
      cbarstatuserrcnt++ ; }
    if (degenActive && dy_ddegenset[j] < dy_lp->degen)
    { cbarj = cbar[j] ;
      if ((flgon(statj,vstatNBLB) && cbarj < -dy_tols->dfeas) ||
	  (flgon(statj,vstatNBUB) && cbarj > dy_tols->dfeas) ||
	  flgon(statj,vstatNBFR|vstatSB))
      { if (lvl >= 2)
	{ warn(347,rtnnme,dy_sys->nme,
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		 dy_prtvstat(statj),j,cbarj,dy_tols->dfeas) ; }
	cbarstatuserrcnt++ ; } } }
/*
  Summary. Print if there's any error.
*/
  if (lvl >= 1 && yerrcnt+cbarerrcnt+ycbarerrcnt+cbarstatuserrcnt > 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  %s: [%s]: (%s)%d: ",rtnnme,
	        dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		dy_lp->tot.iters) ;
    if (yerrcnt > 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "%d y value errors (%g)", yerrcnt,yerrtot) ; }
    if (cbarerrcnt > 0)
    { if (yerrcnt > 0) dyio_outfmt(dy_logchn,dy_gtxecho,", ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%d cbar value errors (%g)",
		  cbarerrcnt,cbarerrtot) ; }
    if (ycbarerrcnt > 0)
    { if (yerrcnt+cbarerrcnt > 0) dyio_outfmt(dy_logchn,dy_gtxecho,", ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "%d y/cbar agreement errors",ycbarerrcnt) ; }
    if (cbarstatuserrcnt > 0)
    { if (yerrcnt+cbarerrcnt+ycbarerrcnt > 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,", ") ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "%d cbar/status agreement errors",cbarstatuserrcnt) ; }
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
/*
  Clean up and we're done.
*/
  FREE(cbar) ;
  FREE(y) ;

  return ; }




bool dy_chkdysys (consys_struct *orig_sys)

/*
  This routine is a paranoid check routine. It checks that dy_sys is really
  the portion of orig_sys specified in dy_origvars and dy_origcons (and that
  this agrees with dy_actvars and dy_actcons).

  Parameters:
    orig_sys:	the original constraint system
  
  Returns: TRUE if the systems seem to correspond, FALSE otherwise.
*/

{ int ovndx,ocndx,avndx,acndx,opkndx,apkndx,cnt ;
  double oaij,aaij,rhscorr,*obj,scaledtol ;
  pkvec_struct *opkcon,*apkcon,*apkvar ;
  flags vstatus ;
  bool retval ;
  char onmebuf[128],ovecnme[128] ;
  const char *rtnnme = "dy_chkdysys" ;

/*
  The first thing we'll do is run through dy_origvars and dy_origcons,
  count up the number of active entries, and check that this agrees with
  the size of dy_sys. At the same time, cross-check with dy_actvars and
  dy_actcons. While we're walking the columns, check that dy_sys and orig_sys
  agree on upper/lower bounds and objective function coefficients. We have
  to be a bit careful. For bounds, the initial check for equality handles the
  case where the bound is +/-inf. For the objective, clearly the phase I
  objective isn't going to match, so make sure we're looking at the phase II
  objective.
*/
  cnt = 0 ;
  if (dy_lp->p1obj.installed == TRUE)
    obj = dy_lp->p1obj.p2obj ;
  else
    obj = dy_sys->obj ;
  for (ovndx = 1 ; ovndx <= orig_sys->varcnt ; ovndx++)
  { if (dy_origvars[ovndx] <= 0)
    { vstatus = (flags) -dy_origvars[ovndx] ;
      if (!(vstatus == vstatNBFR || vstatus == vstatNBFX ||
	    vstatus == vstatNBUB || vstatus == vstatNBLB))
      { errmsg(433,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "inactive",consys_nme(orig_sys,'v',ovndx,TRUE,NULL),ovndx,
	       dy_prtvstat(dy_origvars[ovndx])) ; 
	return (FALSE) ; } }
    else
    { cnt++ ;
      avndx = dy_origvars[ovndx] ;
      if (ovndx != dy_actvars[avndx])
      { errmsg(356,rtnnme,orig_sys->nme,"variable",ovndx,avndx,avndx,
	       dy_actvars[avndx]) ;
	return (FALSE) ; }
      if (!(dy_sys->vub[avndx] == orig_sys->vub[ovndx] ||
	    withintol(dy_sys->vub[avndx],orig_sys->vub[ovndx],dy_tols->zero)))
      { (void) consys_nme(orig_sys,'v',ovndx,TRUE,onmebuf) ;
	strcpy(ovecnme,consys_assocnme(orig_sys,CONSYS_VUB)) ;
	errmsg(367,rtnnme,dy_sys->nme,consys_assocnme(dy_sys,CONSYS_VUB),
	       consys_nme(dy_sys,'v',avndx,TRUE,NULL),avndx,
	       dy_sys->vub[avndx],orig_sys->vub[ovndx],ovecnme,onmebuf,ovndx,
	       dy_sys->vub[avndx]-orig_sys->vub[ovndx],dy_tols->zero) ;
	return (FALSE) ; }
      if (!(dy_sys->vlb[avndx] == orig_sys->vlb[ovndx] ||
	    withintol(dy_sys->vlb[avndx],orig_sys->vlb[ovndx],dy_tols->zero)))
      { (void) consys_nme(orig_sys,'v',ovndx,TRUE,onmebuf) ;
	strcpy(ovecnme,consys_assocnme(orig_sys,CONSYS_VLB)) ;
	errmsg(367,rtnnme,dy_sys->nme,consys_assocnme(dy_sys,CONSYS_VLB),
	       consys_nme(dy_sys,'v',avndx,TRUE,NULL),avndx,
	       dy_sys->vlb[avndx],orig_sys->vlb[ovndx],ovecnme,onmebuf,ovndx,
	       dy_sys->vlb[avndx]-orig_sys->vlb[ovndx],dy_tols->zero) ;
	return (FALSE) ; }
      if (!withintol(obj[avndx],orig_sys->obj[ovndx],dy_tols->cost))
      { (void) consys_nme(orig_sys,'v',ovndx,TRUE,onmebuf) ;
	strcpy(ovecnme,consys_assocnme(orig_sys,CONSYS_OBJ)) ;
 	errmsg(367,rtnnme,dy_sys->nme,consys_assocnme(dy_sys,CONSYS_OBJ),
	       consys_nme(dy_sys,'v',avndx,TRUE,NULL),avndx,
	       obj[avndx],orig_sys->obj[ovndx],ovecnme,onmebuf,ovndx,
	       obj[avndx]-orig_sys->obj[ovndx],dy_tols->cost) ;
	return (FALSE) ; } } }
  if (cnt != dy_sys->archvcnt)
  { errmsg(361,rtnnme,dy_sys->nme,"variable",dy_sys->archvcnt,cnt) ;
    return (FALSE) ; }

  cnt = 0 ;
  for (ocndx = 1 ; ocndx <= orig_sys->concnt ; ocndx++)
  { if (dy_origcons[ocndx] > 0)
    { cnt++ ;
      acndx = dy_origcons[ocndx] ;
      if (ocndx != dy_actcons[acndx])
      { errmsg(356,rtnnme,orig_sys->nme,"constraint",ocndx,acndx,acndx,
	       dy_actcons[acndx]) ;
	return (FALSE) ; } } }
  if (cnt != dy_sys->concnt)
  { errmsg(361,rtnnme,dy_sys->nme,"constraint",dy_sys->concnt,cnt) ;
    return (FALSE) ; }
/*
  Now check each active constraint to make sure that it contains exactly the
  coefficients it should. In the process, check the logical for each constraint
  and make sure that it has exactly one coefficient and the proper bounds.
*/
  retval = TRUE ;
  apkcon = NULL ;
  apkvar = NULL ;
  opkcon = NULL ;
  for (acndx = 1 ; acndx <= dy_sys->concnt && retval == TRUE ; acndx++)
  { if (consys_getrow_pk(dy_sys,acndx,&apkcon) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,"constraint",
	     consys_nme(dy_sys,'c',acndx,TRUE,NULL),acndx) ;
      retval = FALSE ;
      break ; }
    ocndx = dy_actcons[acndx] ;
    if (consys_getrow_pk(orig_sys,ocndx,&opkcon) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,"constraint",
	     consys_nme(orig_sys,'c',ocndx,TRUE,NULL),ocndx) ;
      retval = FALSE ;
      break ; }
    rhscorr = 0 ;
    cnt = 0 ;
    for (opkndx = 0 ; opkndx < opkcon->cnt ; opkndx++)
    { ovndx = opkcon->coeffs[opkndx].ndx ;
      oaij = opkcon->coeffs[opkndx].val ;
      if (dy_origvars[ovndx] < 0)
      { cnt++ ;
	vstatus = (flags) -dy_origvars[ovndx] ;
	switch (vstatus)
	{ case vstatNBLB:
	  { rhscorr += oaij*orig_sys->vlb[ovndx] ;
	    break ; }
	  case vstatNBUB:
	  case vstatNBFX:
	  { rhscorr += oaij*orig_sys->vub[ovndx] ;
	    break ; } } }
      else
      { avndx = dy_origvars[ovndx] ;
	for (apkndx = 0 ; apkndx < apkcon->cnt ; apkndx++)
	{ if (apkcon->coeffs[apkndx].ndx == avndx)
	  { aaij = apkcon->coeffs[apkndx].val ;
	    break ; } }
	if (apkndx >= apkcon->cnt)
	{ errmsg(362,rtnnme,dy_sys->nme,
		 consys_nme(dy_sys,'v',avndx,FALSE,NULL),avndx,
		 consys_nme(dy_sys,'c',acndx,FALSE,NULL),acndx) ;
	  retval = FALSE ;
	  break ; }
	if (aaij != oaij)
	{ errmsg(364,rtnnme,dy_sys->nme,
		 consys_nme(dy_sys,'c',acndx,FALSE,NULL),
		 consys_nme(dy_sys,'v',avndx,FALSE,NULL),
		 acndx,avndx,aaij,
		 consys_nme(orig_sys,'c',ocndx,FALSE,NULL),
		 consys_nme(orig_sys,'v',ovndx,FALSE,NULL),
		 ocndx,ovndx,oaij,orig_sys->nme) ;
	  retval = FALSE ;
	  break ; } } }
    if (retval == FALSE) break ;
    if (opkcon->cnt-cnt != apkcon->cnt-1)
    { errmsg(364,rtnnme,dy_sys->nme,consys_nme(dy_sys,'c',acndx,FALSE,NULL),
	     acndx,opkcon->cnt-cnt,apkcon->cnt-1) ;
      retval = FALSE ;
      break ; }
/*
  We have exactly the right architectural coefficients. What about the
  logical? First, the corresponding entry in actvars should be 0 (logicals
  don't exist in orig_sys). Find its coefficient (there should be exactly one
  in the column), then check that the logical has the right bounds. Recall
  that dylp isn't expecting to see anything except <= and = constraints, with
  range constraints handled as a <= constraint with an upper bounded slack.
*/
    if (dy_actvars[acndx] != -INT_MAX)
    { errmsg(360,rtnnme,dy_sys->nme,acndx,dy_actvars[acndx]) ;
      retval = FALSE ; }
    if (consys_getcol_pk(dy_sys,acndx,&apkvar) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,"variable",
	     consys_nme(dy_sys,'v',acndx,TRUE,NULL),acndx) ;
      retval = FALSE ;
      break ; }
    if (apkvar->cnt != 1)
    { errmsg(402,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,consys_nme(dy_sys,'v',acndx,FALSE,NULL),
	     acndx,apkvar->cnt) ;
      retval = FALSE ;
      break ; }
    for (apkndx = 0 ; apkndx < apkcon->cnt ; apkndx++)
    { if (apkcon->coeffs[apkndx].ndx == acndx) break ; }
    if (apkndx >= apkcon->cnt)
    { errmsg(365,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',acndx,FALSE,NULL),
	     acndx,consys_nme(dy_sys,'c',acndx,FALSE,NULL),acndx) ;
      retval = FALSE ;
      break ; }
    if (dy_sys->vlb[acndx] != 0)
    { errmsg(426,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',acndx,FALSE,NULL),
	     acndx,consys_prtcontyp(dy_sys->ctyp[acndx]),
	     consys_nme(dy_sys,'c',acndx,FALSE,NULL),acndx,"lower",
	     dy_sys->vlb[acndx]) ;
      retval = FALSE ; }
    switch (dy_sys->ctyp[acndx])
    { case contypLE:
      { break ; }
      case contypEQ:
      { if (dy_sys->vub[acndx] != 0)
	{ errmsg(426,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',acndx,FALSE,NULL),
		 acndx,consys_prtcontyp(dy_sys->ctyp[acndx]),
		 consys_nme(dy_sys,'c',acndx,FALSE,NULL),acndx,"upper",
		 dy_sys->vub[acndx]) ;
	  retval = FALSE ; }
	break ; }
      case contypRNG:
      { if (!withintol(dy_sys->vub[acndx],
		       dy_sys->rhs[acndx]-dy_sys->rhslow[acndx],dy_tols->zero))
	{ errmsg(426,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',acndx,FALSE,NULL),
		 acndx,consys_prtcontyp(dy_sys->ctyp[acndx]),
		 consys_nme(dy_sys,'c',acndx,FALSE,NULL),acndx,"upper",
		 dy_sys->vub[acndx]) ;
	  retval = FALSE ; }
	break ; }
      default:
      { errmsg(463,rtnnme,dy_sys->nme,consys_prtcontyp(dy_sys->ctyp[acndx])) ;
	retval = FALSE ;
	break ; } }
    scaledtol = dy_tols->zero*(1+maxx(fabs(dy_sys->rhs[acndx]),fabs(rhscorr))) ;
    if (!withintol(dy_sys->rhs[acndx],orig_sys->rhs[ocndx]-rhscorr,scaledtol))
    { if (withintol(dy_sys->rhs[acndx],
		    orig_sys->rhs[ocndx]-rhscorr,10*scaledtol))
      { warn(366,rtnnme,dy_sys->nme,acndx,dy_sys->rhs[acndx],
	     consys_nme(dy_sys,'c',acndx,FALSE,NULL),
	     orig_sys->rhs[ocndx]-rhscorr,
	     fabs(dy_sys->rhs[acndx]-(orig_sys->rhs[ocndx]-rhscorr)),
	     scaledtol,ocndx,
	     consys_nme(orig_sys,'c',ocndx,FALSE,NULL),orig_sys->rhs[ocndx],
	     -rhscorr) ; }
      else
      { errmsg(366,rtnnme,dy_sys->nme,acndx,dy_sys->rhs[acndx],
	       consys_nme(dy_sys,'c',acndx,FALSE,NULL),
	       orig_sys->rhs[ocndx]-rhscorr,
	       fabs(dy_sys->rhs[acndx]-(orig_sys->rhs[ocndx]-rhscorr)),
	       scaledtol,ocndx,
	       consys_nme(orig_sys,'c',ocndx,FALSE,NULL),orig_sys->rhs[ocndx],
	       -rhscorr) ;
	retval = FALSE ;
	break ; } } }
/*
  Clean up and return.
*/
  if (apkcon != NULL) pkvec_free(apkcon) ;
  if (apkvar != NULL) pkvec_free(apkvar) ;
  if (opkcon != NULL) pkvec_free(opkcon) ;

  return (retval) ; }


#endif /* PARANOIA */

bool dy_dupbasis (int dst_basissze, basis_struct **p_dst_basis,
		  basis_struct *src_basis, int dst_statussze,
		  flags **p_dst_status, int src_statuslen, flags *src_status)

/*
  This routine duplicates the basis and status passed in src_*, leaving a
  copy in the basis and status pointed to by dst_*. Note that we're working
  with the external form of the basis and status --- this routine is intended
  as a convenience for the client code.
 
  If either of p_dst_basis or p_dst_status is null, no copy is made of the
  basis or status, respectively. It's an error if both are null.

  We never shrink an existing vector. If a basis is supplied and it is
  already longer than the requested size, the length is not changed. If
  max(requested size, existing size) < source size, the requested size is
  boosted to be sufficient.

  Parameters:
    dst_basissze: the desired allocated size for the new basis structure
    p_dst_basis: (i) an empty basis structure (space will be allocated as
		     required)
		  (o) duplicate of src_basis
    src_basis:	basis structure to be duplicated
    dst_statussze: the desired allocated size for the new status vector
    p_dst_status: (i) an empty status vector (space will be allocated as
		      required)
    src_statuslen: number of entries in src_status
    src_status:	status vector to be copied
  
  Returns: TRUE if the copy is successful, FALSE otherwise (FALSE is possible
	   only if we're paranoid).
*/

{ basis_struct *dst_basis ;
  bool want_basis,want_status ;


# ifdef PARANOIA
  const char *rtnnme = "dy_dupbasis" ;

  if (p_dst_basis == NULL && p_dst_status == NULL)
  { errmsg(2,rtnnme,"&dst_basis and &dst_status") ;
    return (FALSE) ; }
  if (p_dst_basis != NULL && src_basis == NULL )
  { errmsg(403,rtnnme,"basis") ;
    return (FALSE) ; }
  if (p_dst_status != NULL && src_status == NULL)
  { errmsg(403,rtnnme,"status") ;
    return (FALSE) ; }
# endif

  if (p_dst_basis != NULL)
  { want_basis = TRUE ; }
  else
  { want_basis = FALSE ; }
  if (p_dst_status != NULL)
  { want_status = TRUE ; }
  else
  { want_status = FALSE ; }
/*
  If the user's supplied a basis_struct, check it for size against src_basis
  and the requested allocated size, and force a reallocation of the basis
  vector if necessary. Then make the copy.
*/
  if (want_basis == TRUE)
  { if (*p_dst_basis != NULL)
    { dst_basis = *p_dst_basis ;
      if (dst_basissze < dst_basis->len) dst_basissze = dst_basis->len ;
      if ((dst_basis->len < src_basis->len ||
	   dst_basis->len < dst_basissze) && dst_basis->el != NULL)
      { FREE(dst_basis->el) ;
	dst_basis->el = NULL ; } }
    else
    { dst_basis = (basis_struct *) CALLOC(1,sizeof(basis_struct)) ;
      *p_dst_basis = dst_basis ; }

    if (dst_basissze < src_basis->len)
    {
#     ifdef PARANOIA
      warn(404,rtnnme,"basis",dst_basissze,src_basis->len) ;
#     endif
      dst_basissze = src_basis->len ; }
    if (dst_basis->el == NULL)
    { dst_basis->el = (basisel_struct *)
		      MALLOC((dst_basissze+1)*sizeof(basisel_struct)) ; }

    memcpy(dst_basis->el,src_basis->el,
	   (src_basis->len+1)*sizeof(basisel_struct)) ;
    dst_basis->len = src_basis->len ; }
/*
  And similarly for status. If the client has supplied a status vector, we
  assume it's big enough.
*/
  if (want_status == TRUE)
  { if (*p_dst_status == NULL)
    {
      if (dst_statussze < src_statuslen)
      {
#       ifdef PARANOIA
	warn(404,rtnnme,"status",dst_statussze,src_statuslen) ;
#       endif
	dst_statussze = src_statuslen ; }
      *p_dst_status = (flags *) MALLOC((dst_statussze+1)*sizeof(flags)) ; }
    memcpy(*p_dst_status,src_status,(src_statuslen+1)*sizeof(flags)) ; }
  
  return (TRUE) ; }



static void build_soln (lpprob_struct *orig_lp)

/*
  This routine is responsible for constructing a solution that can be returned
  to the client. The solution consists of the status of all architectural
  variables, a basis vector comprised of (constraint, variable) index pairs,
  the values of the basic variables, and the values of the dual variables.

  Note that this routine will execute even for error returns, just in case
  the client wants to do a postmortem.

  For the status, the normal dylp flag codes are used for nonbasic variables,
  and the negative of the basis position is used for basic variables. Since
  we don't keep explicit values for nonbasic variables, they must be at bound.
  If dylp's pivoting works as advertised, there will be no nonbasic free or
  superbasic variables in an optimal solution, so this won't be a problem.
  (There is the possibility they will occur in an infeasible or unbounded
  problem. C'est la vie.)

  Basic logical variables, which don't exist outside of dylp, are reported
  in the basis vector as the negative of the index of their associated
  constraint.

  Parameters:
    orig_lp:	The LP problem passed by the client.

  Returns: undefined
*/

{ int ovndx,xjndx,xjpos,ubndndx,orig_ubndndx ;
  flags xjstatus ;
  consys_struct *orig_sys ;
  const char *rtnnme = "build_soln" ;

  /* scaling.c */
  extern void dy_unscale_soln(double *x, double *y) ;

/*
  Grab the necessary space, if the user hasn't provided it or the space
  provided is too small. The basis, x, and y arrays have to be big enough to
  hold information for each constraint in the active system, but the status
  and active variable arrays have to be big enough to hold information for
  each variable in the original system. Note that we allocate to colsze or
  rowsze, even if the indicated size is larger than necessary. We have to do
  this to maintain consistency in allocated size between preallocated vectors
  and vectors we allocate here.

  Note that the basis_struct has a length field, but this is the actual length
  of the basis, not the allocated length. The allocated length is assumed to
  match rowsze.
*/
  orig_sys = orig_lp->consys ;
  if (orig_lp->colsze < orig_sys->varcnt)
  { if (orig_lp->status != NULL)
    { FREE(orig_lp->status) ;
      orig_lp->status = NULL ; }
    if (orig_lp->actvars != NULL)
    { FREE(orig_lp->actvars) ;
      orig_lp->actvars = NULL ; }
    orig_lp->colsze = orig_sys->varcnt ; }
  if (orig_lp->status == NULL)
    orig_lp->status = (flags *) CALLOC((orig_lp->colsze+1),sizeof(flags)) ;
  if (flgon(orig_lp->ctlopts,lpctlACTVARSOUT))
  { if (orig_lp->actvars == NULL)
      orig_lp->actvars = (bool *) CALLOC((orig_lp->colsze+1),sizeof(bool)) ; }

  if (orig_lp->rowsze < dy_sys->concnt)
  { if (orig_lp->x != NULL)
    { FREE(orig_lp->x) ;
      orig_lp->x = NULL ; }
    if (orig_lp->y != NULL)
    { FREE(orig_lp->y) ;
      orig_lp->y = NULL ; }
    if (orig_lp->basis != NULL)
    { if (orig_lp->basis->el != NULL)
      { FREE(orig_lp->basis->el) ;
	orig_lp->basis->el = NULL ; } }
    orig_lp->rowsze = dy_sys->concnt ; }
  if (orig_lp->basis == NULL)
    orig_lp->basis = (basis_struct *) CALLOC(sizeof(basis_struct),1) ;
  if (orig_lp->basis->el == NULL)
    orig_lp->basis->el =
	(basisel_struct *) CALLOC((orig_lp->rowsze+1),sizeof(basisel_struct)) ;
  if (orig_lp->x == NULL)
    orig_lp->x = (double *) CALLOC((orig_lp->rowsze+1),sizeof(double)) ;
  if (orig_lp->y == NULL)
    orig_lp->y = (double *) CALLOC((orig_lp->rowsze+1),sizeof(double)) ;

/*
  Fill in the simple stuff --- objective, return code, and iterations. For an
  optimal solution, we'll recalculate the objective on the spot, as it's just
  been incrementally updated through primal II or dual II. If we stopped due
  to hitting an iteration limit, we also want a valid solution, but it's a
  bit more complicated (we may need to flush the P1 objective and recalculate
  duals and reduced costs). A small iteration limit is a common tactic for
  clients to limit computational effort; they still expect a valid objective.

  When the solution is unbounded, the objective is overloaded to hold the
  index of the variable that was discovered to be unbounded (which we need to
  convert to the orig_sys frame, with the usual convention that a logical is
  represented as the negative of the index of it's constraint).
*/
  switch (dy_lp->lpret)
  { case lpOPTIMAL:
    { orig_lp->obj = dy_calcobj() ;
      break ; }
    case lpITERLIM:
    { if (dy_lp->p1obj.installed == TRUE)
      { if (dy_swapobjs(dyPRIMAL2) == FALSE)
	{ errmsg(318,rtnnme,
		 dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,"remove") ;
	  dy_lp->lpret = lpFATAL ; }
	else
	{ dy_calcduals() ;
	  if (dy_calccbar() == FALSE)
	  { errmsg(384,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters) ;
	    dy_lp->lpret = lpFATAL ; } } }
      if (dy_lp->lpret == lpFATAL)
      { orig_lp->obj = dy_tols->inf ; }
      break ; }
    case lpUNBOUNDED:
    { ubndndx = abs(dy_lp->ubnd.ndx) ;
      if (ubndndx > dy_sys->concnt)
      { orig_ubndndx = dy_actvars[ubndndx] ; }
      else
      { orig_ubndndx = -dy_actcons[ubndndx] ; }
      orig_lp->obj = (double) orig_ubndndx ;
      break ; }
    case lpINFEAS:
    { orig_lp->obj = dy_lp->infeas ;
      break ; }
    default:
    { orig_lp->obj = dy_tols->inf ;
      break ; } }
  orig_lp->lpret = dy_lp->lpret ;
  orig_lp->iters = dy_lp->tot.pivs ;
  
/*
  Now fill in the answer. First we walk dy_origvars. For each nonbasic
  architectural variable (active and inactive), we fill in the status entry.
  For each active basic architectural variable, we place the negative of the
  basis position in the status entry and fill in the basis information, then
  fill in the x and y values. Remember that there are no logicals in orig_sys.
  If the user is asking for an active variables vector, fill that in too.
*/
  for (ovndx = 1 ; ovndx <= orig_sys->varcnt ; ovndx++)
  { if (dy_origvars[ovndx] < 0)
    { xjstatus = (flags)(-dy_origvars[ovndx]) ;
      orig_lp->status[ovndx] = xjstatus ;
      if (flgon(orig_lp->ctlopts,lpctlACTVARSOUT))
	orig_lp->actvars[ovndx] = FALSE ; }
    else
    { if (flgon(orig_lp->ctlopts,lpctlACTVARSOUT))
	orig_lp->actvars[ovndx] = TRUE ; 
      xjndx = dy_origvars[ovndx] ;
      xjstatus = dy_status[xjndx] ;
      if (flgon(xjstatus,vstatNONBASIC|vstatEXOTIC))
      { if (flgon(xjstatus,vstatNBFR|vstatSB))
	{ if (dy_lp->lpret == lpOPTIMAL)
	    errmsg(359,rtnnme,dy_sys->nme,
		   consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,
		   dy_prtvstat(xjstatus)) ;
#         ifndef DYLP_NDEBUG
	  else
	    warn(359,rtnnme,dy_sys->nme,
		 consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,
		 dy_prtvstat(xjstatus)) ;
#         endif
	}
	orig_lp->status[ovndx] = xjstatus ; }
      else
      { xjpos = dy_var2basis[xjndx] ;
	orig_lp->status[ovndx] = -xjpos ;
	orig_lp->basis->el[xjpos].cndx = dy_actcons[xjpos] ;
	orig_lp->basis->el[xjpos].vndx = ovndx ;
	orig_lp->x[xjpos] = dy_x[xjndx] ;
	orig_lp->y[xjpos] = dy_y[xjpos] ; } } }
/*
  Now, we walk the basis, to pick up any entries occupied by basic logicals.
*/
  for (xjpos = 1 ; xjpos <= dy_sys->concnt ; xjpos++)
  { if (dy_basis[xjpos] <= dy_sys->concnt)
    { orig_lp->basis->el[xjpos].cndx = dy_actcons[xjpos] ;
      orig_lp->basis->el[xjpos].vndx = -dy_actcons[dy_basis[xjpos]] ;
      orig_lp->x[xjpos] = dy_xbasic[xjpos] ;
      orig_lp->y[xjpos] = dy_y[xjpos] ; } }
/*
  Unscale the solution.
*/
  dy_unscale_soln(orig_lp->x,orig_lp->y) ;
/*
  Set the basis length and we're out of here.
*/
  orig_lp->basis->len = dy_sys->concnt ;

  return ; }



/*
  Handy macro for use by dy_finishup.
*/

#define FREE_AND_CLEAR(zz_var_zz) \
  if (zz_var_zz != NULL) { FREE(zz_var_zz) ; zz_var_zz = NULL ; }

void dy_finishup (lpprob_struct *orig_lp, dyphase_enum phase)

/*
  This routine assembles the final answer from the active and inactive portions
  of the constraint system, and cleans up the working structures of dylp.

  Parameters:
    orig_lp:	The LP problem passed by the client.
    phase:	dylp's phase

  Returns: undefined
*/

{ consys_struct *orig_sys ;

  extern bool dy_retained ;		/* dylp.c */

/*
  Why are we here?

  It may be that the user called dylp strictly to free the data structures;
  this is indicated by orig_lp->phase == dyDONE, phase == dyINV, and the
  lpctlONLYFREE flag.

  It may be that dylp needs to free up previous data structures to do a warm
  or cold start. This is indicated by orig_lp->phase == dyINIT and phase ==
  dyINIT and orig_lp->lpret == dyINV.

  Otherwise, if orig_lp->phase is dyDONE,  we're cleanly done and we'll try
  to build an answer as best we can.  `Cleanly done' covers error and
  non-error returns.

  Anything else is a serious error, and we won't even attempt to build a
  solution, just free the data structures.
*/
  if (orig_lp->phase == dyDONE && phase == dyINV &&
      flgon(orig_lp->ctlopts,lpctlONLYFREE))
  { clrflg(orig_lp->ctlopts,lpctlACTVARSOUT) ;
    dy_initlclsystem(orig_lp,TRUE) ; }
  else
  if (orig_lp->phase == dyINIT && phase == dyINIT && orig_lp->lpret == dyINV)
  { dy_initlclsystem(orig_lp,TRUE) ; }
  else
  if (orig_lp->phase == dyDONE)
  { build_soln(orig_lp) ; }
  else
  { clrflg(orig_lp->ctlopts,lpctlACTVARSOUT) ; }
/*
  Free up the working structures --- the constraint system and the various
  local vectors, the LP control structure, and the option and tolerance
  structures. We also need to detach origvars and origcons from orig_sys, so
  that we don't keep accumulating attachment structures. The code is
  defensive, setting all pointers back to NULL so that we don't get in
  trouble in an environment (like COIN) where multiple cleanup calls are
  possible from independent objects.

  The exception is if we see the lpctlNOFREE flag and we've come to one of
  the mathematically valid results --- optimal, unbounded, or infeasible.
  In this case, we can support reoptimisation with a hot start. Add to these
  the case of lpITERLIM, which is a valid state of affairs and occurs fairly
  commonly in branch & cut applications (e.g., strong branching).
*/
  if (flgoff(orig_lp->ctlopts,lpctlNOFREE) ||
      (!(orig_lp->lpret == lpOPTIMAL || orig_lp->lpret == lpUNBOUNDED ||
	 orig_lp->lpret == lpINFEAS || orig_lp->lpret == lpITERLIM)))
  { orig_sys = orig_lp->consys ;
    if (dy_origvars != NULL)
    { (void) consys_detach(orig_sys,(void **) &dy_origvars,TRUE) ;
      FREE(dy_origvars) ;
      dy_origvars = NULL ; }
    if (dy_origcons != NULL)
    { (void) consys_detach(orig_sys,(void **) &dy_origcons,TRUE) ;
      FREE(dy_origcons) ;
      dy_origcons = NULL ; }
    dy_freelclsystem(orig_lp,TRUE) ;
    if (dy_lp != NULL) (void) dy_swapobjs(dyDONE) ;
    if (dy_sys != NULL)
    { consys_free(dy_sys) ;
      dy_sys = NULL ; }
    FREE_AND_CLEAR(dy_actvars)
    FREE_AND_CLEAR(dy_actcons)
    FREE_AND_CLEAR(dy_basis)
    FREE_AND_CLEAR(dy_var2basis)
    FREE_AND_CLEAR(dy_status)
    FREE_AND_CLEAR(dy_x)
    FREE_AND_CLEAR(dy_xbasic)
    FREE_AND_CLEAR(dy_y)
    FREE_AND_CLEAR(dy_frame)
    FREE_AND_CLEAR(dy_gamma)
    FREE_AND_CLEAR(dy_cbar)
    FREE_AND_CLEAR(dy_rho)
    FREE_AND_CLEAR(dy_brkout)
    FREE_AND_CLEAR(dy_degenset)
    FREE_AND_CLEAR(dy_ddegenset)
    dy_freepivrej() ;
    FREE_AND_CLEAR(dy_lp)
    FREE_AND_CLEAR(dy_tols)
    FREE_AND_CLEAR(dy_opts)
    clrflg(orig_lp->ctlopts,lpctlDYVALID) ;
    dy_retained = FALSE ; }
  else
  { setflg(orig_lp->ctlopts,lpctlDYVALID) ;
    dy_freelclsystem(orig_lp,FALSE) ;
    dy_retained = TRUE ;}

  return ; }

#undef FREE_AND_CLEAR



void dy_freesoln (lpprob_struct *lpprob)

/*
  Simple utility to free the solution structures hanging off an lp problem
  structure. Does not free the constraint system, or the lpprob_struct.

  Parameters:
    lpprob:	the lp problem
  
  Returns: undefined
*/

{ if (lpprob->basis != NULL)
  { if (lpprob->basis->el != NULL) FREE(lpprob->basis->el) ;
    FREE(lpprob->basis) ;
    lpprob->basis = NULL ; }
  if (lpprob->status != NULL)
  { FREE(lpprob->status) ;
    lpprob->status = NULL ; }
  if (lpprob->x != NULL)
  { FREE(lpprob->x) ;
    lpprob->x = NULL ; }
  if (lpprob->y != NULL)
  { FREE(lpprob->y) ;
    lpprob->y = NULL ; }
  if (lpprob->actvars != NULL)
  { FREE(lpprob->actvars) ;
    clrflg(lpprob->ctlopts,lpctlACTVARSIN|lpctlACTVARSOUT) ;
    lpprob->actvars = NULL ; }
  lpprob->colsze = 0 ;
  lpprob->rowsze = 0 ;

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

# ifdef PARANOIA
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

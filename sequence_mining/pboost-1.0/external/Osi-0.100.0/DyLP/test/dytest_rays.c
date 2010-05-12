/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2008 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This file contains routines to test dylp's ray routines: dy_primalRays and
  dy_dualRays.
*/

#include "dylp.h"

extern ioid dy_logchn ;
extern bool dy_gtxecho ;



int dytest_primalRays (int *p_numRays,
		       lpprob_struct *main_lp, lptols_struct *main_lptols,
		       lpopts_struct *main_lpopts)

/*
  This routine checks the primal rays returned by dy_primalRays. For a ray r
  and a constraint ax <= b, the test is that dot(a,r) <= 0. For a constraint
  ax >= b, the test is dot(a,r) >= 0.

  It's up to the calling routine to determine if the number of rays is as
  expected. In particular, it's not an error if dy_primalRays returns fewer
  rays than requested.  If dy_primalRays returns zero rays, this is treated
  as the degenerate case of `all rays pass' and the routine will return
  TRUE.

  Parameters:
    p_numRays:	 (i) the number of rays to request from dy_primalRays
		 (o) the number of rays returned by dy_primalRays
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if all rays returned tested as valid rays, error count otherwise.
*/

{ int m,n,i,k ;
  consys_struct *sys ;

  double **rays ;
  int reqRays,rcvRays ;
  double *rayk ;
  double aidotrayk ;
  bool error ;
  int errcnt ;

  char *rtnnme = "dytest_primalRays" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  n = sys->varcnt ;
  m = sys->concnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.rays >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: checking primal rays using %s (%d x %d).\n",
		rtnnme,sys->nme,m,n) ; }
# endif
/*
  Ask for the requested number of rays.
*/
  rays = NULL ;
  reqRays = *p_numRays ;
  rcvRays = reqRays ;
  *p_numRays = 0 ;
  if (dy_primalRays(main_lp,&rcvRays,&rays) == FALSE)
  { errmsg(955,rtnnme,sys->nme,"primal") ;
    if (rays != NULL)
    { for (k = 0 ; k < rcvRays ; k++)
      { if (rays[k] != NULL) FREE(rays[k]) ; }
      FREE(rays) ; }
    return (1) ; }
  *p_numRays = rcvRays ;
# ifndef DYLP_NDEBUG
  if (main_lpopts->print.rays >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    requested %d rays, received %d.",
		reqRays,rcvRays) ; }
# endif
/*
  Now test each ray. Check first that we actually have a nonzero ray, then
  check dot(a,r) <= 0 (for ax >= b, dot(a,r) >= 0).
*/
  errcnt = 0 ;
  for (k = 0 ; k < rcvRays ; k++)
  { rayk = rays[k] ;
    if (rayk == NULL)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  ERROR: ray %d is NULL.",k) ;
      errcnt++ ;
      continue ; }
    aidotrayk = exvec_1norm(rayk,n) ;
    if (fabs(aidotrayk) <= 0.0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  ERROR: ray %d is zero.",k) ;
      FREE(rayk) ;
      rays[k] = NULL ;
      continue ; }
    for (i = 1 ; i <= m ; i++)
    { aidotrayk = consys_dotrow(sys,i,rayk) ;
      error = FALSE ;
      if (sys->ctyp[i] == contypGE)
      { if (aidotrayk < -main_lptols->zero)
	{ error = TRUE ; } }
      else
      { if (aidotrayk > main_lptols->zero)
	{ error = TRUE ; } }
      if (error == TRUE)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: a<%d> dot ray<%d> = %g ; should be %s 0.",
		    i,k,aidotrayk,((sys->ctyp[i] == contypGE)?">=":"<=")) ;
	errcnt++ ; } }
    FREE(rayk) ;
    rays[k] = NULL ; }
/*
  We're done. Clean up and go home.
*/
  FREE(rays) ;
  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	  "\n%s: found %d errors in %d rays testing Ar <= 0.\n",
	   rtnnme,errcnt,rcvRays) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: pass Ar <= 0.\n",rtnnme) ; }

  return (errcnt) ; }



int dytest_dualRays (int *p_numRays,
		     lpprob_struct *main_lp, lptols_struct *main_lptols,
		     lpopts_struct *main_lpopts)

/*
  This routine checks the dual rays returned by dy_dualRays. For a ray r
  and a dual constraint ya >= c, the test is that dot(r,a) >= 0.

  It's up to the calling routine to determine if the number of rays is as
  expected. In particular, it's not an error if dy_dualRays returns fewer
  rays than requested.  If dy_dualRays returns zero rays, this is treated
  as the degenerate case of `all rays pass' and the routine will return
  TRUE.

  Parameters:
    p_numRays:	 (i) the number of rays to request from dy_dualRays
		 (o) the number of rays returned by dy_dualRays
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if all rays returned tested as valid rays, error count otherwise.
*/

{ int m,n,j,k ;
  consys_struct *sys ;

  double **rays ;
  int reqRays,rcvRays ;
  double *rayk ;
  double aidotrayk ;
  bool error ;
  int errcnt ;

  char *rtnnme = "dytest_dualRays" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  n = sys->varcnt ;
  m = sys->concnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.rays >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: checking dual rays using %s (%d x %d).\n",
		rtnnme,sys->nme,m,n) ; }
# endif
/*
  Ask for the requested number of rays.
*/
  rays = NULL ;
  reqRays = *p_numRays ;
  rcvRays = reqRays ;
  *p_numRays = 0 ;
  if (dy_dualRays(main_lp,&rcvRays,&rays) == FALSE)
  { errmsg(955,rtnnme,sys->nme,"dual") ;
    if (rays != NULL)
    { for (k = 0 ; k < rcvRays ; k++)
      { if (rays[k] != NULL) FREE(rays[k]) ; }
      FREE(rays) ; }
    return (1) ; }
  *p_numRays = rcvRays ;
# ifndef DYLP_NDEBUG
  if (main_lpopts->print.rays >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    requested %d rays, received %d.",
		reqRays,rcvRays) ; }
# endif
/*
  Now test each ray. Check first that we actually have a nonzero ray, then
  check dot(r,a) >= 0. Note that since the dual constraints we're testing are
  implicit, there's no need to consider the constraint type.
*/
  errcnt = 0 ;
  for (k = 0 ; k < rcvRays ; k++)
  { rayk = rays[k] ;
    if (rayk == NULL)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  ERROR: ray %d is NULL.",k) ;
      errcnt++ ;
      continue ; }
    aidotrayk = exvec_1norm(rayk,m) ;
    if (fabs(aidotrayk) <= 0.0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  ERROR: ray %d is zero.",k) ;
      FREE(rayk) ;
      rays[k] = NULL ;
      continue ; }
    for (j = 1 ; j <= n ; j++)
    { aidotrayk = consys_dotcol(sys,j,rayk) ;
      error = FALSE ;
      if (aidotrayk < -main_lptols->zero)
      { error = TRUE ; }
      if (error == TRUE)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: ray<%d> dot a<%d> = %g ; should be >= 0.",
		    k,j,aidotrayk) ;
	errcnt++ ; } }
    FREE(rayk) ;
    rays[k] = NULL ; }
/*
  We're done. Clean up and go home.
*/
  FREE(rays) ;
  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	  "\n%s: found %d errors in %d rays testing rA >= 0.\n",
	   rtnnme,errcnt,rcvRays) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: pass rA >= 0.\n",rtnnme) ; }

  return (errcnt) ; }


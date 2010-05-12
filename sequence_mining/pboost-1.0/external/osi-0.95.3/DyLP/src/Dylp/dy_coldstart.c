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
  This file contains routines to populate the active constraint system and
  select an initial basis for an lp. Both of these activities are discussed in
  more detail below. For here, it suffices to say that dylp offers a great deal
  of flexibility in establishing the initial constraint system and selecting
  the initial basis.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_coldstart.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_coldstart.c 94 2006-06-29 23:06:51Z lou $" ;



/*
  Cold start routines.

  dy_coldstart and subsidiary routines select a set of constraints to populate
  the initial active system. There are a number of options. But first a bit of
  general explanation.

  Ideally, we'd like to be able to guess the constraints that will be tight at
  optimum and load exactly those constraints. If we can also guess the right
  basic variables, we have the answer. It's not that easy, of course. See the
  comments in the second half of the file for an explanation of how dylp
  chooses the initial set of basic variables.

  Back to the subject at hand. Recall that we're dealing with equalities and
  <= constraints. Also recall that we're minimising. Then the normals a<i> to
  the inequality constraints point out of the feasible region, and at optimum
  the normal c to the objective will point into the feasible region. Hence
  ideal `alignment' with the objective occurs when angle between a<i> and c is
  180 degrees. When the angle is 0 degrees, the constraint is forming the far
  side of the polytope.

  For convenience, call the inequalities with 180 <= angle(a<i>,c) < 90 the
  `near' group, those with angle(a<i>,c) = 90 the `perp' group, and those
  with 90 < angle(a<i>,c) <= 0 the `far' group.  The implementation creates a
  sorted list of inequalities, in nonincreasing order using angle(a<i>,c),
  and marks the boundaries between the near, perp, and far groups.

  To populate the active system, equalities are loaded first. All equalities
  are loaded, always.  Empty constraints are never loaded. For the
  inequalities, it's possible to load one or two specified angular ranges,
  with a specified sampling rate. This seems to allow as much flexibility as
  is useful.

  NOTE: The routines that construct the initial basis use knowledge of the
	load order, and knowledge of how consys augments columns, to scan the
	columns in the proper order when selecting architectural variables
	for basis positions. If you change one, you must change both. See
	ib_archvselect.

  Note that specifying a range of 180 - 0 degrees with sampling rate 1.0
  (i.e., activate all the original constraints) is >not< quite the same as
  the fullsys option. 100% activation still allows initial variable
  deactivation (hence multiple variable activation/deactivation phases) and
  final variable and constraint deactivation; fullsys suppresses these.

  Recall that angle(a<i>,c) = (180/pi)*arccos(dot(a<i>,c)/(||a<i>||||c||))
*/


/*
  Utility structures

  angle_struct

  Field		Definition
  -----		----------
  ndx		constraint index
  angle		angle between a<ndx> and c

  ineq_struct

  Field		Definition
  -----		----------
  cnt		Number of inequalities
  perp		index in angles of first inequality with
		angle(a<i>,c) = 90 degrees
  far		index in angles of first inequality with
		angle(a<i>,c) < 90 degrees
  angles	inequalities, sorted by angle (indexed from 0)
*/

typedef struct { int ndx ;
		 double angle ; } angle_struct ;

typedef struct { int cnt ;
		 int perp ;
		 int far ;
		 angle_struct *angles ; } ineq_struct ;


static int near_perp_far (const void *elem1, const void *elem2)

{ const angle_struct *c1,*c2 ;
  c1 = (const angle_struct *) elem1 ;
  c2 = (const angle_struct *) elem2 ;

  if (c1->angle > c2->angle)
    return (-1) ;
  else
  if (c1->angle < c2->angle)
    return (1) ;
  else
    return (0) ; }


static bool cold_sortcons (consys_struct *orig_sys,
			   int **p_eqs, ineq_struct **p_ineqs)

/*
  This routine separates the constraints into equalities and inequalities,
  dropping empty constraints. Depending on options, various additional work is
  performed:
    * If fullsys is not specified, the inequalities are sorted according to
      their angle from the objective function.
    * If full statistics are enabled, constraint angles are loaded into the
      statistics structure and a histogram is compiled.
    * If debug printing is enabled, the requested information is dumped to the
      log.

  Parameters:
    orig_sys:	The original constraint system
    p_eqs:	(i) empty array of int; assumed to be sufficiently large;
		    allocated if NULL
		(o) eqs[0] set to number of equalities
		    eqs[1 .. eqs[0]] set to indices of equality constraints
    p_ineqs:	(i) empty ineq_struct; assumed to be sufficiently large;
		    allocated if NULL
		(o) filled in with inequality information as described in
		    comments at head of file

  Return value: TRUE if all goes as planned, FALSE on error.
*/

{ int i,ndx,m,n,eqcnt,ineqcnt,nearcnt,perpcnt,farcnt ;
  double cnorm,ainorm,aidotc,pi180,anglei ;
  double *c ;
  int *eqs ;
  ineq_struct *ineqs ;
  angle_struct *angles ;
  contyp_enum *ctyp ;

  bool retval,need_angles ;

  pkvec_struct *aj ;

  const char *rtnnme = "cold_sortcons" ;

# ifdef DYLP_STATISTICS
  int k ;
# endif

# ifdef PARANOIA
  if (orig_sys == NULL)
  { errmsg(2,rtnnme,"orig_sys") ;
    return (FALSE) ; }
  if (p_eqs == NULL)
  { errmsg(2,rtnnme,"&eqs") ;
    return (FALSE) ; }
  if (p_ineqs == NULL)
  { errmsg(2,rtnnme,"&ineqs") ;
    return (FALSE) ; }
# endif

/*
  Allocate the equality and inequality data structures, if the client didn't
  supply them.
*/
  if (*p_eqs == NULL)
  { eqs = (int *) MALLOC((orig_sys->concnt+1)*sizeof(int)) ; }
  else
  { eqs = *p_eqs ; }
  if (*p_ineqs == NULL)
  { ineqs = CALLOC(1,sizeof(ineq_struct)) ;
    ineqs->angles =
	(angle_struct *) MALLOC(orig_sys->concnt*sizeof(angle_struct)) ; }
  else
  { ineqs = *p_ineqs ; }
# ifdef PARANOIA
  if (ineqs->angles == NULL)
  { errmsg(2,rtnnme,"angle array") ;
    return (FALSE) ; }
# endif
  angles = ineqs->angles ;

  m = orig_sys->concnt ;
  n = orig_sys->varcnt ;
  retval = TRUE ;
/*
  Scan orig_sys, sorting the constraints into eqs and ineqs->angles and
  discarding empty constraints.
*/
  aj = pkvec_new(0) ;
  eqcnt = 0 ;
  ineqcnt = 0 ;
  ctyp = orig_sys->ctyp ;

  for (i = 1 ; i <= m ; i++)
  { if (consys_getrow_pk(orig_sys,i,&aj) == FALSE)
    { errmsg(122,rtnnme,orig_sys->nme,"row",
	     consys_nme(orig_sys,'c',i,TRUE,NULL),i) ;
      retval = FALSE ;
      break ; }
    if (aj->cnt != 0)
    { if (ctyp[i] == contypEQ)
      { eqs[++eqcnt] = i ; }
      else
      { angles[ineqcnt++].ndx = i ; } } }

  eqs[0] = eqcnt ;
  ineqs->cnt = ineqcnt ;
  if (aj != NULL) pkvec_free(aj) ;
  if (retval == FALSE)
  { if (*p_eqs == NULL && eqs != NULL) FREE(eqs) ;
    if (*p_ineqs == NULL && ineqs != NULL)
    { if (ineqs->angles != NULL) FREE(angles) ;
      FREE(ineqs) ; }
    return (FALSE) ; }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    found %d equalities, %d inequalities",
	        eqcnt,ineqcnt) ;
    if (eqcnt+ineqcnt < m)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", discarded %d empty constraints",
		  (m-eqcnt-ineqcnt)) ; }
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
# endif
/*
  Now, how much more work do we need to do? Decide if we need to calculate
  angles.
*/
  need_angles = TRUE ;
  if (dy_opts->fullsys == TRUE) need_angles = FALSE ;
# ifdef DYLP_STATISTICS
  need_angles = TRUE ;
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 2) need_angles = TRUE ;
# endif
/*
  If we need the angles, get down to it. Calculate
      angle(a<i>,c) = (180/pi)*arccos(dot(a<i>,c)/(||a<i>||||c||))
  for each inequality, then sort the list of inequalities.
*/
  if (need_angles == TRUE)
  { c = orig_sys->obj ;
    nearcnt = 0 ;
    perpcnt = 0 ;
    farcnt = 0 ;
    cnorm = exvec_2norm(c,n) ;
    pi180 = 180/acos(-1.0) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t||c|| = %.4f",cnorm) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n\tConstraint\t\t||a||\t(c/||c||).(a/||a||)\tangle") ; }
#   endif

    for (ndx = 0 ; ndx < ineqcnt ; ndx++)
    { i = angles[ndx].ndx ;
      aidotc = consys_dotrow(orig_sys,i,c) ;
      if (withintol(aidotc,0,dy_tols->zero))
      { 
#	ifndef DYLP_NDEBUG
	ainorm = consys_2normrow(orig_sys,i) ;
#	endif
	anglei = 90 ; }
      else
      { ainorm = consys_2normrow(orig_sys,i) ;
	anglei = pi180*acos(aidotc/(ainorm*cnorm)) ; }
      angles[ndx].angle = anglei ;
      if (anglei > 90)
      { nearcnt++ ; }
      else
      if (anglei < 90)
      { farcnt++ ; }
      else
      { perpcnt++ ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t%-10s (%3d) %12.4f %18.10f%15.6f",
		    consys_nme(orig_sys,'c',i,FALSE,NULL),i,
		    ainorm,aidotc/(ainorm*cnorm),anglei) ; }
#     endif
    }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    inequality partition %d near, %d perp, %d far.",
		  nearcnt,perpcnt,farcnt) ; }
#   endif
    ineqs->perp = nearcnt ;
    ineqs->far = nearcnt+perpcnt ;
    qsort(&angles[0],ineqcnt,sizeof(angle_struct),near_perp_far) ; }

# ifdef DYLP_STATISTICS
/*
  If we're collecting full statistics, record the angle for each inequality
  and calculate a histogram (36 5 degree bins and a dedicated 90 degree
  bin).  This could have been folded into the previous loop, but it's more
  clear as a separate loop. Bins below 90 degrees include their lower bound,
  bins above 90 degrees include their upper bound.
*/
  if (dy_stats != NULL)
  { for (ndx = 0 ; ndx < ineqcnt ; ndx++)
    { i = angles[ndx].ndx ;
      anglei = angles[ndx].angle ;
      dy_stats->cons.angle[i] = anglei ;
      if (anglei == 90)
      { dy_stats->angle.hist[DYSTATS_HISTBINS-1]++ ; }
      else
      { if (anglei < 90)
	{ k = (int) floor(anglei/5) ; }
	else
	{ k = (int) ceil(anglei/5)-1 ; }
	k = minn(k,(DYSTATS_HISTBINS-1)) ;
	dy_stats->angle.hist[k]++ ; }
      if (anglei > dy_stats->angle.max)
      { dy_stats->angle.max = (float) anglei ; }
      if (anglei < dy_stats->angle.min)
      { dy_stats->angle.min = (float) anglei ; } } }
# endif

/*
  Set up the return values and we're done.
*/
  if (*p_eqs == NULL) *p_eqs = eqs ;
  if (*p_ineqs == NULL) *p_ineqs = ineqs ;

  return (TRUE) ; }



static bool cold_createdysys (consys_struct *orig_sys, int eqcnt, int ineqcnt)
/*
  This routine creates the dy_sys constraint system and attaches the
  translation vectors that are used to move between the original (orig_sys)
  and active (dy_sys) constraint systems.
  
  While unrelated to creating dy_sys, it turns out that this routine is a
  convenient place for another pre-loading activity, grooming the bounds on
  variables in orig_sys and assigning an initial status to the inactive
  variables. Where |ub<j> - lb<j>| < dy_tols->pfeas, the bounds are forced to
  equality and the variable is marked as fixed (NBFX). (This sort of thing
  can happen when a client program is tweaking the bounds on variables.) This
  may or may not be visible to the client --- depends on whether dylp has
  made a local copy of the constraint system.  Variables which are not fixed
  are arbitrarily set to the bound that best matches their objective
  coefficient.

  What's a good estimate of maximum size? If the fullsys option is specified,
  just use the original constraint system size. Otherwise, take the attitude
  that some (user) specified fraction of the original constraints and
  variables (active.cons, active.vars) could be active at any one time while
  we're solving the LP. This is purely for efficiency -- if we're undersized,
  the system will automatically expand.

  Given that all equalities are loaded, initcons.frac is interpreted as the
  initial fraction of inequalities that should be loaded. Again, for
  efficiency active.cons should be greater than initcons.frac. dy_sys is
  created with logicals enabled.

  Parameters:
    orig_sys:	the original constraint system
    eqcnt:	the number of equalities in orig_sys
    ineqcnt:	the number of inequalities in orig_sys

  Returns: TRUE if the system is created without error, FALSE otherwise.
*/

{ int j,m_sze,n_sze,flippable ;
  double vlbj,vubj ;
  double *vlb,*vub,*obj ;
  bool infeas ;
  char nmebuf[50] ;

  flags parts = CONSYS_OBJ|CONSYS_VUB|CONSYS_VLB|CONSYS_RHS|CONSYS_RHSLOW|
		CONSYS_VTYP|CONSYS_CTYP,
	opts = CONSYS_LVARS|CONSYS_WRNATT ;

  const char *rtnnme = "cold_createdysys" ;

/*
  Settle the appropriate size for dy_sys, and create the consys structure.
  The variable capacity is set to accommodate logicals.
*/
  dyio_outfxd(nmebuf,-((int) (sizeof(nmebuf)-1)),
	      'l',"%s[actv]",orig_sys->nme) ;
  if (dy_opts->fullsys == TRUE)
  { m_sze = orig_sys->concnt ;
    n_sze = orig_sys->archvcnt+m_sze ; }
  else
  { m_sze = eqcnt+((int) orig_sys->concnt*dy_opts->active.cons) ;
    n_sze = (int) (orig_sys->archvcnt*dy_opts->active.vars+m_sze) ; }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    creating constraint system %s (%d x %d)",
	        nmebuf,m_sze,n_sze) ; }
# endif
  dy_sys = consys_create(nmebuf,parts,opts,m_sze,n_sze,dy_tols->inf) ;
  if (dy_sys == NULL)
  { errmsg(152,rtnnme,nmebuf) ;
    return (FALSE) ; }
/*
  Hang a set of translation vectors onto each system: origcons and origvars
  on orig_sys, and actcons and actvars on dy_sys. dy_origvars is cleared to
  0 as it's attached, and this is important to indicate that the original
  variables have no predefined status.

  Note that we won't leak memory here, even if there's an error. Any attached
  vectors will be deleted when dy_sys is deleted.
*/
  dy_actvars = NULL ;
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(int),(void **) &dy_actvars) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"active -> original variable map") ;
    return (FALSE) ; }
  dy_actcons = NULL ;
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(int),(void **) &dy_actcons) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"active -> original constraint map") ;
    return (FALSE) ; }
  dy_origvars = NULL ;
  if (consys_attach(orig_sys,CONSYS_ROW,
		    sizeof(int),(void **) &dy_origvars) == FALSE)
  { errmsg(100,rtnnme,orig_sys->nme,"original -> active variable map") ;
    return (FALSE) ; }
  dy_origcons = NULL ;
  if (consys_attach(orig_sys,CONSYS_COL,
		    sizeof(int),(void **) &dy_origcons) == FALSE)
  { errmsg(100,rtnnme,orig_sys->nme,"original -> active constraint map") ;
    return (FALSE) ; }
/*
  Assign a status to all variables in orig_sys. There are lots of reasons
  for doing this right off the top:
    * It's not uncommon for people to fix variables in an MPS model. This takes
      them out of consideration right from the start.
    * dylp's intended use is in branch-and-cut codes, which tend to update
      bounds on a regular basis. This can result in bounds which are not
      precisely equal, but are within the primal feasibility tolerance of one
      another. This is a bad situation all around, and needs to be fixed on
      input. Even worse, it's possible we'll be handed bounds which are prima
      facie infeasible --- lower and upper bounds are crossed. This is really
      bad, and dylp does not react well. In this case we need to detect
      infeasibility and report it back to dylp().
    * We can use origvars == 0 as a paranoid check from here on out.
    * If there are no variables with upper and lower bounds (`flippable', in
      dual multipivot) then we might as well turn multipivot off.
  This activity doesn't particularly belong with the creation of dy_sys, but
  this is a convenient place to put it and the information in dy_origvars needs
  to be valid when we begin loading constraints.
*/
  vlb = orig_sys->vlb ;
  vub = orig_sys->vub ;
  obj = orig_sys->obj ;
  flippable = 0 ;
  for (j = 1 ; j <= orig_sys->varcnt ; j++)
  { vlbj = vlb[j] ;
    vubj = vub[j] ;
    if (vlbj > -dy_tols->inf && vubj < dy_tols->inf)
    { if (atbnd(vlbj,vubj) && vlbj != vubj)
      { 
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.setup >= 3)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tForcing equal bound %g for %s (%d)",
		      (vlbj+vubj)/2,consys_nme(orig_sys,'v',j,0,0),j) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\t  original lb = %g, ub = %g, diff = %g, tol = %g",
		      vlbj,vubj,vubj-vlbj,dy_tols->pfeas) ; }
#	endif
	vlb[j] = (vlbj+vubj)/2 ;
	vub[j] = vlb[j] ; }
      if (vlbj == vubj)
      { dy_origvars[j] = -((int) vstatNBFX) ; }
      else
      { flippable++ ;
	if (obj[j] < 0)
	{ dy_origvars[j] = -((int) vstatNBUB) ; }
	else
	{ dy_origvars[j] = -((int) vstatNBLB) ; } } }
    else
    if (vlbj > -dy_tols->inf)
    { dy_origvars[j] = -((int) vstatNBLB) ; }
    else
    if (vubj < dy_tols->inf)
    { dy_origvars[j] = -((int) vstatNBUB) ; }
    else
    { dy_origvars[j] = -((int) vstatNBFR) ; }
    if (vub[j] < vlb[j])
    { dy_lp->lpret = lpINFEAS ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tTrivial infeasibility for %s (%d), lb = %g > ub = %g.",
		    consys_nme(orig_sys,'v',j,0,0),j,vlb[j],vub[j]) ; }
#     endif
    } }
/*
  Disable dual multipivoting? Fairly arbitrarily, give 25% flippable variables
  as the criterion.
*/
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->dmulti.flippable = flippable ;
# endif
  vubj = ((double) flippable)/orig_sys->varcnt ;
  if (vubj < .25 && dy_opts->dpsel.flex == TRUE)
  { dy_opts->dpsel.flex = FALSE ;
    dy_opts->dpsel.strat = 0 ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 2 || dy_opts->print.dual >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	    "\n    %d (%g%%) flippable variables; disabling dual multipivot.",
	    flippable,vubj) ; }
#   endif
  }

  return (TRUE) ; }



static bool cold_loadfull (consys_struct *orig_sys,
			    int *eqs, ineq_struct *ineqs)
/*
  This routine loads the full original constraint system into dy_sys. It's
  just easier to handle this separate from the more complicated business of
  loading a partial system.

  Parameters:
    orig_sys:	the original constraint system
    eqs:	indices of equality constraints; eq[0] contains the number of
		valid entries
    ineqs:	information on inequalities; only the number of inequalities
		and their indices are used here.

  Returns: TRUE if all constraints are loaded successfully, FALSE otherwise.
*/

{ int j,ndx,eqcnt,ineqcnt ;
  bool retval ;

  angle_struct *angles ;

  const char *rtnnme = "cold_loadfull" ;

  eqcnt = eqs[0] ;
  ineqcnt = ineqs->cnt ;
  retval = TRUE ;
/*
  Load any equalities first.
*/
  if (eqcnt > 0)
  { for (ndx = 1 ; ndx <= eqcnt ; ndx++)
    { j = eqs[ndx] ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho, "\n    activating %s %s (%d) ...",
		    consys_prtcontyp(orig_sys->ctyp[j]),
		    consys_nme(orig_sys,'c',j,FALSE,NULL),j) ; }
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL) dy_stats->cons.init[j] = TRUE ;
#     endif
      if (dy_loadcon(orig_sys,j,TRUE,NULL) == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "activate","constraint",
	       consys_nme(orig_sys,'c',j,TRUE,NULL),j) ;
	retval = FALSE ;
	break ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    transferred %d equalities ...",eqcnt) ; }
#   endif
    if (retval == FALSE)
    { return (FALSE) ; } }
/*
  And the inequalities.
*/
  if (ineqcnt > 0)
  { angles = ineqs->angles ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    transferred %d inequalities ...",ineqcnt) ; }
#   endif
    for (ndx = 0 ; ndx < ineqcnt ; ndx++)
    { j = angles[ndx].ndx ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho, "\n    activating %s %s (%d) ...",
		    consys_prtcontyp(orig_sys->ctyp[j]),
		    consys_nme(orig_sys,'c',j,FALSE,NULL),j) ; }
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL) dy_stats->cons.init[j] = TRUE ;
#     endif
      if (dy_loadcon(orig_sys,j,TRUE,NULL) == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "activate","constraint",
	       consys_nme(orig_sys,'c',j,TRUE,NULL),j) ;
	retval = FALSE ;
	break ; } } }

  return (retval) ; }



static bool cold_loadpartial (consys_struct *orig_sys,
			      int *eqs, ineq_struct *ineqs)
/*
  This routine loads some fraction of the original constraint system into
  dy_sys. As explained in the comments at the head of the file, the
  inequalities are sorted by angle(a<i>,c). The client can specify one or two
  angular intervals, and the upper and lower bounds of each interval can be
  closed or open. Both intervals are sampled at the specified fractional rate.

  Parameters:
    orig_sys:	the original constraint system
    eqs:	indices of equality constraints; eq[0] contains the number of
		valid entries
    ineqs:	information on inequalities; only the number of inequalities
		and their indices are used here.

  Returns: TRUE if all constraints are loaded successfully, FALSE otherwise.
*/

{ int j,ndx,eqcnt,ineqcnt,ineq_actvcnt ;
  bool retval ;
  angle_struct *angles ;

  int intcnt,intndx ;
  double albs[2],aubs[2],alb,aub,dblndx,incr ;

  const char *rtnnme = "cold_loadpartial" ;

  eqcnt = eqs[0] ;
  ineqcnt = ineqs->cnt ;
  retval = TRUE ;
/*
  Load any equalities first.
*/
  if (eqcnt > 0)
  { for (ndx = 1 ; ndx <= eqcnt ; ndx++)
    { j = eqs[ndx] ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho, "\n    activating %s %s (%d) ...",
		    consys_prtcontyp(orig_sys->ctyp[j]),
		    consys_nme(orig_sys,'c',j,FALSE,NULL),j) ; }
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL) dy_stats->cons.init[j] = TRUE ;
#     endif
      if (dy_loadcon(orig_sys,j,TRUE,NULL) == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "activate","constraint",
	       consys_nme(orig_sys,'c',j,TRUE,NULL),j) ;
	retval = FALSE ;
	break ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.setup >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    transferred %d equalities ...",eqcnt) ; }
#   endif
    if (retval == FALSE)
    { return (FALSE) ; } }
/*
  If there are no inequalities, we're done.
*/
  if (ineqcnt == 0) return (TRUE) ;
/*
  Establish the number of intervals, angular boundaries, and fractional
  increment. Recall that constraints are sorted nonincreasing from 180
  degrees to 0 degrees. To approximate an open interval, we pull in the
  specified bound by the zero tolerance; for a closed interval, we push it
  out.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    inequalities: sampling %.2f, %c %.4f %.4f %c",
	        dy_opts->initcons.frac,
	        (dy_opts->initcons.i1uopen == TRUE)?'(':'[',
	        dy_opts->initcons.i1u,dy_opts->initcons.i1l,
	        (dy_opts->initcons.i1lopen == TRUE)?')':']') ;
    if (dy_opts->initcons.i2valid == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", %c %.4f %.4f %c",
		  (dy_opts->initcons.i2uopen == TRUE)?'(':'[',
		  dy_opts->initcons.i2u,dy_opts->initcons.i2l,
		  (dy_opts->initcons.i2lopen == TRUE)?')':']') ; } }
# endif
  incr = 1/dy_opts->initcons.frac ;
  if (dy_opts->initcons.i1lopen == TRUE)
  { albs[0] = dy_opts->initcons.i1l+dy_tols->zero ; }
  else
  { albs[0] = dy_opts->initcons.i1l-dy_tols->zero ; }
  if (dy_opts->initcons.i1uopen == TRUE)
  { aubs[0] = dy_opts->initcons.i1u-dy_tols->zero ; }
  else
  { aubs[0] = dy_opts->initcons.i1u+dy_tols->zero ; }
  if (dy_opts->initcons.i2valid == TRUE)
  { if (dy_opts->initcons.i2lopen == TRUE)
    { albs[1] = dy_opts->initcons.i2l+dy_tols->zero ; }
    else
    { albs[1] = dy_opts->initcons.i2l-dy_tols->zero ; }
    if (dy_opts->initcons.i2uopen == TRUE)
    { aubs[1] = dy_opts->initcons.i2u-dy_tols->zero ; }
    else
    { aubs[1] = dy_opts->initcons.i2u+dy_tols->zero ; }
    intcnt = 2 ; }
  else
  { intcnt = 1 ; }
  angles = ineqs->angles ;
  ineq_actvcnt = 0 ;
/*
  The outer loop is just to avoid replicating code for each interval.
*/
  ndx = 0 ;
  for (intndx = 0 ; intndx < intcnt ; intndx++)
  { alb = albs[intndx] ;
    aub = aubs[intndx] ;
/*
  Scan to the upper boundary of the interval. Consider the possibility that
  we might run off the end of angles. An easy and all too common example: the
  intervals exclude 90 degrees, and all inequalities are at 90 degrees to the
  objective. This is an artifact of a modelling style in which the objective
  function is a single variable z, constrained to be equal to the real
  objective cx. Then z appears only in the single equality z - cx = 0.
*/
    while (ndx < ineqcnt && angles[ndx].angle > aub) ndx++ ;
    dblndx = ndx ;
/*
  Now the actual loading loop. We load the constraint, then move ndx to the
  next value according to the sampling fraction.
*/
    while (ndx < ineqcnt && angles[ndx].angle > alb) 
    { j = angles[ndx].ndx ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.setup >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    activating %s %s (%d) %g off 90 ...",
		    consys_prtcontyp(orig_sys->ctyp[j]),
		    consys_nme(orig_sys,'c',j,FALSE,NULL),j,
		    (angles[ndx].angle-90)) ; }
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL) dy_stats->cons.init[j] = TRUE ;
#     endif
      if (dy_loadcon(orig_sys,j,TRUE,NULL) == FALSE)
      { errmsg(430,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       "activate","constraint",
	       consys_nme(orig_sys,'c',j,TRUE,NULL),j) ;
	retval = FALSE ;
	break ; }
      ineq_actvcnt++ ;
      dblndx += incr ;
      ndx = (int) dblndx ; } }
/*
  That's it. A little information and we're out of here.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    transferred %d inequalities ...",ineq_actvcnt) ; }
# endif

  return (retval) ; }



dyret_enum dy_coldstart (consys_struct *orig_sys)

/*
  This routine is responsible for setting up the lp problem that'll be solved
  by dylp. It creates the lp and constraint system structures (dy_lp and
  dy_sys, respectively) and oversees the load of the initial constraint system.

  Parameters:
    orig_sys:	The original constraint system

  Returns: dyrOK if the setup completes without error, dyrFATAL otherwise.
*/

{ int j,n,eqcnt,ineqcnt ;
  double *vlb,*vub,*obj ;
  double vlbj,vubj,objj ; 
  flags statj ;
  int *eqs ;
  ineq_struct *ineqs ;

# ifndef DYLP_NDEBUG
  int nbfxcnt = 0 ;
# endif

  bool retval ;
  
  const char *rtnnme = "dy_coldstart" ;

/*
  To get started, sort the constraints into equalities and inequalities. If
  we're loading a partial active system, the inequalities will also be sorted
  by angle from the objective function.
*/
  eqs = NULL ;
  ineqs = NULL ;
  retval = cold_sortcons(orig_sys,&eqs,&ineqs) ;
  if (retval == FALSE)
  { errmsg(312,rtnnme,
	  orig_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
    if (eqs != NULL) FREE(eqs) ;
    if (ineqs != NULL)
    { if (ineqs->angles != NULL) FREE(ineqs->angles) ;
      FREE(ineqs) ; }
    return (dyrFATAL) ; }
  eqcnt = eqs[0] ;
  ineqcnt = ineqs->cnt ;
/*
  Next, create dy_sys and attach the translation vectors that will allow us to
  move between the original and active systems.
*/
  retval = cold_createdysys(orig_sys,eqcnt,ineqcnt) ;
  if (retval == FALSE)
  { if (eqs != NULL) FREE(eqs) ;
    if (ineqs != NULL)
    { if (ineqs->angles != NULL) FREE(ineqs->angles) ;
      FREE(ineqs) ; }
    return (dyrFATAL) ; }
/*
  Transfer the required constraints. Once this is done,  we're finished with
  the lists of equalities and inequalities.
*/
  if (dy_opts->fullsys == TRUE)
  { retval = cold_loadfull(orig_sys,eqs,ineqs) ; }
  else
  { retval = cold_loadpartial(orig_sys,eqs,ineqs) ; }
  if (eqs != NULL) FREE(eqs) ;
  if (ineqs != NULL)
  { if (ineqs->angles != NULL) FREE(ineqs->angles) ;
    FREE(ineqs) ; }
  if (retval == FALSE)
  { errmsg(313,rtnnme,
	    dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	    orig_sys->nme) ;
    return (dyrFATAL) ; }
/*
  Scan the variables in orig_sys once again and calculate the correction to
  the objective function. (loadcon handled any rhs adjustments for the
  constraints it loaded.)  Inactive free variables are assumed to have value
  0.
*/
  n = orig_sys->varcnt ;
  vlb = orig_sys->vlb ;
  vub = orig_sys->vub ;
  obj = orig_sys->obj ;
  dy_lp->inactzcorr = 0 ;
  for (j = 1 ; j <= n ; j++)
  { if (dy_origvars[j] < 0)
    { statj = (flags) (-dy_origvars[j]) ;
      vlbj = vlb[j] ;
      vubj = vub[j] ;
      objj = obj[j] ;
      switch (statj)
      { case vstatNBLB:
	case vstatNBFX:
	{ dy_lp->inactzcorr += objj*vlbj ;
#	  ifndef DYLP_NDEBUG
	  if (statj == vstatNBFX) nbfxcnt++ ;
#	  endif
	  break ; }
	case vstatNBUB:
	{ dy_lp->inactzcorr += objj*vubj ;
	  break ; } } } }
/*
  Paranoid checks and informational print statements. Apologies for abusing
  previously declared variables for the print loop.
*/
# ifdef PARANOIA
  if (dy_chkdysys(orig_sys) == FALSE) return (dyrFATAL) ;
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.setup >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  system %s has %d constraints, (%d+%d) variables",
	        dy_sys->nme,dy_sys->concnt,dy_sys->archvcnt,dy_sys->logvcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n  %d constraints, %d variables ",
	        orig_sys->concnt-dy_sys->concnt,
	        orig_sys->archvcnt-dy_sys->archvcnt) ;
    if (nbfxcnt > 0) dyio_outfmt(dy_logchn,dy_gtxecho," (%d fixed)",nbfxcnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		" remain inactive in system %s.",orig_sys->nme) ;
    if (dy_opts->print.setup >= 6)
    { vubj = 0 ;
      for (j = 1 ; j <= n ; j++)
      { if (dy_origvars[j] < 0)
	{ statj = (flags)(-dy_origvars[j]) ;
	  switch (statj)
	  { case vstatNBUB:
	    { vlbj = vub[j] ;
	      break ; }
	    case vstatNBLB:
	    { vlbj = vlb[j] ;
	      break ; }
	    case vstatNBFX:
	    case vstatNBFR:
	    { vlbj = 0 ;
	      break ; }
	    default:
	    { errmsg(433,rtnnme,dy_sys->nme,
		     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		     "inactive",consys_nme(orig_sys,'v',j,TRUE,NULL),
		     j,dy_prtvstat(statj)) ;
	      return (dyrFATAL) ; } }
	  if (vlbj != 0)
	  { if (vubj == 0)
	    { dyio_outfmt(dy_logchn,dy_gtxecho,
			  "\n\tinactive variables with nonzero values:") ; }
	    vubj++ ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) = %g, status %s",
		        consys_nme(orig_sys,'v',j,FALSE,NULL),j,vlbj,
		        dy_prtvstat(statj)) ; } } }
      if (vubj == 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tall inactive variables are zero.") ; } } }
# endif

  return (dyrOK) ; }




/*
  Routines to select an initial basis.

  dylp offers three types of starting basis:

    logical:    Basis positions occupied by inequalities are covered with
		slack variables.  Basis positions occupied by equalities are
		covered with artificial variables.  This is the standard
		starting basis.

    slack:	Basis positions occupied by inequalities are covered with
		slack variables. For basis positions occupied by equalities,
		an attempt is made to select architectural variables before
		falling back on artificial variables.

    architectural: For all basis positions, an attempt is made to cover the
		position with an architectural variable, before falling back on
		a slack or artificial. The selection process is performed in a
		way that will use an architectural to cover an equality before
		considering inequalities.
  
  As Bixby notes in an ORSA JOC article, crash is often taken to mean trying
  to guess the optimal basis. That's part of what's happening here, but not
  all. Architectural is chosen as the default basis type on the grounds that
  the default algorithms used in dy_coldstart and subsidiary routines are
  trying to choose constraints that will be tight at optimum. So it makes
  sense to build an initial basis with as few basic logicals as possible.
  But just as important, dy_crash and subsidiary routines try to select
  variables with an eye for numerical stability and ease of factorization.

  The algorithms used here to construct the slack and architectural basis
  styles had their roots in Bixby's article:

    Bixby, R., "Implementing the Simplex Method: The Initial Basis",
    ORSA J. on Computing, 4(3), Summer, 1992, pp. 267-284.
*/
/*
  Some utility declarations used in building the initial basis.
  
  The ibrank structure holds information useful in selecting architectural
  variable for inclusion in the initial basis. The function ib_archvcomp is
  used by qsort to sort the list of eligible candidates.

  ibrank_struct

  bndcnt and nonzero are used for the initial sort of the architecturals.

  If the matrix has been scaled, ajmax = 1 by definition for every row and
  column. But if the matrix hasn't been scaled, it's useful to know the
  maximum.

  Field		Definition
  -----		----------
  ndx		The index of the variable
  bndcnt	The number of bounds
  nonzero	The number of nonzero coefficients in the column
  ajmax		The largest value in the column
*/

typedef struct {
  int ndx ;
  int bndcnt ;
  int nonzero ;
  double ajmax ; } ibrank_struct ;


static int ib_archvcomp (const void *elem1, const void *elem2)
/*
  This function is called by qsort to compare architectural variables.
  The comparison is calculated using the following criteria:
    * Number of bounds: free variables (no finite bounds) are best, followed
      by variables with one finite bound, then two finite bounds
    * Number of coefficients in the column: Fewer coefficients are preferable,
      as this will tend to minimize the work during factoring.
*/


{ const ibrank_struct *v1,*v2 ;

  v1 = (const ibrank_struct *) elem1 ;
  v2 = (const ibrank_struct *) elem2 ;

  if (v1->bndcnt < v2->bndcnt)
  { return (-1) ; }
  else
  if (v1->bndcnt > v2->bndcnt)
  { return (1) ; }
  else
  if (v1->nonzero < v2->nonzero)
  { return (-1) ; }
  else
  if (v1->nonzero > v2->nonzero)
  { return (1) ; }
  else
  { return (0) ; } }



static bool ib_archvrank (int *p_cnt, ibrank_struct **p_archvars)

/*
  This routine sorts the architectural variables in preparation for
  constructing the initial basis. It scans the architectural variables,
  compiling information on the eligible candidates in archvars, and returns a
  sorted list.

  The ranking criteria are explained above in the comments for ib_archvcomp
  (the comparison function passed to qsort).

  Parameters:
    p_cnt:	(o) the number of eligible variables returned in p_archvars
    p_archvars:	(i) empty crashrank_struct vector of sufficient size; if NULL,
		    one will be allocated
		(o) list of eligible architectural variables, sorted in order
		    of preference for inclusion in the basis; may be NULL if
		    there are no candidates (extremely unlikely)

  Returns: TRUE if the sort concludes without error, FALSE otherwise
*/

{ int j,m,n,eligible ;
  double *vlb,*vub ;
  double vlbj,vubj ;
  pkvec_struct *aj ;
  ibrank_struct *archvars ;
  bool scaled,retval ;

  const char *rtnnme = "ib_archvrank" ;

# ifndef DYLP_NDEBUG

  int freecnt,onebndcnt,twobndcnt ;

  freecnt = 0 ;
  onebndcnt = 0 ;
  twobndcnt = 0 ;

# endif

  retval = TRUE ;

  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;

/*
  Did the client supply a vector? If not, allocate one.
*/
  if (*p_archvars == NULL)
  { archvars =
	(ibrank_struct *) MALLOC(dy_sys->archvcnt*sizeof(ibrank_struct)) ; }
  else
  { archvars = *p_archvars ; }
/*
  Scan the architecturals, recording the necessary information in archvars.
  Fixed variables should not be loaded, so we should not see them here. If
  the constraint system has not been scaled, recover the maximum coefficient
  for the column.
*/
  n = dy_sys->varcnt ;
  m = dy_sys->logvcnt ;
  scaled = dy_isscaled() ;
  eligible = 0 ;
  aj = pkvec_new(0) ;
  for (j = m+1 ; j <= n ; j++)
  { if (consys_getcol_pk(dy_sys,j,&aj) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,"column",
	     consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
      retval = FALSE ;
      break ; }
    archvars[eligible].ndx = j ;
    archvars[eligible].nonzero = aj->cnt ;
    if (scaled)
    { archvars[eligible].ajmax = 1.0 ; }
    else
    { archvars[eligible].ajmax = consys_infnormcol(dy_sys,j) ; }
    vlbj = vlb[j] ;
    vubj = vub[j] ;
    if (vlbj > -dy_tols->inf && vubj < dy_tols->inf)
    { archvars[eligible].bndcnt = 2 ;
#     ifdef PARANOIA
      if (vlbj == vubj)
      { errmsg(1,rtnnme,__LINE__) ;
	retval = FALSE ;
	break ; }
#     endif
#     ifndef DYLP_NDEBUG
      twobndcnt++ ;
#     endif
    }
    else
    if (vlbj > -dy_tols->inf || vubj < dy_tols->inf)
    { archvars[eligible].bndcnt = 1 ;
#     ifndef DYLP_NDEBUG
      onebndcnt++ ;
#     endif
    }
    else
    { archvars[eligible].bndcnt = 0 ;
#     ifndef DYLP_NDEBUG
      freecnt++ ;
#     endif
    }
    eligible++ ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	    "\n\t%d eligible architecturals: %d free, %d 1-bound, %d 2-bound",
	    eligible,freecnt,onebndcnt,twobndcnt) ; }
    # endif

  if (aj != NULL) pkvec_free(aj) ;
  if (retval == FALSE)
  { if (archvars != NULL) FREE(archvars) ;
    return (FALSE) ; }
/*
  And sort them. If there aren't any candidates, and we supplied the vector
  for archvars, release it.
*/
  if (eligible > 1)
  { qsort(&archvars[0],eligible,sizeof(ibrank_struct),ib_archvcomp) ; }
  else
  if (eligible == 0)
  { if (*p_archvars == NULL)
    { if (archvars != NULL) FREE(archvars) ;
      archvars = NULL ; } }
/*
  Return the sorted list.
*/
  *p_cnt = eligible ;
  *p_archvars = archvars ;

  return (TRUE) ; }




static int ib_archvselect (int cnt, ibrank_struct *vars)
/*
  This routine tries to populate the basis using the ranked architectural
  variables listed in vars. Equalities are examined first, because if we can't
  cover them with architectural variables, we'll have to use artificials.

  The variables in vars are already ranked according to how desireable it is
  to place them in the basis, taking into account bounds and number of
  nonzeros in a column. Here we'll try to allocate the variables to basis
  positions based on the pivot element a<ii>, and a notion that in an ideal
  world we'd like to construct a lower diagonal basis matrix.
  
  Why lower diagonal? Because we can turn it into a diagonal matrix, which is
  trivially easy to invert! To require that the columns we select conform
  strictly to the lower diagonal view, when we pick a pivot a<ii>, we would
  need to require that a<k,i> = 0, k < i.

  In practice, this is pretty restrictive. So we settle for requiring that
  elimination using multiples of rows a<k>, k < i, have minimal impact on our
  chosen pivot a<ii>.  We do this by requiring that coefficients a<k,i>, k <
  i, be small compared to the pivots a<kk>: specifically, a<k,i> < .1
  a<kk>. (If you're having trouble visualising this, you might want to go off
  and do a quick 3x3 example about now.)

  As we're selecting, we also try to pay a bit of attention to numerical
  stability.  If the matrix is scaled, the maximum coefficient in each row
  and column is 1.0, and we simply require |a<ij>| > .9. If the matrix is
  unscaled, well, we settle for |a<ij>/a<max,j>| > .9.  (Arguably we should
  also look at a<p,max>, but if you're really having stability trouble, you
  should consider allowing dylp to scale the matrix.)

  Note that we don't have to make marginal choices here --- we'll fill any
  vacant positions with the logical (slack or artificial) for the constraint.

  Parameters:
    cnt:	the number of variables in vars
    vars:	ranked vector of architectural variables

  Returns: the number of basis positions covered, or -1 if an error occurs.
*/

{ int i,j,m,ndx,pkndx,covered ;
  double ratio ;
  double *estpiv ;
  ibrank_struct *var ;

  bool scaled,select ;
  pkvec_struct *aj ;
  pkcoeff_struct *aij ;

  dyret_enum retval ;

  const char *rtnnme = "ib_archvselect" ;

  retval = 0 ;
  ratio = 0 ;
/*
  Allocate bookkeeping arrays.
*/
  m = dy_sys->concnt ;
  scaled = dy_isscaled() ;

  aj = NULL ;
  estpiv = (double *) MALLOC((m+1)*sizeof(double)) ;
  for (i = 1 ; i <= m ; i++) estpiv[i] = dy_tols->inf ;
  covered = 0 ;
/*
  Open a loop to scan vars and attempt to use each variable.
*/
  for (ndx = 0 ; ndx < cnt ; ndx++)
  { var = &vars[ndx] ;
    j = var->ndx ;
    if (consys_getcol_pk(dy_sys,j,&aj) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,"column",
	     consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
      retval = -1 ;
      break ; }
/*
  Scan the column: We need a coefficient in an uncovered row (dy_basis[i] ==
  0) which satisfies the stability criteria (ratio > .9). The coefficients we
  scan before finding a suitable pivot can't do too much damage to our goal
  of lower diagonal (|aij->val| < .1*estpiv[i]).

  NOTE: Scan order is crucially important here. We need to consider
	equalities before inequalities in order to minimize eventual use of
	artificials to cover equalities. Because of the way dy_sys is loaded
	(equalities first) and the way consys augments columns (new
	coefficients go at the front) we need to scan from back to front in
	order to consider equalities first.
*/
    select = FALSE ;
    for (pkndx = aj->cnt-1 ; pkndx >= 0 ; pkndx--)
    { aij = &aj->coeffs[pkndx] ;
      i = aij->ndx ;
      ratio = fabs(aij->val) ;
      if (!scaled) ratio /= var->ajmax ;
      if (dy_basis[i] == 0 && ratio > .9)
      { select = TRUE ;
	break ; }
      else
      { if (ratio > .1*estpiv[i]) break ; } }
/*
  Did we select this variable? If so, install it. If we have a full basis,
  break out of the search loop.
*/
    if (select == TRUE)
    { dy_basis[i] = j ;
      dy_var2basis[j] = i ;
      estpiv[i] = ratio ;
      covered++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.crash >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  adding %s (%d)",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," to cover %s (%d),",
		    consys_nme(dy_sys,'c',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," |a<%d,%d>/max(a<*,%d>)| = %g.",
		    i,j,j,ratio) ; }
#     endif
      if (covered == m) break ; }
#   ifndef DYLP_NDEBUG
    else
    if (dy_opts->print.crash >= 4)
    { if (pkndx < 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t  rejected %s (%d); lower diag violation at .1;",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    " a<%d,%d> = %g, estpiv<%d> = %g, ratio %g.",
		    i,j,ratio,i,estpiv[i],ratio/estpiv[i]) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t  rejected %s (%d) at .9; no suitable pivots.",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; } }
#   endif
  }
  if (aj != NULL) pkvec_free(aj) ;
  if (estpiv != NULL) FREE(estpiv) ;
  if (retval < 0) return (retval) ;
/*
  We're done. Print a summary and return.
*/

# ifndef DYLP_NDEBUG
  { if (dy_opts->print.crash >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    added %d architectural variables.",covered) ; } }
# endif

  return (covered) ; }



static int ib_slackselect (void)

/*
  This routine will populate the basis slots corresponding to inequalities
  using the logical (slack) variable for the inequality.

  Parameters: none

  Returns: the number of basis slots filled.
*/

{ int i,m,covered ;
  contyp_enum *ctyp ;

/*
  Walk the constraints, installing the logical for each inequality.
*/
  m = dy_sys->concnt ;
  ctyp = dy_sys->ctyp ;
  covered = 0 ;

  for (i = 1 ; i <= m ; i++)
  { if (ctyp[i] != contypEQ && dy_basis[i] == 0)
    { dy_basis[i] = i ;
      dy_var2basis[i] = i ;
      covered++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.crash >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  adding %s (%d)",
		    consys_nme(dy_sys,'v',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," to cover %s (%d).",
		    consys_nme(dy_sys,'c',i,FALSE,NULL),i) ; }
#     endif
    } }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    added %d slack/surplus variables.",covered) ; }
# endif

  return (covered) ; }


static int ib_artifselect (void)

/*
  This routine will populate the basis slots corresponding to equalities using
  the logical (artificial) variable for the equality.

  Parameters: none

  Returns: the number of basis slots filled
*/

{ int i,m,covered ;
  contyp_enum *ctyp ;

/*
  Walk the constraints, installing the logical for each inequality.
*/
  m = dy_sys->concnt ;
  ctyp = dy_sys->ctyp ;
  covered = 0 ;

  for (i = 1 ; i <= dy_sys->concnt ; i++)
  { if (ctyp[i] == contypEQ && dy_basis[i] == 0)
    { dy_basis[i] = i ;
      dy_var2basis[i] = i ;
      covered++ ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.crash >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  adding %s (%d)",
		    consys_nme(dy_sys,'v',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," to cover %s (%d).",
		    consys_nme(dy_sys,'c',i,FALSE,NULL),i) ; }
#     endif
    } }
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    added %d artificial variables.",covered) ; }
# endif

  return (covered) ; }



static bool ib_populatebasis (void)

/*
  This routine is responsible for populating the basis (i.e., selecting basic
  variables for each basis position) using a mixture of architectural and
  logical variables. The type of basis is determined by the coldbasis option,
  as follows:
    ibLOGICAL:	(logical) uses only logical variables (slacks and artificials).
    ibSLACK:	(slack) uses slacks for inequalities, then covers as many
		equalities as possible with architecturals before falling back
		on artificials.
    ibARCH:	(architectural) uses architecturals first, trying to cover
		equalities before inequalities, and falling back on logicals
		for uncovered positions.

  Parameters: none

  Returns: TRUE if a basis is built, FALSE if something goes wrong.
*/

{ int m,basiscnt,archvcnt,slkcnt,artifcnt,iretval ;
  bool bretval ;

  ibrank_struct *archvars ;

  const char *rtnnme = "ib_populatebasis" ;

# ifndef DYLP_NDEBUG
  int i,j ;
# endif

# ifdef PARANOIA
  if (!(dy_opts->coldbasis >= ibLOGICAL && dy_opts->coldbasis <= ibARCH))
  { errmsg(5,rtnnme,"initial basis type",dy_opts->coldbasis) ;
    return (FALSE) ; }
# endif

  m = dy_sys->concnt ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  constructing ") ;
    switch (dy_opts->coldbasis)
    { case ibLOGICAL:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"logical") ;
	break ; }
      case ibSLACK:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"slack") ;
	break ; }
      case ibARCH:
      { dyio_outfmt(dy_logchn,dy_gtxecho,"architectural") ;
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (FALSE) ; } }
    dyio_outfmt(dy_logchn,dy_gtxecho," basis for system %s, %d constraints.",
	        dy_sys->nme,m) ; }
# endif


  basiscnt = 0 ;

/*
  Logical and slack basis types will use slacks to cover off inequalities.
*/
  if (dy_opts->coldbasis == ibLOGICAL || dy_opts->coldbasis == ibSLACK)
  { slkcnt = ib_slackselect() ;
    basiscnt += slkcnt ; }
/*
  Slack and architectural basis types will use architectural variables. For a
  slack basis, we're just trying to cover equalities before falling back on
  artificial variables. For an architectural basis, we're covering equalities
  and inequalities, scanning the columns in a way that will see equalities
  before inequalities (hence preferentially use architecturals to cover
  equalities).

  Start by sorting them into their respective classes by number of bounds,
  then ordering within classes by the rating described in the header.

  Exclude fixed variables from the two-bound group. A simple equality test is
  all we need, as dy_coldstart has already groomed the bounds.
*/
  if (basiscnt < m &&
      (dy_opts->coldbasis == ibSLACK || dy_opts->coldbasis == ibARCH))
  { archvars = NULL ;
    bretval = ib_archvrank(&archvcnt,&archvars) ;
    if (bretval == FALSE)
    { errmsg(305,rtnnme,dy_sys->nme,
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "rank","architectural") ;
      if (archvars != NULL) FREE(archvars) ;
      return (FALSE) ; }
    iretval = ib_archvselect(archvcnt,archvars) ;
    if (archvars != NULL) FREE(archvars) ;
    if (iretval < 0)
    { errmsg(305,rtnnme,dy_sys->nme,
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "select","architectural") ;
      return (FALSE) ; }
    basiscnt += iretval ; }
/*
  Use logicals for any uncovered inequalities in an architectural basis.
*/
  if (basiscnt < m && dy_opts->coldbasis == ibARCH)
  { slkcnt = ib_slackselect() ;
    basiscnt += slkcnt ; }
/*
  If we have uncovered equalities (and that's all that can remain, at this
  point), give in and use artificial variables.
*/
  if (basiscnt < m)
  { artifcnt = ib_artifselect() ;
    basiscnt += artifcnt ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\t    Pos'n Variable           Constraint") ;
    for (i = 1 ; i <= basiscnt ; i++)
    { j = dy_basis[i] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t     %3d  (%3d) %-15s",i,j,
		  consys_nme(dy_sys,'v',j,FALSE,NULL)) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%-15s",
		  consys_nme(dy_sys,'c',i,FALSE,NULL)) ; } }
# endif

/*
  And if we still have uncovered rows, there's something seriously wrong.
*/
  if (basiscnt < m)
  { errmsg(301,rtnnme,basiscnt,m,dy_sys->nme) ;
    return (FALSE) ; }

  return (TRUE) ; }




dyret_enum dy_crash (void)

/*
  This routine coordinates the selection of an initial basis for an lp
  problem. The task of selecting a set of basic variables is hidden away in
  ib_populatebasis and subsidiary routines (see the comments at the head of
  this group of routines). What's visible here is the work associated with
  factoring the basis, establishing status, and calculating primal and dual
  variables and reduced costs.



  Parameters:  dy_lp:      lp problem dy_opts:    lp algorithm options
    dy_tols:    lp algorithm control and tolerances

  Returns: dyrOK if a basis is constructed; originates dyrFATAL on error, and
	   can relay dyrNUMERIC, dyrBSPACE, and dyrSINGULAR originated by
	   dy_factor.
*/

{ int vndx ;
  double *vub,*vlb,*obj ;
  flags calcflgs ;
  dyret_enum retval ;

  const char *rtnnme = "dy_crash" ;

  extern void dy_setfinalstatus(void) ;		/* dy_hotstart.c */

# ifndef DYLP_NDEBUG
  int cndx ;
# endif

# ifdef PARANOIA
  if (dy_lp == NULL)
  { errmsg(2,rtnnme,"dy_lp") ;
    return (dyrFATAL) ; }
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (dyrFATAL) ; }
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,consys_assocnme(dy_sys,CONSYS_MTX)) ;
    return (dyrFATAL) ; }
  if (dy_sys->vlb == NULL)
  { errmsg(2,rtnnme,consys_assocnme(dy_sys,CONSYS_VLB)) ;
    return (dyrFATAL) ; }
  if (dy_sys->obj == NULL)
  { errmsg(2,rtnnme,consys_assocnme(dy_sys,CONSYS_OBJ)) ;
    return (dyrFATAL) ; }
  if (dy_sys->ctyp == NULL)
  { errmsg(2,rtnnme,consys_assocnme(dy_sys,CONSYS_CTYP)) ;
    return (dyrFATAL) ; }
  if (flgoff(dy_sys->opts,CONSYS_LVARS))
  { errmsg(311,rtnnme,dy_sys->nme) ;
    return (dyrFATAL) ; }
  if (dy_opts == NULL)
  { errmsg(2,rtnnme,"dy_opts") ;
    return (dyrFATAL) ; }
  if (dy_tols == NULL)
  { errmsg(2,rtnnme,"dy_tols") ;
    return (dyrFATAL) ; }
# endif
/*
  Unpack a few arrays we'll use frequently, and attach the basis and inverse
  basis vectors to the constraint system. consys_attach will initialise them
  to 0.
*/
  obj = dy_sys->obj ;
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(int),(void **) &dy_basis) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"basis vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(int),(void **) &dy_var2basis) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"inverse basis vector") ;
    return (dyrFATAL) ; }
/*
  Populate the basis.
*/
  if (ib_populatebasis() == FALSE)
  { errmsg(302,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,"populate") ;
    return (dyrFATAL) ; }
/*
  Factor the basis. We don't want any of the primal or dual variables
  calculated just yet. If this fails we're in deep trouble.
*/
  if (dy_sys->concnt > 0)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.crash >= 2)
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tfactoring ...") ;
#   endif
    calcflgs = 0 ;
    retval = dy_factor(&calcflgs) ;
    switch (retval)
    { case dyrOK:
      case dyrPATCHED:
      { break ; }
      default:
      { errmsg(309,rtnnme,dy_sys->nme) ;
	return (retval) ; } } }
/*
  Attach and clear the vectors which will hold the status and values of
  variables, and the reduced cost.
*/
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(flags),(void **) &dy_status) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"status vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(double),(void **) &dy_xbasic) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"basic variable vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(double),(void **) &dy_x) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"primal variable vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_COL,
		    sizeof(double),(void **) &dy_y) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"dual variable vector") ;
    return (dyrFATAL) ; }
  if (consys_attach(dy_sys,CONSYS_ROW,
		    sizeof(double),(void **) &dy_cbar) == FALSE)
  { errmsg(100,rtnnme,dy_sys->nme,"reduced cost vector") ;
    return (dyrFATAL) ; }
/*
  Might as well work toward a dual feasible start. Calculate the duals, then
  the reduced costs, so we can make intelligent decisions about the status of
  the nonbasic variables.
*/
  dy_calcduals() ;
  if (dy_calccbar() == FALSE)
  { errmsg(384,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
    return (dyrFATAL) ; }
/*
  For nonbasic variables, we have to decide the proper status, based on
  number of bounds and the sign of the reduced cost (we're minimising,
  remember). Once we know the status, we can set a value in dy_x.  Nonbasic
  free variables are arbitrarily awarded a value of 0.  For basic variables,
  if we can't say they're fixed or free, assume strictly between bounds until
  we calculate their initial values. We won't bother to set dy_x and
  dy_xbasic until we've done the calculation.
*/

# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\testablishing initial status and reference frame ...") ;
# endif
  for (vndx = 1 ; vndx <= dy_sys->varcnt ; vndx++)
  { if (dy_var2basis[vndx] != 0)
    { if (vlb[vndx] == vub[vndx])
      { dy_status[vndx] = vstatBFX ; }
      else
      if (vlb[vndx] <= -dy_tols->inf && vub[vndx] >= dy_tols->inf)
      { dy_status[vndx] = vstatBFR ; }
      else
      { dy_status[vndx] = vstatB ; } }
    else
    { if (vlb[vndx] > -dy_tols->inf && vub[vndx] < dy_tols->inf)
      { if (vub[vndx] == vlb[vndx])
	{ dy_status[vndx] = vstatNBFX ;
	  dy_x[vndx] = vub[vndx] ; }
	else
	if (dy_cbar[vndx] >= 0)
	{ dy_status[vndx] = vstatNBLB ;
	  dy_x[vndx] = vlb[vndx] ; }
	else
	{ dy_status[vndx] = vstatNBUB ;
	  dy_x[vndx] = vub[vndx] ; } }
      else
      if (vlb[vndx] > -dy_tols->inf)
      { dy_status[vndx] = vstatNBLB ;
	dy_x[vndx] = vlb[vndx] ; }
      else
      if (vub[vndx] < dy_tols->inf)
      { dy_status[vndx] = vstatNBUB ;
	dy_x[vndx] = vub[vndx] ; }
      else
      { dy_status[vndx] = vstatNBFR ;
	dy_x[vndx] = 0 ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.crash >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t  %s (%d) %s",
		  consys_nme(dy_sys,'v',vndx,FALSE,NULL),vndx,
		  dy_prtvstat(dy_status[vndx])) ;
      if (flgon(dy_status[vndx],vstatNONBASIC|vstatNBFR))
	dyio_outfmt(dy_logchn,dy_gtxecho," with value %g.",dy_x[vndx]) ;
      else
	dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
  }
/*
  Ok, status is set. Now it's time to calculate initial values for the basic
  variables and objective.  Once we have values for the basic variables, see
  how bad it looks, in terms of primal infeasibility, and adjust the status
  for any that are pinned against a bound or out of bounds.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.crash >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\tcalculating basic variable values ...") ; }
# endif
  if (dy_calcprimals() == FALSE)
  { errmsg(316,rtnnme,dy_sys->nme) ;
    return (dyrFATAL) ; }
  dy_lp->z = dy_calcobj() ;
  dy_setfinalstatus() ;
/*
  Make the check for primal and/or dual feasibility, and set the initial
  simplex phase accordingly.
*/
  calcflgs = ladPRIMFEAS|ladPFQUIET|ladDUALFEAS|ladDFQUIET ;
  retval = dy_accchk(&calcflgs) ;
  if (retval != dyrOK)
  { errmsg(304,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
    return (retval) ; }
  if (flgoff(calcflgs,ladPRIMFEAS))
  { dy_lp->simplex.next = dyPRIMAL2 ; }
  else
  if (flgoff(calcflgs,ladDUALFEAS))
  { dy_lp->simplex.next = dyDUAL ; }
  else
  { dy_lp->simplex.next = dyPRIMAL1 ; }

# ifndef DYLP_NDEBUG
/*
  Some debug printing, to dump the initial basis & variables.
*/
  if (dy_opts->print.crash >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  Pos'n Constraint\tDual\tPrimal\t\t   Status\tValue") ;
    for (cndx = 1 ; cndx <= dy_sys->concnt ; cndx++)
    { vndx = dy_basis[cndx] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n   %4d  %-13s%7g",cndx,
		  consys_nme(dy_sys,'c',cndx,FALSE,NULL),dy_y[cndx]) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"  (%4d) %-13s %-7s %7g",vndx,
		  consys_nme(dy_sys,'v',vndx,FALSE,NULL),
		  dy_prtvstat(dy_status[vndx]),dy_x[vndx]) ; } }

  if (dy_opts->print.crash >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tinitial objective %g",dy_lp->z) ;
    if (dy_lp->infeascnt != 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", %d infeasible vars, infeas. = %g",
		  dy_lp->infeascnt,dy_lp->infeas) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,", target simplex %s.",
	        dy_prtlpphase(dy_lp->simplex.next,FALSE)) ; }
# endif

  return (dyrOK) ; }

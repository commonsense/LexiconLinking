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
  This file contains the routines specific to the primal simplex algorithm.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_primal.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_primal.c 269 2009-04-02 05:38:19Z lou $" ;



/*
  A word or three on the handling of infeasible variables during phase I. The
  underlying philosophy is that once a variable becomes feasible, it stays
  feasible, hence the number of infeasible variables decreases monotonically.

  To construct the phase I objective, dy_initp1obj makes a list of the
  infeasible variables, puts +/- 1.0 in the objective in the corresponding
  spots, and calculates the duals and reduced costs. Dylp minimises, so we
  use +1 for BUUB variables, -1 for BLLB variables. To install the phase I
  objective, the phase II objective is detached from dy_sys->obj and cached,
  and the phase I objective is attached in its place.

  With each iteration, tweakp1obj scans the list of infeasible variables to
  see what's happened and make adjustments. The possibilities are:

    * No variable became feasible, in which case nothing need be done.

    * Exactly one variable (the leaving variable) became feasible and became
      nonbasic. Its (now nonbasic) objective coefficient is set to 0. This
      change will not affect the duals (y = c<B>inv(B)) or any reduced costs
      except cbar<i> (cbar<N> = c<N> - yA<N>, so we only need to compensate
      for the change to c<i>). I suppose we could check to see if x<i> is now
      an attractive candidate to pivot into the basis, but it seems more
      trouble than it's worth.

    * One or more variables became feasible, and remained basic. In this
      case, all hell breaks loose, as we'll need to change objective
      coefficients for variables that are still basic, and the change will
      ripple everywhere. I think that a reasonably efficient update
      calculation might be possible, but for now the approach is to simply
      recalculate the duals and reduced costs. At least the column norms
      don't change.

  All variables which gain feasibility are removed from the list of
  infeasible variables, and the list is compressed. Eventually, it dwindles
  to nothing and we're feasible (hence optimal) and done with phase I.

  Now, if you believe this heart-warming story, I have a fine bridge you
  might want to purchase. What'll really happen (not always, but often enough
  to hurt) is that accumulating numeric inaccuracy will creep in. At some
  point we'll refactor and find that previously feasible variables have lost
  feasibility. The variables' status will be reset correctly, so that they'll
  be properly evaluated and handled during pivoting, EXCEPT, the objective
  coefficient will be incorrect (i.e., 0).  Rather than dance around trying
  to add to the list of infeasible variables, we'll simply let it dwindle to
  nothing and catch any variables that became infeasible in the meantime with
  the preoptimality check for primal feasibility (which, as a side effect,
  will make sure that the value of infeascnt in dy_lp is correct). If we fail
  the check (by way of some variables having lost feasibility), we run
  dy_initp1obj again and give it another shot.

  In order to implement an initial variable purge, prior to the first run of
  a simplex routine, we need to be able to install the phase I objective
  before entering the main dynamic simplex loop. Hence the dy_lp structure
  contains phase I objective information. We also have to make dy_initp1obj
  available to dylp.
  
  In the normal course of things, commonstart will initialise the phase I
  objective, and dy_finishup will free it. Nonetheless, dy_primal does a
  check prior to calling primal1, in case we've kicked back into phase I due
  to loss of feasibility. primal1 may also call dy_initp1obj if it finds
  formerly feasible variables have lost feasibility.
*/

/*
  Antidegeneracy comes in two strengths:
    * `Anti-degen lite' attempts to choose the leaving variable using a
      heuristic based on alignment of hyperplanes. There are 6 variations
      available. Two of them actually ignore hyperplane alignment. The other
      four try to choose the leaving variable by considering the relative
      alignment of the plane that will become tight to either the objective
      function or the edge direction for the pivot. See further comments in
      the relevant subroutines (primalout & subroutines) in dy_primalpivot.c.

    * If the lite approach fails, the heavy artillery is a strategy based on
      perturbation of the problem.
  
  When faced with numeric ill-conditioning, the idea is that we gradually boost
  the minimum pivot ratio. This will happen when groombasis has to make a major
  status correction, and when we encounter unexpected loss of dual feasibility
  (optimality) or primal feasibility once we think we've obtained either.
*/



static bool forcesuperbasic (void)

/*
  To put it bluntly, this routine is a bandaid. For good or ill, primal I
  knows that it should never see superbasic variables; they're not necessary
  when feasibility isn't at issue. But, we have the following failure
  scenario:  We're in phase II, and a pivot attempt results in a singular
  basis. dylp soldiers on, calling dy_accchk to refactor (which will patch
  the basis) and recalculate the primal and dual variables. But, sad to say,
  the primal feasibility check fails (we've removed some accumulated
  numerical inaccuracy by refactoring, or changed the numeric conditioning of
  the basis) and we revert to primal phase I. Which promptly fails because of
  the superbasic variable left by the patch.

  So, this bandaid scans the nonbasic variables and, if it finds a superbasic
  variable, adjusts it to the best bound, based on the objective coefficient.
  In the event that we adjust a variable, we'll call dy_accchk to make sure
  we have an accurate count of infeasible basic variables going into phase I.

  Parameters: none

  Returns: TRUE if all goes well, FALSE otherwise.
*/


{ int k,supercnt ;
  flags statk,checks ;
  double valk,lbk,ubk ;

/*
  Open up a loop to walk the variables, looking for superbasics.
*/
  supercnt = 0 ;
  for (k = 1 ; k < dy_sys->varcnt ; k++)
  { 
/*
  It's a superbasic variable. Force it to the appropriate bound based on the
  objective coefficient and presence/absence of the bound.
*/
    if (flgon(dy_status[k],vstatSB))
    { supercnt++ ;
      ubk = dy_sys->vub[k] ;
      lbk = dy_sys->vlb[k] ;
      statk = vstatSB ;
      if (ubk < dy_tols->inf && lbk > -dy_tols->inf)
      { if (dy_sys->obj[k] < 0)
	{ setflg(statk,vstatNBUB) ;
	  valk = ubk ; }
	else
	{ setflg(statk,vstatNBLB) ;
	  valk = lbk ; } }
      else
      if (ubk < dy_tols->inf)
      { setflg(statk,vstatNBUB) ;
	valk = ubk ; }
      else
      if (lbk > -dy_tols->inf)
      { setflg(statk,vstatNBLB) ;
	valk = lbk ; }
      else
      { setflg(statk,vstatNBFR) ;
	valk = 0 ; }
      comflg(dy_status[k],statk) ;
      dy_x[k] = valk ; }

#   ifdef DYLP_PARANOIA
/*
  If we're paranoid, run a status check while we're at it. This << must >>
  follow the code that forces superbasics. dy_chkstatus won't tolerate BLLB
  or BUUB in phase II, and we're here because of loss of feasbility, hence
  we need to show phase I. But dy_chkstatus won't tolerate SB in phase I, so
  we have to get rid of them before we check.
*/
    if (dy_chkstatus(k) == FALSE) return (FALSE) ;
#   endif
  }
/*
  Did we do anything? If so, use dy_accchk to do a primal feasibility check,
  which will ensure that the values in dy_lp are accurate.
*/
  if (supercnt > 0)
  { checks = ladPRIMFEAS|ladPFQUIET ;
    if (dy_accchk(&checks) != dyrOK) return (FALSE) ; }
  
  return (TRUE) ; }



static dyret_enum preoptimality (dyret_enum lpretval, flags *result)

/*
  This routine does the prep work so that we can have confidence in a report
  of optimality, infeasibility (phase I only), or unboundedness (phase II
  only) by the primal simplex routines. It clears the pivot reject list,
  backs out any restricted subproblems, refactors, recalculates the primal
  and (phase II only) dual variables, and performs an accuracy and
  feasibility check.

  Parameters:
    lpretval:	return code assigned by primal1 or primal2
    result:	(o) loaded with the result flags from dy_accchk
  
  Returns:
    dyrOK:	if all goes smoothly
    dyrPATCHED:	if the only bump is that the basis was patched by dy_factor
    dyrLOSTPFEAS: if the primal feasibility check by dy_accchk fails
    dyrLOSTDFEAS: if the dual feasibility check by dy_accchk fails (phase II
		only)
    dyrFATAL:	if anything else goes wrong

    Also can relay error codes from dy_accchk (dyrACCCHK, and various basis
    factoring errors). Loss of primal feasibility dominates loss of dual
    feasibility.
*/

{ flags checkflags ;
  dyret_enum retval ;

# if defined(DYLP_PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "preoptimality" ;
# endif

# ifndef DYLP_NDEBUG
  int print ;
# endif

# ifdef DYLP_PARANOIA
  if (!(lpretval == dyrOPTIMAL || lpretval == dyrINFEAS ||
	lpretval == dyrPUNT || lpretval == dyrUNBOUND))
  { errmsg(4,rtnnme,"lp return code",dy_prtdyret(lpretval)) ;
    return (dyrFATAL) ; }
# endif
# ifndef DYLP_NDEBUG
  if (dy_lp->phase == dyPRIMAL1)
    print = dy_opts->print.phase1 ;
  else
    print = dy_opts->print.phase2 ;
  
  if (print >= 4)
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: validating %s at iteration (%s)%d.",
	        rtnnme,dy_prtdyret(lpretval),
	        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
# endif
/*
  A little prep work. First, we don't want dy_accchk to take heroic measures
  in the event of loss of primal or dual feasibility, so suppress that.

  In phase I, because we're continually playing games with the objective
  function, and we've defined optimality to be primal feasibility, it's
  possible that we won't have dual feasibility when we come here reporting
  optimality. Also, if we go unbounded in phase I, we may not be primal
  or dual feasible.
  
  If this looks to be straightforward optimality and we've just refactored,
  don't request an initial refactor.
*/
  *result = 0 ;
  checkflags = 0 ;
  setflg(checkflags,ladFACTOR|ladPRIMALCHK|ladPRIMFEAS|ladPFQUIET|
		    ladDUALCHK|ladDUALFEAS|ladDFQUIET) ;
  if (lpretval == dyrOPTIMAL && dy_lp->basis.etas == 0)
    clrflg(checkflags,ladFACTOR) ;
/*
  Start with the easy stuff -- clear the pivot reject list and back out any
  restricted subproblems. If degenout notes accuracy loss, request a refactor.
*/
# ifndef DYLP_NDEBUG
  if (print >= 4)
    dyio_outfmt(dy_logchn,dy_gtxecho,
	    "\n\tclearing pivot rejection and antidegeneracy machinery ... ") ;
# endif
  if (dy_clrpivrej(NULL) != TRUE) return (dyrFATAL) ;
  if (dy_lp->degen > 0)
  { if (dy_degenout(0) == dyrREQCHK) setflg(checkflags,ladFACTOR) ; }
/*
  And now the accuracy checks. Failure here is hard failure --- dy_accchk is
  very persistent, as is dy_factor, and any problems would have been fixed if
  possible.
*/
# ifndef DYLP_NDEBUG
  if (print >= 4)
    dyio_outfmt(dy_logchn,dy_gtxecho,"done.\n\t%schecking accuracy ... ",
	        flgon(checkflags,ladFACTOR)?"refactoring and ":"") ;
# endif
  retval = dy_accchk(&checkflags) ;
  *result = checkflags ;
  if (!(retval == dyrOK || retval == dyrPATCHED))
  {
#   ifndef DYLP_NDEBUG
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%sfailed.",(print >= 5)?"\n\t":" ") ; }
#   endif
    return (retval) ; }
  else
  if (flgon(checkflags,ladPRIMALCHK|ladDUALCHK))
  {
#   ifndef DYLP_NDEBUG
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%sfailed",(print >= 5)?"\n\t":" ") ;
      if (flgon(checkflags,ladPRIMALCHK))
	dyio_outfmt(dy_logchn,dy_gtxecho," primal") ;
      if (flgon(checkflags,ladDUALCHK))
	dyio_outfmt(dy_logchn,dy_gtxecho," dual") ;
      dyio_outfmt(dy_logchn,dy_gtxecho," check(s).") ; }
#   endif
    retval = dyrACCCHK ; }
  else
  if (flgon(checkflags,ladPRIMFEAS|ladDUALFEAS))
  {
#   ifndef DYLP_NDEBUG
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%slost",(print >= 5)?"\n\t":" ") ;
      if (flgon(checkflags,ladPRIMALCHK))
	dyio_outfmt(dy_logchn,dy_gtxecho," primal") ;
      if (flgon(checkflags,ladDUALCHK))
	dyio_outfmt(dy_logchn,dy_gtxecho," dual") ;
      dyio_outfmt(dy_logchn,dy_gtxecho," feasibility.") ; }
#   endif
    if (flgon(checkflags,ladPRIMFEAS))
      retval = dyrLOSTPFEAS ;
    else
      retval = dyrLOSTDFEAS ; }
# ifndef DYLP_NDEBUG
  else
  { if (print >= 4)
      dyio_outfmt(dy_logchn,dy_gtxecho,"%s%s.",(print >= 5)?"\n\t":" ",
		  (retval == dyrOK)?"done":"patched") ; }
# endif

  return (retval) ; }
  


bool dy_swapobjs (dyphase_enum phase)

/*
  This routine handles the allocation, exchange, and deallocation of phase I
  and phase II objective vectors.

  The actions are as follows, depending on the value of the phase parameter:
    dyPRIMAL1: We're headed into primal phase I. If the vectors for infvars
	       and the phase I objective are not yet allocated, do it, and
	       attach p1obj and p2obj as additional pointers to the phase I
	       and II objectives, respectively. The actual swap consists of
	       detaching dy_sys->obj as a pointer to the phase II objective
	       and reattaching it as a pointer to the phase I objective. It
	       may be that the phase I objective is already installed; in
	       that case we have only to check infvars for size.
    
    dyPRIMAL2: We're headed into primal phase II. Detach dy_sys->obj as a
	       pointer to the phase I objective, and reattach it as a pointer
	       to the phase II objective.

    dyDONE:    Detach the lot and free infvars and p1obj.

  Parameters:
    phase:	dyPRIMAL1 to swap the phase I objective into place, dyPRIMAL2
		to swap the phase II objective into place, dyDONE to clean up.
  
  Returns: TRUE if the swap goes ok, FALSE otherwise
*/

{ const char *rtnnme = "dy_swapobjs" ;

# ifdef DYLP_PARANOIA
/*
  A little paranoia. The routine will do the right thing (i.e., nothing) if
  called to install the P2 objective and it's already in place. But chances
  are we're confused if it happens. There are valid reasons to want to 
  reinitialize the P1 objective, so no warning is issued.
*/
  if (!(phase == dyPRIMAL1 || phase == dyPRIMAL2 || phase == dyDONE))
  { errmsg(4,rtnnme,"direction",dy_prtlpphase(phase,FALSE)) ;
    return (FALSE) ; }
  if (dy_lp->p1obj.installed == FALSE && phase == dyPRIMAL2)
  { warn(399,rtnnme,
	 dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,"II") ; }
# endif
/*
  We're here to install the phase I objective. If it's already installed, all
  we need to do is check infvars for size.

  If it's not installed, the first thing to check is whether we need to
  (re)allocate vectors for infvars and the phase I objective and attach
  pointers p1obj and p2obj. Then do the routine part of the swap, detaching
  dy_sys->obj as a pointer to the phase II objective and reattaching it as a
  pointer to the phase I objective.
*/
  if (phase == dyPRIMAL1)
  { if (dy_lp->p1obj.installed == TRUE)
    { if (dy_lp->infeascnt > dy_lp->p1obj.infvars_sze)
      { dy_lp->p1obj.infvars_sze = dy_lp->infeascnt ;
	dy_lp->p1obj.infvars = (int *)
	    REALLOC(dy_lp->p1obj.infvars,dy_lp->infeascnt*sizeof(int)) ; } }
    else
    { if (dy_lp->p1obj.p1obj == NULL)
      { dy_lp->p1obj.infvars = (int *) MALLOC(dy_lp->infeascnt*sizeof(int)) ;
        dy_lp->p1obj.infvars_sze = dy_lp->infeascnt ;
	dy_lp->p1obj.p1obj = NULL ;
	if (consys_attach(dy_sys,CONSYS_OBJ,sizeof(double),
			  (void **) &dy_lp->p1obj.p1obj) == FALSE)
	{ errmsg(100,rtnnme,dy_sys->nme,"&dy_lp->p1obj.p1obj") ;
	  return (FALSE) ; }
	dy_lp->p1obj.p2obj = dy_sys->obj ;
	if (consys_attach(dy_sys,CONSYS_OBJ,sizeof(double),
			  (void **) &dy_lp->p1obj.p2obj) == FALSE)
	{ errmsg(100,rtnnme,dy_sys->nme,"&dy_lp->p1obj.p2obj") ;
	  return (FALSE) ; } }
      else
      { if (dy_lp->infeascnt > dy_lp->p1obj.infvars_sze)
	{ dy_lp->p1obj.infvars_sze = dy_lp->infeascnt ;
	  dy_lp->p1obj.infvars = (int *)
	      REALLOC(dy_lp->p1obj.infvars,dy_lp->infeascnt*sizeof(int)) ; } }

    if (consys_detach(dy_sys,(void **) &dy_sys->obj,FALSE) == FALSE)
    { errmsg(105,rtnnme,dy_sys->nme,"&dy_sys->obj (P2)") ;
      return (FALSE) ; }
    dy_sys->obj = dy_lp->p1obj.p1obj ;
    if (consys_attach(dy_sys,CONSYS_OBJ,
		      sizeof(double),(void **) &dy_sys->obj) == FALSE)
    { errmsg(100,rtnnme,dy_sys->nme,"&dy_sys->obj (P1)") ;
      return (FALSE) ; }
    dy_lp->p1obj.installed = TRUE ; } }
/*
  We're here to remove the phase I objective and reattach the phase II
  objective. Detach dy_sys->obj as a pointer to the phase I objective and
  reattach it as a pointer to the phase II objective.
*/
  else
  if (phase == dyPRIMAL2)
  { if (dy_lp->p1obj.installed == FALSE) return (TRUE) ;
    if (consys_detach(dy_sys,(void **) &dy_sys->obj,FALSE) == FALSE)
    { errmsg(105,rtnnme,dy_sys->nme,"&dy_sys->obj (P1)") ;
      return (FALSE) ; }
    dy_sys->obj = dy_lp->p1obj.p2obj ;
    if (consys_attach(dy_sys,CONSYS_OBJ,
		      sizeof(double),(void **) &dy_sys->obj) == FALSE)
    { errmsg(100,rtnnme,dy_sys->nme,"&dy_sys->obj (P2)") ;
      return (FALSE) ; } ;
    dy_lp->p1obj.installed = FALSE ; }
/*
  We're finishing up and releasing the dylp data structures. We need to free
  infvars, and whichever objective isn't currently installed (the other will
  be freed when the dy_sys constraint system is freed).
*/
  else
  { if (dy_lp->p1obj.infvars != NULL) FREE(dy_lp->p1obj.infvars) ;
    if (dy_lp->p1obj.installed == TRUE)
    { if (dy_lp->p1obj.p2obj != NULL) FREE(dy_lp->p1obj.p2obj) ; }
    else
    { if (dy_lp->p1obj.p1obj != NULL) FREE(dy_lp->p1obj.p1obj) ; } }
  
  return (TRUE) ; }



bool dy_initp1obj (void)

/*
  This routine is responsible for initialising the phase I objective.
  dy_swapobjs takes care of making sure the associated structures are ok.
  
  Initialisation consists of scanning the basis to fill in infvars with the
  indices of infeasible variables, and setting the appropriate values in the
  phase I objective. Once the objective function is established, we'll call
  dy_calcduals and dy_calccbar to calculate the duals and reduced costs.

  Parameters: none
  
  Returns: TRUE if the initialisation goes through without problem, FALSE
	   otherwise.
*/

{ int *infvars,infcnt,xipos,xindx ;
  double *p1obj ;
  const char *rtnnme = "dy_initp1obj" ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.phase1 >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    initialising phase 1 objective and reduced costs.") ; }
# endif

/*
  Call dy_swapobjs to have a look at the infvars and phase I objective
  vectors.  They'll be (re)allocated as necessary, and the phase I objective
  will be installed in place of the phase II objective.  Note that whenever
  we're in here, dy_lp->infeascnt should be correct, but we can't confirm
  that with a paranoid check until after we scan the basis.
*/
  if (dy_swapobjs(dyPRIMAL1) == FALSE)
  { errmsg(318,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	   dy_lp->tot.iters,"install/resize") ;
    return (FALSE) ; }
/*
  Clear the phase I objective to 0. We can track the changes under most
  circumstances, but arbitrary changes due to recovery from loss of accuracy
  or basis singularity would defeat us.
*/
  infvars = dy_lp->p1obj.infvars ;
  p1obj = dy_lp->p1obj.p1obj ;
  infcnt = 0 ;
  memset(p1obj,0,(dy_sys->varcnt+1)*sizeof(double)) ;
/*
  Scan the basis and record the indices of infeasible variables. For each
  variable, set the objective coefficient
  to -1 for variables above their upper bound and +1 for variables below
  their lower bound.
*/
  for (xipos = 1 ; xipos <= dy_sys->concnt ; xipos++)
  { xindx = dy_basis[xipos] ;
    if (flgoff(dy_status[xindx],vstatBLLB|vstatBUUB)) continue ;
    infvars[infcnt++] = xindx ;
    if (flgon(dy_status[xindx],vstatBLLB))
    { dy_sys->obj[xindx] = -1.0 ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase1 >= 7)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t%16s (%3d) = %16.8g < lb = %16.8g, infeas = %16.8g",
		    consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		    dy_xbasic[xipos],dy_sys->vlb[xindx],
		    dy_sys->vlb[xindx]-dy_xbasic[xipos]) ; }
#     endif
    }
    else
    { dy_sys->obj[xindx] = 1.0 ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase1 >= 7)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t%16s (%3d) = %16.8g > ub = %16.8g, infeas = %16.8g",
		    consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		    dy_xbasic[xipos],dy_sys->vub[xindx],
		    dy_xbasic[xipos]-dy_sys->vub[xindx]) ; }
#     endif
    } }
  dy_lp->p1obj.infcnt = infcnt ;

# ifdef DYLP_PARANOIA
/*
  Any time that we're in here, dy_lp->infeascnt should be accurate --- either
  we're just starting up, or we've kicked back here after discovering a lack
  of feasibility. In the latter case, a primal feasibility check should have
  run. If we don't match, there's a problem.
*/
  if (infcnt != dy_lp->infeascnt)
  { errmsg(1,rtnnme,__LINE__) ;
    return (FALSE) ; }
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.phase1 >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n      saw %d infeasible variables, tot. infeas. %g.",
	        infcnt,dy_lp->infeas) ; }
# endif

/*
  Calculate the duals and reduced costs for the objective, and the objective
  itself.
*/
  dy_calcduals() ;
  if (dy_calccbar() == FALSE)
  { errmsg(384,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	   dy_lp->tot.iters) ;
    return (FALSE) ; }
  dy_lp->z = dy_calcobj() ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.phase1 >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\t  recalculated duals and reduced costs.") ; }
# endif
# ifdef DYLP_STATISTICS
/*
  Keep track of the number of pivots between changes in feasibility status.
*/
  if (dy_stats != NULL)
  { int pivcnt ;
    dy_stats->infeas.chgcnt2++ ;
    if (dy_stats->infeas.maxcnt < infcnt) dy_stats->infeas.maxcnt = infcnt ;
    pivcnt = dy_lp->tot.pivs-dy_stats->infeas.prevpiv ;
    dy_stats->infeas.totpivs += pivcnt ;
    if (pivcnt > dy_stats->infeas.maxpivs) dy_stats->infeas.maxpivs = pivcnt ;
    dy_stats->infeas.prevpiv = dy_lp->tot.pivs ; }
# endif

  return (TRUE) ; }



static dyret_enum tweakp1obj (bool *reselect, int candxj)

/*
  This routine is responsible for updating the phase I objective function as
  variables gain feasibility.

  If no variable gains feasibility, there's nothing to do --- the objective
  is unchanged.

  If only a single variable has gained feasibility with this pivot, the best
  case involves two variables, such that x<i> left the basis and gained
  feasibility as x<j> entered.  In this case, a little thought brings the
  realisation that since x<i> is now nonbasic, we can change its objective
  coefficient without disturbing the duals (y = c<B>inv(B)).  Since the
  reduced costs are cbar<N> = c<N>-yA<N>, we only have to adjust cbar<i>
  by -c<i>, before setting c<i> to 0. This doesn't entirely get us
  out of trouble (it could be that x<i> has been selected to reenter, and
  when we rewrite cbar<i> we make it unsuitable).

  If a variable gains feasibility and remains basic, all hell breaks loose.
  It looks to me like a reasonably efficient update operation is possible,
  but it's going to be complex. For the nonce, just recalculate all the
  duals and reduced costs. At least the column norms don't change.

  Fixed variables that move from BLLB to BUUB fall in the `all hell breaks
  loose' category.

  Parameters:
    reselect:	(o) set to TRUE if the reduced costs have been recalculated
		    and a new entering variable should be selected, set to
		    FALSE otherwise
    candxj:	the candidate entering variable x<j'>

  Returns: dyrINFEAS if the problem is still infeasible
	   dyrOPTIMAL if the problem is feasible (thus optimal w.r.t. the
		     phase I objective)
	   dyrFATAL for errors.
*/

{ int *infvars,ndx,maxndx,newfeas,xkndx ;
  flags statk ;
  bool recalccbar ;
  const char *rtnnme = "tweakp1obj" ;

# ifndef DYLP_NDEBUG
  int xkpos ;
  double infeas ;
# endif
# ifdef DYLP_PARANOIA
  double chkz ;
# endif

# ifndef DYLP_NDEBUG
  infeas = 0 ;
  if (dy_opts->print.phase1 >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\t  checking feasibility & tweaking phase 1 objective.") ; }
# endif
# ifdef DYLP_PARANOIA
  chkz = dy_calcobj() ;
  if (!withintol(chkz,dy_lp->z,fabs(1000*dy_tols->cost*(1+fabs(chkz)))))
  { warn(405,rtnnme,
	 dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	 dy_lp->z,chkz,dy_lp->z-chkz,fabs(1000*dy_tols->cost*(1+chkz))) ; }
# endif

/*
  Walk the list of infeasible variables, looking for changes in status and
  adjusting the objective and reduced cost accordingly. Each newly feasible
  variable is compressed out of infvars.

  If the only change from infeasible to feasible is the leaving variable
  x<i>, we can make the adjustment here in the loop. Otherwise, we just note
  that a variable has gained feasibility or changed infeasibility and
  remained basic, and recalculate after we leave the loop.
*/
  newfeas = 0 ;
  maxndx = dy_lp->p1obj.infcnt-1 ;
  *reselect = FALSE ;
  recalccbar = FALSE ;
  infvars = dy_lp->p1obj.infvars ;
  for (ndx = 0 ; ndx <= maxndx ; )
  { xkndx = infvars[ndx] ;
    statk = getflg(dy_status[xkndx],vstatSTATUS) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.phase1 >= 5 && flgon(statk,vstatBLLB|vstatBUUB))
    { xkpos = dy_var2basis[xkndx] ;
      if (dy_opts->print.phase1 >= 7)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t    %4s %16s (%3d) = %16.8g",
		    dy_prtvstat(statk),consys_nme(dy_sys,'v',xkndx,FALSE,NULL),
		    xkndx,dy_xbasic[xkpos]) ; }
      if (flgon(statk,vstatBLLB))
      { infeas += dy_sys->vlb[xkndx]-dy_xbasic[xkpos] ;
	if (dy_opts->print.phase1 >= 7)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," < lb = %16.8g, inf = %16.8g",
		      dy_sys->vlb[xkndx],
		      dy_sys->vlb[xkndx]-dy_xbasic[xkpos]) ; } }
      else
      { infeas += dy_xbasic[xkpos]-dy_sys->vub[xkndx] ;
	if (dy_opts->print.phase1 >= 7)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," > ub = %16.8g, inf = %16.8g",
		      dy_sys->vub[xkndx],
		      dy_xbasic[xkpos]-dy_sys->vub[xkndx]) ; } } }
#   endif
    if (flgon(statk,vstatBLLB))
    { ndx++ ;
      if (dy_sys->obj[xkndx] != -1.0)
      { dy_sys->obj[xkndx] = -1.0 ;
	recalccbar = TRUE ; } }
    else
    if (flgon(statk,vstatBUUB))
    { ndx++ ;
      if (dy_sys->obj[xkndx] != 1.0)
      { dy_sys->obj[xkndx] = 1.0 ;
	recalccbar = TRUE ; } }
    else
    { newfeas++ ;
      infvars[ndx] = infvars[maxndx--] ;
      if (flgon(statk,vstatNONBASIC|vstatNBFR))
      { dy_cbar[xkndx] -= dy_sys->obj[xkndx] ;
	setcleanzero(dy_cbar[xkndx],dy_tols->zero) ;
	dy_lp->z -= dy_sys->obj[xkndx]*dy_x[xkndx] ;
	if (xkndx == candxj) *reselect = TRUE ; }
      else
      { recalccbar = TRUE ; }
      dy_sys->obj[xkndx] = 0.0 ; } }
  dy_lp->p1obj.infcnt = maxndx+1 ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.phase1 >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\t  saw %d infeasible variables, down %d, tot. infeas. %g.",
	        dy_lp->p1obj.infcnt,newfeas,infeas) ; }
  if (dy_opts->print.phase1 >= 2 && *reselect == TRUE)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    reselect; newly feasible %s (%d) selected to enter.",
	        consys_nme(dy_sys,'v',candxj,FALSE,NULL),candxj) ; }
# endif
# ifdef DYLP_PARANOIA
/*
  Are we feeling paranoid? Do a scan to make sure that objective coefficients
  are 0 for all variables that are in bound. (At this point, there may be a
  few newly infeasible variables with coefficients of 0, when they should be
  +/- 1.0, hence we don't check for those.)

  In the simple case where the leaving variable is the only variable to gain
  feasibility, we can test the iteratively updated objective.
*/
  for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
  { statk = dy_status[xkndx] ;
    if (flgoff(statk,vstatBLLB|vstatBUUB) && dy_sys->obj[xkndx] != 0.0)
    { errmsg(389,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
	     dy_prtvstat(statk),xkndx,dy_sys->obj[xkndx]) ;
      return (dyrFATAL) ; } }
  if (newfeas == 1 && recalccbar == FALSE)
  { chkz = dy_calcobj() ;
    if (!withintol(chkz,dy_lp->z,fabs(1000*dy_tols->cost*(1+fabs(chkz)))))
    { warn(405,rtnnme,dy_sys->nme,
	   dy_prtlpphase(dy_lp->phase,TRUE),-dy_lp->tot.iters,
	   dy_lp->z,chkz,dy_lp->z-chkz,
	   fabs(1000*dy_tols->cost*(1+chkz))) ; } }
# endif
# ifdef DYLP_STATISTICS
/*
  Keep track of the number of pivots between changes in feasibility status
  which require recalculation of the duals and reduced costs.
*/
  if (dy_stats != NULL && newfeas != 0)
  { int pivcnt ;
    if (recalccbar == FALSE)
      dy_stats->infeas.chgcnt1++ ;
    else
      dy_stats->infeas.chgcnt2++ ;
    pivcnt = dy_lp->tot.pivs-dy_stats->infeas.prevpiv ;
    dy_stats->infeas.totpivs += pivcnt ;
    if (pivcnt > dy_stats->infeas.maxpivs) dy_stats->infeas.maxpivs = pivcnt ;
    dy_stats->infeas.prevpiv = dy_lp->tot.pivs ; }
# endif
/*
  If dy_lp->p1obj.infcnt == 0, then we've gained feasibility (at least with
  respect to the variables in infvars; see notes at the head of the file).
  Return with an indication of optimality. (Note that we might not really be
  optimal w.r.t. the phase I objective, but that's not important.) We'll clear
  dy_y and cbar, just to keep the math straight.
*/
  if (dy_lp->p1obj.infcnt == 0)
  { memset(dy_y,0,(dy_sys->concnt+1)*sizeof(double)) ;
    memset(dy_cbar,0,(dy_sys->varcnt+1)*sizeof(double)) ;
    return (dyrOPTIMAL) ; }
/*
  If recalccbar is TRUE we need to redo the duals, reduced costs, and
  objective.  And we'll need to reselect the entering variable.   Sigh.
*/
  if (recalccbar == TRUE)
  { dy_calcduals() ;
    if (dy_calccbar() == FALSE)
    { errmsg(384,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters) ;
      return (dyrFATAL) ; }
    dy_lp->z = dy_calcobj() ;
    *reselect = TRUE ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.phase1 >= 5)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\trecalculated duals and reduced costs.",
		  maxndx,infeas) ; }
#   endif
  }
/*
  That's it. Return an indication that we're still infeasible.
*/
  return (dyrINFEAS) ; }



static dyret_enum verifyp1obj (void)

/*
  This routine attempts to verify that the phase I objective is indeed
  correct. The main thing we're looking for is unexpected loss or gain of
  feasibility due to accumulated numeric inaccuracy.  We need to do this when
  we have an indication of optimality in phase I, because of the problems
  outlined at the head of the file. If new variables become infeasible after
  infvars is loaded, tweakp1obj will not assign the proper objective
  coefficients. If the objective is all 0's, there's no motivation to reduce
  infeasibility, eh?

  The quick test is equality of dy_lp->p1obj.infcnt and dy_lp->infeascnt. If
  that fails, we have a problem. Even if we have equal counts, we still have
  to check the individual variables, to make sure we haven't had a pair of
  variables swap roles, or had a fixed variable move from BLLB to BUUB.

  If the objective doesn't verify, it's an indication that we're having
  accuracy problems. Boost the pivot selection tolerances and do a refactor
  on the way out.

  Parameters: none

  Returns: dyrOK if the objective verifies,
	   dyrRESELECT if the objective didn't verify, and
	   dyrFATAL, or error code from dy_accchk, if something goes wrong.
*/

{ int *infvars,ndx,xkndx ;
  double ck ;
  flags statk,checks ;
  dyret_enum retval ;
  bool err ;

  retval = dyrINV ;
/*
  Use dy_accchk to do a primal feasibility check, which will ensure that
  the values in dy_lp are accurate. Because we're only asking for a
  recalculation, dy_accchk can return only dyrOK or dyrFATAL.
*/
  checks = 0 ;
  setflg(checks,ladPRIMFEAS|ladPFQUIET) ;
  retval = dy_accchk(&checks) ;
  if (retval != dyrOK) return (retval) ;

/*
  First compare dy_lp->p1obj.infcnt to dy_lp->infeascnt. If they're unequal,
  it's a cinch the objective is incorrect.

  Otherwise, open up a loop to step through infvars and check each entry.  If
  we're debugging, we'll check them all, but if we're not we break on the
  first inconsistency.
*/
  if (dy_lp->p1obj.infcnt != dy_lp->infeascnt)
  { retval = dyrRESELECT ; }
  else
  { retval = dyrOK ;
    err = FALSE ;
    infvars = dy_lp->p1obj.infvars ;
    for (ndx = 0 ; ndx < dy_lp->p1obj.infcnt ; ndx++)
    { xkndx = infvars[ndx] ;
      statk = getflg(dy_status[xkndx],vstatSTATUS) ;
      ck = dy_sys->obj[xkndx] ;
      switch (statk)
      { case vstatBLLB:
	{ if (ck != -1.0) err = TRUE ;
	  break ; }
        case vstatBUUB:
	{ if (ck != 1.0) err = TRUE ;
	  break ; }
	default:
	{ err = TRUE ;
	  break ; } }
      if (err == TRUE)
      { retval = dyrRESELECT ;
#       ifdef DYLP_NDEBUG
	break ;
#       else
	err = FALSE ;
	if (dy_opts->print.phase1 >= 5)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tphase I c<%d> = %g inconsistent for %s %s (%d) = %g;",
		  xkndx,dy_sys->obj[xkndx],dy_prtvstat(statk),
		  consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,dy_x[xkndx]) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho," lb = %g, ub = %g",
		      dy_sys->vlb[xkndx],dy_sys->vub[xkndx]) ;
	  if (!withinbnds(dy_sys->vlb[xkndx],dy_x[xkndx],dy_sys->vub[xkndx]))
	  { if (flgon(statk,vstatBLLB))
	    { dyio_outfmt(dy_logchn,dy_gtxecho,", infeas = %g.",
			  dy_sys->vlb[xkndx]-dy_x[xkndx]) ; }
	    else
	    { dyio_outfmt(dy_logchn,dy_gtxecho,", infeas = %g.",
			  dy_x[xkndx]-dy_sys->vub[xkndx]) ; } }
	  else
	  { dyio_outchr(dy_logchn,dy_gtxecho,'.') ; } }
#       endif
       } } }
/*
  If the objective doesn't check, it indicates that we've lost numerical
  accuracy somewhere along the way.  Boost the pivot selection tolerances and
  refactor now, for two reasons: We'll get back accuracy, if we can, and
  (more subtle) we'll clean up any inconsistency between value and status
  before initp1obj complains about it.
*/
  if (retval == dyrRESELECT)
  { (void) dy_setpivparms(+1,+1) ;
    checks = 0 ;
    setflg(checks, ladFACTOR|ladPRIMFEAS|ladPFQUIET) ;
    retval = dy_accchk(&checks) ;
    if (retval == dyrOK || retval == dyrPATCHED) retval = dyrRESELECT ; }

  return (retval) ; }



static dyret_enum primal1 (void)

/*
  Phase 1 of the primal simplex. There are two nested loops, an inner pivoting
  loop and an outer control loop.
  
  The pivoting loop is a three step sequence: pivot (dy_primalpivot), deal
  with any problems & routine maintenance (dy_duenna), then check for primal
  feasibility and update the objective (tweakp1obj).  As a side effect of the
  pivot, dy_primalpivot returns the index of the preferred candidate, xjcand,
  to enter on the next pivot. (Except when the pivot is a bound-to-bound swing
  by x<j>; in this case, there's no pricing update, hence no selection of a
  candidate x<j>, and we have to do it here.)

  The biggest difference between phase I and phase II is that in phase I we
  tweak the objective each time a variable gains feasibility. There's a lot
  of action hidden in dy_initp1obj and tweakp1obj.

  Exits from the pivoting loop fall into two basic groups:
    * We need to select a new entering x<j>, but there's no other problem.
    * We have something that looks like a termination condition --- optimality,
      infeasibility, unboundedness, a punt, or some fatal problem.
  
  When we need to select a new x<j>, we circle back to the top of the control
  loop where dy_primalin is called to select a new incoming variable. As long
  as an incoming variable can be selected, and nothing fatal goes wrong, the
  inner pivoting loop resumes.
  
  When we suspect a termination condition, things get more complex. In phase
  I, optimality is defined to be primal feasibility and is detected by
  tweakp1obj. Inability to find a candidate for entry (a.k.a. dual
  feasbility, the normal condition for optimality, used by dy_primalin and
  dy_primalpivot) actually equates to infeasibility. A punt is just
  (infeasible) optimality with the possibility we could make further progress
  if we relaxed the pivot selection tolerance. Unboundedness remains
  unboundedness.

  Because life is not fair (see remarks at head of file), when we suspect
  feasibility the first thing we need to do is check that we're really
  considering all infeasible variables. If necessary, we reset the phase I
  objective and resume pivoting.  Otherwise, it's the normal sequence for
  termination: call preoptimality to refactor and do accuracy and feasibility
  checks, then react accordingly.

  Parameters: none

  Returns: most of the dyret_enum codes.
*/

{ dyret_enum lpretval,scanresult,pivresult,duennaresult,preopresult,
	     p1objresult ;
  int startcol,scan,nextcol,candxj,xjndx,indir,xindx,outdir,optcnt ;
  double cbarj,abarij,delta ;
  bool do_pivots,reselect ;
  flags xjstatus,checks ;
  const char *rtnnme = "primal1" ;

# ifndef DYLP_NDEBUG
  bool uxpfeas,uxnpfeas,uxdfeas,uxndfeas ;
  dyret_enum tmpretval ;
# endif

# ifdef DYLP_PARANOIA
  if (dy_lp->degen != 0)
  { errmsg(317,rtnnme,dy_sys->nme,dy_lp->degen) ;
    return (dyrFATAL) ; }
  if (dy_lp->p1.iters != 0)
  { errmsg(5,rtnnme,"phase 1 iteration count",dy_lp->p1.iters) ;
    return (dyrFATAL) ; }
  if (dy_lp->infeascnt <= 0)
  { errmsg(319,rtnnme,dy_sys->nme) ;
    return (dyrFATAL) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->infeas.prevpiv = dy_lp->tot.pivs ;
# endif

  dy_lp->p1.pivs = 0 ;
  dy_lp->pivok = FALSE ;
  dy_lp->prev_pivok = FALSE ;
  lpretval = dyrINV ;
  if (dy_clrpivrej(NULL) != TRUE) return (dyrFATAL) ;
/*
  Initialise the phase I objective, if it's not already in place.
*/
  if (dy_lp->p1obj.installed == FALSE)
  { if (dy_initp1obj() == FALSE)
    { errmsg(318,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"initialise") ;
      return (dyrFATAL) ; } }
/*
  Open the control loop. The purpose of this outer loop is to allow easy
  recovery from false indications of optimality, infeasibility, or
  unboundedness, as well as problems involving pivot selection. All the
  recovery actions happen at the bottom of the loop, after we fall out of the
  pivoting loop. Initialise a few variables and get going.
*/
  scan = dy_opts->scan ;
  nextcol = 1 ;
  optcnt = 0 ;
  while (lpretval == dyrINV)
  { 
/*
  Scan to select an incoming variable x<j>.  If primalin cannot find a
  candidate, we're infeasible (or perhaps we've punted). We're not primal
  feasible (else tweakp1obj would have seen it).

  If we've punted, call dy_dealWithPunt to free up any potential candidates on
  the pivot rejection list and iterate to try again.
*/
    startcol = nextcol ;
    scanresult = dy_primalin(startcol,scan,&candxj,&nextcol) ;
    switch (scanresult)
    { case dyrOK:
      { do_pivots = TRUE ;
	break ; }
      case dyrPUNT:
      { scanresult = dy_dealWithPunt() ;
	if (scanresult == dyrRESELECT)
	{ continue ; }
	else
	{ do_pivots = FALSE ;
	  lpretval = scanresult ; }
	break ; }
      case dyrOPTIMAL:
      { do_pivots = FALSE ;
	lpretval = dyrINFEAS ;
	break ; }
      default:
      { do_pivots = FALSE ;
	lpretval = scanresult ;
	break ; } }
/*
  Open the pivoting loop. While we have a candidate x<j> for entry, we do a
  three-step: attempt the pivot with x<j> (dy_primalpivot), then check that
  everything went off ok (dy_duenna), and finally check for feasibility and
  revise the objective. As part of updating the PSE pricing information, a
  new x<j> will be selected (but note that changes in feasibility status may
  render this choice obsolete, forcing a reselect).

  There are two types of escapes from this loop:
    * We need to (re)select a candidate x<j> for some reason.
    * We have a suspected termination condition --- optimality, infeasibility
      (or punt), unboundedness, or fatal error.

  The first action in the loop is to decide which way the incoming variable
  is moving, based on the variable's status and the sign of cbar<j>. Arguably
  this could move into dy_primalpivot.
*/
    for (xjndx = candxj ; do_pivots == TRUE ; xjndx = candxj)
    { 
#     ifdef DYLP_PARANOIA
      if (xjndx <= 0)
      { errmsg(386,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,"entering") ;
	return (dyrFATAL) ; }
#     endif
      indir = 0 ;
      xjstatus = dy_status[xjndx] ;
      cbarj = dy_cbar[xjndx] ;
      if (cbarj <= 0 && flgon(xjstatus,vstatNBLB|vstatSB|vstatNBFR))
      { indir = 1 ; }
      else
      if (cbarj > 0 && flgon(xjstatus,vstatNBUB|vstatSB|vstatNBFR))
      { indir = -1 ; }
#     ifdef DYLP_PARANOIA
      else
      { errmsg(1,rtnnme,__LINE__) ;
	return (dyrFATAL) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pricing >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n    (%s)%d: %s (%d), entering %s from %s, price = %g ... ",
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     consys_nme(dy_sys,'v',xjndx,TRUE,NULL),xjndx,
	     (indir < 0)?"decreasing":"increasing",dy_prtvstat(xjstatus),
	     cbarj/sqrt(dy_gamma[xjndx])) ; }
#     endif
/*
  Time to get down to business and attempt the pivot. dy_primalpivot does all
  the heavy lifting --- the actual pivot, plus updates of data structures and
  variables. Under normal conditions, we're looking for one of dyrOK (vanilla
  pivot) or dyrDEGEN (degenerate pivot), indicating a successful pivot and
  selection of a new candidate x<j>. dyrOPTIMAL or dyrPUNT are less common,
  and indicate a successful pivot but failure to find a new x<j>.  Other
  possibilities are dyrUNBOUND, dyrREQCHK (suspected accuracy problems),
  dyrMADPIV (pivot rejection), and various errors, including dyrLOSTPFEAS,
  dyrSINGULAR, dyrBSPACE, and dyrFATAL.
*/
      pivresult = dy_primalpivot(xjndx,indir,
				 &xindx,&outdir,&abarij,&delta,&candxj) ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase1 >= 4)
	dy_logpivot(pivresult,xjndx,indir,cbarj,xindx,outdir,abarij,delta) ;
#     endif
/*
  La Duenna makes sure the proprieties are observed and deals with any
  scandals. In the end, these are the cases to consider:

    * dyrOK: The pivot went through without a problem, or whatever went wrong
      has been dealt with transparently from La Duenna's point of view.  Call
      tweakp1obj to tweak the objective. We might discover we're feasible, in
      which case call it optimal and try for termination.  If tweakp1obj ends
      up resetting the reduced costs, escape to reselect, otherwise we'll go
      to the next pivot with the candidate selected by pseupdate.  There are
      a few circumstances where no candidate will be selected (the pivot was
      a nonbasic swing, or primalout returned REQCHK), and these also require
      a reselect.

      dyrOPTIMAL, dyrPUNT are much like dyrOK (i.e., the pivot went through)
      but we also think we're optimal (i.e., we can't find an entering
      variable). If tweakp1obj says we're feasible, try for optimal
      termination, but if it says we're still infeasible, try for infeasible
      termination. dyrPUNT indicates there are potential pivots flagged with
      the NOPIVOT qualifier, but it's pointless to try them again.

    * dyrRESELECT: Whatever happened requires that we select a new entering
      variable. Most commonly, this indicates that primalout couldn't find a
      numerically stable pivot. Less commonly, new pivot candidates have been
      released from the pivot rejection list in response to a punt.
      Occasionally, something more exotic has happened (e.g., the basis has
      been patched due to singularity). Escape to the outer loop to reselect.

    * dyrUNBOUND: We'll kick this back to the caller if preoptimality confirms
      the condition. We need additional constraints. Remember the unbounded
      column.

    * dyrSWING: The primal variables are moving too far, too fast. dyrSWING
      is a recommendation that we pop out of simplex and try to add
      constraints. But we won't do this unless we've completed a minimum
      number of basis changes (currently hardwired to 10).

    * anything else: These are problems too severe to handle here; they, too,
      get kicked back to the calling routine.  The main reason for
      enumerating return values here is to make sure we catch a code that was
      somehow overlooked --- it'll trigger an error message.

  In general, if we're going to continue the pivoting loop, neither do_pivots
  or lpretval should be changed. If we're just going to escape to the outer
  loop for reselection, we set do_pivots to FALSE and leave lpretval at dyrINV.
  To attempt to end phase I, lpretval must be properly set as well.
*/
      duennaresult = dy_duenna(pivresult,xjndx,xindx,candxj,-1) ;
      switch (duennaresult)
      { case dyrOK:
	case dyrOPTIMAL:
	case dyrPUNT:
	{ p1objresult = tweakp1obj(&reselect,candxj) ;
	  switch (p1objresult)
	  { case dyrINFEAS:
	    { if (duennaresult == dyrOK)
	      { if (candxj <= 0) reselect = TRUE ;
		if (reselect == TRUE || xindx == xjndx) do_pivots = FALSE ; }
	      else
	      { do_pivots = FALSE ;
		if (duennaresult == dyrPUNT)
		  lpretval = dyrPUNT ;
		else
		  lpretval = dyrINFEAS ; }
	      break ; }
	    case dyrOPTIMAL:
	    { lpretval = dyrOPTIMAL ;
	      do_pivots = FALSE ;
	      break ; }
	    case dyrFATAL:
	    { errmsg(318,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters,"tweak") ;
	      do_pivots = FALSE ;
	      lpretval = p1objresult ;
	      break ; }
	    default:
	    { errmsg(7,rtnnme,__LINE__,"tweakp1obj code",p1objresult) ;
	      do_pivots = FALSE ;
	      lpretval = p1objresult ;
	      break ; } }
	  break ; }
	case dyrRESELECT:
	{ do_pivots = FALSE ;
	  break ; }
	case dyrUNBOUND:
	{ do_pivots = FALSE ;
	  lpretval = duennaresult ;
	  dy_lp->ubnd.ndx = xjndx*indir ;
	  dy_lp->ubnd.ratio = 0 ;
	  break ; }
	case dyrSWING:
	{ if (dy_lp->basis.pivs >= 10)
	  { lpretval = duennaresult ;
	    do_pivots = FALSE ; }
	  break ; }
	case dyrACCCHK:
	case dyrSINGULAR:
	case dyrBSPACE:
	case dyrSTALLED:
	case dyrITERLIM:
	case dyrNUMERIC:
	case dyrFATAL:
	{ do_pivots = FALSE ;
	  lpretval = duennaresult ;
	  break ; }
	default:
	{ errmsg(7,rtnnme,__LINE__,"La Duenna return code",
		 (int) duennaresult) ;
	  do_pivots = FALSE ;
	  lpretval = dyrFATAL ;
	  break ; } } }
/*
  End of the pivoting loop. Why are we out here? The simplest case is that
  we just want to reselect --- head back to the top of the loop from here.
*/
    if (lpretval == dyrINV) continue ;
/*
  Do we think we're optimal (feasible), infeasible, or unbounded?  (I.e.,
  we're reporting one of dyrOPTIMAL, dyrINFEAS, dyrUNBOUND, or dyrPUNT.) If
  so, check that we have the objective function correct. If it's out-of-date,
  bring it up-to-date, and head back to the top of the loop to try and
  reselect an entering variable.

  dyrSWING isn't subject to the previous checks --- we just want to pop out and
  try to augment the constraint system.

  verifyp1obj recalculates the primal variables but does not refactor unless
  it concludes the objective is incorrect (thus there is the potential for
  loss of primal feasibility in preoptimality).
*/
    if (lpretval == dyrOPTIMAL || lpretval == dyrINFEAS ||
	lpretval == dyrUNBOUND || lpretval == dyrPUNT)
    { p1objresult = verifyp1obj() ;
      switch (p1objresult)
      { case dyrOK:
	{ break ; }
	case dyrRESELECT:
	{
#         ifndef DYLP_NDEBUG
	  if (dy_opts->print.phase1 >= 1)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		   "\n\tfalse termination (%s) due to inconsistent objective",
		   dy_prtdyret(lpretval)) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,
			" at (%s)%d; rebuilding P1 objective.",
		        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
#         endif
	  if (dy_initp1obj() == TRUE)
	  { lpretval = dyrINV ;
	    continue ; }
	  errmsg(318,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,"reinitialise") ;
	  p1objresult = dyrFATAL ;
	  break ; }
	default:
	{ p1objresult = dyrFATAL ;
	  break ; } }
      if (p1objresult == dyrFATAL) return (dyrFATAL) ; }
/*
  To get to here, we have a termination condition and the objective has been
  verified.  We can take useful action for dyrOPTIMAL, dyrUNBOUND, dyrINFEAS,
  and dyrPUNT; everything else is bounced directly back to the caller.

  For all, the first thing to do is call preoptimality to factor the basis
  and do accuracy and feasibility checks. If we get back dyrOK, dyrPATCHED,
  or dyrLOSTDFEAS, then we have primal feasibility and we're out of here
  regardless of the lpretval code.  (Dual feasibility or boundedness are
  nice, but not required). Force a return code of dyrOPTIMAL.

  If we get back dyrLOSTPFEAS, it's a little more complicated. If lpretval is
  dyrINFEAS or dyrUNBOUND, and the primal & dual feasibility status reported
  by preoptimality agrees, we can return.  If lpretval is dyrPUNT, there are
  variables on the pivot rejection list, but dy_dealWithPunt has already
  concluded none of them can be used. If preoptimality shows dual
  feasibility, none of them would help in any event, and we can return
  dyrINFEAS with a clear conscience. Otherwise return dyrPUNT.

  Anything else means that the feasibility status has taken an unexpected
  turn. The working hypothesis is that we need to be more careful in selecting
  pivots for both factoring and pivoting, which we do by tightening the
  current value and lower bound for the pivot selection parameters.

  There are a number of cases which imply unexpected gain/loss of primal or
  dual feasibility, based on lpretval and the results of preoptimality.
  Messages may be generated for these, if DYLP_NDEBUG allows it.
*/
    if (lpretval == dyrOPTIMAL || lpretval == dyrUNBOUND ||
	lpretval == dyrINFEAS || lpretval == dyrPUNT)
    { optcnt++ ;
#     ifndef DYLP_NDEBUG
      uxpfeas = FALSE ;
      uxnpfeas = FALSE ;
      uxdfeas = FALSE ;
      uxndfeas = FALSE ;
      tmpretval = lpretval ;
#     endif
      preopresult = preoptimality(lpretval,&checks) ;
      switch (preopresult)
      { case dyrOK:
	case dyrPATCHED:
	case dyrLOSTDFEAS:
	{
#         ifndef DYLP_NDEBUG
	  if (lpretval == dyrINFEAS || lpretval == dyrPUNT ||
	      lpretval == dyrUNBOUND)
	    uxpfeas = TRUE ;
	  if (preopresult == dyrLOSTDFEAS)
	  { if (lpretval == dyrINFEAS) uxndfeas = TRUE ; }
	  else
	  { if (lpretval == dyrPUNT || lpretval == dyrUNBOUND)
	      uxdfeas = TRUE ; }
#         endif
	  lpretval = dyrOPTIMAL ;
	  break ; }
	case dyrLOSTPFEAS:
	{ 
#	  ifndef DYLP_NDEBUG
	  if (lpretval == dyrOPTIMAL) uxnpfeas = TRUE ;
#	  endif
	  if (flgoff(checks,ladDUALFEAS))
	  { 
#           ifndef DYLP_NDEBUG
	    if (lpretval == dyrUNBOUND || lpretval == dyrPUNT)
	      uxdfeas = TRUE ;
#	    endif
	    if (!(lpretval == dyrINFEAS || lpretval == dyrPUNT))
	    { (void) dy_setpivparms(+1,+1) ;
	      lpretval = dyrINV ; }
	    else
	    { lpretval = dyrINFEAS ; } }
	  else
	  { if (lpretval == dyrUNBOUND || lpretval == dyrPUNT)
	    { /* no action required */ }
	    else
	    {
#	      ifndef DYLP_NDEBUG
	      if (lpretval == dyrINFEAS)
		uxndfeas = TRUE ;
#	      endif
	      (void) dy_setpivparms(+1,+1) ;
	      lpretval = dyrINV ; } }
	  break ; }
	default:
	{ lpretval = preopresult ;
	  break ; } }
/*
  Ok, we've done the necessary processing. There are two things to deal with
  now.  If DYLP_NDEBUG permits, check for warnings about unexpected primal/dual
  (in)feasibility.  If we're going to repeat (dyrINV), check that we haven't
  exceeded the repetition count.

  The warnings could be compacted, but this way we also get to see if the
  four variables have been set to a nonsense pattern.
*/
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase1 >= 2)
      { if (uxpfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tunexpected primal feasibility at iteration (%s)%d.",
		      dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxnpfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n\tunexpected loss of primal feasibility at iteration (%s)%d.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxdfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tunexpected dual feasibility at iteration (%s)%d.",
		      dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxndfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n\tunexpected loss of dual feasibility at iteration (%s)%d.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; } }
#     endif

      if (lpretval == dyrINV)
      { if (optcnt > 15)
	{ errmsg(387,rtnnme,dy_sys->nme,
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,optcnt) ;
	  lpretval = dyrFATAL ; }
#       ifndef DYLP_NDEBUG
	else
	{ if (dy_opts->print.phase1 >= 2)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		   "\n\t(%s)%d: false termination (%s); resuming pivoting.",
		   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		   dy_prtdyret(tmpretval)) ; } }
#       endif
      } } }
/*
  We've finished the outer loop. If we're optimal (and hence headed for phase
  II of the simplex) we'll reinstall the phase II objective, recalculate the
  objective value, duals, and reduced costs, and reset the reference frame.
  This makes primal I transparent to the rest of the code, w.r.t.  playing
  with the objective.

  But ... no sense in doing any of this if we're infeasible, unbounded, or
  otherwise. In this case, we'll be popping out of simplex entirely, and
  there's no sense trying to anticipate.
*/
  if (dy_clrpivrej(NULL) != TRUE) lpretval = dyrFATAL ;
  if (dy_lp->degen > 0) (void) dy_degenout(0) ;

  if (lpretval == dyrOPTIMAL)
  { if (dy_swapobjs(dyPRIMAL2) == FALSE)
    { errmsg(318,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"remove") ;
      return (dyrFATAL) ; }
    dy_calcduals() ;
    if (dy_calccbar() == FALSE)
    { errmsg(384,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters) ;
      return (dyrFATAL) ; }
    dy_lp->z = dy_calcobj() ;
    dy_pseinit() ; }

# ifndef DYLP_NDEBUG
  if (lpretval == dyrUNBOUND && dy_opts->print.phase1 >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    [%s] (%s)%d: system is unbounded.",
	        dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		dy_lp->tot.iters) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL)
  { dy_stats->p1.iters += dy_lp->p1.iters ;
    dy_stats->p1.pivs += dy_lp->p1.pivs ; }
# endif

  return (lpretval) ; }



static dyret_enum primal2 (void)

/*
  Phase 2 of the primal simplex. The basic action is a simple loop: pivot
  (dy_primalpivot), then see what went wrong (dy_duenna). As a side effect of
  the pivot, dy_primalpivot returns the index of the preferred candidate,
  xjcand, to enter on the next pivot. The work in this routine goes into
  trying to deal with problems of pivot rejection.
  
  If the preferred candidate is rejected, dy_primalin is called to select a
  new incoming variable. As long as an incoming variable can be selected, and
  nothing fatal goes wrong, the inner pivoting loop continues.

  Exit from the inner loop occurs due to optimality, unboundedness, a punt,
  or a fatal error (of which there are many). Optimality and unboundedness
  are self-explanatory. We run a preoptimality check (refactor plus primal and
  dual feasibility checks) to confirm. Loss of dual feasibility causes a
  return to the pivoting loop. Loss of primal feasibility causes a reversion
  to phase I.

  A punt indicates that no incoming variable could be selected but there are
  variables flagged with the NOPIVOT qualifier. This can be indicated
  by dy_primalpivot or by dy_primalin. If preoptimality indicates loss of
  dual feasibility (it ignores NOPIVOT qualifiers when doing the check), we'll
  relax the pivot selection tolerance and try again to select a pivot. The
  tolerance is progressively relaxed until a successful pivot occurs, at
  which point it snaps back to the original tolerance. If we relax to the
  bogus number tolerance before managing a successful pivot, we abort.

  Parameters: none

  Returns: most of the dyret_enum codes.
*/

{ dyret_enum lpretval,scanresult,pivresult,duennaresult,preopresult ;
  int startcol,scan,nextcol,candxj,
      xjndx,indir,xindx,outdir,optcnt ;
  double cbarj,abarij,delta ;
  bool do_pivots ;
  flags xjstatus,checks ;
  const char *rtnnme = "primal2" ;

# ifndef DYLP_NDEBUG
  bool uxnpfeas,uxdfeas,uxndfeas ;
  dyret_enum tmpretval ;
# endif

# ifdef DYLP_PARANOIA
  if (dy_lp->degen != 0)
  { errmsg(317,rtnnme,dy_sys->nme,dy_lp->degen) ;
    return (dyrFATAL) ; }
  if (dy_lp->p2.iters != 0)
  { errmsg(5,rtnnme,"phase 2 iteration count",dy_lp->p2.iters) ;
    return (dyrFATAL) ; }
# endif

  dy_lp->p2.pivs = 0 ;
  dy_lp->pivok = FALSE ;
  dy_lp->prev_pivok = FALSE ;
  lpretval = dyrINV ;
/*
  Open the main loop. The purpose of this outer loop is to allow easy
  recovery from false indications of optimality or unboundedness, as well as
  some errors involving pivot selection. All the action happens at the bottom
  of the loop, after we fall out of the pivoting loop. Initialise a few
  variables and get going.
*/
  if (dy_clrpivrej(NULL) != TRUE) return (dyrFATAL) ;
  optcnt = 0 ;
  scan = dy_opts->scan ;
  nextcol = 1 ;
  while (lpretval == dyrINV)
  { 
/*
  The first thing we need to do is a scan to select an incoming variable
  x<j>.  If primalin cannot find a candidate, we're optimal, or perhaps we've
  punted.  Optimality falls into the default case.  If we've punted, call
  dy_dealWithPunt to free up any potential candidates on the pivot rejection
  list and iterate to try again.
*/
    startcol = nextcol ;
    scanresult = dy_primalin(startcol,scan,&candxj,&nextcol) ;
    switch (scanresult)
    { case dyrOK:
      { do_pivots = TRUE ;
	break ; }
      case dyrPUNT:
      { scanresult = dy_dealWithPunt() ;
	if (scanresult == dyrRESELECT)
	{ continue ; }
	else
	{ do_pivots = FALSE ;
	  lpretval = scanresult ; }
	break ; }
      default:
      { do_pivots = FALSE ;
	lpretval = scanresult ;
	break ; } }
/*
  Open the loop that executes pivots. While we have a candidate x<j> for
  entry, we do a two-step: attempt the pivot with x<j> (dy_primalpivot),
  then check that everything went off ok (dy_duenna). As part of updating
  the PSE pricing information, a new x<j> will be selected, and we repeat
  the loop.

  There are two valid escapes from this loop:
    * We find optimality (no candidates for entry) or unboundedness (a
      direction of recession).
    * la Duenna aborts the problem (due to a problem it's not able to
      fix).

  The first action in the loop is to decide which way the incoming variable
  is moving, based on the variable's status and the sign of cbar<j>. Arguably
  this should move into dy_primalpivot.
*/
    for (xjndx = candxj ; do_pivots == TRUE ; xjndx = candxj)
    { 
#     ifdef DYLP_PARANOIA
      if (xjndx <= 0)
      { errmsg(386,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,"entering") ;
	return (dyrFATAL) ; }
#     endif
      indir = 0 ;
      xjstatus = dy_status[xjndx] ;
      cbarj = dy_cbar[xjndx] ;
      if (cbarj <= 0 && flgon(xjstatus,vstatNBLB|vstatSB|vstatNBFR))
      { indir = 1 ; }
      else
      if (cbarj > 0 && flgon(xjstatus,vstatNBUB|vstatSB|vstatNBFR))
      { indir = -1 ; }
#     ifdef DYLP_PARANOIA
      else
      { errmsg(1,rtnnme,__LINE__) ;
	return (dyrFATAL) ; }
#     endif
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pricing >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n    (%s)%d: %s (%d), entering %s from %s, price = %g ... ",
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     consys_nme(dy_sys,'v',xjndx,TRUE,NULL),xjndx,
	     (indir < 0)?"decreasing":"increasing",dy_prtvstat(xjstatus),
	     cbarj/sqrt(dy_gamma[xjndx])) ; }
#     endif
/*
  Time to get down to business and attempt the pivot. dy_primalpivot does all
  the heavy lifting --- the actual pivot, plus updates of data structures and
  variables. Under normal conditions, we're looking for one of dyrOK (vanilla
  pivot) or dyrDEGEN (degenerate pivot), indicating a successful pivot and
  selection of a new candidate x<j>. dyrOPTIMAL or dyrPUNT are less common,
  and indicate a successful pivot but failure to find a new x<j>.  Other
  possibilities are dyrUNBOUND, dyrREQCHK (suspected accuracy problems),
  dyrMADPIV (pivot rejection), and various errors, including dyrLOSTPFEAS,
  dyrSINGULAR, dyrBSPACE, and dyrFATAL.
*/
      pivresult = dy_primalpivot(xjndx,indir,
				 &xindx,&outdir,&abarij,&delta,&candxj) ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase2 >= 4)
	dy_logpivot(pivresult,xjndx,indir,cbarj,xindx,outdir,abarij,delta) ;
#     endif
/*
  La Duenna makes sure the proprieties are observed and deals with any
  scandals. In the end, there are four cases to distinguish:

    * dyrOK: The pivot went through without a problem, or whatever went wrong
      has been dealt with transparently from La Duenna's point of view. We'll
      go on to the next pivot using the candidate selected by dy_primalpivot.
      There are a few circumstances where no candidate will be selected (the
      pivot was a nonbasic swing, or primalout returned REQCHK), and these
      require we force a reselect.

    * dyrOPTIMAL, dyrPUNT, dyrUNBOUND: We'll want to escape the
      pivoting loop and run preoptimality. Depending on what it reports,
      we'll return to the caller or resume pivoting. Unbounded gets a case
      of its own so we can remember the unbounded column. Punt means there
      are potential pivots marked with the NOPIVOT qualifier.
      
    * dyrRESELECT: Whatever happened requires that we select a new entering
      variable. Most commonly, this indicates that primalout couldn't find a
      numerically stable pivot. Less commonly, new pivot candidates have been
      released from the pivot rejection list in response to a punt.
      Occasionally, something more exotic has happened (e.g., the basis has
      been patched due to singularity). Escape to the outer loop to reselect.

    * dyrSWING: The primal variables are moving too far, too fast. dyrSWING
      is a recommendation that we pop out of simplex and try to add
      constraints. But we won't do this unless we've completed a minimum
      number of basis changes (currently hardwired to 10).

    * Anything else: These are errors too severe to handle here; they get
      kicked back to the calling routine.

  The main reason for enumerating return values here is to make sure we catch
  a code that was somehow overlooked --- it'll trigger an error message.
  In general, if we're going to continue the pivoting loop, neither do_pivots
  or lpretval should be changed. If we're going to escape to the outer loop,
  lpretval must be properly set and do_pivots should be set to FALSE.
*/
      duennaresult = dy_duenna(pivresult,xjndx,xindx,candxj,-1) ;
      switch (duennaresult)
      { case dyrOK:
	{ if (candxj <= 0 || xjndx == xindx) do_pivots = FALSE ;
	  break ; }
	case dyrRESELECT:
	{ do_pivots = FALSE ;
	  break ; }
	case dyrOPTIMAL:
	case dyrPUNT:
	{ do_pivots = FALSE ;
	  lpretval = duennaresult ;
	  break ; }
	case dyrUNBOUND:
	{ do_pivots = FALSE ;
	  lpretval = duennaresult ;
	  dy_lp->ubnd.ndx = xjndx*indir ;
	  dy_lp->ubnd.ratio = 0 ;
	  break ; }
	case dyrSWING:
	{ if (dy_lp->basis.pivs >= 10)
	  { lpretval = duennaresult ;
	    do_pivots = FALSE ; }
	  break ; }
	case dyrACCCHK:
	case dyrSINGULAR:
	case dyrBSPACE:
	case dyrSTALLED:
	case dyrITERLIM:
	case dyrLOSTPFEAS:
	case dyrNUMERIC:
	case dyrFATAL:
	{ do_pivots = FALSE ;
	  lpretval = duennaresult ;
	  break ; }
	default:
	{ errmsg(7,rtnnme,__LINE__,"La Duenna return code",
		 (int) duennaresult) ;
	  do_pivots = FALSE ;
	  lpretval = dyrFATAL ;
	  break ; } } }
/*
  Why are we here? The simplest case is that we just need to select a new
  entering variable --- head back to the top of the loop from here.
*/
    if (lpretval == dyrINV) continue ;
/*
  Do we think we're optimal or unbounded? (I.e., we're reporting dyrOPTIMAL,
  dyrPUNT, or dyrUNBOUND.) If so, the first thing to do is call preoptimality
  to refactor and do accuracy and feasibility checks. If we have primal and
  dual feasibility (dyrOK, dyrPATCHED), we can return dyrOPTIMAL with a clear
  conscience.  If preoptimality reports primal and dual feasibility for
  dyrUNBOUND or dyrPUNT, well, unexpected primal feasibility is always a
  pleasure, and we'll play along.

  dyrLOSTDFEAS is what we'd expect for dyrUNBOUND, so this is also a clean
  termination.

  dyrLOSTDFEAS with dyrPUNT says that there are variables on the pivot
  rejection list, but dy_dealWithPunt has already concluded nothing can be
  done to make them useable. Return dyrPUNT.

  Anything else means that the feasibility status has taken an unexpected
  turn. The working hypothesis is that we need to be more careful in
  selecting pivots for both factoring and pivoting, which we do by tightening
  the current value and lower bound for the pivot selection parameters.  If
  we've lost primal feasibility, revert to phase I (dy_primal will take care
  of tightening the pivot selection parameters). If we've only lost dual
  feasibility, head back to the top of the loop to resume pivoting.
*/
    if (lpretval == dyrOPTIMAL ||
	lpretval == dyrUNBOUND || lpretval == dyrPUNT)
    { optcnt++ ;
#     ifndef DYLP_NDEBUG
      uxnpfeas = FALSE ;
      uxdfeas = FALSE ;
      uxndfeas = FALSE ;
      tmpretval = lpretval ;
#     endif
      preopresult = preoptimality(lpretval,&checks) ;
      switch (preopresult)
      { case dyrOK:
	case dyrPATCHED:
	{ 
#	  ifndef DYLP_NDEBUG
	  if (lpretval == dyrUNBOUND || lpretval == dyrPUNT) uxdfeas = TRUE ;
#	  endif
	  lpretval = dyrOPTIMAL ;
	  break ; }
	case dyrLOSTPFEAS:
	{ 
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.phase2 >= 1)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
			"\n  (%s)%d: lost primal feasibility.",
			dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
	  uxnpfeas = TRUE ;
	  if (flgon(checks,ladDUALFEAS))
	  { if (lpretval == dyrOPTIMAL) uxndfeas = TRUE ; }
	  else
	  { if (lpretval != dyrOPTIMAL) uxdfeas = TRUE ; }
#	  endif
	  lpretval = dyrLOSTPFEAS ;
	  break ; }
	case dyrLOSTDFEAS:
	{ if (lpretval == dyrUNBOUND || lpretval == dyrPUNT)
	  { /* no action required */ }
	  else
	  if (lpretval == dyrOPTIMAL)
	  { lpretval = dyrINV ;
	    (void) dy_setpivparms(+1,+1) ;
#	    ifndef DYLP_NDEBUG
	    uxndfeas = TRUE ;
#	    endif
	  }
	  break ; }
	default:
	{ lpretval = preopresult ;
	  break ; } }
/*
  Ok, we've done the necessary processing. There are two things to deal with
  now.  If DYLP_NDEBUG permits, check for warnings about unexpected primal/dual
  (in)feasibility.  If we're going to repeat (dyrINV), check that we haven't
  exceeded the repetition count.

  The warnings could be compacted, but this way we also get to see if the
  four variables have been set to a nonsense pattern.
*/
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase2 >= 2)
      { if (uxnpfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n\tunexpected loss of primal feasibility at iteration (%s)%d.",
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxdfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		 "\n\tunexpected dual feasibility at iteration (%s)%d.",
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxndfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n\tunexpected loss of dual feasibility at iteration (%s)%d.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; } }
#     endif

      if (lpretval == dyrINV)
      { if (optcnt > 15)
	{ errmsg(387,rtnnme,dy_sys->nme,
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,optcnt) ;
	  lpretval = dyrFATAL ; }
#       ifndef DYLP_NDEBUG
	else
	{ if (dy_opts->print.phase2 >= 2)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		   "\n\tfalse termination (%s) at (%s)%d; resuming pivoting.",
		   dy_prtdyret(tmpretval),dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters) ; } }
#       endif
      } } }
/*
  We've finished the outer loop. Clean up before we leave.
*/
  if (dy_clrpivrej(NULL) != TRUE) lpretval = dyrFATAL ;
  if (dy_lp->degen > 0) (void) dy_degenout(0) ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.phase2 >= 2)
  { if (lpretval == dyrUNBOUND)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: system %s is unbounded.",
		  rtnnme,dy_sys->nme) ; }
    else
    if (lpretval == dyrSWING)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: system %s is pseudo-unbounded.",
		  rtnnme,dy_sys->nme) ; } }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL)
  { dy_stats->p2.iters += dy_lp->p2.iters ;
    dy_stats->p2.pivs += dy_lp->p2.pivs ; }
# endif

  return (lpretval) ; }



lpret_enum dy_primal (void)

/*
  This is the driver routine for the primal simplex algorithm. It expects that
  the problem will come in with an initial basis. The algorithm is two-phase,
  running phase I to attain feasibility, if needed, then running phase II to
  obtain optimality.

  Loss of feasibility in phase II will cause reentry into phase I. Presently
  there's a hardcoded limit of 10 occurrences before the code aborts.

  Parameters: none

  Returns: any of lpret_enum codes.
*/

{ lpret_enum retval ;
  dyret_enum dyret ;
  int lostfeascnt ;
  const char *rtnnme = "dy_primal" ;

  retval = lpINV ;
  dy_lp->lpret = lpINV ;
  (void) dy_setpivparms(-100,-100) ;
  (void) dy_setpivparms(+1,+1) ;
  dy_lp->basis.pivs = 0 ;
/*
  We open a loop here to allow for the possibility of returning to phase I if
  we loose feasibility in phase II. There's a limit on how many times we'll
  allow this, though.

  Run phase I if we need it.
*/
  for (lostfeascnt = 0 ; lostfeascnt < 10  ; lostfeascnt++)
  { if (dy_lp->infeas > 0.0)
    { dy_lp->phase = dyPRIMAL1 ;
      dy_lp->p1.iters = 0 ;
      dyret = primal1() ;
#     ifndef DYLP_NDEBUG
      if ((dy_opts->print.phase1 >= 2) ||
	  (dy_opts->print.phase1 >= 1 && dyret != dyrOPTIMAL))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  (%s)%d: primal phase I ended, %d pivots, status %s.",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		    dy_lp->p1.pivs,dy_prtdyret(dyret)) ; }
      if (dy_opts->print.major >= 1 && dyret == dyrOPTIMAL)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\n%s (%s): entering phase %s, iter %d.",
		    "dylp",dy_sys->nme,dy_prtlpphase(dyPRIMAL2,FALSE),
		    dy_lp->tot.iters) ; }
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL && dyret == dyrOPTIMAL)
      { dy_stats->phasecnts[dyPRIMAL2]++ ; }
#     endif
    }
    else
    { dyret = dyrOPTIMAL ; }
    if (dyret != dyrOPTIMAL) break ;
/*
  Phase I succeeded, so we can go on to phase II. Call primal2 to do the work.
*/
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.phase2 >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: entering primal phase II, z = %g",
		  rtnnme,dy_lp->z) ;
      if (dy_opts->print.phase2 >= 4)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", dual active yb = %g",dy_calcdualobj()) ; }
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
    dy_lp->phase = dyPRIMAL2 ;
    dy_lp->simplex.active = dyPRIMAL2 ;
    dy_lp->p2.iters = 0 ;
    dyret = primal2() ;
/*
  What do we have? Optimality is best, but unboundedness is also possible.
  Loss of feasibility takes us back to phase I again. Anything else is an
  error.

  Only loss of feasibility iterates the loop. Tighten up the minimum pivot
  selection parameters, perhaps we can avoid this next time around. Force
  dyPRIMAL1 as the phase for the benefit of paranoid checks on status inside
  forcesuperbasic.
*/
    if (dyret == dyrLOSTPFEAS)
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.phase2 >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  (%s)%d: lost feasibility by %g after %d pivots; ",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		    dy_lp->infeas,dy_lp->p2.pivs) ;
	if (lostfeascnt+1 < 10)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"returning to phase I for try %d.",
		      lostfeascnt+2) ;
	else
	  dyio_outfmt(dy_logchn,dy_gtxecho,"aborting after %d tries.",
		      lostfeascnt+1) ; }
		
#     endif
      if (lostfeascnt+1 < 10)
      { dy_lp->phase = dyPRIMAL1 ;
	if (forcesuperbasic() == FALSE)
	{ errmsg(391,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters) ;
	  dyret = dyrFATAL ;
	  break ; }
	(void) dy_setpivparms(0,+1) ; }
      continue ; }
    else
    {
#     ifndef DYLP_NDEBUG
      if ((dy_opts->print.phase2 >= 2) ||
	  (dy_opts->print.phase2 >= 1 && dyret != dyrOPTIMAL))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n  (%s)%d: primal phase II ended, %d pivots, status %s.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       dy_lp->p2.pivs,dy_prtdyret(dyret)) ; }
#     endif
      break ; } }
/*
  That's it. Why are we here?
    * If all went well, we  have phase II optimality.
    * We could have infeasibility (phase I) or unboundedness (phase I or II).
    * We could have exceeded the occurrence limit for loss of feasibility.
    * We may have reached an iteration limit.
    * Something major went wrong.
  Take a minute to translate the dyret_enum code into an lpret_enum code and
  the do a bunch of printing.
*/
  retval = dyret2lpret(dyret) ;
# ifndef DYLP_NDEBUG
  if (retval == lpOPTIMAL || retval == lpUNBOUNDED || retval == lpINFEAS)
  { if ((dy_lp->phase == dyPRIMAL1 && dy_opts->print.phase1 >= 2) ||
	(dy_lp->phase == dyPRIMAL2 && dy_opts->print.phase2 >= 2))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    (%s)%d: ",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%s ended, %d pivots, ",
		  dy_prtlpphase(dy_lp->phase,FALSE),dy_lp->tot.pivs) ;
      if (retval == lpOPTIMAL)
	dyio_outfmt(dy_logchn,dy_gtxecho,"z<opt> = %g.",dy_lp->z) ;
      else
      if (retval == lpINFEAS)
	dyio_outfmt(dy_logchn,dy_gtxecho,"infeas = %g.",dy_lp->infeas) ;
      else
	dyio_outfmt(dy_logchn,dy_gtxecho,"unbounded.") ; } }
  else
  if (retval == lpLOSTFEAS)
  { if (dy_opts->print.phase2 >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	    "\n  (%s)%d: primal simplex aborted; lost feasibility %d times.",
	    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	    lostfeascnt-1) ; } }
  else
  if (retval == lpITERLIM)
  { if ((dy_lp->phase == dyPRIMAL1 && dy_opts->print.phase1 >= 1) ||
	(dy_lp->phase == dyPRIMAL2 && dy_opts->print.phase2 >= 1))
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n  (%s)%d: primal simplex terminated; iteration limit (%d).",
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     dy_opts->iterlim) ; } }
  else
  { if ((dy_lp->phase == dyPRIMAL1 && dy_opts->print.phase1 >= 1) ||
	(dy_lp->phase == dyPRIMAL2 && dy_opts->print.phase2 >= 1))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  (%s)%d: ",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"%s failed, status %s after %d pivots.",
		  dy_prtlpphase(dy_lp->phase,FALSE),dy_prtdyret(dyret),
		  dy_lp->tot.pivs) ; } }
/*
  For curiousity's sake, try the traditional certificate of optimality:
  dual objective == yb == cx == primal objective.
*/
  if (dy_opts->print.phase2 >= 4 && retval == lpOPTIMAL)
  { double dualobj,primalobj ;

    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: comparing dual and primal objectives.",rtnnme) ;
    dualobj = dy_calcdualobj() ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tdual objective yb = %g.",dualobj) ;
    
    primalobj = dy_calcobj() ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\tprimal objective cx = %g.",primalobj) ;

    if (!withintol(dualobj,primalobj,dy_tols->dchk))
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tWHOOPS! yb - cx = %g - %g = %g > %g.",
		  dualobj,primalobj,dualobj-primalobj,dy_tols->dchk) ;
    else
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tobjectives match.") ; }
# endif

/*
  A last bit of substance --- set the return code & we're out of here.
*/

  dy_lp->lpret = retval ;

  return (retval) ; }

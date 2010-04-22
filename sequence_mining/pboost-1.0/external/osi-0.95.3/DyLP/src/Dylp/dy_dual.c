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
  This file contains the routines specific to dylp's dual simplex algorithm.
  It's a pretty straightforward implementation of dual phase II, with dual
  steepest edge (DSE) pricing.

  It's difficult to give good, clear explanations for this code from the
  viewpoint of running a simplex algorithm on the dual problem. The root of
  the difficulty lies with the handling of bounded primal variables. Suppose
  we have m constraints and n variables in the problem. If you follow the
  usual logic of dual simplex, you see that what's happening is something
  like this:  When a bounded variable is tight against an upper or lower
  bound constraint, we really need the dual variable that should be
  associated with that constraint, because it's nonzero, basic, and
  contributes a term y<m+j>(-l<j>) or y<m+n+j>u<j> to the objective function.
  But the constraints (either -x<j> <= -l<j> or x<j> <= u<j>) are not
  explicitly present in the constraint system, so we'll never see y<m+j> or
  y<m+n+j> in y = c<B>inv(B).

  If you look in a text, they all assume 0 <= x <= +inf when explaining the
  dual algorithm. Hence the contribution for y<m+j> is 0, all y<m+n+j> don't
  exist, and the whole problem is sort of conveniently swept under the rug.
  Just try and find a detailed explanation for bounded variables, l <= x <=
  u, where l != 0 and u != +inf.

  So, why does revised dual simplex work at all running off the primal
  tableaux?

  The first point to make is that the dual objective will, in fact, be
  wrong.  The necessary dual variables are missing from y = c<B>inv(B), which
  contains only those duals associated with the constraints present in the
  primal matrix, A. And the necessary lower and upper bounds are missing from
  the rhs vector, b. The values of the missing duals can be extracted from
  the reduced costs, if desired, since the reduced costs cbar<j> are the
  negative of the values of (all) the basic dual variables, and from this the
  objective can be calculated correctly.

  Anyway, if you slog through the algebra, what you eventually find is that
  the dual surplus variable, call it sigma<j>, associated with x<j> through
  the dual constraint dot(y,a<j>) - sigma<j> = c<j>,  gets pressed into
  service as a surrogate for y<m+j> and y<m+n+j>. You could say that the dual
  simplex algorithm decides which one it'll be, as part of the rules of the
  algorithm. In the case of sigma<j> acting as surrogate for y<m+n+j>,
  there's the added complication that the sign is wrong in the primal
  tableaux, so that where we would have y<m+n+j> >= 0 if the upper bound
  constraint was explicit, we have sigma<j> <= 0 when it's handled
  implicitly. This works out rather neatly (if coincidentally), in the sense
  that (for a max primal) we want a positive reduced cost at optimality for a
  variable at upper bound.

  Generally speaking, one gets a negative dual for variables out at an upper
  bound. The other place where this can occur is when a range constraint
  results in an upper-bounded slack.

  To really understand dual simplex with bounded variables as simplex working
  on the dual problem, you should work the math. When you're done, you'll
  appreciate why texts never get into details. And with that in mind, I'm
  going to comment the code in terms of the results that are required in the
  primal problem.

  One last thing -- the usual comment in texts is that you read the dual
  variables as the negative tranpose of the primal variables. This assumes
  that the primal/dual pair is

    max cx			min yb
	Ax <= b			    yA >= c
    l <= x <= u			    y >= 0
 
  Here, the primal is min cx. The proper way to handle the conversion is to
  ask what would happen if the primal was max -cx. The result is that the
  primal reduced costs are the correct values of the dual variables (i.e.,
  no sign inversions required). When we're hunting for the leaving dual
  variable, we're working on finding the minimum delta_y using
    y<k> = ybar<k> - delta_y<i>(-beta<i>N)
  where ybar<i> is a reduced cost and beta<i> is row k of inv(B). The value
  of the entering dual will be -ybar<j>/abar<i,j> = cbar<j>. The comments
  in dy_dualin sort of gloss over why we're looking for the minimum delta_y
  and concentrate on explaining the operation in terms of getting the correct
  sign for cbar<j>.

  In the context of dylp, the dual algorithm is strictly subordinate, used to
  reoptimise after the addition of constraints. More, the assumption is that
  if the primal simplex uncovers unboundedness, we'll add constraints there
  and return to primal phase I if necessary. The dual algorithm therefore
  assumes that a primal optimal but infeasible solution is available. This
  simplifies things, in that we don't need a stage I for the dual.

  Anti-degeneracy in the dual isn't as strong as in the primal --- only the
  `anti-degen lite' heuristic is implemented, using alignment of the dual
  constraints (columns of the constraint matrix) with the rhs vector.

  When numeric ill-conditioning is uncovered, we attempt to deal with it by
  boosting the minimum pivot tolerances. This can happen if groombasis has
  to correct a major status error, or if dual2 thinks it's reached optimality
  and then unexpectedly looses dual or primal feasibility. We keep trying
  with the dual if only primal feasibility is lost. If we've lost dual
  feasibility, we revert to the primal simplex, as there's no dual phase I.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_dual.c	4.7	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_dual.c 94 2006-06-29 23:06:51Z lou $" ;



static dyret_enum preoptimality (dyret_enum lpretval, flags *result)

/*
  This routine does the prep work so that we can have confidence in a report
  of optimality by the dual simplex routines. It clears the pivot reject
  list, refactors, recalculates the primal and dual variables, and performs
  primal and dual accuracy and feasibility checks.

  No heroic measures are taken, on loss of primal or dual feasibility, as
  the overlying algorithms will deal with it.

  Parameters:
    lpretval:	the lp return code assigned by dual2
    result:	(o) loaded with the result flags from dy_accchk

  Returns:
    dyrOK:	if all goes smoothly
    dyrPATCHED: if the only bump is that the basis was patched by dy_factor
    dyrLOSTDFEAS: if the dual feasibility check by dy_accchk fails
    dyrLOSTPFEAS: if the primal feasibility check by dy_accchk fails.
    dyrFATAL:   if anything else goes wrong

  Also can relay error codes from dy_accchk (dyrACCCHK, and various basis
  factoring errors). Loss of dual feasibility dominates loss of primal
  feasibility.
*/

{ flags checkflags ;
  dyret_enum retval ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "preoptimality" ;
# endif

# ifdef PARANOIA
  if (!(lpretval == dyrOPTIMAL || lpretval == dyrUNBOUND ||
	lpretval == dyrPUNT))
  { errmsg(4,rtnnme,"lp return code",dy_prtdyret(lpretval)) ;
    return (dyrFATAL) ; }
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.dual >= 4)
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: validating optimality at iteration %d.",
	        rtnnme,dy_lp->tot.iters) ;
# endif
/*
  A little prep work. If the dual simplex has returned lpUNBOUNDED (primal
  infeasible) or lpPUNT, we still want the primal feasibility check to
  calculate total infeasibility. We don't want heroic measures taken to try
  and deal with the problem. Don't ask for a refactor if this looks like
  uncomplicated optimality.
*/
  *result = 0 ;
  checkflags = 0 ;
  setflg(checkflags,ladFACTOR|ladPRIMALCHK|ladPRIMFEAS|ladPFQUIET|
		    ladDUALCHK|ladDUALFEAS|ladDFQUIET) ;
  if (lpretval == dyrOPTIMAL && dy_lp->basis.etas == 0)
    clrflg(checkflags,ladFACTOR) ;
/*
  Start with the easy stuff -- clear the pivot reject list and back out any
  restricted subproblem. If dualdegenout notes accuracy loss, request a
  refactor.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.dual >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\tclearing pivot rejection machinery ... ") ; }
# endif
  if (dy_clrpivrej(NULL) != TRUE) return (dyrFATAL) ;
  if (dy_lp->degen > 0)
  { if (dy_dualdegenout(0) == dyrREQCHK) setflg(checkflags,ladFACTOR) ; }
/*
  And now the accuracy checks.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.dual >= 4)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"done.\n\t%schecking accuracy ... ",
	        flgon(checkflags,ladFACTOR)?"refactoring and ":"") ; }
# endif
  retval = dy_accchk(&checkflags) ;
  *result = checkflags ;
  if (!(retval == dyrOK || retval == dyrPATCHED))
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.dual >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%sfailed.",
		  (dy_opts->print.dual >= 5)?"\n\t":" ") ; }
#   endif
    return (retval) ; }
  else
  if (flgon(checkflags,ladPRIMALCHK|ladDUALCHK))
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.dual >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%sfailed",
		  (dy_opts->print.dual >= 5)?"\n\t":" ") ;
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
    if (dy_opts->print.dual >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"lost") ;
      if (flgon(checkflags,ladPRIMFEAS))
	dyio_outfmt(dy_logchn,dy_gtxecho," primal") ;
      if (flgon(checkflags,ladDUALFEAS))
	dyio_outfmt(dy_logchn,dy_gtxecho," dual") ;
      dyio_outfmt(dy_logchn,dy_gtxecho," feasibility.") ; }
#   endif
    if (flgon(checkflags,ladDUALFEAS))
      retval = dyrLOSTDFEAS ;
    else
      retval = dyrLOSTPFEAS ; }
# ifndef DYLP_NDEBUG
  else
  { if (dy_opts->print.dual >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%s%s.",
		  (dy_opts->print.dual >= 5)?"\n\t":" ",
		  (retval == dyrOK)?"done":"patched") ; } }
# endif

  return (retval) ; }
  


static dyret_enum dual2 (void)

/*
  Phase 2 of the dual simplex. The basic action is a simple loop: pivot
  (dy_dualpivot), then see what went wrong (dy_duenna). As a side effect of
  the pivot, dy_dualpivot returns the index of the preferred candidate,
  xicand, to leave on the next pivot. The work in this routine goes into
  trying to deal with problems of pivot rejection.
  
  If the preferred candidate is rejected, dy_dualout is called to select a
  new leaving variable. As long as a leaving variable can be selected, and
  nothing fatal goes wrong, the inner pivoting loop continues.

  Exit from the inner loop occurs due to optimality (which occurs when we
  achieve primal feasibility), (dual) unboundedness (primal infeasibility), a
  punt, or a fatal error (of which there are many). Optimality and
  unboundedness are self-explanatory. We run a preoptimality check (refactor
  plus primal and dual feasibility checks) to confirm. Loss of primal
  feasibility causes a return to the pivoting loop. Loss of dual feasibility
  causes a reversion to the primal simplex.

  A punt indicates that no leaving variable could be selected but there are
  infeasible variables flagged with the NOPIVOT qualifier. This can be indicated
  by dy_dualpivot or by dy_dualout. If preoptimality indicates loss of
  primal feasibility (it ignores NOPIVOT qualifiers when doing the check), we'll
  relax the pivot selection tolerance and try again to select a pivot. The
  tolerance is progressively relaxed until a successful pivot occurs, at
  which point it snaps back to the original tolerance. If we relax to the
  bogus number tolerance before managing a successful pivot, we abort.

  Parameters: none

  Returns: appropriate dual lpret_enum code.
*/

{ int candxi,xindx,outdir,xjndx,indir,optcnt ;
  double cbarj,abarij,delta ;
  flags xistatus,checks ;
  bool do_pivots ;
  dyret_enum lpretval,outresult,pivresult,duennaresult,preopresult ;

  const int successiveDinf = 40 ;

  const char *rtnnme = "dual2" ;

# ifndef DYLP_NDEBUG
  int xipos ;
  bool uxpfeas,uxnpfeas,uxndfeas ;
  dyret_enum tmpretval ;
# endif

# ifdef PARANOIA
  if (dy_lp->degen != 0)
  { errmsg(317,rtnnme,dy_sys->nme,dy_lp->degen) ;
    return (dyrFATAL) ; }
  if (dy_lp->d2.iters != 0)
  { errmsg(5,rtnnme,"dual iteration count",dy_lp->d2.iters) ;
    return (dyrFATAL) ; }
# endif

  dy_lp->d2.pivs = 0 ;
  dy_lp->pivok = FALSE ;
  dy_lp->prev_pivok = FALSE ;
  lpretval = dyrINV ;
  dy_lp->basis.dinf = 0 ;
/*
  Do a little initialisation, then open the outer loop. It's purpose is to
  recover from false terminations (optimality or unboundedness which
  disappeared on refactoring) and pivot selection problems. All the action
  occurs at the bottom of the loop, after we fall out of the inner pivoting
  loop.
*/
  if (dy_clrpivrej(NULL) != TRUE) return(dyrFATAL) ;
  optcnt = 0 ;
  while (lpretval == dyrINV)
  {
/*
  Call dy_dualout to locate a leaving variable x<i> using DSE pricing. From a
  primal viewpoint, we're looking for a variable x<i> with the greatest
  normalised infeasibility. If there are no infeasible variables, we have
  optimality.  The other possibility is that all infeasible variables are on
  the pivot reject list; this is indicated by a return value of dyrPUNT.
  Optimality falls into the default case. Call dy_dealWithPunt to free up any
  potential candidates on the pivot rejection list and iterate to try again.
*/
    outresult = dy_dualout(&candxi) ;
    switch (outresult)
    { case dyrOK:
      { do_pivots = TRUE ;
	break ; }
      case dyrPUNT:
      { outresult = dy_dealWithPunt() ;
	if (outresult == dyrRESELECT)
	{ continue ; }
	else
	{ do_pivots = FALSE ;
	  lpretval = outresult ; }
	break ; }
      default:
      { do_pivots = FALSE ;
	lpretval = outresult ;
	break ; } }
#   ifdef PARANOIA
    if (candxi <= 0 && outresult == dyrOK)
    { dyio_outfmt(dy_logchn,TRUE,
		  "\ndualout: %s(%d) INVALID LEAVING VAR = %d, outresult = %s",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  xindx,dy_prtdyret(outresult)) ;  }
#   endif
/*
  Open the loop that executes pivots. While we have a candidate x<i> to
  leave, we do a two-step: attempt the pivot with x<i> (dy_dualpivot),
  then check that everything went off ok (dy_duenna). As part of updating
  the DSE pricing information, a new x<i> will be selected, and we repeat
  the loop.

  There are two valid escapes from this loop:
    * We find optimality (no candidates to leave) or unboundedness (a
      direction of recession).
    * la Duenna aborts the problem (due to a problem it's not able to
      fix).

  The first action in the loop is to decide which way the leaving variable
  is moving, based on its status. Arguably this should move into dy_dualpivot.
*/
    for (xindx = candxi ; do_pivots == TRUE ; xindx = candxi)
    { 
#     ifdef PARANOIA
      if (xindx <= 0)
      { dyio_outfmt(dy_logchn,TRUE,
		    "\nloop: %s(%d) INVALID LEAVING VAR = %d, outresult = %s",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		    xindx,dy_prtdyret(outresult)) ; 
	errmsg(386,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,"leaving") ;
	return (dyrFATAL) ; }
#     endif
      xistatus = dy_status[xindx] ;
      if (flgon(xistatus,vstatBLLB))
	outdir = 1 ;
      else
	outdir = -1 ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pricing >= 2)
      { xipos = dy_var2basis[xindx] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n    (%s)%d: %s (%d) = %g, ",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		    consys_nme(dy_sys,'v',xindx,TRUE,NULL),xindx,
		    dy_xbasic[xipos]) ;
	if (outdir < 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		 "decreasing and leaving at ub = %g, price = %g.",
		 dy_sys->vub[xindx],
		 (dy_xbasic[xipos]-dy_sys->vub[xindx])/sqrt(dy_rho[xipos])) ; }
	else
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	       "increasing and leaving at lb = %g, price = %g.",
	       dy_sys->vlb[xindx],
	       (dy_sys->vlb[xindx]-dy_xbasic[xipos])/sqrt(dy_rho[xipos])) ; } }
#     endif
/*
  So far so good. Call dy_dualpivot to do the heavy lifting. It will select
  an incoming variable x<j>, pivot the basis representation, update the
  primal and dual variables, the DSE information, and related data
  structures. As a side effect of the DSE updates, a new leaving variable
  x<i> should be selected. We're hoping for a return code of dyrOK or
  dyrDEGEN, indicating a successful (perhaps degenerate) pivot and selection
  of a new leaving variable. dyrOPTIMAL or dyrPUNT indicate a successful
  pivot but failure to select a new x<i>. Other possibilities are dyrUNBOUND
  (primal infeasible), dyrREQCHK (suspected accuracy problems), dyrMADPIV
  (numerically unstable pivot),  and various errors, including dyrLOSTDFEAS,
  dyrSINGULAR, dyrBSPACE, and dyrFATAL.
*/
      xjndx = -1 ;
      pivresult = dy_dualpivot(xindx,outdir,
			       &xjndx,&indir,&cbarj,&abarij,&delta,&candxi) ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.dual >= 4)
	dy_logpivot(pivresult,xjndx,indir,cbarj,xindx,outdir,abarij,delta) ;
#     endif
/*
  La Duenna makes sure the proprieties are observed and deals with any
  scandals. In the end, there are four cases to distinguish:

    * dyrOK: The pivot went through without a problem, or whatever went
      wrong has been fixed transparently from La Duenna's point of view.
      It may be that no candidate was selected because dualin returned REQCHK,
      and in this case we need to force a reselect.

    * dyrOPTIMAL, dyrUNBOUND, dyrPUNT: We'll want to end the pivoting loop
      and run preoptimality. Depending on what it has to say, we'll either
      return or resume pivoting. Unbounded gets its own case so we can
      remember the unbounded row. Punt means that there are potential pivots,
      but they're flagged with the NOPIVOT qualifier.

    * dyrRESELECT: Whatever happened requires that we select a new leaving
      variable. Most commonly, this indicates that dualin couldn't find a
      numerically stable pivot. Less commonly, new pivot candidates have been
      released from the pivot rejection list in response to a punt.
      Occasionally, something more exotic has happened (e.g., the basis has
      been patched due to singularity). Escape to the outer loop to reselect.

    * dyrSWING: The primal variables are moving too far, too fast. dyrSWING
      is a recommendation that we pop out of simplex and try to add
      constraints. But we won't do this unless we've completed a minimum
      number of basis changes (currently hardwired to 10).

    * dyrLOSTDFEAS: The dual simplex has a reasonable tolerance for slight
      loss of feasibility. But if the infeasibility persists through N
      successive refactors, better give it up. Currently, N is hardwired to
      ten.

    * Anything else: errors too severe to handle here. Kick back and let the
      caller deal with it.
  
  The reason for enumerating return values here is to make sure that we haven't
  overlooked some return value. This way, it triggers an error.
*/
	duennaresult = dy_duenna(pivresult,xjndx,xindx,-1,candxi) ;
	switch (duennaresult)
	{ case dyrOK:
	  { if (candxi <= 0) do_pivots = FALSE ;
	    break ; }
	  case dyrRESELECT:
	  { do_pivots = FALSE ;
	    break ; }
	  case dyrOPTIMAL:
	  case dyrPUNT:
	  { lpretval = duennaresult ;
	    do_pivots = FALSE ;
	    break ; }
	  case dyrUNBOUND:
	  { lpretval = dyrUNBOUND ;
	    dy_lp->ubnd.ndx = xindx*outdir ;
	    dy_lp->ubnd.ratio = 0 ;
	    do_pivots = FALSE ;
	    break ; }
	  case dyrSWING:
	  { if (dy_lp->basis.pivs >= 10)
	    { lpretval = duennaresult ;
	      do_pivots = FALSE ; }
	    else
	    if (candxi <= 0)
	    { do_pivots = FALSE ; }
	    break ; }
	  case dyrLOSTDFEAS:
	  {
#	    ifndef DYLP_NDEBUG
	    if (dy_opts->print.dual >= 3 &&
		dy_lp->basis.dinf < successiveDinf)
	    { dyio_outfmt(dy_logchn,dy_gtxecho,
		     "\n(%s)%d: dual infeasible for %d successive refactors.",
		     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		     dy_lp->basis.dinf) ; }
	    if (dy_opts->print.dual >= 1 &&
	    	dy_lp->basis.dinf >= successiveDinf)
	    { dyio_outfmt(dy_logchn,dy_gtxecho,
			  "\n  (%s)%d: dual infeasibility for %d successive",
			  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
			  successiveDinf) ;
	      dyio_outfmt(dy_logchn,dy_gtxecho," refactors; aborting.") ; }
#	    endif
	    if (dy_lp->basis.dinf >= successiveDinf)
	    { lpretval = duennaresult ;
	      do_pivots = FALSE ; }
	    else
	    if (candxi <= 0)
	    { do_pivots = FALSE ; }
	    break ; }
	  case dyrACCCHK:
	  case dyrSINGULAR:
	  case dyrBSPACE:
	  case dyrSTALLED:
	  case dyrITERLIM:
	  case dyrNUMERIC:
	  case dyrFATAL:
	  { lpretval = duennaresult ;
	    do_pivots = FALSE ;
	    break ; }
	  default:
	  { errmsg(7,rtnnme,__LINE__,"La Duenna return code",
		   (int) duennaresult) ;
	    do_pivots = FALSE ;
	    lpretval = dyrFATAL ;
	    break ; } } }
/*
  End of the pivoting loop.  Why are we here? The simplest case is that we
  just need to select a new leaving variable --- head back to the top of the
  loop from here.
*/
  if (lpretval == dyrINV) continue ;
/*
  Do we think we're optimal or (dual) unbounded? (I.e., we're reporting
  dyrOPTIMAL, dyrPUNT, or dyrUNBOUND.) If so, the first thing to do is call
  preoptimality to refactor and do accuracy and feasibility checks. If the
  feasibility status agrees with what we're expecting, we can return. Optimal
  should be primal and dual feasible, while unbounded will be dual feasible
  only.

  If preoptimality reports primal and dual feasibility for dyrUNBOUND or
  dyrPUNT, well, unexpected primal feasibility is always a pleasure.

  dyrLOSTPFEAS with dyrPUNT says that there are variables on the pivot
  rejection list, but dy_dealWithPunt has already concluded nothing can be
  done to make them useable. Return dyrPUNT.

  Anything else means that the feasibility status has taken an unexpected
  turn. The working hypothesis is that we need to be more careful in
  selecting pivots for both factoring and pivoting, which we do by tightening
  the current value and lower bound for the pivot selection parameters.  If
  we've lost dual feasibility, bail out (eventually we'll revert to primal
  phase I).  If we've only lost primal feasibility, head back to the top of
  the loop to resume pivoting.
*/
    if (lpretval == dyrOPTIMAL ||
	lpretval == dyrUNBOUND || lpretval == dyrPUNT)
    { optcnt++ ;
#     ifndef DYLP_NDEBUG
      uxpfeas = FALSE ;
      uxnpfeas = FALSE ;
      uxndfeas = FALSE ;
      tmpretval = lpretval ;
#     endif
      preopresult = preoptimality(lpretval,&checks) ;
      switch (preopresult)
      { case dyrOK:
	case dyrPATCHED:
	{ 
#         ifndef DYLP_NDEBUG
	  if (lpretval == dyrUNBOUND || lpretval == dyrPUNT) uxpfeas = TRUE ;
#         endif
	  lpretval = dyrOPTIMAL ;
	  break ; }
	case dyrLOSTDFEAS:
 	{
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.dual >= 1)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
			"\n  (%s)%d: lost dual feasibility.",
		        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
	  uxndfeas = TRUE ;
	  if (flgon(checks,ladPRIMFEAS))
	  { if (lpretval == dyrOPTIMAL) uxnpfeas = TRUE ; }
	  else
	  { if (lpretval != dyrOPTIMAL) uxpfeas = TRUE ; }
#	  endif
	  lpretval = dyrLOSTDFEAS ;
	  break ; }
	case dyrLOSTPFEAS:
	{ if (lpretval == dyrUNBOUND || lpretval == dyrPUNT)
	  { /* nothing to be done */ }
	  else
	  if (lpretval == dyrOPTIMAL)
	  { lpretval = dyrINV ;
	    (void) dy_setpivparms(+1,+1) ;
#	    ifndef DYLP_NDEBUG
	    uxnpfeas = TRUE ;
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
      if (dy_opts->print.dual >= 2)
      { if (uxndfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n\tunexpected loss of dual feasibility at iteration (%s)%d.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxpfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tunexpected primal feasibility at iteration (%s)%d.",
		      dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
        if (uxnpfeas == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n\tunexpected loss of primal feasibility at iteration (%s)%d.",
	     dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; } }
#     endif

      if (lpretval == dyrINV)
      { if (optcnt > 15)
	{ errmsg(387,rtnnme,dy_sys->nme,
		 dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,optcnt) ;
	  lpretval = dyrFATAL ; }
#       ifndef DYLP_NDEBUG
	else
	{ if (dy_opts->print.dual >= 2)
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
  if (dy_lp->degen > 0) (void) dy_dualdegenout(0) ;

# ifndef DYLP_NDEBUG
  if (lpretval == dyrUNBOUND && dy_opts->print.dual >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: system %s is (dual) unbounded.",
	        rtnnme,dy_sys->nme) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL)
  { dy_stats->d2.iters += dy_lp->d2.iters ;
    dy_stats->d2.pivs += dy_lp->d2.pivs ; }
# endif

  return (lpretval) ; }



lpret_enum dy_dual (void)

/*
  This is the driver routine for the dual simplex algorithm. It expects that
  the problem will arrive with an initial dual feasible basis and solution
  (primal and dual variables). Since there's no provision for a phase I in
  this dual simplex implementation, dy_dual is simpler than its counterpart
  dy_primal.

  Parameters: none

  Returns: lpret_enum code
*/

{ lpret_enum retval ;
  dyret_enum dyret ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "dy_dual" ;
# endif

  retval = lpINV ;
  dy_lp->phase = dyDUAL ;
  dy_lp->lpret = lpINV ;
  (void) dy_setpivparms(-100,-100) ;
  (void) dy_setpivparms(+1,+1) ;
  dy_lp->basis.pivs = 0 ;

/*
  Call dual2 to take us to optimality. We'll also make the usual conversion
  of dual unbounded => primal infeasible at this point. (Remember that we don't
  have to worry about dual infeasible => primal unbounded.)
*/
  dy_lp->z = dy_calcobj() ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.dual >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: entering dual phase II, z = %g",
	        rtnnme,dy_lp->z) ;
    if (dy_opts->print.dual >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  ", dual active yb = %g",dy_calcdualobj()) ; }
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
# endif

  dy_lp->d2.iters = 0 ;
  dyret = dual2() ;
  if (dyret == dyrUNBOUND) dyret = dyrINFEAS ;
/*
  What do we have? Optimality is best, but infeasibility or a punt are also
  possible.  Anything else is regarded as a fatal error. (But dylp will try
  primal I in the event of stalling or loss of dual feasibility.)
*/
  retval = dyret2lpret(dyret) ;

# ifndef DYLP_NDEBUG
  if (retval == lpOPTIMAL || retval == lpINFEAS || retval == lpSWING)
  { if (dy_opts->print.dual >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    (%s)%d: dual phase II ended, %d pivots, ",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dy_lp->d2.pivs) ;
      if (retval == lpOPTIMAL)
	dyio_outfmt(dy_logchn,dy_gtxecho,"optimal z = %g.",dy_lp->z) ;
      else
      if (retval == lpINFEAS)
	dyio_outfmt(dy_logchn,dy_gtxecho,"infeas = %g.",dy_lp->infeas) ;
      else
	dyio_outfmt(dy_logchn,dy_gtxecho,"swing %e for %s (%d).",
		    dy_lp->ubnd.ratio,
		    consys_nme(dy_sys,'v',dy_lp->ubnd.ndx,FALSE,NULL),
		    dy_lp->ubnd.ndx) ; } }
  else
  if (retval == lpITERLIM)
  { if (dy_opts->print.dual >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  (%s)%d: dual simplex terminated; iteration limit (%d).",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dy_opts->iterlim) ; } }
  else
  { if (dy_opts->print.dual >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n  (%s)%d: dual phase II failed, status %s after %d pivots.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	       dy_prtlpret(retval),dy_lp->d2.pivs) ; } }
/*
  For curiousity's sake, try the traditional certificate of optimality:
  dual objective == yb == cx == primal objective.
*/
  if (dy_opts->print.dual >= 4 && retval == lpOPTIMAL)
  { double dualobj,primalobj ;

    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: comparing dual and primal objectives.",rtnnme) ;
    dualobj = dy_calcdualobj() ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tdual objective yb = %g.",dualobj) ;
    
    primalobj = dy_calcobj() ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n\tprimal objective cx = %g.",primalobj) ;

    if (!withintol(dualobj,primalobj,dy_tols->dchk))
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tWHOOPS! yb - cx = %g - %g = %g > %g.",
		  dualobj,primalobj,dualobj-primalobj,dy_tols->dchk) ; }
    else
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tobjectives match.") ; }
# endif

  dy_lp->lpret = retval ;

  return (retval) ; }

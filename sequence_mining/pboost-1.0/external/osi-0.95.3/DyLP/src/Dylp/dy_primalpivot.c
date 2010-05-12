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
  This file contains the routines which select entering and leaving variables
  for a primal pivot, and the primal pivot routine.
*/

/*
  A few words on degeneracy and the antidegeneracy mechanism. There's a
  pathological problem which I don't see any way to avoid. Suppose
      |x<i> - bnd<i>| > dy_tols->pfeas
  but
      delta<i> = (x<i> - bnd<i>)/abar<ij> < dy_tols->zero.
  We have a `dirty' degeneracy, and there doesn't seem to be any way to handle
  this gracefully.

  The method chosen to deal with the problem follows the outlines of the
  handling of bogus values. The first time primalout scans one of these
  gems, it'll call for a refactor unless one has just been done. If we can't
  get rid of the problem that way, the solution adopted is to set x<i> to
  whatever bound it's supposed to leave at, take a degenerate pivot, and
  call for an immediate refactor. Results will vary.

  Another problem to watch out for is the case where the perturbation used
  for a restricted subproblem is overly large. dy_degenin attempts to scale
  the perturbation so that it's at most about 1/1000 of the bound the
  variable is currently pinned at, but this isn't always possible. A
  pathological combination of small and large values of abar<i,*> in the
  pivot column can still allow a too-large delta, so that we get a false
  indication of a breakout.  To guard against this, the code detects cases
  where there are no intervening pivots between dy_degenin and dy_degenout
  and reduces the perturbation by a factor of 10 on each occurrence. It's not
  necessary to actually count pivots.  Detection occurs because we never
  escape the loop in dy_primalpivot which selects the leaving variable. The
  sequence goes like this:

    1) Selection of the leaving variable detects degeneracy. degenin is called
       to form a restricted subproblem, and degen_cyclecnt is incremented. The
       loop iterates.

    2) Selection of a leaving variable detects an apparent breakout. degenout
       is called to remove the restricted subproblem. The loop iterates.

    3) If we really had a breakout, selection of a leaving variable would now
       result in a nonzero delta. But if the breakout is false, we'll be back
       at step 1). We need to escape the selection loop in order to reset
       degen_cyclecnt, so we can do this at most three times (currently
       hardcoded) before the selection loop gives up and takes a degenerate
       pivot.

  At the extreme end of the scale, it can happen that it's not possible to
  perturb a variable at all because the bounds are too close together (a
  pathological problem, to be sure, but it occurs when a program (such as an
  ILP branch-and-bound code) is setting the bounds). What happens is that
  when dy_degenin attempts to scale the perturbation based on the difference
  between the bounds, the perturbation ends up smaller than the tolerance
  around the bound. Such variables are flagged with a vstatNOPER qualifier
  and are treated much like fixed variables for purposes of pivoting (i.e.,
  preferentially pivot this variable out of the basis, and don't trigger a
  new level of antidegeneracy). Though it might not be immediately obvious,
  this same problem (perturbation too small) can show up as a result of the
  perturbation reduction described in the previous paragraph.
*/

/*
  As a `lite' alternative to installing a restricted subproblem, bonsaiG
  can look at hyperplane alignment for the hyperplane made tight by the
  leaving variable. For x<i'> leaving, if it's the slack for a<i> it makes
  a<i> tight. If it's some architectural, it makes a bound hyperplane tight.
  The procedure is pretty obvious except for the AlignEdge and AlignObj
  strategies. Consider AlignObj.
  
  Suppose we have incumbent leaving variable x<i'>, which makes hyperplane
  a<i> tight, and candidate x<k'>, which will make hyperplane a<k> tight, and
  delta<i> = delta<k> (most likely both are 0). Let h<k> be dot(-c,a<k>),
  h<i> be dot(-c,a<i>). If h<k> > 0, then a<k> will tend to shut down motion
  in the direction -c, since it'll become tight at this vertex. (In the
  extreme where -c and a<k> go in precisely the same direction, it's the
  proverbial brick wall.) If h<k> < 0, then a<k> opens in the direction of -c
  (alternatively, -c moves away from a<k>).  Summarised, if only one of h<i>
  and h<k> are strictly greater than 0, you have to go with it, because it
  shuts down the other one. Otherwise, choose on relative magnitude.

  For AlignEdge, it's all the same except that h<*> = dot(eta<j>,a<*>). If
  we have h<i> == h<k>, we can try again to resolve the tie using alignment
  with the objective.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_primalpivot.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_primalpivot.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  Define this symbol to enable a thorough check of the updates to cbar and
  gamma, but be aware that it'll almost certainly trigger fatal errors if
  the LP is numerically ill-conditioned.

  #define CHECK_PSE_UPDATES
*/

#if defined(DYLP_STATISTICS) || !defined(DYLP_NDEBUG)

/*
  Pivot counting structure for the antidegeneracy mechanism. This structure
  is sufficient for simple debugging; more complicated stats are collected
  when DYLP_STATISTICS is defined. DYSTATS_MAXDEGEN is defined in dylp.h.

  Field		Definition
  -----		----------
  iterin	iterin[i] is the value of tot.pivs when degeneracy level i was
		entered.
*/

typedef struct { int iterin[DYSTATS_MAXDEGEN] ; } degenstats_struct ;
static degenstats_struct degenstats ;

#endif		/* DYLP_STATISTICS || !DYLP_NDEBUG */

/*
  Anti-cycling variable for the anti-degeneracy mechanism. degen_cyclecnt
  tracks the number of cycles and controls the reduction in the perturbation.
*/

static int degen_cyclecnt ;




dyret_enum dy_degenout (int level)

/*
  This routine backs out all restricted subproblems to the level given by the
  parameter. It copies the original value of x<k> from dy_x into dy_xbasic,
  properly sets the status, adjusts dy_degenset, and resets dy_lp->degen.

  All variables involved in a restricted subproblem were at bound when they
  were collected into the subproblem. Numeric inaccuracy can cause drift over
  the course of pivoting.  Likewise, for variables not part of the degenerate
  set, dy_x should equal dy_xbasic unless there's drift due to inaccuracy.
  This is tested as part of the accuracy checks and the antidegeneracy
  mechanism will be backed out if a problem is detected.  Here, the same
  tests are just paranoia. But, this is the rationale behind resetting
  dy_xbasic to dy_x for variables not in the restricted subproblem, when we
  back out to level 0, and for coping with the possibility of status other
  than BLB and BUB for variables in the restricted subproblem.

  Parameter:
    level:	The target level for removal of restricted subproblems.
  
  Returns: dyrOK if the restoration goes without problem, dyrREQCHK if
	   there's been too much numerical drift since we began the
	   degenerate subproblem.
*/

{ int xkpos,xkndx ;
  double val,*vub,*vlb ;
  flags statk,qualk ;
  dyret_enum retval ;

# ifdef PARANOIA
  const char *rtnnme = "dy_degenout" ;
# endif

# ifdef DYLP_STATISTICS
  int curlvl,curpivs ;
# endif

# ifndef DYLP_NDEBUG
  if (dy_opts->print.degen >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
       "\n    (%s)%d: antidegeneracy dropping to level %d after %d pivots.",
       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,level,
       dy_lp->tot.pivs-degenstats.iterin[dy_lp->degen]) ; }
# endif

# ifdef DYLP_STATISTICS
/*
  Record the iteration counts. This needs to be a loop because we can peel
  off multiple levels.
*/
  if (dy_stats != NULL)
  for (curlvl = dy_lp->degen ; curlvl > level ; curlvl--)
  { if (curlvl < DYSTATS_MAXDEGEN)
    { curpivs = dy_lp->tot.pivs-degenstats.iterin[curlvl] ;
      dy_stats->pdegen[curlvl].totpivs += curpivs ;
      dy_stats->pdegen[curlvl].avgpivs =
	((float) dy_stats->pdegen[curlvl].totpivs)/
					dy_stats->pdegen[curlvl].cnt ;
      if (curpivs > dy_stats->pdegen[curlvl].maxpivs)
	dy_stats->pdegen[curlvl].maxpivs = curpivs ; } }
# endif

  retval = dyrOK ;
  vub = dy_sys->vub ;
  vlb = dy_sys->vlb ;
/*
  Back out restricted subproblems to the level specified by level. Keep in
  mind that all variables initially sucked into a restricted subproblem were
  basic at bound, but the variables currently occupying those slots could be
  free variables.
*/
  for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
  { xkndx = dy_basis[xkpos] ;
    val = dy_x[xkndx] ;
    if (dy_degenset[xkpos] > level)
    { dy_degenset[xkpos] = level ;
      if (level == 0) clrflg(dy_status[xkndx],vstatNOPER) ;
      statk = dy_status[xkndx] ;
      if (flgon(statk,vstatBFR))
      { dy_xbasic[xkpos] = dy_x[xkndx] ; }
      else
      if (flgon(statk,vstatBFX))
      { dy_xbasic[xkpos] = vub[xkndx] ; }
      else
      { qualk = getflg(dy_status[xkndx],vstatQUALS) ;
	if (atbnd(val,vub[xkndx]))
	{ statk = vstatBUB ;
	  val = vub[xkndx] ;
	  dy_x[xkndx] = val ; }
	else
	if (atbnd(val,vlb[xkndx]))
	{ statk = vstatBLB ;
	  val = vlb[xkndx] ;
	  dy_x[xkndx] = val ; }
	else
	{ retval = dyrREQCHK ;
	  if (val < vlb[xkndx])
	    statk = vstatBLLB ;
	  else
	  if (val > vub[xkndx])
	    statk = vstatBUUB ;
	  else
	    statk = vstatB ;
#         ifdef PARANOIA
	  if (withintol(val,vub[xkndx],
			  1000*dy_tols->pfeas*(1+fabs(vub[xkndx]))))
	  { warn(342,rtnnme,dy_sys->nme,
		 consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,val,
		 vlb[xkndx],vub[xkndx],"bogosity",fabs(val-vub[xkndx]),
		 1000*dy_tols->pfeas*(1+fabs(vub[xkndx]))) ; }
	  else
	  if (withintol(val,vlb[xkndx],
			  1000*dy_tols->zero*(1+fabs(vlb[xkndx]))))
	  { warn(342,rtnnme,dy_sys->nme,
		 consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,val,
		 vlb[xkndx],vub[xkndx],"bogosity",fabs(val-vlb[xkndx]),
		 1000*dy_tols->pfeas*(1+fabs(vlb[xkndx]))) ; }
	  else
	  { if (fabs(val-vlb[xkndx]) < fabs(val-vub[xkndx]))
	    { errmsg(342,rtnnme,dy_sys->nme,
		     consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,val,
		     vlb[xkndx],vub[xkndx],"violation",fabs(val-vlb[xkndx]),
		     dy_tols->pfeas*(1+fabs(vlb[xkndx]))) ; }
	    else
	    { errmsg(342,rtnnme,dy_sys->nme,
		     consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,val,
		     vlb[xkndx],vub[xkndx],"violation",fabs(val-vub[xkndx]),
		     dy_tols->pfeas*(1+fabs(vub[xkndx]))) ; } }
#         endif
	}
	setflg(statk,qualk) ;
	dy_status[xkndx] = statk ;
	dy_xbasic[xkpos] = val ; }
#     ifndef DYLP_NDEBUG
      if ((dy_opts->print.degen >= 4 &&
		    flgoff(dy_status[xkndx],vstatBLB|vstatBFX|vstatBUB)) ||
	  dy_opts->print.degen >= 5)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t%s (%d) restored to %g, status %s",
		    consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		    val,dy_prtvstat(dy_status[xkndx])) ;
	if (flgoff(dy_status[xkndx],vstatBLB|vstatBFX|vstatBUB|vstatBFR))
	{ if (fabs(val-vlb[xkndx]) < fabs(val-vub[xkndx]))
	    dyio_outfmt(dy_logchn,dy_gtxecho,", accum. error %g (tol. %g)",
		        fabs(val-vlb[xkndx]),
		        dy_tols->zero*(1.0+fabs(vlb[xkndx]))) ;
	  else
	    dyio_outfmt(dy_logchn,dy_gtxecho,", accum. error %g (tol. %g)",
		        fabs(val-vub[xkndx]),
		        dy_tols->zero*(1.0+fabs(vub[xkndx]))) ; } }
#     endif
    }
    else
    if (level == 0)
    { if (!atbnd(val,dy_xbasic[xkpos]))
      { retval = dyrREQCHK ;
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.degen >= 4)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\t%s (%d) unperturbed, accum. error %g (tol. %g)",
		      consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		      fabs(val-dy_xbasic[xkpos]),
		      dy_tols->zero*(1.0+fabs(dy_xbasic[xkpos]))) ; }
#	endif
      }
      dy_xbasic[xkpos] = val ; } }
    
  dy_lp->degen = level ;
  
  return (retval) ; }



static void dy_degenin (void)

/*
  This routine forms a new restricted subproblem, increasing the degeneracy
  level kept in dy_lp->degen. An initial base perturbation is calculated so
  that the maximum possible perturbation fraction, perturb = (base)*(concnt),
  is less than 1.0e-3. The base is then increased, if necessary, so that the
  minimum perturbation exceeds the basic feasibility tolerance. For each
  variable, the actual perturbation is calculated as
      perturb = base*xkpos*min((1+bnd),.5*(vub-vlb)).
  If we don't scale the perturbation by the size of the bound, we're in
  trouble at the upper end of the scale because the perturbation won't be big
  enough to move out of the tolerance range around the bound. Adding 1 to the
  relevant bound gets us out of trouble around zero. Taking into account the
  range vub-vlb saves us when both bounds are small, or when for some reason
  the range has been seriously reduced. It can happen that the scaled
  perturbation is just too small --- the routine uses the bogus number
  tolerance --- in which case the variable is flagged with a vstatNOPER
  qualifier, and will be treated essentially the same as a fixed variable for
  pivoting purposes (i.e., preferentially move it out of the basis).

  The routine should not be called if degeneracy isn't present, and will do a
  paranoid check to make sure that's true.

  Parameters: none

  Returns: undefined
*/

{ int xkpos,xkndx,oldlvl ;
  double base,perturb,xk,ubk,lbk,toobig ;
  flags xkstatus ;


# if defined(PARANOID) || defined(DYLP_STATISTICS) || !defined(DYLP_NDEBUG)
  int degencnt ;

  degencnt = 0 ;
# endif

/*
  We want the base perturbation to be such that (concnt)*(base) <= 1.0e-3.
  But ... if we're in here again because the previous perturbation was too
  large (degen_cyclecnt > 0), decrease it by a factor of 10 for each repeat.
  Balance that against the notion that if the perturbation is too small
  (less than pfeas), we simply can't see it.
*/
  base = pow(10.0,(-3-ceil(log10(dy_sys->concnt))-degen_cyclecnt)) ;
  while (base <= dy_tols->pfeas) base *= 10 ;
  oldlvl = dy_lp->degen++ ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.degen >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    (%s)%d: antidegeneracy increasing to level %d",
	        dy_prtlpphase(dy_lp->phase,TRUE),
		dy_lp->tot.iters,dy_lp->degen) ;
    if (degen_cyclecnt > 0)
      dyio_outfmt(dy_logchn,dy_gtxecho,", cycle %d",degen_cyclecnt) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,", base perturbation %g%s",
	        base,(dy_opts->print.degen >= 4)?":":".") ; }
# endif
# if defined(DYLP_STATISTICS) || !defined(DYLP_NDEBUG)
  if (dy_lp->degen < DYSTATS_MAXDEGEN)
    degenstats.iterin[dy_lp->degen] = dy_lp->tot.pivs ;
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL && dy_lp->degen < DYSTATS_MAXDEGEN)
  { if (dy_stats->pdegen[0].cnt < dy_lp->degen)
      dy_stats->pdegen[0].cnt = dy_lp->degen ;
    dy_stats->pdegen[dy_lp->degen].cnt++ ; }
# endif

  for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
  { if (dy_degenset[xkpos] != oldlvl) continue ;
    xkndx = dy_basis[xkpos] ;
    xkstatus = dy_status[xkndx] ;
    if (flgoff(xkstatus,vstatBLB|vstatBFX|vstatBUB)) continue ;
    ubk = dy_sys->vub[xkndx] ;
    lbk = dy_sys->vlb[xkndx] ;
    xk = dy_xbasic[xkpos] ;
    toobig = .001*(ubk-lbk) ;
    dy_degenset[xkpos] = dy_lp->degen ;
    switch (xkstatus)
    { case vstatBLB:
      { dy_brkout[xkpos] = 1 ;
	perturb = base*xkpos*(1+fabs(lbk)) ;
	xk += perturb ;
	if (perturb < toobig && !atbnd(xk,lbk))
	{ dy_xbasic[xkpos] = xk ;
	  dy_status[xkndx] = vstatB ; }
	else
	  setflg(dy_status[xkndx],vstatNOPER) ;
	break ; }
      case vstatBUB:
      { dy_brkout[xkpos] = -1 ;
	perturb = base*xkpos*(1+fabs(ubk)) ;
	xk -= perturb ;
	if (perturb < toobig && !atbnd(xk,ubk))
	{ dy_xbasic[xkpos] = xk ;
	  dy_status[xkndx] = vstatB ; }
	else
	  setflg(dy_status[xkndx],vstatNOPER) ;
	break ; }
      case vstatBFX:
      { dy_brkout[xkpos] = 0 ;
	setflg(dy_status[xkndx],vstatNOPER) ;
	break ; } }
#   if defined(PARANOID) || defined(DYLP_STATISTICS)
    degencnt++ ;
#   endif
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.degen >= 5 ||
	(dy_opts->print.degen >= 4 &&
	 flgon(dy_status[xkndx],vstatNOPER) &&
	 flgoff(dy_status[xkndx],vstatBFX)))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s %s (%d) in pos'n %d ",
		  dy_prtvstat(dy_status[xkndx]),
		  consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,xkpos) ;
      if (flgon(dy_status[xkndx],vstatNOPER))
	dyio_outfmt(dy_logchn,dy_gtxecho,"unperturbed.") ;
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,"perturbed from %g (%s) to %g",
		    dy_x[xkndx],dy_prtvstat(xkstatus),dy_xbasic[xkpos]) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,", breakout %s.",
		    (dy_brkout[xkpos] == 1)?"up":"down") ; } }
#   endif
  }

# ifdef PARANOID
  if (degencnt <= 0)
  { errmsg(327,rtnnme,dy_sys->nme) ;
    return ; }
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.degen >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"%s%d variables.",
	        (dy_opts->print.degen <= 4)?", ":"\n\ttotal ",degencnt) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL && dy_lp->degen < DYSTATS_MAXDEGEN)
  { if (dy_stats->pdegen[dy_lp->degen].maxsiz < degencnt)
      dy_stats->pdegen[dy_lp->degen].maxsiz = degencnt ;
    xkndx = dy_stats->pdegen[dy_lp->degen].cnt-1 ;
    perturb = dy_stats->pdegen[dy_lp->degen].avgsiz ;
    dy_stats->pdegen[dy_lp->degen].avgsiz =
			    (float) ((perturb*xkndx+degencnt)/(xkndx+1)) ; }
# endif

  return ; }



static bool pricexk (int k,
		     int *p_j, double *p_ncbarj, bool *p_pivreject)

/*
  This routine encapsulates the code used to decide if the variable x<k>
  prices out more favourably than the current candidate for entry x<j>. It's
  pulled out as a subroutine for readability and to make sure that the same
  rules are used in each location where variables are priced.

  The first part of the routine deals with various a priori reasons for
  disqualification:
    * x<j> is SB and x<k> is not
    * the sign of cbar<k> is wrong for minimisation, given stat<k>
    * x<j> is NBLB or NBUB and cbar<k> == 0
    * x<k> is flagged as NOPIVOT
  SB and NBFR variables get past the cbar<k> == 0 test because we want them
  to be basic in the final answer, so that we have a valid basic solution.

  Given that x<k> passes the above tests and is thus a feasible candidate,
  we're interested in whether it's the best candidate. The tests are
    1) stat<k> == SB and stat<j> != SB
    2) |ncbar<k>| > |ncbar<j>|
    3) stat<k> == NBFR and stat<j> != NBFR

  Parameters:
    k:		index of x<k>, the variable to be priced
    p_j:	(o) if x<k> should supplant x<j>, xjndx will be updated,
		otherwise it'll be unchanged
    p_ncbarj:	(i) |cbar<j>|/||abar<j>||, where x<j> is the current candidate
		(o) if x<k> should supplant x<j>, ncbarj will be updated,
		otherwise it'll be unchanged
    p_pivreject: (o) will be set to TRUE if x<k> could have been a contender
		except for being flagged with the NOPIVOT status qualifier;
		otherwise, it'll be unchanged.
  
  Returns: TRUE if x<k> replaced x<j>, FALSE otherwise.
*/

{ int j ;
  double ncbarj,cbark,ncbark ;
  flags statk,statj ;

  j = *p_j ;
  ncbarj = *p_ncbarj ;
  if (j != 0)
    statj = dy_status[j] ;
  else
    statj = 0 ;

  statk = dy_status[k] ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pricing >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tpricing %s (%d), status %s; ",
	        consys_nme(dy_sys,'v',k,FALSE,NULL),k,
	        dy_prtvstat(statk)) ; }
# endif

/*
  A quick status check. If x<j> is SB and x<k> is not, x<k> loses and we
  can skip the effort of actually pricing it.
*/
  if (flgoff(statk,vstatSB) && flgon(statj,vstatSB))
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << reject (vstatSB) >>") ;
#   endif
    return (FALSE) ; }
/*
  Price x<k> as ncbar<k> = dot(c,h<k>)/||h<k>||, where dot(c,h<k>) is held in
  cbar<k> and ||h<k>||^2 is held in gamma<k>. Check that cbar<k> has the
  right sign and that ncbar<k> is big enough to consider. For SB and NBFR
  variables, we're potentially interested even if ncbar<k> is 0.
*/
  cbark = dy_cbar[k] ;
  if ((cbark < 0 && flgon(statk,vstatNBUB)) ||
      (cbark > 0 && flgon(statk,vstatNBLB)))
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << reject (incompatible status) >>") ;
#   endif
    return (FALSE) ; }
  if (flgoff(statk,vstatSB|vstatNBFR))
  { if (withintol(cbark,0,dy_tols->dfeas))
    {
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pricing >= 3)
	dyio_outfmt(dy_logchn,dy_gtxecho," << reject (zero) >>") ;
#     endif
      return (FALSE) ; } }
  ncbark = fabs(cbark)/sqrt(dy_gamma[k]) ;
  setcleanzero(ncbark,dy_tols->cost) ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pricing >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "cbar<k> = %g, ||h<k>|| = %g, |cbar<k>|/||h<k>|| = %g.",
	        cbark,sqrt(dy_gamma[k]),ncbark) ; }
# endif
/*
  x<k> could enter. Reject if it's flagged with the NOPIVOT qualifier
  and note that we've done so.
*/
  if (flgon(statk,vstatNOPIVOT))
  { *p_pivreject = TRUE ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << reject (vstatNOPIVOT) >>") ;
#   endif
    return (FALSE) ; }
/*
  x<k> is suitable. The only question remaining is whether it's more suitable
  than x<j>. The criteria, in order, are:
    1) stat<k> == vstatSB
    2) ncbar<k> > ncbar<j>
    3) stat<k> == vstatNBFR and stat<j> != vstatNBFR
*/
  if (flgon(statk,vstatSB) && flgoff(statj,vstatSB))
  { j = k ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << accept (vstatSB) >>") ;
#   endif
  }
  else
  if (ncbark > ncbarj)
  { j = k ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << accept (ncbar) >>") ;
#   endif
  }
  else
  if (withintol(ncbark,ncbarj,dy_tols->dfeas) &&
      flgon(statk,vstatNBFR) && flgoff(statj,vstatNBFR))
  { j = k ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << accept (vstatNBFR) >>") ;
#   endif
  }
  else
  { 
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pricing >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho," << reject (inferior) >>") ;
#   endif
  }

  if (j != *p_j)
  { *p_j = j ;
    *p_ncbarj = ncbark ;
    return (TRUE) ; }
  else
  { return (FALSE) ; } }



dyret_enum dy_primalin (int startcol, int scan, int *xjndx, int *nextcol)

/*
  In the normal course of events with PSE pricing, selection of the next
  incoming variable is made during the update of cbar and gamma. But we
  need a backup to select a new candidate when the preselected candidate is
  rejected, and we need a way to select an initial pivot. Hence this routine.
 
  dy_primalin prices columns to come up with an entering candidate.  It uses
  a sort of partial pricing, scanning a block of columns of size scan and
  returning the column index of the best candidate.  The routine is
  persistent, in the sense that it will scan until a candidate is discovered
  or all columns have been scanned. A return value of dyrOPTIMAL means that
  all columns were scanned without finding a variable suitable for pivoting
  into the basis. The last column scanned is returned in lastcol.  (Why all
  this fuss? Historical inertia --- a previous version of the code used
  full-blown partial pricing.)

  If |cbar<j>| < bogus*dfeas, x<j> is not considered for entry unless we've
  just refactored. The assumed sequence of events is that primalin returns
  dyrOPTIMAL with xjndx == 0, triggering preoptimality, which refactors, finds
  loss of dual feasibility, and we're back here ready to consider tiny
  reduced costs.

  If we see cbar<j> of the correct sign, but x<j> is flagged as NOPIVOT, we
  return dyrPUNT. Roughly the same sequence will occur (preoptimality &
  loss of dual feasibility) but with the important side-effect of relaxing
  the pivot selection tolerance (so that some of the x<j> flagged as
  NOPIVOT may become useable).

  Parameters:
    startcol:	The first column to be priced by this call.
    scan:	The number of columns to be priced.
    xjndx:	(o) Index of the candidate x<j>.
    nextcol:	(o) The next column to be priced.

  Returns: dyrOK if a candidate x<j> is found without problem,
	   dyrOPTIMAL if no candidate is selected,
	   dyrPUNT if no candidate was selected but there are potential
		   x<j> on the pivot reject list,
	   dyrFATAL otherwise
*/

{ int xkndx,scanned,total_scanned,scan_blk,this_blk ;
  flags xkstatus ;
  bool pivreject ;
  double ncbarj ;
  dyret_enum retval ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "dy_primalin" ;
# endif

# ifdef PARANOIA
  if (dy_cbar == NULL)
  { errmsg(101,rtnnme,dy_sys->nme,"dy_cbar") ;
    return (dyrFATAL) ; }
  if (startcol <= 0 || startcol > dy_sys->varcnt)
  { errmsg(102,rtnnme,dy_sys->nme,"column",startcol,1,dy_sys->varcnt) ;
    return (dyrFATAL) ; }
# endif
/*
  Set up for the scan. We'll be looking for the largest normalised reduced cost
  |cbar<j>|/||abar<j>|| over nonbasic columns k, so 0 is bad.
*/
  *xjndx = 0 ;
  ncbarj = -dy_tols->inf ;
  xkndx = startcol ;
  scan_blk = minn(scan,dy_sys->varcnt) ;
  retval = dyrINV ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pricing >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n%s: pricing %d columns from %d for %d candidate.",
	        rtnnme,scan_blk,xkndx,1) ; }
# endif
/*
  Open up a pair of loops to start pricing columns. The outer loop gives
  persistence -- scan until we find something or have scanned all the
  columns.  The inner loop steps through the columns.  Make sure the size of
  the scan_blk is not more than the number of columns remaining to be
  scanned, nor larger than the distance to the end of the array.
*/
  pivreject = FALSE ;
  for (total_scanned = 0 ;
       (total_scanned < scan_blk) ||
	 (total_scanned < dy_sys->varcnt && *xjndx == 0) ;
       total_scanned += scanned)
  { this_blk = minn(scan_blk,dy_sys->varcnt-total_scanned) ;
    this_blk = minn(this_blk,dy_sys->varcnt-xkndx+1) ;
    for (scanned = 0 ; scanned < this_blk ; scanned++,xkndx++)
    { xkstatus = dy_status[xkndx] ;
#     ifdef PARANOIA
      if (dy_chkstatus(xkndx) == FALSE) return (dyrFATAL) ;
#     endif
/*
  We're not interested in basic variables or nonbasic fixed variables, so
  step right over them. (Artificials will fall into the NBFX class). We also
  skip variables on the rejected pivot list, flagged as not eligible for
  entry with the NOPIVOT qualifier.
*/
      if (flgon(xkstatus,vstatBASIC|vstatNBFX))
      { 
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.pricing >= 3)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tpricing %s (%d), status %s; << status >>",
		      consys_nme(dy_sys,'v',xkndx,TRUE,NULL),xkndx,
		      dy_prtvstat(xkstatus)) ; }
#       endif
	continue ; }
/*
  Price x<k> and replace the current candidate x<j> if appropriate.
*/
      (void) pricexk(xkndx,xjndx,&ncbarj,&pivreject) ; }
/*
  End of loop for this_blk. We need to wrap xkndx here, after falling out of
  the this_blk loop, to make sure nextcol returns with the proper value when
  the scan ends precisely on the last column.
*/
    if (xkndx > dy_sys->varcnt) xkndx = 1 ; }
/*
  If we're here, then the scan went ok. As a paranoid check, if we didn't
  find any candidates to enter we should have scanned all the columns, hence
  total_scanned should equal dy_sys->varcnt and xkndx should equal startcol.
*/
# ifdef PARANOIA
  if (*xjndx == 0)
  { if (total_scanned != dy_sys->varcnt || xkndx != startcol)
    { errmsg(1,rtnnme,__LINE__) ;
      return (dyrFATAL) ; } }
# endif
/*
  What's the proper return value? If we have a candidate, return dyrOK.
  
  There are three possible reasons for finding no candidate (*xjndx == 0):
   * We have potential pivots on the reject list:
     pivreject == TRUE. We want to return dyrPUNT; see comments at head of
     routine.
   * We're optimal (phase II) or infeasible (phase I):
     pivreject == FALSE. dyrOPTIMAL is the proper return value.
   * We saw some cbar<k> with the correct sign, but they were bogus numbers:
     pivreject == FALSE: dyrOPTIMAL is still the correct return code; see
     comments at head of routine.
*/
  if (*xjndx == 0)
  { if (pivreject == TRUE)
      retval = dyrPUNT ;
    else
      retval = dyrOPTIMAL ; }
  else
  { retval = dyrOK ; }
  
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pricing >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n%s: (%s)%d: scanned %d columns %d to %d, selected %d",
	        rtnnme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	        total_scanned,startcol,(xkndx-1 < 1)?dy_sys->varcnt:(xkndx-1),
	        (*xjndx == 0)?0:1) ;
    if (dy_opts->print.pricing >= 2 && *xjndx != 0)
    { dyio_outchr(dy_logchn,dy_gtxecho,':') ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) scaled reduced cost %g.",
		  consys_nme(dy_sys,'v',*xjndx,TRUE,NULL),*xjndx,ncbarj) ; }
    else
    if (retval == dyrPUNT)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  ",\n\tall suitable x<j> on rejected pivot list.") ; }
    else
    { dyio_outchr(dy_logchn,dy_gtxecho,'.') ; } }
# endif

/*
  That's it. Report the next column to be scanned and return.
*/
  *nextcol = xkndx ;

  return (retval) ; }



static double cdothyper (int xkndx, int dirk)

/*
  This routine is used in the context of primal antidegen lite, where we're
  trying to choose the leaving variable using hyperplane alignment with
  the objective. The choice of leaving variable determines the hyperplane
  which will become tight on completion of the pivot.

  The routine calculates the dot product of -c with a constraint normal a<k>,
  normalised by ||a<k>||. ||c|| doesn't change over the various candidates,
  so we don't need to include it in the normalisation for comparison purposes.

  The only thing to watch is that we get the right direction for the
  hyperplane. See pdirdothyper for extensive comments that explain how
  dirk is applied.

  Parameters:
    xkndx:	index of x<k>, the candidate to leave
    dirk:	the direction of change; -1 to decrease, +1 to increase

  Returns: -dirk*dot(-c,a<k>)/||a<k>||, as described above, or NaN if something
	   goes awry.
*/

{ double dotprod ;

/*
  If x<k> is a slack, we want the coefficients of the explicit constraint.
*/
  if (xkndx <= dy_sys->concnt)
  { dotprod = consys_dotrow(dy_sys,xkndx,dy_sys->obj) ;
    dotprod = dirk*dotprod/consys_2normrow(dy_sys,xkndx) ; }
/*
  If x<k> is an architectural, we have a bound constraint, coefficient 1.0
  for x<k> only.
*/
  else
  { dotprod = dirk*(-dy_sys->obj[xkndx]) ; }

  setcleanzero(dotprod,dy_tols->zero) ;
  
  return (dotprod) ; }



static double pdirdothyper (int xjndx, double *abarj, int dirj,
			    int xkndx, int dirk)

/*
  This routine is used in the context of primal antidegen lite, where we're
  trying to choose the leaving variable x<i> using hyperplane alignment.  The
  choice of entering variable x<j> determined the edge eta<j> to be traversed.
  The choice of leaving variable will determine the hyperplane a<i>x <= b<i>
  which will become tight on completion of the pivot.
  
  This routine is called with a candidate leaving variable x<k>, and we're
  interested in the alignment of eta<j> with the normal a<k>, which will be
  dot(eta<j>,a<k>)/(||eta<j>||)(||a<k>||). Since ||eta<j>|| won't change
  in the course of selecting a leaving variable, we don't bother with it
  here. If x<j> is decreasing, we need to multiply eta<j> by -1 to get the
  direction of motion right.

  When x<k> is a slack, the constraint coming tight is (nominally) a<k>x <=
  b<k>.  If x<k> is in fact an upper bounded slack increasing to its upper
  bound, we need to multiply a<k> by -1, since the constraint is really a<k>x
  >= blow<i>. We also need to be careful in phase I that we compensate when
  we're approaching a constraint from the wrong side. Here's a table:

      Direction		bound	    constraint		    correction
      decreasing	  lb	    ax <= b	 		+1
      decreasing	  ub	    ax >= blow, wrong side	+1
      increasing	  ub	    ax >= blow			-1
      increasing	  lb	    ax <= b, wrong side		-1

  The `wrong side' cases occur in phase I. The bottom line is multiplication
  by -dir<k> will do it. Here too, the dot product boils down to -abar<kj>.

  When x<k> is an architectural variable, the constraint coming tight is
  implicit --- either the upper bound constraint x<k> <= ub<k> or the lower
  bound constraint x<k> >= l<k>. It's easiest to do this as a separate
  case. The constraint normal is e<k>, and the norm is 1.0. The dot product
  simply selects -abar<kj> from eta<j>.  A little care is required to make
  sure we always put the constraint into <= form, particularly when we
  approach a bound from the wrong side in phase I, but it boils down to
  multiplication by dir<k> to get it right.

  The paranoid calculation actually does the dot product for explicit
  constraints as a check, retrieving the coefficients and then translating
  between basis pos'n and variable index, since abar<j> is indexed by basis
  pos'n while a<k> uses variable indices.

  It's a small but ugly truth that the dynamic nature of the constraint system
  allows for us to encounter `empty' constraints, with no active variables
  except for the slack. Hence we need to add 1 to ||a<k>|| to protect against
  division by zero.

  Parameters:
    xjndx:	index of x<j>, the entering variable
    abarj:	inv(B)a<j>, part of the desired direction of motion for
		this pivot
    dirj:	direction of motion for x<j>, -1 if it's decreasing from its
		upper bound, +1 if it's increasing from its lower bound
    xkndx:	index of x<k>, the candidate to leave
    dirk:	proposed direction of motion for x<k>, -1 to decrease,
		+1 to increase

  Returns: dir<j>*(-dir<k>)*dot(eta<j>,a<k>)/||a<k>||, or NaN if something
	   goes awry.
*/

{ double dotprod,normak,abarkj ;

# ifdef PARANOIA
  int pkndx,xqndx,xqpos ;
  double xqcoeff ;
  pkvec_struct *ak ;
  pkcoeff_struct *akq ;

  const char *rtnnme = "pdirdothyper" ;
# endif

  abarkj = abarj[dy_var2basis[xkndx]] ;
/*
  Start with the case where x<k> is a slack, and we're dealing with an
  explicit constraint.
*/
  if (xkndx <= dy_sys->concnt)
  { 
#   ifdef PARANOIA
/*
  Retrieve the constraint normal a<k> and do the dot product. The procedure
  is to walk the (sparse) row, translating each column coefficient to a basis
  pos'n to find the proper value in eta<j>. Exclude the slack variable (we're
  working with the inequality, not the equality that results from adding the
  slack). Multiply the accumulated result by dir<j> to account for the
  entering variable's direction of motion, and -dir<k> to account for the
  incoming variable's motion (hence the type of constraint coming tight).
  Check that we're equal to abar<kj>.
*/
    ak = NULL ;
    if (consys_getrow_pk(dy_sys,xkndx,&ak) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,
	     "row",consys_nme(dy_sys,'c',xkndx,TRUE,NULL),xkndx) ;
      if (ak != NULL) pkvec_free(ak) ;
      return (quiet_nan(0)) ; }
    dotprod = 0 ;
    for (pkndx = 0, akq = ak->coeffs ; pkndx < ak->cnt ; pkndx++, akq++)
    { xqndx = akq->ndx ;
      if (xqndx == xkndx) continue ;
      xqcoeff = akq->val ;
      xqpos = dy_var2basis[xqndx] ;
      if (xqpos > 0)
      { dotprod -= xqcoeff*abarj[xqpos] ; }
      else
      if (xqndx == xjndx)
      { dotprod += xqcoeff ; } }
    pkvec_free(ak) ;
    if (!withintol(dotprod,abarkj,dy_tols->zero))
    { errmsg(401,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters+1,dotprod,xjndx,xkndx,xkndx,xjndx,-abarkj,
	       fabs(dotprod-abarkj),dy_tols->zero) ;
      return (quiet_nan(0)) ; }
#   endif

    normak = consys_2normrow(dy_sys,xkndx) ;
    dotprod = dirj*(-dirk)*(-abarkj)/(normak+1) ; }
/*
  Architecturals need a bound constraint. As per the analysis at the head of
  the routine, multiplication by dir<k> is necessary to properly convert to
  a <= constraint under all circumstances.
*/
  else
  { dotprod = dirj*dirk*(-abarkj) ; }
  
  setcleanzero(dotprod,dy_tols->zero) ;
  return (dotprod) ; }




static dyret_enum primalout (int xjndx, int indir,
			     double *abarj, double maxabarj,
			     int *xindx, int *outdir, double *deltaj)

/*
  This routine selects the leaving variable x<i> given the entering variable
  x<j> and the direction of movement (indir) as x<j> enters the basis. The
  calculation is the standard limit test, starting from the expression that
  x<B> = inv(B)b - inv(B)a<j>x<j>. We're interesting in finding out the limit
  on the change in x<j> (delta<j>) imposed by lb<j>, ub<j>, and by the upper
  and lower bounds on the basic variables x<k>.  The minimum delta<j> over all
  these determines the leaving variable x<i>.  x<j> can move from one bound
  to the other and be both the entering and leaving variable.

  This is complicated slightly if we're working a perturbed subproblem due to
  degeneracy. In this case, we only consider the basic variables associated
  with the constraints involved in the degeneracy.

  If two variables have equal delta<j> (we're about to pivot into a
  degenerate vertex) the tie is broken based on the stability of the pivot
  (default), or according to one of the alignment-based antidegen lite
  strategies (which see). If we find an x<k> with a delta of 0 (we're already
  at a degenerate vertex) the default is to take it (and abort the scan)
  unless the pivot looks bogus. Again, one of the antidegen lite strategies
  can be specified as an option.

  Where the tie is between a basic variable x<k> vs. a bound-to-bound swing
  of x<j>, x<j> always wins, because we avoid pivoting the basis.

  Note that during phase I, if we're out-of-bound, we calculate the limit on
  delta<j> based on the maximum move (across the near bound to the far bound,
  when it exists, otherwise just to the near bound). The wisdom of this is
  debatable. It offers the possibility of many variables gaining feasibilty
  in a single pivot, but there are real costs associated with maintaining
  PSE pricing information, and hypothetical costs associated with whether this
  is actually good in terms of the phase I objective (c<j> for a variable
  should go to 0 when it hits feasibility, so moving to the far bound might
  not look favorable if we were to recalculate the reduced cost) and in
  terms of the phase II objective (difficult to quantify).

  Parameters:
    xjndx:	Index of the entering variable x<j>.
    indir:	Direction of motion of x<j>.
		    -1: decreasing from upper bound
		     1: increasing from lower bound
    abarj:	Ftran'd column inv(B)a<j> associated with x<j>.
    maxabarj:	MAX{i}(abar<ij>)
    xindx:	(o) Index of the leaving variable x<i>. Also valid for
		    return code dyrLOSTPFEAS (in which case it is the index
		    of the variable where feasibility loss was discovered)
		    and dyrREQCHK (in which case it is the index of the
		    variable whose pivot a<ij> was declared bogus).
    outdir:	(o) Direction of motion of x<i>, coded as:
		    -1: decreasing to lower bound
		     1: increasing to upper bound
    deltaj:	(o) Absolute value of the allowable change in x<j>.

  
  Returns: dyret_enum code, as follows:
    dyrOK:	a strictly basic leaving variable was successfully selected
		(this includes dirty degeneracy)
    dyrDEGEN:	a basic at bound leaving variable is selected; the pivot will
		be (cleanly) degenerate
    dyrMADPIV:	the pivot coefficient abar<ij> would be numerically unstable
    dyrREQCHK:	a possibly bogus abar<ij> was selected as the pivot, and
		refactoring seems wise before trying to use it
		(basis.etas > 1 is the criteria)
    dyrUNBOUND:	the problem is unbounded
    dyrLOSTPFEAS: primal feasibility has been lost
    dyrFATAL:	fatal confusion
*/

{ int xkpos,xkndx,outk,degencnt ;
  flags xkstatus ;
  double abarij,ratioij,aligni,deltak,abarkj,ratiokj,bndk,alignk ;
  bool newxi ;
  dyret_enum retval ;
  const char *rtnnme = "primalout" ;

  retval = dyrINV ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: selecting leaving variable, iteration %d",
	        rtnnme,dy_lp->tot.iters+1) ;
    if (dy_opts->print.pivoting >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tentering variable %s (%d) %s",
		  consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,
		  (indir > 0)?"increasing":"decreasing") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
	    "\n\tPos'n\tVariable\tValue\tabar<k,j>\tbound\tdelta\tDisp") ; } }
# endif
/*
  Start off by assuming the entering variable x<j> will also be the leaving
  variable. Calculate the limits on delta<j> imposed by moving to ub<j> or
  lb<j>. (Using dy_x allows uniform handling of normal nonbasic and
  superbasic variables.) If delta<j> is infinite, the relevant bound doesn't
  exist, and we'll need to find our leaving variable among the basic
  variables.

  If the allowable delta is less than or equal to 0, we have serious
  confusion. If delta < 0, the variable is nonbasic and outside its bounds.
  If delta == 0, the variable is trying to reenter at the same bound it left
  with, and it shouldn't have been chosen by primalin as a candidate for
  entry.
*/
  *xindx = xjndx ;
  *outdir = indir ;
  abarij = 1.0 ;
  if (indir == -1)
  { *deltaj = dy_x[xjndx]-dy_sys->vlb[xjndx] ; }
  else
  { *deltaj = dy_sys->vub[xjndx]-dy_x[xjndx] ; }
# ifdef PARANOIA
  if (*deltaj < -dy_tols->zero)
  { errmsg(325,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	   dy_lp->tot.iters+1,consys_nme(dy_sys,'v',xjndx,FALSE,NULL),
	   xjndx,dy_x[xjndx],dy_sys->vlb[xjndx],dy_sys->vub[xjndx],
	   dy_prtvstat(dy_status[xjndx])) ;
    return (dyrFATAL) ; }
  if (withintol(*deltaj,0,dy_tols->zero))
  { errmsg(326,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	   dy_lp->tot.iters+1,consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,
	   dy_x[xjndx],(indir == 1)?"ub":"lb",dy_prtvstat(dy_status[xjndx]),
	   (indir == 1)?"ub":"lb") ;
    return (dyrFATAL) ; }
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 3 && *deltaj < dy_tols->inf)
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\t  n/a\t%-8s (%d)\t%7g\t%9g\t%7g\t%7g\tleaving at %s",
	        consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,
	        (indir > 0)?dy_sys->vlb[xjndx]:dy_sys->vub[xjndx],1.0,
	        (indir > 0)?dy_sys->vub[xjndx]:dy_sys->vlb[xjndx],*deltaj,
	        (*outdir > 0)?"ub":"lb") ;
# endif
/*
  Open a loop to step through the basic variables. We'll keep at it until
  we've scanned them all, barring something extraordinary.
  
  For each variable x<k>, we
  * Check if it's eligible in the presence of degeneracy (i.e., the variable
    is associated with a constraint in the current degenerate set).
  * Check that abar<kj> is nonzero.
  * Check for basic free -- these never leave.
  * Flip the sign of abar<kj> if x<j> is entering and decreasing. This gives
    two cases in the next step, instead of four.
*/
  newxi = FALSE ;
  ratioij = quiet_nan(0) ;
  aligni = quiet_nan(0) ;
  degencnt = 0 ;
  alignk = 0 ;
  bndk = quiet_nan(0) ;
  for (xkpos = 1 ; xkpos <= dy_sys->concnt && retval == dyrINV ; xkpos++)
  { if (dy_lp->degen > 0 && dy_degenset[xkpos] != dy_lp->degen) continue ;
    abarkj = abarj[xkpos] ;
    if (withintol(abarkj,0,dy_tols->zero)) continue ;
    xkndx = dy_basis[xkpos] ;
    xkstatus = dy_status[xkndx] ;
    if (flgon(xkstatus,vstatBFR)) continue ;
    if (indir < 0) abarkj = -abarkj ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivoting >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%5d\t%-8s (%d)\t%g\t%g\t",xkpos,
		  consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		  dy_xbasic[xkpos],(indir < 0)?-abarkj:abarkj) ;
#   endif
/*
  Figure out delta<k> = |(bnd<k> - x<k>)/abar<kj>|. The hard part is deciding
  on the proper limiting bound. Remember that if abar<kj> > 0, we're moving
  x<k> toward its lower bound, and if abar<kj> < 0, we're moving x<k> toward
  its upper bound. If we're in phase I dealing with an infeasible variable,
  look for the far bound first, and fall back to the near bound if necessary.
  If we're degenerate (BLB and x<k> decreasing, BUB and x<k> increasing, or
  BFX), declare delta<k> to be 0, period.
*/
    if (abarkj > 0.0)
    { if (flgon(xkstatus,vstatBLLB)) continue ;
      outk = -1 ;
      if (flgon(xkstatus,vstatBFX|vstatBLB))
      { deltak = 0.0 ;
#       ifndef DYLP_NDEBUG
	bndk = dy_sys->vlb[xkndx] ;
#       endif
      }
      else
      { bndk = dy_sys->vlb[xkndx] ;
	if (bndk <= -dy_tols->inf)
	{ if (dy_lp->phase == dyPRIMAL1 && flgon(xkstatus,vstatBUUB))
	    bndk = dy_sys->vub[xkndx] ;
	  else
	    continue ; }
	deltak = dy_xbasic[xkpos]-bndk ; } }
    else
    { if (flgon(xkstatus,vstatBUUB)) continue ;
      outk = 1 ;
      if (flgon(xkstatus,vstatBFX|vstatBUB))
      { deltak = 0.0 ;
#       ifndef DYLP_NDEBUG
	bndk = dy_sys->vlb[xkndx] ;
#       endif
      }
      else
      { bndk = dy_sys->vub[xkndx] ;
	if (bndk >= dy_tols->inf)
	{ if (dy_lp->phase == dyPRIMAL1 && flgon(xkstatus,vstatBLLB))
	    bndk = dy_sys->vlb[xkndx] ;
	  else
	    continue ; }
	deltak = bndk-dy_xbasic[xkpos] ; } }
/*
  To make it to here, delta<k> should be finite and positive. If delta<k> <
  -dy_tols->pfeas*(1+fabs(bndk)), we've lost feasibility and it isn't
  reflected in the status. In phase I, it just means we've lost accuracy. In
  phase II, it's more serious --- we'll have to revert to phase I. In any
  event, we're done with searching. For values between loss of feasibility
  and 0, force delta<k> to 0.
*/
    setcleanzero(deltak,dy_tols->zero) ;
    if (deltak < 0.0)
    { if (deltak < -dy_tols->pfeas*(1+fabs(bndk)))
      {
#       ifndef DYLP_NDEBUG
	if (dy_opts->print.phase2 >= 1 || dy_opts->print.phase1 >= 1)
	  warn(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters+1,consys_nme(dy_sys,'v',xkndx,FALSE,NULL),
	       xkndx,dy_prtvstat(xkstatus),dy_sys->vlb[xkndx],
	       dy_xbasic[xkpos],dy_sys->vub[xkndx],-deltak,
	       dy_tols->pfeas*(1+fabs(bndk))) ;
#       endif
	retval = dyrLOSTPFEAS ;
	*xindx = xkndx ;
	*outdir = outk ;
	*deltaj = deltak ;
	continue ; }
      deltak = 0.0 ; }
/*
  See how much we can really move. Even with a reasonable delta<k>, we can
  still end up with something that looks like degeneracy if abar<kj> is very
  large.
*/
    deltak = fabs(deltak/abarkj) ;
    setcleanzero(deltak,dy_tols->zero) ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivoting >= 3)
      dyio_outfmt(dy_logchn,dy_gtxecho,"%g\t%g\t",bndk,deltak) ;
#   endif
/*
  We have delta<k>. Now, is x<k> a better candidate to leave than the current
  incumbent x<i>? If delta<k> is really smaller, there's no contest.
*/
    ratiokj = quiet_nan(0) ;
    if (deltak < *deltaj)
    { newxi = TRUE ;
      ratiokj = dy_chkpiv(abarkj,maxabarj) ;
      degencnt = 0 ; }
/*
  If it's a tie, the first preference is to stick with x<j> (and avoid pivoting
  the basis altogether). Only if abar<kj> looks good do we go to further
  tie-breakers.
  
  The next choice is to pivot out x<k> if it has status BFX and the pivot is
  tolerable. Once out, it'll never reenter and can be purged from the active
  problem. Still, go for the best pivot value between BFX variables.

  After that, we resort to whatever tiebreaking strategy is in effect. The
  default action is Pivot: scan all basic variables, and break ties with
  |abar<kj>|. See the top of the file for additional comments about AlignObj
  and AlignEdge.

  We do not want a toleranced comparison here --- small differences in
  delta<j> multiplied by large abar<kj> can result in loss of feasibility.
*/
    else
    if (deltak == *deltaj && *xindx != xjndx)
    { ratiokj = dy_chkpiv(abarkj,maxabarj) ;
      if (ratiokj >= 1.0)
      { if (flgon(xkstatus,vstatBFX))
	{ if (!(flgon(dy_status[*xindx],vstatBFX) &&
		fabs(abarij) >= fabs(abarkj)))
	    newxi = TRUE ; }
	else
	{ switch (dy_opts->degenlite)
	  { case 0: /* pivotabort */
	    { if (ratiokj > ratioij) newxi = TRUE ;
	      break ; }
	    case 1: /* pivot */
	    { if (ratiokj > ratioij) newxi = TRUE ;
	      degencnt++ ;
	      break ; }
	    case 2: /* alignobj */
	    case 3: /* alignedge */
	    { if (dy_opts->degenlite == 2)
	      { if (degencnt == 0) aligni = cdothyper(*xindx,*outdir) ;
		alignk = cdothyper(xkndx,outk) ; }
	      else
	      { if (degencnt == 0)
		  aligni =  pdirdothyper(xjndx,abarj,indir,*xindx,*outdir) ;
		alignk =  pdirdothyper(xjndx,abarj,indir,xkndx,outk) ; }
	      degencnt++ ;
	      if (aligni > 0 && alignk <= 0)
	      { /* keep x<i> */ }
	      else
	      if (aligni <= 0 && alignk > 0)
	      { newxi = TRUE ; }
	      else
	      if (fabs(aligni) > fabs(alignk))
	      { newxi = TRUE ; }
	      else
	      if (aligni == alignk)
	      { if (ratiokj > ratioij) newxi = TRUE ; }
	      break ; }
	    case 4: /* perpobj */
	    case 5: /* perpedge */
	    { if (dy_opts->degenlite == 4)
	      { if (degencnt == 0) aligni = cdothyper(*xindx,*outdir) ;
		alignk = cdothyper(xkndx,outk) ; }
	      else
	      { if (degencnt == 0)
		  aligni =  pdirdothyper(xjndx,abarj,indir,*xindx,*outdir) ;
		alignk =  pdirdothyper(xjndx,abarj,indir,xkndx,outk) ; }
	      degencnt++ ;
	      if (alignk > aligni) newxi = TRUE ;
	      break ; }
	  } } } }
/*
  Is x<k> better? If so, make x<k> the new leaving variable x<i>. If the
  user's choice of antidegen lite option is pivotabort, maybe we can skip
  the rest of the scan.
*/
    if (newxi == TRUE)
    { *deltaj = deltak ;
      *xindx = xkndx ;
      *outdir = outk ;
      abarij = indir*abarkj ;
      ratioij = ratiokj ;
      aligni = alignk ;
      newxi = FALSE ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivoting >= 3)
      { if (*xindx == xkndx)
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\tleaving at %s",
		      (bndk == dy_sys->vub[xkndx])?"ub":"lb") ; } 
#     endif
      if (dy_opts->degenlite == 0 && deltak == 0 && ratioij >= 1.0) break ; } }
/*
  Why are we here? We could have broken out of the loop due to detection of
  loss of feasibility (dyrLOSTPFEAS), we could have broken out because we've
  found a BFX variable to pivot out, or we could have finished the scan (in
  the latter two cases, retval should still be dyrINV).  If delta<j> is still
  infinite, we're unbounded. A finite and nonzero delta<j> and a numerically
  stable abar<ij> is what we hope for (dyrOK). We distinguish between clean
  degeneracy (x<i> at bound) vs.  `dirty' degeneracy (x<i> not at bound, but
  delta<j> = 0). The first gets dyrDEGEN, the second dyrOK (which will not
  trigger the antidegeneracy mechanism). If we're about to pivot on a bogus
  number, ask for a refactor first (dyrREQCHK). A numerically unstable pivot
  returns dyrMADPIV. But note that if we're doing a bound-to-bound pivot on a
  nonbasic, we're automatically ok --- neither degeneracy or a numerically
  unstable pivot is possible.

  Note that both the antidegeneracy machinery and the fact that we're running
  with a partial constraint system can lead to apparent unboundedness, even
  in primal phase I.
*/
  switch (retval)
  { case dyrINV:
    { if (*deltaj < dy_tols->inf)
      { if (*xindx != xjndx)
	{ if (ratioij >= 1.0)
	  { if (dy_lp->basis.etas > 1 &&
		withintol(abarij,0,dy_tols->bogus*dy_tols->zero))
	    { retval = dyrREQCHK ;
#  	      ifndef DYLP_NDEBUG
	      warn(381,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters+1,"abar",*xindx,xjndx,abarij,
		   dy_tols->bogus*dy_tols->zero,
		   dy_tols->bogus*dy_tols->zero-fabs(abarij)) ;
#  	      endif
	    }
	    else
	    if (*deltaj == 0 && degencnt > 0 &&
		flgon(dy_status[*xindx],vstatBUB|vstatBLB|vstatBFX))
	    { retval = dyrDEGEN ; }
	    else
	    { retval = dyrOK ; } }
	  else
	  { retval = dyrMADPIV ; } }
	else
	{ retval = dyrOK ; } }
      else
      { *xindx = -1 ;
#       ifndef DYLP_NDEBUG
	if (dy_lp->degen == 0 &&
	    (dy_opts->print.phase1 >= 2 || dy_opts->print.phase2 >= 2))
	  warn(324,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters+1) ;
#       endif
	retval = dyrUNBOUND ; }
      break ; }
    case dyrLOSTPFEAS:
    { break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      return (dyrFATAL) ; } }
/*
  We're done, except perhaps for printing some information.
*/
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting == 1)
    dyio_outfmt(dy_logchn,dy_gtxecho,"...") ;
  if ((retval == dyrOK || retval == dyrDEGEN) && dy_opts->print.pivoting >= 1)
  { if (xjndx != *xindx)
    { xkpos = dy_var2basis[*xindx] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    selected %s (%d) = %g to leave pos'n %d at",
		  consys_nme(dy_sys,'v',*xindx,FALSE,NULL),*xindx,
		  dy_xbasic[xkpos],xkpos) ;
      if (*outdir > 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho," %s = %g, ",
		    (dy_status[*xindx] != vstatBLLB)?"ub":"lb",
		    (dy_status[*xindx] != vstatBLLB)?
			    dy_sys->vub[*xindx]:dy_sys->vlb[*xindx]) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho," %s = %g, ",
		    (dy_status[*xindx] != vstatBUUB)?"lb":"ub",
		    (dy_status[*xindx] != vstatBUUB)?
			    dy_sys->vlb[*xindx]:dy_sys->vub[*xindx]) ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "abar<%d,%d> = %g, ",xjndx,*xindx,abarij) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    selected %s (%d) = %g to change to %s = %g, ",
		  consys_nme(dy_sys,'v',*xindx,FALSE,NULL),*xindx,dy_x[*xindx],
		  (*outdir > 0)?"ub":"lb",
		  (*outdir > 0)?dy_sys->vub[*xindx]:dy_sys->vlb[*xindx]) ; }
    if (retval == dyrOK)
      dyio_outfmt(dy_logchn,dy_gtxecho,"delta = %g.",*deltaj) ;
    else
      dyio_outfmt(dy_logchn,dy_gtxecho,"degenerate.") ; }
  if (retval == dyrDEGEN &&
      (dy_opts->print.phase1 >= 3 || dy_opts->print.phase2 >= 3))
  { xkpos = dy_var2basis[*xindx] ;
    xkstatus = dy_status[*xindx] ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n   (%s)%d %s pos'n %d, %s %s (%d) = %g = %s = %g",
		dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		dy_prtdyret(retval),xkpos,dy_prtvstat(xkstatus),
		consys_nme(dy_sys,'v',*xindx,FALSE,NULL),*xindx,
		dy_xbasic[xkpos],(*outdir > 0)?"ub":"lb",
		(*outdir > 0)?dy_sys->vub[*xindx]:dy_sys->vlb[*xindx]) ;
    if (dy_opts->degenlite >= 2 && dy_opts->degenlite <= 5)
    { if (degencnt > 0)
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", align = %g, deg = %d.",aligni,degencnt) ; }
    else
    { dyio_outchr(dy_logchn,dy_gtxecho,'.') ; } }
# endif

  return (retval) ; }



static dyret_enum primalupdate (int xjndx, int indir,
				int xindx, int outdir,
				double *abarj, double delta, double *betai)

/*
  This routine is responsible for updating the various data structures which
  hold the basis, status, and primal and dual variable values.

  Nondegenerate pivots will have delta > 0, and require a full scan of the
  basic variables to update their values. In the process, the value of the
  leaving variable will be updated.
  
  `Dirty' degeneracy gives us delta == 0, but requires some special handling
  to force the leaving variable to bound.
  
  Other pivots with delta == 0 require no updates to the basic variables.

  There's a fair bit of attention given here to controlling roundoff, with a
  (hopefully) intelligent choice of tolerances so that we get clean values
  when x approaches 0 or a bound. Generally, the rule is to try and scale the
  snap interval by the larger of the absolute value of the target value or
  the distance travelled to get there.

  If we decide we've seen a bogus value, we can't just bail out. That leaves
  the status information completely inconsistent and we'll fail a subsequent
  status check. For the purposes of this routine, a value is bogus if either
  |value| < tols.bogus or |lb-value| < tols.bogus or |ub-value| < tols.bogus.
  (The latter two come from the notion of trying to snap variables to bound.)

  Parameters:
    xjndx:	index of the entering variable x<j>
    indir:	the direction of change of the entering variable
    xindx:	index of the leaving variable x<i>
    outdir:	the direction of change of the outgoing variable
    abarj:	the ftran'd column abar<j> = inv(B)a<j>
    delta:	the amount of change of the entering variable; always >= 0.
    betai:	the ith row of inv(B); used for updating dual variables.

  NOTE that cbar<j>, abar<j>, and beta<i> are all calculated prior to the
  pivot!

  Returns: dyrOK if all goes well,
	   dyrSWING if a variable grew too much,
	   dyrREQCHK if a bogus value is calculated or the pivot is
		     dirty degenerate, and
	   dyrFATAL if paranoid checks fail.
*/

{ int xkpos,xkndx,xipos ;
  flags xkstatus,quals ;
  dyret_enum retval ;
  double val,deltak,ubk,lbk,eps0,epsu,epsl,cbarj,abarij ;
  bool dirtyz,swing ;
  double swingratio,maxswing ;
  int swingndx ;
  const char *rtnnme = "primalupdate" ;

# ifndef DYLP_NDEBUG
  int print ;
# endif

  retval = dyrOK ;
  dirtyz = FALSE ;
  swing = FALSE ;
  maxswing = 0 ;
  swingndx = -1 ;
  xipos = dy_var2basis[xindx] ;
  abarij = abarj[xipos] ;
  cbarj = dy_cbar[xjndx] ;

# ifndef DYLP_NDEBUG
  if (dy_lp->phase == dyPRIMAL1)
    print = dy_opts->print.phase1 ;
  else
    print = dy_opts->print.phase2 ;

  if (print >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: updating at iteration %d:",rtnnme,
	        dy_lp->tot.iters+1) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n\t%s (%d) entering pos'n %d from %s%g, delta %g, cbarj %g.",
	        consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,xipos,
	        (dy_status[xjndx] == vstatSB)?"":((indir == 1)?"lb ":"ub "),
	        dy_x[xjndx],(indir == 1)?delta:-delta,dy_cbar[xjndx]) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) = %g leaving at ",
	        consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		dy_xbasic[xipos]) ;
    if (outdir == 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%s %g, pivot %g.",
		  flgon(dy_status[xindx],vstatBLLB)?"lb":"ub",
		  flgon(dy_status[xindx],vstatBLLB)?
		    dy_sys->vlb[xindx]:dy_sys->vub[xindx],abarij) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,"%s %g, pivot %g.",
		  flgon(dy_status[xindx],vstatBUUB)?"ub":"lb",
		  flgon(dy_status[xindx],vstatBUUB)?
		    dy_sys->vub[xindx]:dy_sys->vlb[xindx],abarij) ; } }
# endif


/*
  Update the value and status of the basic variables to reflect the change in
  x<j>. The calculation is straightforward, from the formulas:
    z = c<B>inv(B) + (c<j> - c<B>inv(B)a<j>)*delta = z<old> + cbar<j>*delta
    x<B> = inv(B)b - (inv(B)a<j>)*delta = x<B,old> - abar<j>*delta
  Note that while the antidegeneracy mechanism is active, we're really doing
  degenerate pivots in the original, unperturbed problem, so we shouldn't
  change the objective or any variables not part of the restricted problem.
  Nor do we update dy_x for the variables in the restricted subproblem -- dy_x
  is holding their original values for when we back out the perturbation. As
  we're updating the basic variables, collect the 1-norm to scale the primal
  feasibility tolerance.

  Paranoid checks:
    * We should never attempt to change the value of a fixed variable. But,
      we'll occasionally fudge this in order to promote a sane pivot over
      mad pivots. In this case, the change should never exceed pfeas.
    * In phase I, infeasible variables should not overshoot their opposite
      bound (see comments in primalout), and feasible variables should not
      loose feasibility.
    * In phase II, variables should always remain feasible.
*/
  if (delta > 0.0)
  { if (indir == -1) delta = -delta ;
    if (dy_lp->degen == 0)
    { dy_lp->z += cbarj*delta ;
      setcleanzero(dy_lp->z,dy_tols->zero) ; }
    for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
    { if (abarj[xkpos] != 0 && dy_degenset[xkpos] == dy_lp->degen)
      { deltak = abarj[xkpos]*delta ;
	if (withintol(deltak,0,dy_tols->zero)) continue ;
	xkndx = dy_basis[xkpos] ;
	xkstatus = getflg(dy_status[xkndx],vstatSTATUS) ;
	quals = getflg(dy_status[xkndx],vstatQUALS) ;
	eps0 = dy_tols->zero ;
	ubk = dy_sys->vub[xkndx] ;
	if (ubk < dy_tols->inf)
	{ epsu = dy_tols->pfeas*(1.0+fabs(ubk)) ; }
	else
	  epsu = 0 ;
	lbk = dy_sys->vlb[xkndx] ;
	if (-dy_tols->inf < lbk)
	{ epsl = dy_tols->pfeas*(1.0+fabs(lbk)) ; }
	else
	  epsl = 0 ;
	val = dy_xbasic[xkpos]-deltak ;
	setcleanzero(val,eps0) ;
	if (val != 0.0 &&
	    fabs(val) < eps0*dy_tols->bogus && dy_lp->basis.etas > 1)
	{ retval = dyrREQCHK ;
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.pivoting >= 1)
	    warn(374,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,"x",xkndx,fabs(val),
		 eps0*dy_tols->bogus,eps0*dy_tols->bogus-val) ;
#	  endif
	}
#       ifdef PARANOIA
	if (flgon(xkstatus,vstatBFX) && fabs(deltak) > dy_tols->pfeas)
	{ errmsg(345,rtnnme,dy_sys->nme,
		 consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,val,
		 dy_lp->tot.iters+1,delta) ;
	  return (dyrFATAL) ; }
	if (dy_lp->phase == dyPRIMAL1)
	{ if ((flgon(xkstatus,vstatBLLB) && val > ubk+epsu) ||
	      (flgon(xkstatus,vstatBUUB) && val < lbk-epsl))
	  { errmsg(344,rtnnme,dy_sys->nme,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(xkstatus),dy_lp->tot.iters+1,val,lbk,ubk,delta) ;
		return (dyrFATAL) ; }
	  if (flgon(xkstatus,vstatBLB|vstatB|vstatBUB) &&
	      (val < lbk-dy_tols->bogus*epsl || val > ubk+dy_tols->bogus*epsu))
	  { if (val < lbk-dy_tols->bogus*epsl)
	      errmsg(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters+1,
		     consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		     dy_prtvstat(xkstatus),lbk,val,ubk,
		     (lbk-dy_tols->bogus*epsl)-val,dy_tols->bogus*epsl) ;
	    else
	      errmsg(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters+1,
		     consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		     dy_prtvstat(xkstatus),lbk,val,ubk,
		     val-(ubk+dy_tols->bogus*epsu),dy_tols->bogus*epsu) ;
	    return (dyrFATAL) ; } }
	else
	if (dy_lp->phase == dyPRIMAL2 &&
	    (val < lbk-dy_tols->bogus*epsl || val > ubk+dy_tols->bogus*epsu))
	{ if (val < lbk-epsl)
	    errmsg(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters+1,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(xkstatus),lbk,val,ubk,
		   (lbk-dy_tols->bogus*epsl)-val,dy_tols->bogus*epsl) ;
	  else
	    errmsg(323,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters+1,
		   consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		   dy_prtvstat(xkstatus),lbk,val,ubk,
		   val-(ubk+dy_tols->bogus*epsu),dy_tols->bogus*epsu) ;
	    return (dyrFATAL) ; }
#       endif

	switch (xkstatus)
	{ case vstatB:
	  case vstatBLB:
	  case vstatBUB:
	  { if (atbnd(val,ubk))
	    { dy_status[xkndx] = vstatBUB ; }
	    else
	    if (atbnd(val,lbk))
	    { dy_status[xkndx] = vstatBLB ; }
	    else
	    { dy_status[xkndx] = vstatB ; }
	    break ; }
	  case vstatBLLB:
	  { if (belowbnd(val,lbk))
	    { /* do nothing */ }
	    else
	    if (atbnd(val,lbk))
	    { if (lbk == ubk)
		dy_status[xkndx] = vstatBFX ;
	      else
		dy_status[xkndx] = vstatBLB ; }
	    else
	    if (belowbnd(val,ubk))
	    { dy_status[xkndx] = vstatB ; }
	    else
	    { dy_status[xkndx] = vstatBUB ; }
	    break ; }
	  case vstatBUUB:
	  { if (abovebnd(val,ubk))
	    { /* do nothing */ }
	    else
	    if (atbnd(val,ubk))
	    { if (lbk == ubk)
		dy_status[xkndx] = vstatBFX ;
	      else
		dy_status[xkndx] = vstatBUB ; }
	    else
	    if (abovebnd(val,lbk))
	    { dy_status[xkndx] = vstatB ; }
	    else
	    { dy_status[xkndx] = vstatBLB ; }
	    break ; }
	  case vstatBFR:
	  case vstatBFX:
	  { break ; }
	  default:
	  { errmsg(1,rtnnme,__LINE__) ;
	    return (dyrFATAL) ; } }
	setflg(dy_status[xkndx],quals) ;
/*
  Check for bogus values, within the bogosity tolerance of a bound but not
  close enough to snap to it.
*/
        if (flgoff(dy_status[xkndx],vstatBFX|vstatBLB|vstatBUB))
	{ if (fabs(val-lbk) < epsl*dy_tols->bogus && dy_lp->basis.etas > 1) 
	  { retval = dyrREQCHK ;
#	    ifndef DYLP_NDEBUG
	    if (dy_opts->print.pivoting >= 1)
	      warn(375,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,"x",xkndx,"lb",xkndx,
		   val,lbk,val-lbk,epsl*dy_tols->bogus,
		   epsl*dy_tols->bogus-(val-lbk)) ;
#	    endif
	  }
	  else
	  if (fabs(ubk-val) < epsu*dy_tols->bogus && dy_lp->basis.etas > 1) 
	  { retval = dyrREQCHK ;
#	    ifndef DYLP_NDEBUG
	    if (dy_opts->print.pivoting >= 1)
	      warn(375,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		   dy_lp->tot.iters,"ub",xkndx,"x",xkndx,
		   ubk,val,ubk-val,epsu*dy_tols->bogus,
		   epsu*dy_tols->bogus-(ubk-val)) ;
#	    endif
	  } }
	swingratio = (fabs(val)+1)/(fabs(dy_xbasic[xkpos])+1) ;
	if (swingratio > dy_tols->swing)
	{ swing = TRUE ;
	  if (swingratio > maxswing)
	  { maxswing = swingratio ;
	    swingndx = xkndx ; } }
	dy_xbasic[xkpos] = val ;
	if (dy_lp->degen == 0) dy_x[xkndx] = val ;
#       ifdef PARANOIA
/*
  Check that x<i> has acquired the proper status after the update of basic
  variables.
*/
	if ((xkndx == xindx) &&
	    !flgon(dy_status[xindx],vstatBLB|vstatBFX|vstatBUB))
	{ if (fabs(ubk-val) < fabs(lbk-val))
	  { if (fabs(ubk-val) < 100*epsu)
	    { warn(357,rtnnme,dy_sys->nme,
		   consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		   dy_prtvstat(dy_status[xindx]),"ub",ubk,val,val-ubk,epsu) ; }
	    else
	    { errmsg(357,rtnnme,dy_sys->nme,
		     consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		     dy_prtvstat(dy_status[xindx]),"ub",ubk,val,val-ubk,epsu) ;
	      return (dyrFATAL) ; } }
	  else
	  { if (fabs(lbk-val) < 100*epsl)
	    { warn(357,rtnnme,dy_sys->nme,
		   consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		   dy_prtvstat(dy_status[xindx]),"lb",lbk,val,lbk-val,epsl) ; }
	    else
	    { errmsg(357,rtnnme,dy_sys->nme,
		     consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
		     dy_prtvstat(dy_status[xindx]),"lb",lbk,val,lbk-val,epsl) ;
	      return (dyrFATAL) ; } } }
#       endif

      } } }
/*
  This next clause is for `dirty' degeneracy. See the notes at the top of
  the file. In effect, we clean up by forcing the variable to the appropriate
  bound. A little care is required to decide what the appropriate bound really
  is when in phase I. Set a flag if we need to recalculate the objective
  function (after completing the update bookkeeping for x<i> and x<j>).
*/
  else
  if (xjndx != xindx && flgoff(dy_status[xindx],vstatBFX|vstatBLB|vstatBUB))
  { if (dy_lp->degen == 0) dirtyz = TRUE ;
    xkstatus = getflg(dy_status[xindx],vstatSTATUS) ;
    lbk = dy_sys->vlb[xindx] ;
    ubk = dy_sys->vub[xindx] ;
    if (lbk == ubk)
    { val = lbk ;
      xkstatus = vstatBFX ; }
    else
    if (outdir < 0)
    { val = lbk ;
      if (val <= -dy_tols->inf)
      { if (dy_lp->phase == dyPRIMAL1 && flgon(xkstatus,vstatBUUB))
	{ val = ubk ;
	  xkstatus = vstatBUB ; }
	else
	{ errmsg(382,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,consys_nme(dy_sys,'v',xindx,FALSE,NULL),
		 xindx,"lb",val) ;
	  return (dyrFATAL) ; } }
      else
      { xkstatus = vstatBLB ; } }
    else
    { val = ubk ;
      if (val >= dy_tols->inf)
      { if (dy_lp->phase == dyPRIMAL1 && flgon(xkstatus,vstatBLLB))
	{ val = lbk ;
	  xkstatus = vstatBLB ; }
	else
	{ errmsg(382,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,consys_nme(dy_sys,'v',xindx,FALSE,NULL),
		 xindx,"ub",val) ;
	  return (dyrFATAL) ; } }
      else
      { xkstatus = vstatBUB ; } }

    dy_xbasic[xipos] = val ;
    dy_status[xindx] = 0 ;
    setflg(dy_status[xindx],xkstatus) ;
    if (dy_lp->degen == 0) dy_x[xindx] = val ;

    retval = dyrREQCHK ;

#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivoting >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	   "\n      %s (%d) = %g, %s, leaving at %s, dirty degenerate pivot.",
	   consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,dy_xbasic[xipos],
	   dy_prtvstat(dy_status[xindx]),(outdir < 0)?"lb":"ub") ; }
#   endif
  }

/*
  Deal with the entering and leaving variables.

  In the case where the entering variable x<j> and the leaving variable x<i>
  are different variables, the value at entry is obtained from dy_x, rather
  than going to the upper or lower bound vectors. This lets us handle
  superbasics and variables at bound uniformly. If delta != 0, the final
  status for x<j> should be B or BFR (if x<j> was NBFR). But it's not that
  simple --- a very small or very large delta, combined with a relatively
  large pfeas tolerance, can leave us at either of BLB or BUB, regardless of
  where we started from.

  If the pivot is degenerate, we simply convert from nonbasic to the
  equivalent basic status.  As above, if the antidegeneracy mechanism is
  active, we don't update dy_x for x<j>. On the other hand, if it is active,
  we have to update the breakout entry for this basis position.  The proper
  value is the entry direction, since in the original, nondegenerate problem,
  the entering variable remains at the bound it entered from.

  We updated all basic variables in the previous loop, including the leaving
  variable x<i>, so it will have status BLB, BUB, or BFX as we move it to the
  nonbasic partition here. We'll set dy_x from the bounds vector as a check
  on accumulated inaccuracy due to incremental modification.

  When x<j> == x<i> is both the entering and leaving variable, we need only
  change its status and value. Note that if a superbasic is selected to enter
  and then driven to bound, it will also be selected as the leaving variable,
  so if it remains in the basis, it'll have status B.
*/
  if (xjndx != xindx)
  { dy_var2basis[xjndx] = xipos ;
    dy_var2basis[xindx] = 0 ;
    dy_basis[xipos] = xjndx ;
    xkstatus = getflg(dy_status[xindx],vstatSTATUS) ;
    if (xkstatus == vstatBLB)
    { dy_status[xindx] = vstatNBLB ;
      dy_x[xindx] = dy_sys->vlb[xindx] ; }
    else
    if (xkstatus == vstatBUB)
    { dy_status[xindx] = vstatNBUB ;
      dy_x[xindx] = dy_sys->vub[xindx] ; }
    else
    { dy_status[xindx] = vstatNBFX ;
      dy_x[xindx] = dy_sys->vlb[xindx] ; }
    if (delta != 0)
    { val = dy_x[xjndx]+delta ;
      swingratio = (fabs(val)+1)/(fabs(dy_x[xjndx])+1) ;
      if (swingratio > dy_tols->swing)
      { swing = TRUE ;
	if (swingratio > maxswing)
	{ maxswing = swingratio ;
	  swingndx = xjndx ; } }
      setcleanzero(val,dy_tols->zero*(1.0+fabs(delta))) ;
      switch (dy_status[xjndx])
      { case vstatNBLB:
	case vstatNBUB:
        case vstatSB:
	{ if (atbnd(val,dy_sys->vub[xjndx]))
	    dy_status[xjndx] = vstatBUB ;
	  else
	  if (atbnd(val,dy_sys->vlb[xjndx]))
	    dy_status[xjndx] = vstatBLB ;
	  else
	    dy_status[xjndx] = vstatB ;
	  break ; }
        case vstatNBFR:
	{ dy_status[xjndx] = vstatBFR ;
	  break ; }
#       ifdef PARANOIA
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (dyrFATAL) ; }
#       endif
      } }
    else
    { val = dy_x[xjndx] ;
      switch (dy_status[xjndx])
      { case vstatNBLB:
	{ dy_status[xjndx] = vstatBLB ;
	  break ; }
	case vstatNBUB:
	{ dy_status[xjndx] = vstatBUB ;
	  break ; }
	case vstatSB:
	{ dy_status[xjndx] = vstatB ;
	  break ; }
	case vstatNBFR:
	{ dy_status[xjndx] = vstatBFR ;
	  break ; }
#       ifdef PARANOIA
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (dyrFATAL) ; }
#       endif
      } }
    if (dy_lp->degen > 0)
      dy_brkout[xipos] = indir ;
    else
      dy_x[xjndx] = val ;
    dy_xbasic[xipos] = val ; }
  else
  { if (outdir == 1)
    { dy_status[xindx] = vstatNBUB ;
      dy_x[xindx] = dy_sys->vub[xindx] ; }
    else
    { dy_status[xindx] = vstatNBLB ; 
      dy_x[xindx] = dy_sys->vlb[xindx] ; } }
/*
  Do we need to recalculate the objective, as a result of dirty degeneracy?
*/
  if (dirtyz == TRUE)
  { dy_lp->z = dy_calcobj() ; }
# ifdef PARANOIA
  else
  { val = dy_calcobj() ;
    if (fabs(val-dy_lp->z) > fabs(.001*(1+val)))
    { warn(405,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	   dy_lp->tot.iters+1,dy_lp->z,val,fabs(dy_lp->z-val),
	   fabs(.001*val)) ; } }
# endif
/*
  Do we need to update the duals? y = c<B>inv(B), so an update is required
  only if we actually pivoted (x<j> != x<i>) and thus changed c<B>. We'll
  update even during primal I, but keep in mind that c<B> might be revised in
  tweakp1obj if a variable other than x<i> gained feasibility with this
  pivot.
*/
  if (xjndx != xindx)
  { if (fabs(cbarj) > dy_tols->cost)
    { for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      { deltak = cbarj*betai[xkpos] ;
	deltak = deltak/abarij ;
	val = dy_y[xkpos]+deltak ;
	setcleanzero(val,dy_tols->cost) ;
	if (val != 0.0 &&
	    dy_lp->basis.etas > 1 && fabs(val) < dy_tols->cost*dy_tols->bogus)
	{ retval = dyrREQCHK ;
#         ifndef DYLP_NDEBUG
	  if (dy_opts->print.pivoting >= 1)
	    warn(374,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,"y",xkpos,fabs(val),
		 dy_tols->cost*dy_tols->bogus,
		 dy_tols->cost*dy_tols->bogus-val) ;
#         endif
	}
	dy_y[xkpos] = val ; } } }
/*
  Decide on a return value. Swing overrides the others, as it'll cause us to
  pop out of simplex. (But if there are no loadable constraints, then let's
  not, eh?)
*/
  if (swing == TRUE)
  { if (dy_lp->sys.loadablecons == TRUE)
    { retval = dyrSWING ; }
    dy_lp->ubnd.ndx = swingndx ;
    dy_lp->ubnd.ratio = maxswing ;
#   ifndef DYLP_NDEBUG
    if (print >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    Pseudo-unbounded: growth %e for %s (%d)",
		  dy_lp->ubnd.ratio,
		  consys_nme(dy_sys,'v',dy_lp->ubnd.ndx,FALSE,NULL),
		  dy_lp->ubnd.ndx) ; }
#   endif
  }
/*
  That's it, except for some informational printing.
*/
# ifndef DYLP_NDEBUG
  if (print >= 5)
  { bool first,all ;
  
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\trevised objective %g.",dy_lp->z) ;
#   ifdef PARANOIA
    if (dy_lp->phase == dyPRIMAL2)
    { deltak = dy_calcobj() ;
      if (!atbnd(deltak,dy_lp->z))
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tWHOOPS! updated obj - true obj = %g - %g = %g > %g",
		    dy_lp->z,deltak,dy_lp->z-deltak,dy_tols->dchk) ; }
#   endif
    if (print >= 6)
      all = TRUE ;
    else
      all = FALSE ;
    first = TRUE ;
    for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      if (abarj[xkpos] != 0 || all == TRUE)
      { if (first == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%sprimal variables:",
		      (all == TRUE)?"":"revised ") ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n%8s%20s%16s%16s%16s %s","pos'n","var (ndx)",
		      "lb","val","ub","status") ;
	  first = FALSE ; }
	xkndx = dy_basis[xkpos] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n%8d%14s (%3d)%16.8g%16.8g%16.8g %s",xkpos,
		    consys_nme(dy_sys,'v',xkndx,FALSE,NULL),xkndx,
		    dy_sys->vlb[xkndx],dy_xbasic[xkpos],dy_sys->vub[xkndx],
		    dy_prtvstat(dy_status[xkndx])) ; }
    if (first == TRUE)
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tno change to primal variables.") ;
    if (print >= 7)
    { if (xindx != xjndx)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    dual variables, cbar tolerance %g",
		    dy_tols->dfeas) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n%8s%20s%16s","pos'n","constraint","val") ;
	for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n%8d%20s%16.8g",xkpos,
		      consys_nme(dy_sys,'c',xkpos,FALSE,NULL),dy_y[xkpos]) ; } }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    no change to dual variables.") ; } } }
# endif

  return (retval) ; }



#ifdef CHECK_PSE_UPDATES

static bool check_pse_update (int xkndx, double u_cbark, double u_gammak)

/*
  This routine checks x<k> for consistent status, then does one of the
  following:
    * For nonbasic variables it does a `from scratch' calculation of cbar<k>
      and gamma<k> to check the accuracy of the PSE update calculations.
    * For basic variables, it checks that inv(B)a<k> is a unit vector with
      a 1 in the basis position occupied by x<k>.
  
  The routine should not be called for a variable with status NBFX.

  The calculations are:

    cbar<k> = c<k> - c<B>(inv(B)a<k>)

    gamma~<k> = ||abar~<k>||^2, where abar<k> = inv(B)a<k> and abar~<k>
    is abar<k> with non-reference-frame entries removed.

  Parameters:
    xkndx:	index for column
    u_cbark:	updated cbar<k>
    u_gammak:	updated gamma<k>

  Returns: TRUE if the updated values agree with values calculated from
	   first principles, FALSE otherwise.
*/

{ int xipos,xindx,xkpos ;
  double *abark,cbark,gammak ;
  bool retval ;

  const char *rtnnme = "check_pse_update" ;

/*
  Make sure we're ok with the status of the variable.
*/
  if (dy_chkstatus(xkndx) == FALSE) return (FALSE) ;
/*
  The next thing we want to do is extract the column and FTRAN it.
*/
  abark = NULL ;
  if (consys_getcol_ex(dy_sys,xkndx,&abark) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,
	   "column",consys_nme(dy_sys,'v',xkndx,TRUE,NULL),xkndx) ;
    if (abark != NULL) FREE(abark) ;
    return (FALSE) ; }
  dy_ftran(abark,FALSE) ;
  retval = TRUE ;
/*
  Do the appropriate check. For x<k> basic, check that inv(B)a<k> is a unit
  vector with a 1 in the basis position occupied by x<k>.
*/
  if (flgon(dy_status[xkndx],vstatBASIC))
  { xkpos = dy_var2basis[xkndx] ;
    for (xipos = 1 ; xipos <= dy_sys->concnt ; xipos++)
    { if (xipos == xkpos)
      { if (!withintol(abark[xipos],1.0,dy_tols->zero))
	{ errmsg(385,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,xipos,xkndx,abark[xipos],1.0,
		 abark[xipos]-1.0,dy_tols->zero) ;
	  retval = FALSE ; } }
      else
      { if (!withintol(abark[xipos],0.0,dy_tols->zero))
	{ errmsg(385,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,xipos,xkndx,abark[xipos],0.0,
		 abark[xipos],dy_tols->zero) ;
	  retval = FALSE ; } } } }
/*
  For nonbasic variables, calculate the reduced cost, c<k> - c<B>abar<k>.
  Calculate the projected column norm using only those elements in the
  reference frame. gamma<k> must be at least 1, since every nonbasic variable
  should be in the reference frame.
*/
  else
  { cbark = dy_sys->obj[xkndx] ;
    gammak = 1.0 ;
    for (xipos = 1 ; xipos <= dy_sys->concnt ; xipos++)
    { xindx = dy_basis[xipos] ;
      cbark -= dy_sys->obj[xindx]*abark[xipos] ;
      if (dy_frame[xindx] == TRUE)
	gammak += abark[xipos]*abark[xipos] ; }
    if (!withintol(cbark,u_cbark,dy_tols->reframe*(1+fabs(cbark))))
    { errmsg(388,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"cbar",xkndx,u_cbark,cbark,fabs(u_cbark-cbark),
	     dy_tols->reframe*(1+fabs(cbark))) ;
      retval = FALSE ; }
    if (!withintol(gammak,u_gammak,dy_tols->reframe*(1+fabs(gammak))))
    { errmsg(388,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"gamma",xkndx,u_gammak,gammak,
	     fabs(u_gammak-gammak),dy_tols->reframe*(1+fabs(gammak))) ;
	retval = FALSE ; } }
  
  if (abark != NULL) FREE(abark) ;

  return (retval) ; }

#endif



static dyret_enum pseupdate (int xjndx, int xindx, int *candxj,
			     double *abarj, double *v, double *betai)

/*
  This routine updates the reduced cost vector dy_cbar and the column norm
  vector dy_gamma. We're doing projected steepest edge, and vectors tagged
  with `~' include only entries corresponding to variables in the reference
  frame.
  
  The update formulas are as follows:

    cbar'<i> = -cbar<j>/abar<ij>
    cbar'<k> = cbar<k> - cbar<j>*(abar<ik>/abar<ij>)	k != i

    gamma'<i> = gamma<j>/abar<ij>^2
    gamma'<k> = gamma<k> - 2*(abar<ik>/abar<ij>)*dot(a<k>,v) +
			(abar<ik>/abar<ij>)^2*gamma<j>		k != i
  
  It'd be easy to collapse the update formula for cbar'<k> to that of
  cbar'<i>, since cbar<i> = 0 and abar<ii> = 1 for x<i> basic. But, the
  algebra to collapse the formulas for gamma isn't so transparent, and we can
  skip a dot product, so they are kept separate. Having introduced the special
  case, we might as well use it for cbar'<i> too, since there's some advantage
  in robustness (we don't have to worry about values of gamma<k> or cbar<k>
  for basic variables, since they get reset when the become nonbasic).

  We use betai to calculate the values abar<ik> for each column. If the
  leaving variable x<i> is not already a member of the reference framework,
  it's added, and gamma<i> is boosted by 1.

  pseupdate is also responsible for deciding if a reference frame reset is in
  order, and carrying it out if necessary. The test is that the iteratively
  updated column norm gamma<j> be within a percentage (dy_tols->reframe) of
  the exact norm ||abar~<j>||^2.

  If the LP is numerically illconditioned, the gamma updates can begin to
  drift, and this may not be picked up in a timely manner by the reframe test.
  The one thing we watch out for when calculating any gamma update is a
  value < 1, and reset the update to 1 when this happens. The cbar updates
  don't seem to give any trouble, numerically speaking.

  NOTE that the basis has been pivoted by the time this routine is called, and
  the arrays dy_basis, dy_var2basis, and dy_status reflect this.

  Parameters:
    xjndx:	index of the entering variable
    xindx:	index of the leaving variable
    candxj:	(o) the index of the variable chosen to enter on the next
		pivot
    abarj:	inv(B)a<j>, calculated prior to pivoting the basis
    v:		abar~<j>inv(B), calculated prior to pivoting the basis
    betai:	e<i>inv(B), row i of the basis inverse, calculated prior to
		pivoting the basis

  Returns: dyrOK if the update proceeds without error and a new candidate
		 x<j> is selected;
	   dyrOPTIMAL if the update proceeds without error and no x<j>
		 is selected;
	   dyrPUNT if the update proceeds without error but all candidate
		 x<j> were flagged with the NOPIVOT qualifier;
	   dyrFATAL if there's a problem (only if we're paranoid).
*/

{ int xkpos,xkndx,xipos ;
  double abarij,cbarj,gammaj,abarik,alphak,cbark,gammak,akdotv,candcbarj ;
  flags xkstatus ;
  bool reset,pivreject ;
  dyret_enum retval ;

# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "pseupdate" ;
# endif

/*
  Do a little prep, pulling out common values and setting initial values.
  Remember that x<j> is now occupying the basis position vacated by x<i>.
  Set the reduced cost for x<j> to 0, since it's now basic.
*/
  xipos = dy_var2basis[xjndx] ;
  abarij = abarj[xipos] ;
  cbarj = dy_cbar[xjndx] ;
  dy_cbar[xjndx] = 0 ;
  retval = dyrINV ;
  candcbarj = -dy_tols->inf ;
  *candxj = 0 ;
/*
  Do we need to reset the frame of reference? The test is that the iteratively
  updated norm gamma<j> is within dy_tols->reframe percent of the exact value
  ||abar~<j>||^2+1. We need to be careful here --- we're contemplating norms
  as of prior to the pivot, so we need the previous basis image.  Working on
  the ``make the common case fast'' theory, even if a reset is needed we'll
  proceed with the updates as if nothing has happened, then rewrite dy_frame
  and dy_gamma at the end.
*/
  dy_basis[xipos] = xindx ;
  gammaj = 1 ;
  for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
    if (dy_frame[dy_basis[xkpos]] == TRUE) gammaj += abarj[xkpos]*abarj[xkpos] ;
  if (!withintol(dy_gamma[xjndx],gammaj,dy_tols->reframe*gammaj))
  { reset = TRUE ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivoting >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  %s: (%s)%d: resetting reference frame; trigger %s (%d)",
		  dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\texact gamma<j> = %g, approx = %g, error = %g, tol = %g.",
		  gammaj,dy_gamma[xjndx],fabs(gammaj-dy_gamma[xjndx]),
		  dy_tols->reframe*gammaj) ; }
#   endif
  }
  else
  { reset = FALSE ; }
  dy_basis[xipos] = xjndx ;
  dy_gamma[xjndx] = gammaj ;
/*
  Open a loop to walk the nonbasic variables, updating the reduced costs and
  column norms. (But note that we don't bother with nonbasic fixed variables,
  which will never pivot in.)

  The first thing we do, once we enter the loop, is calculate alpha<k>. If
  it's zero, then no update will be needed.
*/
  pivreject = FALSE ;
  for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
  { xkstatus = dy_status[xkndx] ;
    if (flgon(xkstatus,vstatBASIC|vstatNBFX))
    {
#     ifdef CHECK_PSE_UPDATES
      if (flgon(xkstatus,vstatBASIC))
	if (check_pse_update(xkndx,0,0) == FALSE) return (dyrFATAL) ;
#     endif
      continue ; }
    abarik = consys_dotcol(dy_sys,xkndx,betai) ;
#   ifdef PARANOIA
    if (isnan(abarik) == TRUE)
    { errmsg(320,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"beta<i>",xkndx,"PSE update") ;
      return (dyrFATAL) ; }
#   endif
    alphak = abarik/abarij ;
    setcleanzero(alphak,dy_tols->zero) ;
#   ifdef PARANOIA
/*
  Since x<i> was basic when we extracted betai, abar<i,i> should be 1.0.
  While we don't check for it explicitly, alpha<i> should not be 0 --- it
  could only happen if |1/abar<i,j>| < dy_tols->zero, and you have to wonder
  why we're here in that case. So we're guaranteed to execute the updates on
  cbar<i> and gamma<i> (which don't depend on alpha<i>, in spite of the
  appearance that they could be skipped due to alpha<i> = 0).

  Arguably we should return a fatal error when this check fails, but it can
  be dealt with as an accuracy problem. Scaling can tighten tols.zero, which
  we don't want here, so hardwire the default value of 1.0e-11.
*/
    if (xindx == xkndx)
    { if (!withintol(abarik,1.0,dy_tols->bogus*1.0e-11))
      { if (!withintol(abarik,1.0,.001))
	{ errmsg(385,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,xindx,xkndx,abarik,1.0,abarik-1.0,.001) ; }
	else
	{ warn(385,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	       dy_lp->tot.iters,xindx,xkndx,abarik,1.0,abarik-1.0,
	       dy_tols->bogus*dy_tols->zero) ; } } }
#   endif
/*
  If alpha<k> is nonzero, we have actual update work ahead. For handy
  reference, the formulas are:

    cbar'<i> = -cbar<j>/abar<ij>
    cbar'<k> = cbar<k> - cbar<j>*alpha<k>

    gamma'<i> = gamma<j>/abar<ij>^2
    gamma'<k> = gamma<k> - 2*alpha<k>*dot(a<k>,v) + alpha<k>^2*gamma<j>	
  
  Also, we need to add x<i> to the reference frame, if it's not already in
  it (including bumping gamma<i> by +1). While it may look like
  we're risking the update of gamma<i>, we're not --- as mentioned above,
  alpha<i> should not be 0.

  As a pragmatic solution to numerical error, if gamma<k> comes up less than
  1, set it to 1.
*/
    if (alphak != 0)
    { if (xkndx != xindx)
      { cbark = dy_cbar[xkndx]-cbarj*alphak ;
	akdotv = consys_dotcol(dy_sys,xkndx,v) ;
#       ifdef PARANOIA
	if (isnan(akdotv) == TRUE)
	{ errmsg(320,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		     dy_lp->tot.iters,"v",xkndx,"PSE update") ;
	  return (dyrFATAL) ; }
#       endif
	gammak = dy_gamma[xkndx]-alphak*(2*akdotv-alphak*gammaj) ; }
      else
      { cbark = -cbarj/abarij ;
	gammak = gammaj/(abarij*abarij) ;
	if (dy_frame[xkndx] != TRUE)
	{ dy_frame[xkndx] = TRUE ;
	  gammak += 1 ; } }
      setcleanzero(cbark,dy_tols->cost) ;
      dy_cbar[xkndx] = cbark ;
      if (gammak < 1.0) gammak = 1.0 ;
      dy_gamma[xkndx] = gammak ;
    }
#   ifdef CHECK_PSE_UPDATES
    if (check_pse_update(xkndx,dy_cbar[xkndx],dy_gamma[xkndx]) == FALSE)
      return (dyrFATAL) ;
#   endif
/*
  Updates are finished.  Price x<k> and replace the current candidate x<j> if
  appropriate.
*/
    (void) pricexk(xkndx,candxj,&candcbarj,&pivreject) ; }
/*
  That's the end of the PSE update & pricing loop.  Did we decide up at the
  top that we need a reference frame reset? If so, do it now, before we get
  into setting the return value.
*/
  if (reset == TRUE)
  { memset(dy_frame,0,(dy_sys->varcnt+1)*sizeof(bool)) ;
    memset(dy_gamma,0,(dy_sys->varcnt+1)*sizeof(double)) ;
    for (xkndx = 1 ; xkndx <= dy_sys->varcnt ; xkndx++)
    { if (flgon(dy_status[xkndx],vstatNONBASIC|vstatEXOTIC))
      { dy_frame[xkndx] = TRUE ;
	dy_gamma[xkndx] = 1.0 ; } } }
/*
  What's the proper return value? If we've found a candidate x<j>, return
  dyrOK.
  
  There are three possible reasons for finding no candidate (candxj == 0):
   * We have potential pivots on the reject list:
     pivreject == TRUE. We want to return dyrPUNT; see comments at head of
     dy_primalin.
   * We're optimal (phase II) or infeasible (phase I):
     pivreject == FALSE. dyrOPTIMAL is the proper return value.
   * We saw some cbar<k> with the correct sign, but they were bogus numbers:
     pivreject == FALSE. dyrOPTIMAL is still the correct return code; see
     comments at the head of dy_primalin.
*/
  if (*candxj == 0)
  { if (pivreject == TRUE)
      retval = dyrPUNT ;
    else
      retval = dyrOPTIMAL ; }
  else
  { retval = dyrOK ; }
 
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pricing >= 2)
  { if (*candxj != 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: (%s)%d: selected %s (%d), PSE price %g.",
		  rtnnme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  consys_nme(dy_sys,'v',*candxj,TRUE,NULL),
		  *candxj,candcbarj) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: (%s)%d: no suitable candidates.",rtnnme,
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; } }
  if (dy_opts->print.pricing >= 1)
  { if (retval == dyrPUNT)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: (%s)%d: all suitable x<j> on rejected pivot list.",
		  rtnnme,dy_prtlpphase(dy_lp->phase,TRUE),
		  dy_lp->tot.iters) ; } }
# endif

  return (retval) ; }



dyret_enum dy_primalpivot (int xjndx, int indir,
			   int *p_xindx, int *p_outdir,
			   double *p_abarij, double *p_delta, int *p_xjcand)

/*
  This routine executes a primal pivot. It calculates abar<j> = inv(B)a<j>
  (i.e., it ftran's the pivot column) and calls primalout to determine the
  leaving variable. It then calls dy_pivot (which in turn calls inv_update) to
  update the basis representation, and makes the necessary changes in the
  dylp data structures. In the course of updating primal and dual variables,
  we also update PSE pricing structures and choose a candidate to enter on
  the next pivot.

  There is a fairly elaborate antidegeneracy algorithm implemented here. When
  degeneracy is detected, a restricted problem is formed, composed of only
  the subset of constraints involved in the degeneracy. This problem is
  perturbed (radically!) and used in place of the original problem until a
  direction of recession is found (i.e., a pivot which results in a move away
  from the degenerate vertex). In the presence of finite upper and lower
  bounds, it's a bit tricky to distinguish a valid nondegenerate pivot from
  movement due to the perturbations, but we get by with some careful
  bookkeeping. The reference below discusses the algorithm in the standard
  (but overly simple) context of lower bounds of 0 and upper bounds of
  infinity.

  Note that there's an implicit assumption that the upper and lower bound
  on a variable will not be equal unless the variable's status is indicated
  as fixed (and hence it's ineligible for entering). Note also that while
  the antidegeneracy algorithm is active, the perturbation is applied to
  dy_xbasic but entries in dy_x are kept unperturbed so that we can quickly
  back out the perturbation.

  Reference:
    Ryan, D., Osborne, M., "On the Solution of Highly Degenerate Linear
    Programmes", Mathematical Programming, v.41, pp. 385-392, 1988.

  Parameters:
    xjndx:	The index of the entering variable.
    indir:	The direction of movement; +1 if x<j> is increasing, -1 if
		x<j> is decreasing.
    p_xindx:	(o) Index of the leaving variable x<i>. For return code
		    dyrLOSTPFEAS, set to the index of the variable where
		    primal feasibility loss was discovered
    p_outdir:	(o) returns the direction of movement for x<i>; +1 if it
		    increased and left at its upper bound, -1 if it decreased
		    and left at its lower bound.
    p_abarij:	(o) the pivot element abar<i,j>
    p_delta:	(o) the amount of change to x<j>
    p_xjcand:	(o) Index of the candidate entering variable for the next
		    pivot.

    The four output values are also valid for return codes dyrSINGULAR and
    dyrBSPACE, as they are all determined before inv_update is asked to pivot
    the basis.

  Returns: dyret_enum code, as follows:
    successful pivots:
    dyrOK:	The pivot completed successfully and a new candidate x<j>
		was selected.
    dyrDEGEN:	(primalout) As dyrOK, but the pivot was degenerate.
    dyrOPTIMAL:	(pseupdate) The pivot completed successfully, but no candidate
		x<j> could be found.
    dyrPUNT:    (pseupdate) The pivot completed successfully, but no
		candidate x<j> could be selected because all candidates were
		flagged with the NOPIVOT qualifier.
    dyrREQCHK:	(primalupdate) The pivot completed successfully, but a bogus
		number was calculated, or the pivot was dirty degenerate.

    unsuccessful (aborted) pivots:
    dyrMADPIV:	The pivot coefficient was judged (numerically) unstable
		(primalout, dy_pivot).
    dyrREQCHK:	The pivot coefficient is a bogus number (primalout), or
		there's been too much numerical drift while the antidegeneracy
		mechanism was active (degenout).
    dyrUNBOUND:	The problem is unbounded (primalout).
    dyrLOSTPFEAS: Primal feasibility has been lost (primalout).
    dyrSINGULAR: The pivot resulted in a singular basis (dy_pivot).
    dyrBSPACE:	basis package ran out of room to work (dy_pivot).
    dyrFATAL:	Fatal confusion (data structure error, internal confusion,
		etc.) (various sources)
*/

{ int xipos,xindx,xkpos,outdir ;
  double *abarj,*v,*betai,maxabarj,abarij,delta ;
  dyret_enum retval,outretval,pseretval ;
  bool reselect ;
  const char *rtnnme = "dy_primalpivot" ;

  extern dyret_enum primmultiout(int j, int indir, double *abarj,
				 double maxabarj, int *p_xindx,
				 int *p_outdir, double *p_deltaj) ;

/*
  Force invalid return values in case we abort before they're set normally.
  (This also avoids spurious `read from unitialized' errors.)
*/
  retval = dyrINV ;
  *p_xindx = -1 ;
  *p_outdir = 0 ;
  *p_xjcand = -1 ;
  *p_abarij = quiet_nan(0) ;
  *p_delta = quiet_nan(0) ;
/*
  First we do some prep work. Retrieve and ftran column a<j>. We need its max
  to use later to check for pivot stability.
*/
  abarj = NULL ;
  if (consys_getcol_ex(dy_sys,xjndx,&abarj) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,
	   "column",consys_nme(dy_sys,'v',xjndx,TRUE,NULL),xjndx) ;
    if (abarj != NULL) FREE(abarj) ;
    return (dyrFATAL) ; }

# ifndef DYLP_NDEBUG
/*
  Print the column, if the user's interested. There should be no dirty zeroes,
  and the print will expose them if they occur.
*/
  if (dy_opts->print.pivoting >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: x<%d> (%s) entering, status %s, %s from %g, ",
		rtnnme,xjndx,consys_nme(dy_sys,'v',xjndx,FALSE,NULL),
		dy_prtvstat(dy_status[xjndx]),
		(indir < 0)?"decreasing":"increasing",dy_x[xjndx]) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"lb<j> = %g, ub<j> = %g.",
	        dy_sys->vlb[xjndx],dy_sys->vub[xjndx]) ;
    if (dy_opts->print.pivoting >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    entering column a<%d>:",xjndx) ;
      xkpos = 1 ;
      for (xipos = 1 ; xipos <= dy_sys->concnt ; xipos++)
      { if (abarj[xipos] == 0) continue ;
	xkpos = (xkpos+1)%2 ;
	if (xkpos == 0) dyio_outchr(dy_logchn,dy_gtxecho,'\n') ;
	xindx = dy_basis[xipos] ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\t%ca<%d,%d> = %g",
	    (dy_lp->degen > 0 && dy_lp->degen == dy_degenset[xindx])?'*':'\0',
	    xipos,xjndx,abarj[xipos]) ; } } }
# endif

  dy_ftran(abarj,TRUE) ;
  maxabarj = exvec_infnorm(abarj,dy_sys->concnt,NULL) ;

# ifndef DYLP_NDEBUG
/*
  Print the ftran'd column. Again, there should be no dirty zeroes.
*/
  if (dy_opts->print.pivoting >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n    entering column abar<%d> = inv(B)a<%d>, max %g:",
	        xjndx,xjndx,maxabarj) ;
    xkpos = 1 ;
    for (xipos = 1 ; xipos <= dy_sys->concnt ; xipos++)
    { if (abarj[xipos] == 0) continue ;
      xkpos = (xkpos+1)%2 ;
      if (xkpos == 0) dyio_outchr(dy_logchn,dy_gtxecho,'\n') ;
      xindx = dy_basis[xipos] ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\t%ca<%d,%d> = %g",
             (dy_lp->degen > 0 && dy_lp->degen == dy_degenset[xindx])?'*':'\0',
	     xipos,xjndx,abarj[xipos]) ; } }
# endif
/*
  Open a loop to choose a leaving variable and perform a pivot. The reason we
  need a loop is because of the antidegeneracy algorithm. If it decides that
  this pivot will break us out of the degenerate vertex, we'll back out one
  level of restricted subproblem and repeat the selection process. If, on the
  other hand, primalout comes back with a degenerate pivot, we may want to
  form and perturb a restricted subproblem and reselect the pivot within that
  subproblem. Note that an indication from primalout that the problem is
  unbounded cannot be taken seriously while we're dealing with a restricted
  subproblem. The switch that follows is really spread out by comments, so
  watch for the comment indicating the end.
*/
  reselect = TRUE ;
  degen_cyclecnt = 0 ;
  xipos = -1 ;
  while (reselect)
  { if (dy_opts->ppsel.strat == 0)
    { outretval = primalout(xjndx,indir,abarj,maxabarj,
			    &xindx,&outdir,&delta) ; }
    else
    { outretval = primmultiout(xjndx,indir,abarj,maxabarj,
				 &xindx,&outdir,&delta) ; }
    switch (outretval)
    { 
/*
  primalout returns dyrUNBOUND

  Are we currently coping with degeneracy? If not, the problem is truly
  unbounded and we need to return to some higher level to deal with it.

  If there's a restricted subproblem installed, we've discovered a breakout
  direction from the degenerate vertex, and need to reselect the leaving
  variable after backing out the restricted subproblem. (Presumably we'll
  find a limiting variable in the full problem.)

  dy_degenout returns dyrREQCHK if it notices too much numerical drift of
  the current values in dy_x from the values when the degenerate subproblem
  was formed.
*/
      case dyrUNBOUND:
      { if (dy_lp->degen <= 0)
	{ FREE(abarj) ;
	  return (dyrUNBOUND) ; }
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.pivoting >= 1 || dy_opts->print.degen >= 1)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
	       "\n  (%s)%d: backing out level %d after %d pivots, unbounded.",
	       dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,dy_lp->degen,
	       dy_lp->tot.pivs-degenstats.iterin[dy_lp->degen]) ; }
#	endif
	if (dy_degenout(dy_lp->degen-1) != dyrOK)
	{ outretval = dyrREQCHK ;
	  reselect = FALSE ; }
	break ; }
/*
  primalout returns dyrOK

  Are we currently coping with degeneracy? If not, we have an uncomplicated,
  nondegenerate pivot. Yeah!
  
  If we are at a degenerate vertex, does this pivot break us out?  The trick
  is distinguishing a genuine bounded but nondegenerate pivot from a pivot
  that appears nondegenerate due to the perturbation in the restricted
  subproblem. To do this the code depends on dy_brkout, which specifies the
  required direction -- away from the bound where the variable entered the
  basis. By extension, a pivot where the entering variable swings to its
  opposite bound and leaves is also nondegenerate. Note that xipos is not
  valid when xjndx == xindx (x<j> isn't basic), so the order of the tests for
  breakout is important.

  If we've discovered a breakout direction, we need to reselect the leaving
  variable after backing out the restricted subproblem. (With more variables
  to consider and the perturbation gone, the limiting variable will likely
  be different.)

  dy_degenout returns dyrREQCHK if it notices too much numerical drift of
  the current values in dy_x from the values when the degenerate subproblem
  was formed.
*/
      case dyrOK:
      { if (xjndx != xindx) xipos = dy_var2basis[xindx] ;
	if (dy_lp->degen > 0 && (xjndx == xindx || outdir == dy_brkout[xipos]))
	{
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.pivoting >= 1 || dy_opts->print.degen >= 1)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		        "\n  (%s)%d: backing out level %d after %d pivots, ",
		        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		        dy_lp->degen,
		        dy_lp->tot.pivs-degenstats.iterin[dy_lp->degen]) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,"nondegenerate pivot.") ; }
#	  endif
	  if (dy_degenout(dy_lp->degen-1) != dyrOK)
	  { outretval = dyrREQCHK ;
	    reselect = FALSE ; } }
	else
	  reselect = FALSE ;
	break ; }
/*
  primalout returns dyrDEGEN

  Do we want (and are we allowed) to activate the antidegeneracy mechanism?
  If so, set up and perturb the restricted subproblem and then repeat the
  pivot selection. In order to create a restricted subproblem, opts.degen
  must permit it, and we must have executed opts.degenpivlim successive
  degenerate and nonconstructive pivots.

  The idea is to activate the antidegeneracy algorithm only when we have
  serious degeneracy involving explicit constraints where we can perturb the
  right-hand side (which we accomplish by the equivalent action of
  perturbing the values of the basic variables). To this end, we exclude
  degenerate pivots where:
  * A fixed variable is leaving. A fixed variable is certainly part of the
    cause of degeneracy, and the rules for selecting an incoming variable
    guarantee it won't come back.
  * A free variable is entering. Free variables have no bound, hence can't be
    a cause of degeneracy, and the rules for selecting a leaving variable
    guarantee that a free variable will never leave the basis. And nonbasic
    free variables preclude dual feasibility.
  * A superbasic variable is entering. Since it's not at bound, it won't be
    an immediate cause of degeneracy. And superbasic variables preclude dual
    feasibility.
  * The leaving variable is flagged as `do not perturb'. For all intents and
    purposes, it might as well be fixed, and we need to get it out of the
    basis.
  We'll allow only three attempts at getting the perturbation right, then we
  take the degenerate pivot and be done with it.
*/
      case dyrDEGEN:
      { if (flgon(dy_status[xindx],vstatBFX|vstatNOPER) ||
	    flgon(dy_status[xjndx],vstatNBFR|vstatSB))
	{ reselect = FALSE ;
	  xipos = dy_var2basis[xindx] ;
#         ifndef DYLP_NDEBUG
	  if (dy_opts->print.degen >= 3)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		        "\n      (%s)%d: constructive degenerate pivot.",
		        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s %s (%d) leaving,",
		        dy_prtvstat(dy_status[xindx]),
		        consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho," %s %s (%d) entering.",
		        dy_prtvstat(dy_status[xjndx]),
		        consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx) ; }
#	  endif
	}
	else
        if (dy_opts->degen == TRUE &&
	    dy_opts->degenpivlim < dy_lp->degenpivcnt &&
	    degen_cyclecnt < 3)
	{ 
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.pivoting >= 1 || dy_opts->print.degen >= 1)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		        "\n  (%s)%d: antidegeneracy increasing to level %d.",
		        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		        dy_lp->degen+1) ; }
#	  endif
	  dy_degenin() ;
	  degen_cyclecnt++ ; }
	else
	{ reselect = FALSE ;
	  xipos = dy_var2basis[xindx] ;
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.degen >= 2 && degen_cyclecnt >= 3)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		   "\n    (%s)%d: forced degenerate pivot after %d cycles;",
		   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		   degen_cyclecnt) ;
	    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s %s (%d) leaving.",
		        dy_prtvstat(dy_status[xindx]),
		        consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx) ; }
	  else
	  if (dy_opts->print.degen >= 3)
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
		   "\n      (%s)%d: degenerate pivot, %s %s (%d) leaving.",
		   dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		   dy_prtvstat(dy_status[xindx]),
		   consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx) ; }
#	  endif
	}
	break ; }
/*
  Remaining cases, and end of the pivot selection loop. If primalout returned
  anything other than dyrOK, dyrUNBOUND, or dyrDEGEN, it'll fall through to
  here and we'll punt back to the caller. Possibilities are dyrREQCHK,
  dyrMADPIV, dyrLOSTPFEAS and dyrFATAL. Make the default a trap for internal
  confusion.
*/
      case dyrMADPIV:
      { xipos = dy_var2basis[xindx] ;
	abarij =  abarj[xipos] ;
	(void) dy_addtopivrej(xjndx,dyrMADPIV,abarij,maxabarj) ;
	reselect = FALSE ;
	break ; }
      case dyrREQCHK:
      case dyrLOSTPFEAS:
      { xipos = dy_var2basis[xindx] ;
	reselect = FALSE ;
	break ; }
      case dyrFATAL:
      { FREE(abarj) ;
	return (outretval) ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	FREE(abarj) ;
	return (outretval) ; } } }
/*
  Set the various output parameters, then return for everthing except dyrOK
  and dyrDEGEN (the second case arising because we've chosen not to activate
  the antidegeneracy mechanism). We'll claim 1.0 as the pivot coefficient if
  the pivot is a nonbasic variable swinging bound-to-bound. Note that we can
  wind up out here with outretval = dyrREQCHK and xjndx != xindx = -1 if
  primalout returned dyrUNBOUND, then dy_degenout ran into trouble and
  returned dyrREQCHK.
*/
  *p_xindx = xindx ;
  *p_outdir = outdir ;
  *p_delta = delta ;
  abarij = quiet_nan(0) ;
  if (xjndx != xindx)
  { if (xindx > 0)
    { abarij = abarj[xipos] ;
      *p_abarij = abarij ; } }
  else
  { *p_abarij = 1.0 ; }
  if (!(outretval == dyrOK || outretval == dyrDEGEN))
  { FREE(abarj) ;
    return (outretval) ; }
/*
  The notion is that BFX will never reenter and NBFR will never leave, once the
  pivot is complete. So we're in no danger of cycling.
*/
  if (outretval == dyrOK)
  { dy_lp->degenpivcnt = 0 ; }
  else
  { if (flgon(dy_status[xindx],vstatBFX|vstatNOPER) ||
	flgon(dy_status[xjndx],vstatNBFR|vstatSB))
      dy_lp->degenpivcnt = 0 ;
    else
      dy_lp->degenpivcnt++ ; }
/*
  It looks like the pivot will go through, so get down to business.  Updating
  the PSE and pricing information (dy_cbar and dy_gamma) will require beta<i>
  (row i of the basis inverse), and the vector v = abar~<j>inv(B). We need to
  calculate these prior to the basis change. Unfortunately, we'll then need to
  recalculate inv(B)abar<j> for primal updates.

  After the prep work, attempt the pivot to update the LU factorisation. This
  can fail for three reasons: the pivot element didn't meet the numerical
  stability criteria (but we've checked this already), the pivot produced a
  singular basis, or the basis package ran out of space.

  On the dual side, it was a big win, computationally, to attempt to salvage
  the pivot at this point with a refactor if dy_pivot reported dyrNUMERIC
  (near singularity). It's not clear that this is a problem on the primal side,
  but if I find myself staring at this bit of code, it's worth a shot. All this
  code should be pulled out to a small subroutine if I add recovery.
*/
  if (xjndx != xindx)
  { for (xkpos = 1 ; xkpos <= dy_sys->concnt ; xkpos++)
      if (dy_frame[dy_basis[xkpos]] == FALSE) abarj[xkpos] = 0 ;
    dy_btran(abarj) ;
    v = (double *) MALLOC((dy_sys->concnt+1)*sizeof(double)) ;
    memcpy(v,abarj,(dy_sys->concnt+1)*sizeof(double)) ;
    betai = (double *) CALLOC((dy_sys->concnt+1),sizeof(double)) ;
    betai[xipos] = 1.0 ;
    dy_btran(betai) ;
    if (consys_getcol_ex(dy_sys,xjndx,&abarj) == FALSE)
    { errmsg(122,rtnnme,dy_sys->nme,"column",
	     consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx) ;
      retval = dyrFATAL ; }
    else
    { dy_ftran(abarj,TRUE) ;
      retval = dy_pivot(xipos,abarij,maxabarj) ; } }
  else
  { v = NULL ;
    betai = NULL ;
    retval = dyrOK ; }
/*
  Then update the dylp data structures -- primalupdate does the basis,
  status, and primal and dual variable values, pseupdate does the column
  norms and reduced costs. In the process of updating the reduced costs,
  we'll select a candidate for the next entering variable.

  Note that there's no pseupdate if the pivot was a bound-to-bound move by a
  nonbasic variable (no basis change, no pseupdate).

  Finally, decide on the appropriate return value. Assuming pseupdate doesn't
  run into trouble, dyrREQCHK (primalupdate) wins over dyrDEGEN (primalout)
  which wins over dyrOK (various).
*/
  if (retval == dyrOK)
  { dy_lp->pivok = TRUE ;
    retval = primalupdate(xjndx,indir,xindx,outdir,abarj,delta,betai) ;
    if (retval == dyrOK || retval == dyrREQCHK || retval == dyrSWING)
    { if (xjndx != xindx)
      { pseretval = pseupdate(xjndx,xindx,p_xjcand,abarj,v,betai) ; }
      else
      { pseretval = dyrOK ;
	*p_xjcand = 0 ; } }
    else
    { pseretval = dyrOK ; }
    if (pseretval != dyrOK)
      retval = pseretval ;
    else
    if (retval == dyrOK)
    { if (outretval == dyrDEGEN) retval = dyrDEGEN ; } }
  else
  if (retval == dyrNUMERIC)
  { retval = dyrSINGULAR ; }
/*
  Tidy up and return.
*/
  FREE(abarj) ;
  if (v != NULL) FREE(v) ;
  if (betai != NULL) FREE(betai) ;

  return (retval) ; }


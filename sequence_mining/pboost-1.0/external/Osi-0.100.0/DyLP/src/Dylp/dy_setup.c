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
  This file contains routines to establish options and tolerances for the dylp
  package.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_setup.c	4.7	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_setup.c 269 2009-04-02 05:38:19Z lou $" ;

/*
  dyopts_dflt

  This structure holds the default values assigned prior to processing the
  user's option specifications. A value of -1 indicates that there is no
  default. Either we need to know more about the constraint system before
  choosing a default, or it's something only the user can tell us.
  
  The companion structures dyopts_lb and dyopts_ub give allowable ranges
  (where relevant) and are used to make sure the user doesn't give us nonsense
  values. They're sparsely populated -- it's not possible (or useful) to give
  bounds for a fair number of the values. An (r) indicates a recommended
  value -- the code will gripe about a violation, but won't enforce the limit.

  `Chvatal' is Chvatal, V., Linear Programming, W.H. Freeman, 1983.

  Some specific comments:

  scan:         Controls the number of variables scanned when the incoming
		variable selected during PSE updating is rejected due to
		pivoting problems and we need to scan for alternates.
  iterlim:      The value is set to 10000 <= 5*concnt in the absence of help
		from the user. Comments on pp. 45-46 of Chvatal indicate that
		3*concnt should be sufficient for all but the worst problems.
		A value of 0 means no limit is enforced.
  idlelim:      The value is set to 1000 <= 5*concnt <= 10000 in the absence
		of help from the user. (During dual simplex, it's enforced at
		80% of that, since the dual doesn't have antidegeneracy
		installed.)
  factor:	The default choice of 20 comes from the analysis on
		pp. 111 - 113 of Chvatal, "Linear Programming". The upper
		bound of 100 comes from a passing mention in the XMP
		documentation that factor > 100 shouldn't be used except
		for network problems. Given upgrades since the original version
		of the code, this should be relaxed, but I need to do some
		experiments to see just how far.
  check:	Check actually defaults to factor/2, with an upper limit of
		factor.
  active.*:	The value is for efficiency, not a limit. If it's too small,
		the dylp constraint system will have to do more expansions.
		25% is a pure guess. I'll try and refine it with experience.
  print.*:	The majority of the print controls are set to 0. The output
		is really only of interest if you're looking at some detail
		of the lp implementation.
*/

static
lpopts_struct dyopts_dflt = { cxINITIALLP, /* context */
			      -1,	/* scan */
			      -1,	/* iterlim */
			      -1,	/* idlelim */
			    { -1,	/* dpsel.strat */
			      TRUE,	/* dpsel.flex */
			      FALSE },	/* dpsel.allownopiv */
			    {  1 },	/* ppsel.strat */
			      50,	/* factor */
			      -1,	/* check */
			       1,	/* groom */
			     { 1,	/* con.actlvl */
			       0,	/* con.actlim */
			       0 },	/* con.deactlvl */
			       0,	/* addvar */
			       3,	/* dualadd */
			       5000,	/* coldvars */
			      FALSE,	/* forcecold */
			      FALSE,	/* forcewarm */
			      TRUE,	/* usedual */
			      TRUE,	/* degen */
			      1,	/* degenpivlim */
			      1,	/* degenlite */
			      TRUE,	/* patch */
			      FALSE,	/* fullsys */
			      FALSE,	/* copyorigsys */
			      2,	/* scaling */
			      { .25,	/* active.vars */
				.25 },	/* active.cons */
			      { .5,	/* initcons.frac */
				TRUE,	/* initcons.i1lopen */
				90,	/* initcons.i1l */
				FALSE,	/* initcons.i1uopen */
				180,	/* initcons.i1u */
				TRUE,	/* initcons.i2valid */
				FALSE,	/* initcons.i2lopen */
				0,	/* initcons.i2l */
				TRUE,	/* initcons.i2uopen */
				90 },	/* initcons.i2u */
			      ibLOGICAL, /* coldbasis */
			      { TRUE,	/* finpurge.cons */
				TRUE },	/* finpurge.vars */
			      { FALSE,	/* heroics.d2p */
				FALSE,	/* heroics.p2d */
				FALSE }, /* heroics.flips */
			      { 0,	/* print.major */
			        0,	/* print.scaling */
			        0,	/* print.setup */
				0,	/* print.crash */
				0,	/* print.pricing */
				0,	/* print.pivoting */
				0,	/* print.pivreject */
				0,	/* print.degen */
				1,	/* print.phase1 */
				1,	/* print.phase2 */
				1,	/* print.dual */
				1,	/* print.basis */
				0,	/* print.conmgmt */
				0,	/* print.varmgmt */
				1, 	/* print.force */
				0,	/* print.tableau */
				0, 	/* print.rays */
				0 	/* print.soln */
			      }
			    } ;

static
lpopts_struct dyopts_lb = { cxSINGLELP,	/* context */
			    200,	/* scan */
			    0,		/* iterlim */
			    0,		/* idlelim */
			    { 0,	/* dpsel.strat */
			      FALSE,	/* dpsel.flex */
			      FALSE },	/* dpsel.allownopiv */
			    {  0 },	/* ppsel.strat */
			    1,		/* factor */
			    1,		/* check */
			    0,		/* groom */
			     { 0,	/* con.actlvl */
			       0,	/* con.actlim */
			       0 },	/* con.deactlvl */
			    0,		/* addvar */
			    0,		/* dualadd */
			    0,		/* coldvars */
			    FALSE,	/* forcecold */
			    FALSE,	/* forcewarm */
			    FALSE,	/* usedual */
			    FALSE,	/* degen */
			    0,		/* degenpivlim */
			    0,		/* degenlite */
			    FALSE,	/* patch */
			    FALSE,	/* fullsys */
			    FALSE,	/* copyorigsys */
			    0,		/* scaling */
			    { 0.0,	/* active.vars */
			      0.0 },	/* active.cons */
			    { 0.0,	/* initcons.frac */
			      FALSE,	/* initcons.i1lopen */
			      0,	/* initcons.i1l */
			      FALSE,	/* initcons.i1uopen */
			      0,	/* initcons.i1u */
			      FALSE,	/* initcons.i2valid */
			      FALSE,	/* initcons.i2lopen */
			      0,	/* initcons.i2l */
			      FALSE,	/* initcons.i2uopen */
			      0 },	/* initcons.i2u */
			      ibLOGICAL, /* coldbasis */
			      { FALSE,	/* finpurge.cons */
				FALSE }, /* finpurge.vars */
			    { FALSE,	/* heroics.d2p */
			      FALSE,	/* heroics.p2d */
			      FALSE },	/* heroics.flips */
			    { 0,	/* print.major */
			      0,	/* print.scaling */
			      0,	/* print.setup */
			      0,	/* print.crash */
			      0, 	/* print.pricing */
			      0, 	/* print.pivoting */
			      0,	/* print.pivreject */
			      0,	/* print.degen */
			      0,	/* print.phase1 */
			      0,	/* print.phase2 */
			      0,	/* print.dual */
			      0,	/* print.basis */
			      0,	/* print.conmgmt */
			      0,	/* print.varmgmt */
			      0,	/* print.force */
			      0,	/* print.tableau */
			      0, 	/* print.rays */
			      0 	/* print.soln */
			    }
			  } ;

/*
  Roughly, MAXITERLIM = MAXINT/4, so we can set the overall iteration limit
  to 3*iterlim without getting into integer overflow.
*/

#define MAXITERLIM (int) (((unsigned) 1<<(sizeof(dyopts_ub.iterlim)*8-3))-1)

static
lpopts_struct dyopts_ub = { cxBANDC,	/* context */
			    1000,	/* scan */
			    MAXITERLIM,	/* iterlim */
			    MAXITERLIM,	/* idlelim */
			    { 3,	/* dpsel.strat */
			      TRUE,	/* dpsel.flex */
			      TRUE },	/* dpsel.allownopiv */
			    {  1 },	/* ppsel.strat */
			    100,	/* factor (r) */
			    -1,		/* check == factor */
			     2,		/* groom */
			    {  1,	/* con.actlvl */
			      -1,	/* con.actlim */
			       2 },	/* con.deactlvl */
			    -1,		/* addvar */
			     3,		/* dualadd */
			    100000,	/* coldvars (r) */
			    TRUE,	/* forcecold */
			    TRUE,	/* forcewarm */
			    TRUE,	/* usedual */
			    TRUE,	/* degen */
			    -1,		/* degenpivlim */
			    5,		/* degenlite */
			    TRUE,	/* patch */
			    TRUE,	/* fullsys */
			    TRUE,	/* copyorigsys */
			    2,		/* scaling */
			    { 1.0,	/* active.vars */
			      1.0 },	/* active.cons */
			    { 1.0,	/* initcons.frac */
			      TRUE,	/* initcons.i1lopen */
			      180,	/* initcons.i1l */
			      TRUE,	/* initcons.i1uopen */
			      180,	/* initcons.i1u */
			      TRUE,	/* initcons.i2valid */
			      TRUE,	/* initcons.i2lopen */
			      180,	/* initcons.i2l */
			      TRUE,	/* initcons.i2uopen */
			      180 },	/* initcons.i2u */
			      ibARCH,	/* coldbasis */
			      { TRUE,	/* finpurge.cons */
				TRUE },	/* finpurge.vars */
			    { TRUE,	/* heroics.d2p */
			      TRUE,	/* heroics.p2d */
			      TRUE },	/* heroics.flips */
			    { 1,	/* print.major */
			      2,	/* print.scaling */
			      5,	/* print.setup */
			      4,	/* print.crash */
			      3,	/* print.pricing */
			      5,	/* print.pivoting */
			      2,	/* print.pivreject */
			      5,	/* print.degen */
			      7,	/* print.phase1 */
			      7,	/* print.phase2 */
			      7,	/* print.dual */
			      5,	/* print.basis */
			      5,	/* print.conmgmt */
			      4,	/* print.varmgmt */
			      3,	/* print.force */
			      6,	/* print.tableau */
			      4, 	/* print.rays */
			      4 	/* print.soln */
			    }
			  } ;


/*
  dytols_dflt

  This structure holds the default values assigned prior to processing the
  user's tolerance specifications.

  Some specific comments:

  inf:          Infinity. Dylp can work with either IEEE infinity or a
		`finite' infinity (most often, DBL_MAX). The default is
		HUGE_VAL, which will resolve to IEEE infinity in a
		Sun/Solaris environment.  HUGE_VAL isn't always a
		compile-time constant, so we load it in dy_defaults. If the
		client code has a finite infinity, you surely want to pass
		this in to dylp. Otherwise, dylp will hand back IEEE
		infinity, and finite and infinite infinities do not play well
		together.
  zero:         Defaults to 1.0e-11. Historically I've seen it between
		1.0e-10 and 1.0e-12.  A colleague offers the following rule
		of thumb: ``The zero tolerance should be the product of the
		machine accuracy and the largest number you expect to
		encounter during processing.'' Since the 64 bit IEEE floating
		point format has a 52 bit mantissa, the machine precision is
		2^-52, or about 10^-15. 1.0e-11 is reasonable by this rule.
  pfeas:	This value will be scaled by an amount proportional to the
		1-norm of the basic variables, with a minimum value of zero.
		It's set to (pfeas_scale)*(zero tolerance) initially (right
		at the start of dylp), so that we have something that's valid
		when establishing the basis.
  pfeas_scale:	Allows decoupling of pfeas from zero. Defaults to 10, but
		can be adjusted by the user.
  cost:		The moral equivalent of 0 for things related to the objective
		function, reduced costs, dual variables, etc. Defaults to
		1.0e-10.
  dfeas:	This value will be scaled by an amount proportional to the
		1-norm of the dual variables, with a minimum value of cost.
		There is no default held here.
  dfeas_scale:	As pfeas_scale, for dfeas.
  pivot:        This is the pivot acceptance tolerance, as a fraction of the
		pivot selection tolerance used by the basis package during
		factoring. Defaults to 1.0e-5.
  bogus:        This multiplier is used to detect `bogus' values for reduced
		costs and variable values. Bogus values are values which
		exceed a zero tolerance but are less than bogus*tolerance.
		The idea is that, generally speaking, numbers close to a zero
		tolerance are problematic, and may well be the result of
		accumulated numerical inaccuracy.  dylp attempts to nip this
		problem in the bud by watching for numbers such that tol <
		|val| < bogus*tol for reduced costs and variable values and
		requesting refactorisation when it sees one. Defaults to 1.
		Higher values prompt more refactoring.
  swing:	When (new primal value)/(old primal value) > swing, dylp
		takes it as an indication that the problem needs more
		constraints (pseudo-unbounded).
  reframe:      Multiplier used to trigger a PSE or DSE reference framework
		reset. The default is .1, based on the computational results
		reported for the primal simplex in Forrest & Goldfarb,
		"Steepest Edge Algorithms for Linear Programming",
		Mathematical Programming v.57, pp. 341--374, 1992.
*/

static
lptols_struct dytols_dflt = { 0,		/* inf = HUGE_VAL */
			      1.0e-11,		/* zero */
			      1.0e-5, 		/* pchk */
			      -1, 		/* pfeas */
			      100, 		/* pfeas_scale */
			      1.0e-11,		/* cost	*/
			      1.0e-4,		/* dchk */
			      -1,		/* dfeas */
			      100,		/* dfeas_scale */
			      1.0e-5,		/* pivot */
			      1,		/* bogus */
			      1.0e15,		/* swing */
			      1.0e30,		/* toobig */
			      1.0e-4,		/* purge */
			      .5,		/* purgevar */
			      .1		/* reframe */
			    } ;


void dy_exposeOptDefaults (lpopts_struct **opts_lb,
			   lpopts_struct **opts_dflt,
			   lpopts_struct **opts_ub)
/*
  The sole purpose of this routine is to allow other parts of the code to get
  hold of the defaults without exposing them as global variables. At present,
  used only by dy_options.c.
*/

{ if (opts_lb != NULL) *opts_lb = &dyopts_lb ;
  if (opts_dflt != NULL) *opts_dflt = &dyopts_dflt ;
  if (opts_ub != NULL) *opts_ub = &dyopts_ub ; }

void dy_exposeTolDefaults (lptols_struct **tols_dflt)

/*
  As for dy_exposeOptDefaults --- make default tolerances available to the
  rest of the code.
*/

{ if (tols_dflt != NULL) *tols_dflt = &dytols_dflt ; }



void dy_defaults (lpopts_struct **opts, lptols_struct **tols)

/*
  This routine loads the data structures supplied as parameters with default
  tolerances and options. Typically, a client calls this routine to establish
  base values, then tweaks the options and parameters as desired.

  Parameters:
    opts	(i) lpopts structure; allocated if NULL
		(o) returns an lpopts structure
    tols	(i) lptols structure ; allocated if NULL
		(o) returns an lptols structure

  Returns: undefined
*/

{ 

#ifdef DYLP_PARANOIA
  const char *rtnnme = "dy_defaults" ;

  if (opts == NULL)
  { errmsg(2,rtnnme,"&opts") ;
    return ; }
  if (tols == NULL)
  { errmsg(2,rtnnme,"&tols") ;
    return ; }
#endif

  if (*opts == NULL)
  { (*opts) = (lpopts_struct *) MALLOC(sizeof(lpopts_struct)) ; }
  memcpy(*opts,&dyopts_dflt,sizeof(lpopts_struct)) ;

  if (*tols == NULL)
  { (*tols) = (lptols_struct *) MALLOC(sizeof(lptols_struct)) ; }
  memcpy(*tols,&dytols_dflt,sizeof(lptols_struct)) ;
  (*tols)->inf = HUGE_VAL ;

  return ; }



void dy_checkdefaults (consys_struct *sys,
		       lpopts_struct *opts, lptols_struct *tols)

/*
  This routine looks over various option and tolerance settings with an eye
  toward setting or adjusting values based on the size of the constraint
  system or other options that might be set by the user.

  The default values here are, by and large, grossly larger than required.
  The more outrageous ones are motivated by the Netlib examples.

  Parameters:
    sys:	a constraint system
    opts:	an options structure; may be modified on return
  
  Returns: undefined
*/

{ int scalefactor ;

  if (opts->check < 0) opts->check = opts->factor/2 ;
  if (opts->check <= 0) opts->check = 1 ;

  if (opts->scan < 0)
  { opts->scan = maxx(dyopts_lb.scan,sys->archvcnt/2) ;
    opts->scan = minn(opts->scan,dyopts_ub.scan) ; }

  if (opts->iterlim < 0)
  { opts->iterlim = minn(5*(sys->concnt+sys->varcnt),100000) ;
    opts->iterlim = maxx(opts->iterlim,10000) ; }

  if (opts->idlelim < 0)
  { opts->idlelim = minn(2*(sys->concnt+sys->varcnt),50000) ;
    opts->idlelim = maxx(opts->idlelim,1000) ; }

  if (opts->degenpivlim < 0)
  { opts->degenpivlim = minn(1000,sys->concnt/2) ;
    opts->degenpivlim = maxx(100,opts->degenpivlim) ; }

/*
  If the user has specified a dual pivot strategy, observe it. If not, start
  with strategy 1 (max dual objective improvement).
*/
  if (opts->dpsel.strat >= 0)
  { opts->dpsel.flex = FALSE ; }
  else
  { opts->dpsel.strat = 1 ;
    opts->dpsel.flex = TRUE ; }

  if (opts->fullsys == TRUE)
  { opts->active.vars = 1.0 ;
    opts->active.cons = 1.0 ; }
/*
  Loosen the base primal and dual accuracy check values for larger systems,
  and put a little more distance between the zero tolerances and feasibility
  tolerances.
*/
  scalefactor = ((int) (.5 + log10((double) sys->varcnt))) - 2 ;
  if (scalefactor > 0)
  { tols->pchk *= pow(10.0,(double) scalefactor) ;
    tols->pfeas_scale *= pow(10.0,(double) scalefactor) ; }
  scalefactor = ((int) (.5 + log10((double) sys->concnt))) - 2 ;
  if (scalefactor > 0)
  { tols->dchk *= pow(10.0,(double) scalefactor) ;
    tols->dfeas_scale *= pow(10.0,(double) scalefactor) ; }

/*
  XX_DEBUG_XX

  There's no good way to control this print statement, given the timing and
  purpose of this call. But it's occasionally handy for debugging.

  dyio_outfmt(dy_logchn,TRUE,"\nPTOLS: pzero = %g, pscale = %g, pchk = %g",
	      tols->zero,tols->pfeas_scale,tols->pchk) ;
  dyio_outfmt(dy_logchn,TRUE,"\nDTOLS: dzero = %g, dscale = %g, dchk = %g",
	      tols->cost,tols->dfeas_scale,tols->dchk) ;
*/

  return ; }



void dy_setprintopts (int lvl, lpopts_struct *opts)

/*
  This routine tweaks the lp print level options based on a single integer
  code. It's intended to allow clients of dylp to easily set up overall print
  levels. Just a big case statement.
  
  Level 0 forces dylp to shut up.
  
  Level 1 assumes the normal dylp defaults (phase1, phase2, dual, force, and
  basis print levels set to 1, which allows messages about extraordinary
  events).

  Levels 2 -- 5 provide increasing amounts of information.
  
  At level 1 and above, a specific setting of a dylp value to a higher value
  in the options structure passed in as a parameter will override the value
  here.

  Parameters:
    lvl:	overall print level
    opts:	options structure; for all except lvl = 0, should be preloaded
		with valid values for print options.
  
  Returns: undefined
*/

{ if (lvl < 0) lvl = 0 ;

  switch (lvl)
  { case 0:
    { opts->print.major = 0 ;
      opts->print.scaling = 0 ;
      opts->print.setup = 0 ;
      opts->print.crash = 0 ;
      opts->print.pricing = 0 ;
      opts->print.pivoting = 0 ;
      opts->print.pivreject = 0 ;
      opts->print.degen = 0 ;
      opts->print.phase1 = 0 ;
      opts->print.phase2 = 0 ;
      opts->print.dual = 0 ;
      opts->print.basis = 0 ;
      opts->print.conmgmt = 0 ;
      opts->print.varmgmt = 0 ;
      opts->print.force = 0 ;
      opts->print.tableau = 0 ;
      opts->print.rays = 0 ;
      opts->print.soln = 0 ;
      break ; }
    case 1:
    { opts->print.major = maxx(opts->print.major,dyopts_dflt.print.major) ;
      opts->print.scaling =
	  maxx(opts->print.scaling,dyopts_dflt.print.scaling) ;
      opts->print.setup = maxx(opts->print.setup,dyopts_dflt.print.setup) ;
      opts->print.crash = maxx(opts->print.crash,dyopts_dflt.print.crash) ;
      opts->print.pricing = maxx(opts->print.pricing,
				 dyopts_dflt.print.pricing) ;
      opts->print.pivoting = maxx(opts->print.pivoting,
				  dyopts_dflt.print.pivoting) ;
      opts->print.pivreject = maxx(opts->print.pivreject,
				  dyopts_dflt.print.pivreject) ;
      opts->print.degen = maxx(opts->print.degen,dyopts_dflt.print.degen) ;
      opts->print.phase1 = maxx(opts->print.phase1,dyopts_dflt.print.phase1) ;
      opts->print.phase2 = maxx(opts->print.phase2,dyopts_dflt.print.phase2) ;
      opts->print.dual = maxx(opts->print.dual,dyopts_dflt.print.dual) ;
      opts->print.basis = maxx(opts->print.basis,dyopts_dflt.print.basis) ;
      opts->print.conmgmt = maxx(opts->print.conmgmt,
				 dyopts_dflt.print.conmgmt) ;
      opts->print.varmgmt = maxx(opts->print.varmgmt,
				 dyopts_dflt.print.varmgmt) ;
      opts->print.force = maxx(opts->print.force,dyopts_dflt.print.force) ;
      opts->print.tableau = maxx(opts->print.tableau,
				 dyopts_dflt.print.tableau) ;
      opts->print.rays = maxx(opts->print.rays,
				 dyopts_dflt.print.rays) ;
      opts->print.soln = maxx(opts->print.soln,
				 dyopts_dflt.print.soln) ;
      break ; }
    case 2:
    { opts->print.major = maxx(opts->print.major,1) ;
      opts->print.scaling = maxx(opts->print.scaling,1) ;
      opts->print.setup = maxx(opts->print.setup,1) ;
      opts->print.crash = maxx(opts->print.crash,1) ;
      opts->print.pricing = maxx(opts->print.pricing,0) ;
      opts->print.pivoting = maxx(opts->print.pivoting,0) ;
      opts->print.pivreject = maxx(opts->print.pivreject,0) ;
      opts->print.degen = maxx(opts->print.degen,1) ;
      opts->print.phase1 = maxx(opts->print.phase1,1) ;
      opts->print.phase2 = maxx(opts->print.phase2,1) ;
      opts->print.dual = maxx(opts->print.dual,1) ;
      opts->print.basis = maxx(opts->print.basis,1) ;
      opts->print.conmgmt = maxx(opts->print.conmgmt,1) ;
      opts->print.varmgmt = maxx(opts->print.varmgmt,1) ;
      opts->print.force = maxx(opts->print.force,1) ;
      opts->print.tableau = maxx(opts->print.tableau,1) ;
      opts->print.rays = maxx(opts->print.rays,1) ;
      opts->print.soln = maxx(opts->print.soln,1) ;
      break ; }
    case 3:
    { opts->print.major = maxx(opts->print.major,1) ;
      opts->print.scaling = maxx(opts->print.scaling,2) ;
      opts->print.setup = maxx(opts->print.setup,2) ;
      opts->print.crash = maxx(opts->print.crash,2) ;
      opts->print.pricing = maxx(opts->print.pricing,0) ;
      opts->print.pivoting = maxx(opts->print.pivoting,0) ;
      opts->print.pivreject = maxx(opts->print.pivreject,0) ;
      opts->print.degen = maxx(opts->print.degen,1) ;
      opts->print.phase1 = maxx(opts->print.phase1,3) ;
      opts->print.phase2 = maxx(opts->print.phase2,3) ;
      opts->print.dual = maxx(opts->print.dual,3) ;
      opts->print.basis = maxx(opts->print.basis,2) ;
      opts->print.conmgmt = maxx(opts->print.conmgmt,2) ;
      opts->print.varmgmt = maxx(opts->print.varmgmt,2) ;
      opts->print.force = maxx(opts->print.force,1) ;
      opts->print.tableau = maxx(opts->print.tableau,1) ;
      opts->print.rays = maxx(opts->print.rays,1) ;
      opts->print.soln = maxx(opts->print.soln,1) ;
      break ; }
    case 4:
    { opts->print.major = maxx(opts->print.major,1) ;
      opts->print.scaling = maxx(opts->print.scaling,2) ;
      opts->print.setup = maxx(opts->print.setup,3) ;
      opts->print.crash = maxx(opts->print.crash,3) ;
      opts->print.pricing = maxx(opts->print.pricing,0) ;
      opts->print.pivoting = maxx(opts->print.pivoting,0) ;
      opts->print.pivreject = maxx(opts->print.pivreject,0) ;
      opts->print.degen = maxx(opts->print.degen,2) ;
      opts->print.phase1 = maxx(opts->print.phase1,4) ;
      opts->print.phase2 = maxx(opts->print.phase2,4) ;
      opts->print.dual = maxx(opts->print.dual,4) ;
      opts->print.basis = maxx(opts->print.basis,3) ;
      opts->print.conmgmt = maxx(opts->print.conmgmt,3) ;
      opts->print.varmgmt = maxx(opts->print.varmgmt,2) ;
      opts->print.force = maxx(opts->print.force,2) ;
      opts->print.tableau = maxx(opts->print.tableau,1) ;
      opts->print.rays = maxx(opts->print.rays,3) ;
      opts->print.soln = maxx(opts->print.soln,3) ;
      break ; }
    default:
    { opts->print.major = maxx(opts->print.major,1) ;
      opts->print.scaling = maxx(opts->print.scaling,2) ;
      opts->print.setup = maxx(opts->print.setup,5) ;
      opts->print.crash = maxx(opts->print.crash,4) ;
      opts->print.pricing = maxx(opts->print.pricing,1) ;
      opts->print.pivoting = maxx(opts->print.pivoting,1) ;
      opts->print.pivreject = maxx(opts->print.pivreject,1) ;
      opts->print.degen = maxx(opts->print.degen,2) ;
      opts->print.phase1 = maxx(opts->print.phase1,5) ;
      opts->print.phase2 = maxx(opts->print.phase2,5) ;
      opts->print.dual = maxx(opts->print.dual,5) ;
      opts->print.basis = maxx(opts->print.basis,5) ;
      opts->print.conmgmt = maxx(opts->print.conmgmt,3) ;
      opts->print.varmgmt = maxx(opts->print.varmgmt,2) ;
      opts->print.force = maxx(opts->print.force,3) ;
      opts->print.tableau = maxx(opts->print.tableau,4) ;
      opts->print.rays = maxx(opts->print.rays,4) ;
      opts->print.soln = maxx(opts->print.soln,4) ;
      break ; } }
  
  return ; }


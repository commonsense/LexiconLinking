/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dy_dualmultipivot.c	4.6	10/15/05" ;
static char svnid[] UNUSED = "$Id: dy_dualmultipivot.c 269 2009-04-02 05:38:19Z lou $" ;

/*
  This is an experimental dual pivoting algorithm. The inspiration comes from
  two places:

    * Observing the `dual death spiral', where the dual simplex executes a
      succession of pivots that result in steadily increasing (and infeasible)
      values of the primal variables. Each pivot is driving one variable to
      bound, but at the cost of driving other primal variables farther from
      bound. It's still unclear to me how/why this occurs. It's not simply a
      side effect of working with the partial system --- grow22 is a
      pathological case, and it's all equalities, loaded in full. But it
      leads to the thought ``Wouldn't it be nice to choose the dual pivot to
      reduce primal infeasibility?'' After all, attaining primal feasibility
      is equivalent to reaching optimality, given we start with dual
      feasibility.

    * Looking over the dual pivot selection algorithm in clp. Instead of just
      choosing the variable x<j> associated with the most restrictive delta,
      it examines a set of candidates. If it prefers a candidate with a larger
      delta, it considers the effect of flipping candidates with smaller delta
      to their opposite bound (thus maintaining dual feasibility).

  With the germ of the idea in hand, some poking around on the web turned up
  a paper by Istvan Maros, ``A Generalized Dual Phase-2 Simplex Algorithm'',
  Departmental Technical Report 2001/2, Department of Computing, Imperial
  College, London, January 2001, revised December 2002, ISSN 1469-4174. It
  describes another variation on this idea, and describes some very
  preliminary results, ending with a comment that this needs to be explored
  further.  This material is repeated in Chap. 10 of Maros' book,
  "Computational Techniques of the Simplex Method".

      The general outline of generalised dual pivoting is as follows:

	scan for entry candidates x<j> and sort by delta<j>

	foreach sorted x<j>

	  calculate abar<j> = inv(B)a<j> ;
	  calculate the change in infeasibility if we pivot on x<j>
	   = (change due to previous flips) + (change moving delta<j>) ;
	  calculate the change in infeasibility if we flip x<j>
	   = (change due to previous flips) + (change flipping x<j>) ;

	  if we can't flip x<j>, quit --- it, or a predecessor, is the pivot
	  otherwise, remember the infeasibility if we choose
	  x<j>, and continue to the next candidate
	end foreach

	pivot on the x<j> with the best infeasibility, flipping variables as
	required.

  The foreach loop ends when either 1) we reach a variable that can't be
  flipped, hence must be pivoted, or 2) flipping a variable would drive
  x<i> within bound. (Remember, the value of x<i> is the dual reduced cost.
  When it comes within bound, the reduced cost changes from favourable to
  unfavourable.) The pivot can be the variable that ended the scan, or any
  of its predecessors (modulo numerical stability considerations).

  The development above implies we're tracking the primal infeasibilty of all
  basic variables. And indeed the original thought was to choose the incoming
  primal variable to minimise the maximum infeasibility over the primal basic
  variables.  It turns out that always choosing a pivot this way is not a
  good strategy --- the DYLP tech report says more about this --- but it
  looks to be worth some further experimentation. The default strategy is to
  look solely at the infeasibility of x<i> (the leaving variable) unless
  swings in other primal variables get out of hand. Put another way,
  improving the dual objective is the goal unless swings in primal variable
  values become large enough to threaten numerical stability.

  This is enhanced with a couple of tweaks. The most important of these is the
  notion of promoting a sane (numerically stable) pivot. It'll often happen
  that the front of the candidate list is occupied with variables x<f> with
  delta<f> = 0 (or some tiny value) and abar<if> numerically unstable. But
  precisely because abar<if> is tiny, it might be that we can pass over it to
  another x<q> with a non-zero delta<q> and a much larger abar<iq>. As long
  as (delta<q>/abar<iq>)*abar<if> is sufficiently small that we can tolerate
  the dual infeasibility of cbar<f>, we can take the better pivot on x<q>.
  We avoid a degenerate pivot and we have better numerical stability.
*/

/*
  Structure used to carry around candidates for entry.

  Field		Description
  -----		-----------
  ndx		index of the variable x<k>
  abarik	|abar<ik>|
  ratioik	stability ratio of abar<ik>; larger is better, 1.0 is the
		boundary between stable and unstable
  madpiv	TRUE if abar<ik> failed the pivot stability check
  ddelta	magnitude of dual delta |delta<k>| = |cbar<k>/abar<ik>|
  pivdir	direction of movement if x<k> is the entering primal
		variable; 1 for rising, -1 for falling
  flippable	TRUE if x<k> has an opposite bound (hence can be flipped).
  rev		TRUE if x<k> is actually a reverse pivot (i.e., cbar is
		slightly on the wrong side of 0, and this will bring us
		back to dual feasibility).
  
  piv		If x<k> is used as the entering variable
    delta 	  standard primal delta<k> = delta<i>/abar<ik> when x<k> is
		  the entering variable, where delta<i> is the movement
		  required to drive x<i> out of the basis.
    inf		  change in total primal infeasibility if x<k> enters
    maxinf	  maximum infeasibility of a primal variable after x<k>
		  is pivoted into the basis
  
  flip		If x<k> is flipped to its opposite bound
    delta	  primal delta for movement to opposite bound
    inf		  change in total primal infeasibility if x<k> is flipped
		  to the opposite bound
    maxinf	  maximum infeasibility of a primal variable after x<k> is
		  flipped
*/

typedef struct { int ndx ;
		 double abarik ;
		 double ratioik ;
		 bool madpiv ;
		 double ddelta ;
		 int pivdir ;
		 bool flippable ;
		 bool rev ;
		 struct { double delta ;
			  double inf ;
			  double maxinf ; } piv ;
		 struct { double delta ;
			  double inf ;
			  double maxinf ; } flip ; } dualcand_struct ;

/*
  # ifdef DYLP_PARANOIA
    static int predictiter ;
    static double predicttotinf,predictmaxinf ;
  # endif
*/


static int dualcand_cmp (const void *p_dualcand1, const void *p_dualcand2)
/*
  Comparison function for qsort to sort an array of dualcand_structs. See
  scanForDualInCands for the sort criteria. The primary criterion is dual
  delta; it only gets interesting when the deltas are equal.

  It turns out that sorting by nonincreasing ratio<ik> as the tiebreaker is
  critical. Otherwise, we're trying to bring x<i> to bound using dinky little
  pivots, and in the meantime we're swinging the values of the other basic
  variables off into the ozone.


  Parameters:
    p_dualcand1,2: dualcand_structs to be compared

  Returns: -1 if dualcand1 < dualcand2
	    0 for equality
	    1 if dualcand1 > dualcand2
*/
{ double delta1,delta2,ratio1,ratio2 ;
  bool flip1,flip2,mad1,mad2 ;
  const dualcand_struct *dualcand1,*dualcand2 ;

  dualcand1 = (const dualcand_struct *) p_dualcand1 ;
  dualcand2 = (const dualcand_struct *) p_dualcand2 ;

  delta1 = dualcand1->ddelta ;
  delta2 = dualcand2->ddelta ;
  mad1 = dualcand1->madpiv ;
  mad2 = dualcand2->madpiv ;
/*
  The primary criterion is the size of the delta. See promoteSanePivot for
  an optimization. Unfortunately, we'll have sort problems if we do it here.
*/
  if (delta1 < delta2)
  { return (-1) ; }
  else
  if (delta1 > delta2)
  { return (1) ; }
/*
  Prefer non-mad pivots. These need to be pushed up here lest they be trapped
  behind tiny pivots with equal delta.
*/
  if (mad1 != mad2)
  { if (mad1 == TRUE)
    { return (1) ; }
    else
    { return (-1) ; } }
/*
  Now we have equal deltas. For degenerate candidates (delta == 0) prefer
  flippable; for nondegenerate, prefer nonflippable.
*/
  flip1 = dualcand1->flippable ;
  flip2 = dualcand2->flippable ;
  if (delta1 == 0)
  { if (flip1 != flip2)
    { if (flip1 == TRUE)
	return (-1) ;
      else
	return (1) ; } }
  else
  { if (flip1 != flip2)
    { if (flip1 == FALSE)
	return (-1) ;
      else
	return (1) ; } }
/*
  Equal deltas, and both flippable or nonflippable. Order by nonincreasing
  pivot ratio.
*/
  ratio1 = dualcand1->ratioik ;
  ratio2 = dualcand2->ratioik ;
  if (ratio1 > ratio2)
    return (-1) ;
  else
  if (ratio1 < ratio2)
    return (1) ;
  else
    return (0) ; }

static void promoteSanePivot (dualcand_struct *incands)
/*
  This routine attempts to promote a sane pivot high enough in the list to be
  used. The rationale is this: Mad pivots are tiny. If the delta of a sane
  pivot is small enough, we won't cause dual infeasibility. Suppose we have
  two variables x<s> (sane) and x<m> (mad). If 
    (delta<s>-delta<m>)*a<i,m> < tol
  then it's a good bet we can promote x<s> past x<m> in the pivot list.

  So --- our goal is to walk along incands until we find the first sane
  pivot. If everything up to that point has been flippable, we need do
  nothing.  But if we've passed nonflippable mad pivots, try to promote the
  sane one.

  Parameters:
    incands:	(i) the candidate list;
		    incands[0].ndx is the number of candidates
		(o) the candidate list, possibly rearranged as described

  Returns: undefined
*/

{ int ndx,candcnt,firstnoflip,firstsane ;
  double tol ;
  dualcand_struct sane,insane ;

/*
  Do the initial scan. If we find a sane pivot and havn't run into a
  non-flippable, non-reverse pivot, all is copacetic and we can return
  with no further effort. There may be no sane pivots.
*/
  candcnt = incands[0].ndx ;
  firstnoflip = -1 ;
  firstsane = -1 ;
  for (ndx = 1 ; ndx <= candcnt ; ndx++)
  { if (incands[ndx].madpiv == FALSE)
    { firstsane = ndx ;
      break ; }
    if (firstnoflip < 0 &&
	incands[ndx].flippable == FALSE && incands[ndx].rev == FALSE)
    { firstnoflip = ndx ; } }

  if (firstsane < 0 || (firstnoflip < 0 && firstsane > 0)) return ;
/*
  We have a sane pivot at firstsane and we have mad pivots from
  firstnoflip to firstsane-1. At least the pivot at firstnoflip cannot be
  flipped.

  The tolerance is calculated to split the difference in exponent between the
  zero tolerance and the feasibility tolerance. Reverse pivots can always be
  stepped over.
*/
  tol = log10(dy_tols->dfeas/dy_tols->cost)/2 ;
  tol = dy_tols->cost*pow(10.0,tol) ;
  for (ndx = firstsane ; ndx > firstnoflip ; ndx--)
  { sane = incands[ndx] ;
    insane = incands[ndx-1] ;
    if (insane.rev == TRUE ||
	(sane.ddelta-insane.ddelta)*insane.abarik < tol)
    { incands[ndx-1] = sane ;
      incands[ndx] = insane ; }
    else
    { break ; } }
/*
  That's it. With any luck, we've promoted a sane candidate to a place where it
  can be used.
*/
# ifndef DYLP_NDEBUG
  if ((dy_opts->print.pivoting >= 2 && ndx == firstnoflip) ||
      dy_opts->print.pivoting >= 3)
  { int j ;
    j = incands[firstnoflip].ndx ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t      first no-flip %s (%d) at %d, ",
	        consys_nme(dy_sys,'v',j,FALSE,NULL),j,firstnoflip) ;
    j = incands[firstsane].ndx ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"first sane %s (%d) at %d",
	        consys_nme(dy_sys,'v',j,FALSE,NULL),j,firstsane) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,", promoted to %d.  ",ndx) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL)
  { if (firstsane > 1 && ndx == 1) dy_stats->dmulti.promote++ ; }
# endif

  return ; }



static dyret_enum scanForDualInCands (dualcand_struct *incands, int outdir,
				      double *abari, double maxabari)

/*
  This routine scans the nonbasic variables for candidates to become the
  entering variable; the criteria are outlined below at the head of the
  scan loop. Candidates are stored in the incands array. For each candidate,
  we record the index, magnitude of the pivot and dual delta, whether the
  pivot is unstable, and whether or not the variable is flippable.
 
  Once we've completed the scan, multiple candidates are sorted. The primary
  sort order is nondecreasing dual delta, delta<k>.
    * Within the degenerate pivot group, delta<k> == 0, we can flip variables
      at will, so the secondary sort orders are flippable < nonflippable, then
      nonincreasing ratio<ik>. Paraphrased, flip all we can, and pivot on the
      variable with the most stable pivot. But hope that all degenerate
      variables are flippable and we can move on to ...
    * The nondegenerate group, where we need to be more circumspect. To flip
      a variable, we need to drive the reduced cost past zero to the opposite
      sign. So, to flip a variable with a given delta, we actually need to
      pivot on a variable with a greater delta. Sort nonflippable <
      flippable, then nonincreasing ratio<ik>.  Then we don't worry about
      flipping when we're going to have to pivot, and we'll get the most
      stable pivot.

  Parameters:
    incands:	(i) pointer to array to hold the candidates for entry
		(o) filled and sorted array; in particular
		    incands[0].ndx is the number of candidates
		    incands[1..] are the candidates
    outdir:	direction of motion of the leaving variable x<i>
    abari:	The transformed pivot row, row i of inv(B)N
    maxabari:	The maximum value in the pivot row abar<i>

  Returns: dyrOK if there are candidates and no errors
	   dyrUNBOUND if there are no candidates
	   dyrFATAL on error (only when we're paranoid)
*/

{ int n,k,reject,candcnt,dirk ;
  double abarik,cbark,deltak,ratioik ;
  double *vub,*vlb ;
  flags statk ;
  bool rev ;
  dyret_enum retval ;

# ifdef DYLP_PARANOIA
  const char *rtnnme = "scanForDualInCands" ;
# endif

# if !defined(DYLP_NDEBUG) || defined(DYLP_PARANOIA)
  int print ;

  print = dy_opts->print.pivoting ;	/* suppress print in dy_chkpiv */
  dy_opts->print.pivoting = 0 ;
# endif

# ifdef DYLP_PARANOIA
  if (incands == NULL)
  { errmsg(2,rtnnme,"incands array") ;
    return (dyrFATAL) ; }
# endif

  n = dy_sys->varcnt ;
  vub = dy_sys->vub ;
  vlb = dy_sys->vlb ;
  incands[0].ndx = -1 ;

# ifndef DYLP_NDEBUG
  if (print >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    gathering candidates to enter ... ") ;
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tVariable\t  cbar<k>\tabar<i,k>\t  delta\tDisp") ; } }
# endif
/*
  Open a loop and scan for suitable nonbasic variables.

  Note that we shouldn't see superbasics here, as we're primal infeasible
  while running dual simplex and superbasics should not be created. Nor
  should we see nonbasic free variables; they can't be nonbasic in a dual
  feasible solution. dy_chkstatus enforces this when we're paranoid.
*/
  candcnt = 0 ;
  for (k = 1 ; k <= n ; k++)
  { if (dy_lp->degen > 0 && dy_ddegenset[k] != dy_lp->degen) continue ;
#   ifdef DYLP_PARANOIA
    if (dy_chkstatus(k) == FALSE)
    { incands[0].ndx = -1 ;
      dy_opts->print.pivoting = print ;
      return (dyrFATAL) ; }
#   endif
/*
  Test if x<k> is a candidate to enter:
    * exclude basic and NBFX variables,
    * exclude if abar<ik> == 0,
    * exclude if the sign of abar<ik> is wrong given the direction
      of motion for x<i> and x<k>,
    * flag x<k> if abar<ik> is too small to be stable.

  Case analysis for cbar<i> = -cbar<k>/abar<i,k> yields:
    * x<i> rising to lb<i> and leaving, cbar<i> to be >= 0
      + x<k> nonbasic at l<k>, c<k> >= 0, implies abar<i,k> must be <= 0
      + x<k> nonbasic at u<k>, c<k> <= 0, implies abar<i,k> must be >= 0
    * x<i> dropping to ub<i> and leaving, cbar<i> to be <= 0
      + x<k> nonbasic at l<k>, c<k> >= 0, implies abar<i,k> must be >= 0
      + x<k> nonbasic at u<k>, c<k> <= 0, implies abar<i,k> must be <= 0

  Some fancy footwork for dealing with variables where cbar<k> has the wrong
  sign, but is within the dual feasibility tolerance: Case analysis says that
  we can deal with this by lying about status<k> and reversing the normal
  direction of motion for the entering primal. (For example, if x<j> is NBLB
  but -dfeas < cbar<j> < 0, claim x<j> is NBUB for purposes of evaluating it,
  and note that when x<j> enters it will decrease, becoming BLLB).

  That won't quite get us out of trouble, though. Suppose we have a variable
  x<k> that's NBLB with cbar<j> < 0 and abar<i,k> < 0 (i.e., the pivot is the
  wrong sign for the `slightly infeasible' case). We can safely ignore the
  case where the pivot is the wrong sign when a candidate variable is dual
  feasible. But if it's slightly infeasible, we can't be so cavalier. The
  risk is that we'll decide to do a nondegenerate pivot with a candidate x<q>
  where abar<i,q> < 0 (a standard NBLB candidate, or a slightly infeasible
  NBUB candidate). Then we'll drive cbar<k> further into infeasibility, and
  may well lose dual feasibility. The solution is to declare x<k> a candidate
  for a normal degenerate pivot by forcing cbark to 0 below.
*/
    statk = dy_status[k] ;
    cbark = dy_cbar[k] ;
    if ((flgon(statk,vstatNBUB) && cbark > 0) ||
	(flgon(statk,vstatNBLB) && cbark < 0))
    { comflg(statk,vstatNBUB|vstatNBLB) ;
      rev = TRUE ; }
    else
    { rev = FALSE ; }
    reject = 0 ;
    abarik = abari[k] ;
    if (flgon(statk,vstatBASIC|vstatNBFX))
    { reject = -1 ; }
    else
    if (withintol(abarik,0.0,dy_tols->zero))
    { reject = -2 ; }
    else
    { if (outdir == -1)
      { if ((flgon(statk,vstatNBUB) && abarik > 0) ||
	    (flgon(statk,vstatNBLB) && abarik < 0))
	  { reject = -3 ; } }
      else
      { if ((flgon(statk,vstatNBUB) && abarik < 0) ||
	    (flgon(statk,vstatNBLB) && abarik > 0))
	  { reject = -3 ; } } }
    if (rev == TRUE && reject == -3)
    { cbark = 0 ;
      rev = FALSE ;
      comflg(statk,vstatNBUB|vstatNBLB) ;
      reject = 0 ; }
    if (reject >= 0)
    { ratioik = dy_chkpiv(abarik,maxabari) ;
      if (ratioik < 1.0)
      { reject = 1 ; } }
    else
    { ratioik = quiet_nan(0) ; }
    
#   ifndef DYLP_NDEBUG
    if (print >= 5 || (print >= 4 && reject >= 0))
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%-8s (%d)",
		  consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t%8g",cbark,abarik) ;
      if (print >= 5)
      { switch (reject)
	{ case -1:
	  { dyio_outfmt(dy_logchn,dy_gtxecho,"\t\trejected -- status %s",
		        dy_prtvstat(statk)) ;
	    break ; }
	  case -2:
	  { dyio_outfmt(dy_logchn,dy_gtxecho,"\t\trejected -- zero pivot") ;
	    break ; }
	  case -3:
	  { dyio_outfmt(dy_logchn,dy_gtxecho,
			"\t\trejected -- wrong sign abar<ik>") ;
	    break ; } } } }
#   endif

    if (reject < 0) continue ;
/*
  x<k> is a possible candidate for entry. Record it in incands with its
  associated dual delta. For reverse pivots, negate the deltak to make
  sure it ends up at the very front of the list.
*/
    incands[++candcnt].ndx = k ;
    incands[candcnt].abarik = fabs(abarik) ;
    incands[candcnt].ratioik = ratioik ;
    if (reject == 1)
    { incands[candcnt].madpiv = TRUE ; }
    else
    { incands[candcnt].madpiv = FALSE ; }
    deltak = fabs(cbark/abarik) ;
    if (rev == TRUE)
    { incands[candcnt].rev = TRUE ;
      incands[candcnt].ddelta = -deltak ;
      comflg(statk,vstatNBLB|vstatNBUB) ; }
    else
    { incands[candcnt].rev = FALSE ;
      incands[candcnt].ddelta = deltak ; }
    if (outdir == -1)
    { if (abarik > 0)
      { dirk = 1 ; }
      else
      { dirk = -1 ; } }
    else
    { if (abarik > 0)
      { dirk = -1 ; }
      else
      { dirk = 1 ; } }
    if (rev == FALSE)
    { incands[candcnt].pivdir = dirk ; }
    else
    { incands[candcnt].pivdir = -dirk ; }
    if (flgon(statk,vstatNBUB))
    { if (vlb[k] > -dy_tols->inf)
      { incands[candcnt].flippable = TRUE ;
	incands[candcnt].flip.delta = vlb[k]-vub[k] ; }
      else
      { incands[candcnt].flippable = FALSE ; } }
    else
    { if (vub[k] < dy_tols->inf)
      { incands[candcnt].flippable = TRUE ;
	incands[candcnt].flip.delta = vub[k]-vlb[k] ; }
      else
      { incands[candcnt].flippable = FALSE ; } }

#   ifndef DYLP_NDEBUG
    if (print >= 4)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t accepted",deltak) ;
      if (reject == 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (mad)") ; }
      if (deltak == 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (degen)") ; }
      if (incands[candcnt].flippable == TRUE)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (flip)") ; }
      if (rev == TRUE)
      { dyio_outfmt(dy_logchn,dy_gtxecho," (rev)") ; } }
#   endif
  }
/*
  If we have a list of candidates, sort the list. If we have no candidates,
  we're unbounded.

  The assignments to incands[0].[ddelta,madpiv] are targetted at a quirk of
  qsort: it seems to occasionally sort incands[0]. An optimization of some
  sort, presumably.
*/
  incands[0].ndx = candcnt ;
  incands[0].ddelta = -dy_tols->inf ;
  incands[0].madpiv = FALSE ;
  if (candcnt > 0)
  { if (candcnt > 1)
    { qsort(&incands[1],candcnt,sizeof(dualcand_struct),dualcand_cmp) ; }
    retval = dyrOK ; }
  else
  { retval = dyrUNBOUND ; }
/*
  Now try a little bit of an optimization, to push decent pivots up past
  mad pivots.
*/
# ifndef DYLP_NDEBUG
  dy_opts->print.pivoting = print ;
# endif

  promoteSanePivot(incands) ;

# ifndef DYLP_NDEBUG
  if (print >= 1)
  { if (print >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n   ") ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tVariable\tratio<ik>\tdelta<k>") ;
      for (n = 1 ; n <= candcnt ; n++)
      { k = incands[n].ndx ;
	ratioik = incands[n].ratioik ;
	deltak = incands[n].ddelta ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%-8s (%d)",
		    consys_nme(dy_sys,'v',k,FALSE,NULL),k) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\t%8g\t%8g",ratioik,deltak) ;
	if (incands[n].madpiv == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (mad)") ; }
	if (deltak == 0)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (degen)") ; }
	if (incands[n].flippable == TRUE)
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (flip)") ; }
	statk = dy_status[k] ;
	if ((flgon(statk,vstatNBUB) && incands[n].pivdir == 1) ||
	    (flgon(statk,vstatNBLB) && incands[n].pivdir == -1))
	{ dyio_outfmt(dy_logchn,dy_gtxecho," (rev)") ; }
	if (n > 1)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      " (%g)",deltak-incands[n-1].ddelta) ; } }
      dyio_outfmt(dy_logchn,dy_gtxecho,"\n   ") ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,"%d candidates.",candcnt) ;
    if (retval != dyrOK)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  " Returning %s.",dy_prtdyret(retval)) ; } }
  dy_opts->print.pivoting = print ;
# endif

  return (retval) ; }



static int calcInfChange (dualcand_struct *candk, int i, double *xbasic)

/*
  This routine calculates the change in primal infeasibility due to the change
  in the value of x<k> under two assumptions:
    * x<k> is the pivot variable
    * x<k> will be flipped to its opposite bound

  Parameters:
    candk:	(i) dualcand structure for x<k>
		(o) information for x<k> pivoted or flipped is set
    i:		index of x<i>, the leaving variable
    xbasic:	(i/o) current values of basic variables; updated with each
		pricing to reflect the effect of flipping the candidate.

  Returns:  1 if the calculation completes without error, and x<i> is still
	      infeasible
	    0 if the calculation completes without error, and the flip makes
	      x<i> feasible
	   -1 if there's an error.
*/

{ int m,ipos ;
  flags stati ;
  double *vlb,*vub ;
  double xi,lbi,ubi,deltai ;

  int k ;
  double newxk,lbk,ubk,pivdeltak,flipdeltak,
	 pivinf,maxpivinf,flipinf,maxflipinf ;
  double *abark,abarik ;
  bool flippable ;
  flags statk ;

  int l,lpos ;
  double xl,newxl,lbl,ubl,abarlk,curinf,infl ;

  const char *rtnnme = "calcInfChange" ;

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->dmulti.evals++ ;
# endif

  m = dy_sys->concnt ;
  vub = dy_sys->vub ;
  vlb = dy_sys->vlb ;
/*
  Get x<i>, its status and upper and lower bounds, and calculate delta<i>.
*/
  ipos = dy_var2basis[i] ;
  xi = xbasic[ipos] ;
  stati = dy_status[i] ;
  ubi = vub[i] ;
  lbi = vlb[i] ;
  if (flgon(stati,vstatBUUB))
  { deltai = ubi-xi ; }
  else
  { deltai = lbi-xi ; }
/*
  Get k, its status, and upper and lower bounds.
*/
  k = candk->ndx ;
  statk = dy_status[k] ;
  lbk = vlb[k] ;
  ubk = vub[k] ;
  if (candk->flippable == TRUE)
  { flippable = TRUE ;
    flipdeltak = candk->flip.delta ; }
  else
  { flippable = FALSE ;
    flipdeltak = quiet_nan(0) ; }
/*
  Retrieve a<k>, calculate abar<k>, then calculate the delta<k> for the
  case where x<k> is chosen as the entering variable.
*/
  abark = NULL ;
  if (consys_getcol_ex(dy_sys,k,&abark) == FALSE)
  { errmsg(122,rtnnme,dy_sys->nme,
	   "column",consys_nme(dy_sys,'v',k,TRUE,NULL),k) ;
    if (abark != NULL) FREE(abark) ;
    return (-1) ; }
  dy_ftran(abark,FALSE) ;
  abarik = abark[ipos] ;
  pivdeltak = -deltai/abarik ;
  setcleanzero(pivdeltak,dy_tols->zero) ;
  candk->piv.delta = pivdeltak ;
  pivinf = 0 ;
  if (pivdeltak < 0)
  { /* candk->pivdir = -1 ; */
    newxk = ubk+pivdeltak ;
    if (belowbnd(newxk,lbk))
    { pivinf = lbk-newxk ; } }
  else
  { /* candk->pivdir = 1 ; */
    newxk = lbk+pivdeltak ;
    if (abovebnd(newxk,ubk))
    { pivinf = newxk-ubk ; } }
  maxpivinf = pivinf ;
/*
  Now it's just a straightforward grind.  For each basic variable x<l> s.t.
  abar<lk> != 0, calculate the affect on primal feasibility for x<k> pivoted
  and flipped. abark should be groomed, so we can test for abar<lk> == 0.
  Update xbasic when we calculate the effect of the flip.

  Note that we need to consider unaffected x<l> when calculating the max
  values --- they will hold for both pivot and flip.
*/
  curinf = 0 ;
  flipinf = 0 ;
  maxflipinf = 0 ;
  for (lpos = 1 ; lpos <= m ; lpos++)
  { l = dy_basis[lpos] ;
    xl = xbasic[lpos] ;
    lbl = vlb[l] ;
    ubl = vub[l] ;
    if (abovebnd(xl,ubl))
    { infl = xl-ubl ; }
    else
    if (belowbnd(xl,lbl))
    { infl = lbl-xl ; }
    else
    { infl = 0 ; }
    abarlk = abark[lpos] ;
    if (abarlk == 0)
    { if (infl > maxpivinf) maxpivinf = infl ;
      if (infl > maxflipinf) maxflipinf = infl ; 
      continue ; }
    curinf += infl ;

    newxl = xl-abarlk*pivdeltak ;
    if (abovebnd(newxl,ubl))
    { infl = newxl-ubl ; }
    else
    if (belowbnd(newxl,lbl))
    { infl = lbl-newxl ; }
    else
    { infl = 0 ; }
    pivinf += infl ;
    if (infl > maxpivinf) maxpivinf = infl ;

    if (flippable == TRUE)
    { newxl = xl-abarlk*flipdeltak ;
      if (abovebnd(newxl,ubl))
      { infl = newxl-ubl ; }
      else
      if (belowbnd(newxl,lbl))
      { infl = lbl-newxl ; }
      else
      { infl = 0 ; }
      xbasic[lpos] = newxl ;
      flipinf += infl ;
      if (infl > maxflipinf) maxflipinf = infl ; } }
  
  candk->piv.inf = -curinf+pivinf ;
  candk->piv.maxinf = maxpivinf ;
  if (flippable == TRUE)
  { candk->flip.inf = -curinf+flipinf ;
    candk->flip.maxinf = maxflipinf ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t%s (%d) piv = %g, pivmax = %g",
	        consys_nme(dy_sys,'v',k,FALSE,NULL),k,
	        candk->piv.inf,candk->piv.maxinf) ;
    if (flippable == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho," flip = %g, flipmax = %g",
		  candk->flip.inf,candk->flip.maxinf) ; }
    dyio_outfmt(dy_logchn,dy_gtxecho," pivdelta = %g",pivdeltak) ;
    if (flippable == TRUE)
    { dyio_outfmt(dy_logchn,dy_gtxecho," flipdelta = %g",flipdeltak) ; } }
# endif
/*
  That's it. The final thing we need to do is check on the status of the
  entering variable x<i>. If this flip would change it, we have to disallow
  the flip.
*/
  FREE(abark) ;

  xi = xbasic[ipos] ;
  if ((flgon(stati,vstatBUUB) && !abovebnd(xi,ubi)) ||
      (flgon(stati,vstatBLLB) && !belowbnd(xi,lbi)))
  { return (0) ; }
  else
  { return (1) ; } }



bool selectWithInf (int i, dualcand_struct *incands,
		    int *indices, double *candinf, double *startinf)

/*
  This routine scans the list of candidates x<k>, calculating the primal
  infeasibility over all basic variables if x<k> is flipped and if it's
  chosen as the pivot. This is an expensive calculation, because we'll need
  to calculate inv(B)a<k> for each variable that we price.

  The candidates in incands are assumed to be sorted (see the comments with
  scanForDualInCands), and the calculation of infeasibility is based on the
  notion that all variables preceding the variable under evaluation will be
  flipped.

  In the end, selectWithInf will return up to three indices:
    (1) An index marking the end of a sequence of flips, with no final pivot,
	chosen so as to minimise the maximum primal infeasibility.
    (2) An index marking the end of a sequence of flips, with a final pivot,
	chose so as to minimise the maximum primal infeasibility.
    (3) An index marking the end of the longest possible sequence of flips,
	with a final pivot.
  Any or all of the above can be -1, indicating that there was no qualified
  candidate with the required qualities.

  Parameters:
    i:		The index of x<i>, the entering variable.
    incands:	The list of candidates, indexed from 1. incands[0].ndx is the
		number of candidates.
    indices:	(i) An array with three entries.
		(o) The index for candidate k is in indices[k-1].
    candinf:	(i) An array with three entries.
		(o) The infeasibilities corresponding to the candidates
		    specified in indices.
    startinf	(i) An array with two entries.
		(o) [0] is total infeasibility
		    [1] is maximum infeasibility

  Returns: TRUE if the selection process proceeds without error, FALSE
	   otherwise.
*/

{ int m,j,jpos ;
  double xj,lbj,ubj,infj ;
  double *vlb,*vub ;

  double starttotinf,startmaxinf ;
  double *xbasic ;

  int ndx,candcnt,price_retval ;
  double flipinfk,pivinfk,totinfk ;
  int bestflipcand,bestpivcand,lastpivcand ;
  double bestflipinf,bestpivinf,lastpivinf ;
  bool pivEndsScan,flipEndsScan ;
  dualcand_struct *candk ;

  const char *rtnnme = "selectWithInf" ;

# ifndef DYLP_NDEBUG
  int lastdegen = 0 ;
# endif

  m = dy_sys->concnt ;
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
/*
  We have multiple candidates and at least a chance to flip variables.  To
  properly calculate the change in infeasibility, we'll need to know the
  infeasibility now, and we'll need a copy of xbasic that we can modify as we
  go.
*/
  
  xbasic = (double *) MALLOC((m+1)*sizeof(double)) ;
  starttotinf = 0 ;
  startmaxinf = 0 ;
  for (jpos = 1 ; jpos <= m ; jpos++)
  { xbasic[jpos] = dy_xbasic[jpos] ;
    xj = xbasic[jpos] ;
    j = dy_basis[jpos] ;
    lbj = vlb[j] ;
    ubj = vub[j] ;
    infj = 0 ;
    if (belowbnd(xj,lbj))
    { infj = lbj-xj ; }
    else
    if (abovebnd(xj,ubj))
    { infj = xj-ubj ; }
    if (infj > startmaxinf) startmaxinf = infj ;
    starttotinf += infj ; }
  startinf[0] = starttotinf ;
  startinf[1] = startmaxinf ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n      starting inf tot = %g, max = %g",
	        starttotinf,startmaxinf) ; }
# endif
/*
  # ifdef DYLP_PARANOIA
    if (dy_lp->d2.iters > 0 && predictiter+1 == dy_lp->d2.iters)
    { if (!atbnd(predicttotinf,starttotinf))
      { warn(350,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "total",predicttotinf,predictiter,
	     starttotinf,predicttotinf-starttotinf) ; }
      if (!atbnd(predictmaxinf,startmaxinf))
      { warn(350,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "maximum",predictmaxinf,predictiter,
	     startmaxinf,predictmaxinf-startmaxinf) ; } }
  # endif
*/

/*
  For better or worse, we'll stick with the sort order while deciding what to
  do.  calcInfChange will price a candidate x<k>, determining the change in
  primal infeasibility when we flip the variable and when we choose it as the
  pivot. This is expensive, because we have to calculate inv(B)a<k> to do the
  calculation.

  The predicted maximum primal infeasibility when we use a variable x<k> as a
  pivot (piv.maxinf), includes the change due to flipping all candidates prior
  to x<k> in the candidate list.

  The first nonflippable variable stops the scan --- we have to use it to
  pivot. Within the set of degenerate candidates, we can choose to do only
  flips. Once we get to nondegenerate candidates, we need to do a pivot with
  a dual delta greater than previous candidates in order to flip the previous
  candidates.

  We're perfectly happy to flip variables with mad pivots, but we won't use
  them as pivot candidates.
*/
  bestpivcand = -1 ;
  bestflipcand = -1 ;
  lastpivcand = -1 ;
  bestpivinf = dy_tols->inf ;
  bestflipinf = dy_tols->inf ;
  lastpivinf = quiet_nan(0) ;
  pivEndsScan = FALSE ;
  flipEndsScan = FALSE ;
  candcnt = incands[0].ndx ;
  totinfk = 0 ;
/*
  Scan the degenerate pivots. If we run across a nonflippable, nonreverse
  variable, we've found our entering variable. It's also possible that the
  cumulative effect of flips up to and including cand<k> will make x<i>
  feasible or drive it right past the opposite bound (price_retval == 0). In
  that case, we stop with cand<k>, and can choose to flip it or pivot on it.
  Note that for reverse pivots, the primal variable will not change bound
  (we're recovering from slight dual infeasibility).
*/
  for (ndx = 1,candk = &incands[1] ;
       candk->ddelta <= 0 && ndx <= candcnt ;
       ndx++,candk++)
  { if (candk->rev == FALSE)
    { price_retval = calcInfChange(candk,i,xbasic) ; }
    else
    { price_retval = 1 ; }
    if (price_retval < 0)
    { FREE(xbasic) ;
      errmsg(348,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters+1,
	     consys_nme(dy_sys,'v',candk->ndx,FALSE,NULL),candk->ndx) ;
      return (FALSE) ; }
    flipinfk = candk->flip.maxinf ;
    pivinfk = candk->piv.maxinf ;
    totinfk += candk->flip.inf ;

# if 0
  predicttotinf = 0 ;
  for (jpos = 1 ; jpos <= m ; jpos++)
  { xj = xbasic[jpos] ;
    j = dy_basis[jpos] ;
    lbj = vlb[j] ;
    ubj = vub[j] ;
    if (belowbnd(xj,lbj))
    { predicttotinf += lbj-xj ; }
    else
    if (abovebnd(xj,ubj))
    { predicttotinf += xj-ubj ; } }
  if (dy_opts->print.pivoting >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	     "\n      debug piv = %g, flip = %g, tot = %g, infeasibility = %g",
	     pivinfk,flipinfk,totinfk,predicttotinf) ; }
# endif

    if (candk->madpiv == FALSE)
    { lastpivcand = ndx ;
      lastpivinf = pivinfk ;
      if (pivinfk < bestpivinf)
      { bestpivcand = ndx ;
	bestpivinf = pivinfk ; } }
    if (candk->flippable == FALSE && candk->rev == FALSE)
    { pivEndsScan = TRUE ;
      break ; }
    if (flipinfk < bestflipinf)
    { bestflipcand = ndx ;
      bestflipinf = flipinfk ; }
    if (price_retval == 0)
    { flipEndsScan = TRUE ;
      if (candk->madpiv == FALSE)
      { lastpivcand = ndx ;
	lastpivinf = pivinfk ; }
      break ; } }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 2 && ndx > 1)
  { if (pivEndsScan == TRUE || flipEndsScan == TRUE)
    { jpos = ndx ; }
    else
    { jpos = ndx-1 ; }
    lastdegen = jpos ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"\n      after %d degen",jpos) ;
    if (bestflipcand > 0)
    { j = incands[bestflipcand].ndx ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  ", best flip #%d, %s (%d) = %g",bestflipcand,
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j,bestflipinf) ; }
    if (bestpivcand > 0)
    { j = incands[bestpivcand].ndx ;
      dyio_outfmt(dy_logchn,dy_gtxecho,", best piv #%d, %s (%d) = %g",
		  bestpivcand,
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j,bestpivinf) ; }
    if (lastpivcand > 0)
    { j = incands[lastpivcand].ndx ;
      dyio_outfmt(dy_logchn,dy_gtxecho,", last piv #%d, %s (%d) = %g",
		  lastpivcand,
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j,lastpivinf) ; }
    if (bestflipcand < 0 && bestpivcand < 0 && lastpivcand < 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,", nothing") ; }
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
# endif

/*
  If we're not done, continue to scan the nondegenerate pivots. We no longer
  have the choice of just flipping, so we're only searching for a pivot
  candidate. If we get a return value of 0 from calcInfChange, it means we have
  to pivot at that variable.
*/
  if (pivEndsScan == FALSE && flipEndsScan == FALSE)
  { for ( ; ndx <= candcnt ; ndx++,candk++)
    { price_retval = calcInfChange(candk,i,xbasic) ;
      if (price_retval < 0)
      { FREE(xbasic) ;
	errmsg(348,rtnnme,
	       dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters+1,
	       consys_nme(dy_sys,'v',candk->ndx,FALSE,NULL),candk->ndx) ;
	return (FALSE) ; }
      flipinfk = candk->flip.maxinf ;
      pivinfk = candk->piv.maxinf ;
      totinfk += candk->flip.inf ;

# if 0
  bestinf = 0 ;
  for (jpos = 1 ; jpos <= m ; jpos++)
  { xj = xbasic[jpos] ;
    j = dy_basis[jpos] ;
    lbj = vlb[j] ;
    ubj = vub[j] ;
    if (belowbnd(xj,lbj))
    { bestinf += lbj-xj ; }
    else
    if (abovebnd(xj,ubj))
    { bestinf += xj-ubj ; } }
  if (dy_opts->print.pivoting >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	   "\n      debug piv = %g, flip = %g, tot = %g, infeasibility = %g",
	   pivinfk,flipinfk,starttotinf,bestinf) ; }
# endif

      if (candk->madpiv == FALSE)
      { lastpivcand = ndx ;
	lastpivinf = pivinfk ;
	if (pivinfk < bestpivinf)
	{ bestpivcand = ndx ;
	  bestpivinf = pivinfk ; } }
      if (candk->flippable == FALSE || price_retval == 0)
      { pivEndsScan = TRUE ;
	if (candk->madpiv == FALSE)
	{ lastpivcand = ndx ;
	  lastpivinf = pivinfk ; }
	break ; } }

#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivoting >= 2)
    { if (pivEndsScan == TRUE || flipEndsScan == TRUE)
      { jpos = ndx ; }
      else
      { jpos = ndx-1 ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      after %d nondegen",jpos-lastdegen) ;
      if (bestflipcand > 0)
      { j = incands[bestflipcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", best flip #%d, %s (%d) = %g",bestflipcand,
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,bestflipinf) ; }
      if (bestpivcand > 0)
      { j = incands[bestpivcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", best piv #%d, %s (%d) = %g",bestpivcand,
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,bestpivinf) ; }
      if (lastpivcand > 0)
      { j = incands[lastpivcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", last piv #%d, %s (%d) = %g",lastpivcand,
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,lastpivinf) ; }
      if (bestflipcand < 0 && bestpivcand < 0 && lastpivcand < 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,", nothing") ; }
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
  }

  FREE(xbasic) ;

# ifdef DYLP_PARANOIA
  if ((bestpivcand > 0 && lastpivcand < 0) ||
      (bestpivcand < 0 && lastpivcand > 0))
  { errmsg(1,rtnnme,__LINE__) ;
    return (FALSE) ; }
# endif

/*
  Load up the return values and we're done.
*/
  indices[0] = bestflipcand ;
  indices[1] = bestpivcand ;
  indices[2] = lastpivcand ;

  if (bestflipcand > 0)
  { candinf[0] = bestflipinf ; }
  else
  { candinf[0] = quiet_nan(0) ; }
  if (bestpivcand > 0)
  { candinf[1] = bestpivinf ; }
  else
  { candinf[0] = quiet_nan(0) ; }
  if (lastpivcand > 0)
  { candinf[2] = lastpivinf ; }
  else
  { candinf[0] = quiet_nan(0) ; }

  return (TRUE) ; }



bool selectWithoutInf (int i, double *abari, dualcand_struct *incands,
		       int *indices)

/*
  This routine scans the list of candidates x<k>, calculating the primal
  infeasibility of the leaving variable x<i> if x<k> is flipped and if it's
  chosen as the pivot. This is pretty cheap, because we already have the
  pivot row handy.

  The candidates in incands are assumed to be sorted (see the comments with
  scanForDualInCands), and the calculation of infeasibility is based on the
  notion that all variables preceding the variable under evaluation will be
  flipped (or are reverse pivots, correcting slight dual infeasibility).

  In the end, selectWithoutInf will return up to two indices:
    (1) An index marking the end of the longest possible sequence of flips,
	with no final pivot.
    (2) (Reserved, for compatibility with selectWithInf.)
    (3) An index marking the end of the longest possible sequence of flips,
	with a final pivot.
  Any or all of the above can be -1, indicating that there was no qualified
  candidate with the required qualities.

  Parameters:
    i:		The index of x<i>, the entering variable.
    abari:	The pivot row inv(B)N.
    incands:	The list of candidates, indexed from 1. incands[0].ndx is the
		number of candidates.
    indices:	(i) An array with three entries.
		(o) The index for candidate k is in indices[k-1]. indices[1] is
		    always -1.

  Returns: TRUE if the selection process proceeds without error, FALSE
	   otherwise.
*/

{ int m ;
  double *vlb,*vub ;

  int ipos ;
  double xi,lbi,ubi ;
  flags stati ;

  int k ;
  double lbk,ubk,deltak,abarik ;
  flags statk ;

  double startinf ;

  int ndx,candcnt ;
  int lastflipcand,lastpivcand ;
  bool pivEndsScan,flipEndsScan,flipsOnly ;
  dualcand_struct *candk ;

# ifndef DYLP_NDEBUG
  int j,jndx,lastdegen ;
  lastdegen = 0 ;
# endif
/*
  # ifdef DYLP_PARANOIA
    const char *rtnnme = "selectWithOutInf" ;
  # endif
*/

  m = dy_sys->concnt ;
  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
/*
  We have multiple candidates and at least a chance to flip variables.  To
  properly track the change in infeasibility for x<i>, we'll need to know
  the infeasibility now.
*/
  ipos = dy_var2basis[i] ;
  stati = dy_status[i] ;
  xi = dy_x[i] ;
  ubi = vub[i] ;
  lbi = vlb[i] ;
  if (flgon(stati,vstatBUUB))
  { startinf = xi-ubi ; }
  else
  { startinf = lbi-xi ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n      starting inf<%d> = %g",i,startinf) ; }
# endif
/*
  # ifdef DYLP_PARANOIA
    if (dy_lp->d2.iters == 0)
    { predictinf = 0 ;
      predictiter = 0 ; }
    else
    if (dy_lp->d2.iters > 0 && predictiter+1 == dy_lp->d2.iters)
    { if (!atbnd(predictinf,startinf))
      { warn(350,rtnnme,
	     dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	     "x<i>",predictinf,predictiter,startinf,predictinf-startinf) ; } }
  # endif
*/

/*
  For better or worse, we'll stick with the sort order while deciding what to
  do. Given that we're only tracking the infeasibility of the leaving variable
  x<i>, we walk the candidates until a sequence of flips makes x<i> feasible,
  or until we come across a variable we have to use as a pivot.

  We're perfectly happy to flip variables with mad pivots, but we won't use
  them as pivot candidates.
*/
  lastflipcand = -1 ;
  lastpivcand = -1 ;
  pivEndsScan = FALSE ;
  flipEndsScan = FALSE ;
  candcnt = incands[0].ndx ;
/*
  Scan the degenerate and reverse pivots. The scan ends when we run into a
  nonflippable, non-reverse candidate, or when the accumulated effect of
  flips drives x<i> within bounds or right out the other side of its feasible
  region.

  Any x<k> with a stable pivot (including reverse) is a candidate for
  pivoting. The primal variable associated with a reverse pivot will not
  actually flip (remember, we're recovering from slight dual infeasibility).

  We'll consider a sequence of flips until the accumulated change drives x<i>
  into or through feasibility. If x<i> ends up feasible, then it's possible
  to consider a `pivot' sequence that has only flips, with no final pivot.
  Mad pivots are not a problem for flipping --- the flip simply has little
  effect on x<i>.
*/
  for (ndx = 1,candk = &incands[1] ;
       ndx <= candcnt && candk->ddelta <= 0 ;
       ndx++,candk++)
  { k = candk->ndx ;
    if (candk->madpiv == FALSE)
    { lastpivcand = ndx ; }
    if (candk->flippable == FALSE && candk->rev == FALSE)
    { pivEndsScan = TRUE ;
      break ; }
    if (candk->rev == TRUE) continue ;
    statk = dy_status[k] ;
    lbk = vlb[k] ;
    ubk = vub[k] ;
    if (flgon(statk,vstatNBLB))
    { deltak = ubk-lbk ; }
    else
    { deltak = lbk-ubk ; }
    abarik = abari[k] ;
    xi -= abarik*deltak ;
    if ((flgon(stati,vstatBUUB) && !abovebnd(xi,ubi)) ||
	(flgon(stati,vstatBLLB) && !belowbnd(xi,lbi)))
    { if (withinbnds(lbi,xi,ubi))
      { lastflipcand = ndx ; }
      flipEndsScan = TRUE ;
      break ; } }
/*
  If flips alone can drive us to feasibility, then we'll consider using a
  sequence of flips with no final pivot.

  Note that the condition that x<i> be feasible for a flip-only sequence is
  not enough to prevent cycling. Consider two basic variables x<1> and x<2>
  and a nonbasic variable x<3>. I've seen instances where x<1> is selected to
  leave, and a flip of x<3> drives x<1> feasible and x<2> infeasible. x<2> is
  now selected to leave, and x<3> is selected to flip, which drives x<2>
  feasible and x<i> infeasible. Etc., ad infinitum.
*/
  if (flipEndsScan == TRUE && lastflipcand > 0)
  { flipsOnly = TRUE ; }
  else
  { flipsOnly = FALSE ; }

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 2 && ndx >= 1)
  { if (pivEndsScan == TRUE || flipEndsScan == TRUE)
    { jndx = ndx ; }
    else
    { jndx = ndx-1 ; }
    lastdegen = jndx ;
    if (jndx > 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n      after %d degen",jndx) ;
      if (lastflipcand > 0)
      { j = incands[lastflipcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", last flip #%d, %s (%d)",lastflipcand,
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
      if (lastpivcand > 0)
      { j = incands[lastpivcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,", last piv #%d, %s (%d)",
		    lastpivcand,consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
      if (lastflipcand < 0 && lastpivcand < 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,", nothing") ; }
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; } }
# endif

/*
  If we're not done, continue to scan the nondegenerate pivots. We no longer
  have the choice of just flipping, so we're only searching for a pivot
  candidate. But we still need to track reductions in x<i> due to flips.
*/
  if (pivEndsScan == FALSE && flipEndsScan == FALSE)
  { for ( ; ndx <= candcnt ; ndx++,candk++)
    { k = candk->ndx ;
      if (candk->madpiv == FALSE)
      { lastpivcand = ndx ; }
      if (candk->flippable == FALSE)
      { pivEndsScan = TRUE ;
	break ; }
      statk = dy_status[k] ;
      lbk = vlb[k] ;
      ubk = vub[k] ;
      if (flgon(statk,vstatNBLB))
      { deltak = ubk-lbk ; }
      else
      { deltak = lbk-ubk ; }
      abarik = abari[k] ;
      xi -= abarik*deltak ;
      if ((flgon(stati,vstatBUUB) && !abovebnd(xi,ubi)) ||
	  (flgon(stati,vstatBLLB) && !belowbnd(xi,lbi)))
      { flipEndsScan = TRUE ;
	break ; } }

#   ifndef DYLP_NDEBUG
    if (dy_opts->print.pivoting >= 2)
    { if (pivEndsScan == TRUE || flipEndsScan == TRUE)
      { jndx = ndx ; }
      else
      { jndx = ndx-1 ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      after %d nondegen",jndx-lastdegen) ;
      if (lastflipcand > 0)
      { j = incands[lastflipcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    ", last flip #%d, %s (%d)",lastflipcand,
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
      if (lastpivcand > 0)
      { j = incands[lastpivcand].ndx ;
	dyio_outfmt(dy_logchn,dy_gtxecho,", last piv #%d, %s (%d)",lastpivcand,
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
      if (lastflipcand < 0 && lastpivcand < 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,", nothing") ; }
      dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }
#   endif
  }


/*
  Load up the return values and we're done.
*/
  if (flipsOnly == TRUE)
  { indices[0] = lastflipcand ; }
  else
  { indices[0] = -1 ; }
  indices[1] = -1 ;
  indices[2] = lastpivcand ;

  return (TRUE) ; }




dyret_enum dualmultiin (int i, int outdir,
			int *p_xjndx, int *p_indir,
			double *abari, double maxabari, double **p_abarj)

/*
  Implements the selection algorithm described at the head of the file. If
  there are variables to flip, the routine will flip them before returning
  the final pivot variable. On return, p_xjndx and p_indir will be set just as
  for normal dualin.

  Parameters:
    i:		index of the leaving variable x<i>
    outdir:	direction of movement of x<i>,
		   1 if rising to lower bound
		  -1 if falling to upper bound
    p_xjndx:	index of the variable chosen as pivot
    p_indir:	direction of motion for x<j>; 1 for rising, -1 for falling
    abari:	row i of inv(B)N (the pivot row)
    maxabari:	maximum value in abar<i>
    p_abarj:	(o) used to return abar<j>
  
  Returns: dyret_enum code, as follows:
    dyrOK:	the pivot is (dual) nondegenerate
    dyrDEGEN:	the pivot is (dual) degenerate
    dyrUNBOUND:	the problem is dual unbounded (primal infeasible) (i.e.,
		no incoming variable can be selected)
    dyrMADPIV:	if a<ij> fails the numerical stability test
    dyrRESELECT: if the routine elects to do only flips
    dyrFATAL:	fatal confusion
*/

{ int n,m,candcnt ;
  double *vub,*vlb,*accumj ;
  dyret_enum retval,upd_retval,confirm ;
  bool swing,reqchk,flipOnly,sel_retval ;

  int j ;
  double xj,deltaj,starttotinf,startmaxinf ;
  flags statj ;

  int dpstrat,candndx[4] ;
  double candinf[4],startinf[2] ;

  int ndx,bestpivcand,bestflipcand,lastpivcand,bestcand ;
  double bestpivinf,lastpivinf,bestflipinf,bestinf ;

  dualcand_struct *incands,*candk ;

  const char *rtnnme = "dualmultiin" ;

  /* dy_dualpivot.c */
  dyret_enum dy_confirmDualPivot(int i, int j, double *abari,
				 double maxabari, double **p_abarj) ;
# ifndef DYLP_NDEBUG
  double infj ;
# endif

/*
  Setup
*/
  retval = dyrINV ;
  *p_xjndx = 0 ;
  *p_indir = 0 ;
  n = dy_sys->varcnt ;
  m = dy_sys->concnt ;
  vub = dy_sys->vub ;
  vlb = dy_sys->vlb ;
  incands = (dualcand_struct *) MALLOC((n-m+1)*sizeof(dualcand_struct)) ;
/*
  Determine our strategy:
    1: maximum objective improvement
    2: minimum predicted infeasibility
    3: infeasibility reduction if possible, otherwise maximum objective
       improvement.
  If we move to strategy 2 or 3, loosen dfeas accordingly.

  Chances are good that these breakpoints and toobig will change, and/or
  become parameters. Some experimentation is required.
*/
  if (dy_opts->dpsel.flex == TRUE)
  { if (dy_lp->prim.max > dy_tols->toobig/100)
    { dy_opts->dpsel.strat = 2 ; }
    else
    if (dy_lp->prim.max > dy_tols->toobig/1000)
    { dy_opts->dpsel.strat = 3 ; }
    else
    { dy_opts->dpsel.strat = 1 ; } }
  dpstrat = dy_opts->dpsel.strat ;

# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  (%s)%d: selecting entering variable, strategy %d",
	        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters+1,dpstrat) ; }
# endif
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->dmulti.cnt++ ;
# endif

/*
  Generate a sorted list of candidates to enter. No candidates means we're dual
  unbounded (dyrUNBOUND). Fatal error (dyrFATAL) is possible only if we're
  paranoid.
*/
  retval = scanForDualInCands(incands,outdir,abari,maxabari) ;
  if (retval != dyrOK)
  { *p_xjndx = -1 ;
    FREE(incands) ;
    return (retval) ; }
  candk = &incands[0] ;
  candcnt = candk->ndx ;
  candk++ ;
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->dmulti.cands += candcnt ;
# endif
/*
  If we have only one candidate, or the first candidate isn't flippable,
  we have no choice. Confirm that there's no numerical drift, then we can
  return, saving a large amount of work.
*/
  if (candcnt == 1 || (candk->flippable == FALSE && candk->rev == FALSE))
  { *p_xjndx = candk->ndx ;
    *p_indir = candk->pivdir ;
    if (candk->madpiv == TRUE)
    { retval = dyrMADPIV ; }
    else
    if (candk->ddelta == 0)
    { retval = dyrDEGEN ; }
    else
    { retval = dyrOK ; }
    confirm = dy_confirmDualPivot(i,candk->ndx,abari,maxabari,p_abarj) ;
    if (confirm != dyrOK)
      retval = confirm ;
    FREE(incands) ;
    return (retval) ; }
/*
  We have multiple candidates and at least a chance to flip variables or
  correct slight infeasibilities. The question is how much effort we'll need
  to put into pricing.  If primal magnitudes are reasonable, we can attempt
  to select a sequence considering only the dual objective (i.e., find the
  longest sequence). If the primal magnitudes are over the top, we'd better
  pay attention to primal infeasibility.

  Note that the various infeasibilities are not considered in any way for
  strategy #1. The assignment below (quiet_nan) mollifies access checking and
  guarantees we'll see the effect of improper use.
*/
# ifdef DYLP_STATISTICS
  if (dy_stats != NULL) dy_stats->dmulti.nontrivial++ ;
# endif
  if (dpstrat == 1)
  { sel_retval = selectWithoutInf(i,abari,incands,candndx) ; }
  else
  { sel_retval = selectWithInf(i,incands,candndx,candinf,startinf) ; }
  if (sel_retval == FALSE)
  { FREE(incands) ;
    return (dyrFATAL) ; }
  bestflipcand = candndx[0] ;
  bestpivcand = candndx[1] ;
  lastpivcand = candndx[2] ;
  if (dpstrat != 1)
  { bestflipinf = candinf[0] ;
    bestpivinf = candinf[1] ;
    lastpivinf = candinf[2] ;
    starttotinf = startinf[0] ;
    startmaxinf = startinf[1] ; }
  else
  { bestflipinf = quiet_nan(0) ;
    bestpivinf = quiet_nan(0) ;
    lastpivinf = quiet_nan(0) ;
    starttotinf = quiet_nan(0) ;
    startmaxinf = quiet_nan(0) ; }
# if (!defined(DYLP_NDEBUG) || defined(DYLP_PARANOIA))
  if (dpstrat == 1)
  { if (flgon(dy_status[i],vstatBLLB))
    { starttotinf = vlb[i]-dy_x[i] ; }
    else
    { starttotinf = dy_x[i]-vub[i] ; } }
# endif
/*
  What shall we do?
    * dualmultipiv == 1 forces lastpivcand (max dual objective change)
    * dualmultipiv == 2 forces bestpivcand (min predicted infeasibility)
    * dualmultipiv == 3 will choose bestpivcand only if there's an actual
      reduction in infeasibility; otherwise, it'll choose lastpivcand,
      unless lastpivcand is degenerate, in which case it'll go back to
      bestpivcand.
  Note that bestpivcand > 0 iff lastpivcand > 0.

  If we have no pivot candidate, we'll consider doing only flips, but only if
  it'll actually reduce infeasibility. In the event that we have no pivot
  candidate and can't or choose not to flip, return the first candidate (which
  will be a mad pivot).
*/
  switch (dpstrat)
  { case 1:
    { bestcand = lastpivcand ;
      bestinf = lastpivinf ;
      break ; }
    case 2:
    { bestcand = bestpivcand ;
      bestinf = bestpivinf ;
      break ; }
    case 3:
    { if (bestpivcand > 0)
      { if (bestpivinf < startmaxinf)
	{ bestcand = bestpivcand ;
	  bestinf = bestpivinf ; }
	else
	if (incands[lastpivcand].ddelta > 0)
	{ bestcand = lastpivcand ;
	  bestinf = lastpivinf ; }
	else
	{ bestcand = bestpivcand ;
	  bestinf = bestpivinf ; } }
      else
      { bestcand = -1 ;
	bestinf = quiet_nan(0) ; }
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      FREE(incands) ;
      return (dyrFATAL) ; } }

  flipOnly = FALSE ;
  if (bestcand <= 0)
  { if (dy_opts->dpsel.allownopiv == TRUE && bestflipcand > 0 &&
	(dpstrat == 1 || bestflipinf < startmaxinf))
    { flipOnly = TRUE ;
      bestcand = bestflipcand ;
      bestinf = bestflipinf ; }
    else
    { bestcand = 1 ; } }
  candk = &incands[bestcand] ;

# ifdef DYLP_STATISTICS
    if (dy_stats != NULL && bestcand > 0)
    { dy_stats->dmulti.pivrnks += bestcand ;
      if (dy_stats->dmulti.maxrnk < bestcand)
	dy_stats->dmulti.maxrnk = bestcand ; } 
# endif
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 1)
  { j = candk->ndx ;
    infj = starttotinf ;
    if (dpstrat > 1)
    { for (ndx = 1 ; ndx < bestcand ; ndx++)
      { infj += incands[ndx].flip.inf ; }
      if (flipOnly == TRUE)
	infj += candk->flip.inf ;
      else
	infj += candk->piv.inf ; }
    dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\n  selected %s #%d, %s (%d), est. maxinf =  %g, totinf = %g",
	        (flipOnly == TRUE)?"flip":"pivot",
	        bestcand,consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	        bestinf,infj) ; }
# endif
/*
  # ifdef DYLP_PARANOIA
    if (dy_opts->dpsel.strat >= 2 && incands[bestcand].madpiv == FALSE)
    { predictiter = dy_lp->d2.iters ;
      predicttotinf = starttotinf ;
      for (ndx = 1 ; ndx < bestcand ; ndx++)
      { predicttotinf += incands[ndx].flip.inf ; }
      if (flipOnly == TRUE)
	predicttotinf += candk->flip.inf ;
      else
	predicttotinf += candk->piv.inf ;
      predictmaxinf = bestinf ; }
  # endif
*/
/*
  We've made our choice. But ... if this is a degenerate pivot, and
  antidegeneracy can be activated, let's hold off for a moment.
*/
  if (candk->ddelta == 0)
  { if (dy_opts->degen == TRUE && dy_opts->degenpivlim < dy_lp->degenpivcnt)
    { *p_xjndx = candk->ndx ;
      *p_indir = candk->pivdir ;
      retval = dyrDEGEN ;
      FREE(incands) ;
      return (retval) ; } }
/*
  One last check --- confirm that abar<ij> from abar<i> agrees with abar<ij>
  from abar<j> = inv(B)a<j>.
*/
  confirm = dy_confirmDualPivot(i,candk->ndx,abari,maxabari,p_abarj) ;
  if (confirm != dyrOK)
  { *p_xjndx = candk->ndx ;
    *p_indir = candk->pivdir ;
    FREE(incands) ;
    return (confirm) ; }
/*
  If there are flips to perform, we need to handle
  them now before we can return and let dy_dualpivot handle the actual
  pivot. Rather than doing a separate ftran for each a<j>, accumulate
  delta<j>*a<j> in accum<j> for all variables and then do a single ftran. But
  do update the objective as we process each column.

  Reverse pivots will not flip.
  
  When we're doing only flips, boost bestcand by 1 so that we'll
  flip the candidate specified by bestcand, instead of reserving it for
  pivoting (and undo the increment after the loop ends).
*/
  swing = FALSE ;
  reqchk = FALSE ;
  if (bestcand > 1 || flipOnly == TRUE)
  { if (flipOnly == TRUE) bestcand++ ;
    accumj = (double *) CALLOC((m+1),sizeof(double)) ;
    for (ndx = 1 ; ndx < bestcand ; ndx++)
    { candk = &incands[ndx] ;
      j = candk->ndx ;
      statj = dy_status[j] ;
      if (candk->rev == TRUE)
      { 
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.pivoting >= 3)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\tpassing reverse %s (%d) %s = %g.",
		      consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		      dy_prtvstat(statj),dy_x[j]) ; }
#	endif
	continue ; }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivoting >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tflipping %s (%d) %s = %g to ",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		    dy_prtvstat(statj),dy_x[j]) ; }
#     endif
      comflg(statj,vstatNBUB|vstatNBLB) ;
      if (flgon(statj,vstatNBLB))
      { xj = vlb[j] ; }
      else
      { xj = vub[j] ; }
      dy_status[j] = statj ;
      dy_x[j] = xj ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivoting >= 3)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"%s = %g.",
		    dy_prtvstat(dy_status[j]),dy_x[j]) ; }
#     endif
#     ifdef DYLP_STATISTICS
      if (dy_stats != NULL) dy_stats->dmulti.flips++ ;
#     endif

      deltaj = candk->flip.delta ;
      if (dy_ddegenset[j] == 0)
      { dy_lp->z += dy_cbar[j]*deltaj ; }
      if (consys_mulaccumcol(dy_sys,j,deltaj,accumj) == FALSE)
      { errmsg(122,rtnnme,dy_sys->nme,
	       "column",consys_nme(dy_sys,'v',j,TRUE,NULL),j) ;
	retval = dyrFATAL ;
	break ; } }
    if (retval == dyrFATAL)
    { if (accumj != NULL) FREE(accumj) ;
      FREE(incands) ;
      return (retval) ; }
    if (flipOnly == TRUE) bestcand-- ;
/*
  Ftran the accumulated SUM{j in cands}delta<j>*a<j> and pass it to
  dy_updateprimals to update the basic variables.
*/
    dy_ftran(accumj,FALSE) ;
    upd_retval = dy_updateprimals(candk->ndx,1.0,accumj) ;
    switch (upd_retval)
    { case dyrOK:
      { break ; }
      case dyrSWING:
      { swing = TRUE ;
	break ; }
      case dyrREQCHK:
      { reqchk = TRUE ;
	break ; }
      default:
      { retval = dyrFATAL ;
	break ;  } }
    FREE(accumj) ;
/*
  There's a remote chance that the flips will reduce x<i> exactly to bound.
  (Yes, it's happened.) Catch this at the end of the loop, and claim flips only
  when it happens, abandoning the prospective pivot.
*/
    if (flgoff(dy_status[i],vstatBUUB|vstatBLLB))
    { flipOnly = TRUE ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.pivoting >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  whoa! %s (%d) = %g %s after flips; cancelling pivot.",
		    consys_nme(dy_sys,'v',i,FALSE,NULL),i,dy_x[i],
		    dy_prtvstat(dy_status[i])) ; }
#     endif
    }
    if (retval == dyrFATAL)
    { FREE(incands) ;
      return (retval) ; } }
/*
  If we have a pivot to recommend, set up the return values so we look like
  the usual dualin return. If we don't have a pivot to recommend, return
  dyrRESELECT, to indicate that we should try reselecting the leaving
  variable.

  Note that here we prefer dyrMADPIV over dyrDEGEN --- if antidegeneracy was
  active, we'd have taken the earlier return. We're going to use this pivot,
  degenerate or not, so it's more important to know that it's unstable.
*/
  if (flipOnly == FALSE)
  { candk = &incands[bestcand] ;
    *p_xjndx = candk->ndx ;
    *p_indir = candk->pivdir ;
    if (candk->madpiv == TRUE)
    { retval = dyrMADPIV ; }
    else
    if (candk->ddelta == 0)
    { retval = dyrDEGEN ; }
    else
    { retval = dyrOK ; } }
  else
  { retval = dyrRESELECT ; }

  FREE(incands) ;
  return (retval) ; }


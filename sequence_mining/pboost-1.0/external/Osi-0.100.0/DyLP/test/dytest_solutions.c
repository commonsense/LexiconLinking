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
  This file contains routines to test dylp's solution routines: dy_colPrimals,
  dy_rowPrimals, dy_colDuals, and dy_rowDuals.
*/

#include "dylp.h"

extern ioid dy_logchn ;
extern bool dy_gtxecho ;



int dytest_rowDuals (lpprob_struct *main_lp, lptols_struct *main_lptols,
		     lpopts_struct *main_lpopts)
/*
  This routine checks the dual variables returned by dy_rowDuals. It checks
  that y<i> = c<B>(inv(B))<i>. Columns of the basis inverse are obtained from
  the routine dy_betaj.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if y = c<B>inv(B), error count otherwise.
*/

{ int i,j,k,m,n ;
  consys_struct *sys ;
  flags *status ;

  double *y ;

  double *cB ;
  int *basis2sys ;
  basisel_struct *basisVec ;
  int basisLen ;

  double *betai ;
  double cBdotbetai ;
  int errcnt ;

  char *rtnnme = "dytest_rowDuals" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.soln >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: checking y = c<B>inv(B) using %s (%d x %d).",
		rtnnme,sys->nme,m,n) ;
    if (main_lpopts->print.soln >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "  basis contains %d entries.\n",main_lp->basis->len) ; } }
# endif
/*
  Acquire the row duals and the status vector.
*/
  y = NULL ;
  dy_rowDuals(main_lp,&y) ;
  status = main_lp->status ;
/*
  Make a vector c<B> of objective coefficients in basis order. This is
  considerably easier than creating a basis matrix (as is done for tableau
  testing). By construction, the basic variables for inactive constraints
  are the logicals, which have an objective coefficient of zero, and this
  is how cB and basis2sys are initialised. All that need be done for c<B>
  is to change the entries that are associated with architecturals. For the
  basis, we need to set all entries (logicals can be basic out of natural
  position). Recall that basic logicals are represented by negative indices.
*/
  cB = (double *) MALLOC((m+1)*sizeof(double)) ;
  basis2sys = (int *) MALLOC((m+1)*sizeof(int)) ;
  for (i = 1 ; i <= m ; i++)
  { cB[i] = 0.0 ;
    basis2sys[i] = -i ; }
  basisLen = main_lp->basis->len ;
  basisVec = main_lp->basis->el ;
  for (k = 1 ; k <= basisLen ; k++)
  { i = basisVec[k].cndx ;
    j = basisVec[k].vndx ;
    if (j > 0)
    { cB[i] = sys->obj[j] ;
      basis2sys[i] = j ; }
    else
    { basis2sys[i] = j ; } }
# ifndef DYLP_NDEBUG
  if (main_lpopts->print.soln >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n\tc<B> =") ;
    k = 0 ;
    for (i = 1 ; i <= m ; i++)
    { if (cB[i] != 0)
      { if ((++k)%4 == 0)
	{ k = 0 ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\n\t      ") ; }
	j = basis2sys[i] ;
	dyio_outfmt(dy_logchn,dy_gtxecho," (%d %g  %s %d)",
		    i,cB[i],consys_nme(sys,'v',j,FALSE,NULL),j) ; } } }
#   endif
/*
  Now step through the rows (equivalently, walk the basis) and see if
  y<i> = c<B>beta<j>, where beta<j> is the column of inv(B) such that
  x<j> is basic in pos'n i.
*/
  errcnt = 0 ;
  betai = NULL ;
  for (i = 1 ; i <= m ; i++)
  { j = basis2sys[i] ;
    if (dy_betaj(main_lp,j,&betai) == FALSE)
    { errcnt++ ;
      if (j < 0)
      { j = n-j ; }
      errmsg(952,rtnnme,sys->nme,"column",i,"variable",
	     consys_nme(sys,'v',j,FALSE,NULL),j-n) ;
      continue ; }
    cBdotbetai = 0 ;
    for (k = 1 ; k <= m ; k++)
    { 
      /*
         dyio_outfmt(dy_logchn,dy_gtxecho,
    		  "\n %s (%d) %g * %g",
		  consys_nme(sys,'c',k,FALSE,NULL),k,cB[k],betai[k]) ;
      */
      cBdotbetai += cB[k]*betai[k] ; }
    if (fabs(cBdotbetai-y[i]) > main_lptols->cost)
    { errcnt++ ;
      if (j < 0)
      { j = n-j ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  ERROR: pos'n %d %s (%d) c<B> dot beta<j> = %g; ",
		    i,consys_nme(sys,'v',j,FALSE,NULL),j-n,cBdotbetai) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"expected %g; err %g, tol %g.",
		    y[i],(cBdotbetai-y[i]),main_lptols->zero) ; } }
/*
  Free up space and report the result.
*/
  if (cB != NULL) FREE(cB) ;
  if (basis2sys != NULL) FREE(basis2sys) ;
  if (betai != NULL) FREE(betai) ;
  if (y != NULL) FREE(y) ;

  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: found %d errors testing y = c<B>inv(B).\n",
		rtnnme,errcnt) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: pass y = c<B>inv(B).\n",rtnnme) ; }

  return (errcnt) ; }



int dytest_colDuals (lpprob_struct *main_lp, lptols_struct *main_lptols,
		     lpopts_struct *main_lpopts)
/*
  This routine checks the dual variables returned by dy_colDuals (more
  usually called the reduced costs of the architectural variables).
  
  It checks that cbar<N> = c<N> - yN, where y is the vector of row duals
  returned by dy_rowDuals and N is the set of nonbasic architectural columns
  of A (or the matching index set, as appropriate). It also checks that the
  reduced cost is in agreement with the status.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if cbar<N> = c<N> - yN, error count otherwise.
*/

{ int j,m,n ;
  flags statj ;

  consys_struct *sys ;
  double *obj ;
  double *y ;
  flags *status ;

  double *cbarN ;
  double cbarj ;

  int errcnt ;
  bool staterr ;
  char *errstring ;

  char *rtnnme = "dytest_colDuals" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.soln >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: checking cbar<N> = c<N> - yN using %s (%d x %d).",
		rtnnme,sys->nme,m,n) ; }
# endif
/*
  Acquire the row duals, column duals, status vector, and objective.
*/
  y = NULL ;
  dy_rowDuals(main_lp,&y) ;
  cbarN = NULL ;
  dy_colDuals(main_lp,&cbarN) ;
  status = main_lp->status ;
  obj = sys->obj ;
/*
  Now step through the columns checking that cbar<j> = c<j> = dot(y,a<j>).
  Also check to see that the sign is correct for the status of the variable
  in a minimisation problem.  For status values not listed (vstatSB and any
  of the basic status codes), there's no `correct' sign.
*/
  errcnt = 0 ;
  for (j = 1 ; j <= n ; j++)
  { cbarj = obj[j] - consys_dotcol(sys,j,y) ;
    statj = getflg(status[j],vstatSTATUS) ;
    if (fabs(cbarj-cbarN[j]) > main_lptols->cost)
    { errcnt++ ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\nERROR: col %s (%d) %s cbar<%d> = %g; expected %g;",
		  consys_nme(sys,'v',j,FALSE,NULL),j,dy_prtvstat(statj),
		  j,cbarj,cbarN[j]) ;
      dyio_outfmt(dy_logchn,dy_gtxecho," error %g, tol %g.",
		  fabs(cbarj),main_lptols->cost) ; }
    staterr = FALSE ;
    switch (statj)
    { case vstatNBLB:
      { if (cbarj < -main_lptols->zero)
	{ staterr = TRUE ;
	  errstring = "positive" ; }
	break ; }
      case vstatNBUB:
      { if (cbarj > main_lptols->zero)
	{ staterr = TRUE ;
	  errstring = "negative" ; }
	break ; }
      case vstatNBFR:
      { if (fabs(cbarj) > main_lptols->zero)
	{ staterr = TRUE ;
	  errstring = "zero" ; }
	break ; }
      default:
      { break ; } }
    if (staterr == TRUE)
    { errcnt++ ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\nERROR: col %s (%d) %s cbar<%d> = %g; should be %s.",
		  consys_nme(sys,'v',j,FALSE,NULL),j,dy_prtvstat(statj),
		  j,cbarj,errstring) ; } }
/*
  Free up space and report the result.
*/
  if (y != NULL) FREE(y) ;
  if (cbarN != NULL) FREE(cbarN) ;

  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: found %d errors testing cbar<N> = c<N> - yN.\n",
		rtnnme,errcnt) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: pass cbar<N> = c<N> - yN.\n",rtnnme) ; }

  return (errcnt) ; }




int dytest_colPrimals (lpprob_struct *main_lp, lptols_struct *main_lptols,
		       lpopts_struct *main_lpopts)
/*
  This routine checks the values of the primal architectural variables
  returned by dy_colPrimals.
  
  For basic variables x<B>, the routine checks
    x<B> =  inv(B)b - inv(B)Nx<N>
  To do this, the routine accumulates the values of the basic variables
  during the column scan. When the current column is basic in pos'n i, the
  routine calculates dot(beta<i>,b) and adds it to the total. When the current
  column is nonbasic, the routine calculates abar<j>x<j> and subtracts it
  from the total. Just to make things really annoying, we have to account for
  nonbasic bounded slacks due to range constraints tight at their lower bound
  (which makes the slack nonbasic at its upper bound).

  For a nonbasic variable, the routine checks the value of x<j> against the
  bound specified by the status of x<j>.

  Parameters:
    main_lp:     the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if the values check out, error count otherwise.
*/

{ int i,i_bpos,j,k,m,n ;
  flags statj,stati ;
  double xj,lbj,ubj,betaidotb ;

  consys_struct *sys ;
  flags *status,*logstatus ;
  double *rhs,*rhslow,*vlb,*vub,*betai,*xB,*abarj ;
  contyp_enum *ctyp ;
  basisel_struct *basis ;

  double *x ;

  int berrs,nberrs ;
  bool staterr ;

  char *rtnnme = "dytest_colPrimals" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.soln >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	  "\n%s: checking primal architectural variables using %s (%d x %d).",
	  rtnnme,sys->nme,m,n) ; }
# endif
/*
  Acquire the variable bound and status vectors, the constraint type, rhs,
  and rhslow vectors, the basis vector, and the values of the primal
  architectural variables. Allocate a vector to accumulate x<B>
*/
  x = NULL ;
  dy_colPrimals(main_lp,&x) ;
  basis = main_lp->basis->el ;
  status = main_lp->status ;
  ctyp = sys->ctyp ;
  rhs = sys->rhs ;
  rhslow = sys->rhslow ;
  vlb = sys->vlb ;
  vub = sys->vub ;
  xB = (double *) CALLOC((m+1),sizeof(double)) ;
/*
  Now step through the columns checking the values in x. For a variable basic
  in pos'n i, add dot(beta<i>,b) to the running total for the basic variable.
  
  For a nonbasic variable, confirm that the value, bound, and status agree.
  Then subtract abar<j>x<j> from x<B> if x<j> is at a nonzero bound.

  The only nonbasic status code not explicitly listed is SB (superbasic). This
  really should never appear. The only legitimate reason is that dylp patched
  the basis in primal phase II and then discovered the problem to be unbounded
  before the SB variable could be pivoted back into the basis. This is
  sufficiently exotic to deserve a message.
*/
  berrs = 0 ;
  nberrs = 0 ;
  betai = NULL ;
  abarj = NULL ;
  for (j = 1 ; j <= n ; j++)
  { statj = status[j] ;
    xj = x[j] ;
    if (((int) statj) < 0)
    { k = -((int) statj) ;
      i_bpos = basis[k].cndx ;
      if (dy_betai(main_lp,i_bpos,&betai) == FALSE)
      { berrs++ ;
	errmsg(952,rtnnme,sys->nme,"row",i_bpos,"variable",
	       consys_nme(sys,'v',j,FALSE,NULL),j) ;
	continue ; }
      betaidotb = 0 ;
      for (i = 1 ; i <= m ; i++)
      { betaidotb += betai[i]*rhs[i] ; }
      xB[i_bpos] += betaidotb ; }
    else
    { staterr = FALSE ;
      lbj = vlb[j] ;
      ubj = vub[j] ;
      statj = getflg(statj,vstatSTATUS) ;
      switch (statj)
      { case vstatNBLB:
	{ if (fabs(xj-lbj) > main_lptols->zero)
	  { staterr = TRUE ;
	    betaidotb = lbj ; }
	  break ; }
	case vstatNBUB:
	{ if (fabs(xj-ubj) > main_lptols->zero)
	  { staterr = TRUE ;
	    betaidotb = ubj ; }
	  break ; }
	case vstatNBFX:
	{ if (fabs(xj-lbj) > main_lptols->zero)
	  { staterr = TRUE ;
	    betaidotb = lbj ; }
	  break ; }
	case vstatNBFR:
	{ if (fabs(xj) > main_lptols->zero)
	  { staterr = TRUE ;
	    betaidotb = 0.0 ; }
	  break ; }
	default:
	{ staterr = TRUE ;
	  betaidotb = quiet_nan(42.0L) ;
	  break ; } }
    if (staterr == TRUE)
    { nberrs++ ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\nERROR: %s col %s (%d) = %g; expected %g;",
		  dy_prtvstat(statj),consys_nme(sys,'v',j,FALSE,NULL),j,
		  xj,betaidotb) ;
      dyio_outfmt(dy_logchn,dy_gtxecho," error %g, tol %g.",
		  fabs(xj-betaidotb),main_lptols->zero) ;
      continue ; }
    if (xj == 0.0) continue ;
    if (dy_abarj(main_lp,j,&abarj) == FALSE)
    { nberrs++ ;
      errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	     consys_nme(sys,'v',j,FALSE,NULL),j) ;
      continue ; }
    for (k = 1 ; k <= m ; k++)
    { xB[k] -= abarj[k]*xj ; } } }
/*
  But wait! We're not quite done. We need to account for bounded slacks
  associated with range constraints. If the constraint is tight at its lower
  bound, the slack is nonbasic at its upper bound.
*/
  logstatus = NULL ;
  dy_logStatus(main_lp,&logstatus) ;
  for (i = 1 ; i <= m ; i++)
  { stati = getflg(logstatus[i],vstatSTATUS) ;
    if (ctyp[i] == contypRNG && stati == vstatNBUB)
    { xj = rhs[i]-rhslow[i] ;
      if (dy_abarj(main_lp,-i,&abarj) == FALSE)
      { nberrs++ ;
	errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	       consys_nme(sys,'v',n+i,FALSE,NULL),i) ;
	continue ; }
      for (k = 1 ; k <= m ; k++)
      { xB[k] -= abarj[k]*xj ; } } }
/*
  Scan the variables one more time and check the values of the basic variables.
*/
  for (j = 1 ; j <= n ; j++)
  { statj = status[j] ;
    xj = x[j] ;
    if (((int) statj) < 0)
    { k = -((int) statj) ;
      i_bpos = basis[k].cndx ;
      if (fabs(xj-xB[i_bpos]) > main_lptols->zero)
      { berrs++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\nERROR: basis pos'n %d %s (%d) = %g; expected %g;",
		    i_bpos,consys_nme(sys,'v',j,FALSE,NULL),j,xj,xB[i_bpos]) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," error %g, tol %g.",
		    fabs(xj-xB[i_bpos]),main_lptols->zero) ; } } }
/*
  Free up space and report the result.
*/
  if (logstatus != NULL) FREE(logstatus) ;
  if (abarj != NULL) FREE(abarj) ;
  if (xB != NULL) FREE(xB) ;
  if (betai != NULL) FREE(betai) ;
  if (x != NULL) FREE(x) ;

  if ((berrs+nberrs) != 0)
  { if (berrs != 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: found %d errors testing x<B> = inv(B)b.\n",
		  rtnnme,berrs) ; }
    if (nberrs != 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n%s: found %d errors testing x<N> against bounds & status.\n",
	      rtnnme,nberrs) ; } }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: pass test of primal architectural variable values.\n",
		rtnnme) ; }

  return (berrs+nberrs) ; }



int dytest_rowPrimals (lpprob_struct *main_lp, lptols_struct *main_lptols,
		       lpopts_struct *main_lpopts)
/*
  This routine checks the ind<B> and x<B> vectors returned by dy_rowPrimals.
  It first cross-checks the basis, status and indB arrays, bailing out if the
  cross-checks fail.
  
  Next it checks the values of the basic variables, architectural and logical.

  For basic variables x<B>, the routine checks x<B> =  inv(B)b - inv(B)Nx<N>
    To do this, it first walks the rows of the constraint system and
    initialises x<B> with dot(beta<i>,b). Then it walks the columns and
    accumulates the contributions abar<j>x<j> from nonzero nonbasic
    variables.  Finally, it walks the rows again and subtracts the
    contributions from nonbasic bounded logicals (due to range constraints
    tight at the lower bound).

  Parameters:
    main_lp:     the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if the basic variables validate, error count otherwise.
*/

{ int i,j,k,m,n,i_basis ;
  flags statj,stati ;
  double xj,betaidotb,tol ;

  consys_struct *sys ;
  flags *status,*logstatus ;
  double *rhs,*rhslow,*vlb,*vub,*betai,*xBaccum,*abarj ;
  contyp_enum *ctyp ;
  basisel_struct *basis ;
  int basisLen ;

  double *xB ;
  int *indB ;

  int berrs,nberrs,inderrs ;

  char *rtnnme = "dytest_rowPrimals" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.soln >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	  "\n%s: checking primal basic variables using %s (%d x %d).",
	  rtnnme,sys->nme,m,n) ; }
# endif
/*
  Acquire the variable bound and status vectors, the constraint type, rhs,
  and rhslow vectors, and the basis vector.
*/
  basisLen = main_lp->basis->len ;
  basis = main_lp->basis->el ;
  status = main_lp->status ;
  ctyp = sys->ctyp ;
  rhs = sys->rhs ;
  rhslow = sys->rhslow ;
  vlb = sys->vlb ;
  vub = sys->vub ;
/*
  Call dy_rowPrimals to acquire x<B> (values of basic variables) and ind<B>
  (indices of basic variables).
*/
  xB = NULL ;
  indB = NULL ;
  dy_rowPrimals(main_lp,&xB,&indB) ;
/*
  Validate ind<B>, status, and basis against each other, within the limits of
  each.

  IndB specifies basic variables in row order. Logicals are specified as the
  negative of the row. IndB contains an entry for every constraint. By
  construction, the basic variable for an inactive constraint should be the
  logical for the constraint.

  Basis has one entry for each active constraint. Each entry in basis
  specifies a constraint and a basic variable. Basic logicals are specified
  by the negative of the constraint index. Then for an active constraint i
  and a basis entry k such that basis[k].cndx == i, indB[i] == basis[k].vndx.

  Status only contains information on architecturals. A basic architectural
  is specified as the negative of its entry in the basis vector. Thus
  basis[-status[j]].vndx == j. 
*/
  inderrs = 0 ;
  for (i = 1 ; i <= m ; i++)
  { 
/*
  Scan the basis vector for an entry for this constraint. If it's not present,
  assume the constraint is inactive.
*/
    i_basis = -1 ;
    for (k = 1 ; k <= basisLen ; k++)
    { if (basis[k].cndx == i)
      { i_basis = k ;
	break ; } }
    j = indB[i] ;
/*
  Inactive constraints should specify the associated logical as the basic
  variable.
*/
    if (i_basis < 0)
    { if (j > 0)
      { inderrs++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\nERROR: constraint %s (%d)",
		    consys_nme(sys,'c',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; basis entry = %d; should specify a logical.",j) ; }
      else
      if (-j != i)
      { inderrs++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\nERROR: basis[%d] (%s)",
		    i,consys_nme(sys,'c',i,FALSE,NULL)) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," is %s (%d);",
		    consys_nme(sys,'c',n-j,FALSE,NULL),-j) ;
	dyio_outfmt(dy_logchn,dy_gtxecho," expected %s (%d).",
		    consys_nme(sys,'c',n+i,FALSE,NULL),i) ; } }
/*
  The constraint is active. We should have indB[i] = basis[i_basis].vndx. It
  takes way more work than it should to construct the error message.
*/
    else
    { k = basis[i_basis].vndx ;
      if (j != k)
      { inderrs++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\nERROR: constraint %s (%d)",
		    consys_nme(sys,'c',i,FALSE,NULL),i) ;
	statj = (k < 0)?(n-k):(k) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; basis[%d] specifies %s (%d)",
		    i_basis,consys_nme(sys,'v',statj,FALSE,NULL),k) ;
	statj = (j < 0)?(n-j):(j) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "; indB[%d] specifies %s (%d); they should agree.",
		    i,consys_nme(sys,'v',statj,FALSE,NULL),j) ; }
/*
  If the basic variable k is an architectural, status[k] should agree that it's
  basic and point to the basis vector entry.
*/
      if (k > 0)
      { statj = -((int) status[k]) ;
	if (i_basis != statj)
	{ inderrs++ ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,"\nERROR: constraint %s (%d)",
		      consys_nme(sys,'c',i,FALSE,NULL),i) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "; status[%d] = %d but basis[%d].vndx = %d",
		      k,statj,i_basis,k) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "; they should point to each other.") ; } } } }
  if (inderrs > 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: found %d errors cross-checking basis index vectors.\n",
		rtnnme,inderrs) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"\tTests of basic variable values not performed.\n") ;
    if (xB != NULL) FREE(xB) ;
    if (indB != NULL) FREE(indB) ;
    return (inderrs) ; }
/*
  Now we know the index arrays are correct and we can use them with
  confidence.  Step through the rows, placing the initial component
  dot(beta<i>,b) into each position.
*/
  xBaccum = (double *) CALLOC((m+1),sizeof(double)) ;
  berrs = 0 ;
  betai = NULL ;
  for (i = 1 ; i <= m ; i++)
  { if (dy_betai(main_lp,i,&betai) == FALSE)
    { berrs++ ;
      j = indB[i] ;
      if (j < 0)
      { statj = n-j ; }
      else
      { statj = j ; }
      errmsg(952,rtnnme,sys->nme,"row",i,"basic variable",
	     consys_nme(sys,'v',statj,FALSE,NULL),j) ;
      continue ; }
    betaidotb = 0 ;
    for (k = 1 ; k <= m ; k++)
    { betaidotb += betai[k]*rhs[k] ; }
    xBaccum[i] += betaidotb ; }
/*
  Now step through the columns. Subtract abar<j>x<j> from x<B> if x<j> is at
  a nonzero bound. Anything other than the enumerated status codes is
  extraordinary. vstatSB might be correct if dylp declared unboundedness
  immediately after refactoring in primal phase II, but that's such an unlikely
  coincidence it deserves attention. Anything else is outright wrong.
*/
  nberrs = 0 ;
  abarj = NULL ;
  for (j = 1 ; j <= n ; j++)
  { statj = status[j] ;
    if (((int) statj) < 0) continue ;
    statj = getflg(statj,vstatSTATUS) ;
    switch (statj)
    { case vstatNBLB:
      case vstatNBFX:
      { xj = vlb[j] ;
	break ; }
      case vstatNBUB:
      { xj = vub[j] ;
	break ; }
      case vstatNBFR:
      { xj = 0.0 ;
	break ; }
      default:
      { nberrs++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"\nERROR: constraint %s (%d)",
		    consys_nme(sys,'c',i,FALSE,NULL),i) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"; status of %s (%d) is %s.",
		    consys_nme(sys,'v',j,FALSE,NULL),j,dy_prtvstat(statj)) ;
	xj = 0.0 ;
	break ; } }
    if (xj == 0.0) continue ;
    if (dy_abarj(main_lp,j,&abarj) == FALSE)
    { nberrs++ ;
      errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	     consys_nme(sys,'v',j,FALSE,NULL),j) ;
      continue ; }
    for (k = 1 ; k <= m ; k++)
    { xBaccum[k] -= abarj[k]*xj ; } }
/*
  We're not quite done. We need to account for bounded slacks
  associated with range constraints. If the constraint is tight at its lower
  bound, the slack is nonbasic at its upper bound.
*/
  logstatus = NULL ;
  dy_logStatus(main_lp,&logstatus) ;
  for (i = 1 ; i <= m ; i++)
  { stati = getflg(logstatus[i],vstatSTATUS) ;
    if (ctyp[i] == contypRNG && stati == vstatNBUB)
    { xj = rhs[i]-rhslow[i] ;
      if (dy_abarj(main_lp,-i,&abarj) == FALSE)
      { nberrs++ ;
	errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	       consys_nme(sys,'v',n+i,FALSE,NULL),i) ;
	continue ; }
      for (k = 1 ; k <= m ; k++)
      { xBaccum[k] -= abarj[k]*xj ; } } }
/*
  Scan the rows one more time and check the values of the basic variables.
  Scale this test just a bit so we don't get spurious indications due to
  roundoff. The average of the two values seems safest as a scaling factor.
*/
  for (i = 1 ; i <= m ; i++)
  { tol = ((fabs(xBaccum[i])+fabs(xB[i]))/2)+1 ;
    if (fabs(xBaccum[i]-xB[i]) > tol*main_lptols->zero)
    { berrs++ ;
      j = indB[i] ;
      if (j < 0)
      { statj = n-j ; }
      else
      { statj = j ; }
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\nERROR: basis pos'n %d %s (%d) = %g; expected %g;",
		  i,consys_nme(sys,'v',statj,FALSE,NULL),j,xB[i],xBaccum[i]) ;
      dyio_outfmt(dy_logchn,dy_gtxecho," error %g, tol %g.",
		  fabs(xB[i]-xBaccum[i]),main_lptols->zero) ; } }
/*
  Free up space and report the result.
*/
  if (logstatus != NULL) FREE(logstatus) ;
  if (abarj != NULL) FREE(abarj) ;
  if (xB != NULL) FREE(xB) ;
  if (indB != NULL) FREE(indB) ;
  if (xBaccum != NULL) FREE(xBaccum) ;
  if (betai != NULL) FREE(betai) ;

  if ((berrs+nberrs) != 0)
  { if (berrs != 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: found %d errors testing x<B> = inv(B)b.\n",
		rtnnme,berrs) ; }
    if (nberrs != 0)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n%s: found %d errors attempting to use nonbasic variables.\n",
	      rtnnme,nberrs) ; } }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: pass test of primal basic variable values.\n",
		rtnnme) ; }

  return (berrs+nberrs) ; }


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
  This file contains routines to test dylp's tableau routines: dy_betaj,
  dy_betai, dy_abarj, and dy_abari.
*/

#include "dylp.h"

extern ioid dy_logchn ;
extern bool dy_gtxecho ;



static int compare_basisel (const void *el1, const void *el2)
/*
  Helper routine to sort the basis vector by row index. Since the basis vector
  cannot have equal indices, this routine will never return 0.

  Parameters: a pair of basisel_struct's

  Returns: -1 if el1 < el2
	    0 if el1 = el2
	    1 if el1 > el2
*/
{ const basisel_struct *basisel1,*basisel2 ;

  basisel1 = (const basisel_struct *) el1 ;
  basisel2 = (const basisel_struct *) el2 ;

  if (basisel1->cndx < basisel2->cndx)
  { return (-1) ; }
  else
  { return (1) ; } }



static consys_struct *create_basis (lpprob_struct *main_lp,
				    lptols_struct *main_lptols,
			     	    lpopts_struct *main_lpopts,
			     	    int **p_basis2sys)
/*
  This routine builds a basis matrix to match the basis returned by dylp.
  Very handy for testing dylp's tableau routines.

  If this basis matrix is to be directly useable in calculations involving
  the basis inverse, we need to install rows and columns in the proper
  order.  The row order must correspond to the row order of the external
  system and the basic variable that is associated with row i must be placed
  in column i.
  
  Each entry in the basis vector returned from dylp specifies indices, in the
  external frame of reference, for the constraint and variable associated
  with the basis position. However, the order of basis entries matches the
  constraint order in dylp's active system, which has no predictable
  relationship to the order of constraints in the external system.  In
  addition, if dylp is working with a partial system, not all constraints
  will be included. So we have some work to do.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure
    basis2sys:	 (i) vector to hold translation from basis column order to
		     external system column order; allocated if NULL.
		 (o) completed translation vector

    Note that basis2sys is an attached vector for basis and will be freed when
    the constraint system is freed.

  Returns: pointer to a basis matrix, or NULL if an error occurs during
	   construction
*/

{ int m,n,i,j,k ;
  double infty ;
  consys_struct *sys ;

  consys_struct *basis ;
  basisel_struct *basisVec ;
  int basisLen,lastRow ;

  int *basis2sys ;

  pkvec_struct *ai,*aj ;

  char *rtnnme = "create_basis" ;

# ifndef DYLP_NDEBUG
  int printlvl ;

  printlvl = maxx(main_lpopts->print.tableau,main_lpopts->print.rays) ;
  /* printlvl = 3 ; */
# endif

/*
  Do a little initialisation.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;
  infty = main_lptols->inf ;

# ifndef DYLP_NDEBUG
  if (printlvl >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: generating basis matrix from %s (%d x %d).\n",
		rtnnme,sys->nme,m,n) ; }
# endif
/*
  Construct the basis based on the solution in main_lp. All we need here is the
  bare coefficient matrix. None of the usual attached vectors are required.
*/
  basisLen = main_lp->basis->len ;
# ifndef DYLP_NDEBUG
  if (printlvl >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"  basis contains %d entries.\n",basisLen) ; }
# endif
/*
  Create a constraint system to hold the basis.
*/
  basis = consys_create("basisMtx",0,0,m,m,infty) ;
  if (basis == NULL)
  { errmsg(152,rtnnme,"basis") ;
    return (FALSE) ; }
# ifndef DYLP_NDEBUG
  if (printlvl >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"  created system %s (%d x %d).\n",
		basis->nme,basis->rowsze,basis->colsze) ; }
# endif
/*
  Copy over the row headers for all constraints.
*/
  ai = pkvec_new(0) ;
  for (i = 1 ; i <= m ; i++)
  { if (consys_getrow_pk(sys,i,&ai) == FALSE)
    { errmsg(122,rtnnme,sys->nme,"row",consys_nme(sys,'c',i,FALSE,NULL),i) ;
      if (ai != NULL) pkvec_free(ai) ;
      consys_free(basis) ;
      return (NULL) ; }
    if (consys_addrow_pk(basis,'a',0,ai,0,0,NULL,NULL) == FALSE)
    { errmsg(112,rtnnme,basis->nme,"add row","constraint",
	     consys_nme(sys,'c',i,FALSE,NULL),i) ;
      if (ai != NULL) pkvec_free(ai) ;
      consys_free(basis) ;
      return (NULL) ; }
# ifndef DYLP_NDEBUG
  if (printlvl >= 5)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"  added %s (%d) to basis as row %d.\n",
		consys_nme(sys,'c',i,FALSE,NULL),i,ai->ndx) ; }
# endif
  }
  if (ai != NULL) pkvec_free(ai) ;
/*
  Now install the basic columns. As explained at the head of the routine, we
  have some work to do.

  First, make a copy of the basis vector and sort it by row index, so that
  the basis entries we have match the order in sys.

  Then we can walk the basis and install the columns in order. Basic logicals
  that were part of the basis in dylp's active system are represented by the
  negative of the index of the associated constraint. When we encounter one,
  synthesize a unit column for the logical.

  Except ... the basis will not mention inactive constraints.  Each time we
  encounter a gap in constraint indices, install unit columns as needed.
  We're doing a sort of on-the-fly activation of inactive constraints, using
  the logical as the basic variable.

  With all that, we'll need a translation vector, to take indices in the
  basis frame and translate them to indices in the sys frame.
*/
  if (consys_attach(basis,
		    CONSYS_ROW,sizeof(int),((void **) p_basis2sys)) == FALSE)
  { errmsg(100,rtnnme,basis->nme,"basis2sys") ;
    consys_free(basis) ;
    return (NULL) ; }
  basis2sys = *p_basis2sys ;

  basisVec = CALLOC((basisLen+1),sizeof(basisel_struct)) ;
  memcpy(basisVec,main_lp->basis->el,(basisLen+1)*sizeof(basisel_struct)) ;
  qsort(((void *) &basisVec[1]),
	basisLen,sizeof(basisel_struct),compare_basisel) ;

  aj = pkvec_new(sys->maxcollen) ;
  lastRow = 0 ;
  for (k = 1 ; k <= basisLen ; k++)
  { i = basisVec[k].cndx ;

    while ((i-lastRow) > 1)
    { lastRow++ ;
      aj->coeffs[0].ndx = lastRow ;
      aj->coeffs[0].val = 1.0 ;
      aj->nme = consys_nme(sys,'v',n+lastRow,FALSE,NULL) ;
      aj->cnt = 1 ;
      if (consys_addcol_pk(basis,vartypCON,aj,0.0,0.0,0.0) == FALSE)
      { errmsg(112,rtnnme,
	       basis->nme,"add column","variable",aj->nme,lastRow) ;
	if (aj != NULL) pkvec_free(aj) ;
	if (basisVec != NULL) FREE(basisVec) ;
	consys_free(basis) ;
	return (NULL) ; }
      basis2sys[aj->ndx] = -lastRow ;
#     ifndef DYLP_NDEBUG
      if (printlvl >= 5)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "  fabricated unit column for inactive %s (%d) at %d.\n",
		    aj->nme,lastRow,aj->ndx) ; }
#     endif
    }

    j = basisVec[k].vndx ;
    if (j < 0)
    { aj->coeffs[0].ndx = -j ;
      aj->coeffs[0].val = 1.0 ;
      aj->nme = consys_nme(sys,'v',n-j,FALSE,NULL) ;
      aj->cnt = 1 ; }
    else
    { if (consys_getcol_pk(sys,j,&aj) == FALSE)
      { errmsg(122,rtnnme,sys->nme,"column",
	       consys_nme(sys,'v',j,FALSE,NULL),j) ;
	if (aj != NULL) pkvec_free(ai) ;
	if (basisVec != NULL) FREE(basisVec) ;
	consys_free(basis) ;
	return (NULL) ; } }
    if (consys_addcol_pk(basis,vartypCON,aj,0.0,0.0,0.0) == FALSE)
    { errmsg(112,rtnnme,basis->nme,"add column","variable",aj->nme,abs(j)) ;
      if (aj != NULL) pkvec_free(aj) ;
      if (basisVec != NULL) FREE(basisVec) ;
      consys_free(basis) ;
      return (NULL) ; }
    basis2sys[aj->ndx] = j ;

#   ifndef DYLP_NDEBUG
    if (printlvl >= 5)
    { if (j < 0)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "  fabricated unit column for %s (%d) at %d.\n",
		    aj->nme,-j,aj->ndx) ; }
      else
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "  inserted column for %s (%d) at %d.\n",
		    aj->nme,j,aj->ndx) ; } }
#   endif

    lastRow++ ; }
  if (basisVec != NULL) FREE(basisVec) ;
/*
  There's no guarantee that we've covered all constraints. See if there are
  any left. Those that remain are definitely inactive.
*/
  for (k = lastRow+1 ; k <= m ; k++)
  { aj->coeffs[0].ndx = k ;
    aj->coeffs[0].val = 1.0 ;
    aj->nme = consys_nme(sys,'v',n+k,FALSE,NULL) ;
    aj->cnt = 1 ;
    if (consys_addcol_pk(basis,vartypCON,aj,0.0,0.0,0.0) == FALSE)
    { errmsg(112,rtnnme,
	     basis->nme,"add column","variable",aj->nme,k) ;
      if (aj != NULL) pkvec_free(aj) ;
      consys_free(basis) ;
      return (NULL) ; }
    basis2sys[aj->ndx] = -k ;
#   ifndef DYLP_NDEBUG
    if (printlvl >= 5)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "  fabricated unit column for inactive %s (%d) at %d.\n",
		  aj->nme,k,aj->ndx) ; }
#   endif
  }

  if (aj != NULL) pkvec_free(aj) ;

# ifndef DYLP_NDEBUG
  if (printlvl >= 3)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"  Basis matrix is:\n") ;
    dyio_outfmt(dy_logchn,dy_gtxecho,
		"Pos'n\tConstraint\t  Variable\t  Orig.Col\n") ;
    for (i = 1 ; i <= m ; i++)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "%5d\t%-16s",i,consys_nme(basis,'c',i,FALSE,NULL)) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "  %-16s  %5d\n",consys_nme(basis,'v',i,FALSE,NULL),
		  basis2sys[i]) ; } }
# endif

  return (basis) ; }




int dytest_betaj (lpprob_struct *main_lp, lptols_struct *main_lptols,
		  lpopts_struct *main_lpopts)

/*
  This routine checks the accuracy of the tableau routine dy_betaj (column
  beta<j> of the basis inverse) by testing that B beta<j> = e<j>.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if Binv(B) = I, number of errors otherwise.
*/

{ int m,i,j,k ;
  consys_struct *sys ;

  consys_struct *basis ;

  int *basis2sys ;

  double *betaj ;
  double aidotbetaj,expected ;
  int errcnt ;

  char *rtnnme = "dytest_betaj" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: checking columns of basis inverse using %s (%d x %d).\n",
		rtnnme,sys->nme,sys->concnt,sys->varcnt) ;
    if (main_lpopts->print.tableau >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "  basis contains %d entries.\n",main_lp->basis->len) ; } }
# endif
/*
  Create the basis matrix.
*/
  basis2sys = NULL ;
  basis = create_basis(main_lp,main_lptols,main_lpopts,&basis2sys) ;
  if (basis == NULL || basis2sys == NULL)
  { errmsg(152,rtnnme,"basisMtx") ;
    consys_free(basis) ;
    if (basis2sys != NULL) FREE(basis2sys) ;
    return (1) ; }
/*
  Now that we have a basis matrix with columns in the proper order, matching
  the constraints, we can simply call dy_betaj to obtain columns of the basis
  inverse and call consys_dotrow(i,beta<j>), i = 1, ..., m, and check that we
  have e<j>.
*/
  m = sys->concnt ;
  betaj = NULL ;
  errcnt = 0 ;
  for (k = 1 ; k <= m ; k++)
  { j = basis2sys[k] ;
    if (dy_betaj(main_lp,j,&betaj) == FALSE)
    { errmsg(952,rtnnme,sys->nme,"column",j,"variable",
	     consys_nme(sys,'v',j,FALSE,NULL),j) ;
      errcnt++ ;
      continue ; }
    for (i = 1 ; i <= m ; i++)
    { aidotbetaj = consys_dotrow(basis,i,betaj) ;
      if (i == k)
      { expected = 1.0 ; }
      else
      { expected = 0.0 ; }
      if (fabs(aidotbetaj-expected) > main_lptols->zero)
      { errcnt++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: a<%d> dot beta<%d> = %g ; expected %g; ",
		    i,j,aidotbetaj,expected) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"err %g, tol %g.",
		    (aidotbetaj-expected),main_lptols->zero) ; } } }
/*
  We're done. Do a bit of cleanup.
*/
  if (betaj != 0) FREE(betaj) ;
  consys_free(basis) ;
  if (basis2sys != NULL) FREE(basis2sys) ;

  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: found %d errors testing Binv(B).\n",
	   rtnnme,errcnt) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: pass Binv(B).\n",rtnnme) ; }

  return (errcnt) ; }



int dytest_abarj (lpprob_struct *main_lp, lptols_struct *main_lptols,
		  lpopts_struct *main_lpopts)

/*
  This routine checks the accuracy of the tableau routine dy_abarj, where
  abar<j> = inv(B)a<j>, by testing that B abar<j> = a<j>.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if B(inv(B)A) = A, error count otherwise.
*/

{ int m,n,i,j,k ;
  consys_struct *sys ;

  consys_struct *basis ;

  int *basis2sys ;

  double *abarj,*aj ;
  double aidotabarj ;
  int errcnt ;

  char *rtnnme = "dytest_abarj" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: checking ftran'd columns abar<j> using %s (%d x %d).\n",
		rtnnme,sys->nme,m,n) ;
    if (main_lpopts->print.tableau >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "  basis contains %d entries.\n",main_lp->basis->len) ; } }
# endif
/*
  Create the basis matrix.
*/
  basis2sys = NULL ;
  basis = create_basis(main_lp,main_lptols,main_lpopts,&basis2sys) ;
  if (basis == NULL || basis2sys == NULL)
  { errmsg(152,rtnnme,"basisMtx") ;
    consys_free(basis) ;
    if (basis2sys != NULL) FREE(basis2sys) ;
    return (1) ; }
/*
  Now that we have a basis matrix with columns in the proper order, matching
  the constraints, we can simply call dy_abarj to obtain ftran'd columns
  abar<j> and call consys_dotrow(i,abar<j>), i = 1, ..., m, and check that we
  have a<j>.
*/
  aj = NULL ;
  abarj = NULL ;
  errcnt = 0 ;
  for (j = 1 ; j <= n ; j++)
  { if (dy_abarj(main_lp,j,&abarj) == FALSE)
    { errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	     consys_nme(sys,'v',j,FALSE,NULL),j) ;
      errcnt++ ;
      continue ; }
    if (consys_getcol_ex(sys,j,&aj) == FALSE)
    { errmsg(122,rtnnme,sys->nme,"column",consys_nme(sys,'v',j,FALSE,NULL),j) ;
      errcnt++ ;
      continue ; }
    for (i = 1 ; i <= m ; i++)
    { aidotabarj = consys_dotrow(basis,i,abarj) ;
      if (fabs(aidotabarj-aj[i]) > main_lptols->zero)
      { errcnt++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: a<%d> dot abar<%d> = %g ; expected %g; ",
		    i,j,aidotabarj,aj[i]) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"err %g, tol %g.",
		    (aidotabarj-aj[i]),main_lptols->zero) ; } } }
/*
  And to be really thorough, test the columns associated with logicals.
*/
  memset(aj,0,((size_t) ((m+1)*sizeof(double)))) ;
  for (k = 1 ; k <= m ; k++)
  { if (dy_abarj(main_lp,-k,&abarj) == FALSE)
    { errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	     consys_nme(sys,'v',n+k,FALSE,NULL),k) ;
      errcnt++ ;
      continue ; }
    aj[k] = 1.0 ;
    for (i = 1 ; i <= m ; i++)
    { aidotabarj = consys_dotrow(basis,i,abarj) ;
      if (fabs(aidotabarj-aj[i]) > main_lptols->zero)
      { errcnt++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: a<%d> dot abar<%d> = %g ; expected %g; ",
		    i,-k,aidotabarj,aj[i]) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"err %g, tol %g.",
		    (aidotabarj-aj[i]),main_lptols->zero) ; } }
    aj[k] = 0.0 ; }
/*
  We're done. Do a bit of cleanup.
*/
  if (abarj != 0) FREE(abarj) ;
  if (aj != 0) FREE(aj) ;
  consys_free(basis) ;
  if (basis2sys != NULL) FREE(basis2sys) ;

  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: found %d errors testing B(inv(B)A) = A.\n",
	   rtnnme,errcnt) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: pass B(inv(B)A).\n",rtnnme) ; }
  
  return (errcnt) ; }



int dytest_betai (lpprob_struct *main_lp, lptols_struct *main_lptols,
		  lpopts_struct *main_lpopts)

/*
  This routine checks the accuracy of the tableau routine dy_betai (row
  beta<i> of the basis inverse) by testing that beta<i> B = e<i>.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if inv(B)B = I, error count otherwise.
*/

{ int m,i,j ;
  consys_struct *sys ;

  consys_struct *basis ;

  int *basis2sys ;

  double *betai ;
  double betaidotaj,expected ;
  int errcnt ;

  char *rtnnme = "dytest_betai" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: checking rows of basis inverse using %s (%d x %d).\n",
		rtnnme,sys->nme,sys->concnt,sys->varcnt) ;
    if (main_lpopts->print.tableau >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "  basis contains %d entries.\n",main_lp->basis->len) ; } }
# endif
/*
  Create the basis matrix.
*/
  basis2sys = NULL ;
  basis = create_basis(main_lp,main_lptols,main_lpopts,&basis2sys) ;
  if (basis == NULL || basis2sys == NULL)
  { errmsg(152,rtnnme,"basisMtx") ;
    consys_free(basis) ;
    if (basis2sys != NULL) FREE(basis2sys) ;
    return (1) ; }
/*
  Now that we have a basis matrix with columns in the proper order, matching
  the constraints, we can simply call dy_betai to obtain rows of the basis
  inverse and call consys_dotcol(j,beta<i>), j = 1, ..., m, and check that we
  have e<i>.
*/
  m = sys->concnt ;
  betai = NULL ;
  errcnt = 0 ;
  for (i = 1 ; i <= m ; i++)
  { if (dy_betai(main_lp,i,&betai) == FALSE)
    { errmsg(952,rtnnme,sys->nme,"row",i,"constraint",
	     consys_nme(sys,'c',i,FALSE,NULL),i) ;
      errcnt++ ;
      continue ; }
    for (j = 1 ; j <= m ; j++)
    { betaidotaj = consys_dotcol(basis,j,betai) ;
      if (i == j)
      { expected = 1.0 ; }
      else
      { expected = 0.0 ; }
      if (fabs(betaidotaj-expected) > main_lptols->zero)
      { errcnt++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: beta<%d> dot a<%d> = %g ; expected %g; ",
		    i,j,betaidotaj,expected) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"err %g, tol %g.",
		    (betaidotaj-expected),main_lptols->zero) ; } } }
/*
  We're done. Do a bit of cleanup.
*/
  if (betai != 0) FREE(betai) ;
  consys_free(basis) ;
  if (basis2sys != NULL) FREE(basis2sys) ;

  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: found %d errors testing inv(B)B.\n",
	   rtnnme,errcnt) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s: pass inv(B)B.\n",rtnnme) ; }
  
  return (errcnt) ; }


int dytest_abari (lpprob_struct *main_lp, lptols_struct *main_lptols,
		  lpopts_struct *main_lpopts)

/*
  This routine checks the accuracy of the tableau routine dy_abari, which
  calculates row i of inv(B)[ B N I ], where B is the basic columns, N the
  nonbasic columns, and I is the identity matrix produced by the coefficients
  of logical variables.

  This is an inconvenient calculation to check --- we can't premultiply the
  resulting row by the basis, as we do for all the other routines. So we do the
  next best thing: call dy_abarj and check that abar<ij> matches. It's
  expensive, but hey, this is a test routine. This does imply that dy_abarj
  should be tested first.

  Parameters:
    main_lp:	 the lp problem structure
    main_lptols: the lp tolerance structure
    main_lpopts: the lp options structure

  Returns: 0 if inv(B)B = I, error count otherwise.
*/

{ int m,n,i,j ;
  consys_struct *sys ;

  double *abari,*abarj,*betai ;
  double abarij,expected ;
  int errcnt ;

  char *rtnnme = "dytest_abari" ;

/*
  Do a little initialisation. Mention that we've started.
*/
  sys = main_lp->consys ;
  m = sys->concnt ;
  n = sys->varcnt ;

# ifndef DYLP_NDEBUG
  if (main_lpopts->print.tableau >= 1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"%s: checking rows of inv(B)A using %s (%d x %d).\n",
		rtnnme,sys->nme,m,n) ; }
# endif

/*
  Open a pair of loops to do the testing.
*/
  errcnt = 0 ;
  abari = NULL ;
  betai = NULL ;
  abarj = NULL ;
  for (i = 1 ; i <= m ; i++)
  { if (dy_abari(main_lp,i,&abari,&betai) == FALSE)
    { errmsg(953,rtnnme,sys->nme,"transformed","row",
	     consys_nme(sys,'c',i,FALSE,NULL),i) ;
      errcnt++ ;
      continue ; }
    for (j = 1 ; j <= n ; j++)
    { if (dy_abarj(main_lp,j,&abarj) == FALSE)
      { errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	       consys_nme(sys,'v',j,FALSE,NULL),j) ;
	errcnt++ ;
	continue ; }
      expected = abarj[i] ;
      abarij = abari[j] ;
      if (fabs(abarij-expected) > main_lptols->zero)
      { errcnt++ ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n  ERROR: beta<%d> dot a<%d> = %g ; expected %g; ",
		    i,j,abarij,expected) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,"err %g, tol %g.",
		    (abarij-expected),main_lptols->zero) ; } }
/*
  Now test the columns for the logical variables.
*/
  for (j = 1 ; j <= m ; j++)
  { if (dy_abarj(main_lp,-j,&abarj) == FALSE)
    { errmsg(953,rtnnme,sys->nme,"ftran'd","column",
	     consys_nme(sys,'v',n+j,FALSE,NULL),j) ;
      errcnt++ ;
      continue ; }
    expected = abarj[i] ;
    abarij = betai[j] ;
    if (fabs(abarij-expected) > main_lptols->zero)
    { errcnt++ ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  ERROR: beta<%d> dot a<%d> = %g ; expected %g; ",
		  i,-j,abarij,expected) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"err %g, tol %g.",
		  (abarij-expected),main_lptols->zero) ; } } }
/*
  We're done. Do a bit of cleanup.
*/
  if (abari != 0) FREE(abari) ;
  if (betai != 0) FREE(betai) ;
  if (abarj != 0) FREE(abarj) ;

  if (errcnt != 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	  "\n%s: found %d errors testing e<i>(inv(B)A) against inv(B)a<j>.\n",
	   rtnnme,errcnt) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n%s: pass e<i>(inv(B)A) against inv(B)a<j>.\n",rtnnme) ; }

  return (errcnt) ; }


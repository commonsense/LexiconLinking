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
  This file contains the interface layer between dylp and the glpk routines
  which maintain an LU factorization of the basis inverse.

  One of the nicer aspects of using glpk is that it uses the same convention
  for indexing as dylp: vector[0] is unused, and entries start with vector[1].
*/

#define DYLP_INTERNAL

#include "dylp.h"
#include "glpinv.h"

static char sccsid[] UNUSED = "@(#)dy_basis.c	4.6	11/06/04" ;
static char svnid[] UNUSED = "$Id: dy_basis.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  The basis data structure. glpinv/glpluf does not distinguish between
  allocated capacity and actual basis size. So we cheat --- we allocate the
  size we want, to get vectors of the proper capacity, then reset the size
  fields so that glpinv/glpluf use only what's needed.
*/

static INV *luf_basis = NULL ;
static int luf_capacity = 0 ;

/*
  Basis patch structure.

  This structure is used to record changes to the basis when a singular basis
  is repaired. It's filled in by adjust_basis and then used by adjust_therest.

  Field		Definition
  -----		----------
  pos		the basis position that has changed
  out		the variable removed from the basis
  in		the variable inserted into the basis
*/

typedef struct
{ int pos ;
  int out ;
  int in ; } patch_struct ;



/*
  Pivot selection during factoring

  glpinv/glpluf uses dynamic Markowitz pivot selection with partial (PTPS)
  threshold pivot selection. A pivot for row i must satisfy
  piv_tol*max{j}|a<i,j>|.

  There are seven levels of pivot selection, as follows:

    stable	look
    ------	-------
       .01	  4
       .05	  4
       .1	  4
       .2	  4
       .4	  6
       .8	  8
       .95	  10

  The values for stable are based on guidelines in Suhl & Suhl, Computing
  Sparse LU Factorizations for Large-Scale Linear Programming Bases, ORSA
  Journal on Computing, v2(4), Fall, 1990. glpluf uses a default threshold
  of .1, and tightens to .3, then .7, before giving up. We're moving further
  at both extremes. Note that 1.0 is illegal, according to glpluf.

  glpinv/glpluf also provides for a heuristic limit on the number of pivot
  candidates considered. After k candidates have been examined, it returns
  the best so far. This parameter will affect sparsity (numerical stability
  must be satisfied by all pivots). As we get fussier about numerical
  stability, make glpluf look a bit harder for sparsity.

  What's implemented in practice is a table with the above values, and a
  subroutine which takes a parameter, delta, and raises (delta > 0) or lowers
  (delta < 0) the pivot selection parameters delta steps. Any large value of
  delta serves to slam the tolerances to one extreme or the other. There's
  also a static value, minpivlevel, which defaults to 0 but can be raised to
  deal with persistent accuracy problems.

  The routine dy_factor can tighten the pivot selection parameters if an
  attempt to factor the basis is aborted due to numerical instability. It
  does not loosen them, however --- this decision is left to higher levels.
  (In dylp, dy_accchk can raise and lower the current pivot selection
  parameters. Unexpected loss of accuracy or feasibility during preoptimality
  checks results in a boost to the minimum pivot level which remains in
  effect until the start of a new LP phase.)

  pivtols_struct

  Field		Description
  -----		-----------
  stable	The pivot selection ratio.
  look		The number of candidates considered.
*/

typedef struct { double stable ;
		 int look ; } pivtols_struct ;

static pivtols_struct pivtols[] = {{.01,4},
				   {.05,4},
				   {.1,4},
				   {.2,4},
				   {.4,6},
				   {.8,8},
				   {.95,10}} ;
static int pivlevel,minpivlevel ;

#define MAX_PIVLEVEL (((sizeof pivtols)/(sizeof(pivtols_struct)))-1)



char *dy_prtpivparms (int lvl)

/*
  Utility routine to construct a print string for a set of pivoting
  parameters. It exists just to keep things local to this file.

  Parameters:
    lvl: pivot level to print; if < 0, the current pivot level (as set in
	 luf_basis) is printed
  Returns: pointer to a static character string
*/

{ static char buffer[20] ;
  pivtols_struct pivtol ;

# ifdef PARANOIA
  const char *rtnnme = "dy_prtpivparms" ;

  if (luf_basis == NULL)
  { errmsg(2,rtnnme,"luf_basis") ;
    return ("<<error: no basis!>>") ; }
# endif

  if (lvl < 0 || lvl > MAX_PIVLEVEL)
  { pivtol.stable = luf_basis->luf->piv_tol ;
    pivtol.look = luf_basis->luf->piv_lim ; }
  else
  { pivtol.stable = pivtols[lvl].stable ;
    pivtol.look = pivtols[lvl].look ; }

  dyio_outfxd(buffer,-19,'l',"PTPS(%.2f,%d)",pivtol.stable,pivtol.look) ;

  return (buffer) ; }


bool dy_setpivparms (int curdelta, int mindelta)

/*
  This routine exists to allow other parts of the lp code to adjust the pivot
  regimen used when the basis is factored.  curdelta is the change in the
  current pivot level, mindelta the change in the minimum pivot level.  A
  positive delta increases (tightens) the tolerances by delta steps, and a
  negative delta lowers (loosens) them.

  If either delta is 0, there's no change. A large positive or negative
  value slams the tolerance to the appropriate extreme.


  Parameters:
    curdelta:	change to the current pivot level
    mindelta:	change to the minimum pivot level

  Returns: TRUE if any change actually occurred, FALSE if we were already
	   at the relevant limit(s) (BE CAREFUL interpreting the return value
	   if you're changing both at once)
*/

{ bool minretval,curretval ;

# ifdef PARANOIA
  const char *rtnnme = "dy_setpivparms" ;

  if (luf_basis == NULL)
  { errmsg(2,rtnnme,"luf_basis") ;
    return (FALSE) ; }
# endif

  minretval = FALSE ;
  curretval = FALSE ;
/*
  Adjust the minimum pivot level first. This may imply an adjustment in the
  current pivot level.
*/
  if (mindelta != 0)
  { if ((minpivlevel <= 0 && mindelta < 0) ||
	(minpivlevel >= MAX_PIVLEVEL && mindelta > 0))
    { 
#     ifndef DYLP_NDEBUG
      if ((dy_opts->print.basis >= 3) ||
          (dy_opts->print.basis >= 2 && mindelta > 0))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t    min. pivot ratio unchanged at %s (%d)",
		    dy_prtpivparms(minpivlevel),minpivlevel) ; }
#     endif
    }
    else
    { minretval = TRUE ;
      minpivlevel += mindelta ;
      if (minpivlevel < 0)
	minpivlevel = 0 ;
      else
      if (minpivlevel > MAX_PIVLEVEL)
	minpivlevel = MAX_PIVLEVEL ;
      if (pivlevel < minpivlevel) 
	curdelta = maxx(curdelta,(minpivlevel-pivlevel)) ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.basis >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t    setting min. pivot ratio to %s (%d)",
		    dy_prtpivparms(minpivlevel),minpivlevel) ; }
#     endif
    } }
/*
  Adjust the current pivot level.
*/
  if (curdelta != 0)
  { if ((pivlevel <= minpivlevel && curdelta < 0) ||
	(pivlevel >= MAX_PIVLEVEL && curdelta > 0))
    { 
#     ifndef DYLP_NDEBUG
      if ((dy_opts->print.basis >= 3) ||
          (dy_opts->print.basis >= 2 && mindelta > 0))
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t    cur. pivot ratio unchanged at %s (%d)",
		    dy_prtpivparms(-1),pivlevel) ; }
#     endif
    }
    else
    { curretval = TRUE ;
      pivlevel += curdelta ;
      if (pivlevel < minpivlevel)
	pivlevel = minpivlevel ;
      else
      if (pivlevel > MAX_PIVLEVEL)
	pivlevel = MAX_PIVLEVEL ;
      luf_basis->luf->piv_tol = pivtols[pivlevel].stable ;
      luf_basis->luf->piv_lim =  pivtols[pivlevel].look ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.basis >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\t    setting cur. pivot ratio to %s (%d)",
		    dy_prtpivparms(-1),pivlevel) ; }
#     endif
    } }

  if (curretval == FALSE && minretval == FALSE)
    return (FALSE) ;
  else
    return (TRUE) ; }


double dy_chkpiv (double abarij, double maxabar)

/*
  This routine checks that the pivot element satisfies a stability test. The
  stability test is
    |abarij/maxabar| > dy_tols->pivot*piv_tol

  The motivation for the check is to try and choose reasonably stable pivots
  --- inv_update will do what we tell it, might as well try to choose half
  wisely. Balance that against the desire to make progress toward the
  optimum. Setting dy_tols->pivot to 1 would give the same selection criteria
  as is used during factorization.

  Parameters:
    abarij:	the pivot element (only the absolute value is used)
    maxabar:	(primal) max{i} |abar<i,j>|
		(dual)   max{j} |abar<i,j>|

  Returns: |abarij/(dy_tols->pivot*piv_tol*maxabar)|
	   return value < 1.0 indicates an unstable pivot
*/

{ double ratio,abspiv,stable ;
  
# if defined(PARANOIA) || !defined(DYLP_NDEBUG)
  const char *rtnnme = "dy_chkpiv" ;
# endif

# ifdef PARANOIA
  if (luf_basis == NULL)
  { errmsg(2,rtnnme,"luf_basis") ;
    return (dyrFATAL) ; }
# endif

/*
  In some instances (for example, netlib/grow22), we run into trouble because
  column coefficients inflate to outrageous values --- 1.0e18, for example.
  Take the attitude that any pivot s.t. |abar<ij> > 1.0| should be accepted,
  to give a bit more choice in pivot selection. To do this, lie about the
  stability ratio,
*/
  abspiv = fabs(abarij) ;
  ratio = dy_tols->pivot*luf_basis->luf->piv_tol ;
  stable = ratio*maxabar ;
# ifndef DYLP_NDEBUG
  if (dy_opts->print.pivoting >= 1)
  { if (abspiv/stable < 1.0)
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n%s: %s pivot = %g < %g; column max = %g, ratio = %g.",
		  rtnnme,(abspiv < 1.0)?"rejecting":"tolerating",
		  abarij,stable,maxabar,ratio) ; }
# endif
  
  if (abspiv/stable >= 1.0)
    return (abspiv/stable) ;
  else
  if (abspiv >= 1.0)
    return (1.0) ;
  else
    return (abspiv/stable) ; }




void dy_initbasis (int concnt, int factor, double zero_tol)

/*
  This routine calls the glpk routine inv_create to initialize the basis
  data structures, then sets values for the zero tolerance (eps_tol), pivot
  ratio (piv_tol) and number of candidates examined (piv_lim).

  NOTE: This routine can be (and typically is) called before any of the main
	dylp data structures exist. Be careful what you reference.

  Parameters:
    concnt:	the number of constraints (rows) that the basis representation
		should be capable of handling
    factor:	the planned refactorisation frequency; passed to glpk as the
		basis inverse update capacity (i.e., the limit on the number
		of pivots between refactorisations)
    zero_tol:	zero tolerance; a value of 0.0 uses the glpk default
		(INV->LUF->eps_tol = 1.0e-15).

  Returns: void
*/

{ int sva_size ;
  const char *rtnnme = "dy_initbasis" ;

/*
  Create the basis.
*/
  luf_capacity = concnt ;
  luf_basis = inv_create(luf_capacity,factor) ;
  if (luf_basis == NULL)
  { if (dy_lp == NULL)
    { errmsg(302,rtnnme,"empty","pre-init",0,"create") ; }
    else
    { errmsg(302,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters,"create") ; }
    return ; }
/*
  WARNING: We're going to reach inside glpluf to get it to triple the amount
  of space that it allocates for the sparse vector area. We're doing this by
  triggering the reallocation mechanism built into luf_decomp (called by
  inv_decomp).
*/
  sva_size = luf_basis->luf->sv_size ;
  luf_basis->luf->new_sva = 3*sva_size ;
# ifndef DYLP_NDEBUG
  if (dy_opts != NULL && dy_opts->print.basis >= 2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
	        "\ninitbasis: %s(%d) basis capacity %d, piv lim %d.",
	        dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
	        luf_basis->luf->n,luf_basis->hh_max) ; }
/*
  XX_DEBUG_XX

  There's no good way to control this output, given the timing of the call
  (before dy_opts is initialised), but it's sometimes useful when debugging.
  
  else
  { dyio_outfmt(dy_logchn,TRUE,
	        "\ninitbasis: EXTERN(0) basis capacity %d, piv lim %d.",
	        luf_basis->luf->n,luf_basis->hh_max) ; }
*/
# endif
/*
  Set the initial pivot level to {.01,4}, and allow it to drop to {.01,4}.
*/

  pivlevel = 0 ;
  minpivlevel = 0 ;

  if (zero_tol != 0.0) luf_basis->luf->eps_tol = zero_tol ;
  luf_basis->luf->piv_tol = pivtols[pivlevel].stable ;
  luf_basis->luf->piv_lim =  pivtols[pivlevel].look ;
  luf_basis->luf->max_gro = 1.0e7 ;
/*
  This is the smallest value that can appear on the diagonal of U after a
  pivot update. dylp will (in extremis) drop its pivot selection tolerance
  tols.pivot to 1e-9 (or thereabouts), so upd_tol had better be less or we
  spend a lot of time refactoring. This should probably be adjusted as
  needed, in response to adjustments in tols.pivot, but I need to sit down and
  really think about the math. In the meantime, this seems to be adequate.
*/
  luf_basis->upd_tol = 1.0e-10 ;

  return ; }


void dy_freebasis (void)

/*
   This routine calls inv_delete to clean up the space allocated for the
   basis representation, then free_lib_env to clean up the library.
  
   Parameters: none
   Returns: undefined
*/

{

  /* glplib2.c */
  extern int _glp_free_lib_env(void) ;


  if (luf_basis != NULL)
  { inv_delete(luf_basis) ;
    luf_basis = NULL ; }
  
  (void) _glp_free_lib_env() ;

  return ; }


static void luf_adjustsize (void)
/*
  This routine sets the proper working size in the basis. It also checks that
  the allocated capacity of the basis is sufficient for the current size of
  the constraint system and increases it if necessary.

  Parameters: none

  Returns: undefined
*/

{ double upd_tol,piv_tol,zero_tol,max_gro ;
  int factor,look,capacity ;
# ifndef DYLP_NDEBUG
  int oldcapacity ;
# endif
/*
  Check that we're within the allocated limits, and resize if necessary.
  Preserve the tolerance settings that we may have tweaked.
*/
  if (dy_sys->concnt > luf_capacity)
  { factor = luf_basis->hh_max ;
    upd_tol = luf_basis->upd_tol ;
    zero_tol = luf_basis->luf->eps_tol ;
    piv_tol = luf_basis->luf->piv_tol ;
    look = luf_basis->luf->piv_lim ;
    max_gro = luf_basis->luf->max_gro ;
#   ifndef DYLP_NDEBUG
    oldcapacity = luf_basis->m ;
#   endif
    dy_freebasis() ;
    capacity = (int) (dy_sys->concnt*1.5) ;
    dy_initbasis(capacity,factor,zero_tol) ;
    luf_basis->upd_tol = upd_tol ;
    luf_basis->luf->eps_tol = zero_tol ;
    luf_basis->luf->piv_tol = piv_tol ;
    luf_basis->luf->piv_lim = look ;
    luf_basis->luf->max_gro = max_gro ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.basis >= 2)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    increased basis capacity from %d to %d constraints",
		  oldcapacity,luf_basis->m) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,", piv lim %d.",luf_basis->hh_max) ; }
#   endif
  }
/*
  Now set the working size.
*/
  luf_basis->m = dy_sys->concnt ;
  luf_basis->luf->n = dy_sys->concnt ;

  return ; }




void dy_ftran (double *col, bool save)

/*
  This is a shell that calls inv_ftran to perform inv(B)*col. We need to set
  clean zeros before returning the result.

  Parameter:
    col:	The column vector to be ftran'ed.
    save:	This will be the entering column in the next pivot, and
		the basis maintenance code should remember it (as the implicit
		parameter to the pivot).

  Returns: undefined
*/

{ int isave ;

  if (save == TRUE)
    isave = 1 ;
  else
    isave = 0 ;

  inv_ftran(luf_basis,col,isave) ;

  return ; }


void dy_btran (double *col)

/*
  This is a shell that calls inv_btran to perform col*inv(B). We need to set
  clean zeros before returning the result.

  Parameter:
    col:	The column vector to be btran'ed.

  Returns: undefined
*/

{ inv_btran(luf_basis,col) ;

  return ; }



static void adjust_basis (int *p_patchcnt, patch_struct **p_patches)

/*
  This routine corrects the dylp basis arrays when glpinv/glpluf declares
  the current basis to be singular. glpluf doesn't actually salvage the
  basis --- it just reports the linearly dependent (hence unpivoted) columns
  and corresponding unpivoted rows. Once we've adjusted the basis accordingly,
  we can make another attempt to factor.

  The convention is as follows: luf_basis.rank gives the rank of the basis.
  qq_col[rank+1 .. m] contain the indices of the basic columns that must
  be removed from the basis. pp_row[rank+1 .. m] contain the indices of the
  basic rows that remain unpivoted; we'll put the logicals for these rows into
  the basis.

  Both of the above are expressed in terms of basis positions. For the rows,
  basis position i is equivalent to constraint i; for the columns, we need
  to look up the variable j in basis position i.

  Recognise that in general this is the easiest part of salvaging the
  situation. We record the changes for adjust_therest, which will do the
  remainder of the work after the basis is successfully factored.

  Parameters:
    p_patchcnt:	(o) the number of basis corrections
    p_patches:	(o) patch array recording the basis corrections

  Returns: undefined
*/

{ int *qq_col, *pp_row ;
  int rank,pqndx,i,j,k,pndx ;
  patch_struct *patches ;


#ifdef PARANOIA
  const char *rtnnme = "adjust_basis" ;

  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return ; }
  if (dy_basis == NULL)
  { errmsg(2,rtnnme,"basis") ;
    return ; }
  if (dy_var2basis == NULL)
  { errmsg(2,rtnnme,"var2basis") ;
    return ; }
  if (luf_basis == NULL)
  { errmsg(2,rtnnme,"LUF basis") ;
    return ; }
  if (p_patches == NULL)
  { errmsg(2,rtnnme,"p_patches") ;
    return ; }
#endif

  qq_col = luf_basis->luf->qq_col ;
  pp_row = luf_basis->luf->pp_row ;
  rank = luf_basis->luf->rank ;

  patches =
    (patch_struct *) MALLOC((dy_sys->concnt-rank)*sizeof(patch_struct)) ;

/*
  Walk qq_col, retrieving the basis position that must be corrected. Remove
  the corresponding variable from the basis, and put in its place the slack
  for the basis row.
*/
  for (pqndx = rank+1, pndx = 0 ; pqndx <= dy_sys->concnt ; pqndx++, pndx++)
  { k = qq_col[pqndx] ;
    j = dy_basis[k] ;
    i = pp_row[pqndx] ;
    dy_basis[k] = i ;
    dy_var2basis[j] = 0 ;
    dy_var2basis[i] = k ;
    patches[pndx].pos = k ;
    patches[pndx].out = j ;
    patches[pndx].in = i ;
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.basis >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      pos'n %d (%s (%d)) replacing %s (%d) with %s (%d).",
		  k,consys_nme(dy_sys,'c',k,FALSE,NULL),k,
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		  consys_nme(dy_sys,'v',i,FALSE,NULL),i) ; }
#   endif
  }

  *p_patchcnt = pndx ;
  *p_patches = patches ;
  return ; }



static dyret_enum adjust_therest (int patchcnt, patch_struct *patches)

/*
  We're here because we've successfully patched a singular basis. The patches
  array contains entries of the form <basis pos'n, x<j>, x<i>>, where x<j>
  has just been kicked out of the basis and replaced by x<i>.  The basis and
  var2basis vectors are already corrected (we needed them to complete the
  factorization).  Now we need to adjust other dylp data structures to
  reflect the unexpected change.  The amount of additional work to be done
  depends on the phase of the simplex algorithm.

    dyINIT: We're done. We've just factored the initial basis and none of the
	    other data structures have been initialised. We didn't really
	    need this call, but the code is cleaner this way.

  If we're farther along, we might be in the middle of simplex (dyPRIMAL1,
  dyPRIMAL2, or dyDUAL), or we might be manipulating the constraint system.

  If we're running simplex, the first actions are cleanup: clear the pivot
  reject list and back out any antidegeneracy activity.
    
  Next, set the status of the newly nonbasic variables, consistent with their
  previous status. The general rule is to perturb the solution as little as
  possible. If we're in a primal or dual simplex phase, try to make decisions
  that are compatible with primal or dual feasibility. Two specific points:
    * Superbasic (SB) variables are only created in dyPRIMAL2.
    * Nonbasic free (NBFR) variables imply loss of dual feasibility.

  Once we have nonbasic status set, we can calculate new primals, duals, and
  reduced costs and fine-tune the status of the newly basic variables. If
  we've arrived here from one of the constraint system manipulation phases,
  there will almost certainly be duplication of effort once we return. But
  hey, how often does a basis patch happen, anyway?
  
  If we're in a simplex phase, there's still some work to do to make the
  patch as transparent as possible.

  For dual simplex, we'll check the status of the nonbasic variables and
  try to maintain dual feasibility. This may not be possible. If we do
  maintain dual feasibility, reset the DSE norms.

  For primal simplex, we need to reset the PSE norms.

  Parameters:
    patchcnt:	the number of basis changes
    patches:	array of basis changes

  Returns: dyrOK if the repair proceeds without error, dyrLOSTDFEAS if
	   feasibility is lost in dual phase II, and dyrFATAL if anything
	   else goes wrong.
*/

{ int i,j,pndx ;
  pkvec_struct *aj ;
  flags statj ;
  dyret_enum retval ;
  dyphase_enum phase ;
  double valj,cbarj,*vub,*vlb,*obj ;

  const char *rtnnme = "adjust_therest" ;

# ifndef DYLP_NDEBUG
  flags stati ;
  double vali ;
# endif

# ifdef PARANOIA
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (dyrFATAL) ; }
  if (dy_basis == NULL)
  { errmsg(2,rtnnme,"basis") ;
    return (dyrFATAL) ; }
  if (dy_var2basis == NULL)
  { errmsg(2,rtnnme,"var2basis") ;
    return (dyrFATAL) ; }
  if (patches == NULL)
  { errmsg(2,rtnnme,"patch") ;
    return (dyrFATAL) ; }
# endif

  phase = dy_lp->phase ;

# ifdef PARANOIA
  if (!(phase == dyINIT || phase == dyADDVAR || phase == dyADDCON ||
	phase == dyPRIMAL1 || phase == dyPRIMAL2 || phase == dyDUAL ||
	phase == dyFORCEPRIMAL || phase == dyFORCEDUAL))
  { errmsg(1,rtnnme,__LINE__) ;
    return (dyrFATAL) ; }
  if (!(phase == dyINIT))
  { if (dy_status == NULL)
    { errmsg(2,rtnnme,"status") ;
      return (dyrFATAL) ; }
    if (dy_x == NULL)
    { errmsg(2,rtnnme,"x") ;
      return (dyrFATAL) ; }
    if (dy_xbasic == NULL)
    { errmsg(2,rtnnme,"x<B>") ;
      return (dyrFATAL) ; } }
#endif

  if (phase == dyINIT) return (dyrOK) ;

  vlb = dy_sys->vlb ;
  vub = dy_sys->vub ;
  obj = dy_sys->obj ;
  aj = NULL ;
  retval = dyrOK ;

/*
  If we're in one of the simplex phases, back out any antidegeneracy activity
  and clear the pivot rejection list.  It's easiest to clear the pivot reject
  list ahead of the status modifications so that we don't have to worry about
  the NOPIVOT qualifier when checking status values.
*/
  if (phase == dyPRIMAL1 || phase == dyPRIMAL2 || phase == dyDUAL)
  { if (dy_clrpivrej(NULL) != TRUE) return (dyrFATAL) ;
    if (dy_lp->degen > 0)
    { if (phase == dyDUAL)
      { (void) dy_dualdegenout(0) ; }
      else
      { (void) dy_degenout(0) ; } } }
/*
  Now correct the status for newly nonbasic variables. We need to correct
  dy_x if the status change forces a change in value.  If we end up with a
  NBFR variable, we've lost dual feasibility.

  While we're walking the patches, set the status for x<i> (the newly basic
  variable) to vstatB. No need to be more precise at this point.
*/

  for (pndx = 0 ; pndx < patchcnt ; pndx++)
  { i = patches[pndx].in ;
#   ifndef DYLP_NDEBUG
    stati = dy_status[i] ;
    vali = dy_x[i] ;
#   endif
    dy_status[i] = vstatB ;
    j = patches[pndx].out ;
    statj = dy_status[j] ;
    valj = dy_x[j] ;
    switch (statj)
    { case vstatBLLB:
      { dy_status[j] = vstatNBLB ;
	dy_x[j] = vlb[j] ;
	break ; }
      case vstatBLB:
      { dy_status[j] = vstatNBLB ;
	break ; }
      case vstatB:
      { if (phase == dyPRIMAL2)
	  dy_status[j] = vstatSB ;
	else
	if (valj-vlb[j] < vub[j]-valj)
	{ dy_status[j] = vstatNBLB ;
	  dy_x[j] = vlb[j] ; }
	else
	{ dy_status[j] = vstatNBUB ;
	  dy_x[j] = vub[j] ; }
	break ; }
      case vstatBUB:
      { dy_status[j] = vstatNBUB ;
	break ; }
      case vstatBUUB:
      { dy_status[j] = vstatNBUB ;
	dy_x[j] = vub[j] ;
	break ; }
      case vstatBFX:
      { dy_status[j] = vstatNBFX ;
	break ; }
      case vstatBFR:
      { dy_status[j] = vstatNBFR ;
	if (phase == dyDUAL)
	{ 
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.dual >= 1)
	  { warn(346,rtnnme,
		 dy_sys->nme,dy_prtlpphase(phase,TRUE),dy_lp->tot.iters+1,
		 dy_prtvstat(statj),consys_nme(dy_sys,'v',j,FALSE,NULL),j) ; }
#	  endif
	  retval = dyrLOSTDFEAS ; }
	break ; }
      default:
      { errmsg(380,rtnnme,dy_sys->nme,consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	       dy_prtvstat(statj),"basic") ;
	return (dyrFATAL) ; } }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.basis >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t%s (%d) had status %s, value %g, ",
		  consys_nme(dy_sys,'v',i,FALSE,NULL),i,
		  dy_prtvstat(stati),vali) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"now status %s.",
		  dy_prtvstat(dy_status[i])) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t%s (%d) had status %s, value %g, ",
		  consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		  dy_prtvstat(statj),valj) ;
      dyio_outfmt(dy_logchn,dy_gtxecho,"now status %s, value %g.",
		  dy_prtvstat(dy_status[j]),dy_x[j]) ; }
#   endif
  }

# ifdef PARANOIA
/*
  If paranoid checks are in place, we need agreement between dy_status, dy_x,
  and dy_xbasic, lest dy_calccbar fail. Call dy_calcprimals and
  dy_setbasicstatus to get the basic status right. This is restricted to
  paranoid mode because the proper place to do this is after making
  corrections to nonbasic status for dual feasibility.
*/
  if (dy_calcprimals() == FALSE) return (dyrFATAL) ;
  dy_setbasicstatus() ;
# endif

/*
  Calculate the duals and reduced costs.
*/
  dy_calcduals() ;
  if (dy_calccbar() == FALSE)
  { errmsg(384,rtnnme,
	   dy_sys->nme,dy_prtlpphase(phase,TRUE),dy_lp->tot.iters) ;
    return (dyrFATAL) ; }

/*
  If we're in phase dyDUAL, it's worth a scan to check dual feasibility and
  make adjustments to maintain it, if possible. (retval = dyrLOSTDFEAS says
  we introduced a NBFR variable, in which case we have no hope).

  Open a loop to scan the nonbasic variables. NBFX variables are always dual
  feasible, NBFR variables are never dual feasible.  We're minimising, so
  dual feasibility (primal optimality) is cbarj < 0 && x<j> at upper bound,
  or cbarj > 0 && x<j> at lower bound.  It's important that the zero
  tolerance for cbar<j> here be the same as the one used in dy_dualin when it
  checks for loss of dual feasibility.
*/
  if (phase == dyDUAL && retval != dyrLOSTDFEAS)
  { for (j = 1 ; j <= dy_sys->varcnt ; j++)
    { statj = dy_status[j] ;
      if (flgon(statj,vstatBASIC|vstatNBFX)) continue ;
      if (flgon(statj,vstatNBFR))
      { retval = dyrLOSTDFEAS ;
#	ifndef DYLP_NDEBUG
	cbarj = dy_cbar[j] ;
	if (dy_opts->print.dual >= 1)
	{ warn(347,rtnnme,
	       dy_sys->nme,dy_prtlpphase(phase,TRUE),dy_lp->tot.iters+1,
	       consys_nme(dy_sys,'v',j,FALSE,NULL),j,
	       dy_prtvstat(statj),j,cbarj,dy_tols->dfeas) ; }
#	endif
	break ; }
      cbarj = dy_cbar[j] ;
      if (cbarj < -dy_tols->dfeas && flgoff(statj,vstatNBUB))
      { if (vub[j] >= dy_tols->inf)
	{
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.dual >= 1)
	  { warn(347,rtnnme,
		 dy_sys->nme,dy_prtlpphase(phase,TRUE),dy_lp->tot.iters+1,
		 consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		 dy_prtvstat(statj),j,cbarj,dy_tols->dfeas) ; }
#	  endif
	  retval = dyrLOSTDFEAS ;
	  break ; }
	else
	{ dy_status[j] = vstatNBUB ;
	  dy_x[j] = vub[j] ; } }
      else
      if (cbarj > dy_tols->dfeas && flgoff(statj,vstatNBLB))
      { if (vlb[j] >= dy_tols->inf)
	{
#	  ifndef DYLP_NDEBUG
	  if (dy_opts->print.dual >= 1)
	  { warn(347,rtnnme,
		 dy_sys->nme,dy_prtlpphase(phase,TRUE),dy_lp->tot.iters+1,
		 consys_nme(dy_sys,'v',j,FALSE,NULL),j,
		 dy_prtvstat(statj),j,cbarj,dy_tols->dfeas) ; }
#	  endif
	  retval = dyrLOSTDFEAS ;
	  break ; }
	else
	{ dy_status[j] = vstatNBLB ;
	  dy_x[j] = vlb[j] ; } }
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.basis >= 3 && dy_status[j] != statj)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n\tchanged status of %s (%d) from %s to",
		    consys_nme(dy_sys,'v',j,FALSE,NULL),j,dy_prtvstat(statj)) ;
	dyio_outfmt(dy_logchn,dy_gtxecho,
		    " %s to maintain dual feasibility; cbar = %g.",
		    dy_prtvstat(dy_status[j]),cbarj) ; }
#     endif
    } }

/*
  The dual variables and reduced costs have been recalculated, and we have
  the final status for all nonbasic variables.  Recalculate the primal
  variables and set the status of the basic variables.
*/
  if (dy_calcprimals() == FALSE) return (dyrFATAL) ;
  dy_setbasicstatus() ;
/*
  If we're running primal simplex, reset the PSE reference frame. If we're
  running dual simplex and haven't lost dual feasibility, recalculate the
  basis inverse row norms.
*/
  if (phase == dyPRIMAL1 || phase == dyPRIMAL2)
  { dy_pseinit() ; }
  else
  if (phase == dyDUAL && retval != dyrLOSTDFEAS)
  { dy_dseinit() ; }

  return (retval) ; }



static int factor_loadcol (void *p_consys, int i, int *rndx, double *coeff)

/*
  This routine is used by luf_decomp to load columns of the basis into its
  internal data structure. The requirements are that it load rndx[] and
  coeff[] with the row indices and coefficients, respectively, of column i
  of the basis, returning the number of coefficients.

  Here, we need to look up the index j of the variable that actually occupies
  basis position i, retrieve the column, and load it into rndx and coeff.

  Parameters:
    p_consys:	a consys_struct
    i:		basis index
    rndx:	(o) row indices i of coefficients a<ij>
    coeff:	(o) coefficients a<ij>

  Returns: number of coefficients, or -1 in the event of an error.
*/

{ int j,vecndx,pkndx ;
  double aij ;
  pkvec_struct *aj ;
  consys_struct *consys ;

  const char *rtnnme = "factor_loadcol" ;

# ifdef PARANOIA
  if (p_consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (-1) ; }
  if (rndx == NULL)
  { errmsg(2,rtnnme,"row index vector") ;
    return (-1) ; }
  if (coeff == NULL)
  { errmsg(2,rtnnme,"coefficient") ;
    return (-1) ; }
# endif

  consys = (consys_struct *) p_consys ;
  aj = pkvec_new(consys->maxcollen) ;

/*
  Retrieve the column ...
*/
  j = dy_basis[i] ;
  if (consys_getcol_pk(consys,j,&aj) == FALSE)
  { errmsg(112,rtnnme,dy_sys->nme,"retrieve","column",
	   consys_nme(dy_sys,'v',j,FALSE,NULL),j) ;
    if (aj != NULL) pkvec_free(aj) ;
    return (-1) ; }
/*
  ... and load it into the return vectors rndx and coeff.
*/
  vecndx = 1 ;
  for (pkndx = 0 ; pkndx < aj->cnt ; pkndx++)
  { aij = aj->coeffs[pkndx].val ;
    if (aij != 0.0)
    { rndx[vecndx] = aj->coeffs[pkndx].ndx ;
      coeff[vecndx] = aij ;
      vecndx++ ; } }
/*
  Clean up and return.
*/
  pkvec_free(aj) ;
  return (vecndx-1) ; }




dyret_enum dy_factor (flags *calcflgs)

/*
  This routine orchestrates the LU factorisation of the basis. The glpk
  routines do the grunt work. This routine provides the intelligence.

  If inv_decomp aborts the attempt to factor due to numerical instability, we
  tighten the pivot selection parameters one notch and try again, giving up
  only when no further increase is possible.  The sequence of values for the
  pivot selection parameters are defined in a table at the top of this file.

  If inv_decomp aborts the attempt to factor because the basis is singular,
  we correct the basis with adjust_basis and take another run at factoring.
  In the event that the basis is successfully patched, we have serious work
  to do.  See the comments with adjust_therest for further information. If
  the user has for some reason disabled basis patching, we return
  dyrSINGULAR.

  inv_decomp (actually, luf_decomp) is self-expanding --- if more space is
  needed to hold the factorization, the expansion is handled internally.
  dylp uses ladEXPAND to force basis expansion after a pivot fails due to lack
  of space. In glpk, inv_update will set instructions in the basis structure
  and luf_decomp will handle the expansion, so ladEXPAND is redundant. No
  action need be taken in this routine. It's also not possible to tell if the
  basis has been expanded, so ladEXPAND is not set on output.


  Parameters:
    calcflgs:   (i) ladPRIMALS indicates the primal variables should be
		    recalculated after factoring the basis.
		    ladDUALS indicates the dual variables should be
		    recalculated after factoring the basis.
		    ladEXPAND indicates that the basis should be expanded prior
		    to refactoring.
		(o) flags are set to indicate if the corresponding variables
		    have been recalculated.

  Returns: dyrOK if the basis is factored without incident
	   dyrPATCHED if the basis was singular and has been repaired
	   dyrSINGULAR if the basis was singular and has not been repaired
	   dyrNUMERIC if factoring failed for the strictest pivoting regimen
	   dyrFATAL for other fatal errors

  NOTE: glpinv/glpluf will crash and burn if they encounter what they consider
	to be a fatal error, rather than returning a fatal error code. This
	needs to be addressed at some point. In particular, failure to expand
	the basis, failure to load the basis from the constraint system, and
	various parameter errors fall into this category.
*/

{ int retval,patchcnt ;
  bool try_again,patched ;
  dyret_enum retcode ;
  patch_struct *patches ;

  const char *rtnnme = "dy_factor" ;

#ifdef PARANOIA
  if (dy_sys == NULL)
  { errmsg(2,rtnnme,"dy_sys") ;
    return (dyrFATAL) ; }
  if (dy_basis == NULL)
  { errmsg(2,rtnnme,"basis") ;
    return (dyrFATAL) ; }
#endif

# ifdef DYLP_STATISTICS
  if (dy_stats != NULL)
  { int pivcnt ;
    pivcnt = dy_lp->tot.pivs-dy_stats->factor.prevpiv ;
    dy_stats->factor.avgpivs = dy_stats->factor.avgpivs*dy_stats->factor.cnt ;
    dy_stats->factor.avgpivs += pivcnt ;
    dy_stats->factor.cnt++ ;
    dy_stats->factor.avgpivs /= dy_stats->factor.cnt ;
    if (pivcnt > dy_stats->factor.maxpivs) dy_stats->factor.maxpivs = pivcnt ;
    dy_stats->factor.prevpiv = dy_lp->tot.pivs ; }
# endif

  retcode = dyrINV ;

/*
  Call luf_adjustsize to set the actual size of the basis. If the allocated
  capacity is too small, it will be expanded.
*/
  luf_adjustsize() ;
/*
  Open a loop for factorisation attempts. We'll persist in the face of
  numerical stability problems as long as there's room to tighten the pivot
  selection.

  At present, glpinv/glpluf will crash and burn if they encounter fatal
  problems. The basis load is implicit --- the routine factor_loadcol is
  called from luf_decomp to load up the coefficients.
*/
  try_again = TRUE ;
  patched = FALSE ;
  while (try_again)
  { retval = inv_decomp(luf_basis,dy_sys,factor_loadcol) ;
#   ifndef DYLP_NDEBUG
    if ((retval == 0 && dy_opts->print.basis >= 4) ||
	(retval > 0 && dy_opts->print.basis >= 2))
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    (%s)%d: factored with %s, basis stability %g.",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dy_prtpivparms(-1),luf_basis->min_vrratio) ; }
#   endif
/*
  Deal with the result. A return code of 0 means there were no difficulties;
  1 says the basis was singular and had to be patched before the
  factorisation could be completed. Either is success, and we're done.
*/
    switch (retval)
    { case 0:
      { try_again = FALSE ;
	retcode = dyrOK ;
	break ; }
/*
  Alas, the failures.

  If the problem is a singular basis (retval = 1), fix up the basis structures
  as indicated in the luf_basis structure and try again to factor the basis,
  unless the user has forbidden it.

  If the problem is numerical instability (retval = 2) try to make the pivot
  selection more stringent, and keep trying until we can try no more, at
  which point we'll return numeric instability to the caller.

  What's left is fatal confusion; pass the buck back to the caller.
*/
      case 1:
      { if (dy_opts->patch == FALSE)
	{ errmsg(308,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,dy_prtdyret(dyrSINGULAR)) ;
	  clrflg(*calcflgs,ladPRIMALS|ladDUALS) ;
	  return (dyrSINGULAR) ; }
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.basis >= 2)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n    (%s)%d: attempting to patch singular basis.",
		      dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
#	endif
	adjust_basis(&patchcnt,&patches) ;
	patched = TRUE ;
	break ; }
      case 2:
      { retcode = dyrNUMERIC ;
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.basis >= 2)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n    (%s)%d: factor failed at %s, numerical instability,",
		  dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		  dy_prtpivparms(-1)) ;
	  dyio_outfmt(dy_logchn,dy_gtxecho," max = %g, gro = %g.",
		      luf_basis->luf->big_v,luf_basis->luf->max_gro) ; }
# 	endif
	if (dy_setpivparms(+1,0) == FALSE)
	{ errmsg(307,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
		 dy_lp->tot.iters,dy_prtpivparms(-1)) ;
	  return (retcode) ; }
#	ifndef DYLP_NDEBUG
	if (dy_opts->print.basis >= 2)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n\ttrying again with %s.",
		      dy_prtpivparms(-1)) ; }
#	endif
	break ; }
      default:
      { errmsg(7,rtnnme,__LINE__,"inv_decomp return code",retval) ;
	return (dyrFATAL) ; } }
  }
/*
  If we reach here, we managed to factor the basis.  Reset the count of
  pivots since the last refactor.  If the basis was patched, we have some
  serious cleanup to do, so call adjust_therest to deal with the details.
  Otherwise, turn to the requests to calculate values for the primal and/or
  dual variables.
*/
  dy_lp->basis.etas = 0 ;
  if (patched == TRUE)
  { retcode = adjust_therest(patchcnt,patches) ;
    FREE(patches) ;
    if (retcode == dyrFATAL)
    { errmsg(306,rtnnme,dy_sys->nme,dy_prtlpphase(dy_lp->phase,TRUE),
	     dy_lp->tot.iters) ;
      return (dyrFATAL) ; }
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.basis >= 1)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\t[%s]: compensated for basis correction.",
		  dy_sys->nme) ; }
#   endif
    if (!(dy_lp->phase == dyINIT))
    { setflg(*calcflgs,ladPRIMALS|ladDUALS) ;
      if (retcode == dyrLOSTDFEAS) setflg(*calcflgs,ladDUALFEAS) ; }
    retcode = dyrPATCHED ; }
  else
  { if (flgon(*calcflgs,ladPRIMALS))
    { if (dy_calcprimals() == FALSE)
      { clrflg(*calcflgs,ladPRIMALS) ;
	return (dyrFATAL) ; } }
    if (flgon(*calcflgs,ladDUALS)) dy_calcduals() ; }

  return (retcode) ; }




dyret_enum dy_pivot (int xipos, double abarij, double maxabarj)

/*
  This routine handles a single pivot. It first checks that the pivot element
  satisfies a stability test, then calls inv_update to pivot the basis. We
  can still run into trouble, however, if the pivot results in a singular or
  near-singular basis.
  
  NOTE: There is an implicit argument here that's not immediately obvious.
	inv_update gets the entering column from a cached result set with the
	most recent call to inv_ftran(*,1) (dy_ftran(*,true), if you prefer).
	The underlying assumption is that this is readily available from when
	we ftran'd the entering column to find the leaving variable.

  Parameters:
    xipos:	the basis position of the entering variable
    abarij:	the pivot element (only the absolute value is used)
    maxabarj:	for a primal pivot, max{i} |abar<i,j>|,
		for a dual pivot, max{j} |abar<i,j>|

  Returns:
    dyrOK:	the pivot was accomplished without incident (inv_update)
    dyrMADPIV:	the pivot element abar<i,j> was rejected as numerically
		unstable (dy_chkpiv)
    dyrSINGULAR: the pivot attempt resulted in a structurally singular basis
		(i.e., some diagonal element is zero) (inv_update)
    dyrNUMERIC:	the pivot attempt resulted in a numerically singular (unstable)
		basis (i.e, some diagonal element is too small compared to
		other elements in the associated row and column) (inv_update)
    dyrBSPACE:	glpinv/glpluf ran out of space for the basis representation
		(inv_update)
    dyrFATAL:	internal confusion
*/

{ int retval ;
  double ratio ;
  dyret_enum retcode ;

  const char *rtnnme = "dy_pivot" ;

/*
  Check that the pivot element meets the current criterion for numerical
  stability. Arguably this should have been checked by the caller, but that's
  no excuse for not doing it now.
*/
  ratio = dy_chkpiv(abarij,maxabarj) ;
  if (ratio < 1.0)
  {
#   ifndef DYLP_NDEBUG
    if (dy_opts->print.basis >= 3)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n      %s(%d) pivot aborted; est. pivot stability %g.",
		  dy_prtlpphase(dy_lp->phase,TRUE),
		  dy_lp->tot.iters,rtnnme,ratio) ; }
#   endif
    return (dyrMADPIV) ; }
/*
  Make the call to inv_update, then recode the result.
*/
  retval = inv_update(luf_basis,xipos) ;
# ifndef DYLP_NDEBUG
  if ((retval == 0 && dy_opts->print.basis >= 5) ||
      (retval > 0 && dy_opts->print.basis >= 3))
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		"\n    %s(%d) estimated pivot stability %g; ",
		dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,ratio) ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"measured pivot stability %g.",
	        luf_basis->min_vrratio) ; }
# endif
  switch (retval)
  { case 0:
    { retcode = dyrOK ;
      break ; }
    case 1:
    { retcode = dyrSINGULAR ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.basis >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    %s(%d) singular basis (structural) after pivot.",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
#     endif
      break ; }
    case 2:
    { retcode = dyrNUMERIC ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.basis >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,
		    "\n    %s(%d) singular basis (numeric) after pivot.",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters) ; }
#     endif
      break ; }
    case 3:
    case 4:
    { retcode = dyrBSPACE ;
#     ifndef DYLP_NDEBUG
      if (dy_opts->print.basis >= 2)
      { dyio_outfmt(dy_logchn,dy_gtxecho,"\n    %s(%d) out of space (%s)",
		    dy_prtlpphase(dy_lp->phase,TRUE),dy_lp->tot.iters,
		    (retval == 3)?"eta matrix limit":"sparse vector area") ; }
#     endif
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      retcode = dyrFATAL ;
      break ; } }
  
  return (retcode) ; }


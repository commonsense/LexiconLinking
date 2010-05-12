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
  This file contains i/o routines related to the dylp subroutine library.
*/

#define DYLP_INTERNAL

#include "dylp.h"

static char sccsid[] UNUSED = "@(#)dylp_io.c	4.5	11/06/04" ;
static char svnid[] UNUSED = "$Id: dylp_io.c 94 2006-06-29 23:06:51Z lou $" ;



const char *dy_prtlpret (lpret_enum lpret)

/*
  Generates a print string corresponding to the dylp return codes.

  Parameters:
    lpret:	lp return code

  Returns: print string
*/

{ const char *rtnnme = "dy_prtlpret" ;

  switch (lpret)
  { case lpINV:
    { return ("invalid") ; }
    case lpOPTIMAL:
    { return ("optimal") ; }
    case lpUNBOUNDED:
    { return ("unbounded") ; }
    case lpSWING:
    { return ("pseudo-unbounded") ; }
    case lpINFEAS:
    { return ("infeasible") ; }
    case lpACCCHK:
    { return ("accuracy check") ; }
    case lpSTALLED:
    { return ("stalled") ; }
    case lpITERLIM:
    { return ("iteration limit exceeded") ; }
    case lpNOSPACE:
    { return ("insufficient memory") ; }
    case lpLOSTFEAS:
    { return ("lost feasibility") ; }
    case lpPUNT:
    { return ("punt!") ; }
    case lpFORCEDUAL:
    { return ("force primal -> dual") ; }
    case lpFORCEPRIMAL:
    { return ("force dual -> primal") ; }
    case lpFORCEFULL:
    { return ("force full activation") ; }
    case lpFATAL:
    { return ("unspecified fatal error") ; }
    default:
    { errmsg(5,rtnnme,"lpret",(int) lpret) ;
      return ("nonsense") ; } } }



char *dy_prtvstat (flags status)

/*
  This routine returns a print string corresponding to the status code passed
  as a parameter.

  Parameter:
    status:	variable status code

  Returns: print string
*/

{ flags mystatus ;
  static char buffer[100] ;
  const char *rtnnme = "dy_prtvstat" ;

/*
  If we've been passed a completely empty status, return invalid. But we could
  be asked to print just a qualifier, so be prepared for mystatus to be
  vstatINV once the qualifiers are cleared.
*/
  buffer[0] = '\0' ;
  if (status != vstatINV)
  { mystatus = status ;
    clrflg(mystatus,vstatQUALS) ; }
  else
  { strcpy(buffer,"INV") ;
    return (buffer) ; }

  if (mystatus != vstatINV)
    switch (mystatus)
    { case vstatBFX:
      { strcpy(buffer,"BFX") ;
	break ; }
      case vstatBUUB:
      { strcpy(buffer,"BUUB") ;
	break ; }
      case vstatBUB:
      { strcpy(buffer,"BUB") ;
	break ; }
      case vstatB:
      { strcpy(buffer,"B") ;
	break ; }
      case vstatBLB:
      { strcpy(buffer,"BLB") ;
	break ; }
      case vstatBLLB:
      { strcpy(buffer,"BLLB") ;
	break ; }
      case vstatBFR:
      { strcpy(buffer,"BFR") ;
	break ; }
      case vstatNBFX:
      { strcpy(buffer,"NBFX") ;
	break ; }
      case vstatNBUB:
      { strcpy(buffer,"NBUB") ;
	break ; }
      case vstatNBLB:
      { strcpy(buffer,"NBLB") ;
	break ; }
      case vstatNBFR:
      { strcpy(buffer,"NBFR") ;
	break ; }
      case vstatSB:
      { strcpy(buffer,"SB") ;
	break ; }
      case vstatINV:
      { strcpy(buffer,"INV") ;
	break ; }
      default:
      { errmsg(6,rtnnme,"status",(int) status) ;
        strcpy(buffer,"NONSENSE") ;
	return (buffer) ; } }
/*
  Add any qualifiers.
*/
  if (status != mystatus)
  { strcat(buffer,"(") ;
    if (flgon(status,vstatNOPIVOT)) strcat(buffer,"r") ;
    if (flgon(status,vstatNOPER)) strcat(buffer,"p") ;
    strcat(buffer,")") ; }

  return (buffer) ; }



const char *dy_prtlpphase (dyphase_enum phase, bool abbrv)

/*
  This routine returns a print representation of the lp phase passed as a 
  parameter.

  Parameter:
    phase:	an lp phase code.
    abbrv:	TRUE to get a two-letter abbreviation, FALSE for the long
		version

  Returns: print string
*/

{ const char *rtnnme = "dy_prtlpphase" ;

  switch (phase)
  { case dyINIT:
    { return ((abbrv == TRUE)?"IN":"initialisation") ; }
    case dyPURGEVAR:
    { return ((abbrv == TRUE)?"VD":"variable deactivation") ; }
    case dyGENVAR:
    { return ((abbrv == TRUE)?"VG":"variable generation") ; }
    case dyADDVAR:
    { return ((abbrv == TRUE)?"VA":"variable activation") ; }
    case dyPRIMAL1:
    { return ((abbrv == TRUE)?"P1":"primal phase I") ; }
    case dyPRIMAL2:
    { return ((abbrv == TRUE)?"P2":"primal phase II") ; }
    case dyPURGECON:
    { return ((abbrv == TRUE)?"CD":"constraint deactivation") ; }
    case dyGENCON:
    { return ((abbrv == TRUE)?"CG":"constraint generation") ; }
    case dyADDCON:
    { return ((abbrv == TRUE)?"CA":"constraint activation") ; }
    case dyDUAL:
    { return ((abbrv == TRUE)?"D2":"dual") ; }
    case dyFORCEDUAL:
    { return ((abbrv == TRUE)?"FD":"force dual") ; }
    case dyFORCEPRIMAL:
    { return ((abbrv == TRUE)?"FP":"force primal") ; }
    case dyFORCEFULL:
    { return ((abbrv == TRUE)?"FF":"force full") ; }
    case dyDONE:
    { return ((abbrv == TRUE)?"DN":"done") ; }
    case dyINV:
    { return ((abbrv == TRUE)?"NV":"invalid") ; }
    default:
    { errmsg(6,rtnnme,"lp phase",(int) phase) ;
      return ((abbrv == TRUE)?"??":"nonsense") ; } } }



const char *dy_prtdyret (dyret_enum retcode)

/*
  This routine returns a print representation of the dyret_enum code passed
  as a parameter.

  Parameter:
    retcode:	a dyret_enum return code

  Returns: print string
*/

{ const char *rtnnme = "dy_prtdyret" ;

  switch (retcode)
  { case dyrOK:
    { return ("ok") ; }
    case dyrOPTIMAL:
    { return ("optimal") ; }
    case dyrUNBOUND:
    { return ("unbounded") ; }
    case dyrSWING:
    { return ("pseudo-unbounded") ; }
    case dyrINFEAS:
    { return ("infeasible") ; }
    case dyrREQCHK:
    { return ("request accuracy check") ; }
    case dyrACCCHK:
    { return ("accuracy check failure") ; }
    case dyrLOSTPFEAS:
    { return ("loss of primal feasibility") ; }
    case dyrLOSTDFEAS:
    { return ("loss of dual feasibility") ; }
    case dyrDEGEN:
    { return ("degenerate pivot") ; }
    case dyrRESELECT:
    { if (dy_lp->phase == dyDUAL)
	return ("reselect leaving variable") ;
      else
	return ("reselect entering variable") ; }
    case dyrMADPIV:
    { return ("numerically unstable pivot") ; }
    case dyrPUNT:
    { return ("punt!") ; }
    case dyrPATCHED:
    { return ("basis patched") ; }
    case dyrSINGULAR:
    { return ("basis singular") ; }
    case dyrNUMERIC:
    { return ("ill-conditioned basis") ; }
    case dyrBSPACE:
    { return ("no space for basis") ; }
    case dyrSTALLED:
    { return ("stalled") ; }
    case dyrITERLIM:
    { return ("iteration limit") ; }
    case dyrFATAL:
    { return ("fatal error") ; }
    case dyINV:
    { return ("invalid") ; }
    default:
    { errmsg(6,rtnnme,"dyret_enum code",(int) retcode) ;
      return ("nonsense") ; } } }



void dy_logpivot (dyret_enum result, int xjndx, int indir, double cbarj,
		  int xindx, int outdir, double abarij, double delta)

/*
  This routine prints a standard log line for a pivot.

  Parameters:
    result:	the return code resulting from attempting the pivot
    xjndx:	index of the entering variable x<j>
    indir:	direction of movement of x<j> (1 to increase, -1 to decrease)
    cbarj:	reduced cost cbar<j> for the entering variable
    xindx:	index of the leaving variable x<i>
    outdir:	direction of motion of x<i>
    abarij:	pivot element abar<i,j>
    delta:	amount of change in x<j>

  Returns: undefined
*/

{ bool validin,validout ;
  const char *resstr ;

/*
  logpivot is called from within the dual simplex routine, so we have to
  convert dual unboundedness to primal infeasibility. Pseudo-unboundedness
  (swing) always refers to primal variables.
*/
  validin = TRUE ;
  validout = TRUE ;

  switch (result)
  { case dyrOK:
    { resstr = "(ok)" ;
      break ; }
    case dyrDEGEN:
    { resstr = "(degen)" ;
      break ; }
    case dyrUNBOUND:
    { if (dy_lp->phase == dyDUAL)
      { resstr = "(infea)" ;
	validin = FALSE ; }
      else
      { resstr = "(unbnd)" ;
	validout = FALSE ; }
      break ; }
    case dyrSWING:
    { resstr = "(swing)" ;
      break ; }
    case dyrLOSTPFEAS:
    { resstr = "(!pfea)" ;
      break ; }
    case dyrLOSTDFEAS:
    { resstr = "(!dfea)" ;
      break ; }
    case dyrOPTIMAL:
    { if (dy_lp->phase == dyPRIMAL1)
	resstr = "(infea)" ;
      else
	resstr = "(opt)" ;
      break ; }
    case dyrPUNT:
    { resstr = "(punt!)" ;
      if (xjndx <= 0)
	validin = FALSE ;
      break ; }
    case dyrREQCHK:
    { if (dy_lp->pivok == FALSE)
	resstr = "(chkab)" ;
      else
	resstr = "(chkrq)" ;
      break ; }
    case dyrMADPIV:
    { resstr = "(mad)" ;
      if (xjndx <= 0)
	validin = FALSE ;
      break ; }
    case dyrSINGULAR:
    { resstr = "(sing)" ;
      break ; }
    case dyrBSPACE:
    { resstr = "(nosp)" ;
      break ; }
    case dyrFATAL:
    { resstr = "(fatal)" ;
      break ; }
    case dyrRESELECT:
    { resstr = "(resel)" ;
      if (dy_lp->phase == dyDUAL)
	validout = TRUE ;
      break ; }
    default:
    { resstr = "(huh?)" ;
      result = dyrINV ;
      break ; } }
  dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\n%s%6d %-7s ",dy_prtlpphase(dy_lp->phase,TRUE),
	      dy_lp->tot.iters+1,resstr) ;

  if (result == dyrINV) return ;

  if (validin == TRUE && xjndx > 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho,"In: %s (%d) %s cbarj = %g ;",
	        consys_nme(dy_sys,'v',xjndx,FALSE,NULL),xjndx,
	        (indir == 1)?"inc":"dec",cbarj) ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho,"In: <not selected>") ; }
  
  if (result == dyrFATAL) return ;

  if (result == dyrLOSTPFEAS)
  { dyio_outfmt(dy_logchn,dy_gtxecho,
		" Infeas: %s (%d) = %g, lb = %g, ub = %g",
	        consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
	        dy_xbasic[dy_var2basis[xindx]],
	        dy_sys->vlb[xindx],dy_sys->vub[xindx]) ;
    return ; }

  if (validout == TRUE && xindx > 0)
  { dyio_outfmt(dy_logchn,dy_gtxecho," Out: %s (%d) %s",
	        consys_nme(dy_sys,'v',xindx,FALSE,NULL),xindx,
	        (outdir == 1)?"inc":"dec") ; }
  else
  { dyio_outfmt(dy_logchn,dy_gtxecho," Out: <not selected>") ; }
  
  if (validin == TRUE && validout == TRUE)
  { dyio_outfmt(dy_logchn,dy_gtxecho,", abarij = %g, delta = %g",
	        abarij,(indir == 1)?delta:-delta) ; }

  if (dy_lp->phase == dyDUAL)
  { dyio_outfmt(dy_logchn,dy_gtxecho,", yb = %g.",dy_calcdualobj()) ; }
  else
  if (dy_lp->phase == dyPRIMAL1)
  { dyio_outfmt(dy_logchn,dy_gtxecho,", infeas = %g.",dy_calcpinfeas()) ; }
  else
  if (dy_lp->phase == dyPRIMAL2)
  { dyio_outfmt(dy_logchn,dy_gtxecho,", cx = %g.",dy_calcobj()) ; }
  else
  { dyio_outchr(dy_logchn,dy_gtxecho,'.') ; }

  return ; }



bool dy_dumpcompact (ioid chn, bool echo, lpprob_struct *soln, bool nbzeros)

/*
  This routine prints the compact version of the solution, as returned by
  dylp. The layout of the solution is basic information, followed by nonbasic
  information.

  For each active constraint, we print the basis position; the constraint
  index and name; the basic variable name, index, status and value; and the
  dual variable and value.

  For each nonbasic variable, we print the variable name, index, status, and
  value, and reduced cost. The caller can optionally suppress `common' zeros
  --- e.g., variables with status NBLB, NBFX, or NBUB and value of 0.

  Note that unscaling is not required here. soln->x and soln->y were unscaled
  when the solution was generated, and the client's constraint system
  (soln->consys) is not touched when dylp scales.

  Parameters:
    chn:	file channnel for output
    echo:	TRUE to echo to stdout, FALSE otherwise
    soln:	an lpprob structure, containing a solution as returned by
		dylp.
    nbzeros:	TRUE to print all nonbasic variables with a value of 0,
		FALSE to print only nonbasic variables with nonzero value
		or with status NBFR or SB.
  
  Returns: TRUE if the solution could be printed without error, FALSE
	   otherwise.
*/

{ int vndx,cndx,bpos ;
  double val ;
  bool nononbasic ;
  consys_struct *sys ;
  basis_struct *basis ;
  const char *rtnnme = "dy_dumpcompact" ;

# ifdef PARANOIA
  if (soln == NULL)
  { errmsg(2,rtnnme,"solution") ;
    return (FALSE) ; }
  sys = soln->consys ;
  if (sys == NULL)
  { errmsg(2,rtnnme,"constraint system") ;
    return (FALSE) ; }
# else
  sys = soln->consys ;
# endif

/*
  Begin by printing identifying information about the system and the solution.
  If the phase is anything but dyDONE, we're done too.
*/
  dyio_outfmt(chn,echo,
	      "\n\nSystem: %s\t\t\tfinal status: %s after %d iterations.",
	      sys->nme,dy_prtlpphase(soln->phase,FALSE),soln->iters) ;
  if (soln->phase != dyDONE)
  { dyio_outchr(chn,echo,'\n') ;
    return (TRUE) ; }
/*
  Consider the lp return code. If it's optimal, infeasible, or unbounded, we'll
  continue on to print the solution, otherwise we're done.
*/
  dyio_outfmt(chn,echo,"\n    lp status: %s",dy_prtlpret(soln->lpret)) ;
  switch (soln->lpret)
  { case lpOPTIMAL:
    { dyio_outfmt(chn,echo,"\t\tobjective: %.9g",soln->obj) ;
      break ; }
    case lpINFEAS:
    { dyio_outfmt(chn,echo,"\t\tinfeasibility: %.9g",soln->obj) ;
      break ; }
    case lpUNBOUNDED:
    { if (soln->obj != 0)
      { if (soln->obj < 0)
	{ vndx = abs((int) soln->obj) ;
	  bpos = -1 ; }
	else
	{ vndx = (int) soln->obj ;
	  bpos = 1 ; }
	dyio_outfmt(chn,echo,"\t\tunbounded variable %s (%d) (%s)",
		    consys_nme(sys,'v',vndx,FALSE,NULL),vndx,
		    (bpos < 0)?"decreasing":"increasing") ; }
      break ; }
    default:
    { dyio_outchr(chn,echo,'\n') ;
      return (TRUE) ; } }
/*
  There's a solution to be dumped. Do the basis, duals, and basic variables
  first.
*/
  dyio_outfmt(chn,echo,"\n\nPosn\tConstraint\tDual\t\tPrimal\n") ;
  basis = soln->basis ;
  for (bpos = 1 ; bpos <= basis->len ; bpos++)
  { cndx = basis->el[bpos].cndx ;
    vndx = basis->el[bpos].vndx ;
    if (vndx < 0) vndx = sys->varcnt-vndx ;
    dyio_outfmt(chn,echo,"\n%5d\t(%4d) %-8s\t%12.4g\t(%4d) %-8s %12.7g",
	        bpos,cndx,consys_nme(sys,'c',cndx,FALSE,NULL),soln->y[bpos],
	        vndx,consys_nme(sys,'v',vndx,FALSE,NULL),soln->x[bpos]) ; }
/*
  Now the nonbasic variables. Nonzero values only.
*/
  nononbasic = TRUE ;
  for (vndx = 1 ; vndx <= sys->varcnt ; vndx++)
    if ((int) soln->status[vndx] > 0)
    { if (nononbasic == TRUE)
      { dyio_outfmt(chn,echo,"\n\nNonbasic Primal\n") ;
	nononbasic = FALSE ; }
      switch (soln->status[vndx])
      { case vstatNBLB:
	case vstatNBFX:
	{ val = sys->vlb[vndx] ;
	  if (nbzeros == FALSE && val == 0) continue ;
	  break ; }
	case vstatNBUB:
	{ val = sys->vub[vndx] ;
	  if (nbzeros == FALSE && val == 0) continue ;
	  break ; }
	case vstatNBFR:
	case vstatSB:
	{ val = 0 ;
	  break ; }
	default:
	{ val = quiet_nan(0) ;
	  errmsg(1,rtnnme,__LINE__) ;
	  break ; } }
      dyio_outfmt(chn,echo,"\n(%4d) %-8s %3s %12.7g",vndx,
		  consys_nme(sys,'v',vndx,FALSE,NULL),
		  dy_prtvstat(soln->status[vndx]),val) ; }
  
  if (nononbasic == TRUE)
    dyio_outfmt(chn,echo,"\n\nNo nonbasic architectural variables.\n") ;
  else
    dyio_outchr(chn,echo,'\n') ;

  return (TRUE) ; }


/*
  This file is a portion of the OsiDylp LP distribution.

        Copyright (C) 2004 Lou Hafer

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
  with this program; if not, write to the Free Software Foundation, Inc., 59
  Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
  This file contains i/o routines for the constraint system data structures,
  and related routines for generating name strings for various objects.
*/



#include "dylib_errs.h"
#include "dylib_io.h"
#include "dylib_std.h"
#include "consys.h"

static char sccsid[] UNUSED = "@(#)consys_io.c	4.6	11/11/04" ;
static char svnid[] UNUSED = "$Id: consys_io.c 94 2006-06-29 23:06:51Z lou $" ;



const char *consys_prtcontyp (contyp_enum contyp)

/*
  Utility to print a readable string for a constraint type.

  Parameter:
    contyp:	constraint type

  Returns: appropriate string for the type, or an error string.
*/

{ const char *rtnnme = "consys_prtcontyp" ;

  switch (contyp)
  { case contypLE:
    { return ("<=") ; }
    case contypGE:
    { return (">=") ; }
    case contypEQ:
    { return ("=") ; }
    case contypNB:
    { return ("><") ; }
    case contypRNG:
    { return ("<=>") ; }
    case contypINV:
    { return ("invalid") ; }
    default:
    { errmsg(5,rtnnme,"contyp",(int) contyp) ;
      return ("unrecognised") ; } } }


const char *consys_prtvartyp (vartyp_enum vartyp)

/*
  Utility to print a readable string for a variable type.

  Parameter:
    vartyp:	variable type

  Returns: appropriate string for the type, or an error string.
*/

{ const char *rtnnme = "consys_prtvartyp" ;

  switch (vartyp)
  { case vartypCON:
    { return ("continuous") ; }
    case vartypINT:
    { return ("general integer") ; }
    case vartypBIN:
    { return ("binary") ; }
    case vartypINV:
    { return ("invalid") ; }
    default:
    { errmsg(5,rtnnme,"vartyp",(int) vartyp) ;
      return ("unrecognised") ; } } }



char *consys_assocnme (consys_struct *consys, flags which)

/*
  Utility routine to produce a name for an associated vector. If consys
  is non-NULL, the name is fully qualified (consys.which) using the short form
  of which. If consys is NULL, a longer form of which is returned.

  Parameters:
    consys:	constraint system
    which:	associated vector type (from codes in consys.h)
  
  Returns: a name string for the vector, or a string indicating error.
*/

{ static char nmbuf[128] ;
  int nmlen ;

  if (consys != NULL)
  { nmlen = sizeof(nmbuf)/2 ;
    (void) dyio_outfxd(nmbuf,-nmlen,'l',"%s",consys->nme) ;
    strcat(nmbuf,".") ; }
  else
  { nmbuf[0] = '\0' ; }
  nmlen = strlen(nmbuf) ;

  switch (which)
  { case CONSYS_MTX:
    { strcat(nmbuf,(consys == NULL)?"constraint matrix":"mtx") ;
      break ; }
    case CONSYS_ROWHDR:
    { strcat(nmbuf,(consys == NULL)?"row header array":"rowhdr") ;
      break ; }
    case CONSYS_COLHDR:
    { strcat(nmbuf,(consys == NULL)?"column header array":"colhdr") ;
      break ; }
    case CONSYS_OBJ:
    { strcat(nmbuf,(consys == NULL)?"objective function":"obj") ;
      break ; }
    case CONSYS_VUB:
    { strcat(nmbuf,(consys == NULL)?"variable upper bounds":"vub") ;
      break ; }
    case CONSYS_VLB:
    { strcat(nmbuf,(consys == NULL)?"variable lower bounds":"vlb") ;
      break ; }
    case CONSYS_RHS:
    { strcat(nmbuf,(consys == NULL)?"right-hand-side":"rhs") ;
      break ; }
    case CONSYS_RHSLOW:
    { strcat(nmbuf,(consys == NULL)?"range right-hand-side":"rhslow") ;
      break ; }
    case CONSYS_CUB:
    { strcat(nmbuf,(consys == NULL)?"constraint upper bounds":"cub") ;
      break ; }
    case CONSYS_CLB:
    { strcat(nmbuf,(consys == NULL)?"constraint lower bounds":"clb") ;
      break ; }
    case CONSYS_VTYP:
    { strcat(nmbuf,(consys == NULL)?"variable type":"vtyp") ;
      break ; }
    case CONSYS_CTYP:
    { strcat(nmbuf,(consys == NULL)?"constraint type":"ctyp") ;
      break ; }
    case CONSYS_RSCALE:
    { strcat(nmbuf,(consys == NULL)?"row scaling":"rsc") ;
      break ; }
    case CONSYS_CSCALE:
    { strcat(nmbuf,(consys == NULL)?"column scaling":"csc") ;
      break ; }
    case CONSYS_COL:
    { strcat(nmbuf,(consys == NULL)?"generic column":"col") ;
      break ; }
    case CONSYS_ROW:
    { strcat(nmbuf,(consys == NULL)?"generic row":"row") ;
      break ; }
    default:
    { dyio_outfxd(&nmbuf[nmlen],-26,'l',"<<type error: %#08x>>",(int) which) ;
      break ; } }

  return (nmbuf) ; }



char *consys_lognme (consys_struct *consys, int rowndx, char *clientbuf)

/*
  This is a utility routine that will construct the proper name for a logical
  variable, given the constraint index. Names constructed for logicals are
  guaranteed to be no more than 32 characters, including the final null.

  This routine exists only so that consys_nme and consys_utils:add_logical
  are guaranteed to produce the same names. They are the only two routines
  that should call consys_lognme, and add_logical is the only reason this
  routine isn't static.

  Parameters:
    consys:	the constraint system
    rowndx:	the constraint index
    clientbuf:	a character buffer >= 32 characters, or NULL
  
  Returns: the user's buffer, if supplied, or an internal buffer, containing
	   the name of the logical variable.
*/

{ int len ;
  rowhdr_struct *rowhdr ;
  char *nmebuf ;
  static char ownbuf[32] ;

# ifdef PARANOIA

  const char *rtnnme = "consys_lognme" ;
/*
  The usual, and checks that some necessary associated arrays are present.
  This is an internal routine, so index bounds checks are dropped unless
  we're being paranoid.
*/
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->mtx.rows == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_ROWHDR)) ;
    return (FALSE) ; }
  if (consys->ctyp == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CTYP)) ;
    return (FALSE) ; }
  if (rowndx <= 0 || rowndx > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"row",rowndx,1,consys->concnt) ;
    return (FALSE) ; }
  if (consys->mtx.rows[rowndx] == NULL)
  { errmsg(103,rtnnme,consys->nme,"row",rowndx) ;
    return (FALSE) ; }
# endif

/*
  Use the client's buffer, if supplied, otherwise use the local static buffer.
*/
  if (clientbuf == NULL)
  { nmebuf = ownbuf ; }
  else
  { nmebuf = clientbuf ; }
/*
  Construct a name, based on the type of constraint.
*/
  rowhdr = consys->mtx.rows[rowndx] ;
  len = strlen(rowhdr->nme) ;
  if (len > sizeof(nmebuf)-5) len = sizeof(nmebuf)-5 ;
  strncpy(nmebuf,rowhdr->nme,len) ;
  switch (consys->ctyp[rowndx])
  { case contypLE:
    { strcpy(&nmebuf[len],".slk") ;
      break ; }
    case contypEQ:
    { strcpy(&nmebuf[len],".art") ;
      break ; }
    case contypGE:
    { strcpy(&nmebuf[len],".sur") ;
      break ; }
    case contypRNG:
    { strcpy(&nmebuf[len],".rng") ;
      break ; }
    default:
    { strcpy(&nmebuf[len],".inv") ;
      break ; } }

  return (nmebuf) ; }


const char *consys_nme (consys_struct *consys,
			char cv, int ndx, bool pfx, char *clientbuf)

/*
  Utility routine to retrieve the name of a constraint or variable. If pfx is
  false, the base name is returned. If pfx is true, 'consys->nme.' is added
  as a prefix.

  If the constraint or variable name is stored in a row or column header, the
  stored pointer is returned.

  The name has to be constructed in two cases:
    * when the prefixed form is requested, and
    * when  the name of a logical is requested and logicals aren't enabled.
  In these cases, the name is built in a buffer and a pointer to the buffer
  is returned. If the client doesn't supply a buffer, an internal static
  buffer is used (which will be overwritten the next time it's needed).

  If the user supplies a buffer, that buffer is always used (whether the name
  is constructed or not). The constant CONSYS_MAXBUFLEN is the maximum buffer
  length.

  Parameters:
    consys:	constraint system
    cv:		'c' for constraint, 'v' for variable
    ndx:	the constraint/variable (row/column) index; to request the
		name of the logical for constraint k when logicals aren't
		enabled, use (consys->varcnt+k)
    pfx:	TRUE if the fully qualified name should be generated,
		FALSE for the constraint/variable name only.
    clientbuf:	if non-NULL, the name is constructed and returned in this
		buffer.
  
  Returns: the name of the constraint/variable, or some appropriate string
	   indicating error.
*/

{ static char ourbuf[CONSYS_MAXBUFLEN],ourbuftoo[CONSYS_MAXBUFLEN] ;
  char *nmbuf ;
  const char *rtnbuf ;
  int nmlen,partlen ;

#ifdef PARANOIA

  const char *rtnnme = "consys_nme",
	     *errname = "<<error>>" ;

  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (errname) ; }
  switch (cv)
  { case 'c':
    { if (consys->mtx.rows == NULL)
      { errmsg(101,rtnnme,consys->nme,"row header") ;
	return (errname) ; }
      if (ndx <= 0 || ndx > consys->concnt)
      { errmsg(102,rtnnme,consys->nme,"constraint",ndx,1,consys->concnt) ;
	return (errname) ; }
      if (consys->mtx.rows[ndx] == NULL)
      { errmsg(103,rtnnme,consys->nme,"row",ndx) ;
	return (errname) ; }
      break ; }
    case 'v':
    { if (consys->mtx.cols == NULL)
      { errmsg(101,rtnnme,consys->nme,"column header") ;
	return (errname) ; }
      if (flgon(consys->opts,CONSYS_LVARS))
	nmlen = consys->varcnt ;
      else
	nmlen = consys->varcnt+consys->concnt ;
      if (ndx <= 0 || ndx > nmlen)
      { errmsg(102,rtnnme,consys->nme,"variable",ndx,1,nmlen) ;
	return (errname) ; }
      if (ndx < consys->varcnt && consys->mtx.cols[ndx] == NULL)
      { errmsg(103,rtnnme,consys->nme,"column",ndx) ;
	return (errname) ; }
      break ; }
    default:
    { errmsg(3,rtnnme,"cv",cv) ;
      return (errname) ; } }
#endif

/*
  Can we just return the pointer to the name from the row/column header?
  (Perhaps after copying it into the client's buffer.)
*/
  if (pfx == FALSE && (cv == 'c' || (cv == 'v' && ndx <= consys->varcnt)))
  { if (cv == 'c')
      rtnbuf = consys->mtx.rows[ndx]->nme ;
    else
      rtnbuf = consys->mtx.cols[ndx]->nme ;
    if (clientbuf != NULL)
    { if (strlen(rtnbuf) < CONSYS_MAXBUFLEN)
	strcpy(clientbuf,rtnbuf) ;
      else
      { strncpy(clientbuf,rtnbuf,CONSYS_MAXBUFLEN-1) ;
	clientbuf[CONSYS_MAXBUFLEN-1] = '\0' ; }
      rtnbuf = clientbuf ; } }
/*
  We have to build the name. Not quite so bad as it seems at first glance.
  Figure out what buffer to use, then dump in the prefix, the dump in the
  name.

  The call to consys_lognme down below is made with a private buffer so that
  we don't inadvertently screw up by having a call to consys_nme interfere with
  a previous call to consys_lognme. Sure, it shouldn't happen, but why take
  the chance. The nme field in the row/column header shouldn't be null either,
  but we're going for robustness here -- this routine gets called a lot to
  generate error messages.
*/
  else
  { if (clientbuf == NULL)
      nmbuf = ourbuf ;
    else
      nmbuf = clientbuf ;
    
    if (pfx == TRUE)
    { nmlen = strlen(consys->nme) ;
      if (nmlen > CONSYS_MAXBUFLEN/2-1) nmlen = CONSYS_MAXBUFLEN/2-1 ;
      strncpy(nmbuf,consys->nme,nmlen) ;
      nmbuf[nmlen++] = '.' ; }
    else
    { nmlen = 0 ; }

    switch (cv)
    { case 'c':
      { if (consys->mtx.rows[ndx]->nme == NULL)
	{ strcpy(&nmbuf[nmlen],"<<null>>") ; }
	else
	{ partlen = strlen(consys->mtx.rows[ndx]->nme) ;
	  if (partlen > CONSYS_MAXBUFLEN-nmlen-1)
	    partlen = CONSYS_MAXBUFLEN-nmlen-1 ;
	  strncpy(&nmbuf[nmlen],consys->mtx.rows[ndx]->nme,partlen) ;
	  nmlen += partlen ;
	  nmbuf[nmlen] = '\0' ; }
	break ; }
      case 'v':
      { if (ndx <= consys->varcnt)
	{ if (consys->mtx.cols[ndx]->nme == NULL)
	  { strcpy(&nmbuf[nmlen],"<<null>>") ; }
	  else
	  { partlen = strlen(consys->mtx.cols[ndx]->nme) ;
	    if (partlen > CONSYS_MAXBUFLEN-nmlen-1)
	      partlen = CONSYS_MAXBUFLEN-nmlen-1 ;
	    strncpy(&nmbuf[nmlen],consys->mtx.cols[ndx]->nme,partlen) ;
	    nmlen += partlen ;
	    nmbuf[nmlen] = '\0' ; } }
	else
	{ (void) consys_lognme(consys,ndx-consys->varcnt,ourbuftoo) ;
	  partlen = strlen(ourbuftoo) ;
	  if (partlen > CONSYS_MAXBUFLEN-nmlen-1)
	    partlen = CONSYS_MAXBUFLEN-nmlen-1 ;
	  strncpy(&nmbuf[nmlen],ourbuftoo,partlen) ;
	  nmlen += partlen ;
	  nmbuf[nmlen] = '\0' ; }

	break ; } }
    rtnbuf = nmbuf ; }
  
  return (rtnbuf) ; }



/*
  A pair of utility routines to print constraint bound names and values.
*/

char *consys_conbndnme (char bndlett, int cndx, conbnd_struct *bnd)

/*
  Prints a constraint lower bound name as LB(vndx) (for a bound which is
  finite or has more than one infinite contribution) or LB(vndx\infndx) for a
  bound which has exactly one infinite contribution (infndx is the index of
  the variable contributing the infinity). Analogous for upper bounds.

  The routine is robust in the face of incorrect values for bnd.inf, in the
  sense that it will still produce a printable string. It doesn't actually
  check for errors --- that would require a constraint system as a parameter
  in order to check the limits on bnd.inf.

  Parameters:
    bndlett:    one of 'L' or 'U'
    cndx:	constraint index
    bnd:	constraint bound
  
  Returns: a name string for the bound
*/

{ static char buf[32] ;
  char *bufptr ;

  bufptr = &buf[0] ;
  bufptr += dyio_outfxd(bufptr,-((int) (sizeof(buf)/2-1)),
			'l',"%cB(%d",bndlett,cndx) ;
  if (bnd->inf < 0)
    dyio_outfxd(bufptr,-((int) (sizeof(buf)/2-1)),'l',"\\%d)",-bnd->inf) ;
  else
    dyio_outfxd(bufptr,-1,'l',")") ;
  
  return (&buf[0]) ; }


char *consys_conbndval (conbnd_struct *bnd)

/*
  Prints the constraint bound value, as nn*inf+val, or inf+val, or val,
  depending as the number of infinite contributions is >1, 1, or 0,
  respectively.

  As with consys_conbndnme, the routine is robust in the face of a bogus
  value for bnd.inf.

  Parameter:
    bnd:	constraint bound
  
  Returns: a value string for the bound
*/

{ static char buf[32] ;
  char *bufptr ;

  bufptr = &buf[0] ;
  if (bnd->inf > 0)
    bufptr += dyio_outfxd(bufptr,-((int) (sizeof(buf)/2-1)),
			  'l',"%d*inf+",bnd->inf) ;
  else
  if (bnd->inf < 0)
    bufptr += dyio_outfxd(bufptr,-((int) (sizeof(buf)/2-1)),'l',"inf+") ;
  
  dyio_outfxd(bufptr,-((int) (sizeof(buf)/2-1)),'l',"%g",bnd->bnd) ;

  return (&buf[0]) ; }



void consys_prtcon (ioid chn, bool echo,
		    consys_struct *consys, int i, const char *pfx)

/*
  This routine prints a constraint. Since this could be a lengthy string in
  larger problems, the routine prints directly to the output rather than to a
  string. The print is of the form:

    rhslow <= name (ndx) <= rhs
      coef*var(ndx) +/- coeff*var(ndx) +/- coeff*var(ndx) ...
  
  where coefficients take as many lines as necessary.

  Parameters:
    chn:	i/o channel
    echo:	TRUE to echo to tty, FALSE otherwise.
    consys:	reference constraint system for coefficients in newcon
    i:		index of the constraint to be printed
    pfx:	prefix to be printed at the start of each line
		(typically blank space for indentation)

  Returns: undefined
*/

{ int linecnt,charcnt,ndx ;
  contyp_enum ctypi ;
  pkvec_struct *coni ;
  pkcoeff_struct *ai ;
  char buf[64] ;
  const char *rtnnme = "consys_prtcon",
	     *errstring = "<< !consys_prtcon print error! >>",
	     *dfltpfx = "" ;


# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    dyio_outfmt(chn,echo,errstring) ;
    return ; }
  if (i < 0 || i > consys->concnt)
  { errmsg(102,rtnnme,consys->nme,"constraint",i,1,consys->concnt) ;
    dyio_outfmt(chn,echo,errstring) ;
    return ; }
# endif

  if (pfx == NULL) pfx = dfltpfx ;
/*
  The basics: name, index, constraint type, rhs (and rhslow, if needed).
*/
  ctypi = consys->ctyp[i] ;
  dyio_outfmt(chn,echo,"\n%s",pfx) ;
  if (ctypi == contypRNG)
  { dyio_outfmt(chn,echo,"%g <= ",consys->rhslow[i]) ; }
  dyio_outfmt(chn,echo,"%s (%d) %s %g",consys_nme(consys,'c',i,FALSE,NULL),i,
	 consys_prtcontyp(ctypi),consys->rhs[i]) ;
/*
  Now the coefficients, limiting each line to 80 characters. The initial
  coefficient on a line is indented and printed with an optional `-' followed
  by the coefficient, variable name, and index.  Subsequent coefficients
  consist of +/-, coefficient, variable name, and index. To keep things
  simple, if a coefficient doesn't fit on a line, we throw it back for
  reformatting as the initial coefficient of the next line.
*/
  coni = NULL ;
  if (consys_getrow_pk(consys,i,&coni) == FALSE)
  { errmsg(122,rtnnme,consys->nme,
	   "constraint",consys_nme(consys,'c',i,FALSE,NULL),i) ;
    dyio_outfmt(chn,echo,errstring) ;
    if (coni != NULL) pkvec_free(coni) ;
    return ; }

  ai = coni->coeffs ;
  linecnt = 0 ;
  for (ndx = 0 ; ndx < coni->cnt ; ndx++)
  { if (linecnt == 0)
    { charcnt = dyio_outfxd(&buf[0],-60,'l',"\n%s  % g %s(%d)",pfx,ai[ndx].val,
			    consys_nme(consys,'v',ai[ndx].ndx,FALSE,NULL),
			    ai[ndx].ndx) ; }
    else
    { charcnt = dyio_outfxd(&buf[0],-60,'l'," %+g %s(%d)",ai[ndx].val,
			    consys_nme(consys,'v',ai[ndx].ndx,FALSE,NULL),
			    ai[ndx].ndx) ; }
    if (linecnt+charcnt < 70)
    { dyio_outfmt(chn,echo,"%s",&buf[0]) ;
      linecnt += charcnt ; }
    else
    { ndx-- ;
      linecnt = 0 ; } }
  pkvec_free(coni) ;
  
  return ; }

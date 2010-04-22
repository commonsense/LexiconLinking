/*
  This file is part of the support library  for the OsiDylp LP distribution.

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
  This module contains utility print routines which will dump a bnf data
  structure in a (more-or-less) human readable format.
*/

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)bnfrdrio.c	2.3	09/25/04" ;
static char svnid[] UNUSED = "$Id: dylib_bnfrdrio.c 94 2006-06-29 23:06:51Z lou $" ;

#include "dylib_io.h"
#include "dylib_bnfrdr.h"




static const char *prtbnftype (bnftype_enum type)

/*
  This routine returns a pointer to the string representation of the bnf type.

  Parameter:
    type:	a bnf type

  Returns: string representation.
*/

{ static char badtype[30] ;

  switch (type)
  { case bnfG:
      return ("G") ;
    case bnfNP:
      return ("NP") ;
    case bnfP:
      return ("P") ;
    case bnfT:
      return ("T") ;
    case bnfDS:
      return ("DF") ;
    case bnfDL:
      return ("DB") ;
    case bnfRS:
      return ("RF") ;
    case bnfRL:
      return ("RB") ;
    case bnfI:
      return ("I") ;
    case bnfL:
      return ("L") ;
    default:
    { dyio_outfxd(badtype,-(sizeof(badtype)-1),'l',"bad bnf type (%d)",type) ;
      return (badtype) ; } } }


static const char *prtbnfttype (bnfttype_enum ttype)

/*
  This routine returns a pointer to the string representation of the terminal
  type.

  Parameter:
    ttype:	a terminal type

  Returns: string representation.
*/

{ static char badtype[40] ;

  switch (ttype)
  { case bnfttNIL:
      return ("NIL") ;
    case bnfttN:
      return ("N") ;
    case bnfttID:
      return ("ID") ;
    case bnfttD:
      return ("D") ;
    case bnfttF:
      return ("F") ;
    case bnfttQ:
      return ("Q") ;
    default:
    { dyio_outfxd(badtype,-(sizeof(badtype)-1),
		  'l',"bad terminal type (%d)",ttype) ;
      return (badtype) ; } } }



/*
  prtdefname and prtrefname are basically identical, but we need two to avoid
  playing yet more games with unions and casts. A null def or ref parameter
  is an error, and a null name field is most likely an error (indicating
  incorrect initialisation or overwriting of the structure), but valid
  strings are returned in both cases.
*/


static void prtdefname (ioid chn, bool echo, bnfdef_struct *def)

/*
  Parameters:
    chn:	stream id from openfile.
    echo:	TRUE to echo to stdout
    def:	a bnf reference

  Returns: undefined
*/

{ if (def == NULL)
    dyio_outfmt(chn,echo,"<<null pointer>>") ;
  else
  if (def->name == NULL)
    dyio_outfmt(chn,echo,"unnamed(%#08x)",def) ;
  else
    dyio_outfmt(chn,echo,"%s",def->name) ;
    
  return ; }


static void prtrefname (ioid chn, bool echo, bnfref_struct *ref)

/*
  Parameters:
    chn:	stream id from openfile.
    echo:	TRUE to echo to stdout
    ref:	a bnf reference

  Returns: undefined
*/

{ if (ref == NULL)
    dyio_outfmt(chn,echo,"<<null pointer>>") ;
  else
  if (ref->name == NULL)
    dyio_outfmt(chn,echo,"unnamed(%#08x)",ref) ;
  else
    dyio_outfmt(chn,echo,"%s",ref->name) ;
    
  return ; }


static void prtstring (ioid chn, bool echo, const char *string)

/*
  This routine prints a text string, or * if the string is null.

  Parameters:
    chn:	stream id from openfile.
    echo:	TRUE to echo to stdout
    string:	text string
  
  Returns: undefined
*/

{ if (string != NULL)
    dyio_outfmt(chn,echo,"\"%s\"",string) ;
  else
    dyio_outchr(chn,echo,'*') ; }



static void prtstore (ioid chn, bool echo, bnfref_any ref)

/*
  This routine prints the storage directions in a bnf reference.

  Parameters:
    chn:	stream id from openfile.
    ref:	a bnf reference
  
  Returns: undefined
*/

{ if (flgon(ref.t2->uflgs,bnfstore) == TRUE || ref.t2->type == bnfI)
  { if (flgon(ref.t2->uflgs,bnfatsgn) == TRUE) dyio_outchr(chn,echo,'@') ;
    dyio_outfmt(chn,echo,"%d",ref.t2->offset) ; }
  else
    dyio_outchr(chn,echo,'*') ; }


static void prtlst (ioid chn, bool echo, struct bnfref_type3 *ref)

/*
  This routine prints the list information in a bnf reference.

  Parameter:
    chn:	stream id from openfile.
    ref:	a bnf reference
  
  Returns: undefined
*/

{ if (flgon(ref->uflgs,bnflst) == TRUE)
  { if (flgon(ref->uflgs,bnfstbg) == TRUE) dyio_outchr(chn,echo,'@') ;
    prtrefname(chn,echo,(bnfref_struct *) ref) ; }
  else
    dyio_outchr(chn,echo,'*') ; }


static void prtlbl (ioid chn, bool echo,
		    bnflblsrc_enum code, bnfref_struct *src)

/*
  This routine prints the information in the code and src fields of a label
  definition.

  Parameters:
    chn:	stream id from openfile
    code:	value of nmcd or ndcd field.
    src:	value of nmsrc or ndsrc field.
  
  Returns: undefined
*/

{ switch (code)
  { case bnfncBNF:
    { prtrefname(chn,echo,src) ;
      break ; }
    case bnfncS:
    { dyio_outfmt(chn,echo,"%%%d",(int) src) ;
      break ; }
    case bnfncC:
    { dyio_outchr(chn,echo,'c') ;
      break ; }
    case bnfncN:
    { dyio_outchr(chn,echo,'n') ;
      break ; }
    default:
    { dyio_outfmt(chn,echo,"invalid! (%d)",code) ;
      break ; } } }


static void prtlblsav (ioid chn, bool echo, int flg, bnfLBdef_struct *def)

/*
  This routine prints the save information for a label definition.

  Parameters:
    chn:	stream id from openfile.
    flg:	either bnfsvnd or bnfsvnm.
    def:	a label definition
  
  Returns: undefined
*/

{ if (flgon(def->dflgs,flg) == TRUE)
    switch (flg)
    { case bnfsvnd:
      { dyio_outfmt(chn,echo,"%%%d",def->savnd) ;
	break ; }
      case bnfsvnm:
      { dyio_outfmt(chn,echo,"%%%d",def->savnm) ;
	break ; } }
  else
    dyio_outchr(chn,echo,'*') ; }
  
      



void prtbnfref (ioid chn, bool echo, bnfref_struct *bnfref)

/*
  This routine prints a bnf reference data structure in a (sort of) human
  readable format.

  Parameters:
    chn:	stream id, obtained from openfile.
    ref:	bnf reference.
  
  Returns: undefined
*/

{ bnfref_any ref ;

/*
  Check that ref is non-null - print <null!> if it isn't.
*/
  ref.com = bnfref ;
  if (ref.com == NULL)
  { dyio_outfmt(chn,echo,"<null!>") ;
    return ; }
/*
  Print the basics: type, name, and name of associated definition.
*/
  dyio_outfmt(chn,echo,"<%s,",prtbnftype(ref.com->type)) ;
  prtrefname(chn,echo,ref.com) ;
  dyio_outfmt(chn,echo,"->") ;
  prtdefname(chn,echo,ref.com->defn) ;
/*
  And flesh out with the remainder of the type-specific information.
*/
  switch (ref.com->type)
  { case bnfG:
    { dyio_outchr(chn,echo,',') ;
      prtstore(chn,echo,ref) ;
      dyio_outchr(chn,echo,',') ;
      prtlst(chn,echo,ref.t3) ;
      break ; }
    case bnfNP:
    { dyio_outchr(chn,echo,',') ;
      prtlst(chn,echo,ref.t3) ;
      break ; }
    case bnfP:
    { dyio_outchr(chn,echo,',') ;
      if (flgon(ref.P->uflgs,bnfsv) == TRUE)
        dyio_outfmt(chn,echo,"%%%d",ref.P->offset) ;
      else
	prtstore(chn,echo,ref) ;
      dyio_outchr(chn,echo,',') ;
      prtlst(chn,echo,ref.t3) ;
      break ; }
    case bnfT:
    { dyio_outchr(chn,echo,',') ;
      if (flgon(ref.T->uflgs,bnfflt) == TRUE)
        dyio_outfmt(chn,echo,"flt,") ;
      if (flgon(ref.T->uflgs,bnfdbl) == TRUE)
        dyio_outfmt(chn,echo,"dbl,") ;
      if (flgon(ref.T->uflgs,bnfcs) == TRUE)
        dyio_outfmt(chn,echo,"cs,") ;
      if (flgon(ref.T->uflgs,bnfmin) == TRUE)
        dyio_outfmt(chn,echo,"min,") ;
      if (flgon(ref.T->uflgs,bnfexact) == TRUE)
        dyio_outfmt(chn,echo,"exact,") ;
      prtstore(chn,echo,ref) ;
      break ; }
    case bnfI:
    { dyio_outchr(chn,echo,',') ;
      prtstore(chn,echo,ref) ;
      break ; }
    case bnfDS:
    case bnfDL:
    case bnfRS:
    case bnfRL:
    case bnfL:
    { break ; }
    default:
    { break ; } }
  dyio_outchr(chn,echo,'>') ; }



void prtbnfdef (ioid chn, bool echo, bnfdef_struct *bnfdef)

/*
  This routine prints a bnf definition in (sort of) human readable format.

  Parameters:
    chn:	stream id as returned from openfile
    def:	a bnf definition
  
  Returns: undefined
*/

{ bnfdef_any def ;
  bnfref_struct ***altrefs,**comprefs ;
  int altnum,altndx,compnum,compndx ;
  const char *txm1 = "|\n\t" ;

/*
  Suppress compiler warning.
*/
  comprefs = NULL ;
/*
  Check that bnfdef is non-null - print <null!> if it isn't.
*/
  def.com = bnfdef ;
  if (def.com == NULL)
  { dyio_outfmt(chn,echo,"<null!>") ;
    return ; }
/*
  Print the definition header first: type and name.
*/
  dyio_outfmt(chn,echo,"\n<%s,",prtbnftype(def.com->type)) ;
  prtdefname(chn,echo,def.com) ;
  if (def.com->type == bnfG)
  { dyio_outfmt(chn,echo,",%d,",def.G->size) ;
    if (def.G->link >= 0)
      dyio_outfmt(chn,echo,"%d",def.G->link) ;
    else
      dyio_outchr(chn,echo,'*') ; }
  dyio_outfmt(chn,echo,"> ::= ") ;
/*
  And now the body of the definition.
*/
  switch (def.com->type)
  { case bnfG:
    { comprefs = def.G->comps ;
      if (comprefs != NULL)
        compnum = (int) *comprefs++ ;
      else
	compnum = 0 ;
      for (compndx = 0 ; compndx < compnum ; compndx++)
	prtbnfref(chn,echo,*comprefs++) ;
      break ; }
    case bnfNP:
    { altrefs = def.NP->alts ;
      if (altrefs != NULL)
      { altnum = (int) *altrefs++ ;
	comprefs = *altrefs++ ; }
      else
	altnum = 0 ;
      for (altndx = 0 ; altndx < altnum ; altndx++, comprefs = *altrefs++)
      { if (comprefs != NULL)
	{ compnum = (int) *comprefs++ ;
	  for (compndx = 0 ; compndx < compnum ; compndx++)
	    prtbnfref(chn,echo,*comprefs++) ; }
	else
	  dyio_outfmt(chn,echo,"<null alternative! (%d)>",altndx+1) ;
	if (altndx < altnum-1) dyio_outfmt(chn,echo,txm1) ; }
      break ; }
    case bnfP:
    { altrefs = def.P->alts ;
      if (altrefs != NULL)
      { altnum = (int) *altrefs++ ;
	comprefs = *altrefs++ ; }
      else
	altnum = 0 ;
      for (altndx = 0 ; altndx < altnum ; altndx++, comprefs = *altrefs++)
      { if (comprefs != NULL)
	{ compnum = (int) *comprefs++ ;
	  for (compndx = 0 ; compndx < compnum ; compndx++)
	    prtbnfref(chn,echo,*comprefs++) ; }
	else
	  dyio_outfmt(chn,echo,"<null alternative! (%d)>",altndx+1) ;
	if (altndx < altnum-1) dyio_outfmt(chn,echo,txm1) ; }
      break ; }
    case bnfT:
    { dyio_outfmt(chn,echo,"<%s",prtbnfttype(def.T->ttype)) ;
      switch (def.T->ttype)
      { case bnfttNIL:
	{ break ; }
	case bnfttN:
	{ dyio_outfmt(chn,echo,"(%d),",def.T->parm1) ;
	  prtstring(chn,echo,def.T->val) ;
	  break ; }
	case bnfttID:
	{ prtstring(chn,echo,def.T->val) ;
	  break ; }
	case bnfttD:
	{ prtstring(chn,echo,def.T->val) ;
	  break ; }
	case bnfttF:
	{ dyio_outfmt(chn,echo,"(%d),",def.T->parm1) ;
	  prtstring(chn,echo,def.T->val) ;
	  break ; }
	case bnfttQ:
	{ dyio_outchr(chn,echo,'(') ;
	  if (def.T->qschr < ' ')
	    dyio_outfmt(chn,echo,"%#02x",def.T->qschr) ;
	  else
	    dyio_outchr(chn,echo,def.T->qschr) ;
	  dyio_outchr(chn,echo,',') ;
	  if (def.T->qechr < ' ')
	    dyio_outfmt(chn,echo,"%#02x",def.T->qechr) ;
	  else
	    dyio_outchr(chn,echo,def.T->qechr) ;
	  dyio_outfmt(chn,echo,"),") ;
	  prtstring(chn,echo,def.T->val) ;
	  break ; } }
      dyio_outchr(chn,echo,'>') ;
      break ; }
    case bnfDS:
    case bnfDL:
    { dyio_outchr(chn,echo,'<') ;
      if (flgon(def.LB->dflgs,bnfadv) == TRUE)
	dyio_outfmt(chn,echo,"a,") ;
      else
	dyio_outfmt(chn,echo,"b,") ;
      prtlbl(chn,echo,def.LB->nmcd,def.LB->nmsrc) ;
      dyio_outchr(chn,echo,',') ;
      prtlblsav(chn,echo,bnfsvnm,def.LB) ;
      dyio_outchr(chn,echo,',') ;
      prtlbl(chn,echo,def.LB->ndcd,def.LB->ndsrc) ;
      dyio_outfmt(chn,echo,"(%d)",def.LB->offset) ;
      dyio_outchr(chn,echo,',') ;
      prtlblsav(chn,echo,bnfsvnd,def.LB) ;
      dyio_outchr(chn,echo,'>') ;
      break ; }
    case bnfRS:
    case bnfRL:
    { dyio_outchr(chn,echo,'<') ;
      if (flgon(def.LB->dflgs,bnfadv) == TRUE)
	dyio_outfmt(chn,echo,"a,") ;
      else
	dyio_outfmt(chn,echo,"b,") ;
      prtlbl(chn,echo,def.LB->nmcd,def.LB->nmsrc) ;
      dyio_outfmt(chn,echo,"(%d)",def.LB->offset) ;
      dyio_outchr(chn,echo,',') ;
      prtlblsav(chn,echo,bnfsvnm,def.LB) ;
      dyio_outchr(chn,echo,',') ;
      prtlbl(chn,echo,def.LB->ndcd,def.LB->ndsrc) ;
      dyio_outfmt(chn,echo,"(%d)",def.LB->offset2) ;
      dyio_outchr(chn,echo,',') ;
      prtlblsav(chn,echo,bnfsvnd,def.LB) ;
      dyio_outchr(chn,echo,'>') ;
      break ; }
    case bnfI:
    { dyio_outfmt(chn,echo,"<%d>",def.I->ival) ;
      break ; }
    case bnfL:
    { dyio_outchr(chn,echo,'<') ;
      if (flgon(def.L->dflgs,bnfsvnm) == TRUE)
	dyio_outfmt(chn,echo,"%%%d",(int) def.L->txt) ;
      else
	prtstring(chn,echo,def.L->txt) ;
      dyio_outchr(chn,echo,'>') ;
      break ; } } }

#ifndef DYLP_NDEBUG


void printtab (ioid chn, bool echo,
	       int nestlvl, bool numlvl, bool tablvl)

/*
  Utility routine to print optional line numbers and tab leadering for
  debugging trace lines.

  Parameters:
    chn:	i/o id for trace output
    echo:	TRUE to echo trace to stdout
    nestlvl:	nesting level of the trace line
    numlvl:	TRUE to print the nesting level at the left of the line
    tablvl:	TRUE to show nesting level with indentation
*/

{ int tabcnt ;

  if (numlvl == TRUE)
  { dyio_outfmt(chn,echo,"%2d: ",nestlvl) ;
    nestlvl-- ; }
/*
  Change the field with in the dyio_outfmt statement to change the amount of
  indentation per level.
*/
  if (tablvl == TRUE)
  { for (tabcnt = nestlvl ; tabcnt > 0 ; tabcnt--)
    { dyio_outfmt(chn,echo,"%4s"," ") ; } }

  return ; }

#endif /* DYLP_NDEBUG */

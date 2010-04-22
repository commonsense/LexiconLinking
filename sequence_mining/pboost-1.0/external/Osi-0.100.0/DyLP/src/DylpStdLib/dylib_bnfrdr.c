/*
  This file is part of the support library for the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This module contains a set of routines which will parse input text and
  construct a data structure according to a bnf-like specification describing
  the input syntax and the data structure to be built.

  The routines do a top-down parse. The major deficiencies in the package at
  the present are that it cannot do arrays and it cannot back up over
  constructs that modify the user's data structure or generate labels. These
  seem restrictive at first glance but in practice have not been enough of a
  problem for me to expend the effort to fix them.
*/

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)bnfrdr.c	3.9	09/25/04" ;
static char svnid[] UNUSED = "$Id: dylib_bnfrdr.c 245 2008-07-08 13:40:02Z lou $" ;

#include "dylib_io.h"
#include "dylib_errs.h"
#include "dylib_strrtns.h"
#include "dylib_bnfrdr.h"

/*
  curnde is always the node created by the currently active generator
  reference; cndesze is its size. newnde is the node returned by the most
  recent successfully parsed generator reference. curtxt is the text string
  being formed by the currently active primitive, newtxt is the most recently
  formed text string.

  NOTE: while dolist is active, it manipulates curnde and curtxt so that they
	are the previous list element and the cumulative text collected by the
	list. This is transparent on either side of dolist and need not be
	taken into account when writing bnfs.
*/

static void *curnde,*newnde ;
static int cndesze ;
static char *curtxt,*newtxt ;

/*
  bnfchn is the input channel ; it is global to the module and set by parse so
  that we don't have to pass it around as a parameter.
*/

static ioid bnfchn ;

/*
  savedtxt is a scratchpad array for keeping around useful text strings
  constructed while parsing the input.
*/

#define SAVE_SLOTS 10
static const char *savedtxt[SAVE_SLOTS+1] ;

#ifndef DYLP_NDEBUG

/*
  Debugging aids  -- compiled in and ready for use unless DYLP_NDEBUG is defined.
  See comments on activating debugging prints in bnfrdr.h.

  nestlvl is the nesting level of bnf rule; if set to -1, no tabs are used.
  numlvl causes the nesting level to be printed at the left margin.
  tablvl causes indentation proportional to the nesting level.
  warnzlbl causes a warning to be issued if a socket/label evaluates to 0.

  printtab is over in bnfrdrio.c; it does indentation of the debug trace.
*/

extern void printtab(ioid dbgchn, bool dbgecho,
		     int nestlvl, bool numlvl, bool tablvl) ;

  static int debug = 0 ;
  static ioid dbgchn = IOID_NOSTRM ;
  static bool dbgecho = TRUE ;
  static int nestlvl ;
  static bool warnzlbl = TRUE, numlvl = TRUE, tablvl = TRUE ;

#endif /* DYLP_NDEBUG */




/*
  The reader maintains lists of sockets and labels, one list for defined and
  one list for undefined of each type. Entries are made on the defined lists
  by all socket and label definitions. Entries are made on the undefined
  lists by any socket or label reference where one of the socket or label is
  undefined.

  For entries on the defined lists, the fields are interpreted as follows:

  Field		Description
  -----		-----------
  lblnxt	Next entry in this list.
  lblnmtxt	The name of the socket/label.
  lblval	The value of the socket/label.
  lbladv	TRUE if this socket/label has been advanced (redefined),
		FALSE otherwise.

  For entries on the undefined lists, the fields are interpreted as
  follows:

  Field		Description
  -----		-----------
  lblnxt	Next entry in this list.
  lblnmtxt	The name of the socket/label.
  lblval	The value of the defined label/socket from the reference.

*/

typedef struct deflbl_struct_tag { struct deflbl_struct_tag *lblnxt ;
				   const char *lblnmtxt ;
				   void *lblval ;
				   bool lbladv ; } deflbl_struct ;

typedef struct udeflbl_struct_tag { struct udeflbl_struct_tag *lblnxt ;
				    const char *lblnmtxt ;
				    void *lblval ; } udeflbl_struct ;

/*
  The list heads.
*/

static deflbl_struct *blbllst,*flbllst ;
static udeflbl_struct *ublbllst,*uflbllst ;



/*
  Routines to access the saved text array. The major reason these are routines
  is to keep the error checking and string management in one place. They're
  just big enough that I don't want to insert the code inline each time.
*/

static void strenter (int ndx, const char *txt)

/*
  This routine is used to insert text into the saved text array.

  Parameters:
    ndx:	index in the saved text array
    txt:	string to be saved
  
  Returns: undefined
*/

{ const char *rtnnme = "strenter" ;

/*
  Check out the parameters.
*/
  if (ndx < 0 || ndx > SAVE_SLOTS)
  { errmsg(40,rtnnme,ndx,SAVE_SLOTS) ;
    return ; }
  if (txt == NULL)
  { errmsg(2,rtnnme,"txt") ;
    return ; }
/*
  Clean out the slot, if necessary, and store the new pointer.
*/
  if (savedtxt[ndx] != NULL) STRFREE(savedtxt[ndx]) ;
  savedtxt[ndx] = STRALLOC(txt) ;

  return ; }


static const char *strretrv (int ndx)

/*
  This routine is used to retrieve text from the stored text array.

  Parameter:
    ndx:	index in the saved text array.

  Returns: The character string stored in savedtxt[ndx], or NULL if there is a
	   problem.
*/

{ const char *rtnnme = "strretrv" ;

/*
  Check out the parameter.
*/
  if (ndx < 0 || ndx > SAVE_SLOTS)
  { errmsg(40,rtnnme,ndx,SAVE_SLOTS) ;
    return (NULL) ; }
/*
  Fetch back the text.
*/
  if (savedtxt[ndx] == NULL) errmsg(51,rtnnme,ndx) ;

  return (savedtxt[ndx]) ; }




/*
  A few useful macros. offset_in_range takes an offset and a size and checks
  whether the offset specified is in range for storing an object of the given
  size in curnde. make_label generates a pointer, and make_socket makes a
  socket (pointer to pointer) for storing things.
*/

#define offset_in_range(qq_off,qq_sze) \
	(((qq_off) >= 0 && (qq_off) <= cndesze - (qq_sze))?TRUE:FALSE)

#define max_offset(qq_sze) \
	(cndesze - (qq_sze))

#define make_label(qq_base,qq_offset) \
	((void *) ((ptrdiff_t) (qq_base) + (qq_offset)))

#define make_socket(qq_base,qq_offset) \
	((void **) ((ptrdiff_t) (qq_base) + (qq_offset)))



static deflbl_struct *finddlbl (deflbl_struct **lst, const char *txt)

/*
  This routine searches for the label whose name is txt in the defined label
  list headed by lst. It does a simple linear search. It assumes that txt has
  not necessarily been entered in the literal tree, hence an actual string
  comparison is performed (case-independent).

  Parameters:
    lst:	head of the label list
    txt:	name of the label

  Returns: The entry for the label, or NULL if none is found.
*/

{ deflbl_struct *lblent ;
  const char *rtnnme = "finddlbl" ;

/*
  Check out the parameters.
*/
  if (lst == NULL)
  { errmsg(2,rtnnme,"label list") ;
    return (NULL) ; }
  if (txt == NULL)
  { errmsg(2,rtnnme,"label name") ;
    return (NULL) ; }
/*
  Search the list for the named entry.
*/
  for (lblent = *lst ; lblent != NULL ; lblent = lblent->lblnxt)
    if (cistrcmp(lblent->lblnmtxt,txt) == 0) return (lblent) ;

  return (NULL) ; }



static void lblresolve (udeflbl_struct **lst, deflbl_struct *dlbl)

/*
  This routine searches an undefined label list for entries which can be
  resolved by the (presumably newly) defined label dlbl. The routine assumes
  that the text strings it considers will have been entered in the literal tree
  and thus compares pointers to test for equality.

  Parameters:
    lst:	one of the undefined label lists.
    dlbl:	a defined label
  
  Returns: undefined
*/

{ udeflbl_struct **ulbl2,*ulbl1 ;
  const char *nmtxt ;
  void **socket ;
  const char *rtnnme = "lblresolve" ;

/*
  Check out the parameters.
*/
  if (lst == NULL)
  { errmsg(2,rtnnme,"label list") ;
    return ; }
  if (dlbl == NULL)
  { errmsg(2,rtnnme,"defined label") ;
    return ; }
/*
  Search through the list, resolving any references we find.
*/
  nmtxt = dlbl->lblnmtxt ;
  for (ulbl2 = lst, ulbl1 = *ulbl2 ;
       ulbl1 != NULL ;
       ulbl2 = &ulbl1->lblnxt, ulbl1 = *ulbl2)
    if (ulbl1->lblnmtxt == nmtxt)
    { *ulbl2 = ulbl1->lblnxt ;
      if (lst == &uflbllst)
      { socket = make_socket(dlbl->lblval,0) ;
	*socket = ulbl1->lblval ; }
      else
      { socket = make_socket(ulbl1->lblval,0) ;
	*socket = dlbl->lblval ;
      FREE((char *) ulbl1) ; } }

  return ; }



void rdrinit (void)

/*
  This routine initializes the bnf reader.

  Parameters:
    zerolbl:	set to TRUE if the reader should issue a warning when a
		label evaluates to 0

  Returns: undefined
*/

{ int ndx ;

  curnde = NULL ;
  cndesze = 0 ;
  newnde = NULL ;
  curtxt = NULL ;
  newtxt = NULL ;
  for (ndx = 0 ; ndx < SAVE_SLOTS ; ndx++) savedtxt[ndx] = NULL ;
  flbllst = NULL ;
  blbllst = NULL ;
  uflbllst = NULL ;
  ublbllst = NULL ;
  
  return ; }



void rdrclear (void)

/*
  This routine clears the bnf reader by releasing the label lists and savedtxt
  array, and newtxt. It also makes two consistency checks: At the end of the
  parse, it should be true that curnde and curtxt are null.

  Parameters: none.

  Returns: undefined.
*/

{ int ndx ;
  deflbl_struct *dlbl,*nxtdlbl ;
  udeflbl_struct *udlbl,*nxtudlbl ;
  const char *rtnnme = "rdrclear" ;

  for (dlbl = blbllst ; dlbl != NULL ; dlbl = nxtdlbl)
  { nxtdlbl = dlbl->lblnxt ;
    STRFREE(dlbl->lblnmtxt) ;
    FREE((char *) dlbl) ; }
  blbllst = NULL ;
  for (dlbl = flbllst ; dlbl != NULL ; dlbl = nxtdlbl)
  { nxtdlbl = dlbl->lblnxt ;
    STRFREE(dlbl->lblnmtxt) ;
    FREE((char *) dlbl) ; }
  flbllst = NULL ;
  for (udlbl = ublbllst ; udlbl != NULL ; udlbl = nxtudlbl)
  { nxtudlbl = udlbl->lblnxt ;
    STRFREE(udlbl->lblnmtxt) ;
    FREE((char *) udlbl) ; }
  ublbllst = NULL ;
  for (udlbl = uflbllst ; udlbl != NULL ; udlbl = nxtudlbl)
  { nxtudlbl = udlbl->lblnxt ;
    STRFREE(udlbl->lblnmtxt) ;
    FREE((char *) udlbl) ; }
  uflbllst = NULL ;

  for (ndx = 0 ; ndx < SAVE_SLOTS ; ndx++)
  { if (savedtxt[ndx] != NULL) STRFREE(savedtxt[ndx]) ;
    savedtxt[ndx] = NULL ; }

  if (newtxt != NULL)
  { FREE(newtxt) ;
    newtxt = NULL ; }

  if (curnde != NULL) warn(71,rtnnme,"curnde") ;
  if (curtxt != NULL) warn(71,rtnnme,"curtxt") ;
  
  return ; }



/*
  The heart of the package, the routines which actually deal with the bnf
  constructs and build the user's data strucure.  There is one routine to
  handle each type of bnf construct, plus a routine which handles lists.
*/

static bool dogenerator(bnfGref_struct *ref),
	    dononprimitive(bnfNPref_struct *ref),
	    doprimitive(bnfNPref_struct *ref),
	    doterminal(bnfTref_struct *ref),
	    doimmediate(bnfIref_struct *ref),
	    doliteral(bnfLref_struct *ref),
	    dolabel(bnfLBref_struct *ref),
	    doreference(bnfLBref_struct *ref),
	    dolist(bnfref_any ref) ;



static bool dogenerator (bnfGref_struct *ref)

/*
  This routine handles constructs of type generator. It checks to make sure
  that the reference and the definition are well-formed, then constructs the
  node as specified in the definition. It then enters a loop to parse the
  components of the body of the definition. Failure of the parse results in
  an error return.

  Parameter:
    ref:	reference to a generator definition

  Result: TRUE if the parse succeeds, FALSE if it fails.
*/

{ void *savcnde,**socket ;
  int savcndesze,compnum,compndx ;
  bnfGdef_struct *def ;
  bnfref_struct **comprefs ;
  bnfref_any compref ;
  bool success ;
  const char *rtnnme = "dogenerator" ;

/*
  First test the reference to make sure that it references a generator.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfGdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfG)
  { errmsg(36,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref->uflgs,bnfdebug))
    if (debug++ == 0)
    { dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
      nestlvl = 0 ; }
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    dyio_outfmt(dbgchn,dbgecho," ::=\n") ; }
# endif

/*
  Now look at the storage spec. If we're not in a list, we check that any
  storage request is within curnde. It's an error if curnde is null. If
  the spec passes, we form the socket address. If we're in a list, then we
  assume that dolist is handling the linking and no storage is done here.
*/
  if (flgon(ref->uflgs,bnfstore) == TRUE && flgoff(ref->uflgs,bnflst) == TRUE)
  { if (curnde == NULL)
    { errmsg(68,rtnnme) ;
      return (FALSE) ; }
    if (offset_in_range(ref->offset,sizeof(void *)) == FALSE)
    { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(void *))) ;
      return (FALSE) ; }
    socket = make_socket(curnde,ref->offset) ; }
  else
    socket = NULL ;
/*
  Now check that the definition has a valid size. If it passes, stack the old
  curnde and make the new one. Clearing the new node is not just being nice -
  it allows us to sidestep the issue of just what sort of NIL thing we're
  dealing with down in doterminal.
*/
  if (def->size <= 0)
  { errmsg(31,rtnnme,def->size) ;
    return (FALSE) ; }
  savcnde = curnde ;
  savcndesze = cndesze ;
  cndesze = def->size ;
  curnde = (void *) MALLOC(cndesze) ;
  memset((void *) curnde,0,cndesze) ;
/*
  Set up for the parse loop by getting the component array and extracting the
  number of components.
*/
  comprefs = def->comps ;
  if (comprefs != NULL)
  { compnum = addrToInt(*comprefs++) ;
    compref.com = *comprefs++ ; }
  else
    compnum = 0 ;
/*
  The main loop to handle the components. 
*/
  for (compndx = 0, success = TRUE ;
       compndx < compnum && success == TRUE ;
       compndx++, compref.com = *comprefs++)
  { if (compref.com == NULL)
    { errmsg(32,rtnnme,compndx+1,compnum) ;
      success = FALSE ;
      break ; }
    switch (compref.com->type)
    { case bnfG:
      { if (flgon(compref.G->uflgs,bnflst) == TRUE)
	  success = dolist(compref) ;
	else
	  success = dogenerator(compref.G) ;
	break ; }
      case bnfNP:
      { if (flgon(compref.NP->uflgs,bnflst) == TRUE)
	  success = dolist(compref) ;
	else
	  success = dononprimitive(compref.NP) ;
	break ; }
      case bnfP:
      { if (flgon(compref.P->uflgs,bnflst) == TRUE)
	  success = dolist(compref) ;
	else
	  success = doprimitive(compref.P) ;
	break ; }
      case bnfT:
      { success = doterminal(compref.T) ;
	break ; }
      case bnfI:
      { success = doimmediate(compref.I) ;
	break ; }
      case bnfL:
      { errmsg(34,rtnnme,"literal","generator") ;
	success = FALSE ;
	break ; }
      case bnfDS:
      case bnfDL:
      { success = dolabel(compref.LB) ;
	break ; }
      case bnfRS:
      case bnfRL:
      { success = doreference(compref.LB) ;
	break ; }
      default:
      { errmsg(35,rtnnme,compref.com->type) ;
	success = FALSE ;
	break ; } } }

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { if (success == FALSE)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      dyio_outfmt(dbgchn,dbgecho,"-- fail @ %d of %d --",compndx,compnum) ; }
    dyio_outchr(dbgchn,dbgecho,'\n') ;
    nestlvl-- ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
        dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

/*
  If we've successfully parsed the components of the generator definition,
  update newnde and store a pointer to the new node if requested.
*/
  if (success == TRUE)
  { if (socket != NULL) *socket = curnde ;
    newnde = curnde ; }
/*
  Restore the old curnde and then return.
*/
  curnde = savcnde ;
  cndesze = savcndesze ;

  return (success) ; }



static bool dononprimitive (bnfNPref_struct *ref)

/*
  This routine handles constructs of type non-primitive. It first checks to
  make sure that the reference is well-formed. Next, it creates a copy of
  curnde and marks the current position in the input stream for error recovery
  purposes. The actual parsing is handled by two nested loops. The outer loop
  sequences through the alternative parses, while the inner loop parses the
  components of each alternative.

  Parameter:
    ref:	reference to a non-primitive definition

  Returns: TRUE if the non-primitive is successfully parsed, FALSE otherwise.
*/

{ void *savcnde = NULL ; 
  bnfNPdef_struct *def ;
  bnfref_struct ***altrefs,**comprefs ;
  bnfref_any compref ;
  int altnum,altndx,compnum,compndx ;
  long marker ;
  bool success ;
  const char *rtnnme = "dononprimitive" ;

/*
  First test the reference to make sure that it references a non-primitive.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfNPdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfNP)
  { errmsg(38,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref->uflgs,bnfdebug))
    if (debug++ == 0)
    { dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
      nestlvl = 0 ; }
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    dyio_outfmt(dbgchn,dbgecho," ::=\n") ; }
# endif

/*
  Prepare for failure of an alternative parse. Mark the input and create a copy
  of curnde, so we can recover if a parse fails.
*/
  marker = dyio_mark(bnfchn) ;
  if (curnde != NULL)
  { if (cndesze <= 0)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; }
    savcnde = (void *) MALLOC(cndesze) ;
    memcpy((void *) savcnde,(void *) curnde,cndesze) ; }
/*
  Set up for the loop to do alternative parses by getting the alternative array
  and extracting the number of alternatives.
*/
  altrefs = def->alts ;
  if (altrefs != NULL)
  { altnum = addrToInt(*altrefs++) ;
    comprefs = *altrefs++ ; }
  else
  { altnum = 0 ;
    comprefs = NULL ; }
/*
  The outer loop to handle the alternatives. 
*/
  for (altndx = 0 ;
       altndx < altnum ;
       altndx++, comprefs = *altrefs++)
/*
  Set up for the parse loop by getting the component array and extracting the
  number of components. It's an error if an alternative has no components.
*/
  { 

#   ifndef DYLP_NDEBUG
    if (debug > 0)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      dyio_outfmt(dbgchn,dbgecho,
		  "[ alternative %d of %d ]\n",altndx+1,altnum) ; }
#   endif

    if (comprefs == NULL)
    { errmsg(37,rtnnme,altndx+1,altnum) ;
      break ; }
    else
    { compnum = addrToInt(*comprefs++) ;
      compref.com = *comprefs++ ; }
/*
  The inner loop to handle the components of an alternative. 
*/
    for (compndx = 0, success = TRUE ;
	 compndx < compnum && success == TRUE ;
	 compndx++, compref.com = *comprefs++)
    { if (compref.com == NULL)
      { errmsg(32,rtnnme,compndx+1,compnum) ;
	success = FALSE ;
	break ; }
      switch (compref.com->type)
      { case bnfG:
	{ if (flgon(compref.G->uflgs,bnflst) == TRUE)
	    success = dolist(compref) ;
	  else
	    success = dogenerator(compref.G) ;
	  break ; }
	case bnfNP:
	{ if (flgon(compref.NP->uflgs,bnflst) == TRUE)
	    success = dolist(compref) ;
	  else
	    success = dononprimitive(compref.NP) ;
	  break ; }
	case bnfP:
	{ if (flgon(compref.P->uflgs,bnflst) == TRUE)
	    success = dolist(compref) ;
	  else
	    success = doprimitive(compref.P) ;
	  break ; }
	case bnfT:
	{ success = doterminal(compref.T) ;
	  break ; }
	case bnfI:
	{ success = doimmediate(compref.I) ;
	  break ; }
	case bnfL:
	{ errmsg(34,rtnnme,"literal","non-primitive") ;
	  success = FALSE ;
	  break ; }
	case bnfDS:
	case bnfDL:
	{ success = dolabel(compref.LB) ;
	  break ; }
	case bnfRS:
	case bnfRL:
	{ success = doreference(compref.LB) ;
	  break ; }
	default:
	{ errmsg(35,rtnnme,compref.com->type) ;
	  success = FALSE ;
	  break ; } } }

#   ifndef DYLP_NDEBUG
    if (debug > 0)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      if (success == TRUE)
      { dyio_outfmt(dbgchn,dbgecho,"[ pass %d ]",altndx+1) ;
	nestlvl-- ;
	if (flgon(ref->uflgs,bnfdebug))
	  if (--debug == 0)
	    dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
      else
        dyio_outfmt(dbgchn,dbgecho,
		    "[ fail %d @ %d of %d ]",altndx+1,compndx,compnum) ;
      dyio_outchr(dbgchn,dbgecho,'\n') ; }
#   endif

/*
  End loop to process the components of an alternative. If we've successfully
  parsed the alternative, clean up and return. Otherwise, back up the input
  and restore curnde.
*/
    if (success == TRUE)
    { if (curnde != NULL) FREE((char *) savcnde) ;
      return (TRUE) ; }
    else
    { if (curnde != NULL) memcpy((void *) curnde,(void *) savcnde,cndesze) ;
      dyio_backup(bnfchn,marker) ; } }
/*
  To get here, all the alternative parses have failed. Clean up and return
  failure.
*/
  if (curnde != NULL)
  { memcpy((void *) curnde,(void *) savcnde,cndesze) ;
    FREE(savcnde) ; }
# ifndef DYLP_NDEBUG
  if (debug > 0)
  { nestlvl-- ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif
  
  return (FALSE) ; }



static bool doprimitive (bnfNPref_struct *ref)

/*
  This routine handles constructs of type primitive. It first checks to make
  sure that the reference is well-formed. Next, it creates a copy of curnde
  and marks the current position in the input stream for error recovery
  purposes. The actual parsing is handled by two nested loops. The outer loop
  sequences through the alternative parses, while the inner loop parses the
  components of each alternative. As each alternative is parsed, a text string
  is constructed which is the concatenation of the text strings returned by
  all terminals, literals, and primitives parsed for that alternative. If the
  parse is successful, a pointer to the string may be stored as indicated in
  the reference.

  Parameter:
    ref:	reference to a primitive definition

  Returns: TRUE if the primitive is successfully parsed, FALSE otherwise.
*/

{ void *savcnde = NULL ;
  void **socket ;
  char *savctxt,*lcltxt ;
  bnfPdef_struct *def ;
  bnfref_struct ***altrefs,**comprefs ;
  bnfref_any compref ;
  int altnum,altndx,compnum,compndx ;
  long marker ;
  bool success ;
  const char *rtnnme = "doprimitive" ;

/*
  First test the reference to make sure that it references a primitive.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfPdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfP)
  { errmsg(39,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref->uflgs,bnfdebug))
    if (debug++ == 0)
    { dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
      nestlvl = 0 ; }
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    dyio_outfmt(dbgchn,dbgecho," ::=\n") ; }
# endif

/*
  Prepare for failure of an alternative parse. Mark the input and create a copy
  of curnde, so we can recover if a parse fails. Prepare for character string
  collection by stacking curtxt and creating a null string. curtxt is left NULL
  here and only made non-null while we're dealing with a terminal, literal, or
  primitive in the parse loop.
*/
  marker = dyio_mark(bnfchn) ;
  if (curnde != NULL)
  { if (cndesze <= 0)
    { errmsg(1,rtnnme,__LINE__) ;
      return (FALSE) ; }
    savcnde = (void *) MALLOC(cndesze) ;
    memcpy((void *) savcnde,(void *) curnde,cndesze) ; }
  savctxt = curtxt ;
  curtxt = NULL ;
  lcltxt = (char *) MALLOC(sizeof(char)) ;
  *lcltxt = '\0' ;
/*
  Set up for the loop to do alternative parses by getting the alternative array
  and extracting the number of alternatives.
*/
  altrefs = def->alts ;
  if (altrefs != NULL)
  { altnum = addrToInt(*altrefs++) ;
    comprefs = *altrefs++ ; }
  else
  { altnum = 0 ;
    comprefs = NULL ; }
/*
  The outer loop to handle the alternatives. 
*/
  for (altndx = 0 ;
       altndx < altnum ;
       altndx++, comprefs = *altrefs++)
/*
  Set up for the parse loop by getting the component array and extracting the
  number of components. It's an error if an alternative has no components.
*/
  {

#   ifndef DYLP_NDEBUG
    if (debug > 0)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      dyio_outfmt(dbgchn,dbgecho,
		  "[ alternative %d of %d ]\n",altndx+1,altnum) ; }
#   endif

    if (comprefs == NULL)
    { errmsg(37,rtnnme,altndx+1,altnum) ;
      break ; }
    else
    { compnum = addrToInt(*comprefs++) ;
      compref.com = *comprefs++ ; }
/*
  The inner loop to handle the components of an alternative. 
*/
    for (compndx = 0, success = TRUE ;
	 compndx < compnum && success == TRUE ;
	 compndx++, compref.com = *comprefs++)
    { if (compref.com == NULL)
      { errmsg(32,rtnnme,compndx+1,compnum) ;
	success = FALSE ;
	break ; }
      switch (compref.com->type)
      { case bnfG:
	{ if (flgon(compref.G->uflgs,bnflst) == TRUE)
	    success = dolist(compref) ;
	  else
	    success = dogenerator(compref.G) ;
	  break ; }
	case bnfNP:
	{ if (flgon(compref.NP->uflgs,bnflst) == TRUE)
	    success = dolist(compref) ;
	  else
	    success = dononprimitive(compref.NP) ;
	  break ; }
	case bnfP:
	{ curtxt = lcltxt ;
	  if (flgon(compref.P->uflgs,bnflst) == TRUE)
	    success = dolist(compref) ;
	  else
	    success = doprimitive(compref.P) ;
	  lcltxt = curtxt ;
	  curtxt = NULL ;
	  break ; }
	case bnfT:
	{ curtxt = lcltxt ;
	  success = doterminal(compref.T) ;
	  lcltxt = curtxt ;
	  curtxt = NULL ;
	  break ; }
	case bnfI:
	{ success = doimmediate(compref.I) ;
	  break ; }
	case bnfL:
	{ curtxt = lcltxt ;
	  success = doliteral(compref.L) ;
	  lcltxt = curtxt ;
	  curtxt = NULL ;
	  break ; }
	case bnfDS:
	case bnfDL:
	{ success = dolabel(compref.LB) ;
	  break ; }
	case bnfRS:
	case bnfRL:
	{ success = doreference(compref.LB) ;
	  break ; }
	default:
	{ errmsg(35,rtnnme,compref.com->type) ;
	  success = FALSE ;
	  break ; } } }

#   ifndef DYLP_NDEBUG
    if (debug > 0)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      if (success == TRUE)
      { dyio_outfmt(dbgchn,dbgecho,"[ pass %d ]",altndx+1) ;
	nestlvl-- ;
	if (flgon(ref->uflgs,bnfdebug))
	  if (--debug == 0)
	    dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
      else
        dyio_outfmt(dbgchn,dbgecho,
		    "[ fail %d @ %d of %d ]",altndx+1,compndx,compnum) ;
      dyio_outchr(dbgchn,dbgecho,'\n') ; }
#   endif

/*
  End loop to process the components of an alternative. If we've successfully
  parsed the alternative, there are three actions to do. First, set newtxt to
  the collected text and restore the old curtxt. Next, store the new text as
  indicated in the reference. Finally, if curtxt is non-null, concatenate
  newtxt to it. 
  For an unsuccussful parse, restore curnde, throw away the collected text,
  back up the input, and try again.
*/
  if (success == TRUE)
  { if (newtxt != NULL) FREE(newtxt) ;
    newtxt = lcltxt ;
    curtxt = savctxt ;
    if (flgon(ref->uflgs,bnfstore) == TRUE &&
	flgoff(ref->uflgs,bnflst) == TRUE)
    { if (flgon(ref->uflgs,bnfsv) == TRUE)
      { strenter(ref->offset,newtxt) ; }
      else
      { if (curnde == NULL)
	{ errmsg(68,rtnnme) ;
	  newtxt = NULL ;
	  break ; }
        socket = make_socket(curnde,ref->offset) ;
	if (flgon(ref->uflgs,bnfatsgn) == TRUE)
	{ if (offset_in_range(ref->offset,sizeof(char *)) == FALSE)
	  { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(char *))) ;
	    newtxt = NULL ;
	    break ; }
	  *((const char **) socket) = STRALLOC(newtxt) ; }
	else
	{ if (offset_in_range(ref->offset,strlen(newtxt)+1) == FALSE)
	  { errmsg(30,rtnnme,ref->offset,max_offset(strlen(newtxt)+1)) ;
	    newtxt = NULL ;
	    break ; }
	  (void) strcpy((char *) socket,newtxt) ; } } }
    if (curtxt != NULL && strlen(newtxt) > 0)
    { savctxt = (char *) MALLOC(strlen(curtxt)+strlen(newtxt)+1) ;
      (void) strcpy(savctxt,curtxt) ;
      (void) strcat(savctxt,newtxt) ;
      FREE(curtxt) ;
      curtxt = savctxt ; }
    if (curnde != NULL) FREE((char *) savcnde) ;
    return (TRUE) ; }
  else
  { if (curnde != NULL) memcpy((void *) curnde,(void *) savcnde,cndesze) ;
    *lcltxt = '\0' ;
    dyio_backup(bnfchn,marker) ; } }
/*
  To get here, all the alternative parses have failed or an error was detected
  which caused the alternatives loop to abort. Clean up and return
  failure.
*/
  if (curnde != NULL)
  { memcpy((void *) curnde,(void *) savcnde,cndesze) ;
    FREE(savcnde) ; }
  FREE(lcltxt) ;
  curtxt = savctxt ;

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { nestlvl-- ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif
  
  return (FALSE) ; }



static int scanbinary (const char *string, int *intp)

/*
  Quick and dirty routine to convert a binary number from string form to
  internal representation. Has the same behaviour as sscanf.

  Parameters:
    string:	string to be scanned
    intp:	pointer to location where an integer can be stored.
  
  Returns: 1 if the number is successfully converted, 0 otherwise.
*/

{ const char *ptr ;
  int val ;

  if (string == NULL) return 0 ;

  val = 0 ;
  for (ptr = string ; *ptr != '\0' ; ptr++)
  { if (*ptr != '0' && *ptr != '1')
      return (0) ;
    val = (val<<1)+*ptr-'0' ; }

  *intp = val ;
  return (1) ; }



bool doterminal (bnfTref_struct *ref)

/*
  This routine handles constructs of type terminal. It checks to be sure the
  reference and definition are well-formed, then scans the required terminal
  from the input stream. The terminal scanned must be the type requested and,
  if an expected string is provided, the terminal must match this string. The
  value is then converted to internal form if necessary and stored as directed
  in the reference.

  Parameter:
    ref:	reference to a terminal definition

  Returns: TRUE if the terminal is successfully parsed, FALSE otherwise.
*/

{ bnfTdef_struct *def ;
  lex_struct *lex = NULL ;
  void **socket ;
  char *lcltxt ;
  int cnt ;
  bool success ;
  const char *rtnnme = "doterminal" ;
  static lex_struct lex_nil = {DY_LCNIL,NULL} ;

/*
  From stdio.h
*/
  extern int sscanf(const char *str, const char *format, ...) ;

/*
  First test the reference to make sure that it references a terminal.
  We'll put off testing for correct storage spec 'til further down when
  we sort out just what will be stored.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfTdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfT)
  { errmsg(41,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref->uflgs,bnfdebug))
    if (debug++ == 0)
      dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    nestlvl-- ; }
# endif

/*
  Call the appropriate scanner to scan the terminal from the input. We check
  here to make sure that we got the class of thing that we expected and abort
  if we didn't.
*/
  success = TRUE ;
  switch (def->ttype)
  { case bnfttNIL:
    { lex = &lex_nil ;
      break ; }
    case bnfttN:
    { lex = dyio_scanlex(bnfchn) ;
      if (lex->class != DY_LCNUM) success = FALSE ;
      break ; }
    case bnfttID:
    { lex = dyio_scanlex(bnfchn) ;
      if (lex->class != DY_LCID) success = FALSE ;
      break ; }
    case bnfttD:
    { lex = dyio_scanlex(bnfchn) ;
      if (lex->class != DY_LCDEL) success = FALSE ;
      break ; }
    case bnfttF:
    { lex = dyio_scanstr(bnfchn,DY_LCFS,def->parm1,'\0','\0') ;
      if (lex->class != DY_LCFS) success = FALSE ;
      break ; }
    case bnfttQ:
    { lex = dyio_scanstr(bnfchn,DY_LCQS,0,def->qschr,def->qechr) ;
      if (lex->class != DY_LCQS && lex->class != DY_LCNIL) success = FALSE ;
      break ; }
    default:
    { errmsg(42,rtnnme,def->ttype) ;
      success = FALSE ;
      break ; } }

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { if (lex->string != NULL)
      dyio_outfmt(dbgchn,dbgecho," = \"%s\"\n",lex->string) ;
    else
      dyio_outfmt(dbgchn,dbgecho," = nil\n") ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

  if (success == FALSE) return (FALSE) ;
  if (def->val != NULL)
  { if (lex->class == DY_LCNIL) return (FALSE) ;
    if (flgon(ref->uflgs,bnfcs) == TRUE)
    { if (flgon(ref->uflgs,bnfmin) == TRUE)
      { if (mstrcmp(lex->string,def->val) != 0) return (FALSE) ; }
      else
      { if (strcmp(lex->string,def->val) != 0) return (FALSE) ; } }
    else
    { if (flgon(ref->uflgs,bnfmin) == TRUE)
      { if (cimstrcmp(lex->string,def->val) != 0) return (FALSE) ; }
      else
      { if (cistrcmp(lex->string,def->val) != 0) return (FALSE) ; } } }
/*
  Deal with the terminal as instructed by the storage instructions. We can
  either store a pointer to the text string or convert the string to some
  internal representation and store it. We don't have to worry about just
  what sort of thing any null object represents, since the whole node is
  zeroed when it is created.
*/
  if (flgon(ref->uflgs,bnfstore) == TRUE)
  { socket = make_socket(curnde,ref->offset) ;
    if (flgon(ref->uflgs,bnfatsgn) == TRUE)
    { if (offset_in_range(ref->offset,sizeof(char *)) == FALSE)
      { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(char *))) ;
	return (FALSE) ; }
      if (lex->class != DY_LCNIL)
      { *((const char **) socket) = STRALLOC(lex->string) ; } }
    else
      switch (def->ttype)
      { case bnfttNIL:
	{ break ; }
	case bnfttN:
	{ if (flgon(ref->uflgs,bnfflt) == TRUE)
	  { if (offset_in_range(ref->offset,sizeof(float)) == FALSE)
	    { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(float))) ;
	      return (FALSE) ; }
	    if (def->parm1 != 10)
	    { errmsg(46,rtnnme,def->parm1,"float") ;
	      return (FALSE) ; }
	    cnt = sscanf(lex->string,"%f",(float *) socket) ; }
	  else
	  if (flgon(ref->uflgs,bnfdbl) == TRUE)
	  { if (offset_in_range(ref->offset,sizeof(double)) == FALSE)
	    { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(double))) ;
	      return (FALSE) ; }
	    if (def->parm1 != 10)
	    { errmsg(46,rtnnme,def->parm1,"double") ;
	      return (FALSE) ; }
	    cnt = sscanf(lex->string,"%lf",(double *) socket) ; }
	  else
	  { if (offset_in_range(ref->offset,sizeof(int)) == FALSE)
	    { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(int))) ;
	      return (FALSE) ; }
	    switch (def->parm1)
	    { case 10:
	      { cnt = sscanf(lex->string,"%d",(int *) socket) ;
		break ; }
	      case 8:
	      { if (*lex->string == '#')
		  cnt = sscanf(lex->string+1,"%o",(unsigned int *) socket) ;
		else
		  cnt = sscanf(lex->string,"%o",(unsigned int *) socket) ;
		break ; }
	      case 2:
	      { cnt = scanbinary(lex->string,(int *) socket) ;
		break ; }
	      default:
	      { errmsg(46,rtnnme,def->parm1,"int") ;
		return (FALSE) ; } } }
	  if (cnt != 1)
	  { errmsg(45,rtnnme,lex->string) ;
	    return (FALSE) ; }
	  break ; }
	case bnfttID:
	{ if (offset_in_range(ref->offset,strlen(lex->string)+1) == FALSE)
	  { errmsg(30,rtnnme,ref->offset,max_offset(strlen(lex->string)+1)) ;
	    return (FALSE) ; }
	  (void) strcpy((char *) socket,lex->string) ;
	  break ; }
	case bnfttD:
	{ if (offset_in_range(ref->offset,sizeof(char)) == FALSE)
	  { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(char))) ;
	    return (FALSE) ; }
	  *((char *) socket) = *lex->string ;
	  break ; }
	case bnfttF:
	{ if (flgon(ref->uflgs,bnfexact) == TRUE)
	  { if (offset_in_range(ref->offset,strlen(lex->string)) == FALSE)
	    { errmsg(30,rtnnme,ref->offset,max_offset(strlen(lex->string))) ;
	      return (FALSE) ; }
	    (void) strncpy((char *) socket,lex->string,def->parm1) ; }
	  else
	  { if (offset_in_range(ref->offset,strlen(lex->string)+1) == FALSE)
	    { errmsg(30,rtnnme,ref->offset,max_offset(strlen(lex->string)+1)) ;
	      return (FALSE) ; }
	    (void) strcpy((char *) socket,lex->string) ; }
	  break ; }
	case bnfttQ:
	{ if (offset_in_range(ref->offset,strlen(lex->string)+1) == FALSE)
	  { errmsg(30,rtnnme,ref->offset,max_offset(strlen(lex->string)+1)) ;
	    return (FALSE) ; }
	  if (lex->class != DY_LCNIL) (void) strcpy((char *) socket,lex->string) ;
	  break ; } } }
/*
  The last thing to do is see if there is an active primitive collecting text,
  and if so concatenate our string onto it.
*/
  if (curtxt != NULL && lex->class != DY_LCNIL)
  { lcltxt = (char *) MALLOC(strlen(curtxt)+strlen(lex->string)+1) ;
    (void) strcpy(lcltxt,curtxt) ;
    (void) strcat(lcltxt,lex->string) ;
    FREE(curtxt) ;
    curtxt = lcltxt ; }
  
  return (TRUE) ; }



bool doimmediate (bnfIref_struct *ref)

/*
  This routine handles constructs of type immediate. It checks to make sure
  that the reference is well-formed, then stores the value specified in the
  definition into the field specified in the reference.

  Parameter:
    ref:	reference to an immediate definition.

  Returns: TRUE if the value is successfully stored, FALSE otherwise.
*/

{ void **socket ;
  bnfIdef_struct *def ;
  const char *rtnnme = "doimmediate" ;

/*
  First test the reference to make sure that it references an immediate and
  that the storage offset is in bounds.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfIdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfI)
  { errmsg(47,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref->uflgs,bnfdebug))
    if (debug++ == 0)
      dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    dyio_outfmt(dbgchn,dbgecho,"\n") ;
    nestlvl-- ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

  if (offset_in_range(ref->offset,sizeof(int)) == FALSE)
  { errmsg(30,rtnnme,ref->offset,max_offset(sizeof(int))) ;
    return (FALSE) ; }
/*
  Do the store.
*/
  socket = make_socket(curnde,ref->offset) ;
  *((int *) socket) = def->ival ;

  return (TRUE) ; }



bool doliteral (bnfLref_struct *ref)

/*
  This routine handles constructs of type literal. It checks to make sure that
  the reference is well-formed, then concatenates the string requested onto
  curtxt.

  Parameter:
    ref:	reference to a literal definition.

  Returns: TRUE if the string is successfully added to curtxt, FALSE otherwise.
*/

{ bnfLdef_struct *def ;
  const char *txt ;
  char *txt2 ;
  const char *rtnnme = "doliteral" ;

/*
  First test the reference to make sure that it references a literal. Also
  check on how the literal is to be obtained. Finally, make sure that we're
  in an active primitive.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfLdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfL)
  { errmsg(48,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref->uflgs,bnfdebug))
    if (debug++ == 0)
      dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    nestlvl-- ; }
# endif

  if (flgon(def->dflgs,bnfsvnm) == TRUE)
    txt = strretrv(addrToInt(def->txt)) ;
  else
    txt = def->txt ;

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { if (txt != NULL)
      dyio_outfmt(dbgchn,dbgecho," = \"%s\"\n",txt) ;
    else
      dyio_outfmt(dbgchn,dbgecho," = nil\n") ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

  if (txt == NULL)
  { errmsg(49,rtnnme) ;
    return (FALSE) ; }
  if (curtxt == NULL)
  { errmsg(50,rtnnme) ;
    return (FALSE) ; }
/*
  It looks OK, so go ahead and get the text and concatenate it to curtxt.
*/
  txt2 = (char *) MALLOC(strlen(curtxt)+strlen(txt)+1) ;
  (void) strcpy(txt2,curtxt) ;
  (void) strcat(txt2,txt) ;
  FREE(curtxt) ;
  curtxt = txt2 ;

  return (TRUE) ; }



bool dolabel (bnfLBref_struct *ref)

/*
  This routine handles socket and label definitions. It first checks to make
  sure defid is a label definition. Next the name field is processed to form
  the string which will be the name of the label, and the node field is
  processed to obtain the value of the label. The offset field from the
  definition is added to the value obtained from the node field to come up
  with the final value of the label. The routine then searches for the label
  in the appropriate defined label list. If it is not defined, it is added to
  the list. If it is defined, the routine checks to see if this definition is
  a redefinition (advance) of the label. If it is, the label is redefined. If
  the current definition is not flagged as a redefinition but the existing
  label is, the label is not redefined. If neither the old or new labels are
  flagged as advances the old label is redefined. The appropriate undefined
  label list is searched if the label is a new definition or the base
  definition of a currently existing label flagged as an advance.

  Parameter:
    ref:	reference to a label definition.
  
  Returns: TRUE if the label definition was processed without error, FALSE
	   otherwise.
*/

{ bnfLBdef_struct *def ;
  const char *nmtxt,*ndtxt ;
  deflbl_struct *lblnde ;
  void *lblval ;
  const char *rtnnme = "dolabel" ;

/*
  First test the reference to make sure that it references a label definition.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfLBdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfDS && def->type != bnfDL)
  { errmsg(52,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon (ref->uflgs,bnfdebug))
    if (debug++ == 0)
      dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    dyio_outchr(dbgchn,dbgecho,'\n') ;
    nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ name ]") ; }
# endif

/*
  Process the name parameters to get the name of the label.
*/
  switch (def->nmcd)
  { case bnfncBNF:
    { if (def->nmsrc == NULL)
      { errmsg(53,rtnnme,"name") ;
	return (FALSE) ; }
      if (doprimitive((bnfPref_struct *) def->nmsrc) == FALSE)
      { errmsg(56,rtnnme,"name") ;
	return (FALSE) ; }
      nmtxt = newtxt ;
      if (flgon(def->dflgs,bnfsvnm) == TRUE) strenter(def->savnm,nmtxt) ;
      break ; }
    case bnfncS:
    { nmtxt = strretrv(addrToInt(def->nmsrc)) ;
      break ; }
    default:
    { errmsg(55,rtnnme,def->nmcd,"name") ;
      return (FALSE) ; } }
  if (nmtxt == NULL)
  { errmsg(57,rtnnme,"name") ;
    return (FALSE) ; }
  nmtxt = STRALLOC(nmtxt) ;

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ \"%s\" ]\n",nmtxt) ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ value ]\n") ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ; }
# endif

/*
  Now process the value parameters to get the label value.
*/
  switch (def->ndcd)
  { case bnfncBNF:
    { if (def->ndsrc == NULL)
      { errmsg(53,rtnnme,"value") ;
	STRFREE(nmtxt) ;
	return (FALSE) ; }
      if (doprimitive((bnfPref_struct *) def->ndsrc) == FALSE)
      { errmsg(56,rtnnme,"value") ;
	STRFREE(nmtxt) ;
	return (FALSE) ; }
      ndtxt = newtxt ;
      lblnde = finddlbl(&blbllst,ndtxt) ;
      if (lblnde == NULL)
      { if (ndtxt == NULL)
	  errmsg(57,rtnnme,"value") ;
	else
	  errmsg(54,rtnnme,ndtxt) ;
	STRFREE(nmtxt) ;
	return (FALSE) ; }
      if (flgon(def->dflgs,bnfsvnd) == TRUE) strenter(def->savnd,ndtxt) ;
      lblval = lblnde->lblval ;
      break ; }
    case bnfncC:
    { lblval = curnde ;
      break ; }
    case bnfncN:
    { lblval = newnde ;
      break ; }
    case bnfncS:
    { ndtxt = strretrv(addrToInt(def->ndsrc)) ;
      lblnde = finddlbl(&blbllst,ndtxt) ;
      if (lblnde == NULL)
      { if (ndtxt == NULL)
	  errmsg(57,rtnnme,"value") ;
	else
	  errmsg(54,rtnnme,ndtxt) ;
	STRFREE(nmtxt) ;
	return (FALSE) ; }
      lblval = lblnde->lblval ;
      break ; }
    default:
    { errmsg(55,rtnnme,def->nmcd,"value") ;
      STRFREE(nmtxt) ;
      return (FALSE) ; } }
/*
  Form the full label from the value plus the offset. Warn the user if it
  evaluates to NULL.
*/
  lblval = make_label(lblval,def->offset) ;

# ifndef DYLP_NDEBUG
  if (lblval == NULL && warnzlbl == TRUE)
    warn(58,rtnnme,nmtxt) ;
  if (debug > 0)
  { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ value: %#1x ]\n",nmtxt) ;
    nestlvl -= 2 ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

/*
  Now search for the label in the appropriate defined label list.
*/
  if (def->type == bnfDS)
    lblnde = finddlbl(&flbllst,nmtxt) ;
  else
    lblnde = finddlbl(&blbllst,nmtxt) ;
/*
  Enter in the appropriate list if the label was previously undefined. Redefine
  the present entry if it is not an advance, or if the label we're now defining
  is flagged to be an advance of the present entry.
*/
  if (lblnde == NULL)
  { lblnde = (deflbl_struct *) MALLOC(sizeof(deflbl_struct)) ;
    if (def->type == bnfDS)
    { lblnde->lblnxt = flbllst ;
      flbllst = lblnde ; }
    else
    { lblnde->lblnxt = blbllst ;
      blbllst = lblnde ; }
    lblnde->lbladv = flgon(def->dflgs,bnfadv) ;
    lblnde->lblnmtxt = nmtxt ;
    lblnde->lblval = lblval ; }
  else
  { if (flgon(def->dflgs,bnfadv) == TRUE || lblnde->lbladv == FALSE)
    { lblnde->lblval = lblval ;
      lblnde->lbladv = flgon(def->dflgs,bnfadv) ; } }
/*
  Now check the appropriate undefined label lists to see if we can resolve any
  entries.
*/
  if (flgoff(def->dflgs,bnfadv) == TRUE)
  { if (def->type == bnfDS)
      lblresolve(&uflbllst,lblnde) ;
    else
      lblresolve(&ublbllst,lblnde) ; }

  return (TRUE) ; }



bool doreference (bnfLBref_struct *ref)

/*
  This routine handles references to sockets and labels. Each is handled in
  basically the same manner, with the exception that for a socket reference,
  the label to be stored in the socket must be defined, and for a label
  reference, the socket which will hold the label must be defined. The
  routine first checks to make sure that defid is a socket/label reference.
  It then processes nmsrc to obtain the socket and ndsrc to obtain the label.
  It completes the reference by making the necessary assignment (if both the
  socket and the label are defined) or by making an entry on the appropriate
  undefined list.

  Parameter:
    ref:	reference to a label reference
  
  Returns: TRUE if the actions described above are completed successfully,
	   FALSE otherwise.
*/

{ bnfLBdef_struct *def ;
  const char *nmtxt = NULL ;
  const char *ndtxt = NULL ;
  deflbl_struct *lblnde ;
  udeflbl_struct *ulblnde ;
  void *val,**socket ;
  bool socket_valid,val_valid ;
  const char *rtnnme = "doreference" ;

/*
  First test the reference to make sure that it references a label reference.
*/
  if (ref == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  def = (bnfLBdef_struct *) ref->defn ;
  if (def == NULL)
  { errmsg(33,rtnnme) ;
    return (FALSE) ; }
  if (def->type != bnfRS && def->type != bnfRL)
  { errmsg(52,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon (ref->uflgs,bnfdebug))
    if (debug++ == 0)
      dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    prtbnfref(dbgchn,dbgecho,(bnfref_struct *) ref) ;
    dyio_outchr(dbgchn,dbgecho,'\n') ;
    nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ socket ]") ; }
# endif

/*
  Now process the nm parameters to get the value of the socket.
*/
  socket_valid = TRUE ;
  switch (def->nmcd)
  { case bnfncBNF:
    { if (def->nmsrc == NULL)
      { errmsg(59,rtnnme,"socket") ;
	return (FALSE) ; }
      if (doprimitive((bnfPref_struct *) def->nmsrc) == FALSE)
      { errmsg(60,rtnnme,"socket") ;
	return (FALSE) ; }
      nmtxt = newtxt ;
      if (nmtxt == NULL)
      { errmsg(61,rtnnme,"socket") ;
	return (FALSE) ; }
      lblnde = finddlbl(&flbllst,nmtxt) ;
      if (lblnde == NULL)
      { socket = NULL ;
	socket_valid = FALSE ; }
      else
	socket = make_socket(lblnde->lblval,0) ;
      if (flgon(def->dflgs,bnfsvnm) == TRUE) strenter(def->savnm,nmtxt) ;
      break ; }
    case bnfncC:
    { socket = make_socket(curnde,def->offset) ;
      break ; }
    case bnfncN:
    { socket = make_socket(newnde,def->offset) ;
      break ; }
    case bnfncS:
    { nmtxt = strretrv(addrToInt(def->nmsrc)) ;
      if (nmtxt == NULL)
      { errmsg(61,rtnnme,"socket") ;
	return (FALSE) ; }
      lblnde = finddlbl(&flbllst,nmtxt) ;
      if (lblnde == NULL)
      { socket = NULL ;
	socket_valid = FALSE ; }
      else
	socket = make_socket(lblnde->lblval,0) ;
      break ; }
    default:
    { errmsg(64,rtnnme,def->nmcd,"socket") ;
      return (FALSE) ; } }
  if (socket_valid == TRUE && socket == NULL)
  { errmsg(62,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ \"%s\" = %#1x ]\n",nmtxt) ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ label ]\n") ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ; }
# endif

/*
  And do the same for the nd parameters to get the label to be stored in the
  socket.
*/
  val_valid = TRUE ;
  switch (def->ndcd)
  { case bnfncBNF:
    { if (def->ndsrc == NULL)
      { errmsg(59,rtnnme,"value") ;
	return (FALSE) ; }
      if (doprimitive((bnfPref_struct *) def->ndsrc) == FALSE)
      { errmsg(60,rtnnme,"value") ;
	return (FALSE) ; }
      ndtxt = newtxt ;
      if (ndtxt == NULL)
      { errmsg(61,rtnnme,"value") ;
	return (FALSE) ; }
      lblnde = finddlbl(&blbllst,ndtxt) ;
      if (lblnde == NULL)
      { val = NULL ;
	val_valid = FALSE ; }
      else
	val = make_label(lblnde->lblval,0) ;
      if (flgon(def->dflgs,bnfsvnd) == TRUE) strenter(def->savnd,ndtxt) ;
      break ; }
    case bnfncC:
    { val = make_label(curnde,def->offset2) ;
      break ; }
    case bnfncN:
    { val = make_label(newnde,def->offset2) ;
      break ; }
    case bnfncS:
    { ndtxt = strretrv(addrToInt(def->ndsrc)) ;
      if (ndtxt == NULL)
      { errmsg(61,rtnnme,"value") ;
	return (FALSE) ; }
      lblnde = finddlbl(&blbllst,ndtxt) ;
      if (lblnde == NULL)
      { val = NULL ;
	val_valid = FALSE ; }
      else
	val = make_label(lblnde->lblval,0) ;
      break ; }
    default:
    { errmsg(64,rtnnme,def->ndcd,"value") ;
      return (FALSE) ; } }

# ifndef DYLP_NDEBUG
  if (warnzlbl == TRUE && val == NULL && val_valid == TRUE)
    warn(65,rtnnme) ;
  if (debug > 0)
  { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ \"%s\" = %#1x ]\n",nmtxt) ;
    nestlvl -= 2 ;
    if (flgon(ref->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

/*
  Complete the reference if both the socket and the label are defined, or make
  an entry in the proper undefined list.
*/
  if (def->type == bnfRS)
  { if (val_valid  == FALSE)
    { errmsg(63,rtnnme) ;
      return (FALSE) ; }
    if (socket_valid == TRUE)
      *socket = val ;
    else
    { ulblnde = (udeflbl_struct *) MALLOC(sizeof(udeflbl_struct)) ;
      ulblnde->lblnxt = uflbllst ;
      uflbllst = ulblnde ;
      ulblnde->lblnmtxt = STRALLOC(nmtxt) ;
      ulblnde->lblval = val ; } }
  else
  { if (socket_valid  == FALSE)
    { errmsg(66,rtnnme) ;
      return (FALSE) ; }
    if (val_valid == TRUE)
      *socket = val ;
    else
    { ulblnde = (udeflbl_struct *) MALLOC(sizeof(udeflbl_struct)) ;
      ulblnde->lblnxt = ublbllst ;
      ublbllst = ulblnde ;
      ulblnde->lblnmtxt = STRALLOC(ndtxt) ;
      ulblnde->lblval = (void *) socket ; } }

  return (TRUE) ; }



bool dolist (bnfref_any ref)

/*
  This routine handles the -LIST construct (applicable to generators,
  non-primitives, and primitives). It first checks to be sure that the
  definition and separator specified by ref are of the correct type. If
  both check out, the routine enters a loop, doing first the body of the
  definition, then the separator. The loop terminates successfully when the
  separator parse fails. It terminates unsuccessfully when the main parse
  fails.

  Parameter:
    ref:	bnf reference flagged as a list construct (bnflst)

  Returns: TRUE if the list parses successfully, FALSE otherwise.
*/

{ bnfref_struct *sepref ;
  bnfdef_any def ;
  void **socket ;
  bool firsttime,success,sepsuccess ;
  void *savcnde = NULL ;
  void *firstnde = NULL ;
  char *savctxt = NULL ;
  char *lclsavtxt = NULL ;
  long marker ;
  const char *rtnnme = "dolist" ;

/*
  Consistency checks. First check that we're dealing with the proper bnf
  reference types for the separator. Errors in body type will be caught below
  as part of the loop setup.
*/
  if (ref.t3 == NULL)
  { errmsg(2,rtnnme,"bnf ref") ;
    return (FALSE) ; }
  sepref = ref.t3->sep ;
  if (sepref->type != bnfP && sepref->type != bnfT)
  { errmsg(67,rtnnme) ;
    return (FALSE) ; }

# ifndef DYLP_NDEBUG
  if (flgon(ref.t3->uflgs,bnfdebug))
    if (debug++ == 0)
    { dyio_outfmt(dbgchn,dbgecho,"\n\n>>>>>> trace begins >>>>>>\n") ;
      nestlvl = 0 ; }
  if (debug > 0)
  { nestlvl++ ;
    printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ begin list ]\n") ; }
# endif

/*
  Set up for the loop which will process the list. Nothing needs to be done for
  non-primitives. For generators, we want to check out the reference storage
  spec and the link storage spec. We put curnde out of the way so we can link
  the list up properly. For primitives, we can check the reference storage
  spec now only if we're storing a pointer to the string. A request for storage
  of the actual string is checked after the end of the parsing loop.  We also
  need to put curtxt out of the way so that we can properly collect the text
  from all primitives in the loop.
*/
  def.com = ref.t3->defn ;
  switch (def.com->type)
  { case bnfG:
    { if (flgon(ref.G->uflgs,bnfstore) == TRUE)
      { if (curnde == NULL)
	{ errmsg(68,rtnnme) ;
	  return (FALSE) ; }
	if (offset_in_range(ref.G->offset,sizeof(void *)) == FALSE)
	{ errmsg(30,rtnnme,ref.G->offset,max_offset(sizeof(void *))) ;
	  return (FALSE) ; } }
      if (def.G->link < 0 ||
	  def.G->link > def.G->size - sizeof(void *))
      { errmsg(69,rtnnme,def.G->link,def.G->size-sizeof(void *)) ;
	return (FALSE) ; }
      savcnde = curnde ;
      break ; }
    case bnfP:
    { if (flgon(ref.P->uflgs,bnfstore) == TRUE)
      { if (curnde == NULL)
	{ errmsg(68,rtnnme) ;
	  return (FALSE) ; }
	if (flgon(ref.P->uflgs,bnfatsgn) == TRUE)
	  if (offset_in_range(ref.P->offset,sizeof(char *)) == FALSE)
	  { errmsg(30,rtnnme,ref.P->offset,max_offset(sizeof(char *))) ;
	    return (FALSE) ; } }
      savctxt = curtxt ;
      curtxt = (char *) MALLOC(sizeof(char)) ;
      *curtxt = '\0' ;
      break ; }
    case bnfNP:
    { break ; }
    default:
    { errmsg(70,rtnnme,def.com->type) ;
      return (FALSE) ; } }
/*
  Now we come to the processing loop. First we try to parse a list element
  according to the body bnf. If this fails, the loop aborts and we return an
  error. If an element is parsed, then an attempt is made to parse the
  separator. If this succeeds, the loop prepares for the next iteration and
  repeats. If the separator parse fails, it is taken to indicate the end of the
  list and we fall out of the loop to the post-processing.
*/
  firsttime = TRUE ;
  while (1)
  {

#   ifndef DYLP_NDEBUG
    if (debug > 0)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      dyio_outfmt(dbgchn,dbgecho,"[ body ]\n") ; }
#   endif

    switch (def.com->type)
    { case bnfG:
      { success = dogenerator(ref.G) ;
	break ; }
      case bnfNP:
      { success = dononprimitive(ref.NP) ;
	break ; }
      case bnfP:
      { success = doprimitive(ref.P) ;
	break ; }
      default:
      { errmsg(1,rtnnme,__LINE__) ;
	return (FALSE) ; } }
    if (success == FALSE) break ;

    if (def.com->type == bnfG)
    { if (firsttime == TRUE)
      { firstnde = newnde ;
	firsttime = FALSE ; }
      else
      { socket = make_socket(curnde,def.G->link) ;
	*socket = newnde ; }
      curnde = newnde ; }

    marker = dyio_mark(bnfchn) ;
    if (flgon(ref.t3->uflgs,bnfstbg) == TRUE && def.com->type == bnfP)
    { lclsavtxt = curtxt ;
      curtxt = NULL ; }

#   ifndef DYLP_NDEBUG
    if (debug > 0)
    { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
      dyio_outfmt(dbgchn,dbgecho,"[ separator ]\n") ; }
#   endif

    if (sepref->type == bnfP)
      sepsuccess = doprimitive((bnfPref_struct *) sepref) ;
    else
      sepsuccess = doterminal((bnfTref_struct *) sepref) ;
    if (sepsuccess == FALSE)
    { if (sepref->type == bnfT) dyio_backup(bnfchn,marker) ;
      if (flgon(ref.t3->uflgs,bnfstbg) == TRUE && def.com->type == bnfP)
	curtxt = lclsavtxt ;
      break ; }

    if (flgon(ref.t3->uflgs,bnfstbg) == TRUE)
    { dyio_backup(bnfchn,marker) ;
      if (def.com->type == bnfP) curtxt = lclsavtxt ; } }
/*
  End of the loop to parse the list. Now we clean up and return. For
  generators, we have to restore curnde and link the head of the list onto it.
  For primitives, we need to look at storing the collected text, then restore
  curtxt and concatenate the collected text to it, if necessary.
*/
  switch (def.com->type)
  { case bnfG:
    { if (success == TRUE)
      { newnde = firstnde ;
	curnde = savcnde ;
	if (flgon(ref.G->uflgs,bnfstore) == TRUE)
	{ socket = make_socket(curnde,ref.G->offset) ;
	  *socket = firstnde ; } }
      else
      { curnde = savcnde ; }
      break ; }
    case bnfP:
    { if (success == TRUE)
      { if (newtxt != NULL) FREE(newtxt) ;
	newtxt = curtxt ;
	curtxt = savctxt ;
	if (flgon(ref.P->uflgs,bnfstore) == TRUE)
	{ if (flgon(ref.P->uflgs,bnfsv) == TRUE)
	  { strenter(ref.P->offset,newtxt) ; }
	  else
	  { socket = make_socket(curnde,ref.P->offset) ;
	    if ( flgon(ref.P->uflgs,bnfatsgn) == TRUE)
	    { *((const char **) socket) = STRALLOC(newtxt) ; }
	    else
	    { if (offset_in_range(ref.P->offset,strlen(newtxt)+1) == FALSE)
	      { errmsg(30,rtnnme,ref.P->offset,max_offset(strlen(newtxt)+1)) ;
		FREE(newtxt) ;
		newtxt = NULL ;
		return (FALSE) ; }
	      (void) strcpy((char *) socket,newtxt) ; } } }
	if (curtxt != NULL)
	{ lclsavtxt = (char *) MALLOC(strlen(curtxt)+strlen(newtxt)+1) ;
	  (void) strcpy(lclsavtxt,curtxt) ;
	  (void) strcat(lclsavtxt,newtxt) ;
	  FREE(curtxt) ;
	  curtxt = lclsavtxt ; } }
      else
      { FREE(curtxt) ;
	curtxt = savctxt ; }
      break ; }
    default:
    { break ; } }

# ifndef DYLP_NDEBUG
  if (debug > 0)
  { printtab(dbgchn,dbgecho,nestlvl,numlvl,tablvl) ;
    dyio_outfmt(dbgchn,dbgecho,"[ end list ]\n") ;
    nestlvl-- ;
    if (flgon(ref.t3->uflgs,bnfdebug))
      if (--debug == 0)
	dyio_outfmt(dbgchn,dbgecho,"<<<<<< trace ends <<<<<<\n\n") ; }
# endif

  return (success) ; }



bool parse (ioid chn, struct bnfref_type3 *bnfid, parse_any *result)

/*
  This routine is used to access the bnf reader to parse input.

  Parameters:
    chn:	Id of the input stream, obtained from openfile
    bnfid:	Bnf to be used to parse the input. The top-level construct
		must be either a generator, non-primitive, or primitive
		(which allows us to get away with bnfref_type3 as the common
		type, otherwise we'd be forced to void).
    result:	Will be assigned the data structure or character string
		built by the parse.

  Returns: TRUE if the parse succeeds, FALSE otherwise. Note that the value
	   built by the parse is returned in result.
*/

{ bool success ;
  bnfref_any ref ;
  const char *rtnnme = "parse" ;

/*
  Make sure we have a valid bnf reference, and some place to put the result for
  generators and primitives.
*/
  if (bnfid == NULL)
  { errmsg(2,rtnnme,"bnf") ;
    return (FALSE) ; }
  ref.t3 = bnfid ;
  if (ref.com->type != bnfG &&
      ref.com->type != bnfNP &&
      ref.com->type != bnfP)
  { errmsg(43,rtnnme) ;
    return (FALSE) ; }
  if (ref.com->type != bnfNP)
    if (result == NULL)
    { errmsg(2,rtnnme,"result") ;
      return (FALSE) ; }
/*
  Now set the input channel, call the appropriate routine to do the parse, and
  return the result.
*/
  bnfchn = chn ;
  switch (ref.com->type)
  { case bnfG:
    { if (flgon(ref.G->uflgs,bnflst) == TRUE)
	success = dolist(ref) ;
      else
	success = dogenerator(ref.G) ;
      if (success == TRUE) result->g = newnde ;
      break ; }
    case bnfNP:
    { if (flgon(ref.NP->uflgs,bnflst) == TRUE)
	success = dolist(ref) ;
      else
	success = dononprimitive(ref.NP) ;
      break ; }
    case bnfP:
    { if (flgon(ref.P->uflgs,bnflst) == TRUE)
	success = dolist(ref) ;
      else
	success = doprimitive(ref.P) ;
      if (success == TRUE) result->c = newtxt ;
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      success = FALSE ;
      break ; } }
  
  return (success) ; }

# ifndef DYLP_NDEBUG


void bnfdbgctl (ioid dbgchn_p, bool dbgecho_p,
		bool warnzlbl_p, bool numlvl_p, bool tablvl_p)

/*
  By default the debugging trace warns about labels that evaluate to NULL,
  and (once triggered by the bnfdebug flag in a reference) prints trace
  lines prefixed with the nesting level and proportionally indented. The output
  will go to stdout only. This routine allows the behaviour to be changed.

  Parameters:
    dbgchn:	i/o id for trace output (default IOID_NOSTRM)
    dbgecho:	TRUE to echo trace output to stdout (default TRUE)
    warnzero:	TRUE to produce warning messages for null labels (default TRUE)
    numlvl:	TRUE to prefix trace lines with the nesting level (default
    		TRUE)
    tablvl:	TRUE to indent trace lines proportional to the nesting level
		(default TRUE)
*/

{ dbgchn = dbgchn_p ;
  dbgecho = dbgecho_p ;
  warnzlbl = warnzlbl_p ;
  numlvl = numlvl_p ;
  tablvl = tablvl_p ; }

# endif /* DYLP_NDEBUG */

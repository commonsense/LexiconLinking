#ifndef _DYLIB_BNFRDR_H
#define _DYLIB_BNFRDR_H

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

#include "dylib_io.h"

/*
  sccs: @(#)bnfrdr.h	3.5	09/01/99
  svn/cvs: "$Id: dylib_bnfrdr.h 71 2006-06-09 04:21:15Z andreasw $" ;

  This file contains the structure definitions required to use the bnf reader
  package. The bnf reader depends heavily on the following two assumptions:

  * There is some pointer type into which any other pointer can be cast and
    recovered.
  
  * An int can be cast into a pointer and recovered. This is used to prevent
    the complexity of the bnf data structure from getting out of hand, but
    could be avoided at the expense of a substantial increase in its size and
    awkwardness of use.

  The basic scheme is something like this. At the bottom, we have a number of
  terminal constructs: immediates, literals, terminals, and labels of various
  flavours. Above these are the three non-terminal constructs: primitives,
  non-primitives, and generators. The non-terminals have bodies which are made
  up of references to other terminal or non-terminal constructs. Generator
  bodies have only one parse; primitive and non-primitive bodies can have a
  number of alternative parses. A picture is probably in order here:

  definition      body

  ----------      ----------	  ----------      ----------
  |        | ---> |        | ---> |  ref   | ---> |  defn  |
  ----------      ----------      ----------      ----------
		  |  ....  |
		  ----------	  ----------      ----------
		  |        | ---> |  ref   | ---> |  defn  |
		  ----------	  ----------      ----------

  A definition contains a pointer to its body, which is an array of pointers to
  references. Each reference points to a definition. Essentially, a reference
  specifies how the second definition is to be used in the context of the body
  of the first definition.

  The bnf reader has the capability to create arbitrary links in the data
  structure it's building. Some terminology will make things easier:

  		|  ....  |
  		----------      ----------
      socket -->| label -+----->|        |
  		----------      ----------
  		|  ....  |      |  ....  |

  The value of a label is the address of something. To make a link, the
  value of the label has to be stored.  The place where it is stored is
  called a socket. The value of a socket is the address of the field where
  the value of the label is stored. When it defines a socket, the bnf reader
  associates a name with an address; similarly for a label. Both socket and
  label references cause a label to be stored in a socket; the difference
  between the two lies in which of the socket or label can be undefined when
  the reference is processed.

  When it's not important to distinguish between sockets and labels, the
  documentation uses label to include both.
  
  To write a bnf, you use the set of macros defined at the end of
  the file. For a detailed explanation of the sorts of things that can be
  specified with the bnf, the user should take a look at the supplementary
  documentation. The structures and code will make a lot more sense afterward.
*/



/*
  Definitions of enum types used as codes in the bnf structures which follow.
*/

/*
  bnftype_enum codes the type of the bnf definition.

  Value		Description
  -----		-----------
  bnfG		Generator definition
  bnfNP		Non-primitive definition
  bnfP		Primitive definition
  bnfT		Terminal definition
  bnfDS		Socket definition definition
  bnfDL		Label definition definition
  bnfRS		Socket reference definition 
  bnfRL		Label reference definition
  bnfI		Immediate value definition
  bnfL		Literal definition
*/

typedef enum {bnfG,bnfNP,bnfP,bnfT,bnfDS,
	      bnfDL,bnfRS,bnfRL,bnfI,bnfL} bnftype_enum ;


/*
  bnfttype_enum codes the type of lexeme expected by a terminal.

  Value		Description
  -----		-----------
  bnfttNIL	the null lexeme
  bnfttN	number
  bnfttID	identifier
  bnfttD	delimiter
  bnfttF	fixed-length string
  bnfttQ	quoted string
*/

typedef enum {bnfttNIL,bnfttN,bnfttID,bnfttD,bnfttF,bnfttQ} bnfttype_enum ;


/*
  bnflblsrc_enum codes the way in which text strings used for label names
  are obtained.

  Value		Description
  -----		-----------
  bnfncBNF      A bnf is supplied which will produce a text string. If this
		code appears in the context of a name, the string will be the
		name of the label. If it appears in the context of a value,
		the string will be used as a label name and the value
		associated with the name will become the value of the label
		being defined.
  bnfncS	An index in the saved text array is supplied. The string
		retrieved is interpreted as for bnfncBNF.
  bnfncC	The value of curnde is used as the socket/label value. This
		code is not valid in the context of a name.
  bnfncN	The value of newnde is used as the socket/label value. This
		code is not valid in the context of a name.
*/

typedef enum {bnfncBNF,bnfncS,bnfncC,bnfncN} bnflblsrc_enum ;



/*
  Flag definitions used in bnf definitions.

  Flag		Description
  ----		-----------
  bnfadv	Indicates the redefinition of a previously defined label. The
		usual context for use is to redefine (advance) a label which
		is the link in a linked list.
  bnfsvnd	Save the text string developed by the nd part of a label
		definition definition or label reference definition.
  bnfsvnm	Save the text string developed by the nm part of a label
		definition definition or label reference definition. This
		flag is also used in literal definitions to indicate that text
		should be retrieved from the saved text array.
*/

#define bnfadv	1<<0
#define bnfsvnd	1<<1
#define bnfsvnm	1<<2


/*
  Flag definitions used in bnf references.

  Flag		Description
  ----		-----------
  bnflst	The definition referenced describes one element of a list of
		indefinite length.
  bnfstore	The value produced by the referenced bnf will be stored
		somehow, and the offset field should be valid.
  bnfatsgn	Store a pointer to the character string produced by the
		referenced bnf, rather than the string itself.
  bnfstbg	The bnf referenced as the separator between list elements is
		really the beginning of the next list element. (Hence we'll
		have to back up over it once we recognize it.)
  bnfflt	A float number is expected here.
  bnfdbl	A double number is expected here.
  bnfcs		Forces a case-sensitive comparison of the string read for a
		terminal with the value specified in the terminal definition.
  bnfmin	Requests a minimum-length comparison - as long as the string
		parsed for the terminal matches the value specified in the
		terminal definition up to the end of the parsed string, the
		comparison succeeds.
  bnfsv		Used in primitives to indicate that the string is to be stored
		in the savedtxt array. The offset should be a valid savedtxt
		index in this case.
  bnfexact	(bnfttF only) used to prevent the addition of the null
		terminator at the end of a character string when the string
		is stored directly in a field (must be specified to store a
		single char in a field of size sizeof(char))
  bnfdebug	Debugging should be activated for this reference and all
		bnf rules nested within it.
*/

#define bnflst		1<<0
#define bnfstore	1<<1
#define bnfatsgn	1<<2
#define bnfstbg		1<<3
#define bnfflt		1<<4
#define bnfcs		1<<5
#define bnfmin		1<<6
#define bnfsv		1<<7
#define bnfexact	1<<8
#define bnfdebug	1<<9
#define bnfdbl		1<<10



/*
  Data structures used for bnf definitions. There are three types of things
  here: individual structures for the various definition types, a common
  structure which consists only of the fields common to all of the individual
  structures, and a pointer union which is handy when walking around in a bnf.

  Just to keep the explanation in hand a bit, let's define components and
  alternatives. The body of a bnf definition consists of alternatives
  (alternative parses), each of which has a number of components. A component
  array is an array of pointers to bnf reference structures, each of which in
  turn references a bnf definition structure. An alternative array is an array
  of pointers to component arrays. I know this is ugly and involves a lot of
  dereferencing but it seems to be the only way to handle the variable lengths
  involved. 

  NOTE: To keep things from getting completely out of hand, the first entry in
	a component or alternative array specifies the number of pointers that
	follow. This is one of the (ab)uses of the int - pointer - int cast.
*/

/*
  The common portion.

  Field		Description
  -----		-----------
  type		Type code identifying what sort of definition this is.
  name		The name of the rule (derived from the C variable name;
		see the macros gdef, npdef, etc.)
*/

#define bnfdef_common bnftype_enum type ; \
		      const char *name ;

typedef struct { bnfdef_common } bnfdef_struct ;


/*
  Data structure for a generator definition. Generators cause the creation of a
  node in the data structure being built for the user. For simplicity, they may
  not have alternative parses, but since they can reference non-primitives no
  flexibility is lost.

  Field		Description
  -----		-----------
  bnfdef_common	As above.
  size		Size (in bytes) of the node to be created.
  link		Offset (in bytes) from the base of the node created by the
		generator to the field used as a link field when this node is
		in a linked list.
  comps		Pointer to a component array.
*/

typedef struct { bnfdef_common
		 int size ;
		 int link ;
		 struct bnfref_struct_tag **comps ; } bnfGdef_struct ;


/*
  Data structure for a non-primitive definition. Non-primitives are simply a
  device for defining alternative parses. They don't directly create anything.

  Field		Description
  -----		-----------
  bnfdef_common	As above.
  alts		Pointer to an alternative array.
*/

typedef struct {bnfdef_common
		struct bnfref_struct_tag ***alts ; } bnfNPdef_struct ;
		


/*
  Data structure for a primitive definition. The distinction between a
  primitive and a non-primitive is that a primitive constructs a string which
  is the concatenation of the strings returned by the bnf's referenced in the
  primitive's body. The data structure is identical to that for non-primitives.
*/

typedef bnfNPdef_struct bnfPdef_struct ;


/*
  Data structure for a terminal. Terminals are used to specify specific things
  to be obtained from the input stream. The various parameters required to
  describe a terminal should really be mushed into a union, but then the bnf
  data structure would have to be built dynamically, since unions can't be
  initialized.

  Field		Description
  -----		-----------
  bnfdef_common	As above.
  ttype		Code identifying the type of terminal to be obtained.
  qschr		Starting character for a quoted string.
  qechr		Ending character for a quoted string.
  parm1		Overloaded field, interpreted as follows:
		numbers: specifies the radix
		fixed-length strings: specifies the string length
  val		Expected value of the string obtained from the input stream.
		(This test is applied before the string is converted to the
		 internal form appropriate for whatever is specified in ttype.)
*/

typedef struct { bnfdef_common
		 bnfttype_enum ttype ;
		 char qschr ;
		 char qechr ;
		 int parm1 ;
		 const char *val ; } bnfTdef_struct ;


/*
  Data structure for an immediate value. Immediates are used to jam a code into
  the data structure being built.

  Field		Description
  -----		-----------
  bnfdef_common	As above.
  ival		Integer value.
*/

typedef struct {bnfdef_common
		int ival ; } bnfIdef_struct ;


/*
  Data structure for a literal. Literals are used to insert characters into the
  input stream. (Handy for generating label names, for instance.)

  Field		Description
  -----		-----------
  bnfdef_common	As above.
  dflgs		Flags.
  txt		The string to be inserted. This field is also used to index
		into the saved text array by casting it to an int.
*/

typedef struct { bnfdef_common
		 flags dflgs ;
		 char *txt ; } bnfLdef_struct ;


/*
  Last but not least, the data structure used to define socket/label
  definitions and references. (Definitions, mind you - there is another
  structure to reference socket/label definitions and references.) A
  socket/label definition associates of a name (a text string) with a value
  (almost always an address).  A socket/label reference specifies a socket
  and a label. The label is inserted into the socket. Fields prefixed by nm
  are the name in a socket/label definition and the socket in a socket/label
  reference. Fields prefixed by nd are the value in a socket/label definition
  and the label in a socket/label reference.

  Field		Description
  -----		-----------
  bnfdef_common	As above.
  dflgs		Flags.
  nmcd		Specifies how name/socket will be obtained.
  ndcd		Specifies how value/label will be obtained.
  savnm		Specifies location in saved text array where string associated
		with nm will be stored.
  nmsrc		Pointer to bnf which will produce string for nm, or cast into
		an int and used as a location in the saved text array.
  savnd		Specifies location in saved text array where string associated
		with nd will be stored.
  ndsrc		Pointer to bnf which will produce string for nd, or cast into
		an int and used as a location in the saved text array.
  offset	Correction (in bytes) to socket/label value (socket/label
		definitions) or socket (socket/label references).
  offset2	Correction (in bytes) to label (socket/label references).
*/

typedef struct { bnfdef_common
		 flags dflgs ;
		 bnflblsrc_enum nmcd ;
		 bnflblsrc_enum ndcd ;
		 int savnm ;
		 struct bnfref_struct_tag *nmsrc ;
		 int savnd ;
		 struct bnfref_struct_tag *ndsrc ;
		 int offset ;
		 int offset2 ; } bnfLBdef_struct ;


/*
  And finally, the handy union of pointers promised back at the start. We
  really should be using this in the bnf reference structure declarations,
  rather than (bnfdef_struct *), but since references and definitions are
  mutually recursive we get into ugliness. There's also the point that we
  want to be able to create bnfs at compile time and you can't initialize
  unions.
*/

typedef union { bnfdef_struct *com ;
		bnfGdef_struct *G ;
		bnfNPdef_struct *NP ;
		bnfPdef_struct *P ;
		bnfTdef_struct *T ;
		bnfIdef_struct *I ;
		bnfLdef_struct *L ;
		bnfLBdef_struct *LB ; } bnfdef_any ;



/*
  Now, on to the data structures used to reference bnf definitions. Recall if
  you will the introductory comments about component and alternative arrays and
  the general setup of the bnf data structure. We have the same three types of
  data structures here as for bnf definitions.
*/

/*
  The common portion. It includes a type code, a name, usage flags, and a
  pointer to the bnf definition.

  Field		Description
  -----		-----------
  type		Type code identifying what sort of definition this reference
		points to.
  name		The name of the reference (derived from the C variable name;
		see the macros qref, npref, pref, etc.)
  uflgs		Usage flags.
  defn		Pointer to a bnf definition structure.
*/

#define bnfref_common bnftype_enum type ; \
		      const char *name ; \
		      bnfdef_struct *defn ; \
		      flags uflgs ;

typedef struct bnfref_struct_tag { bnfref_common } bnfref_struct ;


/*
  References to labels of all flavours and to literals require only the
  common fields. The only reason we need the uflgs field is for the bnfdebug
  flag.
*/

typedef bnfref_struct bnfLBref_struct ;
typedef bnfref_struct bnfLref_struct ;


/*
  References to terminals and immediates require an offset for storage.

  Field		Description
  -----		-----------
  bnfref_common	As above.
  offset	Offset (in bytes) into current node to the field where the
		value produced by the referenced bnf will be stored.
*/

struct bnfref_type2 { bnfref_common
		      int offset ; } ;

typedef struct bnfref_type2 bnfTref_struct ;
typedef struct bnfref_type2 bnfIref_struct ;


/*
  References to generators, non-primitives, and primitives can be in lists and
  require a separator specification in addition to the offset. Non-primitives
  do not make use of the offset field.

  Field		Description
  -----		-----------
  bnfref_common	As above.
  offset	Offset (in bytes) into current node to the field where the
		value produced by the referenced bnf will be stored.
  sep		A reference to a bnf definition describing the separator
		between list elements in the input stream.
*/

struct bnfref_type3 { bnfref_common
		      int offset ;
		      bnfref_struct *sep ; } ;

typedef struct bnfref_type3 bnfGref_struct ;
typedef struct bnfref_type3 bnfNPref_struct ;
typedef struct bnfref_type3 bnfPref_struct ;


/*
  And the handy union pointer type. Same general comments as for the
  declaration of bnfdef_any.
*/

typedef union { bnfref_struct *com ;
		struct bnfref_type1 *t1 ;
		struct bnfref_type2 *t2 ;
		struct bnfref_type3 *t3 ;
	        bnfGref_struct *G ;
		bnfNPref_struct *NP ;
		bnfPref_struct *P ;
		bnfTref_struct *T ;
		bnfIref_struct *I ;
		bnfLref_struct *L ;
		bnfLBref_struct *LB ; } bnfref_any ;



/*
  The macros that make defining the bnf data structures marginally
  less painful.
*/

/*
  Macros to help with constructing field offsets. NULLP is specially designed
  to produce a NULL value when used as &NULLP. This is required for some of the
  macros where one must fill the field with either the address of a
  bnfref_struct or the value NULL. By this device we avoid having to make the
  user aware of when and when not to use &. mkoff simply produces the offset of
  a given field in a structure type.
*/

#define NULLP (*((char *) 0))
#define mksav(qqoff) (*((char *) qqoff))
#define mkoff(qqtype,qqfield) (&((qqtype *) 0)->qqfield)

/*
  Macros for alternative and component lists. These just generate the headers;
  the actual lists have to be typed out, as:

  althd(arule_alts) = { altcnt(3),
			mkaref(arule_alt1), mkaref(arule_alt2),
			mkaref(arule_alt3) } ;

  comphd(arule_alt1) = { compcnt(2),
			 mkcref(brule_ref), mkcref(crule_ref) } ;

  where brule_ref and crule_ref are bnf references (most likely constructed
  using the gref, npref, etc. macros).
*/

#define althd(qqnme) bnfref_struct **qqnme[]
#define altcnt(qqcnt) (bnfref_struct **) (qqcnt)
#define mkaref(qqref) (bnfref_struct **) (qqref)

#define comphd(qqnme) bnfref_struct *qqnme[]
#define compcnt(qqcnt) (bnfref_struct *) (qqcnt)
#define mkcref(qqref) (bnfref_struct *) (&qqref)

/*
  Macros to initialise bnf definitions. Note the use of the ANSI C
  'stringisation' operator, '#', to get a text string for the name. For
  non-ANSI implementations, replacing #qqnme with "qqnme" usually works (but
  not all non-ANSI preprocessor implementations will see the macro parameter
  inside a string, and ANSI C explicitly disallows it).
*/

#define gdef(qqnme,qqsze,qqlnk,qqcomps) \
bnfGdef_struct qqnme = { bnfG, #qqnme, (int) (qqsze), (int) (qqlnk), \
			 (bnfref_struct **) qqcomps }

#define npdef(qqnme,qqalts) \
bnfNPdef_struct qqnme = { bnfNP, #qqnme, (bnfref_struct ***) qqalts }

#define pdef(qqnme,qqalts) \
bnfPdef_struct qqnme = { bnfP, #qqnme, (bnfref_struct ***) qqalts }

#define tdef(qqnme,qqttype,qqparm,qqval) \
bnfTdef_struct qqnme = { bnfT, #qqnme, qqttype, '\0', '\0', \
			 (int) (qqparm), (const char *) (qqval) }

#define tqdef(qqnme,qqschr,qqechr,qqval) \
bnfTdef_struct qqnme = { bnfT, #qqnme, bnfttQ, (char) qqschr, (char) qqechr,\
			 0, (char *) (qqval) }

#define dfdef(qqnme,qqdflgs,qqnmcd,qqnm,qqsavnm,qqndcd,qqnd,qqsavnd,qqoff) \
bnfLBdef_struct qqnme = { bnfDS, #qqnme, (flags) (qqdflgs), qqnmcd, qqndcd, \
			  (int) (qqsavnm), (bnfref_struct *) &qqnm, \
			  (int) (qqsavnd), (bnfref_struct *) &qqnd, \
			  (int) (qqoff), 0 }

#define dbdef(qqnme,qqdflgs,qqnmcd,qqnm,qqsavnm,qqndcd,qqnd,qqsavnd,qqoff) \
bnfLBdef_struct qqnme = { bnfDL, #qqnme, (flags) (qqdflgs), qqnmcd, qqndcd, \
			  (int) (qqsavnm), (bnfref_struct *) &qqnm, \
			  (int) (qqsavnd), (bnfref_struct *) &qqnd, \
			  (int) (qqoff), 0 }

#define rfdef(qqnme,qqdflgs,qqnmcd,qqnm,qqsavnm,qqoff,qqndcd,qqnd,qqsavnd,qqoff2) \
bnfLBdef_struct qqnme = { bnfRS, #qqnme, (flags) (qqdflgs), qqnmcd, qqndcd, \
			  (int) (qqsavnm), (bnfref_struct *) &qqnm, \
			  (int) (qqsavnd), (bnfref_struct *) &qqnd, \
			  (int) (qqoff), (int) (qqoff2) }

#define rbdef(qqnme,qqdflgs,qqnmcd,qqnm,qqsavnm,qqoff,qqndcd,qqnd,qqsavnd,qqoff2) \
bnfLBdef_struct qqnme = { bnfRL, #qqnme, (flags) (qqdflgs), qqnmcd, qqndcd, \
			  (int) (qqsavnm), (bnfref_struct *) &qqnm, \
			  (int) (qqsavnd), (bnfref_struct *) &qqnd, \
			  (int) (qqoff), (int) (qqoff2) }

#define idef(qqnme,qqval) \
bnfIdef_struct qqnme = { bnfI, #qqnme, (int) (qqval) }

#define ldef(qqnme,qqdflgs,qqtxt) \
bnfLdef_struct qqnme = { bnfL, #qqnme, (flags) (qqdflgs), (char *) (qqtxt) }



#define gref(qqnme,qqref,qquflgs,qqoff,qqsep) \
bnfGref_struct qqnme = { bnfG, #qqnme, (bnfdef_struct *) &qqref, \
			 (flags) (qquflgs), (int) (qqoff), \
			 (bnfref_struct *) &qqsep }

#define npref(qqnme,qqref,qquflgs,qqsep) \
bnfNPref_struct qqnme = { bnfNP, #qqnme, (bnfdef_struct *) &qqref, \
			 (flags) (qquflgs), (int) 0, (bnfref_struct *) &qqsep }

#define pref(qqnme,qqref,qquflgs,qqoff,qqsep) \
bnfPref_struct qqnme = { bnfP, #qqnme, (bnfdef_struct *) &qqref, \
			 (flags) (qquflgs), (int) (qqoff), \
			 (bnfref_struct *) &qqsep }

#define tref(qqnme,qqref,qquflgs,qqoff) \
bnfTref_struct qqnme = { bnfT, #qqnme, (bnfdef_struct *) &qqref, \
			 (flags) qquflgs, (int) qqoff }

#define dfref(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfDS, #qqnme, (bnfdef_struct *) &qqref, (flags) 0 }

#define dbref(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfDL, #qqnme, (bnfdef_struct *) &qqref, (flags) 0 }

#define rfref(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfRS, #qqnme, (bnfdef_struct *) &qqref, (flags) 0 }

#define rbref(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfRL, #qqnme, (bnfdef_struct *) &qqref, (flags) 0 }

#define iref(qqnme,qqref,qqoff) \
bnfIref_struct qqnme = { bnfI, #qqnme, (bnfdef_struct *) &qqref, \
			(flags) 0, (int) qqoff }

#define lref(qqnme,qqref) \
bnfLref_struct qqnme = { bnfL, #qqnme, (bnfdef_struct *) &qqref, (flags) 0 }

#ifndef DYLP_NDEBUG

/*
  This set of definitions sets the bnfdebug flag, but doesn't add a separate
  uflgs parameter (we don't want to lead the user to think any of the others
  are valid).
*/

#define dfrefdbg(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfDS, #qqnme, (bnfdef_struct *) &qqref, \
			  (flags) bnfdebug }

#define dbrefdbg(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfDL, #qqnme, (bnfdef_struct *) &qqref, \
			  (flags) bnfdebug }

#define rfrefdbg(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfRS, #qqnme, (bnfdef_struct *) &qqref, \
			  (flags) bnfdebug }

#define rbrefdbg(qqnme,qqref) \
bnfLBref_struct qqnme = { bnfRL, #qqnme, (bnfdef_struct *) &qqref, \
			  (flags) bnfdebug }

#define lrefdbg(qqnme,qqref) \
bnfLref_struct qqnme = { bnfL, #qqnme, (bnfdef_struct *) &qqref, \
			  (flags) bnfdebug }

#endif /* DYLP_NDEBUG */



/*
  Last, but not least, some declarations to allow the use of the bnf reader.
  rdrinit and rdrclear initialize and clear the reader; they should bracket
  related groups of calls. parse is the main parsing routine. The union type
  parse_any is the appropriate thing to hold the result.
*/

typedef union { void *g ;
		char *c ; } parse_any ;

extern void rdrinit(void),rdrclear(void) ;
extern bool parse(ioid chn, struct bnfref_type3 *bnfid, parse_any *result) ;

#ifndef DYLP_NDEBUG
/*
  The control routine for the bnf debugging trace output. See the comments
  in bnfrdr.c for the proper use of the parameters.
*/
  
extern void bnfdbgctl(ioid dbgchn, bool dbgecho, bool warnzlbl, bool numlvl,
		      bool tablvl) ;
#else
#define bnfdbgctl(dgbchn,dbgecho,warnzlbl,numlvl,tablvl)
#endif

/*
  Utility print routines from bnfrdrio.c.
*/

extern void prtbnfref(ioid chn, bool echo, bnfref_struct *ref),
	    prtbnfdef(ioid chn, bool echo, bnfdef_struct *def) ;

#endif /* _DYLIB_BNFRDR_H */

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
  The routines in this file implement a package for storage and allocation of
  text strings. The package stores all strings in a hash table, merging
  identical strings. A count is maintained of the number of outstanding
  references to the string. There are two access routines, stralloc for
  inserting strings and strfree for releasing them.

  I/O is handled entirely by stdio facilities unless the malloc debugging
  macros (loustd.h) are compiled in, in which case the io library is needed.
*/

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)littab.c	1.4	09/25/04" ;
static char svnid[] UNUSED = "$Id: dylib_littab.c 71 2006-06-09 04:21:15Z andreasw $" ;

#include <stdio.h>
#include "dylib_hash.h"


/*
  Definition of litent data structure.

  Field		Definition
  -----		----------
  refs		Reference count for this string
  text		Pointer to the string
*/

typedef struct litent_internal { int refs ;
				 char *text ; } litent ;

/*
  ANSI C specifies that littable will be initialised to NULL; this is
  important.
*/

#define LITTABLESIZE 2039

static hel *littable[LITTABLESIZE]  ;



const char *stralloc (const char *string)

/*
  This routine is called by the user to allocate a permanent copy of string.
  In reality, it will search the hash table for an entry corresponding to
  string. If it finds one, it will increment the reference count and return a
  pointer to the string stored in the entry. If there is no existing entry,
  one is created, a copy of the string is made, and a pointer to this copy is
  returned.

  Parameters:
    string:	text string

  Returns: pointer to a usable copy of string, or NULL in the event of an
	   error.
*/

{ litent *lit ;
  const char *rtnnme = "stralloc" ;

  if (string == NULL)
  { fprintf(stderr,"\n%s: null string parameter!\n",rtnnme) ;
    return (NULL) ; }

  lit = (litent *) lookup(string,littable,LITTABLESIZE) ;

  if (lit != NULL)
  { lit->refs++ ;
    return (lit->text) ; } 

  lit = (litent *) MALLOC(sizeof(litent)) ;
  lit->text = MALLOC(strlen(string)+1) ;
  strcpy(lit->text,string) ;
  lit->refs = 1 ;
  if (enter(lit->text,littable,LITTABLESIZE,(char *) lit) == NULL)
  { fprintf(stderr,"\n%s: couldn't enter string \"%s\" in literal table!\n",
	    rtnnme,string) ;
    FREE(lit->text) ;
    FREE(lit) ;
    return (NULL) ; }

  return (lit->text) ; }



bool strfree (const char *string)

/*
  This routine is called by the user to "free" a string. In reality, it will
  search the hash table for an entry corresponding to string. If it finds one,
  it will decrement the reference count. It is an error if there is no entry
  for string, which implies that at some point a pointer to the string was
  acquired without consulting this package.

  Parameter:
    string:	text string

  Returns: TRUE if the string was successfully "freed", FALSE otherwise.
*/

{ litent *lit ;
  const char *rtnnme = "strfree" ;

  if (string == NULL)
  { fprintf(stderr,"\n%s: null string parameter!\n",rtnnme) ;
    return (FALSE) ; }

  lit = (litent *) lookup(string,littable,LITTABLESIZE) ;

  if (lit == NULL)
  { fprintf(stderr,"\n%s: no entry for string \"%s\" in literal table!\n",
	    rtnnme,string) ;
    return (FALSE) ; }
  
  if (--lit->refs == 0)
  { if (erase(lit->text,littable,LITTABLESIZE) == NULL)
    { fprintf(stderr,"\n%s: confusion deleting entry for string \"%s\"!\n",
	      rtnnme,lit->text) ;
      return (FALSE) ; }
    FREE(lit->text) ;
    FREE(lit) ; }
  
  return (TRUE) ; }

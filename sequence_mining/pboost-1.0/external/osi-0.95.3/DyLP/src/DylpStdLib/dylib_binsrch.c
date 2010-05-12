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

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)binsrch.c	1.4	09/25/04" ;
static char svnid[] UNUSED = "$Id: dylib_binsrch.c 71 2006-06-09 04:21:15Z andreasw $" ;

#include "dylib_strrtns.h"
#include "dylib_keytab.h"



int find (char *word, keytab_entry keytab[], int numkeys)     

/*
  Find looks up a word in a keyword table. It uses binary search and is
  insensitive to case.

  Parameters:
    word:	pointer to the word to be looked up
    keytab:	pointer to the array of keywords
    numkeys:	length of the keyword table

  Return value: 
    Normal     token number associated with the keyword
    Otherwise  -1 if the keyword cannot be found
	       
*/

{ int high,low,mid,place ;

  low = 0 ;
  high = numkeys-1 ;
  while (low <= high) 
  { mid = (high+low)/2 ;
    place = cistrcmp(word,keytab[mid].keyword) ;
    switch (place)
    { case -1:
        high = mid-1 ;
        break ;
      case 0:
	return (keytab[mid].token) ;
      case 1:
	low = mid+1 ;
	break ; } }
  return (-1) ; }



int ambig (char *word, keytab_entry keytab[], int numkeys)

/*
  This routine looks up a word in a keyword table using a "shortest unique
  match" algorithm. If the minimum number of characters to be matched is
  greater than the length of the keyword, it's taken to be a request for an
  exact match (think of it as requiring the terminating '\0' to be matched).

  Parameters:
    word:	word to be looked up
    keytab:	array of keywords
    numkeys:	length of the keyword array

  Returns:
    Normal:	token number associated with the keyword
    Otherwise:	-1 if the keyword cannot be found
		-2 if the keyword is ambiguous
*/

{ int high,low,mid,place ;

  low = 0 ;
  high = numkeys-1 ;
  while (low <= high)
  { mid = (high+low)/2 ;
    place = cimstrcmp(word,keytab[mid].keyword) ;
    switch(place)
    { case -1:
      { high = mid-1 ;
	break ; }
      case 0:
      { if (strlen(word) < keytab[mid].min-1)
	  return (-2) ;
	else
	if (strlen(keytab[mid].keyword) >= keytab[mid].min)
	  return (keytab[mid].token) ;
	else
	if (strlen(keytab[mid].keyword) == strlen(word))
	  return (keytab[mid].token) ; }
      case 1:
      { low = mid+1 ;
	break ; } } }
  return (-1) ; }

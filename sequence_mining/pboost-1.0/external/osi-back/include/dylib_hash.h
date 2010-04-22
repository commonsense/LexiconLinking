#ifndef _DYLIB_HASH_H
#define _DYLIB_HASH_H
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

/*                Hash Table Structures                 */

/*
  In order that the hashing routines can be used for a number of different
  tables, they do not have any knowledge of the thing being hashed. All that
  is maintained is the association between a key and some generic object.

  @(#)hash.h	1.3 06/22/04
  svn/cvs: $Id: dylib_hash.h 71 2006-06-09 04:21:15Z andreasw $
*/

/* 
  The basic hash table entry structure
  
  field		description
  -----		-----------
  next		next entry at this bucket
  key		hash key (character string)
  ent		structure associated with this entry
*/

typedef struct hel_tag { struct hel_tag *next ;
			 const char *key ;
			 void *ent ; } hel ;

/* Hash table interface routines */

extern void *lookup(const char *key, hel *hashtab[], int size),
            *search(const char *key, hel *hashtab[], int size, bool init),
	    *enter(const char *key, hel *hashtab[], int size, void *entry),
	    *erase(const char *key, hel *hashtab[], int size) ;

#endif /* _DYLIB_HASH_H */

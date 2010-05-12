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
  This module contains routines for maintaining hash tables. Each entry has
  a character string key and a pointer to whatever data the user chooses to
  keep.
  
  The module is a standalone unit. I/O is handled entirely by stdio
  facilities unless the malloc debugging macros (dylib_std.h) are compiled in,
  in which case the io library is needed.
*/

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)hash.c	1.4	09/25/04" ;
static char svnid[] UNUSED = "$Id: dylib_hash.c 148 2007-06-09 03:15:30Z lou $" ;

#include <stdio.h>
#include "dylib_hash.h"



static int hash (const char *key, int size)

/*
  Function to return an index into the hash table.  The function simply adds up
  the values of the char's in the string.

  Parameter:
    key:	pointer to a character string
    size:	size of the hash table

  Return value: integer hash value of the string
*/

{ int hashval ;

  for (hashval = 0 ; *key != '\0' ; hashval += *key++) ;
  return (hashval%size) ; }



void *lookup (const char *key, hel *hashtab[], int size)

/*
  Function to search a hash table for the first occurence of a particular key.
  The hash index is calculated to get to the proper bucket, then the linked 
  list is searched sequentially. The insertion strategy insures that the most
  recently installed entry for the key is the one found. Since one might use
  this routine to see if something is in the table, it doesn't complain if the
  entry isn't found.

  Parameter:
    key:	character string to be looked up
    hashtab:	hash table to be searched
    size:	number of buckets in the hash table

  Returns: entry associated with the key, or NULL if nothing is found
*/

{ hel *hshel ;
  const char *rtnnme = "lookup" ;

  if (key == NULL)
  { fprintf(stderr,"\n%s: null key!\n", rtnnme) ;
    return (NULL) ; }
  if (hashtab == NULL)
  { fprintf(stderr,"\n%s: null hashtab!\n", rtnnme) ;
    return (NULL) ; }
  if (size <= 0)
  { fprintf(stderr,"\n%s: hashtab size violates 0 < %d!\n", rtnnme, size) ;
    return (NULL) ; }

  for (hshel = hashtab[hash(key,size)] ; hshel != NULL ; hshel = hshel->next)
    if (strcmp(key,hshel->key) == 0)
      return (hshel->ent) ;

  return (NULL) ; }



void *search (const char *key, hel *hashtab[], int size, bool init)

/*
  Function to search a hash table for all occurences of a particular key.
  If init == TRUE, then the hash index is calculated to get to the proper
  bucket, then the linked list is searched sequentially for the first occurence
  of the key. This position is remembered, and subsequent calls with init set
  to FALSE will continue searching the chain from the point where the last call
  left off.

  Parameter:
    key:	character string to be looked up
    hashtab:	hash table to be searched
    size:	number of buckets in the hash table
    init:	controls search method, as above.

  Returns: entry associated with the key, or NULL if nothing is found (this
	   will occur if the key is not present, or when there are no more
	   entries on the hash chain with the requested key).
*/

{ static hel *hshel ;
  static bool search_hshel_is_valid = FALSE ;
  const char *rtnnme = "search" ;

  if (key == NULL)
  { fprintf(stderr,"\n%s: null key!\n", rtnnme) ;
    search_hshel_is_valid = FALSE ;
    return (NULL) ; }
  if (hashtab == NULL)
  { fprintf(stderr,"\n%s: null hashtab!\n", rtnnme) ;
    search_hshel_is_valid = FALSE ;
    return (NULL) ; }
  if (size <= 0)
  { fprintf(stderr,"\n%s: hashtab size violates 0 < %d!\n", rtnnme, size) ;
    search_hshel_is_valid = FALSE ;
    return (NULL) ; }

  if (init == TRUE)
  { hshel = hashtab[hash(key,size)] ;
    search_hshel_is_valid = TRUE ; }
  else
  { if (search_hshel_is_valid != TRUE)
    { fprintf(stderr,"\n%s: attempt to continue before an init!\n",rtnnme) ;
      return (NULL) ; }
    if (hshel != NULL) hshel = hshel->next ; }

  for ( ; hshel != NULL ; hshel = hshel->next)
    if (strcmp(key,hshel->key) == 0) return (hshel->ent) ;

  search_hshel_is_valid = FALSE ;
  return (NULL) ; }



void *enter (const char *key, hel *hashtab[], int size, void *entry)

/*
  Function to add an entry to a hash table. The addition takes place at
  the head of the list (of the particular hash value). This is to facilitate
  the primary use of the hash tables, which is to return a pointer to the
  active symbol table entry for a name.

  Parameter:
    key:	character string to be used as a hash key
    hashtab:	hash table where entry is to be made
    size:	number of buckets in the hash table
    entry:	the entry to be associated with the key
  
  Returns: entry, or NULL if entry is not successfully entered in the table.
*/

{ hel *new ;
  int hashval ;
  const char *rtnnme = "enter" ;

  if (key == NULL)
  { fprintf(stderr,"\n%s: null key!\n", rtnnme) ;
    return (NULL) ; }
  if (hashtab == NULL)
  { fprintf(stderr,"\n%s: null hashtab!\n", rtnnme) ;
    return (NULL) ; }
  if (size <= 0)
  { fprintf(stderr,"\n%s: hashtab size violates 0 < %d!\n", rtnnme, size) ;
    return (NULL) ; }

  new = (hel *) MALLOC(sizeof(hel)) ;
  hashval = hash(key, size) ;
  new->next = hashtab[hashval] ;
  new->key = key ;
  new->ent = entry ;
  hashtab[hashval] = new ;
  return (entry) ; }




void *erase (const char *key, hel *hashtab[], int size)

/*
  Function used to delete a particular entry from a hash table. It erases the
  first hash element with a matching key. It is the caller's responsibility to
  deal with the entry pointed to by the hash element. It is considered an error
  to try and erase a non-existent entry, and the routine will complain.

  Parameter:
    key:	hash key
    hashtab:	hash table from which entry is to be deleted
    size:	number of buckets in the hash table

  Returns: Value of the entry field corresponding to the erased key, or NULL if
	   the key could not be found.
*/

{ int hashval ;
  hel *hshel1, **hshel2 ;
  void *entry ;
  const char *rtnnme = "erase" ;

  if (key == NULL)
  { fprintf(stderr,"\n%s: null key!\n", rtnnme) ;
    return (NULL) ; }
  if (hashtab == NULL)
  { fprintf(stderr,"\n%s: null hashtab!\n", rtnnme) ;
    return (NULL) ; }
  if (size <= 0)
  { fprintf(stderr,"\n%s: hashtab size violates 0 < %d!\n", rtnnme, size) ;
    return (NULL) ; }

  hashval = hash(key,size) ;
  for (hshel2 = &hashtab[hashval], hshel1 = *hshel2 ;
       hshel1 != NULL ;
       hshel2 = &hshel1->next, hshel1 = *hshel2)
    if (strcmp(key,hshel1->key) == 0)
    { *hshel2 = hshel1->next ;
      entry = hshel1->ent ;
      FREE(hshel1) ;
      return (entry) ; }

  fprintf(stderr,"\n%s: can't locate key %s.\n",rtnnme,key) ;
  return (NULL) ; }

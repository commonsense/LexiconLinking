#ifndef _DYLIB_KEYTAB_H
#define _DYLIB_KEYTAB_H

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
  Data structure for keyword tables searched by find and ambig
  
  @(#)keytab.h	1.2	08/31/99
  svn/cvs: $Id: dylib_keytab.h 148 2007-06-09 03:15:30Z lou $
*/

/*
  Field		Contents
  -----		--------
  keyword	Character string for the keyword.
  min		Minimum number of characters which must be matched before
		cimstrcmp will report a match.
  token		Value returned when the keyword is matched.
*/

typedef struct keytab_entry_internal { const char *keyword ;
				       int min ;
				       int token ; } keytab_entry ;


/*
  binsrch.c
*/

extern int find(char *word, keytab_entry keytab[], int numkeys),
	   ambig(char *word, keytab_entry keytab[], int numkeys) ;

#endif /* _DYLIB_KEYTAB_H */

#ifndef _DYLIB_STRRTNS_H
#define _DYLIB_STRRTNS_H

/*
  This file is part of the support library for the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

#include "dylib_std.h"

/*
  This file contains external definitions for the routines in the string
  package.

  @(#)strrtns.h	1.3	06/22/04
  svn/cvs: $Id: dylib_strrtns.h 148 2007-06-09 03:15:30Z lou $
*/

extern int cistrcmp(const char *str1, const char *str2),	/* strrtns.c */
           cimstrcmp(const char *str1, const char *str2),
           mstrcmp(const char *str1, const char *str2) ;
extern char *strsave(char *original) ;

extern const char *stralloc(const char *string) ;		/* littab.c */
extern bool strfree(const char *string) ;

/*
  Some macros to hide the memory allocation functions. Note that the
  debugging versions of these macros use outfmt from the io library
  and assume the existence of a string, rtnnme (typically the name of
  the current subroutine) that's used to identify the origin of the
  message.
*/

#if (MALLOC_DEBUG == 2)

#include "dylib_io.h"

void *zz_ptr_zz ;
ioid  zz_chn_zz ;

#define STRALLOC(zz_sptr_zz) \
  ( zz_ptr_zz = (void *) stralloc(zz_sptr_zz), \
    outfmt(zz_chn_zz,FALSE,":stralloc: %#08x (%s) in %s.\n", \
	   zz_ptr_zz,zz_ptr_zz,rtnnme), \
    (char *) zz_ptr_zz )

#define STRFREE(zz_fptr_zz) \
  ( outfmt(zz_chn_zz,FALSE,":strfree: %#08x (%s) in %s.\n", \
	   zz_fptr_zz,zz_fptr_zz,rtnnme), \
    strfree(zz_fptr_zz) )

#else

#define STRALLOC(zz_sptr_zz) stralloc(zz_sptr_zz)

#define STRFREE(zz_fptr_zz) strfree(zz_fptr_zz)

#endif


#endif /* _DYLIB_STRRTNS_H */

#ifndef _DYLIB_ERRS_H
#define _DYLIB_ERRS_H

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
  sccs: @(#)errs.h	2.3	03/18/04
  svn/cvs: $Id: dylib_errs.h 71 2006-06-09 04:21:15Z andreasw $
*/

#include "dylib_std.h"
#ifdef _DYLIB_FORTRAN
#include "dylib_fortran.h"
#endif

void errinit(const char *emsgpath, const char *elogpath, bool errecho),
     errterm(void) ;

void errmsg(int errid, ... ),
       warn(int errid, ... ) ;

#ifdef _DYLIB_FORTRAN
void errmsg_(integer *errid, char *ident, ... ) ;
void warn_(integer *errid, char *ident, ... ) ;
#endif

#endif /* _DYLIB_ERRS_H */

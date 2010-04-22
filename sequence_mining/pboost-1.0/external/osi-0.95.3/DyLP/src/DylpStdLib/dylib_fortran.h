#ifndef _DYLIB_FORTRAN_H
#define _DYLIB_FORTRAN_H
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
  @(#)fortran.h	1.1	09/01/99
  svn/cvs: $Id: dylib_fortran.h 71 2006-06-09 04:21:15Z andreasw $
*/

/*
  Common typedefs, definitions, macros, etc., which are handy when constructing
  C code that must talk to Fortran code.

  Off the top, typedefs and defines for the basic equivalences between
  Fortran and C data types. This list isn't complete, but it covers the
  common ones.  (Taken from the Sun Fortran Programmer's Guide.)
*/

typedef short int integer_2 ;
typedef long int integer ;
typedef long int logical ;
typedef float real ;
typedef double double_precision ;

#define TRUEL 1L
#define FALSEL 0L

/*
  A note about string handling in mixed code. C goes by the convention that
  strings are terminated by a '\0' character, and padded with '\0' on the
  rare occasions when padding is necessary. Fortran, on the other hand, keeps
  an explicit (though hidden) length, and pads with ' '. The two forms are
  just not compatible. Take care when passing strings back and forth. Adding
  an explicit null in the Fortran definition of the string is your best bet.
  The output routines in io.c and errs.c expect this, and do not make use of
  the 'hidden' parameter giving the string length.
*/

/*
  Some macros to help with Fortran arrays. 2-D arrays should simply be declared
  as array[rows*cols], and let the macro take care of the rest. This avoids
  complications due to column-major order in Fortran vs. row-major order in C.

  NOTE that these macros assume you are using the Fortran indexing convention,
  which starts at 1.
*/

#define f_chr(zz_ptr,zz_ndx,zz_strsze) (*(zz_ptr+((zz_ndx)-1)*(zz_strsze)))
#define f_arr1(zz_ptr,zz_ndx) (*(zz_ptr+(zz_ndx)-1))
#define f_arr2(zz_ptr,zz_row,zz_col,zz_collen) \
	(*(zz_ptr+((zz_col)-1)*(zz_collen)+((zz_row)-1)))


/*
  These codes are used by the Fortran part of the code to identify arguments
  supplied to errmsg_, warn_, and outfmt_. The Fortran code sees them from
  the common block argcod_, which is initialised in errs.c:errinit (from the
  io library. Do not rearrange the structure declaration without making
  corresponding changes in the Fortran common block.

  The codes ftnargVARNAME and ftnargCONNAME are peculiar to the bonsai MILP
  program (which prompted the development of the Fortran interface for i/o
  and error messages) and likely of little use in other contexts. At best,
  modifications to the routines in errs.c and io.c are required to support
  them.
*/

#define ftnargINTEGER ((integer) 1)
#define ftnargDOUBLE_PRECISION ((integer) 2)
#define ftnargCHARACTER ((integer) 3)
#define ftnargVARNAME ((integer) 4)
#define ftnargCONNAME ((integer) 5)
#define ftnargEND ((integer) 6)

extern struct { integer integer_code ;
		integer double_precision_code ;
		integer character_code ;
		integer varname_code ;
		integer conname_code ;
		integer end_code ; } argcod_ ;

#endif /* _DYLIB_FORTRAN_H */

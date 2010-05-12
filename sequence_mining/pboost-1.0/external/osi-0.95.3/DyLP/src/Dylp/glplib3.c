/* glplib3.c */

/*----------------------------------------------------------------------
-- Copyright (C) 2000, 2001, 2002 Andrew Makhorin <mao@mai2.rcnet.ru>,
--               Department for Applied Informatics, Moscow Aviation
--               Institute, Moscow, Russia. All rights reserved.
--
-- This file is a part of GLPK (GNU Linear Programming Kit).
--
-- GLPK is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2, or (at your option)
-- any later version.
--
-- GLPK is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
-- or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
-- License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with GLPK; see the file COPYING. If not, write to the Free
-- Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
-- 02111-1307, USA.
----------------------------------------------------------------------*/

#ifndef UNUSED
# if defined(_GNU_SOURCE) || defined(__GNUC__)
#   define UNUSED __attribute__((unused))
# else
#   define UNUSED
# endif
#endif

static char sccsid[] UNUSED = "@(#)glplib3.c	1.2	09/25/04" ;
static char svnid[] UNUSED = "$Id: glplib3.c 71 2006-06-09 04:21:15Z andreasw $" ;

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "glplib.h"

/*----------------------------------------------------------------------
-- print - print informative message.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void print(char *fmt, ...);
--
-- *Description*
--
-- The routine print prints an informative message specified by the
-- format control string fmt and optional parameter list. */

void print(const char *fmt, ...)
{     va_list arg;
      /* print an informative message */
      va_start(arg, fmt);
      vfprintf(stdout, fmt, arg);
      va_end(arg);
      fputc('\n', stdout);
      /* return to the calling program */
      return;
}

/*----------------------------------------------------------------------
-- fault - print error message and terminate program execution.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void fault(char *fmt, ...);
--
-- *Description*
--
-- The routine fault prints an error message specified by the format
-- control string fmt and optional parameter list, and then abnormally
-- terminates execution of the program.
--
-- *Returns*
--
-- The routine fault never returns. */

void fault(const char *fmt, ...)
{     va_list arg;
      /* print an error message */
      va_start(arg, fmt);
      vfprintf(stdout, fmt, arg);
      va_end(arg);
      fputc('\n', stdout);
      /* deinitialize library environment */
      free_lib_env();
      /* terminate program execution */
#if 0
      abort();
#else
      exit(3);
#endif
      /* no return */
}

/*----------------------------------------------------------------------
-- insist - check for logical condition.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void insist(int expr);
--
-- *Description*
--
-- The routine insist (implemented as a macro) checks for a logical
-- condition specified by the parameter expr. If the condition is false
-- (i.e. expr is zero), the routine prints an appropriate error message
-- and abnormally terminates the program.
--
-- This routine is a replacement of the standard function assert. */

void _insist(const char *expr, const char *file, int line)
{     /* print an error message */
      fputc('\n', stdout);
      fprintf(stdout, "Assertion failed: %s, file %s, line %d\n",
         expr, file, line);
      /* deinitialize library environment */
      free_lib_env();
      /* terminate program execution */
#if 0
      abort();
#else
      exit(3);
#endif
      /* no return */
}

/*----------------------------------------------------------------------
-- watch - take reading of stop-watch.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- double watch(void);
--
-- *Returns*
--
-- The routine watch returns the processor time in seconds. */

double watch(void)
{     return
         (double)clock() / (double)CLOCKS_PER_SEC;
}

/* eof */

/* glplib3.c */

/*----------------------------------------------------------------------
-- Copyright (C) 2000, 2001, 2002 Andrew Makhorin <mao@mai2.rcnet.ru>,
--               Department for Applied Informatics, Moscow Aviation
--               Institute, Moscow, Russia. All rights reserved.
--
-- This file is a part of GLPK (GNU Linear Programming Kit).
--
-- Licensed under the Common Public License (CPL) by permission of the
-- author for inclusion in the DyLP LP distribution.
----------------------------------------------------------------------*/

#ifndef UNUSED
# if defined(_GNU_SOURCE) || defined(__GNUC__)
#   define UNUSED __attribute__((unused))
# else
#   define UNUSED
# endif
#endif

static char sccsid[] UNUSED = "@(#)glplib3.c	1.2	09/25/04" ;
static char svnid[] UNUSED = "$Id: glplib3.c 148 2007-06-09 03:15:30Z lou $" ;

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef __CYGWIN__
/*
  With --pedantic-errors, cygwin won't compile its own signal.h, which is
  included from time.h
*/
# include <time.h>
#endif
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
#ifndef __CYGWIN__
/*
  As mentioned above, we get in trouble if we include time.h, so this
  function has to go.
*/
double watch(void)
{     return
         (double)clock() / (double)CLOCKS_PER_SEC;
}
#endif

/* eof */

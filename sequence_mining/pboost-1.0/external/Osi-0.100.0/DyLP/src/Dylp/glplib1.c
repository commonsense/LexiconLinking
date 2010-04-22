/* glplib1.c */

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

static char sccsid[] UNUSED = "@(#)glplib1.c	1.2	09/25/04" ;
static char svnid[] UNUSED = "$Id: glplib1.c 148 2007-06-09 03:15:30Z lou $" ;

#include <stddef.h>
#include "glplib.h"

static void *pointer = NULL;
/* a secret place to save a pointer */

/*----------------------------------------------------------------------
-- save_pointer - save a pointer.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void save_pointer(void *ptr);
--
-- *Description*
--
-- The routine save_pointer saves a pointer ptr in some secret place.
--
-- This routine is intended for internal needs and should not be used
-- directly by the application program. */

void save_pointer(void *ptr)
{     pointer = ptr;
      return;
}

/*----------------------------------------------------------------------
-- read_pointer - obtain a pointer.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void *read_pointer(void);
--
-- *Description*
--
-- The routine read_pointer obtains the pointer, which was previously
-- saved by the routine save_pointer in some secret place.
--
-- This routine is intended for internal needs and should not be used
-- directly by the application program.
--
-- *Returns*
--
-- The routine read_pointer returns the pointer saved by the routine
-- save_pointer. If the latter routine has not been called yet, NULL is
-- returned. */

void *read_pointer(void)
{     void *ptr = pointer;
      return ptr;
}

/* eof */

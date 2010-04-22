/* glplib1.c */

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

static char sccsid[] UNUSED = "@(#)glplib1.c	1.2	09/25/04" ;
static char svnid[] UNUSED = "$Id: glplib1.c 71 2006-06-09 04:21:15Z andreasw $" ;

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

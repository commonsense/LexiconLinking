/* glplib2.c */

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

static char sccsid[] UNUSED = "@(#)glplib2.c	1.2	09/25/04" ;
static char svnid[] UNUSED = "$Id: glplib2.c 71 2006-06-09 04:21:15Z andreasw $" ;

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "glplib.h"

/*----------------------------------------------------------------------
-- init_lib_env - initialize library environment.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- int init_lib_env(void);
--
-- *Description*
--
-- The routine init_lib_env creates and initializes library environment
-- block used by other low-level library routines.
--
-- If it is neccessary, this routine is called automatically, therefore
-- the user needn't to call it explicitly.
--
-- *Returns*
--
-- The routine returns one of the following error codes:
--
-- 0 - no errors;
-- 1 - the library environment has been already initialized;
-- 2 - initialization failed. */

int init_lib_env(void)
{     ENV *env;
      /* obtain a pointer to the environmental block */
      env = read_pointer();
      /* check if the environment has been initialized */
      if (env != NULL) return 1;
      /* allocate the environmental block */
      env = malloc(sizeof(ENV));
      /* check if the block has been successfully allocated */
      if (env == NULL) return 2;
      /* initialize the environmental block */
      env->mem_ptr = NULL;
      env->mem_limit = INT_MAX;
      env->mem_total = 0;
      env->mem_tpeak = 0;
      env->mem_count = 0;
      env->mem_cpeak = 0;
      /* save a pointer to the environmental block */
      save_pointer(env);
      /* initialization completed */
      return 0;
}

/*----------------------------------------------------------------------
-- get_env_ptr - obtain a pointer to the environmental block.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- ENV *get_env_ptr(void);
--
-- *Description*
--
-- The routine get_env_ptr obtains and returns a pointer to the library
-- environmental block.
--
-- If the library environment has not been initialized yet, the routine
-- performs initialization. If initialization fails, the routine prints
-- an error message to stderr and abnormally terminates the program.
--
-- *Returns*
--
-- The routine returns a pointer to the environmental block. */

ENV *get_env_ptr(void)
{     ENV *env;
      /* obtain a pointer to the environmental block */
      env = read_pointer();
      /* check if the environment has been initialized */
      if (env == NULL)
      {  /* not initialized yet; perform initialization */
         if (init_lib_env() != 0)
         {  /* initialization failed; print an error message */
            fprintf(stderr, "\nget_env_ptr: initialization failed\n");
            fflush(stderr);
            /* and abnormally terminate the program */
            abort();
         }
         /* initialization completed; obtain a pointer */
         env = read_pointer();
      }
      return env;
}

/*----------------------------------------------------------------------
-- free_lib_env - deinitialize library environment.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- int free_lib_env(void);
--
-- *Description*
--
-- The routine free_lib_env frees all resources (memory blocks, etc.),
-- which was allocated by the library routines and which are currently
-- still in use.
--
-- The user needn't to call this routine until he wishes to explicitly
-- free all the resources.
--
-- *Returns*
--
-- The routine returns one of the following codes:
--
-- 0 - no errors;
-- 1 - the library environment is inactive (not initialized);
-- 2 - deinitialization completed, but the routine was unable to free
--     some resources properly. */

int free_lib_env(void)
{     ENV *env;
      /* obtain a pointer to the environmental block */
      env = read_pointer();
      /* check if the environment has been initialized */
      if (env == NULL) return 1;
      /* free memory blocks, which are still allocated */
      while (env->mem_ptr != NULL)
      {  MEM *this = env->mem_ptr;
         env->mem_ptr = this->next;
         free(this);
      }
      /* free memory allocated to the environmental block */
      free(env);
      /* reset a pointer to the environmental block */
      save_pointer(NULL);
      /* deinitialization completed */
      return 0;
}

/* eof */

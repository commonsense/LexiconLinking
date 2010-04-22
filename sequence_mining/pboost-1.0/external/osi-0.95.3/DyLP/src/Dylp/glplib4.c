/* glplib4.c */

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

static char sccsid[] UNUSED = "@(#)glplib4.c	1.2	09/25/04" ;
static char svnid[] UNUSED = "$Id: glplib4.c 71 2006-06-09 04:21:15Z andreasw $" ;

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "glplib.h"

#define MEM_FLAG 0x20101960
/* a value used as memory block descriptor flag (may be changed if
   necessary) */

/*----------------------------------------------------------------------
-- umalloc - allocate memory block.
--
-- *Synopsis*
--
-- #include "glpset.h"
-- void *umalloc(int size);
--
-- *Description*
--
-- The routine umalloc allocates a memory block of size bytes long.
--
-- Note that being allocated the memory block contains arbitrary data
-- (not binary zeros).
--
-- *Returns*
--
-- The routine umalloc returns a pointer to the allocated memory block.
-- To free this block the routine ufree (not free!) should be used. */

void *umalloc(int size)
{     ENV *env = get_env_ptr();
      MEM *desc;
      int size_of_desc = align_datasize(sizeof(MEM));
      if (size < 1)
         fault("umalloc: invalid size");
      if (size > INT_MAX - size_of_desc)
         fault("umalloc: size too big");
      size += size_of_desc;
      if (size > env->mem_limit - env->mem_total)
         fault("umalloc: no memory available");
      desc = malloc(size);
      if (desc == NULL)
         fault("umalloc: malloc failed");
#if 1
      memset(desc, '?', size);
#endif
      desc->size = size;
      desc->flag = MEM_FLAG;
      desc->prev = NULL;
      desc->next = env->mem_ptr;
      if (desc->next != NULL) desc->next->prev = desc;
      env->mem_ptr = desc;
      env->mem_total += size;
      if (env->mem_tpeak < env->mem_total)
         env->mem_tpeak = env->mem_total;
      env->mem_count++;
      if (env->mem_cpeak < env->mem_count)
         env->mem_cpeak = env->mem_count;
      return (void *)((char *)desc + size_of_desc);
}

/*----------------------------------------------------------------------
-- ucalloc - allocate memory block.
--
-- *Synopsis*
--
-- #include "glpset.h"
-- void *ucalloc(int nmemb, int size);
--
-- *Description*
--
-- The routine ucalloc allocates a memory block of (nmemb*size) bytes
-- long.
--
-- Note that being allocated the memory block contains arbitrary data
-- (not binary zeros).
--
-- *Returns*
--
-- The routine ucalloc returns a pointer to the allocated memory block.
-- To free this block the routine ufree (not free!) should be used. */

void *ucalloc(int nmemb, int size)
{     if (nmemb < 1)
         fault("ucalloc: invalid nmemb");
      if (size < 1)
         fault("ucalloc: invalid size");
      if (nmemb > INT_MAX / size)
         fault("ucalloc: array too big");
      return umalloc(nmemb * size);
}

/*----------------------------------------------------------------------
-- ufree - free memory block.
--
-- *Synopsis*
--
-- #include "glpset.h"
-- void ufree(void *ptr);
--
-- *Description*
--
-- The routine ufree frees the memory block pointed to by ptr and which
-- was previuosly allocated by the routine umalloc or ucalloc. */

void ufree(void *ptr)
{     ENV *env = get_env_ptr();
      MEM *desc;
      int size_of_desc = align_datasize(sizeof(MEM));
      if (ptr == NULL)
         fault("ufree: null pointer");
      desc = (void *)((char *)ptr - size_of_desc);
      if (desc->flag != MEM_FLAG)
         fault("ufree: invalid pointer");
      if (env->mem_total < desc->size || env->mem_count == 0)
         fault("ufree: memory allocation error");
      if (desc->prev == NULL)
         env->mem_ptr = desc->next;
      else
         desc->prev->next = desc->next;
      if (desc->next == NULL)
         ;
      else
         desc->next->prev = desc->prev;
      env->mem_total -= desc->size;
      env->mem_count--;
      memset(desc, '?', size_of_desc);
      free(desc);
      return;
}

/* eof */

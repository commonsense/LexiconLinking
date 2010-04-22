/* glplib.h */

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
/*
  @(#)glplib.h	1.1	10/18/02
  svn/cvs: $Id: glplib.h 158 2007-07-06 01:25:14Z lou $
*/

#ifndef _GLPLIB_H
#define _GLPLIB_H

#define save_pointer          dy_glp_save_pointer
#define read_pointer          dy_glp_read_pointer

#define init_lib_env          dy_glp_init_lib_env
#define get_env_ptr           dy_glp_get_env_ptr
#define free_lib_env          dy_glp_free_lib_env

#define print                 dy_glp_print
#define fault                 dy_glp_fault
#define _insist               dy_glp_insist
#define watch                 dy_glp_watch

#define umalloc               dy_glp_umalloc
#define ucalloc               dy_glp_ucalloc
#define ufree                 dy_glp_ufree

#define create_pool           dy_glp_create_pool
#define get_atom              dy_glp_get_atom
#define free_atom             dy_glp_free_atom
#define get_atomv             dy_glp_get_atomv
#define clear_pool            dy_glp_clear_pool
#define delete_pool           dy_glp_delete_pool

extern void save_pointer(void *ptr);
/* save a pointer */

extern void *read_pointer(void);
/* obtain a pointer */

typedef struct ENV ENV;
typedef struct MEM MEM;

struct ENV
{     /* library environmental block */
      MEM *mem_ptr;
      /* pointer to the linked list of allocated memory blocks */
      int mem_limit;
      /* maximal amount of memory (in bytes) available for dynamic
         allocation */
      int mem_total;
      /* total amount of currently allocated memory (in bytes; is the
         sum of the size fields over all memory block descriptors) */
      int mem_tpeak;
      /* peak value of mem_total */
      int mem_count;
      /* total number of currently allocated memory blocks */
      int mem_cpeak;
      /* peak value of mem_count */
};

extern int init_lib_env(void);
/* initialize library environment */

extern ENV *get_env_ptr(void);
/* obtain a pointer to the environmental block */

extern int free_lib_env(void);
/* deinitialize library environment */

extern void print(const char *fmt, ...);
/* print informative message */

extern void fault(const char *fmt, ...);
/* print error message and terminate program execution */

#define insist(expr) \
((void)((expr) || (_insist(#expr, __FILE__, __LINE__), 1)))

extern void _insist(const char *expr, const char *file, int line);
/* check for logical condition */

extern double watch(void);
/* take reading of stop-watch */

/* some processors need data to be properly aligned; the align_boundary
   macro defines the boundary which should fit for all data types; the
   align_datasize macro allows enlarging size of data item in order the
   immediately following data of any type should be properly aligned */

#define align_boundary sizeof(double)

#define align_datasize(size) \
((((size) + (align_boundary - 1)) / align_boundary) * align_boundary)

struct MEM
{     /* memory block descriptor */
      int size;
      /* size of block (in bytes, including descriptor) */
      int flag;
      /* descriptor flag */
      MEM *prev;
      /* pointer to descriptor of the previous block */
      MEM *next;
      /* pointer to descriptor of the next block */
      /* actual data start here (there may be a "hole" between the next
         field and actual data because of data alignment) */
};

extern void *umalloc(int size);
/* allocate memory block */

extern void *ucalloc(int nmemb, int size);
/* allocate memory block */

extern void ufree(void *ptr);
/* free memory block */

typedef struct POOL POOL;

struct POOL
{     /* memory pool (a set of atoms) */
      int size;
      /* size of each atom in bytes (1 <= size <= 256); if size = 0,
         different atoms may have different sizes */
      void *avail;
      /* pointer to the linked list of free atoms */
      void *link;
      /* pointer to the linked list of allocated blocks (it points to
         the last recently allocated block) */
      int used;
      /* number of bytes used in the last allocated block */
      void *stock;
      /* pointer to the linked list of free blocks */
      int count;
      /* total number of allocated atoms */
};

extern POOL *create_pool(int size);
/* create memory pool */

extern void *get_atom(POOL *pool);
/* allocate atom of fixed size */

extern void free_atom(POOL *pool, void *ptr);
/* free an atom */

extern void *get_atomv(POOL *pool, int size);
/* allocate atom of variable size */

extern void clear_pool(POOL *pool);
/* free all atoms */

extern void delete_pool(POOL *pool);
/* delete memory pool */

#endif

/* eof */

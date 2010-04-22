# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 491 2006-05-19 20:37:52Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = dy_errmsgs.txt

########################################################################
#                               libDylp                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libDylp.la

# List all source files for this library, including headers
libDylp_la_SOURCES = \
	consys.h \
	consys_io.c \
	consys_mathutils.c \
	consys_scaling.c \
	consys_utils.c \
	dy_basis.c \
	dy_bound.c \
	dy_cmdint.c dy_cmdint.h \
	dy_coldstart.c \
	dy_conmgmt.c \
	dy_dual.c \
	dy_dualmultipivot.c \
	dy_dualpivot.c \
	dy_duenna.c \
	dy_force.c \
	dy_hotstart.c \
	dylp.c dylp.h \
	dylp_io.c \
	dylp_utils.c \
	dy_options.c \
	dy_penalty.c \
	dy_pivreject.c \
	dy_primal.c \
	dy_primalmultipivot.c \
	dy_primalpivot.c \
	dy_scaling.c \
	dy_setup.c \
	dy_statistics.c \
	dy_varmgmt.c \
	dy_warmstart.c \
	glpinv.c glpinv.h \
	glplib1.c \
	glplib2.c \
	glplib3.c \
	glplib4.c \
	glplib.h \
	glpluf.c glpluf.h \
	vector.h \
	vector_utils.c

# This is for libtool (on Windows)
libDylp_la_LDFLAGS = $(LT_LDFLAGS)

# We want to have also the objects from the DylpStdLib in this library
libDylp_la_LIBADD = ../DylpStdLib/libDylpStdLib.la

# Since automake is not doing this on its own, we need to declare the
# dependencies to the subdirectory libraries here
libDylp_la_DEPENDENCIES = $(libDylp_la_LIBADD)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../DylpStdLib`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	consys.h \
	dy_cmdint.h \
	dylp.h \
	glpinv.h \
	glplib.h \
	glpluf.h \
	vector.h

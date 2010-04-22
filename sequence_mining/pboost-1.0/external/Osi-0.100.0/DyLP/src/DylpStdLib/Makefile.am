# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 491 2006-05-19 20:37:52Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                            libDylpStdLib                             #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be installed
# since it will be included in the libDylp library
noinst_LTLIBRARIES = libDylpStdLib.la

# List all source files for this library, including headers
libDylpStdLib_la_SOURCES = \
	DylpConfig.h \
	dylib_binsrch.c \
	dylib_bnfrdr.c dylib_bnfrdr.h \
	dylib_bnfrdrio.c \
	dylib_errs.c dylib_errs.h \
	dylib_fortran.h \
	dylib_hash.c dylib_hash.h \
	dylib_io.c dylib_io.h \
	dylib_keytab.h \
	dylib_littab.c \
	dylib_std.h \
	dylib_strrtns.c dylib_strrtns.h

# This is for libtool (on Windows)
libDylpStdLib_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
CFLAGS +=

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../../inc/config_dylp.h \
	dylib_errs.h \
	dylib_hash.h \
	dylib_io.h \
	dylib_std.h \
	dylib_strrtns.h

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = DylpConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)

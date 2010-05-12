# Copyright (C) 2007-2009 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 1252 2009-03-31 21:34:31Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = get.GLPK INSTALL.Glpk glpk.patch

if COIN_BUILD_GLPK
  noinst_LTLIBRARIES = libcoinglpk.la
else
  noinst_LTLIBRARIES =
endif

GLPKDIR = @GLPKDIR@

nodist_libcoinglpk_la_SOURCES = \
$(GLPKDIR)/src/glpapi01.c \
$(GLPKDIR)/src/glpapi02.c \
$(GLPKDIR)/src/glpapi03.c \
$(GLPKDIR)/src/glpapi04.c \
$(GLPKDIR)/src/glpapi05.c \
$(GLPKDIR)/src/glpapi06.c \
$(GLPKDIR)/src/glpapi07.c \
$(GLPKDIR)/src/glpapi08.c \
$(GLPKDIR)/src/glpapi09.c \
$(GLPKDIR)/src/glpapi10.c \
$(GLPKDIR)/src/glpapi11.c \
$(GLPKDIR)/src/glpapi12.c \
$(GLPKDIR)/src/glpapi13.c \
$(GLPKDIR)/src/glpapi14.c \
$(GLPKDIR)/src/glpapi15.c \
$(GLPKDIR)/src/glpapi16.c \
$(GLPKDIR)/src/glpapi17.c \
$(GLPKDIR)/src/glpapi18.c \
$(GLPKDIR)/src/glpapi19.c \
$(GLPKDIR)/src/glpapi20.c \
$(GLPKDIR)/src/glpapi21.c \
$(GLPKDIR)/src/glpavl.c \
$(GLPKDIR)/src/glpbfd.c \
$(GLPKDIR)/src/glpbfx.c \
$(GLPKDIR)/src/glpcpx.c \
$(GLPKDIR)/src/glpdmp.c \
$(GLPKDIR)/src/glpdmx.c \
$(GLPKDIR)/src/glpfhv.c \
$(GLPKDIR)/src/glpgmp.c \
$(GLPKDIR)/src/glphbm.c \
$(GLPKDIR)/src/glpini01.c \
$(GLPKDIR)/src/glpini02.c \
$(GLPKDIR)/src/glpios01.c \
$(GLPKDIR)/src/glpios02.c \
$(GLPKDIR)/src/glpios03.c \
$(GLPKDIR)/src/glpios04.c \
$(GLPKDIR)/src/glpios05.c \
$(GLPKDIR)/src/glpios06.c \
$(GLPKDIR)/src/glpios07.c \
$(GLPKDIR)/src/glpios08.c \
$(GLPKDIR)/src/glpios09.c \
$(GLPKDIR)/src/glpios10.c \
$(GLPKDIR)/src/glpipm.c \
$(GLPKDIR)/src/glpipp01.c \
$(GLPKDIR)/src/glpipp02.c \
$(GLPKDIR)/src/glplib01.c \
$(GLPKDIR)/src/glplib02.c \
$(GLPKDIR)/src/glplib03.c \
$(GLPKDIR)/src/glplib04.c \
$(GLPKDIR)/src/glplib05.c \
$(GLPKDIR)/src/glplib06.c \
$(GLPKDIR)/src/glplib07.c \
$(GLPKDIR)/src/glplib08.c \
$(GLPKDIR)/src/glplib09.c \
$(GLPKDIR)/src/glplib10.c \
$(GLPKDIR)/src/glplib11.c \
$(GLPKDIR)/src/glplib12.c \
$(GLPKDIR)/src/glplpf.c \
$(GLPKDIR)/src/glplpp01.c \
$(GLPKDIR)/src/glplpp02.c \
$(GLPKDIR)/src/glplpx01.c \
$(GLPKDIR)/src/glplpx02.c \
$(GLPKDIR)/src/glplpx03.c \
$(GLPKDIR)/src/glplpx04.c \
$(GLPKDIR)/src/glplpx05.c \
$(GLPKDIR)/src/glpluf.c \
$(GLPKDIR)/src/glplux.c \
$(GLPKDIR)/src/glpmat.c \
$(GLPKDIR)/src/glpmpl01.c \
$(GLPKDIR)/src/glpmpl02.c \
$(GLPKDIR)/src/glpmpl03.c \
$(GLPKDIR)/src/glpmpl04.c \
$(GLPKDIR)/src/glpmpl05.c \
$(GLPKDIR)/src/glpmpl06.c \
$(GLPKDIR)/src/glpmps01.c \
$(GLPKDIR)/src/glpmps02.c \
$(GLPKDIR)/src/glpnet01.c \
$(GLPKDIR)/src/glpnet02.c \
$(GLPKDIR)/src/glpnet03.c \
$(GLPKDIR)/src/glpnet04.c \
$(GLPKDIR)/src/glpnet05.c \
$(GLPKDIR)/src/glpnet06.c \
$(GLPKDIR)/src/glpnet07.c \
$(GLPKDIR)/src/glpqmd.c \
$(GLPKDIR)/src/glprgr.c \
$(GLPKDIR)/src/glprng01.c \
$(GLPKDIR)/src/glprng02.c \
$(GLPKDIR)/src/glpscf.c \
$(GLPKDIR)/src/glpscg.c \
$(GLPKDIR)/src/glpscl.c \
$(GLPKDIR)/src/glpsds.c \
$(GLPKDIR)/src/glpspm.c \
$(GLPKDIR)/src/glpspx01.c \
$(GLPKDIR)/src/glpspx02.c \
$(GLPKDIR)/src/glpsql.c \
$(GLPKDIR)/src/glpssx01.c \
$(GLPKDIR)/src/glpssx02.c

libcoinglpk_la_DEPENDENCIES = $(srcdir)/$(GLPKDIR)/include/config.h

# CoinUtils and SYMPHONY include glpmpl.h and have HAVE_CONFIG_H defined
# so we put a dummy config.h file into the glpk include directory. At least
# a newline is needed to avoid a bug in Sun Solaris Studio, 64-bit optimised
# build.
$(srcdir)/$(GLPKDIR)/include/config.h :
	echo "" > $(srcdir)/$(GLPKDIR)/include/config.h

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/$(GLPKDIR)/include` -I$(GLPKDIR)

libcoinglpk_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = $(MAKEOKFILE)

include ../../BuildTools/Makemain.inc

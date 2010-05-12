# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 34 2006-05-15 19:22:07Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src/DylpStdLib src/Dylp examples

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball. Similarly, we don't want to attempt to build the
# documentation unless specifically requested.

if ALWAYS_FALSE
  SUBDIRS += test
  SUBDIRS += doc doc/Figures doc/TexMF
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	examples/generic.spc \
	examples/greenbeb.spc \
	examples/Makefile.in \
	examples/mpsio.c \
	examples/odsi+dylp.cpp \
	examples/osi_dylp.c \
	examples/osi+dylp.cpp \
	examples/plain \
	examples/README

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	@echo "Building and running dylp unit test."
	@cd test ; $(MAKE) test
	@if test -r ../Osi/test/Makefile; then \
	   echo "Building and running OsiDylp unit test." ; \
	   cd ../Osi/test; \
	   $(MAKE) test; \
	 else \
	   echo "Osi not present; skipping Osi unit test."; \
	 fi

# Autotools does a good job on code, but it doesn't deal well with LaTeX
# documentation. A few rules to help it along. Note that we don't delete
# dylp.ps and dylp.pdf until maintainer-clean; these come prebuilt with the
# distribution.

.PHONY: doc mostlyclean-local clean-local distclean-local \
	maintainer-clean-local
doc:
	@echo "Attempting to build dylp typeset documentation."
	@cd doc ; build_dylpdoc --ps --pdf dylp
mostlyclean-local:
	@cd doc ; if test -x build_dylpdoc ; then build_dylpdoc --clean ; fi

clean-local: mostlyclean-local

distclean-local: clean-local
	@cd doc ; rm -f build_dylpdoc dylpabsdir.tex makefile.dylpdoc
	@cd doc ; rm -f Makefile Figures/Makefile TexMF/Makefile

maintainer-clean-local: distclean-local
	@cd doc ; rm -f dylp.pdf dylp.ps

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

addlibsfile = dylp_addlibs.txt
addlibsdir = $(prefix)/share/doc/coin/$(PACKAGE_NAME)/

install-exec-local: install-doc
	$(install_sh_DATA) $(addlibsfile) $(DESTDIR)$(addlibsdir)/$(addlibsfile)

uninstall-local: uninstall-doc
	rm -f $(DESTDIR)$(addlibsdir)/$(addlibsfile)

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# So aclocal will find the local macros in the m4 directory

ACLOCAL_AMFLAGS = -I m4

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include ../BuildTools/Makemain.inc

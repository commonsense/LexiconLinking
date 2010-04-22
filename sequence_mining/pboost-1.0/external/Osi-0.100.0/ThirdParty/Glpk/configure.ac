# Copyright (C) 2007-2008 Stefan Vigerske and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: configure.ac 1252 2009-03-31 21:34:31Z stefan $

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([ThirdPartyGlpk],[0.0],[http://projects.coin-or.org/BuildTools/newticket])

AC_COPYRIGHT([
Copyright 2007-2009 Stefan Vigerske and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(Makefile.am)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AC_COIN_PROJECTDIR_INIT

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

#turn off warnings and pedantic mode since it seem to make trouble on solaris platform
#and we couldn't fix it anyway
coin_skip_warn_cflags=yes

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Glpk)

# Get the name of the C compiler
AC_COIN_PROG_CC

# Initialize automake
AC_COIN_INIT_AUTO_TOOLS

###########################################################################
#                            COIN_HAS_GLPK                                #
###########################################################################

# This macro checks for a library containing the GLPK library.  It
# checks if the user has provided an argument for the GLPK library,
# and if not, it checks whether the GLPK ThirdParty/GLPK directory has
# been configured.  It adds to ADDLIBS any flags required to link with
# an externally provided GLPK.  It defines the makefile conditional
# and preprocessor macro COIN_HAS_GLPK, if GLPK is available, and it
# defines the makefile conditional COIN_BUILD_GLPK, if GLPK is
# compiled within COIN.

## Check if GLPK is already provided by user
AC_COIN_HAS_GLPK(Glpk)

rm -f $MAKEOKFILE
if test x"$use_thirdpartyglpk" = xbuild ; then
  touch $MAKEOKFILE
  AC_MSG_CHECKING([GLPK version])
  majorver=`awk '/define GLP_MAJOR_VERSION/ {print $3}' ${GLPKINCDIR}/glpk.h`
  minorver=`awk '/define GLP_MINOR_VERSION/ {print $3}' ${GLPKINCDIR}/glpk.h`
  AC_MSG_RESULT([${majorver}.${minorver}])
  if test $majorver != 4 || test $minorver != 37 ; then
    AC_MSG_ERROR(["Incorrect GLPK version. This build system is made for GLPK 4.37. Please update by rerunning get.GLPK"])
  fi
else
  if test x$coin_has_glpk = xtrue ; then
    AC_MSG_NOTICE([The Glpk sources don't need to be compiled])
  else
    AC_MSG_WARN([Glpk source is available but will not be built. No valid external library is detected. Are you sure this is what you want?])
  fi
  MAKEOKFILE=
fi
AC_SUBST(MAKEOKFILE)

GLPKDIR=glpk
AC_SUBST(GLPKDIR)

AC_CONFIG_FILES([Makefile])

AC_COIN_FINALIZE

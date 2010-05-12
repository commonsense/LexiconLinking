# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: configure.ac,v 1.1.2.1 2006/04/19 23:25:12 andreasw Exp $

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([DyLP],[0.1],[coin-dylp@list.coin-or.org])


AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

AC_REVISION(0.9)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Dylp/vector.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(DyLP)

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Check for settings in CXXFLAGS that should move to CPPFLAGS
AC_DYLP_FIX_CPPFLAGS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(DyLP)

#############################################################################
#                       Doing project specific tests                        #
#############################################################################

# Determine the C type that corresponds to the C++ bool type in size

AC_DYLP_EQUIV_FOR_CPP_BOOL

########################################################################

# The floating point functions we'll be looking for seem to live in one of
# the following headers: math.h, ieeefp.h, float.h, and sunmath.h. Check
# for presence. We'll use the results to limit the work below. The coin
# macro COIN_PROJECTDIR_INIT should deal with presence/absence of -lm
# (for libmath) in ADDLIBS. Don't overrule it here.

AC_CHECK_HEADER([math.h],
    [AC_DEFINE([HAVE_MATH_H],[1],[Define to 1 if math.h exists.])
     fp_header_list="math.h" 
     AC_C_ADD_TO_INCLUDES([math.h])],
    [AC_MSG_WARN([It appears that math.h is missing.])
     AC_MSG_WARN([There is no hope of building dylp on this system.])])
AC_CHECK_HEADER([ieeefp.h],
    [AC_DEFINE([HAVE_IEEEFP_H],[1],[Define to 1 if ieeefp.h exists.])
     fp_header_list="ieeefp.h $fp_header_list"
     AC_C_ADD_TO_INCLUDES([ieeefp.h])],
    [])
AC_CHECK_HEADER([float.h],
    [AC_DEFINE([HAVE_FLOAT_H],[1],[Define to 1 if float.h exists.])
     fp_header_list="float.h $fp_header_list"
     AC_C_ADD_TO_INCLUDES([float.h])],
    [])
AC_CHECK_HEADER([sunmath.h],
    [AC_DEFINE([HAVE_SUNMATH_H],[1],[Define to 1 if sunmath.h exists.])
     fp_header_list="sunmath.h $fp_header_list"
     AC_C_ADD_TO_INCLUDES([sunmath.h])
     LIBS="-lsunmath $LIBS"
     ADDLIBS="-lsunmath $ADDLIBS"],
    [])

AC_MSG_NOTICE([Found floating point header files $fp_header_list.])
# AC_MSG_NOTICE([includes are now:])
# AC_MSG_NOTICE([$ac_includes_default])

# Now determine the names in this environment for isfinite() and isnan(). The
# most common variants are finite() and isnan(). We check for others. After
# these macros finish, DYLP_ISFINITE has the name for isfinite(), and
# DYLP_ISNAN has the name for isnan().

AC_DYLP_FIND_ISFINITE
AC_DYLP_FIND_ISNAN

# Check for the presence of a quiet_nan function. To my knowledge, this will
# exist only if sunmath is present, but hey, who knows what lurks out there.

AC_CHECK_FUNC([quiet_nan],
    [AC_DEFINE([DYLP_HAS_QUIET_NAN],[1],
	[Define this symbol if the quiet_nan function exists. This function
	 should return the bit pattern for IEEE quiet NaN.])])

# And determine whether we're big-endian or little-endian. This is necessary to
# define the correct patterns for IEEE infinity and quiet_nan, in the event
# that we can't get them any other way. This will define WORDS_BIGENDIAN if
# we're on a big-endian machine, and do nothing otherwise.

AC_C_BIGENDIAN

# Finally, dylp needs to know where it's error message text file is located.
# This is made horribly complicated by Windows, which insists on backslash
# as the path separator. Keep in mind the sed expression is interpreted by the
# shell and then by sed; the backslashes are halved each time. In the end
# we're left with s,\,\\,g as the actual substitution. This gets us a double
# backslash in DYLP_ERRMSGDIR, which the various compilers can strip down to
# one. Sheesh.

AC_ARG_VAR(DYLP_ERRMSGDIR,
    [Directory containing the DyLP error message text file (dy_errmsgs.txt).
     Default is srcdir/src/Dylp])
AC_SUBST(DYLP_ERRMSGDIR)
if test x"$DYLP_ERRMSGDIR" = x; then
  DYLP_ERRMSGDIR=`$CYGPATH_W $abs_source_dir/src/Dylp`
  DYLP_ERRMSGDIR=`echo $DYLP_ERRMSGDIR | sed -e 's,\\\\,\\\\\\\\,g'`
fi
AC_DEFINE_UNQUOTED([DYLP_ERRMSGDIR],
    ["$DYLP_ERRMSGDIR"],
    [Set to the full path directory name for the location of the error text
     message file dy_errmsgs.txt])

##############################################################################
#		Dylp configuration options				     #
##############################################################################

AC_DYLP_PARANOIA([no])
AC_DYLP_STATISTICS([no])
AC_DYLP_INFO([yes])

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

AC_COIN_VPATH_LINK(examples/plain)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/DylpStdLib/Makefile
                 src/Dylp/Makefile
		 dylp_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([inc/config_dylp.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE

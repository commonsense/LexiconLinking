# Copyright (C) 2006 International Business Machines..
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id: coin.m4 146 2006-11-29 18:27:51Z andreasw $
#
# Author: Andreas Wachter    IBM      2006-04-14

# This file defines the common autoconf macros for COIN
#

# Check requirements
AC_PREREQ(2.59)

###########################################################################
#                           COIN_MAIN_SUBDIRS                             #
###########################################################################

# This macro sets up the recursion into configure scripts into
# subdirectories.  Each possible subdirectory should be provided as a
# new argument to this macro.  The current limit is 10 subdirectories.
# This automatically also checks for the Data subdirectory.

AC_DEFUN([AC_COIN_MAIN_SUBDIRS],
[AC_ARG_VAR([COIN_SKIP_PROJECTS],[Set to the subdirectories of projects that should be skipped in the configuration])

m4_ifvaln([$1],[AC_MSG_CHECKING(whether directory $1 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $1; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$1/configure; then
                  coin_subdirs="$coin_subdirs $1"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($1)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$2],[AC_MSG_CHECKING(whether directory $2 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $2; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$2/configure; then
                  coin_subdirs="$coin_subdirs $2"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($2)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$3],[AC_MSG_CHECKING(whether directory $3 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $3; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$3/configure; then
                  coin_subdirs="$coin_subdirs $3"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($3)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$4],[AC_MSG_CHECKING(whether directory $4 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $4; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$4/configure; then
                  coin_subdirs="$coin_subdirs $4"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($4)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$5],[AC_MSG_CHECKING(whether directory $5 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $5; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$5/configure; then
                  coin_subdirs="$coin_subdirs $5"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($5)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$6],[AC_MSG_CHECKING(whether directory $6 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $6; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$6/configure; then
                  coin_subdirs="$coin_subdirs $6"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($6)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$7],[AC_MSG_CHECKING(whether directory $7 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $7; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$7/configure; then
                  coin_subdirs="$coin_subdirs $7"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($7)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$8],[AC_MSG_CHECKING(whether directory $8 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $8; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$8/configure; then
                  coin_subdirs="$coin_subdirs $8"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($8)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$9],[AC_MSG_CHECKING(whether directory $9 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $9; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$9/configure; then
                  coin_subdirs="$coin_subdirs $9"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($9)
                else
                  AC_MSG_RESULT(no)
                fi])
m4_ifvaln([$10],[AC_MSG_CHECKING(whether directory $10 is available)
                coin_skip=no
                if test x"$COIN_SKIP_PROJECTS" != x; then
                  for dir in $COIN_SKIP_PROJECTS; do
                    if test $dir = $10; then
                      coin_skip=yes
                    fi
                  done
                fi
                if test $coin_skip = yes; then
                  AC_MSG_RESULT(skipping)
                elif test -r $srcdir/$10/configure; then
                  coin_subdirs="$coin_subdirs $10"
                  AC_MSG_RESULT(yes)
                  AC_CONFIG_SUBDIRS($10)
                else
                  AC_MSG_RESULT(no)
                fi])
]) # AC_COIN_MAIN_SUBDIRS

###########################################################################
#                        COIN_THIRDPARTY_SUBDIRS                          #
###########################################################################

# This macro sets up the recursion into the configure script in a
# subdirectory for compilation of third party code.  The first
# argument is just the string that appears in the configure output.
# The second argument is the directory with the configure script, and
# the third one is a file that should exists in that directory.  If
# this file does not exist, we assume that the user has not downloaded
# the code, and we are not going to compile it

AC_DEFUN([AC_COIN_THIRDPARTY_SUBDIRS],
[AC_MSG_CHECKING(whether code for third party package $1 is available)
coin_skip=no
if test x"$COIN_SKIP_PROJECTS" != x; then
  for dir in $COIN_SKIP_PROJECTS; do
    if test $dir = $2; then
      coin_skip=yes
    fi
  done
fi
if test $coin_skip = yes; then
  AC_MSG_RESULT(skipping)
elif test -r $srcdir/$2/$3; then
  coin_subdirs="$coin_subdirs $2"
  AC_MSG_RESULT(yes)
  AC_CONFIG_SUBDIRS($2)
else
  AC_MSG_RESULT(no)
fi
]) # AC_COIN_THIRDPARTY_SUBDIRS

###########################################################################
#                           COIN_CHECK_VPATH                              #
###########################################################################

# This macro sets the variable coin_vpath_config to true if this is a
# VPATH configuration, otherwise it sets it to false.
AC_DEFUN([AC_COIN_CHECK_VPATH],
[AC_MSG_CHECKING(whether this is a VPATH configuration)
if test `cd $srcdir; pwd` != `pwd`; then
  coin_vpath_config=yes;
else
  coin_vpath_config=no;
fi
AC_MSG_RESULT($coin_vpath_config)
]) # AC_COIN_CHECK_VPATH

###########################################################################
#                         COIN_PROJECTDIR_INIT                            #
###########################################################################

# This macro does everything that is required in the early part in the
# configure script, such as defining a few variables.  This should only
# be used in the main directory of a project directory (the one under
# which src is)

AC_DEFUN([AC_COIN_PROJECTDIR_INIT],
[# Initialize the ADDLIBS variable
ADDLIBS='-lm'
AC_SUBST(ADDLIBS)

# Initialize the FADDLIBS variable (which is to be used with a fortran
# compiler and will not include FLIBS)
FADDLIBS=
AC_SUBST(FADDLIBS)

# A useful makefile conditional that is always false
AM_CONDITIONAL(ALWAYS_FALSE, false)

# We set the following variable so that we know later in AC_COIN_FINALIZE
# that we are in a project main directory
coin_projectdir=yes
]) # AC_COIN_PROJECTDIR_INIT

###########################################################################
#                          COIN_DEBUG_COMPILE                             #
###########################################################################

# enable the configure flags --enable-debug and --enable-debug-prjct
# (where prcjt is the name of the project in lower case) and set the
# variable coin_debug_compile to true or false This is used by
# COIN_PROG_CXX, COIN_PROG_CC and COIN_PROG_F77 to determine the
# compilation flags.  This macro also makes the switches
# --with-prjct-verbosity and --with-prjct-checklevel available, which
# define the preprocessor macros COIN_PRJCT_VERBOSITY and
# COIN_PRJCT_CHECKLEVEL to the specified value (default is 0).
#
# The project specific flags are only made available, if one gives the
# name of the project as first argument to this macro.

AC_DEFUN([AC_COIN_DEBUG_COMPILE],
[AC_BEFORE([$0],[AC_COIN_PROG_CXX])dnl
AC_BEFORE([$0],[AC_COIN_PROG_CC])dnl
AC_BEFORE([$0],[AC_COIN_PROG_F77])dnl

AC_MSG_CHECKING([whether we want to compile in debug mode])

AC_ARG_ENABLE([debug],
[AC_HELP_STRING([--enable-debug],
                [compile all projects with debug options tests])],
[case "${enableval}" in
   yes) coin_debug_compile=true
        enable_shared=no
        ;;
   no)  coin_debug_compile=false
        ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug)
        ;;
esac],
[coin_debug_compile=false])

m4_ifvaln([$1],
[AC_ARG_ENABLE(debug-m4_tolower($1),
 [AC_HELP_STRING([--enable-debug-m4_tolower($1)],
                 [compile this project ($1) with debug options])],
 [case "${enableval}" in
    yes) coin_debug_compile=true
         enable_shared=no
         ;;
    no)  coin_debug_compile=false
         ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug-m4_tolower($1))
         ;;
 esac],[:])
]) # m4_ifvaln([$1],

if test $coin_debug_compile = true; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

m4_ifvaln([$1],
[AC_ARG_WITH(m4_tolower($1)-verbosity,
             AC_HELP_STRING([--with-m4_tolower($1)-verbosity],
                            [specify the debug verbosity level for project $1]),
             [m4_tolower(coin_$1_verbosity)=$withval],
             [m4_tolower(coin_$1_verbosity)=0])
 AC_DEFINE_UNQUOTED(m4_toupper(COIN_$1_VERBOSITY),
                    m4_tolower($coin_$1_verbosity),
                    [Define to the debug verbosity level (0 is no output)])

 AC_ARG_WITH(m4_tolower($1)-checklevel,
             AC_HELP_STRING([--with-m4_tolower($1)-checklevel],
                            [specify the sanity check level for project $1]),
             [m4_tolower(coin_$1_checklevel)=$withval],
             [m4_tolower(coin_$1_checklevel)=0])
 AC_DEFINE_UNQUOTED(m4_toupper(COIN_$1_CHECKLEVEL),
                    m4_tolower($coin_$1_checklevel),
                    [Define to the debug sanity check level (0 is no test)])

# We use the following variable to have a string with the upper case
# version of the project name
COIN_PRJCT=m4_toupper($1)

]) # m4_ifvaln([$1],
 
]) # AC_COIN_DEBUG_COMPILE

###########################################################################
#                          COIN_MINGW_LD_FIX                              #
###########################################################################

# This macro is included by any PROG_compiler macro, to set the LD
# environment variable on MinWG to the correct value (link)

AC_DEFUN([AC_COIN_MINGW_LD_FIX],
[case $build in
  *-mingw*)
    if test "${LD+set}" = set; then :; else
      LD=link
    fi
    ;;
esac
])

###########################################################################
#                        COIN_ENABLE_DOSCOMPILE                           #
###########################################################################

# This macro is included by any PROG_compiler macro, to enable the
# --enable-doscompile options which is to be used when one wants to
# compile an executable under Cygwin which also runs directly under
# does (without requiring Cygwin1.dll).  Essentially, if enabled and
# the GNU compilers are used, it switches the --mno-cygwin flag on.

AC_DEFUN([AC_COIN_ENABLE_DOSCOMPILE],
[AC_ARG_ENABLE([doscompile],
[AC_HELP_STRING([--enable-doscompile],
                [Under Cygwin, compile so that executables run under DOS (default: disabled)])],
[if test "$enable_doscompile = yes"; then
  case $build in
    *-cygwin*) ;;
    *) AC_MSG_ERROR([--enable-doscompile option makes only sense under Cygwin]) ;;
  esac
fi],
[enable_doscompile=no])
])

###########################################################################
#                             COIN_PROG_CXX                               #
###########################################################################

# Find the compile command by running AC_PROG_CXX (with compiler names
# for different operating systems) and put it into CXX (unless it was
# given my the user), and find an appropriate value for CXXFLAGS.  It is
# possible to provide additional -D flags in the variable CXXDEFS.

AC_DEFUN([AC_COIN_PROG_CXX],
[AC_REQUIRE([AC_COIN_PROG_CC]) #Let's try if that overcomes configuration problem with VC++ 6.0
AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_LANG_PUSH(C++)

AC_ARG_VAR(CXXDEFS,[Additional -D flags to be used when compiling C++ code.])
AC_ARG_VAR(ADD_CXXFLAGS,[Additional C++ compiler options])
AC_ARG_VAR(DBG_CXXFLAGS,[Debug C++ compiler options])
AC_ARG_VAR(OPT_CXXFLAGS,[Optimize C++ compiler options])

coin_has_cxx=yes

save_cxxflags="$CXXFLAGS"
case $build in
  *-cygwin* | *-mingw*)
             comps="g++ cl" ;;
  *-darwin*) comps="g++ c++ CC" ;;
          *) comps="xlC aCC CC g++ c++ pgCC icpc gpp cxx cc++ cl FCC KCC RCC" ;;
esac

# We delete the cached value, since the test might not have been
# performed with out choise of compilers earlier
$as_unset ac_cv_prog_CXX || test "${ac_cv_prog_CXX+set}" != set || { ac_cv_prog_CXX=; export ac_cv_prog_CXX; }
AC_PROG_CXX([$comps])
CXXFLAGS="$save_cxxflags"

# Check if a project specific CXXFLAGS variable has been set
if test x$COIN_PRJCT != x; then
  eval coin_tmp=\${${COIN_PRJCT}_CXXFLAGS+set}
  if test x$coin_tmp = xset; then
    eval CXXFLAGS=\${${COIN_PRJCT}_CXXFLAGS}
  fi
fi

if test x"$CXXFLAGS" = x; then

# ToDo decide whether we want -DNDEBUG for optimization
  coin_add_cxxflags=
  coin_opt_cxxflags=
  coin_dbg_cxxflags=
  coin_warn_cxxflags=

  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc* | */icpc*)
        ;;
      *)
# ToDo decide about unroll-loops
        coin_opt_cxxflags="-O3 -fomit-frame-pointer"
        coin_add_cxxflags="-pipe"
        coin_dbg_cxxflags="-g"
        coin_warn_cxxflags="-pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion"
        if test "$enable_doscompile" = yes; then
          case $build in
            *-cygwin*)
              CXXFLAGS="-mno-cygwin"
              AC_TRY_LINK([],[int i=0; i++;],
                          [coin_add_cxxflags="-mno-cygwin $coin_add_cxxflags"])
              CXXFLAGS=
              ;;
          esac
        fi
    esac
  fi
  if test -z "$coin_opt_cxxflags"; then
    case $build in
      *-cygwin* | *-mingw*)
        case "$CXX" in
          cl* | */cl*)
            coin_opt_cxxflags='-O2'
            coin_add_cxxflags='-nologo -EHsc -GR -MT'
            coin_dbg_cxxflags='-Yd'
            ;;
        esac
        ;;
      *-linux-*)
        case "$CXX" in
          icpc* | */icpc*)
            coin_opt_cxxflags="-O3 -ip"
            coin_add_cxxflags=""
            coin_dbg_cxxflags="-g"
            # Check if -i_dynamic is necessary (for new glibc library)
            CXXFLAGS=
            AC_TRY_LINK([],[int i=0; i++;],[],
                        [coin_add_cxxflags="-i_dynamic $coin_add_cxxflags"])
            ;;
          pgCC* | */pgCC*)
            coin_opt_cxxflags="-fast"
            coin_add_cxxflags="-Kieee -pc 64"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-ibm-*)
        case "$CXX" in
          xlC* | */xlC* | mpxlC* | */mpxlC*)
            coin_opt_cxxflags="-O3 -qarch=auto -qcache=auto -qtune=auto -qmaxmem=-1"
            coin_add_cxxflags="-bmaxdata:0x80000000 -qrtti=dyna"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-hp-*)
        case "$CXX" in
          aCC* | */aCC* )
            coin_opt_cxxflags="-O"
            coin_add_cxxflags="-AA"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-sun-*)
          coin_opt_cxxflags="-O4"
          coin_dbg_cxxflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_cxx_g" = yes && test -z "$coin_dbg_cxxflags" ; then
    coin_dbg_cxxflags="-g"
  fi

  if test -z "$coin_opt_cxxflags"; then
    # Try if -O option works if nothing else is set
    CXXFLAGS=-O
    AC_TRY_LINK([],[int i=0; i++;],[coin_opt_cxxflags="-O"])
  fi

  # if PM doesn't want the warning messages, take them out
  if test x"$coin_skip_warn_cxxflags" = xyes; then
    coin_warn_cxxflags=
  fi

  if test x${DBG_CXXFLAGS+set} != xset; then
    DBG_CXXFLAGS="$coin_dbg_cxxflags $coin_add_cxxflags $coin_warn_cxxflags"
  fi
  if test x${OPT_CXXFLAGS+set} != xset; then
    OPT_CXXFLAGS="$coin_opt_cxxflags $coin_add_cxxflags -DNDEBUG $coin_warn_cxxflags"
  fi

  DBG_CXXFLAGS="$DBG_CXXFLAGS $ADD_CXXFLAGS $CXXDEFS"
  OPT_CXXFLAGS="$OPT_CXXFLAGS $ADD_CXXFLAGS $CXXDEFS"

  if test "$coin_debug_compile" = "true"; then
    CXXFLAGS="$DBG_CXXFLAGS"
  else
    CXXFLAGS="$OPT_CXXFLAGS"
  fi
else
  CXXFLAGS="$CXXFLAGS $ADD_CXXFLAGS $CXXDEFS"
  if test x${DBG_CXXFLAGS+set} != xset; then
    DBG_CXXFLAGS="$CXXFLAGS"
  fi
  if test x${OPT_CXXFLAGS+set} != xset; then
    OPT_CXXFLAGS="$CXXFLAGS"
  fi
fi

# Try if CXXFLAGS works
save_CXXFLAGS="$CXXFLAGS"
AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
if test -z "$CXXFLAGS"; then
  AC_MSG_WARN([The flags CXXFLAGS="$save_CXXFLAGS" do not work.  I will now just try '-O', but you might want to set CXXFLAGS manually.])
  CXXFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
  if test -z "$CXXFLAGS"; then
    AC_MSG_WARN([This value for CXXFLAGS does not work.  I will continue with empty CXXFLAGS, but you might want to set CXXFLAGS manually.])
  fi
fi

AC_MSG_NOTICE([C++ compiler options are: $CXXFLAGS])

AC_ARG_VAR(MPICXX,[C++ MPI Compiler])
if test x"$MPICXX" = x; then :; else
  AC_MSG_NOTICE([Will use MPI C++ compiler $MPICXX])
  CXX="$MPICXX"
fi

case "$CXX" in
  cl*)
    AC_COIN_MINGW_LD_FIX
    ;;
esac

AC_LANG_POP(C++)
]) # AC_COIN_PROG_CXX


###########################################################################
#                             COIN_CXXLIBS                                #
###########################################################################

# Determine the C++ runtime libraries required for linking a C++ library
# with a Fortran or C compiler.  The result is available in CXXLIBS.

AC_DEFUN([AC_COIN_CXXLIBS],
[AC_REQUIRE([AC_PROG_CXX])dnl
AC_LANG_PUSH(C++)
AC_ARG_VAR(CXXLIBS,[Libraries necessary for linking C++ code with Fortran compiler])
if test -z "$CXXLIBS"; then
  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc* | */icpc*)
        CXXLIBS=""
        ;;
      *)
        CXXLIBS="-lstdc++ -lm" # -lgcc"
        ;;
    esac
  else
    case $build in
     *-mingw32 | *-cygwin-* )
      case "$CXX" in
      cl*)
        CXXLIBS=nothing;;
      esac;;
     *-linux-*)
      case "$CXX" in
      icpc* | */icpc*)
        CXXLIBS=""
             ;;
      pgCC* | */pgCC*)
        CXXLIBS="-lstd -lC -lc"
             ;;
      esac;;
    *-ibm-*)
      CXXLIBS="-lC -lc"
      ;;
    *-hp-*)
      CXXLIBS="-L/opt/aCC/lib -l++ -lstd_v2 -lCsup_v2 -lm -lcl -lc"
      ;;
    *-sun-*)
      CXXLIBS="-lCstd -lCrun"
    esac
  fi
fi
if test -z "$CXXLIBS"; then
  AC_MSG_WARN([Could not automatically determine CXXLIBS (C++ link libraries; necessary if main program is in Fortran of C).])
else
  AC_MSG_NOTICE([Assuming that CXXLIBS is \"$CXXLIBS\".])
fi
if test x"$CXXLIBS" = xnothing; then
  CXXLIBS=
fi
AC_LANG_POP(C++)
]) # AC_COIN_CXXLIBS

###########################################################################
#                           COIN_CHECK_HEADER                             #
###########################################################################

# This macro checks for a header file, but it does so without the
# standard header.  This avoids warning messages like:
#
# configure: WARNING: dlfcn.h: present but cannot be compiled
# configure: WARNING: dlfcn.h:     check for missing prerequisite headers?
# configure: WARNING: dlfcn.h: see the Autoconf documentation
# configure: WARNING: dlfcn.h:     section "Present But Cannot Be Compiled"
# configure: WARNING: dlfcn.h: proceeding with the preprocessor's result
# configure: WARNING: dlfcn.h: in the future, the compiler will take precedence

AC_DEFUN([AC_COIN_CHECK_HEADER],
[if test x"$4" = x; then
  hdr="#include <$1>"
else
  hdr="$4"
fi
AC_CHECK_HEADERS([$1],[$2],[$3],[$hdr])
]) # AC_COIN_CHECK_HEADER

###########################################################################
#                       COIN_CHECK_CXX_CHEADER                             #
###########################################################################

# This macro checks for C header files that are used from C++.  For a give
# stub (e.g., math), it first checks if the C++ library (cmath) is available.
# If it is, it defines HAVE_CMATH (or whatever the stub is).  If it is not
# available, it checks for the old C head (math.h) and defines HAVE_MATH_H
# if that one exists.

AC_DEFUN([AC_COIN_CHECK_CXX_CHEADER],
[AC_LANG_PUSH(C++)
AC_COIN_CHECK_HEADER([c$1],[$2],[$3],[$4])
if test "$ac_cv_header_c$1" != "yes"; then
  AC_COIN_CHECK_HEADER([$1.h],[$2],[$3],[$4])
fi
AC_LANG_POP(C++)
]) # AC_COIN_CHECK_CXX_CHEADER

###########################################################################
#                             COIN_PROG_CC                                #
###########################################################################

# Find the compile command by running AC_PROG_CC (with compiler names
# for different operating systems) and put it into CC (unless it was
# given my the user), and find an appropriate value for CFLAGS.  It is
# possible to provide additional -D flags in the variable CDEFS.

AC_DEFUN([AC_COIN_PROG_CC],
[AC_REQUIRE([AC_COIN_MINGW_LD_FIX])
AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_LANG_PUSH(C)

# For consistency, we set the C compiler to the same value of the C++
# compiler, if the C++ is set, but the C compiler isn't (only for CXX=cl)
if test x"$CXX" != x; then
  case "$CXX" in
    cl*)
      if test x"$CC" = x; then
        CC="$CXX"
        AC_MSG_WARN([C++ compiler name provided as $CXX, but CC is unset. Setting CC to $CXX])
      fi
      ;;
  esac
fi

AC_ARG_VAR(CDEFS,[Additional -D flags to be used when compiling C code.])
AC_ARG_VAR(ADD_CFLAGS,[Additional C compiler options])
AC_ARG_VAR(DBG_CFLAGS,[Debug C compiler options])
AC_ARG_VAR(OPT_CFLAGS,[Optimize C compiler options])

coin_has_cc=yes

save_cflags="$CFLAGS"
case $build in
  *-cygwin* | *-mingw*)
             comps="gcc cl" ;;
  *-linux-*) comps="xlc gcc cc pgcc icc" ;;
  *)         comps="xlc_r xlc cc gcc pgcc icc" ;;
esac

# We delete the cached value, since the test might not have been
# performed with out choise of compilers earlier
$as_unset ac_cv_prog_CC || test "${ac_cv_prog_CC+set}" != set || { ac_cv_prog_CC=; export ac_cv_prog_CC; }
AC_PROG_CC([$comps])
CFLAGS="$save_cflags"

# Check if a project specific CFLAGS variable has been set
if test x$COIN_PRJCT != x; then
  eval coin_tmp=\${${COIN_PRJCT}_CFLAGS+set}
  if test x$coin_tmp = xset; then
    eval CFLAGS=\${${COIN_PRJCT}_CFLAGS}
  fi
fi

if test x"$CFLAGS" = x; then

  coin_add_cflags=
  coin_opt_cflags=
  coin_dbg_cflags=
  coin_warn_cflags=

  if test "$GCC" = "yes"; then
    case "$CC" in
      icc* | */icc*)
        ;;
      *)
        coin_opt_cflags="-O3 -fomit-frame-pointer"
        coin_add_cflags="-pipe"
        coin_dbg_cflags="-g"
        coin_warn_cflags="-pedantic-errors -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall"
        if test "$enable_doscompile" = yes; then
          case $build in
            *-cygwin*)
              CFLAGS="-mno-cygwin"
              AC_TRY_LINK([],[int i=0; i++;],
                          [coin_add_cflags="-mno-cygwin $coin_add_cflags"])
              CFLAGS=
            ;;
          esac
        fi
    esac
  fi
  if test -z "$coin_opt_cflags"; then
    case $build in
      *-cygwin* | *-mingw*)
        case "$CC" in
          cl* | */cl*)
            coin_opt_cflags='-O2'
            coin_add_cflags='-nologo'
            coin_dbg_cflags='-Yd'
            ;;
        esac
        ;;
      *-linux-*)
        case "$CC" in
          icc* | */icc*)
            coin_opt_cflags="-O3 -ip"
            coin_add_cflags=""
            coin_dbg_cflags="-g"
            # Check if -i_dynamic is necessary (for new glibc library)
            CFLAGS=
            AC_TRY_LINK([],[int i=0; i++;],[],
                        [coin_add_cflags="-i_dynamic $coin_add_cflags"])
            ;;
          pgcc* | */pgcc*)
            coin_opt_cflags="-fast"
            coin_add_cflags="-Kieee -pc 64"
            coin_dbg_cflags="-g"
            ;;
        esac
        ;;
      *-ibm-*)
        case "$CC" in
          xlc* | */xlc* | mpxlc* | */mpxlc*)
            coin_opt_cflags="-O3 -qarch=auto -qcache=auto -qtune=auto -qmaxmem=-1"
            coin_add_cflags="-bmaxdata:0x80000000"
            coin_dbg_cflags="-g"
          ;;
        esac
        ;;
      *-hp-*)
        coin_opt_cflags="-O"
        coin_add_cflags="-Ae"
        coin_dbg_cflags="-g"
        ;;
      *-sun-*)
        coin_opt_cflags="-xO4"
        coin_dbg_cflags="-g"
        ;;
      *-sgi-*)
        coin_opt_cflags="-O -OPT:Olimit=0"
        coin_dbg_cflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_cc_g" = yes && test -z "$coin_dbg_cflags" ; then
    coin_dbg_cflags="-g"
  fi

  if test -z "$coin_opt_cflags"; then
    # Try if -O option works if nothing else is set
    CFLAGS="-O"
    AC_TRY_LINK([],[int i=0; i++;],[coin_opt_cflags="-O"])
  fi

  # if PM doesn't want the warning messages, take them out
  if test x"$coin_skip_warn_cflags" = xyes; then
    coin_warn_cflags=
  fi

  if test x${DBG_CFLAGS+set} != xset; then
    DBG_CFLAGS="$coin_dbg_cflags $coin_add_cflags $coin_warn_cflags"
  fi
  if test x${OPT_CFLAGS+set} != xset; then
    OPT_CFLAGS="$coin_opt_cflags $coin_add_cflags -DNDEBUG $coin_warn_cflags"
  fi

  DBG_CFLAGS="$DBG_CFLAGS $ADD_CFLAGS $CDEFS"
  OPT_CFLAGS="$OPT_CFLAGS $ADD_CFLAGS $CDEFS"

  if test "$coin_debug_compile" = "true"; then
    CFLAGS="$DBG_CFLAGS"
  else
    CFLAGS="$OPT_CFLAGS"
  fi
else
  CFLAGS="$CFLAGS $ADD_CFLAGS $CDEFS"
  if test x${DBG_CFLAGS+set} != xset; then
    DBG_CFLAGS="$CFLAGS"
  fi
  if test x${OPT_CFLAGS+set} != xset; then
    OPT_CFLAGS="$CFLAGS"
  fi
fi

# Check if user wants to have additional CFLAGS options
AC_ARG_VAR(ADD_CFLAGS,[Additional C compiler options])
if test x"$ADD_CFLAGS" != x; then
  CFLAGS="$CFLAGS $ADD_CFLAGS"
fi

# Try if CFLAGS works
save_CFLAGS="$CFLAGS"
AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
if test -z "$CFLAGS"; then
  AC_MSG_WARN([The value CFLAGS="$save_CFLAGS" do not work.  I will now just try '-O', but you might want to set CFLAGS manually.])
  CFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
  if test -z "$CFLAGS"; then
    AC_MSG_WARN([This value for CFLAGS does not work.  I will continue with empty CFLAGS, but you might want to set CFLAGS manually.])
  fi
fi

AC_MSG_NOTICE([C compiler options are: $CFLAGS])

AC_ARG_VAR(MPICC,[C MPI Compiler])
if test x"$MPICC" = x; then :; else
  AC_MSG_NOTICE([Will use MPI C compiler $MPICC])
  CC="$MPICC"
fi

# Correct ADDLIBS initialization if we are using the MS compiler
case "$CC" in
  cl*)
    ADDLIBS=
    AC_COIN_MINGW_LD_FIX
    ;;
esac

AC_LANG_POP(C)
]) # AC_COIN_PROG_CC

###########################################################################
#                             COIN_PROG_F77                               #
###########################################################################

# Find the compile command by running AC_PROG_F77 (with compiler names
# for different operating systems) and put it into F77 (unless it was
# given my the user), and find an appropriate value for FFLAGS

AC_DEFUN([AC_COIN_PROG_F77],
[AC_REQUIRE([AC_COIN_MINGW_LD_FIX])
AC_REQUIRE([AC_COIN_ENABLE_DOSCOMPILE])
AC_REQUIRE([AC_COIN_PROG_CC])
AC_LANG_PUSH([Fortran 77])

AC_ARG_VAR(ADD_FFLAGS,[Additional Fortran compiler options])
AC_ARG_VAR(DBG_FFLAGS,[Debug Fortran compiler options])
AC_ARG_VAR(OPT_FFLAGS,[Optimize Fortran compiler options])

coin_has_f77=yes

save_fflags="$FFLAGS"
case $build in
  *-cygwin* | *-mingw*)
             comps="gfortran g77 ifort fl32" ;;
  *)         comps="xlf fort77 gfortran f77 g77 pgf90 pgf77 ifort ifc frt af77" ;;
esac

# We delete the cached value, since the test might not have been
# performed with out choise of compilers earlier
$as_unset ac_cv_prog_F77 || test "${ac_cv_prog_F77+set}" != set || { ac_cv_prog_F77=; export ac_cv_prog_F77; }
AC_PROG_F77($comps)
FFLAGS="$save_fflags"

# Check if a project specific FFLAGS variable has been set
if test x$COIN_PRJCT != x; then
  eval coin_tmp=\${${COIN_PRJCT}_FFLAGS+set}
  if test x$coin_tmp = xset; then
    eval FFLAGS=\${${COIN_PRJCT}_FFLAGS}
  fi
fi

if test x"$FFLAGS" = x; then

  coin_add_fflags=
  coin_opt_fflags=
  coin_dbg_fflags=
  coin_warn_fflags=

  if test "$G77" = "yes"; then
    coin_opt_fflags="-O3 -fomit-frame-pointer"
    coin_add_fflags="-pipe"
    coin_dbg_fflags="-g"
    if test "$enable_doscompile" = yes; then
      case $build in
        *-cygwin*)
          FFLAGS="-mno-cygwin"
          AC_TRY_LINK([],[      write(*,*) 'Hello world'],
                      [coin_add_fflags="-mno-cygwin $coin_add_fflags"])
          FFLAGS=
        ;;
      esac
    fi
  else
    case $build in
      *-cygwin* | *-mingw*)
        case $F77 in
          ifort* | */ifort*)
            coin_opt_fflags='-O3'
            coin_add_fflags='-nologo -MT'
            coin_dbg_fflags='-debug'
          ;;
        esac
        ;;
      *-linux-*)
        case $F77 in
          ifc* | */ifc* | ifort* | */ifort*)
            coin_opt_fflags="-O3 -ip"
            coin_add_fflags="-cm -w90 -w95"
            coin_dbg_fflags="-g -CA -CB -CS"
            # Check if -i_dynamic is necessary (for new glibc library)
            FFLAGS=
            AC_TRY_LINK([],[      write(*,*) 'Hello world'],[],
                        [coin_add_fflags="-i_dynamic $coin_add_fflags"])
            ;;
          pgf77* | */pgf77* | pgf90* | */pgf90*)
            coin_opt_fflags="-fast"
            coin_add_fflags="-Kieee -pc 64"
            coin_dbg_fflags="-g"
          ;;
        esac
        ;;
      *-ibm-*)
        case "$F77" in
          xlf* | */xlf* | mpxlf* | */mpxlf* )
            coin_opt_fflags="-O3 -qarch=auto -qcache=auto -qtune=auto -qmaxmem=-1"
            coin_add_fflags="-bmaxdata:0x80000000"
            coin_dbg_fflags="-g -C"
            ;;
        esac
        ;;
      *-hp-*)
        coin_opt_fflags="+O3"
        coin_add_fflags="+U77"
        coin_dbg_fflags="-C -g"
        ;;
      *-sun-*)
        coin_opt_fflags="-O4"
        coin_dbg_fflags="-g"
        ;;
      *-sgi-*)
        coin_opt_fflags="-O5 -OPT:Olimit=0"
        coin_dbg_fflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_f77_g" = yes && test -z "$coin_dbg_fflags" ; then
    coin_dbg_fflags="-g"
  fi

  if test -z "$coin_opt_fflags"; then
    # Try if -O option works if nothing else is set
    FFLAGS=-O
    AC_TRY_LINK([],[      integer i], [coin_opt_fflags="-O"])
  fi

  # if PM doesn't want the warning messages, take them out
  if test x"$coin_skip_warn_fflags" = xyes; then
    coin_warn_fflags=
  fi

  if test x${DBG_FFLAGS+set} != xset; then
    DBG_FFLAGS="$coin_dbg_fflags $coin_add_fflags $coin_warn_fflags"
  fi
  if test x${OPT_FFLAGS+set} != xset; then
    OPT_FFLAGS="$coin_opt_fflags $coin_add_fflags $coin_warn_fflags"
  fi

  DBG_FFLAGS="$DBG_FFLAGS $ADD_FFLAGS"
  OPT_FFLAGS="$OPT_FFLAGS $ADD_FFLAGS"

  if test "$coin_debug_compile" = "true"; then
    FFLAGS="$DBG_FFLAGS"
  else
    FFLAGS="$OPT_FFLAGS"
  fi
else
  FFLAGS="$FFLAGS $ADD_FFLAGS"
  if test x${DBG_FFLAGS+set} != xset; then
    DBG_FFLAGS="$FFLAGS"
  fi
  if test x${OPT_FFLAGS+set} != xset; then
    OPT_FFLAGS="$FFLAGS"
  fi
fi

# Try if FFLAGS works
AC_TRY_LINK([],[      integer i],[],[FFLAGS=])
if test -z "$FFLAGS"; then
  AC_MSG_WARN([The flags FFLAGS="$FFLAGS" do not work.  I will now just try '-O', but you might want to set FFLAGS manually.])
  FFLAGS='-O'
  AC_TRY_LINK([],[      integer i],[],[FFLAGS=])
  if test -z "$FFLAGS"; then
    AC_MSG_WARN([This value for FFLAGS does not work.  I will continue with empty FFLAGS, but you might want to set FFLAGS manually.])
  fi
fi

AC_MSG_NOTICE([Fortran compiler options are: $FFLAGS])

AC_ARG_VAR(MPIF77,[Fortran MPI Compiler])
if test x"$MPIF77" = x; then :; else
  AC_MSG_NOTICE([Will use MPI Fortran compiler $MPIF77])
  F77="$MPIF77"
fi

case "$F77" in
  ifort*)
    AC_COIN_MINGW_LD_FIX
    ;;
esac

AC_LANG_POP([Fortran 77])
]) # AC_COIN_PROG_F77

###########################################################################
#                           COIN_F77_WRAPPERS                             #
###########################################################################

# Calls autoconfs AC_F77_WRAPPERS and does additional corrections to FLIBS

AC_DEFUN([AC_COIN_F77_WRAPPERS],
[AC_BEFORE([AC_COIN_PROG_F77],[$0])dnl
AC_BEFORE([AC_PROG_F77],[$0])dnl

AC_LANG_PUSH([Fortran 77])

AC_F77_WRAPPERS

# This is to correct a missing exclusion in autoconf 2.59
if test x"$FLIBS" != x; then
  my_flibs=
  for flag in $FLIBS; do
    case flag in
      -lcrt*.o) ;;
             *) my_flibs="$my_flibs $flag" ;;
    esac
  done
  FLIBS="$my_flibs"
fi

case $build in
# The following is a fix to define FLIBS for ifort on Windows
   *-cygwin* | *-mingw*)
     case "$F77" in
       ifort* | */ifort*)
           FLIBS="-link libifcorert.lib $LIBS /NODEFAULTLIB:libc.lib";;
     esac;;
   *-hp-*)
       FLIBS="$FLIBS -lm";;
   *-ibm-*)
       FLIBS=`echo $FLIBS | sed 's/-lc)/-lc/g'` ;;
   *-linux-*)
     case "$F77" in
       pgf77* | */pgf77* | pgf90* | */pgf90*)
# ask linker to go through the archives multiple times
# (the Fortran compiler seems to do that automatically...
         FLIBS="-Wl,--start-group $FLIBS -Wl,--end-group" ;;
     esac
esac

]) # AC_COIN_F77_WRAPPERS


###########################################################################
#                          COIN_INIT_AUTOMAKE                             #
###########################################################################

# This macro calls the regular INIT_AUTOMAKE and MAINTAINER_MODE
# macros, and defines additional variables and makefile conditionals,
# that are used in the maintainer parts of the makfile.  It also
# checks if the correct versions of the autotools are used.
#
# This also defines the AC_SUBST variables:
# abs_source_dir     absolute path to source code for this package
# abs_bin_dir        absolute path to the directory where binaries are
#                    going to be installed (prefix/bin)
# abs_lib_dir        absolute path to the directory where libraries are
#                    going to be installed (prefix/lib)
# abs_include_dir    absolute path to the directory where the header files
#                    are installed (prefix/include)

AC_DEFUN([AC_COIN_INIT_AUTOMAKE],
[AC_REQUIRE([AC_PROG_EGREP])
# Stuff for automake
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

coin_have_externals=no
if test "$enable_maintainer_mode" = yes; then

  # If maintainer mode is chosen, we make sure that the correct versions
  # of the tools are used, and that we know where libtoo.m4 is (to
  # recreate acinclude.m4)

  AC_SUBST(LIBTOOLM4)
  LIBTOOLM4=

  # Check if we have autoconf
  AC_CHECK_PROG([have_autoconf],[autoconf],[yes],[no])
  if test $have_autoconf = no; then
    AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find autoconf in your path.])
  fi

  # Check whether autoconf is the correct version
  correct_version='2.59'
  grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
  AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of autoconf])
  autoconf --version > confauto.out 2>&1
  if $EGREP $grep_version confauto.out >/dev/null 2>&1; then
    AC_MSG_RESULT([yes])
  else
    rm -f confauto.out
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([You don't have the correct version of autoconf as the first one in your path.])
  fi
  rm -f confauto.out

  # Check if the executable autoconf is picked up from the correct location
  AC_MSG_CHECKING([whether autoconf is coming from the correct location])
  autoconf_dir=`which autoconf | sed -e 's=/autoconf=='`
  autoconf_dir=`cd $autoconf_dir; pwd`
  if test x$AUTOTOOLS_DIR = x; then
    want_dir=$HOME/bin
  else
    want_dir=$AUTOTOOLS_DIR/bin
  fi
  if test $autoconf_dir = `cd $want_dir; pwd`; then
    AC_MSG_RESULT([yes])
  else
    rm -f confauto.out
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([The autoconf executable should be picked up from \$HOME/bin or \$AUTOTOOLS_DIR/bin.])
  fi

  # Check if we have automake
  AC_CHECK_PROG([have_automake],[automake],[yes],[no])
  if test $have_automake = no; then
    AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find automake in your path.])
  fi
  
  # Check whether automake is the correct version
  correct_version='1.9.6'
  grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
  AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of automake])
  automake --version > confauto.out 2>&1
  if $EGREP $grep_version confauto.out >/dev/null 2>&1; then
    AC_MSG_RESULT([yes])
  else
    rm -f confauto.out
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([You don't have the correct version of automake as the first one in your path.])
  fi
  rm -f confauto.out

  # Check if the executable automake is picked up from the correct location
  AC_MSG_CHECKING([whether automake is coming from the correct location])
  automake_dir=`which automake | sed -e 's=/automake=='`
  automake_dir=`cd $automake_dir; pwd`
  if test x$AUTOTOOLS_DIR = x; then
    want_dir=$HOME/bin
  else
    want_dir=$AUTOTOOLS_DIR/bin
  fi
  if test $automake_dir = `cd $want_dir; pwd`; then
    AC_MSG_RESULT([yes])
  else
    rm -f confauto.out
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([The automake executable should be picked up from \$HOME/bin or \$AUTOTOOLS_DIR/bin.])
  fi

  # Check if we can find the libtool file
  if test x$AUTOTOOLS_DIR = x; then
    want_dir=$HOME/share
  else
    want_dir=$AUTOTOOLS_DIR/share
  fi
  AC_CHECK_FILE([$want_dir/aclocal/libtool.m4],
                [LIBTOOLM4="$want_dir/aclocal/libtool.m4"],
                [AC_MSG_ERROR([I cannot find the libtool.m4 file.])])

  # Check if this is the correct version of libtool (with escaped dots)
  correct_version='1.5.22'
  grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
  AC_CHECK_FILE([$want_dir/libtool/ltmain.sh],
	        [have_ltmain=yes],
                [have_ltmain=no])
  AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of libtool.])
  if test $have_ltmain = yes; then
    if $EGREP $grep_version $want_dir/libtool/ltmain.sh >/dev/null 2>&1; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([You don't have the correct version of libtool.])
    fi
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([I cannot find the ltmain.sh file.])
  fi  

  # Check if we have an Externals file
  if test -r $srcdir/Externals; then
    coin_have_externals=yes
  fi
  # Check if subversion is installed and understands https
  AC_CHECK_PROG([have_svn],[svn],[yes],[no])
  if test x$have_svn = xyes; then
    AC_MSG_CHECKING([svn understands https])
    svn --version > confauto.out 2>&1
    if $EGREP https confauto.out >/dev/null 2>&1; then
      AC_MSG_RESULT(yes)
    else
      AC_MSG_RESULT(no)
      have_svn=no
    fi
    rm -f confauto.out
  fi

  # Find the location of the BuildTools directory
  BUILDTOOLSDIR=
  if test -r $srcdir/BuildTools/coin.m4; then
    BUILDTOOLSDIR=$srcdir/BuildTools
  else
    if test -r $srcdir/../BuildTools/coin.m4; then
      BUILDTOOLSDIR=$srcdir/../BuildTools
    else
      if test -r $srcdir/../../BuildTools/coin.m4; then
        BUILDTOOLSDIR=$srcdir/../../BuildTools
      else
        AC_MSG_ERROR(Cannot find the BuildTools directory)
      fi
    fi
  fi
  AC_SUBST(BUILDTOOLSDIR)

  # The following variable is set to the name of the directory where
  # the autotool scripts are located
  AC_SUBST(AUX_DIR)
  AUX_DIR=$ac_aux_dir
fi

# helpful variable for the base directory of this package
abs_source_dir=`cd $srcdir; pwd`
AC_SUBST(abs_source_dir)

# Stuff for example Makefiles
if test x$prefix = xNONE; then
  full_prefix=$ac_default_prefix
else
  full_prefix=$prefix
fi
full_prefix=`cd $full_prefix ; pwd`
AC_SUBST(abs_lib_dir)
abs_lib_dir=$full_prefix/lib
AC_SUBST(abs_include_dir)
abs_include_dir=$full_prefix/include
AC_SUBST(abs_bin_dir)
abs_bin_dir=$full_prefix/bin

AM_CONDITIONAL(HAVE_EXTERNALS,
               test $coin_have_externals = yes && test x$have_svn = xyes)
]) # AC_COIN_INIT_AUTOMAKE

###########################################################################
#                         COIN_INIT_AUTO_TOOLS                            #
###########################################################################

# Initialize the auto tools automake and libtool, with all
# modifications we want for COIN packages.
#
# RPATH_FLAGS        link flags for hardcoding path to shared objects

# This is a trick to have this code before AC_COIN_PROG_LIBTOOL
AC_DEFUN([AC_COIN_DISABLE_STATIC],
[
# On Cygwin, building DLLs doesn't work
case $build in
  *-cygwin*)
    coin_disable_shared=yes
    platform=Cygwin
  ;;
  *-mingw*)
    coin_disable_shared=yes
    platform="Msys"
#    case "$CXX" in
#      cl*)
#        coin_disable_shared=yes
#        platform="Msys with cl"
#    ;;
#    esac
  ;;
esac
if test x"$coin_disable_shared" = xyes; then
  if test x"$enable_shared" = xyes; then
    AC_MSG_WARN([On $platform, shared objects are not supported. I'm disabling your choice.])
  fi
  enable_shared=no
fi
# By default, we only want the shared objects to be compiled
AC_DISABLE_STATIC
])

AC_DEFUN([AC_COIN_INIT_AUTO_TOOLS],
[AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_REQUIRE([AC_COIN_DISABLE_STATIC])

# Initialize automake
AC_COIN_INIT_AUTOMAKE

# Stuff for libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects
AC_COIN_RPATH_FLAGS($abs_lib_dir)
]) # AC_COIN_INIT_AUTO_TOOLS

###########################################################################
#                           COIN_PROG_LIBTOOL                             #
###########################################################################

# Setup the libtool stuff together with any modifications to make it
# work on additional platforms

AC_DEFUN([AC_COIN_PROG_LIBTOOL],
[AC_REQUIRE([AC_COIN_DLFCN_H])

# We check for this header here in a non-standard way to avoid warning
# messages
AC_PROG_LIBTOOL

# Fix bugs in libtool script for Windows native compilation:
# - cygpath is not correctly quoted in fix_srcfile_path
# - paths generated for .lib files is not run through cygpath -w


# - lib includes subdirectory information; we want to replace
#
# old_archive_cmds="lib /OUT:\$oldlib\$oldobjs\$old_deplibs"
#
# by
#
# old_archive_cmds="echo \$oldlib | grep .libs >/dev/null; if test \$? = 0; then cd .libs; lib /OUT:\`echo \$oldlib\$oldobjs\$old_deplibs | sed -e s@\.libs/@@g\`; cd .. ; else lib /OUT:\$oldlib\$oldobjs\$old_deplibs ; fi"
#
#          -e 's%old_archive_cmds="lib /OUT:\\\$oldlib\\\$oldobjs\\\$old_deplibs"%old_archive_cmds="echo \\\$oldlib \| grep .libs >/dev/null; if test \\\$? = 0; then cd .libs; lib /OUT:\\\`echo \\\$oldlib\\\$oldobjs\\\$old_deplibs \| sed -e s@\\.libs/@@g\\\`; cd .. ; else lib /OUT:\\\$oldlib\\\$oldobjs\\\$old_deplibs; fi"%' \

# The following was a hack for chaniing @BACKSLASH to \
#          -e 'sYcompile_command=`\$echo "X\$compile_command" | \$Xsed -e '"'"'s%@OUTPUT@%'"'"'"\$output"'"'"'%g'"'"'`Ycompile_command=`\$echo "X\$compile_command" | \$Xsed -e '"'"'s%@OUTPUT@%'"'"'"\$output"'"'"'%g'"'"' | \$Xsed -e '"'"'s%@BACKSLASH@%\\\\\\\\\\\\\\\\%g'"'"'`Y' \

# Correct cygpath for minGW (ToDo!)
case $build in
  *-mingw*)
    CYGPATH_W=echo
    ;;
esac

case $build in
  *-cygwin* | *-mingw*)
  case "$CXX" in
    cl* | */cl*) 
      AC_MSG_NOTICE(Applying patches to libtool for cl compiler)
      sed -e 's|fix_srcfile_path=\"`cygpath -w \"\$srcfile\"`\"|fix_srcfile_path=\"\\\`'"$CYGPATH_W"' \\\"\\$srcfile\\\"\\\`\"|' \
          -e 's|fix_srcfile_path=\"\"|fix_srcfile_path=\"\\\`'"$CYGPATH_W"' \\\"\\$srcfile\\\"\\\`\"|' \
          -e 's%compile_deplibs=\"\$dir/\$old_library \$compile_deplibs\"%compile_deplibs="'\`"$CYGPATH_W"' \$dir/\$old_library | sed -e '"'"'sY\\\\\\\\Y/Yg'"'"\`' \$compile_deplibs\"'% \
          -e 's%compile_deplibs=\"\$dir/\$linklib \$compile_deplibs\"%compile_deplibs="'\`"$CYGPATH_W"' \$dir/\$linklib | sed -e '"'"'sY\\\\\\\\Y/Yg'"'"\`' \$compile_deplibs\"'% \
	  -e 's%lib /OUT:%lib -OUT:%' \
	  -e "s%cygpath -w%$CYGPATH_W%" \
	  -e 's%$AR x \\$f_ex_an_ar_oldlib%bla=\\`lib -nologo -list \\$f_ex_an_ar_oldlib | xargs echo\\`; echo \\$bla; for i in \\$bla; do lib -nologo -extract:\\$i \\$f_ex_an_ar_oldlib; done%' \
	  -e 's/$AR t/lib -nologo -list/' \
	  -e 's%f_ex_an_ar_oldlib="\($?*1*\)"%f_ex_an_ar_oldlib='\`"$CYGPATH_W"' \1`%' \ 
	  -e  's%^archive_cmds=.*%archive_cmds="\\$CC -o \\$lib \\$libobjs \\$compiler_flags \\\\\\`echo \\\\\\"\\$deplibs\\\\\\" | \\$SED -e '"\'"'s/ -lc\\$//'"\'"'\\\\\\` -link -dll~linknames="%' \
      libtool > conftest.bla

      mv conftest.bla libtool
      chmod 755 libtool
      ;;
    *)
      AC_MSG_NOTICE(Applying patches to libtool for GNU compiler)
      sed -e 's|fix_srcfile_path=\"`cygpath -w \"\$srcfile\"`\"|fix_srcfile_path=\"\\\`'"$CYGPATH_W"' \\\"\\$srcfile\\\"\\\`\"|' \
          -e 's|"lib /OUT:\\$oldlib\\$oldobjs\\$old_deplibs"|"\\$AR \\$AR_FLAGS \\$oldlib\\$oldobjs\\$old_deplibs~\\$RANLIB \\$oldlib"|' \
          -e 's|libext="lib"|libext="a"|' \
      libtool > conftest.bla

      mv conftest.bla libtool
      chmod 755 libtool
      ;;
  esac
esac

# ToDo
# For now, don't use the -no-undefined flag, since the Makefiles are
# not yet set up that way.  But we need to fix this, when we want
# to comile DLLs under Windows.
LT_LDFLAGS=
AC_SUBST(LT_LDFLAGS)
]) # AC_COIN_PROG_LIBTOOL

# This is a trick to force the check for the dlfcn header to be done before
# the checks for libtool
AC_DEFUN([AC_COIN_DLFCN_H],
[AC_LANG_PUSH(C)
AC_COIN_CHECK_HEADER([dlfcn.h])
AC_LANG_POP(C)
]) # AC_COIN_DLFCN_H

###########################################################################
#                            COIN_RPATH_FLAGS                             #
###########################################################################

# This macro, in case shared objects are used, defines a variable
# RPATH_FLAGS that can be used by the linker to hardwire the library
# search path for the given directories.  This is useful for example
# Makefiles

AC_DEFUN([AC_COIN_RPATH_FLAGS],
[RPATH_FLAGS=

if test $enable_shared = yes; then
  case $build in
    *-linux-*)
      if test "$GXX" = "yes"; then
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -Wl,--rpath -Wl,$dir"
        done
      fi ;;
    *-darwin*)
        RPATH_FLAGS=nothing ;;
    *-ibm-*)
      case "$CXX" in
      xlC* | */xlC* | mpxlC* | */mpxlC*)
        RPATH_FLAGS=nothing ;;
      esac ;;
    *-hp-*)
        RPATH_FLAGS=nothing ;;
    *-mingw32)
        RPATH_FLAGS=nothing ;;
    *-sun-*)
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -R$dir"
        done
  esac

  if test "$RPATH_FLAGS" = ""; then
    AC_MSG_WARN([Could not automatically determine how to tell the linker about automatic inclusion of the path for shared libraries.  The test examples might not work if you link against shared objects.  You will need to set the LD_LIBRARY_PATH, DYLP_LIBRARY_PATH, or LIBDIR variable manually.])
  fi
  if test "$RPATH_FLAGS" = "nothing"; then
    RPATH_FLAGS=
  fi
fi

AC_SUBST(RPATH_FLAGS)
]) # AC_COIN_RPATH_FLAGS

###########################################################################
#                              COIN_FINALIZE                              #
###########################################################################

# This macro should be called at the very end of the configure.ac file.
# It creates the output files (by using AC_OUTPUT), and might do some other
# things (such as generating links to data files in a VPATH configuration).
# It also prints the "success" message.

AC_DEFUN([AC_COIN_FINALIZE],
[
FADDLIBS="$ADDLIBS"
if test x"$coin_need_flibs" = xyes; then
  ADDLIBS="$ADDLIBS $FLIBS"
fi

# library extension
AC_SUBST(LIBEXT)
case "$CC" in
  cl*) LIBEXT=lib ;;
    *) LIBEXT=a ;;
esac

# Define VPATH_DISTCLEANFILES to be everything that needs to be
# cleaned for distclean in a vpath configuration
AC_SUBST(VPATH_DISTCLEANFILES)
VPATH_DISTCLEANFILES="$coin_vpath_link_files"

AC_OUTPUT

if test x"$coin_vpath_link_files" = x; then : ; else
  lnkcmd=
  if test "$enable_doscompile" = yes; then
    lnkcmd=cp
  fi
  case "$CC" in
    cl* | */cl*)
      lnkcmd=cp ;;
  esac
  if test "$lnkcmd" = cp; then
    AC_MSG_NOTICE(Copying data files for VPATH configuration)
  else
    AC_PROG_LN_S
    AC_MSG_NOTICE(Creating VPATH links for data files)
    lnkcmd="$LN_S"
  fi
  for file in $coin_vpath_link_files; do
    dir=`AS_DIRNAME(["./$file"])`
    if test -d $dir; then : ; else
      AS_MKDIR_P($dir)
    fi
    rm -f $file
    $lnkcmd $abs_source_dir/$file $file
  done
fi

if test x$coin_projectdir = xyes; then
  AC_MSG_NOTICE([Configuration of $PACKAGE_NAME successful])
else
  AC_MSG_NOTICE([Main configuration of $PACKAGE_NAME successful])
fi
]) #AC_COIN_FINALIZE

###########################################################################
#                             COIN_VPATH_LINK                             #
###########################################################################

# This macro makes sure that a symbolic link is created to a file in
# the source code directory tree if we are in a VPATH compilation, and
# if this package is the main package to be installed

AC_DEFUN([AC_COIN_VPATH_LINK],
[AC_REQUIRE([AC_COIN_CHECK_VPATH])
if test $coin_vpath_config = yes; then
  coin_vpath_link_files="$coin_vpath_link_files $1"
fi
]) #AC_COIN_VPATH_LINK

###########################################################################
#                       COIN_ENABLE_GNU_PACKAGES                          #
###########################################################################

# This macro defined the --enable-gnu-packages flag.  This can be used
# to check if a user wants to compile GNU packges (such as readline or
# zlib) into the executable.  By default, GNU packages are disabled.
# This also defines the automake conditional COIN_ENABLE_GNU_PACKAGES

AC_DEFUN([AC_COIN_ENABLE_GNU_PACKAGES],
[AC_ARG_ENABLE([gnu-packages],
               [AC_HELP_STRING([--enable-gnu-packages],
                               [compile with GNU packages (disabled by default)])],
	       [coin_enable_gnu=$enableval],
	       [coin_enable_gnu=no])
]) # AC_COIN_ENABLE_GNU_PACKAGES

###########################################################################
#                           COIN_CHECK_GNU_ZLIB                           #
###########################################################################

# This macro checks for the libz library.

AC_DEFUN([AC_COIN_CHECK_GNU_ZLIB],
[AC_REQUIRE([AC_COIN_ENABLE_GNU_PACKAGES])
AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_BEFORE([$0],[AC_COIN_FINISH])

coin_has_zlib=no
if test $coin_enable_gnu = yes; then
  AC_COIN_CHECK_HEADER([zlib.h],[coin_has_zlib=yes])

  if test $coin_has_zlib = yes; then
    AC_CHECK_LIB([z],[gzopen],
                 [ADDLIBS="-lz $ADDLIBS"],
                 [coin_has_zlib=no])
  fi

  if test $coin_has_zlib = yes; then
    AC_DEFINE([COIN_HAS_ZLIB],[1],[Define to 1 if zlib is available])
  fi
fi

AM_CONDITIONAL(COIN_HAS_ZLIB,test x$coin_has_zlib = xyes)
]) # AC_COIN_CHECK_GNU_ZLIB


###########################################################################
#                          COIN_CHECK_GNU_BZLIB                           #
###########################################################################

# This macro checks for the libbz2 library.

AC_DEFUN([AC_COIN_CHECK_GNU_BZLIB],
[AC_REQUIRE([AC_COIN_ENABLE_GNU_PACKAGES])
AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_BEFORE([$0],[AC_COIN_FINISH])

coin_has_bzlib=no
if test $coin_enable_gnu = yes; then
  AC_COIN_CHECK_HEADER([bzlib.h],[coin_has_bzlib=yes])

  if test $coin_has_bzlib = yes; then
    AC_CHECK_LIB([bz2],[BZ2_bzReadOpen],
                 [ADDLIBS="-lbz2 $ADDLIBS"],
                 [coin_has_bzlib=no])
  fi

  if test $coin_has_bzlib = yes; then
    AC_DEFINE([COIN_HAS_BZLIB],[1],[Define to 1 if bzlib is available])
  fi
fi
]) # AC_COIN_CHECK_GNU_BZLIB


###########################################################################
#                         COIN_CHECK_GNU_READLINE                         #
###########################################################################

# This macro checks for GNU's readline.  It verifies that the header
# readline/readline.h is available, and that the -lreadline library
# contains "readline".  It is assumed that #include <stdio.h> is included
# in the source file before the #include<readline/readline.h>

AC_DEFUN([AC_COIN_CHECK_GNU_READLINE],
[AC_REQUIRE([AC_COIN_ENABLE_GNU_PACKAGES])
AC_BEFORE([AC_COIN_PROG_CXX],[$0])
AC_BEFORE([AC_COIN_PROG_CC],[$0])
AC_BEFORE([AC_COIN_PROG_F77],[$0])
AC_BEFORE([$0],[AC_COIN_FINISH])

coin_has_readline=no
if test $coin_enable_gnu = yes; then
  AC_COIN_CHECK_HEADER([readline/readline.h],
                       [coin_has_readline=yes],[],
                       [#include <stdio.h>])

  coin_save_LIBS="$LIBS"
  LIBS=
  # First we check if tputs and friends are available
  if test $coin_has_readline = yes; then
    AC_SEARCH_LIBS([tputs],[ncurses termcap curses],[],
                   [coin_has_readline=no])
  fi

  # Now we check for readline
  if test $coin_has_readline = yes; then
    AC_CHECK_LIB([readline],[readline],
                 [ADDLIBS="-lreadline $LIBS $ADDLIBS"],
                 [coin_has_readline=no])
  fi

  if test $coin_has_readline = yes; then
    AC_DEFINE([COIN_HAS_READLINE],[1],[Define to 1 if readline is available])
  fi

  LIBS="$coin_save_LIBS"
fi
]) # AC_COIN_CHECK_GNU_READLINE

###########################################################################
#                             COIN_DATA_PATH                              #
###########################################################################

# This macro defines a preprocessor macro with the absolute path to a
# subdirectory of Data.  The argument of this macro is the name of the
# subdirectory (in correct case), and the name of the macro is
# COIN_DATA_DIR_PATH, where dir is replaced by the capitalized name of
# the directory.  The path ends with a separator ("/" for linux and
# '\\' for Windows).  The default value for this path can be
# overwritten with the input variable with the same name
# (COIN_DATA_DIR_PATH).  At this point we chech only for the
# $srcdir/../Data subdirectory.

AC_DEFUN([AC_COIN_DATA_PATH],
[AC_MSG_CHECKING([absolute path to data directory $1])

AC_ARG_VAR(m4_toupper(COIN_DATA_$1_PATH),[Set to absolute path to Data/$1 subdirectory])

if test x"$m4_toupper(COIN_DATA_$1_PATH)" = x; then
  m4_toupper(COIN_DATA_$1_PATH)=`cd $srcdir/../Data/$1; pwd`
fi

# Under Cygwin, use Windows path.  Add separator
case $build in
  *-cygwin*)
    m4_toupper(COIN_DATA_$1_PATH)=`cygwin -w $m4_toupper(COIN_DATA_$1_PATH)`\\
    ;;
  *)
    m4_toupper(COIN_DATA_$1_PATH)="$m4_toupper(COIN_DATA_$1_PATH)/"
    ;;
esac

if test -r $m4_toupper(COIN_DATA_$1_PATH); then
  AC_DEFINE_UNQUOTED(m4_toupper(COIN_DATA_$1_PATH),["$m4_toupper(COIN_DATA_$1_PATH)"],
            [Define to absolute path for Data subdirectory $1])
  AC_MSG_RESULT($m4_toupper(COIN_DATA_$1_PATH))
else
  AC_MSG_ERROR(Directory $m4_toupper(COIN_DATA_$1_PATH) does not exist)
fi
]) # AC_COIN_HAS_DATA

###########################################################################
#                          COIN_EXAMPLE_FILES                             #
###########################################################################

# This macro determines the names of the example files (using the
# argument in an "ls" command) and sets up the variables EXAMPLE_FILES
# and EXAMPLE_CLEAN_FILES.  If this is a VPATH configuration, it also
# creates soft links to the example files.

AC_DEFUN([AC_COIN_EXAMPLE_FILES],
[AC_REQUIRE([AC_COIN_CHECK_GNU_ZLIB])
AC_REQUIRE([AC_COIN_CHECK_VPATH])
files=`cd $srcdir; ls $1`
# We need to do the following loop to make sure that are no newlines
# in the variable
for file in $files; do
  EXAMPLE_FILES="$EXAMPLE_FILES $file"
done
if test $coin_vpath_config = yes; then
  lnkcmd=
  if test "$enable_doscompile" = yes; then
    lnkcmd=cp
  fi
  case "$CC" in
    cl* | */cl*)
      lnkcmd=cp ;;
  esac
  if test "$lnkcmd" = cp; then
    AC_MSG_NOTICE([Copying example files ($1)])
  else
    AC_PROG_LN_S
    AC_MSG_NOTICE([Creating links to the example files ($1)])
    lnkcmd="$LN_S"
  fi
  for file in $EXAMPLE_FILES; do
    rm -f $file
    $lnkcmd $srcdir/$file $file
  done
  EXAMPLE_CLEAN_FILES='$1'
else
  EXAMPLE_CLEAN_FILES=
fi

# In case there are compressed files, we create a variable with the
# uncompressed names
EXAMPLE_UNCOMPRESSED_FILES=
for file in $EXAMPLE_FILES; do
  case $file in
    *.gz)
      EXAMPLE_UNCOMPRESSED_FILES="$EXAMPLE_UNCOMPRESSED_FILES `echo $file | sed -e s/.gz//`"
      ;;
  esac
done

AC_SUBST(EXAMPLE_UNCOMPRESSED_FILES)
AC_SUBST(EXAMPLE_FILES)
AC_SUBST(EXAMPLE_CLEAN_FILES)
]) #AC_COIN_EXAMPLE_FILES

###########################################################################
#                            COIN_HAS_PROJECT                             #
###########################################################################

# This macro sets up usage of a Coin package.  It defines the
# PKGSRCDIR and PKGOBJDIR variables, refering to the main source and
# object directory of the package, respectively.  It also defines
# a COIN_HAS_PKG preprocessor macro and makefile conditional.  The
# argument should be the name (Pkg) of the project (in correct lower
# and upper case)

AC_DEFUN([AC_COIN_HAS_PROJECT],
[AC_MSG_CHECKING([for COIN project $1])

# First check, if the sub-project is actually available (ToDo: allow
# other locations)

m4_tolower(coin_has_$1)=unavailable
if test x"$COIN_SKIP_PROJECTS" != x; then
  for dir in $COIN_SKIP_PROJECTS; do
    if test $dir = $1; then
      m4_tolower(coin_has_$1)=skipping
    fi
  done
fi

if test $m4_tolower(coin_has_$1) != skipping; then
  if test $PACKAGE_TARNAME = m4_tolower($1); then
    m4_tolower(coin_has_$1)=.
  else
    if test -d $srcdir/../$1; then
      m4_tolower(coin_has_$1)=../$1
    fi
  fi
fi

if test $m4_tolower(coin_has_$1) != unavailable &&
   test $m4_tolower(coin_has_$1) != skipping; then
  # Set the #define if the component is available
  AC_DEFINE(m4_toupper(COIN_HAS_$1),[1],[Define to 1 if the $1 package is used])

  # Set the variables for source and object code location
  AC_SUBST(m4_toupper($1SRCDIR))
  m4_toupper($1SRCDIR)=`cd $srcdir/$m4_tolower(coin_has_$1); pwd`
  AC_SUBST(m4_toupper($1OBJDIR))
  m4_toupper($1OBJDIR)=`pwd`/$m4_tolower(coin_has_$1)
fi

  # Define the Makefile conditional
AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),
               [test $m4_tolower(coin_has_$1) != unavailable &&
                test $m4_tolower(coin_has_$1) != skipping])
AC_MSG_RESULT([$m4_tolower(coin_has_$1)])
]) # AC_COIN_HAS

###########################################################################
#                        COIN_HAS_USER_LIBRARY                            #
###########################################################################

# This macro sets up usage of a library with header files.  It defines
# the LBRYINCDIR variable, and it defines COIN_HAS_LBRY preprocessor
# macro and makefile conditional.  The first argument should be the
# full name (LibraryName) of the library, and the second argument (in
# upper case letters) the abbreviation (LBRY).  This macro also
# introduces the configure arguments --with-libraryname-incdir and
# --with-libraryname-lib which have to be both given by a user to use
# this solver to tell the configure script where the include files and
# the library are located.  Those arguments can also be given as
# environement variables LBRYINCDIR and LBRYLIB, but a --with-*
# argument overwrites an environment variable.  If a third argument is
# given, it is assumed that this is the name of a header file that can
# be checked for in the given include directory, and if a fourth
# argument is given, it is assumed to be the name of a C function
# which is given and defined in the library, and a test is done to
# check if that symbol is defined in the library.
# If it possible to disable the check, by specifying
# --disable-libraryname-libcheck - this is a workaround for platforms
# where checks don't work (yet) properly.

AC_DEFUN([AC_COIN_HAS_USER_LIBRARY],
[AC_REQUIRE([AC_COIN_PROJECTDIR_INIT])
AC_MSG_CHECKING(if user provides library for $1)

# Check for header file directory
AC_ARG_WITH(m4_tolower($1)-incdir,
            AC_HELP_STRING([--with-m4_tolower($1)-incdir],
                           [specify the directory with the header files for library $1]),
                           [$2INCDIR=`cd $withval; pwd`])
# Check for library directory
AC_ARG_WITH(m4_tolower($1)-lib,
            AC_HELP_STRING([--with-m4_tolower($1)-lib],
                           [specify the flags to link with the library $1]),
                           [$2LIB=$withval])
# Switch to disable library check if requested
AC_ARG_ENABLE(m4_tolower($1)-libcheck,
              AC_HELP_STRING([--enable-m4_tolower($1)-libcheck],
                             [use disable-m4_tolower($1)-libcheck to skip the link check at configuration time]),
              [m4_tolower($1)_libcheck=$enableval],
              [m4_tolower($1)_libcheck=yes])

if test x"$$2INCDIR" != x || test x"$$2LIB" != x; then
  m4_tolower(coin_has_$2)=true
else
  m4_tolower(coin_has_$2)=false
fi

if test $m4_tolower(coin_has_$2) = true; then
# Check either both arguments or none are given
  if test x"$$2INCDIR" = x || test x"$$2LIB" = x; then
    AC_MSG_ERROR([You need to specify both --with-m4_tolower($1)-incdir and --with-m4_tolower($1)-lib if you want to use library $1])
  fi
  AC_MSG_RESULT(yes)
  # Check if the given header file is there
  m4_ifvaln([$3],[AC_CHECK_FILE([$$2INCDIR/$3],[],
                 [AC_MSG_ERROR([Cannot find file $3 in $$2INCDIR])])])
  # Check if the symbol is provided in the library
  # ToDo: FOR NOW WE ASSUME THAT WE ARE USING THE C++ COMPILER
  m4_ifvaln([$4],[if test x"$m4_tolower($1)_libcheck" != xno; then
                    coin_save_LIBS="$LIBS"
                    LIBS="$$2LIB $ADDLIBS"
		    AC_MSG_CHECKING([whether symbol $4 is available with $2])
# ToDo find out what to do about extern "C"
#                    AC_TRY_LINK([extern "C" {void $4();}],[$4()],
                    AC_TRY_LINK([void $4();],[$4()],
                                [AC_MSG_RESULT(yes)],
			        [AC_MSG_RESULT(no)
                                 AC_MSG_ERROR([Cannot find symbol $4 with $2])])
                    LIBS="$coin_save_LIBS"
                  fi])
  ADDLIBS="$$2LIB $ADDLIBS"
  AC_DEFINE(COIN_HAS_$2,[1],[Define to 1 if the $1 package is used])
else
  AC_MSG_RESULT(no)
fi

AC_SUBST($2INCDIR)
AC_SUBST($2LIB)
AM_CONDITIONAL(COIN_HAS_$2,
               test $m4_tolower(coin_has_$2) = true)
]) #AC_COIN_HAS_SOLVER 

###########################################################################
#                               COIN_HAS_ASL                              #
###########################################################################

# This macro checks if the user has provide arguments that say where
# the precompiled ASL files should be found (with the --with-asldir
# flag).  If this is not the case, we check if the ThirdParty/ASL
# directory has been configured, which indicates that the files will
# be in that directory and can be used.

AC_DEFUN([AC_COIN_HAS_ASL],
[coin_aslobjdir=../ThirdParty/ASL
coin_aslsrcdir=$srcdir/$coin_aslobjdir

# Determine the name of the ASL library
case "$CXX" in
  cl* | */cl*)
    ampllib=amplsolv.lib ;;
  *)
    ampllib=amplsolver.a ;;
esac

AC_ARG_WITH([asldir],
            AC_HELP_STRING([--with-asldir],
                           [specify path to AMPL solver directory (or BUILD for compilation, or "no" for disabling AMPL)]),
            [use_asldir=$withval], [use_asldir=])

if test "$use_asldir" = BUILD; then
  AC_CHECK_FILE([$coin_aslobjdir/Makefile],[],
                [AC_MSG_ERROR([option \"BUILD\" specified for asldir, but directory is not configure (sources missing?)])])
elif test -z "$use_asldir"; then
 # try to find sources - if not given don't compile
  AC_CHECK_FILE([$coin_aslobjdir/Makefile],[use_asldir=BUILD],[use_asldir=no])
elif test "$use_asldir" != "no"; then
  AC_CHECK_FILE([$use_asldir/$ampllib],[],
                [AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but library missing])])
  AC_CHECK_FILE([$use_asldir/asl.h],[],
                [AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but header files are missing])])
  use_asldir=`cd $use_asldir; pwd`
  case $build in
    *-cygwin*) use_asldir=`cygpath -w $use_asldir | sed -e sX\\\\\\\\X/Xg` ;;
  esac
fi

# Variable containing ASL library (including full path)
AC_SUBST(ASLLIB)
# Variable containing flags for including ASL header files
AC_SUBST(ASL_CPPFLAGS)

if test "$use_asldir" = BUILD; then
  coin_aslobjdir=`cd $coin_aslobjdir; pwd`
  ASLLIB=`$CYGPATH_W $coin_aslobjdir/$ampllib | sed -e sX\\\\\\\\X/Xg`
  coin_aslsrcdir=`cd $coin_aslsrcdir; pwd`
  ASL_CPPFLAGS="-I"`$CYGPATH_W $coin_aslobjdir | sed -e sX\\\\\\\\X/Xg`" -I"`$CYGPATH_W $coin_aslsrcdir/solvers | sed -e sX\\\\\\\\X/Xg`
elif test "$use_asldir" != no; then
  ASLLIB=`$CYGPATH_W $use_asldir/$ampllib | sed -e sX\\\\\\\\X/Xg`
  ASL_CPPFLAGS="-I"`$CYGPATH_W $use_asldir | sed -e sX\\\\\\\\X/Xg`
fi

if test "$use_asldir" != no; then
  AC_CHECK_LIB(dl,[dlopen],[ASLLIB="$ASLLIB -ldl"],[])
  coin_has_asl=yes
  AC_DEFINE([COIN_HAS_ASL],[1],
            [If defined, the Ampl Solver Library is available.])
else
  coin_has_asl=no
fi
AM_CONDITIONAL(COIN_HAS_ASL, test $coin_has_asl = yes)
]) # AC_COIN_HAS_ASL

###########################################################################
#                            COIN_TRY_FLINK                               #
###########################################################################

# Auxilliary macro to test if a Fortran function name can be linked,
# given the current settings of LIBS.  We determine from the context, what
# the currently active programming language is, and cast the name accordingly.
# The first argument is the name of the function/subroutine, in small letters,
# the second argument are the actions taken when the test works, and the
# third argument are the actions taken if the test fails.

AC_DEFUN([AC_COIN_TRY_FLINK],
[case $ac_ext in
  f)
    AC_TRY_LINK([],[      call $1],[$2],[$3])
    ;;
  c)
    AC_F77_FUNC($1,cfunc$1)
    if test x"$coin_need_flibs" = xyes; then
      flink_try=no;
    else
      AC_TRY_LINK([void $cfunc$1();],[$cfunc$1()],
                  [flink_try=yes],[flink_try=no])
    fi
    if test $flink_try = yes; then
      $2
    else
      if test x"$FLIBS" != x; then
        flink_save_libs="$LIBS"
        LIBS="$LIBS $FLIBS"
        AC_TRY_LINK([void $cfunc$1();],[$cfunc$1()],
                    [LIBS="$flink_save_libs"
                     $2
                     coin_need_flibs=yes],
                    [LIBS="$flink_save_libs"
                     $3])
      else
        $3
      fi
    fi
    ;;
  cc)
    AC_F77_FUNC($1,cfunc$1)
    if test x"$coin_need_flibs" = xyes; then
      flink_try=no;
    else
      AC_TRY_LINK([extern "C" {void $cfunc$1();}],[$cfunc$1()],
                  [flink_try=yes],[flink_try=no])
    fi
    if test $flink_try = yes; then
      $2
    else
      if test x"$FLIBS" != x; then
        flink_save_libs="$LIBS"
        LIBS="$LIBS $FLIBS"
        AC_TRY_LINK([extern "C" {void $cfunc$1();}],[$cfunc$1()],
                    [LIBS="$flink_save_libs"
                     $2
                     coin_need_flibs=yes],
                    [LIBS="$flink_save_libs"
                     $3])
      else
        $3
      fi
    fi
    ;;
esac
]) # AC_COIN_TRY_FLINK

###########################################################################
#                             COIN_HAS_BLAS                               #
###########################################################################

# This macro checks for a library containing the BLAS library.  It
# tried standard libraries, and if none is found to be working, it
# checks whether the BLAS ThirdParty/Blas directory has been configured.
# It adds to ADDLIBS any flags required to link with an externally provided
# BLAS.  It defines the makefile conditional and preprocessor macro
# COIN_HAS_BLAS, if blas is available, and it defines the makefile conditional
# COIN_BUILD_BLAS, if blas is compiled within COIN.

AC_DEFUN([AC_COIN_HAS_BLAS],
[coin_blasobjdir=../ThirdParty/Blas
coin_blassrcdir=$srcdir/$coin_blasobjdir

AC_ARG_WITH([blas],
            AC_HELP_STRING([--with-blas],
                           [specify BLAS library (or BUILD for compilation)]),
            [use_blas=$withval], [use_blas=])

# Check if user supplied option makes sense
if test x"$use_blas" != x; then
  if test "$use_blas" = "BUILD"; then
    AC_CHECK_FILE([$coin_blasobjdir/Makefile],[],
                  [AC_MSG_ERROR([option \"BUILD\" specified for Blas, but $coin_blasobjdir directory is not configured])])
  else
    AC_MSG_CHECKING([whether user supplied BLASLIB=\"$use_blas\" works])
    LIBS="$use_blas $LIBS"
    ADDLIBS="$use_blas $ADDLIBS"
    AC_COIN_TRY_FLINK([daxpy],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])
                       AC_MSG_ERROR([user supplied BLAS library \"$use_blas\" does not work])])
  fi
else
# Try to autodetect the library for blas based on build system
  case $build in
    *-sgi-*) 
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([whether -lcomplib.sgimath has BLAS])
      LIBS="-lcomplib.sgimath $LIBS"
      AC_COIN_TRY_FLINK([daxpy],
                        [AC_MSG_RESULT([yes])
                         use_blas=-lcomplib.sgimath;
                         ADDLIBS="-lcomplib.sgimath $ADDLIBS"],
                        [AC_MSG_RESULT([no])
                         SAVE_LIBS="$LIBS"])
      ;;
    *-sun-*)
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([whether -xlic_lib=sunperf has BLAS])
      LIBS="-xlic_lib=sunperf $LIBS"
      AC_COIN_TRY_FLINK([daxpy],
                        [AC_MSG_RESULT([yes])
                         use_blas='-xlic_lib=sunperf'
                         ADDLIBS="-xlic_lib=sunperf $ADDLIBS"],
                        [AC_MSG_RESULT([no])
                         LIBS="$SAVE_LIBS"])
      ;;
  esac
  # On cygwin, if enable_doscompile is used, recompile blas because it
  # otherwise links with the cygwin blas which doesn't run under DOS
  if test "$enable_doscompile" != yes; then
    if test -z "$use_blas"; then
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([whether -lblas has BLAS])
      LIBS="-lblas $LIBS"
      AC_COIN_TRY_FLINK([daxpy],
                        [AC_MSG_RESULT([yes])
                         ADDLIBS="-lblas $ADDLIBS"
                         use_blas='-lblas'],
                        [AC_MSG_RESULT([no])
                         LIBS="$SAVE_LIBS"])
    fi
  fi
  if test -z "$use_blas"; then
    AC_CHECK_FILE([$coin_blasobjdir/Makefile],[use_blas=BUILD])
  fi
fi

AM_CONDITIONAL([COIN_HAS_BLAS],[test x"$use_blas" != x])
AM_CONDITIONAL([COIN_BUILD_BLAS],[test "$use_blas" = BUILD])

if test x"$use_blas" = x; then
  coin_has_blas=no
else
  coin_has_blas=yes
  AC_DEFINE([COIN_HAS_BLAS],[1],
            [If defined, the BLAS Library is available.])
fi
]) # AC_COIN_HAS_BLAS

###########################################################################
#                            COIN_HAS_LAPACK                              #
###########################################################################

# This macro checks for a library containing the LAPACK library.  It
# tried standard libraries, and if none is found to be working, it
# checks whether the LAPACK ThirdParty/Lapack directory has been
# configured.  It adds to ADDLIBS any flags required to link with an
# externally provided LAPACK.  It defines the makefile conditional and
# preprocessor macro COIN_HAS_LAPACK, if lapack is available, and it
# defines the makefile conditional COIN_BUILD_LAPACK, if lapack is
# compiled within COIN.

AC_DEFUN([AC_COIN_HAS_LAPACK],
[coin_lapackobjdir=../ThirdParty/Lapack
coin_lapacksrcdir=$srcdir/$coin_lapackobjdir

AC_ARG_WITH([lapack],
            AC_HELP_STRING([--with-lapack],
                           [specify LAPACK library (or BUILD for compilation)]),
            [use_lapack=$withval], [use_lapack=])

# Check if user supplied option makes sense
if test x"$use_lapack" != x; then
  if test "$use_lapack" = "BUILD"; then
    AC_CHECK_FILE([$coin_lapackobjdir/Makefile],[],
                  [AC_MSG_ERROR([option \"BUILD\" specified for Lapack, but $coin_lapackobjdir directory is not configured])])
  else
    AC_MSG_CHECKING([whether user supplied LAPACKLIB=\"$use_lapack\" works])
    LIBS="$use_lapack $LIBS"
    ADDLIBS="$use_lapack $ADDLIBS"
    AC_COIN_TRY_FLINK([dsyev],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])
                       AC_MSG_ERROR([user supplied LAPACK library \"$use_lapack\" does not work])])
  fi
else
  if test x$coin_has_blas = xyes; then
    # First try to see if LAPACK is already available with BLAS library
    AC_MSG_CHECKING([whether LAPACK is already available with BLAS library])
    AC_COIN_TRY_FLINK([dsyev],
                      [AC_MSG_RESULT([yes]); use_lapack=ok],
                      [AC_MSG_RESULT([no])])
  fi
  if test -z "$use_lapack"; then
    # Try to autodetect the library for lapack based on build system
    case $build in
      *-sgi-*) 
        SAVE_LIBS="$LIBS"
        AC_MSG_CHECKING([whether -lcomplib.sgimath has LAPACK])
        LIBS="-lcomplib.sgimath $LIBS"
        AC_COIN_TRY_FLINK([dsyev],
                          [AC_MSG_RESULT([yes])
                           use_lapack=-lcomplib.sgimath;
                           ADDLIBS="-lcomplib.sgimath $ADDLIBS"],
                          [AC_MSG_RESULT([no])
                           SAVE_LIBS="$LIBS"])
        ;;
      *-sun-*)
        SAVE_LIBS="$LIBS"
        AC_MSG_CHECKING([whether -xlic_lib=sunperf has LAPACK])
        LIBS="-xlic_lib=sunperf $LIBS"
        AC_COIN_TRY_FLINK([dsyev],
                          [AC_MSG_RESULT([yes])
                           use_lapack='-xlic_lib=sunperf'
                           ADDLIBS="-xlic_lib=sunperf $ADDLIBS"],
                          [AC_MSG_RESULT([no])
                           LIBS="$SAVE_LIBS"])
        ;;
    esac
  fi
  # On cygwin, if enable_doscompile is used, recompile lapack because it
  # otherwise links with the cygwin lapack which doesn't run under DOS
  if test "$enable_doscompile" != yes; then
    if test -z "$use_lapack"; then
      SAVE_LIBS="$LIBS"
      AC_MSG_CHECKING([whether -llapack has LAPACK])
      LIBS="-llapack $LIBS"
      AC_COIN_TRY_FLINK([dsyev],
                        [AC_MSG_RESULT([yes])
                         ADDLIBS="-llapack $ADDLIBS"
                         use_lapack='-llapack'],
                        [AC_MSG_RESULT([no])
                         LIBS="$SAVE_LIBS"])
    fi
  fi
  if test -z "$use_lapack"; then
    AC_CHECK_FILE([$coin_lapackobjdir/Makefile],[use_lapack=BUILD])
  fi
fi

AM_CONDITIONAL([COIN_HAS_LAPACK],[test x"$use_lapack" != x])
AM_CONDITIONAL([COIN_BUILD_LAPACK],[test "$use_lapack" = BUILD])

if test x"$use_lapack" = x; then
  coin_has_lapack=no
else
  coin_has_lapack=yes
  AC_DEFINE([COIN_HAS_LAPACK],[1],
            [If defined, the LAPACK Library is available.])
fi
]) # AC_COIN_HAS_LAPACK

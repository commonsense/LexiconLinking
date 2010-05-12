
# This macro looks for flags in CXXFLAGS and/or CFLAGS that really should be
# in CPPFLAGS because they affect the preprocessor's search path.
# See Section 4.7.1 in the autoconf manual.

# Current tests are:
#   -mno-cygwin	Sort of a special case, needs to be in preprocess, compile,
#		and link.

AC_DEFUN([AC_DYLP_FIX_CPPFLAGS],
[ case "$CXXFLAGS $CFLAGS" in
    *-mno-cygwin*)
      CXXFLAGS=`echo $CXXFLAGS | sed -e 's/-mno-cygwin//'`
      CFLAGS=`echo $CFLAGS | sed -e 's/-mno-cygwin//'`
      CPPFLAGS="-mno-cygwin $CPPFLAGS"
      LDFLAGS="-mno-cygwin $LDFLAGS"
      ;;
  esac
])

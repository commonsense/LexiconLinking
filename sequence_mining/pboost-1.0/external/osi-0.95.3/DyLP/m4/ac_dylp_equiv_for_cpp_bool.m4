
# AC_DYLP_EQUIV_FOR_CPP_BOOL
# ------------------------------------------------------------------------
# Calls to _AC_COMPUTE_INT lifted from the definition of AC_CHECK_SIZEOF in
# autoconf/types.m4. Use of _AC_COMPUTE_INT automatically deals with cross
# compilation (no mean feat; check the macro in autoconf/general.m4!). The
# use of (long) (sizeof (bool)) as the expression is a workaround for some
# HP compiler bug. See full comment in AC_CHECK_SIZEOF.
# ------------------------------------------------------------------------
AC_DEFUN([AC_DYLP_EQUIV_FOR_CPP_BOOL],
[ AC_MSG_NOTICE([Determining C type equivalent for C++ bool.])
  AC_LANG_PUSH(C++)
  _AC_COMPUTE_INT([(long) (sizeof (bool))],
      [ac_cv_sizeof_cpp_bool],
      [AC_INCLUDES_DEFAULT([$3])],
      [AC_MSG_FAILURE([cannot compute sizeof (bool) for C++])])
  AC_LANG_POP(C++)
# Force a particular value to test the code below.
# ac_cv_sizeof_cpp_bool=8
  AC_MSG_NOTICE([C++ bool is $ac_cv_sizeof_cpp_bool bytes.])

  AC_LANG_PUSH(C)
  dylp_booltype="no"
  _AC_COMPUTE_INT([(long) (sizeof (char))],
      [ac_cv_sizeof_c_bool],
      [AC_INCLUDES_DEFAULT([$3])],
      [AC_MSG_FAILURE([cannot compute sizeof (char) for C])])
  if test $ac_cv_sizeof_cpp_bool = $ac_cv_sizeof_c_bool; then
    dylp_booltype="char"
  fi
  if test $dylp_booltype = "no"; then
    _AC_COMPUTE_INT([(long) (sizeof (int))],
	[ac_cv_sizeof_c_bool],
	[AC_INCLUDES_DEFAULT([$3])],
	[AC_MSG_FAILURE([cannot compute sizeof (int) for C])])
    if test $ac_cv_sizeof_cpp_bool = $ac_cv_sizeof_c_bool; then
      dylp_booltype="int"
    fi
  fi
  if test $dylp_booltype = "no"; then
    _AC_COMPUTE_INT([(long) (sizeof (short int))],
	[ac_cv_sizeof_c_bool],
	[AC_INCLUDES_DEFAULT([$3])],
	[AC_MSG_FAILURE([cannot compute sizeof (short int) for C])])
    if test $ac_cv_sizeof_cpp_bool = $ac_cv_sizeof_c_bool; then
      dylp_booltype="short int"
    fi
  fi
  if test $dylp_booltype = "no"; then
    _AC_COMPUTE_INT([(long) (sizeof (long int))],
	[ac_cv_sizeof_c_bool],
	[AC_INCLUDES_DEFAULT([$3])],
	[AC_MSG_FAILURE([cannot compute sizeof (long int) for C])])
    if test $ac_cv_sizeof_cpp_bool = $ac_cv_sizeof_c_bool; then
      dylp_booltype="long int"
    fi
  fi
  if test $dylp_booltype = "no"; then
    dylp_booltype="char"
    AC_MSG_WARN([Cannot determine C type to match C++ bool. Defaulting to
		 char. Dylp will compile, but will certainly crash when
		 executed.])
  fi
  AC_DEFINE_UNQUOTED([BOOL],[$dylp_booltype],
      [Define to the C type whose size in bytes matches the size in bytes
       of the the C++ bool type.])
  AC_LANG_POP(C)

  AC_MSG_NOTICE([C $dylp_booltype will be used as bool by dylp.])
])


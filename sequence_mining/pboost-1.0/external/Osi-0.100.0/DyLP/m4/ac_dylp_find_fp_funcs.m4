
# AC_DYLP_FIND_ISFINITE
# ------------------------------------------------------
# Determines the name of the finite() function in this environment. This is the
# function that recognises whether an IEEE floating point value is finite.  The
# variable ac_name_of_finite will be set to the proper name on return and
# DYLP_ISFINITE will be defined to the same value.
# ------------------------------------------------------
AC_DEFUN([AC_DYLP_FIND_ISFINITE],
[ AC_MSG_NOTICE([Checking for proper name for isfinite().])
  ac_name_of_isfinite="unavailable"
  AC_CHECK_FUNC([finite],[ac_name_of_isfinite=finite])
  if test "$ac_name_of_isfinite" = "unavailable"; then
    AC_CHECK_FUNC([_finite],[ac_name_of_isfinite=_finite])
  fi
  if test "$ac_name_of_isfinite" = "unavailable"; then
    AC_CHECK_FUNC([isfinite],[ac_name_of_isfinite=isfinite])
  fi
  if test "$ac_name_of_isfinite" = "unavailable"; then
    AC_MSG_WARN([Cannot find a C function to check if an IEEE floating point
		 value is finite. There is no hope of building dylp on this
		 system.])
  fi
  AC_DEFINE_UNQUOTED([DYLP_ISFINITE],
      [$ac_name_of_isfinite],
      [Define to be the name of the C function used to check that an IEEE
       floating point value is finite.])
  AC_MSG_NOTICE([Using $ac_name_of_isfinite as isfinite().])
])

# AC_DYLP_FIND_ISNAN
# ------------------------------------------------------
# Determines the name of the isnan() function in this environment. This is the
# function that recognises whether an IEEE floating point value is NaN.  The
# variable ac_name_of_isnan will be set to the proper name on return and
# DYLP_ISNAN will be defined to the same value.
# ------------------------------------------------------
AC_DEFUN([AC_DYLP_FIND_ISNAN],
[ AC_MSG_NOTICE([Checking for proper name for isnan().])
  ac_name_of_isnan="unavailable"
  AC_CHECK_FUNC([isnan],[ac_name_of_isnan=isnan])
  if test "$ac_name_of_isnan" = "unavailable"; then
    AC_CHECK_FUNC([_isnan],[ac_name_of_isnan=_isnan])
  fi
  if test "$ac_name_of_isnan" = "unavailable"; then
    AC_MSG_WARN([Cannot find a C function to check if an IEEE floating point
		 value is NaN. There is no hope of building dylp on this
		 system.])
  fi
  AC_DEFINE_UNQUOTED([DYLP_ISNAN],
      [$ac_name_of_isnan],
      [Define to be the name of the C function used to check that an IEEE
       floating point value is NaN.])
  AC_MSG_NOTICE([Using $ac_name_of_isnan as isnan().])
])


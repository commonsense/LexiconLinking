
# AC_DYLP_PARANOIA(dflt)
# ------------------------------------------------------------------
# Processes the paranoia option.
# ------------------------------------------------------------------
AC_DEFUN([AC_DYLP_PARANOIA],
[ AC_ARG_ENABLE([dylp-paranoia],
      AS_HELP_STRING([--enable-dylp-paranoia],
          [Enable dylp's paranoid checks (default=$1)]),
      [dylp_paranoia=$enableval],
      [dylp_paranoia=$1])
  if test "$dylp_paranoia" = "yes"; then
    AC_DEFINE([DYLP_PARANOIA],[1],
	[Define this variable to enable dylp's paranoid checks.])
    AC_MSG_NOTICE([Dylp paranoid checks enabled.])
  else
    AC_MSG_NOTICE([Dylp paranoid checks disabled.])
  fi
])

# AC_DYLP_STATISTICS(dflt)
# ------------------------------------------------------------------
# Processes the statistics option.
# ------------------------------------------------------------------
AC_DEFUN([AC_DYLP_STATISTICS],
[ AC_ARG_ENABLE([dylp-stats],
      AS_HELP_STRING([--enable-dylp-stats],
          [Enable dylp's statistics collection features (default=$1)]),
      [dylp_stats=$enableval],
      [dylp_stats=$1])
  if test "$dylp_stats" = "yes"; then
    AC_DEFINE([DYLP_STATISTICS],[1],
	[Define this variable to enable dylp's statistics collection
	 features.])
    AC_MSG_NOTICE([Dylp statistics collection enabled.])
  else
    AC_MSG_NOTICE([Dylp statistics collection disabled.])
  fi
])

# AC_DYLP_INFO(dflt)
# ------------------------------------------------------------------
# Processes the information printing (info) option.
# ------------------------------------------------------------------
AC_DEFUN([AC_DYLP_INFO],
[ AC_ARG_ENABLE([dylp-info],
      AS_HELP_STRING([--enable-dylp-info],
          [Enable dylp's informational printing features (default=$1)]),
      [dylp_info=$enableval],
      [dylp_info=$1])
  if test "$dylp_info" = "no"; then
    AC_DEFINE([DYLP_NDEBUG],[1],
	[Define this variable to disable dylp's informational printing
	 features.])
    AC_MSG_NOTICE([Dylp informational printing disabled.])
  else
    AC_MSG_NOTICE([Dylp informational printing enabled.])
  fi
])


# AC_C_ADD_TO_INCLUDES (header-file)
# ------------------------------------------------------
# Augments the list of default include files with the specified header.
# ------------------------------------------------------
AC_DEFUN([AC_C_ADD_TO_INCLUDES],
[ ac_includes_default="\
$ac_includes_default
@%:@include <$1>"
])

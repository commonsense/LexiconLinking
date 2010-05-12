/*
 * Include file for the configuration of Vol.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN packages or third party code) are set
 * here.  The project maintainer needs to remember to update this file
 * and choose reasonable defines.  A user can modify the default
 * setting by editing this file here.
 *
 */

#ifndef __DYLPCONFIG_H__

#ifdef HAVE_CONFIG_H
#include "config_dylp.h"

/* undefine macros that could conflict with those in other config.h
   files */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#else /* HAVE_CONFIG_H */

/* include the COIN-wide system specific configure header */
#include "configall_system.h"

/***************************************************************************/
/*             HERE DEFINE THE CONFIGURATION SPECIFIC MACROS               */
/***************************************************************************/

/* Define to the C type corresponding to the C++ bool type */
#define BOOL char

/* Define to the debug sanity check level (0 is no test) */
#define COIN_DYLP_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_DYLP_VERBOSITY 0

/* Define to 1 if the Osi package is used */
/* #define COIN_HAS_OSI 1 */

/* Define to 1 if the DyLP package is used */
#define COIN_HAS_DYLP 1

/* Set to the full path directory name for the location of the error text
   message file dy_errmsgs.txt */
#define DYLP_ERRMSGDIR "/usr/local/share"

/* Define to 1 if DyLP's fpchecks outputs -D_GNU_SOURCE */
/* #define _GNU_SOURCE 1 */

/* Define to 4321 if DyLP's fpchecks outputs -D__BIG_ENDIAN=4321 */
/* #define __BIG_ENDIAN */

/* Define to 1 if DyLP's fpchecks outputs -D__DYLP_BROKEN_FPCLASS */
#define __DYLP_BROKEN_FPCLASS

/* Define to 1 if DyLP's fpchecks outputs -D__DYLP_FIX_HUGE_VAL */
/* #define __DYLP_FIX_HUGE_VAL */

/* Define to 1 if DyLP's fpchecks outputs -D__DYLP_SUN */
/* #define __DYLP_SUN */

/* Define to 1 if DyLP's fpchecks outputs -D__DYLP_SUNWspro */
/* #define __DYLP_SUNWspro */

/* Define to 1234 if DyLP's fpchecks outputs -D__LITTLE_ENDIAN=1234 */
/* #define __LITTLE_ENDIAN */

#endif /* HAVE_CONFIG_H */

#endif /*__DYLPCONFIG_H__*/

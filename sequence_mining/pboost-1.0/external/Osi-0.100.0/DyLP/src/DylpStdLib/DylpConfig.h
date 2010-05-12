/*
  This file is part of the support library for the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).

  Include file for the configuration of Dylp.

  On systems where the code is configured with the configure script (i.e.,
  compilation is always done with HAVE_CONFIG_H defined), this header file
  includes the automatically generated header file config_dylp.h, then
  undefines macros that might configure with other ProjConfig.h files.

  On systems that are compiled in other ways (e.g., with the Developer
  Studio), the header file configall_system.h is included to define those
  macros that depend on the operating system and the compiler, followed by
  the defines used for configuration of dylp. A user can modify the default
  settings by editing this file.
*/

#ifndef __DYLPCONFIG_H__
#define __DYLPCONFIG_H__

#ifdef HAVE_CONFIG_H
#include "config_dylp.h"

/*
  Undefine macros that could conflict with those in other config.h files
*/

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#else /* HAVE_CONFIG_H */

/*
  Include the COIN-wide system specific configure header (usually in
  BuildTools/inc).
*/

#include "configall_system.h"

/*
  Defines specific to dylp.
*/

/*
  Define to the C type corresponding to the C++ bool type. `char' is
  correct on many systems. The next most likely choice is int.
*/
#define BOOL char

/*
  Define to the debug sanity check level (0 is no test)
*/
#define COIN_DYLP_CHECKLEVEL 0

/*
  But dylp was developed long before COIN came into being, so if you really
  want the paranoid checks, define DYLP_PARANOIA. The value isn't important.
*/
/* #define DYLP_PARANOIA 1 */

/*
  Define to the debug verbosity level (0 is no output)
*/
#define COIN_DYLP_VERBOSITY 0
/*
  But dylp was developed long before COIN came into being, so if you
  want informational printing, DO NOT define DYLP_NDEBUG. The value isn't
  important.
*/
/* #undef DYLP_NDEBUG 1 */

/*
  Define this variable to enable dylp's statistics collection features.
*/
#define DYLP_STATISTICS 1

/*
  Define to 1 if the DyLP package is available.
*/
#define COIN_HAS_DYLP 1

/*
  Set to the full path directory name for the location of the error text
   message file dy_errmsgs.txt. This file is distributed with dylp source and
   not normally installed elsewhere. An absolute path to DyLP/src/Dylp/ is
   appropriate. The string should end with a directory separator ("/" or "\",
   depending on your system). The surrounding quotes are part of the
   definition. There is no good default; the value given here will work from
   the examples directory, on a windows system, which seems the most likely
   environment to be using this part of DylpConfig.h.
*/
/* #define DYLP_ERRMSGDIR "..\\src\\Dylp\\" */

/*
  Define this symbol if your system is `big-endian', i.e., the most significant
  byte of a multibyte quantity is stored in the lowest byte address. Intel x86
  systems are little-endian. SPARC and Motorola are big-endian.
*/
/* #define WORDS_BIGENDIAN 1 */

/*
  Define this symbol if the quiet_nan function exists. This function should
  return the bit pattern for IEEE quiet NaN.
*/
/* #define DYLP_HAS_QUIET_NAN 1 */

/*
  Define to be the name of the C function used to check that an IEEE floating
  point value is finite. Common possibilities are finite, _finite, and
  isfinite.
*/
#define DYLP_ISFINITE finite

/*
  Define to be the name of the C function used to check that an IEEE floating
  point value is NaN. Common possibilities are isnan and _isnan.
*/
#define DYLP_ISNAN isnan

/*
  Define to 1 if sunmath.h exists. As you might guess, define this only on a
  Sun/Solaris system. And really, if you're building on Sun, why are you
  using this part of the configuration file? Run configure!
*/
/* #define HAVE_SUNMATH_H 1 */

#endif /* HAVE_CONFIG_H */

#endif /*__DYLPCONFIG_H__*/

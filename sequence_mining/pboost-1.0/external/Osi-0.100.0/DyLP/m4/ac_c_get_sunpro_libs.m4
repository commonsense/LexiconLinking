
# AC_DYLP_GET_SUNSTUDIO_LIBDIRS
# ------------------------------------------------------

# Determines the correct set of Sun Studio libraries to add to the final link
# command for libDylp. Why do we need this? Glad you asked. Libtool, in its
# wisdom, elects to issue a bare `ld' command for the final link, instead of
# using cc as the linker. Hence ld does not get the set of library directories
# that cc would normally pass to ld with a -YP option. Which means it doesn't
# find libsunmath, which it needs for quiet_nan. I really should just get rid
# of quiet_nan, but at this point it's an obsession. To get cc to cough up the
# value of the YP spec, we can run it on an empty file with the -xdryrun flag.
# We need to actually do the run because the set of directories will vary
# depending on command line flags. And, just to make life interesting, libtool
# does not grok -YP and happily discards it, if you try to use it directly. So
# we need to convert it into a bunch of -L specs.

# At the end of this macro, SUNSTUDIOLIBDIRS will be set to a sequence of -L
# specs matching the -YP spec.

# ------------------------------------------------------

AC_DEFUN([AC_DYLP_GET_SUNSTUDIO_LIBDIRS],
[ AC_MSG_NOTICE([Determining Sun Studio library directories ... ])
  cat /dev/null > dylp_ac_studiodirs.c

# The autoconf convention seems to be to put output in a file. Avoids problems
# stuffing too-long strings into variables, I suppose.

  $CC -xdryrun $CFLAGS $CPPFLAGS dylp_ac_studiodirs.c 2> dylp_ac_studiodirs.c

# Now find the -YP option. Remember that autoconf does not take kindly to [] in
# macros, so we need to use quadrigraphs.

  SUNSTUDIO_LIBDIRS=`cat dylp_ac_studiodirs.c | \
	 sed -n -e 's/.*-Y@<:@^,@:>@*,\(@<:@^ "@:>@*\).*/-L\1/p' | \
	 sed -n -e 's/:/ -L/gp'`
  rm -f dylp_ac_studiodirs.c

  AC_MSG_NOTICE(["$SUNSTUDIO_LIBDIRS"])
])

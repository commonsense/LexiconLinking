# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 26 2006-05-31 22:26:19Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                List files that should be distributed                 #
########################################################################

EXTRA_DIST = $(EXAMPLE_FILES)

########################################################################
#           Extra targets for uncompressing the .gz files              #
########################################################################

# If configure has been run with --enable-gnu-packages, we assume that
# compressed files can be read and that we don't have to uncompress.
# Otherwise, we uncompress the files
if COIN_HAS_ZLIB
uncompress: skipunzip
else
uncompress: unzip
endif

# This target still leaves the original compressed files around
unzip: $(EXAMPLE_UNCOMPRESSED_FILES)

$(EXAMPLE_UNCOMPRESSED_FILES):
	gzip -d -c $@.gz > $@

skipunzip:
	echo "Skipping decompression (package compiled with --enable-gnu-packages)"

.PHONY: uncompress unzip skipunzip

########################################################################
#                  List files that should be cleaned                   #
########################################################################

CLEANFILES = $(EXAMPLE_UNCOMPRESSED_FILES)

DISTCLEANFILES = $(EXAMPLE_CLEAN_FILES)

include ../../BuildTools/Makemain.inc

# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# $Id: Makefile.in 726 2006-04-17 04:16:00Z andreasw $

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# To compile other examples, either changed the following line, or
# add the argument DRIVER=problem_name to make
DRIVER = driver

# CHANGEME: This should be the name of your executable
EXE = $(DRIVER)

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  $(DRIVER).o 

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

# CHANGEME: Directory to the sources for the (example) problem definition
# files
SRCDIR = .


##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile the      #
#  COIN package.                                                         #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -fomit-frame-pointer -pipe -DNDEBUG -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas  

# additional C++ Compiler options for linking
CXXLINKFLAGS = 

# C Compiler command
CC = cc

# C Compiler options
CFLAGS = -O3 -fomit-frame-pointer -pipe -DNDEBUG -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall -Wno-unknown-pragmas  

# Directory with COIN header files
COININCDIR = /Users/dustin/Documents/Projects/sequence-mining/pboost-1.0/external/Osi-0.100.0/include/coin

# Directory with COIN libraries
COINLIBDIR = /Users/dustin/Documents/Projects/sequence-mining/pboost-1.0/external/Osi-0.100.0/lib

# Libraries necessary to link with Clp
LIBS = -L$(COINLIBDIR) -lClp -lCoinUtils \
	-lm  `cat /Users/dustin/Documents/Projects/sequence-mining/pboost-1.0/external/Osi-0.100.0/lib/../share/doc/coin/CoinUtils/coinutils_addlibs.txt`

# Necessary Include dirs (we use the CYGPATH_W variables to allow
# compilation with Windows compilers)
INCL =  -I`$(CYGPATH_W) $(COININCDIR)` $(ADDINCFLAGS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

# Here we list all possible generated objects or executables to delete them
CLEANFILES = \
	addBits.o addBits \
	addColumns.o addColumns \
	addRows.o addRows \
	decompose.o decompose \
	defaults.o defaults \
	driver2.o driver2 \
	driver.o driver \
	driverC.o driverC \
	dualCuts.o dualCuts \
	ekk.o ekk \
	ekk_interface.o ekk_interface \
	hello.o hello \
	makeDual.o makeDual \
	minimum.o minimum \
	network.o network \
	piece.o piece \
	rowColumn.o rowColumn \
	sprint2.o sprint2 \
	sprint.o sprint \
	testBarrier.o testBarrier \
	testBasis.o testBasis \
	testGub2.o testGub2 \
	testGub.o testGub \
	testQP.o testQP \
	useVolume.o useVolume

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(CLEANFILES)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`

.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.c.obj:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`

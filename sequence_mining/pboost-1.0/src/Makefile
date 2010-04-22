
CPP=g++

INC_BOOST=-I/usr/local/include/boost-1_33_1
LIB_BOOST=-lboost_program_options-1_33_1 -lboost_serialization-1_33_1 
include Makefile.options

# General prequisites:
#    libboost-dev
#    libboost-serialization-dev

# 32 bit shared compilation
CPPFLAGS=-Wall -g -O3 -fomit-frame-pointer \
	-funroll-loops -finline-functions -Wno-deprecated
#CPPFLAGS=-g -ggdb -Wall -Wno-deprecated
INCLUDE=-I. -I../external/osi-0.95.3/include -I../external/osi-0.95.3/Osi/src/OsiClp  -I../external/osi-0.95.3/Osi/src -I/usr/local/include/boost-1_33_1
INCLUDE_PSPAN=-I. -I../external/gmm-3.0/include -I../external/osi-0.95.3/Osi/src -I../external/osi-0.95.3/Osi/src/OsiClp -L/usr/local/lib $(INC_BOOST)
LIBS= -L/usr/local/lib	-L../external/osi-0.95.3/lib -L../external/osi-0.95.3/Osi/lib $(LIB_BOOST) -lOsiClp -lOsi -lClp -lm -lCoinUtils
LIBS_PTEST=-L/usr/local/lib $(LIB_BOOST)
LIBS_PSPAN=-L/usr/local/lib $(LIB_BOOST)

# 64 bit static compilation (cluster nodes do not have the desired libraries)
CPP_64=g++-4.2
CPPFLAGS_64=-Wall -static -m64 -march=athlon64 -O3 -fomit-frame-pointer \
	-funroll-loops -finline-functions -Wno-deprecated -L/usr/lib64
#CPPFLAGS_64=-Wall -static -m64 -g -ggdb -Wno-deprecated
INCLUDE_64=-I. -I../external/boost-1.33.1/include -I../external/osi-0.95.3-amd64/include
LIBS_64=-L../external/boost-1.33.1/lib -lboost_program_options-gcc \
	-lboost_serialization-gcc \
	-ldl -L../external/osi-0.95.3-amd64/lib -lOsiClp -lOsi -lClp -lCoinUtils
LIBS_64_PTEST=-L../external/boost-1.33.1/lib -lboost_serialization-gcc \
	-lboost_program_options-gcc
LIBS_64_PSPAN=-L../external/boost-1.33.1/lib -lboost_serialization-gcc \
	-lboost_program_options-gcc

# pboost sources
SOURCES=Sequence.cpp PrefixSpan.cpp PrefixSpanK.cpp \
	PrefixSpanLPBoost.cpp PrefixSpanLPBoost15.cpp \
	LPBoostOSI.cpp LPBoost15OSI.cpp DDAG.cpp
SOURCES_INC=Sequence.h PrefixSpan.h PrefixSpanK.h \
	PrefixSpanLPBoost.h PrefixSpanLPBoost15.h \
	LPBoostOSI.h LPBoost15OSI.h BestPatternList.h ExplorationGuide2.h DDAG.h \
	Classifier.h

# pspan sources
SOURCES_PSPAN=Sequence.cpp PrefixSpan.cpp PrefixSpanK.cpp PrefixSpanDSPCA.cpp
SOURCES_PSPAN_INC=Sequence.h PrefixSpan.h PrefixSpanK.h PrefixSpanDSPCA.h \
	BestPatternList.h ExplorationGuide2.h


# Matlab MEX wrapper
MEX=${MATLABROOT}/bin/mex
MEXFLAGS=-O -I.

MEX_SOURCES=Sequence.cpp PrefixSpan.cpp PrefixSpanLPBoost.cpp PrefixSpanLPBoost15.cpp
MEX_SOURCES_INC=Sequence.h PrefixSpan.h PrefixSpanLPBoost.h PrefixSpanLPBoost15.h \
	BestPatternList.h ExplorationGuide2.h


### Targets

all:	pboost ptest pspan
# pspan removed, a rebuild is not always required
# SequenceTest removed, outdated

all64:	pboost.amd64 ptest.amd64

SequenceTest:	SequenceTest.cpp $(SOURCES) $(SOURCES_INC) Makefile
	$(CPP) $(CPPFLAGS) -o SequenceTest SequenceTest.cpp
		$(SOURCES) $(INCLUDE) $(LIBS)

PrefixSpanTest:	PrefixSpanTest.cpp $(SOURCES) $(SOURCES_INC) Makefile
	$(CPP) $(CPPFLAGS) -o PrefixSpanTest PrefixSpanTest.cpp \
		$(SOURCES) $(INCLUDE) $(LIBS)

pspan:	pspan.cpp $(SOURCES_PSPAN) $(SOURCES_PSPAN_INC) Makefile
	$(CPP) $(CPPFLAGS) -o pspan pspan.cpp \
		$(SOURCES_PSPAN) $(INCLUDE_PSPAN) $(LIBS_PSPAN)

pboost:	pboost.cpp $(SOURCES) $(SOURCES_INC) Makefile
	$(CPP) $(CPPFLAGS) -o pboost pboost.cpp \
		$(SOURCES) $(INCLUDE) $(LIBS)

ptest:	ptest.cpp DDAG.cpp Sequence.cpp Sequence.h Makefile Classifier.h
	$(CPP) $(CPPFLAGS) -o ptest ptest.cpp Sequence.cpp DDAG.cpp \
		$(INCLUDE) $(LIBS_PTEST)

static:	pspan.x86 pboost.x86 ptest.x86 pspan.amd64 pboost.amd64 ptest.amd64

# static x86 targets
pspan.x86:	pspan.cpp $(SOURCES_PSPAN) $(SOURCES_PSPAN_INC) Makefile
	$(CPP) -static $(CPPFLAGS) -o pspan.x86 pspan.cpp \
		$(SOURCES_PSPAN) $(INCLUDE_PSPAN) $(LIBS_PSPAN)
	mv pspan.x86 ../static/

pboost.x86:	pboost.cpp $(SOURCES) $(SOURCES_INC) Makefile
	$(CPP) -static $(CPPFLAGS) -o pboost.x86 pboost.cpp \
		$(SOURCES) $(INCLUDE) $(LIBS)
	mv pboost.x86 ../static/

ptest.x86:	ptest.cpp DDAG.cpp Sequence.cpp Sequence.h Makefile
	$(CPP) -static $(CPPFLAGS) -o ptest.x86 ptest.cpp Sequence.cpp DDAG.cpp \
		$(INCLUDE) $(LIBS_PTEST)
	mv ptest.x86 ../static/


# static amd64 targets
pspan.amd64:	pspan.cpp $(SOURCES_PSPAN) $(SOURCES_PSPAN_INC) Makefile
	$(CPP_64) $(CPPFLAGS_64) -o pspan.amd64 pspan.cpp \
		$(SOURCES_PSPAN) $(INCLUDE_64) $(INCLUDE_PSPAN) $(LIBS_64_PSPAN)
	mv pspan.amd64 ../static/

pboost.amd64:	pboost.cpp $(SOURCES) $(SOURCES_INC) Makefile
	$(CPP_64) $(CPPFLAGS_64) -o pboost.amd64 pboost.cpp \
		$(SOURCES) $(INCLUDE_64) $(LIBS_64)
	mv pboost.amd64 ../static/

ptest.amd64:	ptest.cpp DDAG.cpp Sequence.cpp Sequence.h Makefile
	$(CPP_64) $(CPPFLAGS_64) -o ptest.amd64 ptest.cpp Sequence.cpp DDAG.cpp \
		$(INCLUDE_64) $(LIBS_64_PTEST)
	mv ptest.amd64 ../static/


# Matlab mex wrapper (may be obsolete by now?)
mex:	mexprefixspan

mexprefixspan:	mexprefixspan.cpp $(MEX_SOURCES) $(MEX_SOURCES_INC) Makefile
	${MEX} ${MEXFLAGS} mexprefixspan.cpp ${MEX_SOURCES}
	cp mexprefixspan.mex* ../gboost/

clean:
	rm -f SequenceTest PrefixSpanTest pspan ptest pboost pboost.amd64



BJAM=./tools/build/jam_src/bin.macosxx86/bjam
TOOLSET=darwin
BJAM_CONFIG= -sPYTHON_ROOT=/opt/local/Library/Frameworks/Python.framework/Versions/2.6 -sPYTHON_VERSION=2.6
PREFIX=/usr/local
EPREFIX=/usr/local
LIBDIR=/usr/local/lib
INCLUDEDIR=/usr/local/include
LIBS=

all: .dummy
	@echo "$(BJAM) $(BJAM_CONFIG) -sTOOLS=$(TOOLSET) $(LIBS)"
	@$(BJAM) $(BJAM_CONFIG) -sTOOLS=$(TOOLSET) $(LIBS) || \
	echo "Not all Boost libraries built properly."

clean: .dummy
	rm -rf bin

distclean: clean
	rm -rf Makefile config.log

check: .dummy
	@cd status && ../$(BJAM) $(BJAM_CONFIG) -sTOOLS=$(TOOLSET) test || echo "Some Boost regression tests failed. This is normal."

install: .dummy
	@echo "$(BJAM) $(BJAM_CONFIG) --prefix=$(PREFIX) --exec-prefix=$(EPREFIX) --libdir=$(LIBDIR) --includedir=$(INCLUDEDIR) -sTOOLS=$(TOOLSET) $(LIBS) install"
	@$(BJAM) $(BJAM_CONFIG) --prefix=$(PREFIX) --exec-prefix=$(EPREFIX) --libdir=$(LIBDIR) --includedir=$(INCLUDEDIR) -sTOOLS=$(TOOLSET) $(LIBS) install || echo "Not all Boost libraries built properly."

.dummy:


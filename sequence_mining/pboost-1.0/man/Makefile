
GROFF=groff

all:	pspan.txt pboost.txt ptest.txt \
	pspan.ps pboost.ps ptest.ps \
	pspan.pdf pboost.pdf ptest.pdf

%.ps:    %.1
	$(GROFF) -t -e -mandoc -Tps $< > $@

%.txt:    %.1
	$(GROFF) -man -Tascii -P-cbuof $< > $@

%.pdf:    %.ps
	ps2pdf $< $@

clean:
	rm -f *.ps *.pdf *.txt



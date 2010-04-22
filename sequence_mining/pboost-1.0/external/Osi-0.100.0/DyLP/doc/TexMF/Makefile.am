# Copyright (C) 2009 Lou Hafer
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Lou Hafer		SFU	2009.04.15

AUTOMAKE_OPTIONS = foreign

# Make sure all the supporting TeX files get distributed. The only thing in
# this subtree that's directly part of the dylp documenation is dylp.bib. All
# the rest is supporting definition files for formats and fonts. This won't
# change, so there's no sense in putting Makefile.am's in the subdirectories.

EXTRA_DIST = $(dylpdoc_BIBTEX) $(dylpdoc_LATEXDEFS) \
	     $(dylpdoc_FONTTFM) $(dylpdoc_FONTVF)

dylpdoc_BIBTEX = bibtex/bst/louplain.bst \
		 bibtex/bib/dylp.bib

dylpdoc_LATEXDEFS = tex/latex/loustandard.sty \
		    tex/latex/loubookman.sty \
		    tex/latex/psnfss \
		    tex/latex/psnfss/omspbk.fd \
		    tex/latex/psnfss/ot1pag.fd \
		    tex/latex/psnfss/t1pbk.fd \
		    tex/latex/psnfss/omxpbk.fd \
		    tex/latex/psnfss/omspplcm.fd \
		    tex/latex/psnfss/t1bch.fd \
		    tex/latex/psnfss/ot1pplcm.fd \
		    tex/latex/psnfss/upbk.fd \
		    tex/latex/psnfss/ot1pbk.fd \
		    tex/latex/psnfss/upcr.fd \
		    tex/latex/psnfss/omlpbk.fd \
		    tex/latex/psnfss/t1pag.fd \
		    tex/latex/psnfss/upsy.fd \
		    tex/latex/psnfss/omlpplcm.fd \
		    tex/latex/psnfss/omxpplcm.fd \
		    tex/latex/psnfss/upag.fd \
		    tex/latex/loumath.sty \
		    tex/latex/codedocn.sty

# Arguably these should be in your local TeX distribution. But they're small.

dylpdoc_FONTTFM = fonts/tfm/pbkd8t.tfm
		  fonts/tfm/omselau.tfm
		  fonts/tfm/pbkdo8t.tfm
		  fonts/tfm/pbkdu8t.tfm
		  fonts/tfm/pbkli7y.tfm
		  fonts/tfm/omselabo.tfm
		  fonts/tfm/pbkdi8t.tfm
		  fonts/tfm/pbkd7y.tfm
		  fonts/tfm/omselabu8t.tfm
		  fonts/tfm/pbkdi7t.tfm
		  fonts/tfm/zpplcmb.tfm
		  fonts/tfm/omselabo8t.tfm
		  fonts/tfm/pbkli7t.tfm
		  fonts/tfm/pbkl7t.tfm
		  fonts/tfm/pbkl7v.tfm
		  fonts/tfm/zpplcmbm.tfm
		  fonts/tfm/omselao.tfm
		  fonts/tfm/zpplcmr.tfm
		  fonts/tfm/zpplcmry.tfm
		  fonts/tfm/omselabu.tfm
		  fonts/tfm/zpplcmby.tfm
		  fonts/tfm/zbkl8t.tfm
		  fonts/tfm/pbkdu8r.tfm
		  fonts/tfm/pbkli7m.tfm
		  fonts/tfm/zpplcmrv.tfm
		  fonts/tfm/zbkd8t.tfm
		  fonts/tfm/omselau8t.tfm
		  fonts/tfm/pbkd7t.tfm
		  fonts/tfm/omselao8t.tfm
		  fonts/tfm/pbklu8r.tfm
		  fonts/tfm/psyro.tfm
		  fonts/tfm/pbklu8t.tfm
		  fonts/tfm/pbkli8t.tfm
		  fonts/tfm/pbkdi7m.tfm
		  fonts/tfm/pbkl7y.tfm
		  fonts/tfm/zpplcmrm.tfm
		  fonts/tfm/pbkdi7y.tfm
		  fonts/tfm/pbklo8t.tfm
		  fonts/tfm/pbkl8t.tfm

# Virtual font definitions. These are the core of replacing the Computer Modern
# fonts with something that looks better --- Bookman is more to my taste.

dylpdoc_FONTVF = fonts/vf/zbkl8t.vf
		 fonts/vf/pbkl7y.vf
		 fonts/vf/pbkdo8t.vf
		 fonts/vf/zpplcmb.vf
		 fonts/vf/pbkdu8t.vf
		 fonts/vf/zpplcmby.vf
		 fonts/vf/pbkli8t.vf
		 fonts/vf/zbkd8t.vf
		 fonts/vf/pbkdi7t.vf
		 fonts/vf/zpplcmrv.vf
		 fonts/vf/zpplcmry.vf
		 fonts/vf/pbkdi7m.vf
		 fonts/vf/pbkdi7y.vf
		 fonts/vf/pbkd8t.vf
		 fonts/vf/pbkl8t.vf
		 fonts/vf/pbkd7t.vf
		 fonts/vf/pbkl7t.vf
		 fonts/vf/zpplcmrm.vf
		 fonts/vf/zpplcmr.vf
		 fonts/vf/pbkli7t.vf
		 fonts/vf/zpplcmbm.vf
		 fonts/vf/pbkl7v.vf
		 fonts/vf/pbkli7m.vf
		 fonts/vf/pbklo8t.vf
		 fonts/vf/pbkd7y.vf
		 fonts/vf/pbklu8t.vf
		 fonts/vf/pbkli7y.vf
		 fonts/vf/pbkdi8t.vf

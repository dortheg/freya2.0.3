#
# Included by Makefile in parent directory
#
# Doug Wright

#....location to store final documents
DOCDIR := ../../doc


#....keep up to date pdf file in the doc directory
$(DOCDIR)/$(file).pdf: $(file).pdf
	cp $(file).pdf $(DOCDIR)
#
#....rule for making pdf files
#
%.pdf:%.tex
	latex $* && latex $*  # do twice in case cross-references need to be updated
        #.....can not use dvipdf on ilx because the default is A4 paper and
        #     can not be over-ruled on the command line
	dvips -t letter -q -f $* | gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$*.pdf -c .setpdfwrite -
#
#
#....create pdf file directly from dvi (available if tetex installed)
#      Use this if you use rotating.sty, regular dvipdf doesn't do it right
#
#	latex $*
#	dvipdfm $*

#
#....create gif files from eps figures (suitable for Powerpoint)
#      Don't run this on ilx, much smaller files on MacOSX 
#      the Mac uses AFPL ghostscript instead of Gnu
#
eps := $(wildcard eps/*.eps)
gif: $(eps:.eps=.gif)

%.gif: %.eps
	convert -density 144 $*.eps $*.gif

#
#....clean up tex junk
#
clean:
	rm *.aux *.dvi *.log

.PHONY: default pdf gif clean

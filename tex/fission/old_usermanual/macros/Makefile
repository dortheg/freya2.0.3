#
# Included by Makefile in parent directory
#
# Doug Wright

doc := $(file:.tex=)

export TEXINPUTS := $(TEXINPUTS):macros

default: pdf

#
#....create pdf file
#
pdf:
	latex $(doc)
        #.....can not use dvipdf on ilx because the default is A4 paper and
        #     can not be over-ruled on the command line
	dvips -t letter -q -f $(doc) | gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$(doc).pdf -c .setpdfwrite -

#
#....create pdf file directly from dvi (available if tetex installed)
#      Use this if you use rotating.sty, regular dvipdf doesn't do it right
#
pdfm:
	latex $(doc)
	dvipdfm $(doc)

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

.PHONY: default pdf pdfm eps gif clean

#
# make        Run latex and generate new pdf document
# make gif    Generate gifs of everything in eps directory
# make clean  Delete the LaTeX junk files (rm *.aux *.dvi *.log)
#
# Doug Wright

file = ieee_abstract

include ../macros/Makefile

default: $(DOCDIR)/$(file).pdf

common := ieee_abstract_bib.tex

ieee_abstract.pdf:       ieee_abstract.tex $(common)

#....make sure that LaTex finds our special macros in the local area
export TEXINPUTS := $(TEXINPUTS):../macros

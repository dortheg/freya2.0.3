#
# make             Run pdflatex and generate new pdf documents if needed
#
# pdflatex fission         make the general usermanual (by hand)
# pdflatex fission_geant   make the documentation for geant (by hand)
#
# Doug Wright

default: fission.pdf fission_geant.pdf ../../doc/fission.pdf

common := physics.tex photofission.tex fission_geant_bib.tex implementation.tex snippet.tex arg.tex

fission.pdf:       fission.tex $(common) usermanual.tex geant4.tex mcnpx.tex interface.tex cog.tex

../../doc/fission.pdf:	fission.pdf
	cp fission.pdf ../../doc

fission_geant.pdf: fission_geant.tex $(common) 

%.pdf: %.tex
	pdflatex $*; pdflatex $*

set title "total photon energy versus fission fragment mass, thermal neutrons on Pu-239"

set style line 1 linecolor 1 linetype 10 linewidth 3 pointtype 6

set xlabel "fission fragment mass (amu)"
set ylabel "total photon energy (MeV)"

plot "../build/totpenergy_vs_Af.res" using 1:3:5 with errorbars ls 1 notitle

set out "totpenergy_vs_Af.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "totpenergy_vs_Af.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


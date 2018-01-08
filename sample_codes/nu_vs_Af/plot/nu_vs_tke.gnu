set title "Nubar vs total kinetic energy of fission fragments, thermal neutrons on Pu-239"

set style line 1 linecolor 1 linetype 10 linewidth 3 pointtype 6

set xlabel "total kinetic energy of fission fragments (MeV)"
set ylabel "mean number of neutrons"
set xrange [100:210]

plot "../build/nu_vs_TKE.res" using 1:6:8 with errorbars ls 1 notitle

set out "nu_vs_tke.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "nu_vs_tke.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


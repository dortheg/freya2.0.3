set title "Mean number of neutrons per fission fragment, thermal neutrons on Pu-239"

set style line 1 linecolor 1 linetype 10 linewidth 3 pointtype 6

set xlabel "fission fragment mass A"
set ylabel "mean number of neutrons"
set xrange [60:190]

plot "../build/nu_vs_Af.res" using 1:3:5 with errorbars ls 1 notitle

set out "nu_vs_Af.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "nu_vs_Af.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


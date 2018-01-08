set title "Neutron-neutron correlation versus angle, thermal neutrons on Pu-239"

set style line 1 linecolor 1 linetype 10 linewidth 3 pointtype 6

set xlabel "direction cosine between fission neutrons"
set ylabel "n-n correlation (arbitrary units)"
set xrange [-1:+1]

plot "../build/nnCorr_vs_mu.res" using 1:5:7 with errorbars ls 1 notitle

set out "nnCorr.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "nnCorr.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


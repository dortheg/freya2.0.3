set title "Direction cosine between fission neutrons, 2 MeV neutrons on Pu-239"

set style line 1 linecolor 1 linetype 10 linewidth 3

set xlabel "mu"
set ylabel "probability"

# set format y "%.1te%+2T"
set format y "%.1t 10^{%+2T}"

plot "../build/angular_correlation.res" using ($1+$3)/2:5 with histeps ls 1 notitle

set out "angular_correlation.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "angular_correlation.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


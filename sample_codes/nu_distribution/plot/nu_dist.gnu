set title "Probability of nu neutrons per fission, 1 MeV neutrons on U-235"

set style line 1 linecolor 1 linetype 10 linewidth 3

set xlabel "nu"
set ylabel "probability"
# set logscale y

# set format y "%.1te%+2T"
# set format y "10^{%+2T}"

plot "../build/nu_dist.res" using 1:3 with histeps ls 1 notitle

set out "nu_dist.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "nu_dist.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


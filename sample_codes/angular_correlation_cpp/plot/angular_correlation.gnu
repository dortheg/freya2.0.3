set title "Direction cosine between fission neutrons, spontaneous fission of Pu-240"

set style line 1 linecolor 1 linetype 10 linewidth 3
set style line 2 linecolor 2 linetype 10 linewidth 3
set style line 3 linecolor 3 linetype 10 linewidth 3
set style line 4 linecolor 4 linetype 10 linewidth 3

set xlabel "mu"
set ylabel "probability"

# set format y "%.1te%+2T"
set format y "%.1t 10^{%+2T}"

plot "../build/angular_correlation.res.1.1" using ($1+$3)/2:5 with histeps ls 1 title "x=1.1", \
     "../build/angular_correlation.res.1.2" using ($1+$3)/2:5 with histeps ls 2 title "x=1.2", \
     "../build/angular_correlation.res.1.3" using ($1+$3)/2:5 with histeps ls 3 title "x=1.3", \
     "../build/angular_correlation.res.1.4" using ($1+$3)/2:5 with histeps ls 4 title "x=1.4"

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

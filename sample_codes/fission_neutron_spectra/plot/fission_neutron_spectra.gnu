set title "Fission spectra for 1.85 MeV neutrons on Pu-239"

set style line 1 linecolor 1 linetype 10 linewidth 1
set style line 2 linecolor 2 linetype 10 linewidth 1
set style line 3 linecolor 3 linetype 10 linewidth 1
set style line 4 linecolor 4 linetype 10 linewidth 1
set style line 5 linecolor 5 linetype 10 linewidth 1
set style line 6 linecolor 6 linetype 10 linewidth 1
set style line 7 linecolor 7 linetype 10 linewidth 1
set style line 8 linecolor 8 linetype 10 linewidth 1

set xlabel "energy [MeV]"
set ylabel "probability"
set logscale x
set logscale y

# set format y "%.1te%+2T"
set format y "10^{%+2T}"

set key left top

plot \
"../build/nu_dist_n=1.res" using 1:5 with histeps ls 1 title "1 neutron emitted", \
"../build/nu_dist_n=2.res" using 1:5 with histeps ls 2 title "2 neutron emitted", \
"../build/nu_dist_n=3.res" using 1:5 with histeps ls 3 title "3 neutron emitted", \
"../build/nu_dist_n=4.res" using 1:5 with histeps ls 4 title "4 neutron emitted", \
"../build/nu_dist_n=5.res" using 1:5 with histeps ls 5 title "5 neutron emitted"

# "../build/nu_dist_n=6.res" using 1:5 with histeps ls 6 title "6 neutron emitted"
# "../build/nu_dist_n=7.res" using 1:5 with histeps ls 7 title "7 neutron emitted"
# "../build/nu_dist_n=8.res" using 1:5 with histeps ls 8 title "8 neutron emitted"

set out "neutron_spectra_vs_multiplicity.eps"
set terminal postscript enhanced font 'sans,18' color eps
replot
set terminal pop
replot
set out "neutron_spectra_vs_multiplicity.pdf"
set terminal pdf enhanced font 'sans,12' color
replot
set terminal pop
replot


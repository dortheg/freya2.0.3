if(!exists("iso")) iso = 92235
if(!exists("engindex")) engindex = 0
if(!exists("FreyaEng")) FreyaEng = 0.53
if(!exists("engCons")) engCons = 2.53e-8

set xrange [1e-1:100]
set yrange [1e-7:1e-2]
set xlabel "Total fission neutron and gamma-ray energy [MeV]"
set logs x
set logs y

if(!exists("message")) message = 1; print "Press e to change energy"

plot \
"../regr/testFreyaTotal.res" using 3:((column(1) == iso && column(2) == FreyaEng && column(4) != 0.)?column(4):1/0) with points title sprintf("Total neutron and gamma-ray energy, FREYA, isotope = %g FreyaEng = %g MeV", iso, FreyaEng), \
"../regr/testEngConsTotal.goldref" using 3:((column(1) == iso && column(2) == engCons && column(4) != 0.)?column(4):1/0) with points title sprintf("Total neutron and gamma-ray energy, energy conservation, isotope = %g engCons = %g MeV", iso, engCons)

bind e "engindex=(engindex<1)?engindex+1:0; engCons=(engindex==0)?2.53e-8:(engindex==1)? 1.:1.; load 'specFreyaCorrelatedTotal.gnuplot'"


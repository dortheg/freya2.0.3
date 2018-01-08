if(!exists("iso")) iso = 92235
if(!exists("engindex")) engindex = 0
if(!exists("FreyaEng")) FreyaEng = 0.53
if(!exists("engWatt")) engWatt = 0.
if(!exists("scalingn")) scalingn = 1.77
if(!exists("scalingp")) scalingp = 3.72

set logs x
set logs y

if(!exists("message")) message = 1; print "Press e to change energy"

plot \
"../regr/testspec.goldref" using 4:((column(1) == iso && column(3) == engWatt && "n" eq strcol(2))?column(5):1/0) with linespoints title sprintf("Original Watt spectrum, isotope = %g engWatt = %g MeV", iso, engWatt), \
"../regr/testFreyaSpec.goldref" using 3:((column(1) == iso && column(2) == FreyaEng)?scalingn*column(4):1/0) with linespoints title sprintf("Gold ref neutron spectrum, FREYA, isotope = %g FreyaEng = %g MeV", iso, FreyaEng), \
"../regr/testFreyaSpec.res" using 3:((column(1) == iso && column(2) == FreyaEng)?scalingn*column(4):1/0) with linespoints title sprintf("Neutron spectrum, FREYA, isotope = %g FreyaEng = %g MeV", iso, FreyaEng), \
"../regr/testspec.goldref" using 4:((column(1) == iso && column(3) == engWatt && "p" eq strcol(2))?column(5):1/0) with linespoints title sprintf("Original fission gamma spectrum, isotope = %g engWatt = %g MeV", iso, engWatt), \
"../regr/testFreyaSpec.goldref" using 3:((column(1) == iso && column(2) == FreyaEng)?scalingp*column(5):1/0) with linespoints title sprintf("Gold ref photon spectrum, FREYA, isotope = %g FreyaEng = %g MeV", iso, FreyaEng), \
"../regr/testFreyaSpec.res" using 3:((column(1) == iso && column(2) == FreyaEng)?scalingp*column(5):1/0) with linespoints title sprintf("Photon spectrum, FREYA, isotope = %g FreyaEng = %g MeV", iso, FreyaEng)

bind e "engindex=(engindex<1)?engindex+1:0; engWatt=(engindex==0)?0.:(engindex==1)?1.:1.; load 'testspecFreya.gnuplot'"


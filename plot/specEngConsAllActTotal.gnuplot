if(!exists("iso")) iso = 90231
if(!exists("idx")) idx = 0
if(!exists("engindex")) engindex = 0
if(!exists("engConsEng")) engConsEng = 25.3E-9

set logs x
set logs y

if(!exists("message")) message = 1; print "Press m to change isotope\nPress e to change the energy"

plot \
"../regr/testEngConsAllActSpecTotal.res" using 3:((column(1) == iso && column(2) == engConsEng && column(4) != 0.)?column(4):1/0) with linespoints title sprintf("Total fission neutron energy, energy conservation, all actinides, isotope = %g engConsEng = %g MeV", iso, engConsEng), \
"../regr/testEngConsSpecTotal.res.zeroCorrelation" using 3:((column(1) == iso && column(2) == engConsEng && column(4) != 0.)?column(4):1/0) with linespoints title sprintf("Total fission neutron energy, independent energy sampling, isotope = %g engConsEng = %g MeV", iso, engConsEng), \
"../regr/testEngConsAllActSpecTotal.res" using 3:((column(1) == iso && column(2) == engConsEng && column(5) != 0.)?column(5):1/0) with linespoints title sprintf("Total photon energy, energy conservation, all actinides, isotope = %g engConsEng = %g MeV", iso, engConsEng), \
"../regr/testEngConsSpecTotal.res.zeroCorrelation" using 3:((column(1) == iso && column(2) == engConsEng && column(5) != 0.)?column(5):1/0) with linespoints title sprintf("Total photon energy, independent energy sampling, isotope = %g engConsEng = %g MeV", iso, engConsEng)

bind m "idx=(idx==5)?0:idx+1; iso=(idx==0)?90231:(idx==1)?91233:(idx==2)?92235:(idx==3)?92238:(idx==4)?94239:(idx==5)?98252:90231; load 'specEngConsAllActTotal.gnuplot'"

bind e "engindex=(engindex<10)?engindex+1:0; engConsEng=(engindex==0)?25.3E-9:(engindex==1)?1.:(engindex==2)?2.:(engindex==3)?3.:(engindex==4)?4.:(engindex==5)?5.:(engindex==6)?6.:(engindex==7)?7.:(engindex==8)?8.:(engindex==9)?9.:(engindex==10)?10.:10.; engWatt=(engindex==0)?0.:(engindex==1)?1.:(engindex==2)?2.:(engindex==3)?3.:(engindex==4)?4.:(engindex==5)?5.:(engindex==6)?6.:(engindex==7)?7.:(engindex==8)?8.:(engindex==9)?9.:(engindex==10)?10.:10.; load 'specEngConsAllActTotal.gnuplot'"


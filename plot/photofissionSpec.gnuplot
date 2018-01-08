if(!exists("iso")) iso = 92232
if(!exists("eng")) eng = 0
if(!exists("idx")) idx = 0
if(!exists("energy")) energy = 1.00002E-11
if(!exists("engindex")) engindex = 0

set logs x
set logs y

if(!exists("message")) message = 1; print "press n to change the energy of the sampled distribution\npress m to change isotope\nPress e to change the energy of the U-232 data from the endfb.VII library"

plot "../regr/testPhotofissionSpec.res" using 4:((column(1) == iso && int(column(3)) == eng)?column(5):1/0) with linespoints title sprintf("isotope = %g eng = %g MeV", iso, eng), \
"U-232.txt" using 3:(0.1*column(4)) index engindex with linespoints title sprintf("U-232 energy = %g MeV", energy)


bind n "eng=(eng==0)?1:(eng==1)?2:(eng==2)?3:(eng==3)?4:(eng==4)?5:(eng==5)?6:(eng==6)?7:(eng==7)?8:(eng==8)?9:(eng==9)?10:(eng==10)?0:0; load 'photofissionSpec.gnuplot'"

bind m "eng=0; inc=(idx==0)?1:(idx==8)?-1:inc; idx=idx+inc; iso=(idx==0)?92232:(idx==1)?92233:(idx==2)?92234:(idx==3)?92235:(idx==4)?92236:(idx==5)?92238:(idx==6)?94239:(idx==7)?94241:(idx==8)?94240:(idx==9)?98252:98252; load 'photofissionSpec.gnuplot'"
bind e "engindex=(engindex<6)?engindex+1:0; energy=(engindex==0)?1.00002E-11:(engindex==1)?0.0010:(engindex==2)?0.01:(engindex==3)?0.1:(engindex==4)?1.:(engindex==5)?10.:(engindex==6)?20.:20.; load 'photofissionSpec.gnuplot'"

if(!exists("idx")) idx = 0
if(!exists("isotope")) isotope = 90231
if(!exists("eng")) eng = 0
if(!exists("isoindex")) isoindex = 0

if(!exists("message")) message=1; print "Press m to change isotope\nPress n to change energy\n"

plot "../regr/testEngConsAllActNudist.res" using 5:((column(1) == isotope && int(column(3)) == eng)?column(6):1/0) with linespoints title sprintf("isotope (%g) = %g energy = %g MeV", isoindex, isotope, eng)

bind n "eng=(eng==0)?1:(eng==1)?2:(eng==2)?3:(eng==3)?4:(eng==4)?5:(eng==5)?6:(eng==6)?7:(eng==7)?8:(eng==8)?9:(eng==9)?10:(eng==10)?0:0; load 'nudistEngConsAllAct.gnuplot'"

bind m "idx=(idx==5)?0:idx+1; isotope=(idx==0)?90231:(idx==1)?91233:(idx==2)?92235:(idx==3)?92238:(idx==4)?94239:(idx==5)?98252:90231; load 'nudistEngConsAllAct.gnuplot'"


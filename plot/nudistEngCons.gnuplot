if(!exists("idx")) idx = 0
if(!exists("isotope")) isotope = 92235
if(!exists("eng")) eng = 0
if(!exists("isoindex")) isoindex = 0

if(!exists("message")) message=1; print "Press m to change isotope\nPress n to change energy\n"

plot "../regr/testEngConsNudist.res" using 5:((column(1) == isotope && int(column(3)) == eng)?column(6):1/0) with linespoints title sprintf("isotope (%g) = %g energy = %g MeV", isoindex, isotope, eng)

bind n "eng=(eng==0)?1:(eng==1)?2:(eng==2)?3:(eng==3)?4:(eng==4)?5:(eng==5)?6:(eng==6)?7:(eng==7)?8:(eng==8)?9:(eng==9)?10:(eng==10)?0:0; load 'nudistEngCons.gnuplot'"

bind m "idx=(idx==2)?0:idx+1;isotope=(idx==0)?92235:(idx==1)?92238:(idx==2)?94239:92235; load 'nudistEngCons.gnuplot'"

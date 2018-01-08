if(!exists("isotope")) isotope = 92232
if(!exists("nudist")) nudist = 0
if(!exists("eng")) eng = 0
if(!exists("isoindex")) isoindex = 0

if(!exists("message")) message=1; print "Press m to change isotope\nPress n to change energy\nPress , to change nudist option"

plot "../regr/testNuDist.res" using 5:((column(1) == isotope && int(column(3)) == eng && int(column(4)) == nudist)?column(6):1/0) with linespoints title sprintf("isotope (%g) = %g energy = %g MeV, nudist option = %g", isoindex, isotope, eng, nudist)

bind n "eng=(eng==0)?1:(eng==1)?2:(eng==2)?3:(eng==3)?4:(eng==4)?5:(eng==5)?6:(eng==6)?7:(eng==7)?8:(eng==8)?9:(eng==9)?10:(eng==10)?0:0; load 'nudist.gnuplot'"

bind , "nudist=(nudist==0)?1:(nudist==1)?2:(nudist==2)?3:(nudist==3)?0:0; load 'nudist.gnuplot'"

bind m "eng=0; isoindex=(isoindex<9)?isoindex+1:0; isotope=(isoindex==0)?92232:(isoindex==1)?92233:(isoindex==2)?92234:(isoindex==3)?92235:(isoindex==4)?92236:(isoindex==5)?92238:(isoindex==6)?94239:(isoindex==7)?94241:(isoindex==8)?94240:(isoindex==9)?98252:98252; load 'nudist.gnuplot'"

if(!exists("isotope")) isotope = 92238
if(!exists("isoindex")) isoindex = 0
if(!exists("isoinc")) isoinc = 1

if(!exists("message")) message=0; print "Press m to change isotope"

plot "../regr/testSpNuDist.res" using 3:((column(1) == isotope)?column(4):1/0) with linespoints title sprintf("isotope(%g) = %g", isoindex, isotope)

bind m "isoindex=(isoindex<18)?isoindex+isoinc:0; isotope=(isoindex == 0)?90232:(isoindex == 1)?92232:(isoindex == 2)?92233:(isoindex == 3)?92234:(isoindex == 4)?92235:(isoindex == 5)?92236:(isoindex == 6)?92238:(isoindex == 7)?93237:(isoindex == 8)?94236:(isoindex == 9)?94238:(isoindex == 10)?94239:(isoindex == 11)?94240:(isoindex == 12)?94241:(isoindex == 13)?94242:(isoindex == 14)?95241:(isoindex == 15)?96242:(isoindex == 16)?96244:(isoindex == 17)?97249:(isoindex == 18)?98252:98252; load 'spnudist.gnuplot'"

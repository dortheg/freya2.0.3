Watt(a,b,c,x) = c*sqrt(pi*b/4/a)*exp(b/4/a)/a*exp(-1.*x*a)*sinh(sqrt(b*x))

if(!exists("isotope")) isotope = 90232
if(!exists("isoindex")) isoindex = 0
if(!exists("isoinc")) isoinc = 1

set logs x
set logs y

if(!exists("message")) message = 1; print "press m to change isotope"

plot "../regr/testSpspec.res" using 3:((column(1) == isotope)?column(4):1/0) with linespoints title sprintf("isotope (%g) = %g", isoindex, isotope), \
"../regr/testSpspec.res" using 3:(Watt(1.12082e+00, 3.72278e+00, .0093, $3)) title "Analytical Watt fision spectrum for spontaneous fission in U-232"

bind m "isoindex=(isoindex<18)?isoindex+isoinc:0; isotope=(isoindex == 0)?90232:(isoindex == 1)?92232:(isoindex == 2)?92233:(isoindex == 3)?92234:(isoindex == 4)?92235:(isoindex == 5)?92236:(isoindex == 6)?92238:(isoindex == 7)?93237:(isoindex == 8)?94236:(isoindex == 9)?94238:(isoindex == 10)?94239:(isoindex == 11)?94240:(isoindex == 12)?94241:(isoindex == 13)?94242:(isoindex == 14)?95241:(isoindex == 15)?96242:(isoindex == 16)?96244:(isoindex == 17)?97249:(isoindex == 18)?98252:98252; load 'spspec.gnuplot'"

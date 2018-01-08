a(x)=a2+a1*x+a0*x**2

fit a(x) 'a_vs_eng.txt' using 1:2 via a0,a1,a2

print "a0=", a0, " a1=", a1, " a2=", a2

plot "a_vs_eng.txt" using 1:2 with linespoints title "data", "a_vs_eng.txt" using 1:(a(column(1))) with lines title "fit"


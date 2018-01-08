Watt(a,b,c,x) = c*sqrt(pi*b/4/a)*exp(b/4/a)/a*exp(-a*x)*sinh(sqrt(b*x))

b=1.

if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 0 via a,c 
if(!exists("i")) a0 = a; b0 = b; c0 = c
if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 1 via a,c 
if(!exists("i")) a1 = a; b1 = b; c1 = c
if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 2 via a,c 
if(!exists("i")) a2 = a; b2 = b; c2 = c
if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 3 via a,c 
if(!exists("i")) a3 = a; b3 = b; c3 = c
if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 4 via a,c 
if(!exists("i")) a4 = a; b4 = b; c4 = c
if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 5 via a,c 
if(!exists("i")) a5 = a; b5 = b; c5 = c
if(!exists("i")) fit Watt(a,b,c,x) 'U-232.txt' using 3:4 index 6 via a,c 
if(!exists("i")) a6 = a; b6 = b; c6 = c
if(!exists("i")) print a0, "\t", b0, "\t", c0
if(!exists("i")) print a1, "\t", b1, "\t", c1
if(!exists("i")) print a2, "\t", b2, "\t", c2
if(!exists("i")) print a3, "\t", b3, "\t", c3
if(!exists("i")) print a4, "\t", b4, "\t", c4
if(!exists("i")) print a5, "\t", b5, "\t", c5
if(!exists("i")) print a6, "\t", b6, "\t", c6

if(!exists("i")) i=0
if (i==0) a=a0; b=b0; c=c0
if (i==1) a=a1; b=b1; c=c1
if (i==2) a=a2; b=b2; c=c2
if (i==3) a=a3; b=b3; c=c3
if (i==4) a=a4; b=b4; c=c4
if (i==5) a=a5; b=b5; c=c5
if (i==6) a=a6; b=b6; c=c6

plot "U-232.txt" using 3:4 index i title "data ".int(i), "U-232.txt" using 3:(Watt(a,b,c,column(3))) index i title "Watt ".int(i)

bind > "a=a*1.01; replot"
bind < "a=a/1.01; replot"

bind n "inc=(i==0)?1:(i==6)?-1:inc; i=i+inc; load 'U-232.gnuplot'"

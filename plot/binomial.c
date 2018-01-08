#include <math.h>

#define maxp 40

double logbin(int n, int N) {
   double sum=0.;
   int i;

   if (n==0) return 0;
   for(i=2; i<=N; i++) sum += log((double) i);
   for(i=2; i<=n; i++) sum -= log((double) i);
   for(i=2; i<=N-n; i++) sum -= log((double) i);
   return sum;
}

int main() {
   int i=0;
   double p=0.765156;
   double q=1.-p;
   double logcoeff[maxp];
   double logpi[maxp];

   printf("p=%e, q=%e\n", p, q);
   for (i=0; i<=maxp; i++) {
      logcoeff[i] = logbin(i, i+25);
      printf("log[binomial(%d,%d)]=%20.14e\n", i, i+25, logcoeff[i]);
   }
   for (i=0; i<=maxp; i++) {
      logpi[i] = logcoeff[i]+26.*log(p)+i*log(q);
      printf("pi( %d ) = %e\n", i, exp(logpi[i]));
   }

   return 0;
}

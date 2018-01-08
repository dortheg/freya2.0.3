#define iterations 3000000
#define nbins 100

#include <stdio.h>
#include <stdlib.h>
#include "Fission.h"

void init(void);
FILE* openfile(char* name);
void output(int* hist);

int main() {
   int isotope = 94239;
   double energy_MeV = 2.;
   double nubar = 3.163;
   double time = 0.;

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   int i, hist[nbins];
   for (i=0; i<nbins; i++) hist[i] = 0.;

   init();
   for (i=0; i<iterations; i++) {
      genfissevt_(&isotope, &time, &nubar, &energy_MeV);
      int errorlength=maxerrorlength;
      getfreya_errors_(&errorlength, &errors[0]);
      if (errorlength>1) {
        printf("%s\n",errors);
        exit(1);
      }
      int nneutrons = getnnu_();
      int n1;
      for(n1=0; n1<nneutrons; n1++) {
         double u1 = getndircosu_(&n1), v1 = getndircosv_(&n1), w1 = getndircosw_(&n1);
         int n2;
         for(n2=n1+1; n2<nneutrons; n2++) {
            double u2 = getndircosu_(&n2), v2 = getndircosv_(&n2), w2 = getndircosw_(&n2);
            double scalar_product = u1*u2+v1*v2+w1*w2;
        
            int bin_index = (int) (nbins*(scalar_product+1)/2);
            hist[bin_index]++;
         }
      }
   }
   output(hist);
}

void init(void) {
   unsigned short int s[3] = {1234, 5678, 9012};
   int i, three = 3;
   seed48(s);
   setcorrel_(&three);
   return;
}

FILE* openfile(char* name) {
   FILE* fp = fopen(name, "w");
   if (fp == (FILE *) 0) fprintf(stderr, "Could not open %s for writing", name);
   return fp;
}

void output(int* hist) {
   char filename [1024];
   sprintf(filename, "angular_correlation.res");
   FILE* fp = openfile(filename);

   int i;
   unsigned int sum=0;
   for (i=0; i<nbins; i++) sum += hist[i];
   for (i=0; i<nbins; i++) fprintf(fp, "%e - %e : %e\n", -1+2.*i/nbins, -1+2.*(i+1)/nbins, 1.*hist[i]/sum);

   fclose(fp);
   return;
}


#define iterations 300000
#define nbins 100

#include <stdio.h>
#include "fissionEvent.h"

void init(void);
FILE* openfile(char* name);
void output(int* hist);

int main() {
   bool spontfiss=false;
   int isotope = 94239;
   double energy_MeV = 2.;
   double nubar = 3.163;
   double time = 0.;

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   int hist[nbins];
   for (int i=0; i<nbins; i++) hist[i] = 0.;

   init();
   for (int i=0; i<iterations; i++) {
      fissionEvent* fe = new fissionEvent(isotope, time, nubar, energy_MeV, (spontfiss)?0:1);
      int errorlength=maxerrorlength;
      fe->getFREYAerrors(&errorlength, &errors[0]);
      if (errorlength>1) {
         printf("%s\n",errors);
         exit(1);
      }
      int nneutrons = fe->getNeutronNu();
      for(int n1=0; n1<nneutrons; n1++) {
         double u1 = fe->getNeutronDircosu(n1), v1 = fe->getNeutronDircosv(n1), w1 = fe->getNeutronDircosw(n1);
         for(int n2=n1+1; n2<nneutrons; n2++) {
            double u2 = fe->getNeutronDircosu(n2), v2 = fe->getNeutronDircosv(n2), w2 = fe->getNeutronDircosw(n2);
            double scalar_product = u1*u2+v1*v2+w1*w2;
        
            int bin_index = (int) (nbins*(scalar_product+1)/2);
            hist[bin_index]++;
         }
      }
      delete fe;
   }
   output(hist);
}

void init(void) {
   unsigned short int s[3] = {1234, 5678, 9012};
   int i;
   seed48(s);
   fissionEvent::setCorrelationOption(3);
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

   unsigned int sum=0;
   for (int i=0; i<nbins; i++) sum += hist[i];
   for (int i=0; i<nbins; i++) fprintf(fp, "%e - %e : %e\n", -1+2.*i/nbins, -1+2.*(i+1)/nbins, 1.*hist[i]/sum);

   fclose(fp);
   return;
}


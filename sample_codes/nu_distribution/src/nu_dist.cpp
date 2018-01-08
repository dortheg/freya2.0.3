// #define iterations 30000000
#define iterations 3000000
#define nbins 15

#include <stdio.h>
#include "fissionEvent.h"

void init(void);
FILE* openfile(char* name);
void output(int* hist);

int main() {
   bool spontaneous_fission=true;
   int isotope = 98252;
   double energy_MeV = 2.1;
   double nubar = 2.523670;
   double time = 0.;

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   int hist[nbins];
   for (int i=0; i<nbins; i++) hist[i] = 0.;

   init();
   for (int i=0; i<iterations; i++) {
      fissionEvent* fe = new fissionEvent(isotope, time, nubar, energy_MeV, (spontaneous_fission)?0:1);
      int errorlength=maxerrorlength;
      fe->getFREYAerrors(&errorlength, &errors[0]);
      if (errorlength>1) {
         printf("%s\n",errors);
         exit(1);
      }
      int nneutrons = fe->getNeutronNu();
      hist[nneutrons]++;
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
   sprintf(filename, "nu_dist.res");
   FILE* fp = openfile(filename);

   unsigned int sum=0;
   for (int i=0; i<nbins; i++) sum += hist[i];
   for (int i=0; i<nbins; i++) fprintf(fp, "%d : %10.8f\n", i, 1.*hist[i]/sum);

   for (int i=0; i<nbins; i++) 
     printf("nu[%d]=%g\n", i, 1.*hist[i]/sum);

   double nu_bar=0;
   for (int i=1; i<nbins; i++) nu_bar += 1.*i*hist[i]/sum;
   printf("nu_bar=%g\n", nu_bar);
   double nu_2=0;
   for (int i=2; i<nbins; i++) nu_2 += 0.5*i*(i-1)*hist[i]/sum;
   printf("nu2=%g\n", nu_2);
   double nu_3=0;
   for (int i=3; i<nbins; i++) nu_3 += 1./6*i*(i-1)*(i-2)*hist[i]/sum;
   printf("nu3=%g\n", nu_3);
   double nu_4=0;
   for (int i=4; i<nbins; i++) nu_4 += 1./6/4*i*(i-1)*(i-2)*hist[i]/sum;
   printf("nu4=%g\n", nu_4);
   double nu_5=0;
   for (int i=5; i<nbins; i++) nu_5 += 1./6/4/5*i*(i-1)*(i-2)*hist[i]/sum;
   printf("nu5=%g\n", nu_5);

   printf("D2=%g\n", nu_2/nu_bar);
   printf("D3=%g\n", nu_3/nu_bar);
   printf("D4=%g\n", nu_4/nu_bar);
   printf("D5=%g\n", nu_5/nu_bar);

   fclose(fp);
   return;
}


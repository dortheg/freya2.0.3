#define iterations 30000000
#define maxneutrons 8

#define minEnergy 1e-4 // MeV
#define maxEnergy 1e2  // MeV
#define nbinsperspectrum 600

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "fissionEvent.h"

void init(void);
FILE* openfile(char* name);
void output(double f, int hist[][nbinsperspectrum]);

int main() {
   int isotope = 94239;
   double energy_MeV = 1.85;
   double nubar = 3.00888;
   double time = 0.;
   double f = pow(maxEnergy/minEnergy,1./nbinsperspectrum);
   double logf = log(f);

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   int hist[maxneutrons+1][nbinsperspectrum];
   for (int i=0; i<maxneutrons+1; i++)
      for (int j=0; j<nbinsperspectrum; j++)
         hist[i][j] = 0.;

   init();
   for (int i=0; i<iterations; i++) {
      fissionEvent* fe = new fissionEvent(isotope, time, nubar, energy_MeV, 1);
      int errorlength=maxerrorlength;
      fe->getFREYAerrors(&errorlength, &errors[0]);
      if (errorlength>1) {
         printf("%s\n",errors);
         exit(1);
      }
      int nneutrons = fe->getNeutronNu();
      for(int n1=0; n1<nneutrons; n1++) {
         double eng = fe->getNeutronEnergy(n1);
         int energy_bin_index = (int) (log(eng/minEnergy)/logf);
         if (energy_bin_index>=0 && energy_bin_index<nbinsperspectrum)
           hist[nneutrons][energy_bin_index]++;
      }
      delete fe;
   }
   output(f, hist);
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

void output(double f, int hist[][nbinsperspectrum]) {

   std::cout << "# neutrons\taverage energy [MeV]"
             << std::endl;
   for (int i=1; i<maxneutrons+1; i++) {
     char filename [1024];
     sprintf(filename, "nu_dist_n=%d.res", i);
     FILE* fp = openfile(filename);
     unsigned int sum=0;
     for (int j=0; j<nbinsperspectrum; j++) sum += hist[i][j];
     for (int j=0; j<nbinsperspectrum; j++) fprintf(fp, "%e - %e : %f\n", minEnergy*pow(f,j), minEnergy*pow(f,j+1), 1.*hist[i][j]/sum);
     double ave_energy = 0.;
     for (int j=0; j<nbinsperspectrum; j++) ave_energy += 1.*hist[i][j]/sum*minEnergy*pow(pow(f,j)*pow(f,j+1),0.5);
     std::cout << i << "\t\t" << ave_energy
               << std::endl;
     fclose(fp);
   }

   return;
}


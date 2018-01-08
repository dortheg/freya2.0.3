#define iterations 30000000
#define maxneutrons 15

// Mannhart energy structure for fission neutron spectrum histogram
#define nbinsperspectrum 71

#include <stdio.h>
#include <math.h>
#include "fissionEvent.h"

void init(void);
FILE* openfile(char* name);
void output(double* binuppererg, int hist[][nbinsperspectrum]);

int main() {
   bool spontaneousfission=true;
   int isotope = 98252;
   double energy_MeV = 1.00;
   double nubar = 3.00888;
   double time = 0.;

   double* binuppererg = new double [nbinsperspectrum];
   binuppererg[0]=0.015;
   binuppererg[1]=0.035;
   binuppererg[2]=0.055;
   binuppererg[3]=0.075;
   binuppererg[4]=0.095;
   binuppererg[5]=0.115;
   binuppererg[6]=0.135;
   binuppererg[7]=0.165;
   binuppererg[8]=0.195;
   binuppererg[9]=0.225;
   binuppererg[10]=0.255;
   binuppererg[11]=0.305;
   binuppererg[12]=0.355;
   binuppererg[13]=0.405;
   binuppererg[14]=0.455;
   binuppererg[15]=0.505;
   binuppererg[16]=0.555;
   binuppererg[17]=0.605;
   binuppererg[18]=0.655;
   binuppererg[19]=0.705;
   binuppererg[20]=0.755;
   binuppererg[21]=0.805;
   binuppererg[22]=0.855;
   binuppererg[23]=0.905;
   binuppererg[24]=0.955;
   binuppererg[25]=1.050;
   binuppererg[26]=1.150;
   binuppererg[27]=1.250;
   binuppererg[28]=1.350;
   binuppererg[29]=1.450;
   binuppererg[30]=1.550;
   binuppererg[31]=1.650;
   binuppererg[32]=1.750;
   binuppererg[33]=1.850;
   binuppererg[34]=1.950;
   binuppererg[35]=2.150;
   binuppererg[36]=2.350;
   binuppererg[37]=2.550;
   binuppererg[38]=2.750;
   binuppererg[39]=2.950;
   binuppererg[40]=3.250;
   binuppererg[41]=3.550;
   binuppererg[42]=3.850;
   binuppererg[43]=4.150;
   binuppererg[44]=4.450;
   binuppererg[45]=4.750;
   binuppererg[46]=5.050;
   binuppererg[47]=5.550;
   binuppererg[48]=6.050;
   binuppererg[49]=6.550;
   binuppererg[50]=7.050;
   binuppererg[51]=7.550;
   binuppererg[52]=8.050;
   binuppererg[53]=8.550;
   binuppererg[54]=9.050;
   binuppererg[55]=9.550;
   binuppererg[56]=10.050;
   binuppererg[57]=10.550;
   binuppererg[58]=11.050;
   binuppererg[59]=11.550;
   binuppererg[60]=12.050;
   binuppererg[61]=12.550;
   binuppererg[62]=13.050;
   binuppererg[63]=13.550;
   binuppererg[64]=14.050;
   binuppererg[65]=14.600;
   binuppererg[66]=15.900;
   binuppererg[67]=16.900;
   binuppererg[68]=17.900;
   binuppererg[69]=19.100;
   binuppererg[70]=20.000;

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   int hist[maxneutrons+1][nbinsperspectrum];
   for (int i=0; i<maxneutrons+1; i++)
      for (int j=0; j<nbinsperspectrum; j++)
         hist[i][j] = 0.;

   init();
   for (int i=0; i<iterations; i++) {
      fissionEvent* fe = new fissionEvent(isotope, time, nubar, energy_MeV, (spontaneousfission==true)?0:1);
      int errorlength=maxerrorlength;
      fe->getFREYAerrors(&errorlength, &errors[0]);
      if (errorlength>1) {
         printf("%s\n",errors);
         exit(1);
      }
      int nneutrons = fe->getNeutronNu();
      for(int n1=0; n1<nneutrons; n1++) {
         double eng = fe->getNeutronEnergy(n1);
         int energy_bin_index = 0;
         for (int i=0; i<nbinsperspectrum; i++) {
            if (eng < binuppererg[i]) {
               energy_bin_index=i;
               hist[nneutrons][energy_bin_index]++;
               break;
            }
         }
         if (energy_bin_index==0 && eng > binuppererg[nbinsperspectrum-1]) {
            printf("Neutron energy (%f MeV) out of bounds\n", eng);
         }
      }
      delete fe;
   }
   output(binuppererg, hist);
   delete binuppererg;
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

void output(double* binuppererg, int hist[][nbinsperspectrum]) {

   printf("# neutrons\taverage energy [MeV]\n");
   for (int i=1; i<maxneutrons+1; i++) {
     char filename [1024];
     sprintf(filename, "nu_dist_n=%d.res", i);
     FILE* fp = openfile(filename);
     int sum=0;
     for (int j=0; j<nbinsperspectrum; j++) sum += hist[i][j];
     double ave_energy = 0.;
     if (sum != 0) {
        for (int j=0; j<nbinsperspectrum; j++) {
          double E0=(j==0)?0:binuppererg[j-1];
          double E1=binuppererg[j];
          fprintf(fp, "%e - %e : %f\n", E0, E1, 1.*hist[i][j]);
          // fprintf(fp, "%e - %e : %f\n", E0, E1, 1.*hist[i][j]/sum/(E1-E0));
          ave_energy += 1.*hist[i][j]/sum*pow(E0*E1,0.5);
        }
     }
     printf("%i\t\t%f\n", i, ave_energy);
     fclose(fp);
   }

   return;
}


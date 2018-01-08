#define iterations 100000
#define nisotopes 6
#define nenergies 11
#define specmax 660
#define nmax 20
#define pmax 40
#define eng0 1.e-9
#define f 1.04

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Fission.h"

void init(int* nnu, int* pnu);
void initspec(int* nspec, int* pspec);
FILE* openfile(char* name);
void specoutput(FILE* fp, int isotope, double energy, int* nspec, int* pspec, int nneutrons, int nphotons);

int isotopes[nisotopes] = {
  90231,
  91233, 
  92235, 92238,
  94239,
  98252
};

double nubar[nisotopes] = {
  2.00, // placeholder
  2.00, // placeholder
  2.41, 2.3,
  2.88,
  2.00 // placeholder
}; // nubar is not really needed

double energies[nenergies] = {
  25.3e-9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

void init(int* nnu, int* pnu) {
   /* Initialization */
   unsigned short int s[3] = {1234, 5678, 9012};
   int zero = 0;
   int one = 1;
   int two = 2;
   int i;

   setcorrel_(&two);
   setnudist_(&one);
//   setcorrel_(&zero);
   seed48(s);
   for (i=0; i<nmax; i++) nnu[i] = 0.;
   for (i=0; i<pmax; i++) pnu[i] = 0.;
}

void initspec(int* nspec, int* pspec) {
   /* Initialization */
   int i;

   for (i=0; i<specmax; i++) nspec[i] = 0.;
   for (i=0; i<specmax; i++) pspec[i] = 0.;
}

FILE* openfile(char* name) {
   FILE* tmpfp;
   tmpfp = fopen(name, "w");
   if (tmpfp == (FILE *) 0) fprintf(stderr, "Could not open %s for writing", name);
   return tmpfp;
}

void output(FILE* fp, int isotope, int* nnu, int* pnu, double eng, int terrell) {
   int i;

   for (i=0; i<nmax; i++) {
      if (nnu[i] != 0) fprintf(fp, "%d n %e %d %d %e\n", isotope, eng, terrell, i, 1.*nnu[i]/iterations);
   }
   fprintf(fp, "\n\n");
   for (i=0; i<pmax; i++) {
      if (pnu[i] != 0) fprintf(fp, "%d p %e %d %d %e\n", isotope, eng, terrell, i, 1.*pnu[i]/iterations);
   }
   fprintf(fp, "\n\n");
}

void specoutput(FILE* fp, int isotope, double energy, int* nspec, int* pspec, int nneutrons, int nphotons) {
   int i;

   for (i=0; i<specmax; i++) {
      if (nspec[i] != 0 || pspec[i] != 0) 
         fprintf(fp, "%d %e %e %e %e\n", isotope, energy, eng0*pow((double) f, (double) i+.5), 1.*nspec[i]/nneutrons/(eng0*pow((double) f, (double) i)), 1.*pspec[i]/nphotons/(eng0*pow((double) f, (double) i)));
   }
   fprintf(fp, "\n\n");
   return;
}

int main ()
{
   FILE* fpnudist;
   FILE* fpspec;
   FILE* fpspecTotal;
   FILE* fpTotal;

   char nudistname [1024];
   char specname [1024];
   char specTotalname [1024];
   char Totalname [1024];

   int i, k, isotopeindex, engindex, nindex, pindex, index;
   double time = 0.;
   double eng; // MeV
   int nneutrons, nphotons;
   int nfissions;
   double nEng, pEng;
   int nspec[specmax], pspec[specmax];
   int nspecTotal[specmax], pspecTotal[specmax], specTotal[specmax];
   int nnu[nmax], pnu[pmax];

   sprintf(nudistname, "testEngConsAllActNudist.res");
   sprintf(specname, "testEngConsAllActSpec.res");
   sprintf(specTotalname, "testEngConsAllActSpecTotal.res");
   sprintf(Totalname, "testEngConsAllActTotal.res");
   fpnudist = openfile(nudistname);
   fpspec = openfile(specname);
   fpspecTotal = openfile(specTotalname);
   fpTotal = openfile(Totalname);

   for (isotopeindex=0; isotopeindex<nisotopes; isotopeindex++) {
      for (engindex=0; engindex<nenergies; engindex++) {
         eng = energies[engindex];
         init(nnu, pnu);
         initspec(nspec, pspec);
         initspec(nspecTotal, pspecTotal);
         initspec(nspecTotal, specTotal);
         nneutrons=0, nphotons = 0;
         nfissions=0;
         for (i=0; i<iterations; i++) {
            nEng=0.;
            pEng=0.;
            genfissevt_(&isotopes[isotopeindex], &time, &nubar[isotopeindex], &eng);
            nneutrons += getnnu_();
            nphotons += getpnu_();
            nfissions++;
            nnu[getnnu_()]++;
            pnu[getpnu_()]++;
            for(k=0; k<getnnu_(); k++) {
               nindex = (int) (log(getneng_(&k)/(double) eng0)/log((double) f));
               nindex = (nindex<0)?0:nindex;
               nspec[nindex]++;
               nEng += getneng_(&k);
            }
            if (0 != getnnu_()) {
               nindex = (int) (log(nEng/(double) eng0)/log((double) f));
               nindex = (nindex<0)?0:nindex;
               nspecTotal[nindex]++;
            }

            for(k=0; k<getpnu_(); k++) {
               pindex = (int) (log(getpeng_(&k)/(double) eng0)/log((double) f));
               pindex = (pindex<0)?0:pindex;
               pspec[pindex]++;
               pEng += getpeng_(&k);
            }
            if (0 != getpnu_()) {
               pindex = (int) (log(pEng/(double) eng0)/log((double) f));
               pindex = (pindex<0)?0:pindex;
               pspecTotal[pindex]++;
            }

            if (0 != getnnu_() && 0 != getpnu_()) {
               index = (int) (log((nEng+pEng)/(double) eng0)/log((double) f));
               index = (index<0)?0:index;
               specTotal[index]++;
            }
         }
         output(fpnudist, isotopes[isotopeindex], nnu, pnu, eng, 0);
         specoutput(fpspec, isotopes[isotopeindex], eng, nspec, pspec, nneutrons, nphotons);
         specoutput(fpspecTotal, isotopes[isotopeindex], eng, nspecTotal, pspecTotal, nfissions, nfissions);
         specoutput(fpTotal, isotopes[isotopeindex], eng, specTotal, specTotal, nfissions, nfissions);
      }
   }

   fclose(fpnudist);
   fclose(fpspec);
   fclose(fpspecTotal);
   fclose(fpTotal);
}

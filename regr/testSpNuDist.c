#define iterations 100000
#define nisotopes 19
#define nmax 20
#define pmax 40
#define nspecmax 400
#define pspecmax 400
#define neng0 1.e-8
#define peng0 1.e-8
#define f 1.1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Fission.h"

int isotopes[nisotopes] = {
  90232,
  92232, 92233, 92234, 92235, 92236, 92238,
  93237,
  94236, 94238, 94239, 94240, 94241, 94242,
  95241,
  96242, 96244,
  97249,
  98252
};

FILE *fp;
FILE *fpspec;

void init(int* nnu, int* pnu) {
   /* Initialization */
   int i;

   unsigned short int s[3] = {1234, 5678, 9012};

   seed48(s);
   for (i=0; i<nmax; i++) nnu[i] = 0.;
   for (i=0; i<pmax; i++) pnu[i] = 0.;
}

void initspec(int* nspec, int* pspec) {
   /* Initialization */
   int i;

   for (i=0; i<nspecmax; i++) nspec[i] = 0.;
   for (i=0; i<pspecmax; i++) pspec[i] = 0.;
}

void output(int isotope, int* nnu, int* pnu) {
   int i;

   if (fp == 0) {
      fp = fopen("testSpNuDist.res", "w");
      if (fp == (FILE *) 0) fprintf(stderr, "Could not open testSpNuDist.res for writing");
   }
   for (i=0; i<nmax; i++) {
      if (nnu[i] != 0) fprintf(fp, "%d n %d %e\n", isotope, i, 1.*nnu[i]/iterations);
   }
   fprintf(fp, "\n\n");
   for (i=0; i<pmax; i++) {
      if (pnu[i] != 0) fprintf(fp, "%d p %d %e\n", isotope, i, 1.*pnu[i]/iterations);
   }
   fprintf(fp, "\n\n");
}

void specoutput(int isotope, int* nspec, int* pspec, int nneutrons, int nphotons) {
   int i;

   if (fpspec == 0) {
      fpspec = fopen("testSpspec.res", "w");
      if (fpspec == (FILE *) 0) fprintf(stderr, "Could not open testSpspec.res for writing");
   }
   for (i=0; i<nspecmax; i++) {
      if (nspec[i] != 0) fprintf(fpspec, "%d n %e %e\n", isotope, neng0*pow((double) f, (double) i+.5), 1.*nspec[i]/nneutrons/(neng0*pow((double) f, (double) i)));
   }
   fprintf(fpspec, "\n\n");
   for (i=0; i<pspecmax; i++) {
      if (pspec[i] != 0) fprintf(fpspec, "%d p %e %e\n", isotope, peng0*pow((double) f, (double) i+.5), 1.*pspec[i]/nphotons/(peng0*pow((double) f, (double) i)));
   }
   fprintf(fpspec, "\n\n");
}

int main ()
{
   int i, k, isotopeindex, nindex, pindex;
   int zero = 0, one = 1;

   double time = 0.;
   int nnu[nmax], pnu[pmax];
   int nneutrons, nphotons;
   int nspec[nspecmax], pspec[pspecmax];

   for (isotopeindex=0; isotopeindex<nisotopes; isotopeindex++) {
      init(nnu, pnu);
      initspec(nspec, pspec);
      nneutrons=0, nphotons = 0;
      for (i=0; i<iterations; i++) {
         genspfissevt_(&isotopes[isotopeindex], &time);
         nneutrons += getnnu_();
         nphotons += getpnu_();
         nnu[getnnu_()]++;
         pnu[getpnu_()]++;
         for(k=0; k<getnnu_(); k++) {
            nindex = (int) (log(getneng_(&k)/(double) neng0)/log((double) f));
            nindex = (nindex<0)?0:nindex;
            nspec[nindex]++;
         }
         for(k=0; k<getpnu_(); k++) {
            pindex = (int) (log(getpeng_(&k)/(double) peng0)/log((double) f));
            pindex = (pindex<0)?0:pindex;
            pspec[pindex]++;
         }
      }
      output(isotopes[isotopeindex], nnu, pnu);
      specoutput(isotopes[isotopeindex], nspec, pspec, nneutrons, nphotons);
   }

   init(nnu, pnu);
   setcf252_(&one, &zero);
   for (i=0; i<iterations; i++) {
      genspfissevt_(&isotopes[nisotopes-1], &time);
      nnu[getnnu_()]++;
      pnu[getpnu_()]++;
   }
   output(isotopes[nisotopes-1], nnu, pnu);
   fclose(fp);
   fclose(fpspec);
}

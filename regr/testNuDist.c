#define iterations 30000
#define nisotopes 10
#define nmax 20
#define pmax 40
#define nspecmax 400
#define pspecmax 400
#define maxneng 10.
#define maxpeng 10.
#define neng0 1.e-8
#define peng0 1.e-8
#define f 1.1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Fission.h"

int isotopes[nisotopes] = {
 92232,
 92233,
 92234,
 92235,
 92236,
 92238,
 94239,
 94241,
 94240,
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

void output(int isotope, int* nnu, int* pnu, double eng, int terrell) {
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

void specoutput(int isotope, int* nspec, int* pspec, double eng, int nneutrons, int nphotons) {
   int i;

   for (i=0; i<nspecmax; i++) {
      if (nspec[i] != 0) fprintf(fpspec, "%d n %e %e %e\n", isotope, eng, neng0*pow((double) f, (double) i+.5), 1.*nspec[i]/nneutrons/(neng0*pow((double) f, (double) i)));
   }
   fprintf(fpspec, "\n\n");
   for (i=0; i<pspecmax; i++) {
      if (pspec[i] != 0) fprintf(fpspec, "%d p %e %e %e\n", isotope, eng, peng0*pow((double) f, (double) i+.5), 1.*pspec[i]/nphotons/(peng0*pow((double) f, (double) i)));
   }
   fprintf(fpspec, "\n\n");
}

int main ()
{
   int i, j, k, isotopeindex, nindex, pindex;
   int number[4] = {0, 1, 2, 3};
   double nubar = 2.4305631;
   double eng;

   double time = 0.;
   int nnu[nmax], pnu[pmax];
   int nneutrons, nphotons;
   int nspec[nspecmax], pspec[pspecmax];

   if (fp == 0) {
      fp = fopen("testNuDist.res", "w");
      if (fp == (FILE *) 0) fprintf(stderr, "Could not open testNuDist.res for writing");
   }
   if (fpspec == 0) {
      fpspec = fopen("testspec.res", "w");
      if (fpspec == (FILE *) 0) fprintf(stderr, "Could not open testspec.res for writing");
   }
   for (i=0; i<4; i++) {
      setnudist_(&number[i]);
      for (isotopeindex=0; isotopeindex<nisotopes; isotopeindex++) {
         eng = 0.;
         if (isotopes[isotopeindex] == 92232 ||
             isotopes[isotopeindex] == 92235 ||
             isotopes[isotopeindex] == 92238 ||
             isotopes[isotopeindex] == 94239) {
            for (eng = 0.; eng <=10.; eng=eng+1.) {
               fprintf(fp, "# energy : %d MeV\n", (int) eng);
               if (i == 0) fprintf(fpspec, "# energy : %d MeV\n", (int) eng);
               init(nnu, pnu);
               initspec(nspec, pspec);
               nneutrons=0, nphotons = 0;
               for (j=0; j<iterations; j++) {
                  genfissevt_(&isotopes[isotopeindex], &time, &nubar, &eng);
                  nneutrons += getnnu_();
                  nphotons += getpnu_();
                  nnu[getnnu_()]++;
                  pnu[getpnu_()]++;
                  if (i == 0) {
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
               }
               output(isotopes[isotopeindex], nnu, pnu, eng, number[i]);
               if (i == 0) specoutput(isotopes[isotopeindex], nspec, pspec, eng, nneutrons, nphotons);
            }
         } else {
            init(nnu, pnu);
            initspec(nspec, pspec);
            nneutrons=0, nphotons = 0;
            for (j=0; j<iterations; j++) {
               genfissevt_(&isotopes[isotopeindex], &time, &nubar, &eng);
               nneutrons += getnnu_();
               nphotons += getpnu_();
               nnu[getnnu_()]++;
               pnu[getpnu_()]++;
               if (i == 0) {
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
            }
            output(isotopes[isotopeindex], nnu, pnu, 0., number[i]);
            if (i == 0) specoutput(isotopes[isotopeindex], nspec, pspec, 0., nneutrons, nphotons);
         }
      }
   }

   fclose(fp);
   fclose(fpspec);
}

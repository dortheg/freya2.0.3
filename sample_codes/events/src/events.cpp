#define maxA 200
#define mMax 50       /* The maximum number of ejectiles per fission generated by FREYA */

#include <stdio.h>
#include <iostream>
#include <math.h>
#include "stdlib.h"
#include "Fission.h"

using namespace std;

extern "C" {
   extern int msfreya_setup_c_();
   extern int msfreya_event_c_(int,double,double,double*,int*,int*,double*,int*,int*,double*,int*,double*,int*,double*);
   extern int msfreya_getids_c_(int*,int*,int*);
   extern int msfreya_getniso_c_(int *,int *);
   extern int msfreya_getzas_c_(int *,int *);
   extern double msfreya_sepn_c_(int, int, int);
   extern double msfreya_gsmassn_c_(int, int);
   extern int msfreya_geterrors_c_(char *, int *);
   extern int msfreya_reseterrorflag_c_();
   extern int msfreya_errorflagset_c_();
   extern int msfreya_usehostrng_c_();
}

void init(void);
void initFREYA(int& nisosf, int& nisoif, int& niso,
               int** ZAs, int** fistypes);
bool FREYA_event(FILE* fp, int Z, int A, int fissionindex, double ePart, 
                 int fissiontype, int*& ZAs, int*& fistypes, int niso
                );
FILE* openfile(char* name);
void output_compound(FILE* fp, int Z, int A, double energy_MeV, int niterations);
void output_ff(FILE* fp, int fissionindex, int Z, int A, double exc_erg,
               int nmultff1, int gmultff1, double PP [5]);
void output_secondaries(FILE* fp, int ptypes [mMax], double particles [4*3*mMax], 
                        int npart2skip);
void readinput(int& Z, int& A, double& E, int& fissiontype, int& iterations, char outputfilename [1024]);

int main() {
   int iterations=10000;        // Number of fission events to be generated
   double energy_MeV = 25.3e-9; // thermal
   int Z = 94;
   int A = 239;
   char outputfilename [1024];
   sprintf(outputfilename, "history.res");
   int fissiontype = 1; // 0: spontaneous fission
                        // 1: neutron-induced fission

   int nisosf = 0; // Number of spontaneous fission isotopes
   int nisoif = 0; // Number of induced fission isotopes
   int niso = 0;   // Number of fission isotopes

   int** ZAs;      // ZA's of fission isotopes
   int** fistypes; // types of fission [spontaneous (0), induced (1)
   ZAs = new int*;
   fistypes = new int*;

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   init();
   initFREYA(nisosf, nisoif, niso, ZAs, fistypes);
   niso=nisosf+nisoif;

   // read in Z, A, energy_MeV, fission type, number of iterations, output file name
   readinput(Z, A, energy_MeV, fissiontype, iterations, outputfilename);

   FILE* fp = openfile(outputfilename);
   // fprintf(fp, "This is the output record from %g MeV n + %d%d -> f (%d events):\n\n", energy_MeV, Z, A, iterations);
   output_compound(fp, Z, A+1, (fissiontype==0)?0.:energy_MeV, iterations);
   
   for (int i=0; i<iterations; i++) {
      if (!FREYA_event(fp, Z, A, i, energy_MeV, fissiontype, *ZAs, *fistypes, niso)) {
         int errorlength=maxerrorlength;
         msfreya_geterrors_c_(&errors[0], &errorlength);
         if (errorlength>1) {
            printf("%s\n",errors);
            exit(1);
         }
      }
   }
   fprintf(fp, "    0    0    0\n");
   fclose(fp);
}

FILE* openfile(char* name) {
   FILE* fp = fopen(name, "w");
   if (fp == (FILE *) 0) fprintf(stderr, "Could not open %s for writing\n", name);
   return fp;
}

void init(void) {
   unsigned short int s[3] = {1234, 5678, 9012};
   int i;
   seed48(s);
   // msfreya_usehostrng_c_();
   return;
}

void initFREYA(int& nisosf, int& nisoif, int& niso,
               int** ZAs, int** fistypes) {

   int maxerrorlength=10000;
   char errors[maxerrorlength];

   msfreya_reseterrorflag_c_();
   msfreya_setup_c_();
   if (msfreya_errorflagset_c_()==1) {
      int errorlength=maxerrorlength;
      msfreya_geterrors_c_(&errors[0], &errorlength);
      if (errorlength>1) {
         printf("%s\n",errors);
         exit(1);
      }
   }
   msfreya_getniso_c_(&nisosf, &nisoif);
   niso=nisosf+nisoif;

   // allocate memory to store ZA's for spontaneous and neutron-induced
   // fissions
   *ZAs = new int [niso];
   *fistypes = new int [niso];

   // Populate ZAs and fistypes
   msfreya_getzas_c_(&(*ZAs[0]),&(*fistypes[0]));
}

bool FREYA_event(FILE* fp, int Z, int A, int fissionindex, double ePart, 
                 int fissiontype, int*& ZAs, int*& fistypes, int niso
                ) {
   int isotope = 1000*Z+A;
   // if the compound nucleus is ZA, the original nucleus was
   //   ZA for photofission
   //   Z(A-1) for neutron-induced fission
   // treat photofission as if it were neutron-induced fission
   if (fissiontype==2) isotope--;
   
   // Find the index of the fission/isotope
   bool foundfission=false;
   int iKm1=0;
   for (iKm1=0; iKm1<niso; iKm1++)
      if (isotope == ZAs[iKm1] && ((fissiontype==0) == (fistypes[iKm1]==0))) {
         foundfission=true;
         break;
      }
   if (!foundfission) {
      fprintf(stderr, "ABORT: fission type %d not supported for isotope %d\n", fissiontype, isotope);
      exit(1);
   }

   int iK=iKm1+1; // FORTRAN indexing
   int freyaA=isotope-1000*Z;
   // watch out! in freya, the A for induced fission is the A of the 
   // compound nucleus (for induced fission, add 1 neutron to the nucleus)
   freyaA+=(fissiontype==0)?0:1;
   msfreya_reseterrorflag_c_();

   // Compute nucleus excitation energy for this event
   double eps0;
   double En;
   switch (fissiontype) {
      case 0:
         // spontaneous fission
         eps0 = 0.;
         En=0.;
         break;
      case 1:
         // neutron-induced fission
      case 2:
         // photon-induced fission
         double sepni;
         sepni = msfreya_sepn_c_(iK,Z,freyaA);
         if (msfreya_errorflagset_c_()==1) return false;

         if (fissiontype==1) {
            // neutron-induced fission
            eps0 = sepni+ePart;
            En=ePart;
         } else if (fissiontype==2) {
            // photon-induced fission
            eps0 = ePart;
            En=ePart-sepni;
            if (En<0) En=0.;
         }
         break;
      default:
         fprintf(stderr, "ABORT: fission type %d not supported\n", fissiontype);
         exit(1);
         break;
   }

   // ...generate fission event
   // declare those, msfreya_event needs them
   double V0[3]; // velocity of the initial nucleus
   for (int i=0; i<3; i++) V0[i]=0; // nucleus at rest

   double P0[5]; // excited energy, momentum and kinetic energy
                // of nucleus before interaction
   double P1[5]; // excited energy, momentum and kinetic energy
                // of fission fragment 1
   double P2[5]; // excited energy, momentum and kinetic energy
                // of fission fragment 2
   int Z1, A1;  // Charge & mass number of fission fragment 1
   int Z2, A2;  // Charge & mass number of fission fragment 2

   double W0=msfreya_gsmassn_c_(Z, freyaA);  // ground-state mass of nucleus
   if (msfreya_errorflagset_c_()==1) return false;
   
   double ndir [3];              // neutron direction ((0,0,0) forces isotropic)
   for (int i=0; i<3; i++)
      ndir[i]=0.;

   P0[0]=W0+eps0;                // Rest energy of init nucleus
   double g0=1.0;                // gamma0
   P0[4]=g0*P0[0];               // Total energy of init nucleus
   for (int i=0; i<3; i++)
      P0[i+1]=P0[4]*V0[i];       // Momentum of initial nucleus

   int mult;                     // Number of particles emitted
   int nmultff0,gmultff0;        // Number of pre-fission neutrons/photons emitted
   int nmultff1,nmultff2;        // Number of neutrons emitted by fission fragments
   int gmultff1,gmultff2;        // Number of photons emitted by fission fragments
   double particles [4*3*mMax];  // their momentum and kinetic energy
   int ptypes [3*mMax];          // their type: 0(g) & 1(n)
   int ptypes0 [mMax];           // pre-fission ejectile types
   int ptypes1 [mMax];           // types of 1st fission fragment ejectiles
   int ptypes2 [mMax];           // types of 2nd fission fragment ejectiles
   
   msfreya_event_c_(iK,En,eps0,&(P0[0]),&Z1,&A1,&(P1[0]),&Z2,&A2,&(P2[0]),&mult,&(particles[0]),&(ptypes[0]),&(ndir[0]));
   if (msfreya_errorflagset_c_()==1) return false;

   msfreya_getids_c_(&(ptypes0[0]),&(ptypes1[0]),&(ptypes2[0]));

   // count number of pre-fission neutrons/photons emitted
   nmultff0=0;
   gmultff0=0;
   // total number of pre-ejectiles
   int npart0=0;
   for (int i=0; i<mMax; i++) {
     if (0 == ptypes0[i]) {
       gmultff0=gmultff0+1;
     } else if (1 == ptypes0[i]) {
       nmultff0=nmultff0+1;
     } else if (-1 == ptypes0[i]) {
       npart0=i;
       break;
     }
   }

   // count number of neutrons/photons for fission fragment 1
   nmultff1=0;
   gmultff1=0;
   // total number of secondaries from fission fragment 1
   int npart1=0;
   for (int i=0; i<mMax; i++) {
     if (0 == ptypes1[i]) {
       gmultff1=gmultff1+1;
     } else if (1 == ptypes1[i]) {
       nmultff1=nmultff1+1;
     } else if (-1 == ptypes1[i]) {
       npart1=i;
       break;
     }
   }

   // count number of neutrons/photons for fission fragment 2
   nmultff2=0;
   gmultff2=0;
   // total number of secondaries from fission fragment 2
   int npart2=0;
   for (int i=0; i<mMax; i++) {
     if (0 == ptypes2[i]) {
       gmultff2=gmultff2+1;
     } else if (1 == ptypes2[i]) {
       nmultff2=nmultff2+1;
     } else if (-1 == ptypes2[i]) {
       npart2=i;
       break;
     }
   }
   //....print results for pre-fission neutrons
   output_ff(fp, fissionindex+1, Z, A+1-nmultff0, eps0, nmultff0, gmultff0, P0);
   output_secondaries(fp, ptypes0, particles, 0);

   //....print results for fission fragment #1
   double W0_ff1=msfreya_gsmassn_c_(Z1, A1-nmultff1);  // ground-state mass of ff #1
   if (msfreya_errorflagset_c_()==1) return false;
   output_ff(fp, fissionindex+1, Z1, A1, P1[0]-W0_ff1, nmultff1, gmultff1, P1);
   output_secondaries(fp, ptypes1, particles, npart0);

   //....print results for fission fragment #2
   double W0_ff2=msfreya_gsmassn_c_(Z2, A2-nmultff2);  // ground-state mass of ff #2
   if (msfreya_errorflagset_c_()==1) return false;
   output_ff(fp, fissionindex+1, Z2, A2, P2[0]-W0_ff2, nmultff2, gmultff2, P2);
   output_secondaries(fp, ptypes2, particles, npart0+npart1);

   return true;
}

void output_compound(FILE* fp, int Z, int A, double energy_MeV, int niterations) {
   fprintf(fp, "%5d%5d%10.3f:%8d events\n", Z, A, energy_MeV, niterations);
   return;
}

void output_ff(FILE* fp, int fissionindex, int Z, int A, double exc_erg,
               int nmultff, int gmultff, double PP [5]) {
   double s = pow(PP[1],2)+pow(PP[2],2)+pow(PP[3],2);
   double ke_ff=0.5*s/PP[4];

   fprintf(fp, "%8d%5d%5d%10.3f%5d%5d\n", fissionindex, Z, A, exc_erg, nmultff, gmultff);
   if (0 == s) {
     fprintf(fp, "%7.3f%10.3f%10.3f%10.3f\n", ke_ff, 0., 0., 0.);
   } else {
     s = sqrt(s);
     fprintf(fp, "%7.3f%10.3f%10.3f%10.3f\n", ke_ff, PP[1]/s, PP[2]/s, PP[3]/s);
   }
   return;
}

void output_secondaries(FILE* fp, int ptypes [mMax], double particles [4*3*mMax],
                        int npart2skip) {
   const double wn=939.5651828; // neutron mass in MeV
   int count;
   double u, v, w, s, ke;
   //....print u,v,w,energy for neutrons
   count=0;
   for (int i=npart2skip; i<mMax; i++) {
     if (1 == ptypes[i-npart2skip]) {
       u=particles[i*4];
       v=particles[i*4+1];
       w=particles[i*4+2];
       s=pow(u,2)+pow(v,2)+pow(w,2);
       ke=0.5*s/wn;
       s=sqrt(s);
       fprintf(fp, "%7.3f%7.3f%7.3f%7.3f ", ke, u/s, v/s, w/s);
       count++;
     } else if (-1 == ptypes[i]) break;
   }
   if (count>0) fprintf(fp, "\n");

   //....print u,v,w,energy for photons
   count=0;
   for (int i=npart2skip; i<mMax; i++) {
     if (0 == ptypes[i-npart2skip]) {
       u=particles[i*4];
       v=particles[i*4+1];
       w=particles[i*4+2];
       s=pow(u,2)+pow(v,2)+pow(w,2);
       s=sqrt(s);
       fprintf(fp, "%7.3f%7.3f%7.3f%7.3f ", s, u/s, v/s, w/s);
       count++;
     } else if (-1 == ptypes[i-npart2skip]) break;
   }
   if (count>0) fprintf(fp, "\n");
   return;
}

void readinput(int& Z, int& A, double& E, int& fissiontype, int& iterations, char outputfilename [1024]) {
  cout << "Value of Z: ";
  cin >> Z;
  cout << "Value of A: ";
  cin >> A;
  cout << "Value of energy (in units of MeV, -1 for spontaneous fission): ";
  cin >> E;
  cout << "Number of iterations: ";
  cin >> iterations;
  cout << "Output file name: ";
  cin >> outputfilename;
  fissiontype=1;
  if (-1==E) fissiontype=0;
  if (0==fissiontype)
    cout << iterations << " spontaneous fissions of " << Z << A << endl;
  else
    cout << iterations << " neutron-induced fissions of " << Z << A << " (E=" << E << " MeV)" << endl;
  cout << "output file name: " << outputfilename << endl;
}

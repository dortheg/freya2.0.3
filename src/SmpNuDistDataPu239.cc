/* 
Copyright (c) 2006-2016 Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory 
UCRL-CODE-224807.

All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

o Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.

o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the documentation and/or other materials provided with the distribution.

o Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE. 

2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights. 

3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
*/


#include <math.h>
#include "fissionEvent.h"

int fissionEvent::SmpNuDistDataPu239(double erg) {

/*
  Description
    Sample Number of Neutrons from fission in Pu-239 using 
    Zucker and Holden's tabulated data for Pu-239
*/

/*
  Input
    erg      - incident neutron energy
  Output
    SmpNuDistDataPu239  - sampled multiplicity
    
*/
 
  double pnu[8]={0., 0., 0., 0., 0., 0., 0., 0.};
  double cpnu;
  double eng;
  double r;

/* 
  Check if energy is within the range of experimental values
*/
  if (erg > 10) eng=10.;
  else eng=erg;

  r=fisslibrng();
/*
  Pu-239 nu distribution
*/
  if (eng <= 5.0) pnu[0]=0.0108826e0-0.00207694e0*eng-6.5e-4*pow(eng,2)+4.023e-4*pow(eng,3)-7.93e-5*pow(eng,4)+5.53666667e-6*pow(eng,5);
       
  if (eng > 5 && eng <= 10) pnu[0]=0.078606e0-5.17531e-2*eng+1.42034e-2*pow(eng,2)-1.96292e-3*pow(eng,3)+1.34512e-4*pow(eng,4)-3.63416e-6*pow(eng,5);
  if (r <= pnu[0]) return 0;

  if (eng <= 5.0) pnu[1]=0.0994916e0-0.01979542e0*eng-0.00236583e0*pow(eng,2)+0.0020581e0*pow(eng,3)-4.14016667e-4*pow(eng,4)+2.85666667e-5*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[1]=0.10052e0-2.61361e-2*eng+3.78355e-3*pow(eng,2)-3.70667e-4*pow(eng,3)+1.95458e-5*pow(eng,4)-3.87499e-7*pow(eng,5);
  cpnu=pnu[0]+pnu[1];
  if (r <= cpnu) return 1;

  if (eng <= 5.0) pnu[2]=0.2748898e0-0.01565248e0*eng-0.00749681e0*pow(eng,2)+0.00217121e0*pow(eng,3)-3.13041667e-4*pow(eng,4)+1.88183333e-5*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[2]=0.282487e0-0.0261342e0*eng-1.16895e-3*pow(eng,2)+1.9888e-4*pow(eng,3)-6.41257e-6*pow(eng,4)+1.02502e-7*pow(eng,5);
  cpnu=cpnu+pnu[2];
  if (r <= cpnu) return 2;

  if (eng <= 5.0) pnu[3]=0.3269196e0+0.00428312e0*eng-0.00189322e0*pow(eng,2)-4.31925001e-4*pow(eng,3)+1.18466667e-4*pow(eng,4)-9.04166668e-6*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[3]=0.329058e0+4.023e-3*eng-3.06402e-3*pow(eng,2)+2.2628e-4*pow(eng,3)-1.50875e-5*pow(eng,4)+4.39168e-7*pow(eng,5);
  cpnu=cpnu+pnu[3];
  if (r <= cpnu) return 3;

  if (eng <= 5.0) pnu[4]=0.2046061e0+0.02633899e0*eng+0.0041514e0*pow(eng,2)-0.00275542e0*pow(eng,3)+5.0325e-4*pow(eng,4)-3.32158333e-5*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[4]=0.18992e0+4.55188e-2*eng-7.06316e-3*pow(eng,2)+7.29916e-4*pow(eng,3)-4.71791e-5*pow(eng,4)+1.185e-6*pow(eng,5);
  cpnu=cpnu+pnu[4];
  if (r <= cpnu) return 4;

  if (eng <= 5.0) pnu[5]=0.0726834e0+0.00116043e0*eng+0.007572e0*pow(eng,2)-0.00161972e0*pow(eng,3)+2.3545e-4*pow(eng,4)-1.546e-5*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[5]=0.0779212e0-1.35849e-3*eng+6.68583e-3*pow(eng,2)-7.98649e-4*pow(eng,3)+4.88625e-5*pow(eng,4)-1.54167e-6*pow(eng,5);
  cpnu=cpnu+pnu[5];
  if (r <= cpnu) return 5;

  if (eng <= 5.0) pnu[6]=0.0097282e0+0.00494589e0*eng+0.00115294e0*pow(eng,2)-3.25191667e-4*pow(eng,3)+6.00083333e-5*pow(eng,4)-3.745e-6*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[6]=7.85432e-3+7.33182e-3*eng-2.03705e-4*pow(eng,2)+8.73787e-5*pow(eng,3)-4.24164e-6*pow(eng,4)+2.37499e-7*pow(eng,5);
  cpnu=cpnu+pnu[6];
  if (r <= cpnu) return 6;

  if (eng <= 5.0) pnu[7]=6.301e-4+1.10666667e-4*eng+4.28016667e-4*pow(eng,2)+1.12041667e-5*pow(eng,3)-4.31666667e-6*pow(eng,4)+3.29166667e-7*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[7]=1.5323e-3-7.91857e-4*eng+8.01017e-4*pow(eng,2)-6.82833e-5*pow(eng,3)+4.38333e-6*pow(eng,4)-6.0e-8*pow(eng,5);
  cpnu=cpnu+pnu[7];
  if (r <= cpnu) return 7;
  else return 8;
}

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

int fissionEvent::SmpNuDistDataU238(double erg) {

/*
  Description
    Sample Number of Neutrons from fission in U-238 using 
      Zucker and Holden's tabulated data for U-238
*/

/*
  Input
    erg      - incident neutron energy
  Output
    SmpNuDistDataU238  - sampled multiplicity
    
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
  U-238 nu distribution
*/
  if (eng <= 5.0) pnu[0]=0.0396484e0-1.14202e-2*eng+1.94627e-3*pow(eng,2)-2.95412e-4*pow(eng,3)+2.98333e-5*pow(eng,4)-1.31417e-6*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[0]=0.0360013e0-8.06662e-3*eng+6.59461e-4*pow(eng,2)-3.54123e-5*pow(eng,3)+2.03749e-6*pow(eng,4)-5.91663e-8*pow(eng,5);
  if (r <= pnu[0]) return 0;

  if (eng <= 5.0) pnu[1]=0.252954e0-5.17151e-2*eng+2.84558e-3*pow(eng,2)+2.93563e-4*pow(eng,3)-5.99833e-5*pow(eng,4)+3.34417e-6*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[1]=0.259093e0-5.7936e-2*eng+5.50734e-3*pow(eng,2)-3.09607e-4*pow(eng,3)+1.20957e-5*pow(eng,4)-2.49997e-7*pow(eng,5);
  cpnu=pnu[0]+pnu[1];
  if (r <= cpnu) return 1;

  pnu[2]=0.29395353e0+0.01098908e0*eng-0.00565976e0*pow(eng,2)+3.14515399e-4*pow(eng,3)-5.66793415e-6*pow(eng,4)+1.54070513e-7*pow(eng,5);
  cpnu=cpnu+pnu[2];
  if (r <= cpnu) return 2;

  if (eng <= 5.0) pnu[3]=0.264447e0+3.02825e-2*eng-3.12762e-3*pow(eng,2)-1.5875e-4*pow(eng,3)+4.91667e-5*pow(eng,4)-3.38667e-6*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[3]=0.262024e0+3.37134e-2*eng-5.01711e-3*pow(eng,2)+3.58761e-4*pow(eng,3)-2.17959e-5*pow(eng,4)+5.10834e-7*pow(eng,5);
  cpnu=cpnu+pnu[3];
  if (r <= cpnu) return 3;

  if (eng <= 5) pnu[4]=0.111176e0+1.66321e-2*eng+2.56307e-3*pow(eng,2)-2.17754e-4*pow(eng,3)-5.96667e-6*pow(eng,4)+7.44167e-7*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[4]=0.107859e0+1.88862e-2*eng+2.07521e-3*pow(eng,2)-2.08099e-4*pow(eng,3)+3.23745e-6*pow(eng,4)-1.24999e-7*pow(eng,5);
  cpnu=cpnu+pnu[4];
  if (r <= cpnu) return 4;

  if (eng <= 5.0) pnu[5]=0.0312261e0+4.12932e-3*eng+9.18413e-4*pow(eng,2)+4.36542e-5*pow(eng,3)-5.9125e-6*pow(eng,4)+3.20833e-7*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[5]=0.0324008e0+3.04772e-3*eng+1.3327e-3*pow(eng,2)-3.96916e-5*pow(eng,3)+2.94583e-6*pow(eng,4)-7.66666e-8*pow(eng,5);
  cpnu=cpnu+pnu[5];
  if (r <= cpnu) return 5;

  if (eng <= 5.0) pnu[6]=5.9347e-3+9.80023e-4*eng+4.24667e-4*pow(eng,2)+3.04458e-5*pow(eng,3)-6.46667e-6*pow(eng,4)+4.30833e-7*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[6]=6.5421e-3+3.07834e-4*eng+7.39454e-4*pow(eng,2)-4.70459e-5*pow(eng,3)+3.44583e-6*pow(eng,4)-8.91667e-8*pow(eng,5);
  cpnu=cpnu+pnu[6];
  if (r <= cpnu) return 6;

  if (eng <= 5.0) pnu[7]=5.436e-4+1.3756e-4*eng-5.0e-7*pow(eng,2)+1.35917e-5*pow(eng,3)+5.0e-7*pow(eng,4)-5.16667e-8*pow(eng,5);
  if (eng > 5 && eng <= 10) pnu[7]=9.212e-4-1.57585e-4*eng+8.41126e-5*pow(eng,2)+4.14166e-6*pow(eng,3)+5.37501e-7*pow(eng,4)-6.66668e-9*pow(eng,5);
  cpnu=cpnu+pnu[7];
  if (r <= cpnu) return 7;
  else return 8;
}

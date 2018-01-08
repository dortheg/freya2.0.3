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


#define nZAEngConsAllActP 12 /* 12 isotopes for photons in same paper */

#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "fissionEvent.h"

double fissionEvent::getTotEngPEnergyConsAllActinides(double ePart, int iso) {

/*
  Description
    Determine the average total energy of prompt photons emitted by fission.

    The data in this function comes from the paper "Energy-Dependent 
    Fission Q Values Generalized for All Actinides" by R. Vogt. The paper 
    gives the total energy for all emitted prompt fission photons. Data 
    in the paper is given for all actinides, major and minor, in the 
    Evaluated Nuclear Data Library 2008 release, ENDL2008.
    The average total prompt fission photon energy is given as a 
    function of the incident neutron energy by the formula:
      PEng = coeffPhoton[2]+coeffPhoton[1]*ePart+coeffPhoton[0]*ePart^2
*/

/*
  Input
    ePart     - energy of incoming particle
    iso       - isotope
  Output
    getTotEngPEnergyConsAllActinides - average total energy of prompt fission
                                       gamma-rays
*/

   static int ZAEngConsAllActP [nZAEngConsAllActP]= {
                   92232, // U-232
                   92233, // U-233
                   92234, // U-234
                   92235, // U-235
                   92236, // U-236
                   92237, // U-237
                   92238, // U-238
                   92240, // U-240
                   92241, // U-241
                   94239, // Pu-239
                   98252, // Cf-252
                   0      // generic
                            };

   static double coeffPhoton [nZAEngConsAllActP][3] = {
                  {0.000182, 0.0255,  7.256},   // U-232
                  {0.000182, 0.0255,  7.256},   // U-233
                  {0.000182, 0.0255,  7.256},   // U-234
                  {-0.00474, 0.2295,  7.284},   // U-235
                  {0.000182, 0.0255,  7.256},   // U-236
                  {0.000182, 0.0255,  7.256},   // U-237
                  {-1.22e-7, 0.01607, 6.658},   // U-238
                  {0.000182, 0.0255,  7.256},   // U-240
                  {0.000182, 0.0255,  7.256},   // U-241
                  {-0.009878,0.4249,  6.857},   // Pu-239
                  {0.,       0.01831, 6.44186}, // Cf-252
                  {7.238e-8, 0.01693, 6.95}     // generic, all others
                            };
   double muEp;

   int i;

// Find photon parameters for isotope
   int isoindexp=-1;
   for (i=0; isoindexp == -1 && i<nZAEngConsAllActP-1; i++) {
      if (iso == ZAEngConsAllActP[i]) {
         isoindexp = i;
         break;
      }
   }
   if (isoindexp == -1) isoindexp = nZAEngConsAllActP-1; // revert to generic if nothing else
   
   muEp = coeffPhoton[isoindexp][2] + ePart*(coeffPhoton[isoindexp][1] + ePart*coeffPhoton[isoindexp][0]);

   return muEp;
}

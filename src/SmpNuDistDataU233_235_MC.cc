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

int fissionEvent::SmpNuDistDataU233_235_MC(double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in U-233 and U-235 using 
    (a) Gwin, Spencer and Ingle tabulated data at thermal 
        energies (0 MeV),
    (b) Zucker and Holden's tabulated data for U-235 at 1 MeV and 
        higher.
    The 11 P(nu) distributions are given as a function of nubar, 
    the average number of neutrons from induced fission for the 
    11 different energies (0 to 10 MeV), based on the U-235 data 
    above.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    SmpNuDistDataU233_235_MC  - sampled multiplicity
    
*/

  static double U235nu [11] [8] = {
     {.0291000, .1660000, .3362000, .3074000, .1333000, .0259000, .0021000, .0002000},
     {.0237898, .1555525, .3216515, .3150433, .1444732, .0356013, .0034339, .0004546},
     {.0183989, .1384891, .3062123, .3217566, .1628673, .0455972, .0055694, .0011093},
     {.0141460, .1194839, .2883075, .3266568, .1836014, .0569113, .0089426, .0019504},
     {.0115208, .1032624, .2716849, .3283426, .2021206, .0674456, .0128924, .0027307},
     {.0078498, .0802010, .2456595, .3308175, .2291646, .0836912, .0187016, .0039148},
     {.0046272, .0563321, .2132296, .3290407, .2599806, .1045974, .0265604, .0056322},
     {.0024659, .0360957, .1788634, .3210507, .2892537, .1282576, .0360887, .0079244},
     {.0012702, .0216090, .1472227, .3083032, .3123950, .1522540, .0462449, .0107009},
     {.0007288, .0134879, .1231200, .2949390, .3258251, .1731879, .0551737, .0135376},
     {.0004373, .0080115, .1002329, .2779283, .3342611, .1966100, .0650090, .0175099}
    };
  static double U235nubar [11] = {
      2.4370000, 
      2.5236700, 
      2.6368200, 
      2.7623400, 
      2.8738400, 
      3.0386999, 
      3.2316099, 
      3.4272800, 
      3.6041900, 
      3.7395900, 
      3.8749800
    };
  double fraction, r, cum;
  int engind, nu;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= U235nubar[0] && nubar <= U235nubar[10]) {
/*
     Use Zucker and Holden Data
*/
     engind = 1;
     while (nubar > U235nubar[engind]){ engind++;}
     fraction = (nubar-U235nubar[engind-1])/(U235nubar[engind]-U235nubar[engind-1]);
     if(fisslibrng() > fraction) engind--;

     r = fisslibrng();
     nu = 0;
     cum = U235nu[engind][0];
     while (r > cum && nu < 7){ 
       nu++;
       cum += U235nu[engind][nu];
     }
     return nu;
  } else {
/*
     Use Terrell's formula
*/
     return (int) SmpTerrell(nubar);
  }
}

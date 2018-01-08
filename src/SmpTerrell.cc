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
#include <string>
#include <sstream>
#include "fissionEvent.h"

#define TWOPI 6.283185307
#define SQRT2 1.414213562
#define BSHIFT -0.43287
#define WIDTH 1.079

double fissionEvent::SmpTerrell(double nubar) {
/*
  Description
    Sample Fission Number from Terrell's modified Gaussian distribution

    method uses Red Cullen's algoritm UCRL-TR-222526
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    SmpTerrell  - sampled multiplicity
    
*/

  double width;
  double temp1, temp2, expo, cshift;
  double rw, theta, sampleg;


  if (nubar < WIDTH) {
    std::ostringstream o;
    o << nubar;
    std::string errMsg = "fission nubar out of range, nubar=" + o.str();
    fissionerr(6, "SmpTerrell", errMsg);
  }

  width = SQRT2 * WIDTH;
  temp1 = nubar + 0.5;
  temp2 = temp1/width;
  temp2 *= temp2;
  expo = exp(-temp2);
  cshift = temp1 + BSHIFT * WIDTH * expo/(1. - expo);

  do {
    rw = sqrt(-log(fisslibrng()));
    theta = TWOPI * fisslibrng();
    sampleg = width * rw * cos(theta) + cshift;
  } while (sampleg < 0.0);

  return floor(sampleg);
}

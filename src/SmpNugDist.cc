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

#define nfissg 40
#define alphanegbin 26

int fissionEvent::SmpNugDist(int isotope, double nubar, double ePart, bool spontaneous) {

/*
  Description
    Sample Number of Photons from neutron induced fission in 
    all isotopes using Tim Valentine's model (negative binomial
    distribution, using nubar as a model parameter)
*/

/*
  Input
    isotope     - isotope
    nubar       - average number of fission neutrons
    ePart       - energy of incident neutron causing fission.
    spontaneous - true for spontaneous fission

  Output
    SmpNugDist - sampled multiplicity
*/
 
  static double logcoeff[nfissg+1] = {
     0.00000000000000e+00,
     3.25809653802149e+00,
     5.86078622346587e+00,
     8.09437844497297e+00,
     1.00753799138395e+01,
     1.18671393830676e+01,
     1.35093671183247e+01,
     1.50291928720691e+01,
     1.64462588918558e+01,
     1.77753948391357e+01,
     1.90281578076311e+01,
     2.02137814732888e+01,
     2.13397927361450e+01,
     2.24124295384099e+01,
     2.34369338549243e+01,
     2.44177631079360e+01,
     2.53587464524005e+01,
     2.62632027266277e+01,
     2.71340310844251e+01,
     2.79737817391769e+01,
     2.87847119553932e+01,
     2.95688309141589e+01,
     3.03279360625106e+01,
     3.10636428574894e+01,
     3.17774093252521e+01,
     3.24705565058120e+01,
     3.31442856005149e+01,
     3.37996924530920e+01,
     3.44377798564689e+01,
     3.50594680730467e+01,
     3.56656038766170e+01,
     3.62569683628670e+01,
     3.68342837279018e+01,
     3.73982191769817e+01,
     3.79493960962713e+01,
     3.84883925970040e+01,
     3.90157475227212e+01,
     3.95319639951220e+01,
     4.00375125617872e+01,
     4.05328339990172e+01,
     4.10183418147990e+01
  };
  int i;
  double cpi[nfissg+1];
  double p, q, nubarg;
  double r;

/* 
  No data is available for induced fission gamma number
  distributions. Sample the negative binomial cumulative 
  probability distribution.
*/
  nubarg = getNubarg(isotope, nubar, ePart, spontaneous);
  p = 1.*alphanegbin/(alphanegbin+nubarg);
  q = 1.-p;
  cpi[0] = exp(logcoeff[0]+26.*log(p));
  for (i=1; i<=nfissg; i++) cpi[i] = cpi[i-1] + exp(logcoeff[i]+26.*log(p)+i*log(q));
  for (i=0; i<=nfissg; i++) cpi[i] = cpi[i]/cpi[nfissg-1];

  r=fisslibrng();

  for(i=0; i<=nfissg; i++) if (r <= cpi[i]) return i;
  return 0;
}

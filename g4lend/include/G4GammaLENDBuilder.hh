//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4GammaLENDBuilder_h
#define G4GammaLENDBuilder_h 1

#include "globals.hh"

#include "G4PhotoFissionProcess.hh"
#include "G4PhotoCaptureProcess.hh"
#include "G4PhotoNuclearProcess.hh"

#include "G4LENDElasticCrossSection.hh"
#include "G4LENDElastic.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDInelasticCrossSection.hh"
#include "G4LENDFission.hh"
#include "G4LENDFissionCrossSection.hh"
#include "G4LENDCapture.hh"
#include "G4LENDCaptureCrossSection.hh"
#include "G4ZeroXS.hh"

class G4GammaLENDBuilder
{
  public: 
    G4GammaLENDBuilder(G4String eva="");
    virtual ~G4GammaLENDBuilder();

  public: 
    virtual void Build(G4PhotoNuclearProcess * aP);
    virtual void Build(G4PhotoFissionProcess * aP);
    virtual void Build(G4PhotoCaptureProcess * aP);

    void SetMinEnergy(G4double aM) 
    {
      theMin=aM;
      theIMin = theMin;
    }
    void SetMinInelasticEnergy(G4double aM) 
    {
      theIMin=aM;
    }
    void SetMaxEnergy(G4double aM) 
    {
      theIMax = aM;
      theMax=aM;
    }
    void SetMaxInelasticEnergy(G4double aM)
    {
      theIMax = aM;
    }


  private:

    G4double theMin;
    G4double theIMin;
    G4double theMax;
    G4double theIMax;

    G4LENDInelastic * theLENDInelastic;
    G4LENDInelasticCrossSection * theLENDInelasticCrossSection;
    G4LENDFission * theLENDFission;
    G4LENDFissionCrossSection * theLENDFissionCrossSection;
    G4LENDCapture * theLENDCapture;
    G4LENDCaptureCrossSection * theLENDCaptureCrossSection;

    G4String evaluation;
};

#endif


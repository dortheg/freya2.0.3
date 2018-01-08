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
#ifndef G4NeutronLLNLFissLENDBuilder_h
#define G4NeutronLLNLFissLENDBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4LENDElasticCrossSection.hh"
#include "G4LENDElastic.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDInelasticCrossSection.hh"
#include "G4LENDFission.hh"
#include "G4LENDFissionCrossSection.hh"
#include "G4LENDCapture.hh"
#include "G4LENDCaptureCrossSection.hh"
#include "G4NeutronLENDBuilder.hh"
#include "G4FissLib_new.hh"

class G4NeutronLLNLFissLENDBuilder : public G4NeutronLENDBuilder
{
  public: 
    G4NeutronLLNLFissLENDBuilder(G4String eva="");
    virtual ~G4NeutronLLNLFissLENDBuilder();

  public: 
    using G4NeutronLENDBuilder::Build;
    virtual void Build(G4HadronFissionProcess * aP);

    void SetMinEnergy(G4double aM)
    {
      theMin=aM;
    }
    void SetMaxEnergy(G4double aM)
    {
      theMax=aM;
    }

  private:

    G4double theMin;
    G4double theMax;

    G4FissLib_new * theFissionModel;
    G4LENDFissionCrossSection * theLENDFissionCrossSection;

    static double rng4llnlfisslib(void);

    G4String evaluation;
};

#endif

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
/// \file hadronic/Hadr04/src/NeutronHPphysics.cc
/// \brief Implementation of the NeutronHPphysics class
//
// $Id: NeutronHPphysics.cc 66587 2012-12-21 11:06:44Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NeutronHPphysics.hh"

#include "NeutronHPMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes

#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPInelastic.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPCapture.hh"

#include "G4HadronFissionProcess.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPFission.hh"

// neutron-induced fission (new implementation)
#ifdef FISSION_NEW
#include "fissionEvent.h"
#include "G4FissLib_new.hh"
#else
#include "G4FissLib.hh"
#endif

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::NeutronHPphysics(const G4String& name)
:  G4VPhysicsConstructor(name), fThermal(false), fNeutronMessenger(0)
{
  fNeutronMessenger = new NeutronHPMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::~NeutronHPphysics()
{
  delete fNeutronMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHPphysics::ConstructProcess()
{
   G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
   
   // process: elastic
   //
   G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
   pManager->AddDiscreteProcess(theElasticProcess);   
   //
   // cross section data set
   G4NeutronHPElasticData* dataSet1a = new G4NeutronHPElasticData();
   G4NeutronHPThermalScatteringData* dataSet1b 
                               = new G4NeutronHPThermalScatteringData();
   theElasticProcess->AddDataSet(dataSet1a);                               
   if (fThermal) theElasticProcess->AddDataSet(dataSet1b);
   //
   // models
   G4NeutronHPElastic*           model1a = new G4NeutronHPElastic();
   G4NeutronHPThermalScattering* model1b = new G4NeutronHPThermalScattering();
  if (fThermal)  model1a->SetMinEnergy(4*eV);
   theElasticProcess->RegisterMe(model1a);    
   if (fThermal) theElasticProcess->RegisterMe(model1b);
   
   // process: inelastic
   //
   G4NeutronInelasticProcess* theInelasticProcess = new G4NeutronInelasticProcess();
   pManager->AddDiscreteProcess(theInelasticProcess);   
   //
   // cross section data set
   G4NeutronHPInelasticData* dataSet2 = new G4NeutronHPInelasticData();
   theInelasticProcess->AddDataSet(dataSet2);                               
   //
   // models
   G4NeutronHPInelastic* model2 = new G4NeutronHPInelastic();
   theInelasticProcess->RegisterMe(model2);    

   // process: nCapture   
   //
   G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess();
   pManager->AddDiscreteProcess(theCaptureProcess);    
   //
   // cross section data set
   G4NeutronHPCaptureData* dataSet3 = new G4NeutronHPCaptureData();
   theCaptureProcess->AddDataSet(dataSet3);                               
   //
   // models
   G4NeutronHPCapture* model3 = new G4NeutronHPCapture();
   theCaptureProcess->RegisterMe(model3);
   
   // process: nFission   
   //
   G4HadronFissionProcess* theFissionProcess = new G4HadronFissionProcess();
   pManager->AddDiscreteProcess(theFissionProcess);    
   //
   // cross section data set
   G4NeutronHPFissionData* dataSet4 = new G4NeutronHPFissionData();
   theFissionProcess->AddDataSet(dataSet4);                               
   //
   // models
#ifdef FISSION_NEW
   // pass the random number generator to class fissionEvent
   fissionEvent::setRNGd(rng4llnlfisslib); 
#ifdef USEFREYA
   G4cout << "NeutronHPphysics: using new version of LLNL Fission Library (overriding built-in version) with FREYA turned on\n";
   fissionEvent::setCorrelationOption(3);
#else
   G4cout << "NeutronHPphysics: using new version of LLNL Fission Library (overriding built-in version) without FREYA\n";
#endif
   G4FissLib_new* theFissionModel = new G4FissLib_new;
#else
   G4cout <<"NeutronHPphysics: using built-in version of LLNL Fission Library\n";
   G4FissLib* theFissionModel = new G4FissLib;
#endif
   theFissionProcess->RegisterMe(theFissionModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double NeutronHPphysics::rng4llnlfisslib(void)
{
  return G4UniformRand();
}

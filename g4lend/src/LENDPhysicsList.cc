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
//
// $Id: LENDPhysicsList.cc,v 1.26 2009/11/15 14:48:40 maire Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LENDPhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4NeutronLLNLFissLENDBuilder.hh"
#include "G4GammaLENDBuilder.hh"
#include "G4FissLib_new.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LENDPhysicsList::LENDPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);
  evaluation="ENDF/BVII.1";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LENDPhysicsList::~LENDPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::ConstructMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

void LENDPhysicsList::ConstructBaryons()
{
  //  barions

  G4BaryonConstructor pBConstructor;
  pBConstructor.ConstructParticle();
  G4IonConstructor pIConstructor;
  pIConstructor.ConstructParticle();

  G4GenericIon::GenericIonDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructDecay();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

#include "G4HadronElasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4PhotoCaptureProcess.hh"
#include "G4PhotoFissionProcess.hh"

#include "G4CascadeInterface.hh"
#include "G4ZeroXS.hh"
#include "G4LENDElasticCrossSection.hh"
#include "G4LENDInelasticCrossSection.hh"
#include "G4LENDCaptureCrossSection.hh"
#include "G4LENDFissionCrossSection.hh"
#include "G4LENDElastic.hh"
#include "G4LENDCapture.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDFission.hh"

#include "G4HadronElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4NeutronRadCapture.hh"
#include "G4LFission.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::ConstructEM()
{
  // get particle iterator to add processes ---------------------
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = 0;
  theParticleIterator = theParticleTable->GetIterator();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "gamma") {

      // setup process for photon
      G4PhotoNuclearProcess* thePhotoNuclearProcess = new G4PhotoNuclearProcess();
      G4PhotoFissionProcess* thePhotoFissionProcess = new G4PhotoFissionProcess("photoFission");
      G4PhotoCaptureProcess* thePhotoCaptureProcess = new G4PhotoCaptureProcess("photoCapture");

      // setup model
      G4CascadeInterface* theGammaReaction = new G4CascadeInterface;
      theGammaReaction->SetMinEnergy(20*MeV);
      thePhotoNuclearProcess->RegisterMe( theGammaReaction );  

      // setup LEND
      G4GammaLENDBuilder LENDbuilder(evaluation.data());
      LENDbuilder.SetMaxEnergy(20.*MeV);
      LENDbuilder.Build(thePhotoNuclearProcess);
      LENDbuilder.Build(thePhotoFissionProcess);
      LENDbuilder.Build(thePhotoCaptureProcess);

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      pmanager->AddDiscreteProcess( thePhotoNuclearProcess );
      pmanager->AddDiscreteProcess( thePhotoFissionProcess );
      pmanager->AddDiscreteProcess( thePhotoCaptureProcess );
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);      

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);
    
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MuMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 4, 4);
             
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {
      //proton  
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);       
     
    } else if( particleName == "alpha" || 
	       particleName == "He3" )     {
      //alpha 
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
     
    } else if( particleName == "GenericIon" ) { 
      //Ions 
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);     
      
    } else if( particleName == "neutron" ) { 

      // setup process for neutron 
      G4HadronElasticProcess* theNeutronElasticProcess = new G4HadronElasticProcess();
      G4HadronCaptureProcess* theNeutronCaptureProcess = new G4HadronCaptureProcess();
      G4NeutronInelasticProcess* theNeutronInelasticProcess = new G4NeutronInelasticProcess();
      G4HadronFissionProcess* theNeutronFissionProcess = new G4HadronFissionProcess();

      //setup model
      G4HadronElastic* theElasticModel1 = new G4HadronElastic;
      theElasticModel1->SetMinEnergy( 19.*MeV );
      theNeutronElasticProcess->RegisterMe( theElasticModel1 );

      G4CascadeInterface* theBertiniForN = new G4CascadeInterface;
      theBertiniForN->SetMinEnergy( 19.*MeV );
      theBertiniForN->SetMaxEnergy( 6.*GeV );
      theNeutronInelasticProcess->RegisterMe(theBertiniForN);

      G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
      theCaptureModel->SetMinEnergy(19*MeV);
      theNeutronCaptureProcess->RegisterMe( theCaptureModel );

      G4FissLib_new* theFissionModel = new G4FissLib_new;
      theFissionModel->SetMinEnergy(15.99*MeV);
      theNeutronFissionProcess->RegisterMe(theFissionModel);

      //setup LEND XS for low energy
      G4NeutronLLNLFissLENDBuilder LENDbuilder(evaluation.data());
      LENDbuilder.SetMaxEnergy(19.*MeV);
      LENDbuilder.Build(theNeutronElasticProcess);
      LENDbuilder.Build(theNeutronInelasticProcess);
      LENDbuilder.Build(theNeutronCaptureProcess);
      LENDbuilder.Build(theNeutronFissionProcess);

      pmanager->AddDiscreteProcess( theNeutronElasticProcess );
      pmanager->AddDiscreteProcess( theNeutronInelasticProcess );
      pmanager->AddDiscreteProcess( theNeutronCaptureProcess );
      pmanager->AddDiscreteProcess( theNeutronFissionProcess );

    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,        -1, 2, 2);        
    }     
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void LENDPhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  // get particle iterator to add processes ---------------------
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = 0;
  theParticleIterator = theParticleTable->GetIterator();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LENDPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "LENDPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


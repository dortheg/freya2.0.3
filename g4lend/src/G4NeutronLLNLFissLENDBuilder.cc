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
#include "G4NeutronLLNLFissLENDBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "fissionEvent.h"

G4NeutronLLNLFissLENDBuilder::
G4NeutronLLNLFissLENDBuilder( G4String eva ) : G4NeutronLENDBuilder(eva)
{
  theFissionModel = 0;
  theLENDFissionCrossSection = 0;
  theMin = 0;
  theMax = 20*MeV;
  evaluation = eva;
}

G4NeutronLLNLFissLENDBuilder::
~G4NeutronLLNLFissLENDBuilder() 
{}

void G4NeutronLLNLFissLENDBuilder::Build(G4HadronFissionProcess * aP) {
  if(theFissionModel == 0) theFissionModel = new G4FissLib_new;
  theFissionModel->SetMinEnergy(theMin);
  theFissionModel->SetMaxEnergy(theMax);
  fissionEvent::setRNGd(rng4llnlfisslib);
#ifdef USEFREYA
  fissionEvent::setCorrelationOption(3);
  G4cout << "G4NeutronLLNLFissLENDBuilder: using new version of LLNL Fission Library (overriding built-in version) with FREYA turned on\n";
#else
  G4cout << "G4NeutronLLNLFissLENDBuilder: using new version of LLNL Fission Library (overriding built-in version) without FREYA\n";
#endif
  // cross section data set
  if(theLENDFissionCrossSection==0) theLENDFissionCrossSection=new G4LENDFissionCrossSection( G4Neutron::Neutron() );
  if ( evaluation != "" ) theLENDFissionCrossSection->ChangeDefaultEvaluation( evaluation );
  // theLENDFissionCrossSection->AllowNaturalAbundanceTarget();
  theLENDFissionCrossSection->AllowAnyCandidateTarget();
  aP->AddDataSet(theLENDFissionCrossSection);
  aP->RegisterMe(theFissionModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double G4NeutronLLNLFissLENDBuilder::rng4llnlfisslib(void)
{
  return G4UniformRand();
}


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
#include "G4GammaLENDBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4GammaLENDBuilder::
G4GammaLENDBuilder( G4String eva ) 
{
  theLENDInelastic = 0;
  theLENDInelasticCrossSection = 0;
  
  theLENDFission = 0;
  theLENDFissionCrossSection = 0;
  
  theLENDCapture = 0;
  theLENDCaptureCrossSection = 0;
  
  theMin = 0;
  theIMin = theMin;
  theMax = 20*MeV;
  theIMax = theMax;
  evaluation = eva;

}

G4GammaLENDBuilder::
~G4GammaLENDBuilder() 
{}

void G4GammaLENDBuilder::
Build(G4PhotoNuclearProcess * aP)
{
  if(theLENDInelastic==0) theLENDInelastic = new G4LENDInelastic( G4Gamma::Gamma() );
  theLENDInelastic->SetMinEnergy(theMin);
  theLENDInelastic->SetMaxEnergy(theMax);

  if ( evaluation != "" ) theLENDInelastic->ChangeDefaultEvaluation( evaluation );
  theLENDInelastic->AllowNaturalAbundanceTarget();
  if(theLENDInelasticCrossSection == 0) theLENDInelasticCrossSection = new G4LENDInelasticCrossSection( G4Gamma::Gamma() );
  if ( evaluation != "" ) theLENDInelasticCrossSection->ChangeDefaultEvaluation( evaluation );
  theLENDInelasticCrossSection->AllowNaturalAbundanceTarget();
  aP->AddDataSet(theLENDInelasticCrossSection);
  aP->RegisterMe(theLENDInelastic);
}

void G4GammaLENDBuilder::
Build(G4PhotoFissionProcess * aP)
{
  if(theLENDFission == 0) theLENDFission = new G4LENDFission( G4Gamma::Gamma() );
  theLENDFission->SetMinEnergy(theMin);
  theLENDFission->SetMaxEnergy(theMax);
  if ( evaluation != "" ) theLENDFission->ChangeDefaultEvaluation( evaluation );
  theLENDFission->AllowNaturalAbundanceTarget();
  if(theLENDFissionCrossSection==0) theLENDFissionCrossSection=new G4LENDFissionCrossSection( G4Gamma::Gamma() );
  if ( evaluation != "" ) theLENDFissionCrossSection->ChangeDefaultEvaluation( evaluation );
  theLENDFissionCrossSection->AllowNaturalAbundanceTarget();
  aP->AddDataSet(new G4ZeroXS);
  aP->AddDataSet(theLENDFissionCrossSection);
  aP->RegisterMe(theLENDFission);
}

void G4GammaLENDBuilder::
Build(G4PhotoCaptureProcess * aP)
{
  if(theLENDCapture==0) theLENDCapture = new G4LENDCapture( G4Gamma::Gamma() );
  theLENDCapture->SetMinEnergy(theMin);
  theLENDCapture->SetMaxEnergy(theMax);
  if ( evaluation != "" ) theLENDCapture->ChangeDefaultEvaluation( evaluation );
  theLENDCapture->AllowNaturalAbundanceTarget();
  if(theLENDCaptureCrossSection==0) theLENDCaptureCrossSection = new G4LENDCaptureCrossSection( G4Gamma::Gamma() );
  if ( evaluation != "" ) theLENDCaptureCrossSection->ChangeDefaultEvaluation( evaluation );
  theLENDCaptureCrossSection->AllowNaturalAbundanceTarget();
  aP->AddDataSet(new G4ZeroXS);
  aP->AddDataSet(theLENDCaptureCrossSection);
  aP->RegisterMe(theLENDCapture);
}

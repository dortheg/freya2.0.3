//******************************************************************************
// PhysicsList.cc
//
// Defines physics processes for this application.
//
// 1.01 JMV, LLNL, JUN-2014:  Second version.
//******************************************************************************
//
#include "globals.hh"
#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "NeutronHPphysics.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include <iomanip>

//----------------------------------------------------------------------------//
PhysicsList::PhysicsList():  G4VModularPhysicsList()
{
  defaultCutValue = 0.01*mm;
  SetVerboseLevel(1);

  // Neutron Physics
  RegisterPhysics( new NeutronHPphysics("neutronHP"));
}

//----------------------------------------------------------------------------//
PhysicsList::~PhysicsList()
{
}

//----------------------------------------------------------------------------//
// Create all particles that can appear in this application.
//----------------------------------------------------------------------------//
void PhysicsList::ConstructParticle()
{
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
}

//----------------------------------------------------------------------------//
void PhysicsList::SetCuts()
{
  // set cut values for gamma and neutron
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "neutron");
 
  //  SetCutsWithDefault();   
  if (verboseLevel > 0) DumpCutValuesTable();  
}

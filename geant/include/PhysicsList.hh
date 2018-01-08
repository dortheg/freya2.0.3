//******************************************************************************
// PhysicsList.hh
//
// This class is a class derived from G4VModularPhysicsList for constructing 
// particles and physical interaction processes.
//
// 1.01 JMV, LLNL, JUN-2014:  Second version.
//******************************************************************************
//
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList();

  protected:
    void ConstructParticle();
    void SetCuts();
};

#endif

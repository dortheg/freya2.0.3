//******************************************************************************
// SponFissIsotope.hh
//
// This class is a class derived from SingleSource for 
// constructing the process used to generate incident particles.
//
// 1.00 JMV, LLNL, OCT-2010: version compatible with Geant 4.9.3.
//******************************************************************************
// 
#ifndef SponFissIsotope_h
#define SponFissIsotope_h 1

#include "SingleSource.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Event.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4LLNLFission.hh"
#ifdef FISSION_NEW
#include "fissionEvent.h"
#else
#include "G4fissionEvent.hh"
#endif

class G4Event;

class SponFissIsotope : public SingleSource
{
  public:
    SponFissIsotope();    
    SponFissIsotope(G4int iso);    
    ~SponFissIsotope();

  public:
    void GeneratePrimaryVertex(G4Event* anEvent);
    // Set the verbosity level.
    void SetVerbosity(G4int verb) {verbosityLevel = verb;};

  private:
    G4int isotope;
    G4ThreeVector particle_polarization;
    G4ParticleDefinition* neutron_definition;
    G4ParticleDefinition* photon_definition;
    G4SPSPosDistribution* posDist;

    // Verbosity
    G4int verbosityLevel;

};

#endif

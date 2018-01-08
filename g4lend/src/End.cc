//
// Chance to do things at the end of the job
//
// Doug Wright, LLNL

#include "End.hh"

#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4PhysicsModelCatalog.hh"

#include "User.hh"
void End(){
    
    double Ekin = 10.*MeV;
    //printf("thick=%f\n",user.thickness);
    
    //.... materials must already be defined (but not necessarily used)
    const G4Material* material = G4Material::GetMaterial("G4_Pb");

    G4EmCalculator emCalculator;
    G4cout
    << "\n"
    << "Ekin (MeV) = " <<  Ekin
    
    << "  gamma atten. length (cm) = "
    << emCalculator.ComputeGammaAttenuationLength (Ekin, material)/cm

    << "  compute cross section per volume (compt) = "
    << emCalculator.ComputeCrossSectionPerVolume(Ekin, G4Gamma::Gamma(),"compt",material)

    << "  (photoFission) = "
    << emCalculator.ComputeCrossSectionPerVolume(Ekin, G4Gamma::Gamma(),"photoFission",material)
    
    // Get is only possible for materials actually used (not just defined)
    //<< "  get cross section per volume (compt) = "
    //<< emCalculator.GetCrossSectionPerVolume(10.*MeV, G4Gamma::Gamma(),"compt",material)
    << "\n" << G4endl;

    //....loop over all defined particles and then printout every process type, subtype, and name
    //    (this is useful reference for processing root file later, also should check that subtypes are unique)
    //
    //    developed from code found in BookForAppliDev.pdf (Example 7.2. An example of G4ParticleGunMessenger.cc.)
    //    and from /Users/wright20/cern/geant4.10.01/source/run/src/G4PhysicsListHelper.cc
    G4ParticleTable* table = G4ParticleTable::GetParticleTable();
    for( int i=0; i<table->entries(); i++)
    {
        G4ParticleDefinition* particle = table->GetParticle(i);
        G4ProcessVector* list = particle->GetProcessManager()->GetProcessList();
        for (int j=0; j<list->size(); j++)
            cout
                << particle->GetParticleName() << " "
                << particle->GetPDGEncoding() << " "
                << ((*list)[j])->GetProcessType() << " "
                << ((*list)[j])->GetProcessSubType() << " "
                << ((*list)[j])->GetProcessName()
                << endl;
        cout << endl;
    }

    //....check contents of G4PhysicsModelCatalog
    G4int entries = G4PhysicsModelCatalog::Entries();
    for (G4int i = 0; i < entries; ++i) {
        G4cout << i << ": " << G4PhysicsModelCatalog::GetModelName(i) << G4endl;
    }
}

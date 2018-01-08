//******************************************************************************
// SponFissIsotope.cc
//
// 1.00 JMV, LLNL, OCT-2010: version compatible with Geant 4.9.3.
//******************************************************************************
//
#include "SponFissIsotope.hh"
#include "G4SystemOfUnits.hh"

//----------------------------------------------------------------------------//
SponFissIsotope::SponFissIsotope()
{
}

//----------------------------------------------------------------------------//
SponFissIsotope::SponFissIsotope(G4int iso)
{
  neutron_definition = G4Neutron::Neutron();
  photon_definition = G4Gamma::Gamma();

  // verbosity
  verbosityLevel = 0;

  isotope = iso;
}

//----------------------------------------------------------------------------//
SponFissIsotope::~SponFissIsotope()
{
}

//----------------------------------------------------------------------------//
void SponFissIsotope::GeneratePrimaryVertex(G4Event* anEvent)
{ 
  // Generate a spontaneous fission using the fission library and emit
  // the neutrons and gamma-rays

  G4double time = GetParticleTime()/second;
#ifdef FISSION_NEW
  fissionEvent* fe = new fissionEvent(isotope, time, -1., 0., 0);
#ifdef USEFREYA
  if (3 == fe->getCorrelationOption()) {
    int err_len = 1000;
    char* error_message = new char[err_len];
    fe->getFREYAerrors(&err_len, error_message);
    if (err_len>1) {
      G4ExceptionDescription ed;
      ed << "Call to new fissionEvent("
         << "isotope=" << isotope << ", "
         << "time=" << time << ", "
         << "nubar=-1." << ", "
         << "eng=0." << ", "
         << "0) failed with error message from FREYA: "
         << G4endl
         << error_message;
      delete [] error_message;
      G4Exception("G4FissionLibrary_new::SampleMult", "freya001", FatalException,
                  ed);
    }
    delete [] error_message;
  }
#endif
#else
  G4fissionEvent* fe = new G4fissionEvent(isotope, time, -1., 0.);
#endif
  
  G4int nPrompt, gPrompt;
  nPrompt = fe->getNeutronNu();
  gPrompt = fe->getPhotonNu();

  if(verbosityLevel > 1) {
    G4cout << " nPrompt: " << nPrompt << G4endl
           << " gPrompt: " << gPrompt << G4endl;
  }

  // Position
  posDist = GetPosDist();
  G4ThreeVector sampled_particle_position = posDist->GenerateOne();

  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(sampled_particle_position, 0.);

  G4double mom, momx, momy, momz, eng;

  if(verbosityLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;

  G4DynamicParticle* it;
  // Build neutrons
  for(G4int i=0; i<nPrompt; i++)
  {
    it = new G4DynamicParticle();
    it->SetDefinition(neutron_definition);
    eng = fe->getNeutronEnergy(i)*MeV;
    it->SetKineticEnergy(eng);
    mom = it->GetTotalMomentum();

    momx = mom*fe->getNeutronDircosu(i);
    momy = mom*fe->getNeutronDircosv(i);
    momz = mom*fe->getNeutronDircosw(i);

    G4PrimaryParticle* particle = new G4PrimaryParticle(
                             neutron_definition,
                             momx, momy, momz,
                             eng);
    particle->SetMass(neutron_definition->GetPDGMass());
    particle->SetCharge(neutron_definition->GetPDGCharge());
    particle->SetPolarization(particle_polarization.x(),
                              particle_polarization.y(),
                              particle_polarization.z());

    if(verbosityLevel > 1){
      G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
      G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
      G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
    }
    vertex->SetPrimary(particle);
    delete it;
  }
  // Build gammas
  for(G4int i=0; i<gPrompt; i++)
  {
    it = new G4DynamicParticle();
    it->SetDefinition(photon_definition);
    eng = fe->getPhotonEnergy(i)*MeV;
    it->SetKineticEnergy(eng);
    mom = it->GetTotalMomentum();

    momx = mom*fe->getPhotonDircosu(i);
    momy = mom*fe->getPhotonDircosv(i);
    momz = mom*fe->getPhotonDircosw(i);

    G4PrimaryParticle* particle = new G4PrimaryParticle(
                             photon_definition,
                             momx, momy, momz,
                             eng);
    particle->SetMass(photon_definition->GetPDGMass());
    particle->SetCharge(photon_definition->GetPDGCharge());
    particle->SetPolarization(particle_polarization.x(),
                              particle_polarization.y(),
                              particle_polarization.z());

    if(verbosityLevel > 1){
      G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
      G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
      G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
    }

    vertex->SetPrimary(particle);
    delete it;
  }
  delete fe;
  vertex->SetT0(time*second);
//  G4cout << "         Time: "<<vertex->GetT0()/second << G4endl;

  anEvent->AddPrimaryVertex( vertex );
  if(verbosityLevel > 1)
    G4cout << " Primary Vetex generated !"<< G4endl;
}

#ifdef FISSION_NEW
#include "G4FissLib_new.hh"
#include "G4FissionLibrary_new.hh"
#include "G4SystemOfUnits.hh"

G4FissLib_new::G4FissLib_new()
 :xSec(0)
{
  SetMinEnergy(0.0);
  SetMaxEnergy(20.*MeV);
  if(!getenv("G4NEUTRONHPDATA")) {
     G4cout << "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files." << G4endl;
     throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
  }
  dirName = getenv("G4NEUTRONHPDATA");
  G4String tString = "/Fission/";
  dirName = dirName + tString;
  numEle = G4Element::GetNumberOfElements();
  theFission = new G4NeutronHPChannel[numEle];

  for (G4int i=0; i<numEle; i++)
  { 
//    G4cout << "G4FissLib_new::G4FissLib_new(): element "<< i << " : " << (*(G4Element::GetElementTable()))[i]->GetZ()<< G4endl;
    if((*(G4Element::GetElementTable()))[i]->GetZ()>89)
    {
      theFission[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      theFission[i].Register(new G4FissionLibrary_new);
    }
  }
}
  
G4FissLib_new::~G4FissLib_new()
{
  delete [] theFission;
}
  
G4HadFinalState*
G4FissLib_new::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus&)
{
  const G4Material* theMaterial = aTrack.GetMaterial();
  G4int n = theMaterial->GetNumberOfElements();
  G4int index = theMaterial->GetElement(0)->GetIndex();

  if (n != 1) {
    xSec = new G4double[n];
    G4double sum = 0;
    G4int i;
    G4int imat;
    const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
    G4double rWeight;    
    G4NeutronHPThermalBoost aThermalE;
    for (i = 0; i < n; i++) {
      imat = theMaterial->GetElement(i)->GetIndex();
      rWeight = NumAtomsPerVolume[i];
      xSec[i] = theFission[imat].GetXsec(aThermalE.GetThermalEnergy(aTrack,
		                                                    theMaterial->GetElement(i),
  						      theMaterial->GetTemperature()));
      xSec[i] *= rWeight;
      sum+=xSec[i];
    }

    G4double random = G4UniformRand();
    G4double running = 0;
    for (i = 0; i < n; i++) {
      running += xSec[i];
      index = theMaterial->GetElement(i)->GetIndex();
      if(random<=running/sum) break;
    }
    delete [] xSec;
  }

  return theFission[index].ApplyYourself(aTrack);
}

const std::pair<G4double, G4double> G4FissLib_new::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus (taken from G4LFission)
  return std::pair<G4double, G4double>(5*perCent,250*GeV);
}
#endif

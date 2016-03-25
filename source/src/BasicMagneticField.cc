//There's a lot of functionality here that I'm not currently using since I just copied this class over from the errorprop example, and combined it with pieces from the B1 example.
//It doesn't hurt however so I'm leaving it in for now.

#include "BasicMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BasicMagneticField::BasicMagneticField()
  : G4UniformMagField(G4ThreeVector())
{
  GetGlobalFieldManager()->SetDetectorField(this);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BasicMagneticField::BasicMagneticField(G4ThreeVector fieldVector)
  : G4UniformMagField(fieldVector)
{
  GetGlobalFieldManager()->SetDetectorField(this);    
  GetGlobalFieldManager()->CreateChordFinder(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BasicMagneticField::SetFieldValue(G4double fieldValue)
{
  G4UniformMagField::SetFieldValue(G4ThreeVector(0,0,fieldValue));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BasicMagneticField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector!=G4ThreeVector(0.,0.,0.))
  { 
    G4UniformMagField::SetFieldValue(fieldVector);
    fieldMgr->SetDetectorField(this);
  } else {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    G4MagneticField* magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BasicMagneticField::~BasicMagneticField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4FieldManager*  BasicMagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}


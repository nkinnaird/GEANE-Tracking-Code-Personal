

#include "BasicPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4PhysicalConstants.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicPrimaryGeneratorAction::BasicPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  particleName("e+")
 // fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleMomentum(2000.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicPrimaryGeneratorAction::~BasicPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  // G4double envSizeXY = 0;
    G4double envSizeY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("TruthPlaneZero");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    // envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeY = fEnvelopeBox->GetYHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    // msg << "The gun will be place at the center.";
    G4Exception("BasicPrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  G4double size = 0.2; 
  // G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 =  envSizeY * (G4UniformRand()-0.5);
  G4double z0 =  size * envSizeZ * (G4UniformRand()-0.5);
  // G4double z0 = -5.*cm; //to fix postion of particle gun for now

  // fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));


  // fParticleGun->SetParticlePosition(G4ThreeVector(-100.,y0,z0)); // uses mm for these values 
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));

  /////////////////////////////////////////////////////////////////////////////////////

  double angle = 28.; //If field is on.
  double angleForNoField = 0.;//15.; //If field is off.


  double randomAngle = 1.*(G4UniformRand()-0.5);
  G4ThreeVector particleStartPos((-4600.*std::sin(angle*pi/180)+z0*std::sin(angle*pi/180))*mm,y0,((4600.-4600.*std::cos(angle*pi/180))+z0*std::cos(angle*pi/180))*mm); 
  // Center point comes from a radius of 4.6 m for 2 GeV positrons, rotated ~2 meters in a circle out corresponding to an angle of 26 degrees (a little before the 25 degree rotated 0 plane), then + some random number multipled by sin and cos in x and z to randomly distribute generated particles properly along the rotated plane. (Since x and z rotate.)
  fParticleGun->SetParticlePosition(particleStartPos);

  G4ThreeVector particleMomDirection(std::cos((angle-angleForNoField+randomAngle/2)*pi/180),std::sin((randomAngle)*pi/180),-std::sin((angle-angleForNoField+randomAngle/2)*pi/180)); // Cos angle and -Sin angle degrees for particles starting a little before the 0 plane to be angled correctly.
  fParticleGun->SetParticleMomentumDirection(particleMomDirection); // Have to divide by 20 for some reason for the angle to not be too large - I think because the momentum vector direction is no longer a unit vector

  G4double particleMomentum = 2250.+(1500.)*(G4UniformRand()-0.5)*MeV; // Energies evenly distributed from 1 GeV to 3 GeV.
  fParticleGun->SetParticleMomentum(particleMomentum);
  /////////////////////////////////////////////////////////////////////////////////////



  // fParticleGun->SetParticlePosition(G4ThreeVector(-100.,0,0)); // uses mm for these values


  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


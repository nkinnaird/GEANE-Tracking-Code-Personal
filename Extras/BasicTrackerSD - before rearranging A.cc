

#include "BasicTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"

#include "CLHEP/Matrix/Matrix.h" // Should replace this with the eigen package at some point.

// #include <Eigen/Dense>


// #include "BasicEventAction.hh"

const double noHit = 900000.; // Unphysical double value to signify the lack of a hit within a tracker plane.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicTrackerSD::BasicTrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicTrackerSD::~BasicTrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new BasicTrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool BasicTrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  

  CLHEP::HepMatrix ZYtoUVcoordinateTransformationMatrix(2,2);

  double rotationAngle = 7.5*pi/180;

  ZYtoUVcoordinateTransformationMatrix[0][0]=std::cos(rotationAngle);
  ZYtoUVcoordinateTransformationMatrix[0][1]=-std::sin(rotationAngle);
  ZYtoUVcoordinateTransformationMatrix[1][0]=std::cos(rotationAngle);
  ZYtoUVcoordinateTransformationMatrix[1][1]=std::sin(rotationAngle);


 // Statements for truth planes.
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  G4double energyThreshold = 500.*MeV;

  G4StepPoint* pPreStepPoint = aStep->GetPreStepPoint();

  if(pPreStepPoint->GetStepStatus() == fGeomBoundary && aStep->GetTrack()->GetTrackID() == 1 && pPreStepPoint->GetKineticEnergy() > energyThreshold)
  { // Boundary check for step just crossing into sensitive detector, trackID check to make sure we're still looking at the primary positron.
    // Energy check to remove low energy tracks that circle through the detector multiple times - the first part of the track will still be included which is fine.

    // G4cout << "Track " << aStep->GetTrack()->GetTrackID() << " Hit at " << aStep->GetPreStepPoint()->GetPosition() << G4endl;
    // G4cout << hcID << G4endl;

    BasicTrackerHit* newHit = new BasicTrackerHit();

    newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
    newHit->SetPos (pPreStepPoint->GetPosition()); 
    newHit->SetMom (pPreStepPoint->GetMomentum()); 
    newHit->SetGlobalTime (pPreStepPoint->GetGlobalTime()); 
    newHit->SetProperTime (pPreStepPoint->GetProperTime()); 
    newHit->SetCopyNo (pPreStepPoint->GetPhysicalVolume()->GetCopyNo());


    if (pPreStepPoint->GetPhysicalVolume()->GetName() == "TruthPlaneZero")
    {
        // newHit->SetUPos (pPreStepPoint->GetPosition().z()); // Set U and V to Y and Z respectively for now.
        // newHit->SetVPos (pPreStepPoint->GetPosition().y());

        newHit->SetUPos (ZYtoUVcoordinateTransformationMatrix[0][0]*pPreStepPoint->GetPosition().z() + ZYtoUVcoordinateTransformationMatrix[0][1]*pPreStepPoint->GetPosition().y());
        newHit->SetVPos (ZYtoUVcoordinateTransformationMatrix[1][0]*pPreStepPoint->GetPosition().z() + ZYtoUVcoordinateTransformationMatrix[1][1]*pPreStepPoint->GetPosition().y());

    }
    else if (pPreStepPoint->GetPhysicalVolume()->GetName() == "TruthPlaneU")
    {
        newHit->SetUPos (ZYtoUVcoordinateTransformationMatrix[0][0]*pPreStepPoint->GetPosition().z() + ZYtoUVcoordinateTransformationMatrix[0][1]*pPreStepPoint->GetPosition().y());
        newHit->SetVPos (noHit);

        // newHit->SetUPos (pPreStepPoint->GetPosition().z());
        // newHit->SetVPos (noHit);
    }
    else if (pPreStepPoint->GetPhysicalVolume()->GetName() == "TruthPlaneV")
    {
        newHit->SetUPos (noHit);
        newHit->SetVPos (ZYtoUVcoordinateTransformationMatrix[1][0]*pPreStepPoint->GetPosition().z() + ZYtoUVcoordinateTransformationMatrix[1][1]*pPreStepPoint->GetPosition().y());
    
        // newHit->SetUPos (noHit);
        // newHit->SetVPos (pPreStepPoint->GetPosition().y());
    }



    fHitsCollection->insert( newHit ); // Don't insert hits for truth values once I turn on straw tracker SD.

    G4cout << G4endl << "The copyNo of the hit plane is: " << pPreStepPoint->GetPhysicalVolume()->GetCopyNo() << G4endl;
  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  


/*
// Statements for straw tubes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  G4StepPoint* pPreStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();

  if(aStep->GetTrack()->GetTrackID() == 1)
  { // Overarching if statement since I only care about the "main" particle and no secondaries.

  if(pPreStepPoint->GetStepStatus() == fGeomBoundary) // Add a Volume name if statement here if I'm going to have 2 sensitive detectors referencing the same class - still might not work depending on how it creates hit collections.
  { // Will probably need a statement like this for the pPostStepPoint to make sure I get the point when the particle is exiting the straw.

  if (aStep->IsFirstStepInVolume())
  { G4cout << "Supposedly FIRST step in volume." << G4endl; }
  else if(aStep->IsLastStepInVolume())
  { G4cout << "Supposedly LAST step in volume." << G4endl; }
  else { G4cout << "Not first or last step, hopefully a middle step." << G4endl; }

  G4ThreeVector chordStart = pPreStepPoint->GetPosition();
  G4ThreeVector chordEnd = pPostStepPoint->GetPosition();


  G4RotationMatrix* strawRotation = pPreStepPoint->GetPhysicalVolume()->GetObjectRotation();
  G4ThreeVector strawTranslationVector = pPreStepPoint->GetPhysicalVolume()->GetObjectTranslation();

  G4ThreeVector topPointVector;
  topPointVector.setY(strawRotation->yz()*50.); // These +- 50's are the z ends of the straw before it is rotated into its proper orientation, and then translated to it's proper place.
  topPointVector.setZ(strawRotation->zz()*50.);  

  G4ThreeVector bottomPointVector;
  bottomPointVector.setY(strawRotation->yz()*-50.);
  bottomPointVector.setZ(strawRotation->zz()*-50.);  

  topPointVector = topPointVector+strawTranslationVector;
  bottomPointVector = bottomPointVector+strawTranslationVector;

  G4ThreeVector wireLineVector = topPointVector - bottomPointVector;
  G4ThreeVector chordLineVector = chordStart - chordEnd;
  G4ThreeVector normalVector = wireLineVector.cross(chordLineVector);
  G4ThreeVector distanceVector = chordStart - strawTranslationVector;
  G4double closestApproachDistance = std::abs(distanceVector.dot(normalVector))/normalVector.mag(); // To get the closest distance of approach, find the normal vector between the two parallel planes formed by the two separate vectors, then project a vector connecting a point on each line onto this normal vector.

  CLHEP::HepMatrix coordinateTransformationMatrix(3,3); // Matrix of basis vectors for UV coordinate system. Made from rotating the Z axis -7.5 degrees about X to get the U axis, and rotating the Y axis about X by -82.5 degrees to get the V axis.
  coordinateTransformationMatrix[0][0]=1.;
  coordinateTransformationMatrix[0][1]=coordinateTransformationMatrix[0][2]=coordinateTransformationMatrix[1][0]=coordinateTransformationMatrix[2][0]=0.;
  coordinateTransformationMatrix[1][1]=std::sin(-7.5*CLHEP::deg);
  coordinateTransformationMatrix[1][2]=std::cos(-82.5*CLHEP::deg);
  coordinateTransformationMatrix[2][1]=std::cos(-7.5*CLHEP::deg);
  coordinateTransformationMatrix[2][2]=-std::sin(-82.5*CLHEP::deg);

  // std::cout << coordinateTransformationMatrix;

  CLHEP::HepMatrix transformInverse(3,3); // Have to multiply the inverse of the basis transformation matrix by the vector in X,Y,Z space to get to X,U,V space.
  transformInverse = coordinateTransformationMatrix.inverse();

  // std::cout << transformInverse;

// Was testing some vector stuff here, can probably delete.
/////////////////////////////////////////////////////////////////////////////////////
//   G4ThreeVector testVector(5,8,-1.05322);
//   G4ThreeVector testVectorNew;
//   CLHEP::Hep3Vector testVector2(0,0,2);

// // (transformInverse) * (testVector2); // Can't get this to work for some reason.

// testVectorNew.setX(transformInverse[0][0]*testVector.x()+transformInverse[0][1]*testVector.y()+transformInverse[0][2]*testVector.z());
// testVectorNew.setY(transformInverse[1][0]*testVector.x()+transformInverse[1][1]*testVector.y()+transformInverse[1][2]*testVector.z());
// testVectorNew.setZ(transformInverse[2][0]*testVector.x()+transformInverse[2][1]*testVector.y()+transformInverse[2][2]*testVector.z());
/////////////////////////////////////////////////////////////////////////////////////


G4ThreeVector strawVectorUV; // Will want to just grab U part of this vector for the U position of where the track went. (With U sensitive dector straws.)
strawVectorUV.setX(transformInverse[0][0]*strawTranslationVector.x()+transformInverse[0][1]*strawTranslationVector.y()+transformInverse[0][2]*strawTranslationVector.z());
strawVectorUV.setY(transformInverse[1][0]*strawTranslationVector.x()+transformInverse[1][1]*strawTranslationVector.y()+transformInverse[1][2]*strawTranslationVector.z());
strawVectorUV.setZ(transformInverse[2][0]*strawTranslationVector.x()+transformInverse[2][1]*strawTranslationVector.y()+transformInverse[2][2]*strawTranslationVector.z());





  G4cout << "Copy number preStepPoint: " <<  pPreStepPoint->GetPhysicalVolume()->GetCopyNo() << G4endl
         << "Copy number postStepPoint: " <<  pPostStepPoint->GetPhysicalVolume()->GetCopyNo() << G4endl
         << "Straw Center Position: " << strawTranslationVector << G4endl
         // << "Angled straw top point: " << topPointVector << G4endl
         // << "Angled straw bottom point: " << bottomPointVector << G4endl
         // << "Chord start point: " << chordStart << G4endl
         // << "Chord end point: " << chordEnd << G4endl
         << "ClosestApproachDistance: " << closestApproachDistance << G4endl
         << "Transformed straw vector: " << strawVectorUV << G4endl
         << "Step length: " << aStep->GetStepLength() << G4endl << G4endl;




    // G4cout << "Track " << aStep->GetTrack()->GetTrackID() << " Hit at " << aStep->GetPreStepPoint()->GetPosition() << G4endl;
    // G4cout << hcID << G4endl;

    BasicTrackerHit* newHit = new BasicTrackerHit();

    newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
    newHit->SetPos (pPreStepPoint->GetPosition()); 
    newHit->SetMom (pPreStepPoint->GetMomentum()); 
    newHit->SetGlobalTime (pPreStepPoint->GetGlobalTime()); 
    newHit->SetProperTime (pPreStepPoint->GetProperTime()); 

    fHitsCollection->insert( newHit );

  }
  
  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/


  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void BasicTrackerSD::EndOfEvent(G4HCofThisEvent*)
{

  // Printing hit collection properties from BasicEventAction.cc, so this method unnecessary for now.

  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

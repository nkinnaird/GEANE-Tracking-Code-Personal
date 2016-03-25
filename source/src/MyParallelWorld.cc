#include "BasicDetectorConstruction.hh"
#include "MyParallelWorld.hh"

#include "BasicTrackerSD.hh"
// #include "G4SDManager.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"


#include "G4Material.hh"

MyParallelWorld::MyParallelWorld(G4String worldName)
	: G4VUserParallelWorld(worldName),
  truthPlaneHalfX(0.01*cm), truthPlaneHalfY(5.5*cm), truthPlaneHalfZ(11.*cm)
	{;}

MyParallelWorld::~MyParallelWorld(){;}

void MyParallelWorld::Construct()
{
	G4VPhysicalVolume* ghostWorld = GetWorld();
	G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  BasicDetectorConstruction* myPhysicalWorld = new BasicDetectorConstruction(); // This here so that I can get the main world physical parameters.
  // Tried to get it through ghostWorld but I'm not sure if that's possible.

	  //Vacuum
	  // G4double a = 1.*g/mole;
	  // G4double density = 1.E-9*g/cm3;
	  // G4double z;
	  // G4Material* Vacuum = new G4Material("Vacuum", z=1., a, density);
    G4Material* Vacuum = G4Material::GetMaterial("Vacuum"); // Should already be defined or shouldn't matter since we're in the ghost world.

	   G4int moduleNum = 8;
     G4int layerNum = 2;

     G4double moduleSeparation = myPhysicalWorld->GetmoduleSeparation();
     G4double layerSeparation = myPhysicalWorld->GetlayerSeparation();
     G4double UVSeparation = myPhysicalWorld->GetUVSeparation();
     // G4double HalfYZLength = myPhysicalWorld->GetHalfYZLength();
     // G4double truthPlaneHalfX = myPhysicalWorld->GetTruthPlaneHalfX();
     // G4double truthPlaneHalfY = myPhysicalWorld->GetTruthPlaneHalfY();
     // G4double truthPlaneHalfZ = myPhysicalWorld->GetTruthPlaneHalfZ();



/////////////////////////////////////////////////////////////////////////////////////
  //------------------------------ 
  // Truth Detector Planes - planes that lie at the center of modules to acquire truth information for traceback
  //------------------------------   

     G4VisAttributes * TruthVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
     TruthVisAtt->SetForceWireframe(true);

     logicTruthPlane = new G4LogicalVolume*[1+moduleNum*layerNum*2];


      G4double Xposition = 0.;
      G4double Yposition = 0.;
      G4double Zposition = 0.;

      // G4RotationMatrix* PlaneRotation = new G4RotationMatrix();
      // PlaneRotation->rotateX(0.*deg);
      // PlaneRotation->rotateY(0.*deg);
      // PlaneRotation->rotateZ(0.*deg);

      G4ThreeVector planeTranslationPoint = G4ThreeVector(Xposition, Yposition, Zposition);
      G4Box* TruthPlaneZero;
      G4Box* TruthPlane;

/////////////////////////////////////////////////////////////////////////////////////
// Create the 0 plane separately.

              double angle = 25.;

              Xposition = (-460.*std::sin(angle*pi/180))*cm; // These for a turned plane ~2 meters out in a circle for a 4.6 m radius for 2 GeV positrons, corresponding to an angle of 25 degrees.
              Zposition = (460.-460.*std::cos(angle*pi/180))*cm;
              planeTranslationPoint.setX(Xposition);
              planeTranslationPoint.setZ(Zposition);
              // PlaneRotation->rotateY(-angle*deg);

              TruthPlaneZero = new G4Box("TruthPlaneZero",truthPlaneHalfX,2.*truthPlaneHalfY,10.*truthPlaneHalfZ); 

              logicTruthPlane[0] = new G4LogicalVolume(TruthPlaneZero,Vacuum,"TruthPlaneZero",0,0,0);
              logicTruthPlane[0]->SetVisAttributes(TruthVisAtt);

              new G4PVPlacement(0, planeTranslationPoint,"TruthPlaneZero",logicTruthPlane[0],ghostWorld,false,0);


/////////////////////////////////////////////////////////////////////////////////////
           G4String planeName;

    for (int i = 1; i < moduleNum+1; ++i)
     {
      for (int k = 0; k < 2; ++k)
      { // loop for U or V planes
       for (int j = 0; j < layerNum; ++j)
       { 

          if (k == 0)
          { planeName = "TruthPlaneU"; }
          else if (k==1)
          { planeName = "TruthPlaneV"; }
          
        
           Xposition = (i-1) * moduleSeparation + j * layerSeparation + k * UVSeparation + truthPlaneHalfX; // + truthPlaneHalfX so the front face of the box is at the center of the module.
           Yposition=0.;
           Zposition=0.;

           planeTranslationPoint = G4ThreeVector(Xposition, Yposition, Zposition);

           TruthPlane = new G4Box(planeName,truthPlaneHalfX,truthPlaneHalfY,truthPlaneHalfZ);

           logicTruthPlane[(i-1)*4+j+k*2+1] = new G4LogicalVolume(TruthPlane,Vacuum,planeName,0,0,0);
           logicTruthPlane[(i-1)*4+j+k*2+1]->SetVisAttributes(TruthVisAtt);

           new G4PVPlacement(0, planeTranslationPoint,planeName,logicTruthPlane[(i-1)*4+j+k*2+1],ghostWorld,false,(i-1)*4+j+k*2+1);

       }
      }
     }

/////////////////////////////////////////////////////////////////////////////////////







     // // Below for 16 parallel world detector planes, 2 per module, one for U and V each at the center of the 2 U or V layers. 
     // logicTruthPlane = new G4LogicalVolume*[moduleNum*layerNum];

     // for (int i = 0; i < moduleNum; ++i)
     // {
     //   for (int j = 0; j < layerNum; ++j)
     //   {
     //     // for (int copyNo = 0; copyNo < 2; ++copyNo) // Loop for U or V placement.
     //     // {

     //      //NOTE: Keep in mind that the center of the box x position is given by below, so if I want to measure U or V at the center of the straw, I need to offset the box so the edge of the box lies along the center plane of the straws.

     //       G4double Xposition = 5.*cm + i * moduleSeparation + j * UVSeparation + layerSeparation/2. + truthPlaneHalfX;

     //       G4Box* TruthPlane = new G4Box("TruthPlane",0.1*cm,HalfYZLength,HalfYZLength);

     //       logicTruthPlane[(j+(2*i))] = new G4LogicalVolume(TruthPlane,Vacuum,"TruthPlane",0,0,0); // (j+(2*i)) is the pattern to get successive numbers for plane numbers, with 8 modules have 4 planes each.
     //       logicTruthPlane[(j+(2*i))]->SetVisAttributes(TruthVisAtt);

     //       new G4PVPlacement(0, G4ThreeVector(Xposition,0,0),logicTruthPlane[(j+(2*i))],"TruthPlane",worldLogical,false,(j+(2*i)));

     //     // }
     //   }
     // }











       G4VisAttributes* parallelworldVisAtt = new G4VisAttributes(0);
       // G4VisAttributes* parallelworldVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.));
      worldLogical->SetVisAttributes(parallelworldVisAtt);

}

void MyParallelWorld::ConstructSD()
{

  G4String truthTrackerChamberSDname = "Basic/TruthTrackerChamberSD";
  BasicTrackerSD* truthTrackerSD = new BasicTrackerSD(truthTrackerChamberSDname,
                                            "TruthTrackerHitsCollection");

  // Setting truthTrackerSD to all logical volumes with the same name.
    SetSensitiveDetector("TruthPlaneZero", truthTrackerSD, true);
        SetSensitiveDetector("TruthPlaneU", truthTrackerSD, true);
            SetSensitiveDetector("TruthPlaneV", truthTrackerSD, true);


}

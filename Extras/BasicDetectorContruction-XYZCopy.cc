

#include "BasicDetectorConstruction.hh"
// #include "BasicDetectorMessenger.hh"
#include "BasicMagneticField.hh"

#include "BasicTrackerSD.hh"

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BasicDetectorConstruction::BasicDetectorConstruction()
  : G4VUserDetectorConstruction(),
    /*HalfAirGapLength(5.*cm),*/ HalfStrawLength(5.*cm), outerRadiusOfTheMylarStraw(.255*cm), 
    distBtwnWires(.6*cm), moduleSeparation(13.*cm), layerSeparation(.6*cm), UVSeparation(2.*cm),
    outerRadiusOfTheGas(.2535*cm), outerRadiusOfTheWire(.00125*cm), 
    truthPlaneHalfX(0.01*cm), truthPlaneHalfY(5.5*cm), truthPlaneHalfZ(11.*cm),
    HalfYZLength(200.*cm), HalfWorldLength(1000.*cm), 
    fUserLimits(0), fMagField(0)//, fDetectorMessenger(0) 
{

  // create UserLimits
  fUserLimits = new G4UserLimits();

  fMagField = new BasicMagneticField(
                    G4ThreeVector(0.0*kilogauss,14.5*kilogauss,0.0*kilogauss)); // 14.5 kG in Y
  // fDetectorMessenger = new BasicDetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BasicDetectorConstruction::~BasicDetectorConstruction()
{
  delete fMagField;
  // delete fDetectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* BasicDetectorConstruction::Construct()
{
//--------- Material definition ---------

  //Vacuum
  G4double a = 1.*g/mole;
  G4double density = 1.E-16*g/cm3;
  G4double z;
  G4Material* Vacuum = new G4Material("Vacuum", z=1., a, density);
  

  G4NistManager* nistMgr = G4NistManager::Instance();
  // G4Material* air = nistMgr->FindOrBuildMaterial("G4_AIR");
  //Al
  // G4Material* al = nistMgr->FindOrBuildMaterial("G4_Al");
  G4Material* tungsten = nistMgr->FindOrBuildMaterial("G4_W"); // The wire is "gold plated" tungsten, not sure if that matters.

  // G4Material* plastic = nistMgr->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* mylar = nistMgr->FindOrBuildMaterial("G4_MYLAR"); // We call it "aluminized mylar," not sure if that is different...

  G4Material* argon = nistMgr->FindOrBuildMaterial("G4_Ar");
  G4Material* C02 = nistMgr->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    
  /*G4double trackerGasDensity = argon->GetDensity() * .8 + C02->GetDensity() * .2;
  G4cout << "Argon gas density is here :" << argon->GetDensity() << G4endl;  
  G4cout << "Argon chemical formula is here :" << argon->GetChemicalFormula() << G4endl;  
  G4cout << "C02 gas density is here :" << C02->GetDensity() << G4endl;  
  G4cout << "C02 chemical formula is here :" << C02->GetChemicalFormula() << G4endl; 
  G4cout << "Tracker gas density is here :" << trackerGasDensity << G4endl;   */
  // Having trouble with the above to get gas densities - keep getting e+16 numbers, and argon chemical formula
  // doesn't show up even though it apparently finds the material. For now I will plug in densities manually
  // from http://geant4.cern.ch/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/apas08.html.

  // G4double trackerGasDensity = 0.00166201 * .8 + 0.00184212 * .2;
  G4double trackerGasDensity = 1.8223*mg/cm3; // This from a geant4 materials example. http://geant4.web.cern.ch/geant4/UserDocumentation/Doxygen/examples_doc/html_TestEm8/html/DetectorConstruction_8cc_source.html

  G4Material* trackerGas = new G4Material("Gas",trackerGasDensity,2); // 2 For num components.
  // trackerGas->AddMaterial(argon, 80*perCent);
  // trackerGas->AddMaterial(C02, 20*perCent);
  trackerGas->AddMaterial(argon, 0.783);
  trackerGas->AddMaterial(C02, 0.217);



  // Print all the materials defined.
  //
  // G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  
  
   
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 
  //-  G4double HalfWorldLength = fXWorldLength;
  G4cout << " HalfWorldLength " << HalfWorldLength << G4endl;
  
  G4Box* solidWorld= new G4Box("world",HalfWorldLength,HalfYZLength,HalfYZLength);
  G4LogicalVolume* logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World", 0, 0, 0);
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* physiWorld 
    = new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        "World",         // its name
                        logicWorld,      // its logical volume
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0);              // no field specific to volum
                                 
  

  //------------------------------ 
  // Straw Volumes
  //------------------------------   


     G4int numStrawsLayer = 32;
     G4int moduleNum = 8;
     G4int layerNum = 2;
     G4int totalStraws = numStrawsLayer * moduleNum * layerNum;

     G4RotationMatrix* URotation = new G4RotationMatrix();
     URotation->rotateX(82.5*deg);
     URotation->rotateY(0.*deg);
     URotation->rotateZ(0.*deg);

     G4RotationMatrix* VRotation = new G4RotationMatrix();
     VRotation->rotateX(97.5*deg);
     VRotation->rotateY(0.*deg);
     VRotation->rotateZ(0.*deg);

     G4VisAttributes* StrawAtt = new G4VisAttributes(G4Colour(0.,1.,1.));
     StrawAtt->SetForceWireframe(false);

     ////Array of logical volumes to fill with straw components - hopefully this isn't taking up too much memory.
     // It's probably possible to do this in a better way but I'm not sure...
     logicMylarStrawU = new G4LogicalVolume*[totalStraws];
     logicMylarStrawV = new G4LogicalVolume*[totalStraws];

     logicGasStrawU = new G4LogicalVolume*[totalStraws];
     logicGasStrawV = new G4LogicalVolume*[totalStraws];

     logicWireStrawU = new G4LogicalVolume*[totalStraws];
     logicWireStrawV = new G4LogicalVolume*[totalStraws];

  for (G4int i = 0; i < moduleNum; ++i)
  {
    for (G4int j = 0; j < layerNum; ++j)
    {        
      for (G4int copyNo=0; copyNo<numStrawsLayer; copyNo++) 
      {
          
           ///physical geometry of straw components
           G4Tubs* MylarStraw = new G4Tubs("MylarStraw",outerRadiusOfTheGas,outerRadiusOfTheMylarStraw,HalfStrawLength,0.*deg,360.*deg);
           G4Tubs* GasStraw = new G4Tubs("GasStraw",outerRadiusOfTheWire,outerRadiusOfTheGas,HalfStrawLength,0.*deg,360.*deg);
           G4Tubs* WireStraw = new G4Tubs("WireStraw",0,outerRadiusOfTheWire,HalfStrawLength,0.*deg,360.*deg);
           
           ////construction of all u layers//////////

           ///positions first
           G4double Xposition = 5.*cm + i * moduleSeparation + j * layerSeparation;
           G4double Yposition = 0.*cm;
           G4double Zposition = -9.3*cm + copyNo * distBtwnWires + j * layerSeparation/2; //Distance between wires in the same plane of 32. Spacing between wires at ends is 18.6cm, so shift left by half that.
           
           logicMylarStrawU[copyNo+32*j+64*i] = new G4LogicalVolume(MylarStraw,/*Vacuum*/mylar,"MylarStraw",0,0,0); // +32 * i+j up to moduleNum so I don't overwrite previous logic volume elements in arrar
           logicMylarStrawU[copyNo+32*j+64*i]->SetVisAttributes(StrawAtt);

           logicGasStrawU[copyNo+32*j+64*i] = new G4LogicalVolume(GasStraw,/*Vacuum*/trackerGas,"GasStraw",0,0,0); // +32 * i+j up to moduleNum so I don't overwrite previous logic volume elements in arrar
           logicWireStrawU[copyNo+32*j+64*i] = new G4LogicalVolume(WireStraw,/*Vacuum*/tungsten,"WireStraw",0,0,0); // +32 * i+j up to moduleNum so I don't overwrite previous logic volume elements in arrar

           new G4PVPlacement(URotation, G4ThreeVector(Xposition,Yposition,Zposition),"MylarStraw",logicMylarStrawU[copyNo+32*j+64*i],physiWorld,false,copyNo+32*j+64*i);       
           new G4PVPlacement(URotation, G4ThreeVector(Xposition,Yposition,Zposition),"GasStraw",logicGasStrawU[copyNo+32*j+64*i],physiWorld,false,copyNo+32*j+64*i);       
           new G4PVPlacement(URotation, G4ThreeVector(Xposition,Yposition,Zposition),"WireStraw",logicWireStrawU[copyNo+32*j+64*i],physiWorld,false,copyNo+32*j+64*i);   


           ////construction of all v layers///////////
           Xposition = 5.*cm + i * moduleSeparation + j * layerSeparation + UVSeparation;
           // Zposition = copyNo * (0.6*cm) + j * layerSeparation/2; //Would add shift of left or right here.

           logicMylarStrawV[copyNo+32*j+64*i] = new G4LogicalVolume(MylarStraw,/*Vacuum*/mylar,"MylarStraw",0,0,0); 
           logicMylarStrawV[copyNo+32*j+64*i]->SetVisAttributes(StrawAtt);

           logicGasStrawV[copyNo+32*j+64*i] = new G4LogicalVolume(GasStraw,/*Vacuum*/trackerGas,"GasStraw",0,0,0); // +32 * i+j up to moduleNum so I don't overwrite previous logic volume elements in arrar
           logicWireStrawV[copyNo+32*j+64*i] = new G4LogicalVolume(WireStraw,/*Vacuum*/tungsten,"WireStraw",0,0,0); // +32 * i+j up to moduleNum so I don't overwrite previous logic volume elements in arrar


           new G4PVPlacement(VRotation, G4ThreeVector(Xposition,Yposition,Zposition),"MylarStraw",logicMylarStrawV[copyNo+32*j+64*i],physiWorld,false,copyNo+32*j+64*i);            
           new G4PVPlacement(VRotation, G4ThreeVector(Xposition,Yposition,Zposition),"GasStraw",logicGasStrawV[copyNo+32*j+64*i],physiWorld,false,copyNo+32*j+64*i);       
           new G4PVPlacement(VRotation, G4ThreeVector(Xposition,Yposition,Zposition),"WireStraw",logicWireStrawV[copyNo+32*j+64*i],physiWorld,false,copyNo+32*j+64*i);


      }
    }
  }


  //------------------------------ 
  // Truth Detector Planes - planes that lie along wire planes to acquire truth information for traceback
  //------------------------------   

     G4VisAttributes * TruthVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
     TruthVisAtt->SetForceWireframe(true);

     logicTruthPlane = new G4LogicalVolume*[moduleNum*layerNum*2];

     for (int i = 0; i < moduleNum+1; ++i)
     {
           G4double Xposition = 5.*cm + (i-1) * moduleSeparation + (layerSeparation + UVSeparation)/2. + truthPlaneHalfX; // + truthPlaneHalfX so the front face of the box is at the center of the module.
           G4double Yposition = 0.;
           G4double Zposition = 0.;

           G4RotationMatrix* PlaneRotation = new G4RotationMatrix();
           PlaneRotation->rotateX(0.*deg);
           PlaneRotation->rotateY(0.*deg);
           PlaneRotation->rotateZ(0.*deg);

           G4ThreeVector planeTranslationPoint = G4ThreeVector(Xposition, Yposition, Zposition);
           G4Box* TruthPlane;

/////////////////////////////////////////////////////////////////////////////////////

           if (i == 0)// To move the 0 plane.
           { 
              double angle = 25.;

              Xposition = (-460.*std::sin(angle*pi/180))*cm; // These for a turned plane ~2 meters out in a circle for a 4.6 m radius for 2 GeV positrons, corresponding to an angle of 25 degrees.
              Zposition = (460.-460.*std::cos(angle*pi/180))*cm;
              planeTranslationPoint.setX(Xposition);
              planeTranslationPoint.setZ(Zposition);
              // PlaneRotation->rotateY(-angle*deg);

              TruthPlane = new G4Box("TruthPlane",truthPlaneHalfX,2.*truthPlaneHalfY,10.*truthPlaneHalfZ); 

           } 

/////////////////////////////////////////////////////////////////////////////////////

             else { 
            TruthPlane = new G4Box("TruthPlane",truthPlaneHalfX,truthPlaneHalfY,truthPlaneHalfZ); 
             }

           logicTruthPlane[i] = new G4LogicalVolume(TruthPlane,Vacuum,"TruthPlane",0,0,0);
           logicTruthPlane[i]->SetVisAttributes(TruthVisAtt);

           new G4PVPlacement(PlaneRotation, planeTranslationPoint /*G4ThreeVector(Xposition,Yposition,Zposition)*/,"TruthPlane",logicTruthPlane[i],physiWorld,false,i);

     }

     // for (int i = 0; i < moduleNum; ++i)
     // {
     //   for (int j = 0; j < layerNum; ++j)
     //   {
     //     for (int copyNo = 0; copyNo < 2; ++copyNo) // Loop for U or V placement.
     //     {
     //       G4double Xposition = 5.*cm + i * moduleSeparation + j * layerSeparation + copyNo * UVSeparation;

     //       G4Box* TruthPlane = new G4Box("TruthPlane",0.01*cm,HalfYZLength,HalfYZLength);

     //       logicTruthPlane[((2*copyNo)+j+(4*i))] = new G4LogicalVolume(TruthPlane,Vacuum,"TruthPlane",0,0,0); // ((2*copyNo)+j+(4*i)) is the pattern to get successive numbers for plane numbers, with 8 modules have 4 planes each.
     //       logicTruthPlane[((2*copyNo)+j+(4*i))]->SetVisAttributes(TruthVisAtt);

     //       new G4PVPlacement(0, G4ThreeVector(Xposition,0,0),"TruthPlane",logicTruthPlane[((2*copyNo)+j+(4*i))],physiWorld,false,((2*copyNo)+j+(4*i)));

     //     }
     //   }
     // }


  G4VisAttributes* worldVisAtt = new G4VisAttributes(0);
  // G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.));
  logicWorld->SetVisAttributes(worldVisAtt);
  return physiWorld;
}

////////////////////////////////////////
//Try to ignore field part here at first maybe and see if it works, otherwise maybe get rid of rest of field stuff and do it here.

void BasicDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  G4String trackerChamberSDname = "Basic/TrackerChamberSD";


  BasicTrackerSD* aTrackerSD = new BasicTrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");

  // Setting aTrackerSD to all logical volumes with the same name 
  // of "whatever Straw".
  SetSensitiveDetector("TruthPlane", aTrackerSD, true);

/*
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);*/
}

///////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BasicDetectorConstruction::SetMagField(G4double fieldValue)
{
  fMagField->SetFieldValue(fieldValue);
}

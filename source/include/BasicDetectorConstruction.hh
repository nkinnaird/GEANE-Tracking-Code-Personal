

#ifndef BasicDetectorConstruction_hh
#define BasicDetectorConstruction_hh 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "BasicMagneticField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
// class BasicDetectorMessenger;
class G4UserLimits;


//------------------------------------------------------------------------

class BasicDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  BasicDetectorConstruction();
  ~BasicDetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  void SetMagField(G4double);

  G4double GetHalfStrawLength() const {return HalfStrawLength;};
  G4double GetouterRadiusOfMylarStraw() const {return outerRadiusOfTheMylarStraw;};
  G4double GetdistBtwnWires() const {return distBtwnWires;};
  G4double GetmoduleSeparation() const {return moduleSeparation;};
  G4double GetlayerSeparation() const {return layerSeparation;};
  G4double GetUVSeparation() const {return UVSeparation;};
  G4double GetouterRadiusOfTheGas() const {return outerRadiusOfTheGas;};
  G4double GetouterRadiusOfTheWire() const {return outerRadiusOfTheWire;};
  G4double GetHalfYZLength() const {return HalfYZLength;};

  G4double GetTruthPlaneHalfX() const {return truthPlaneHalfX;};
  G4double GetTruthPlaneHalfY() const {return truthPlaneHalfY;};
  G4double GetTruthPlaneHalfZ() const {return truthPlaneHalfZ;};


private:
  // G4double HalfAirGapLength;

  G4double HalfStrawLength;
  G4double outerRadiusOfTheMylarStraw;
  G4double distBtwnWires;
  G4double moduleSeparation;
  G4double layerSeparation;
  G4double UVSeparation;
  G4double outerRadiusOfTheGas;
  G4double outerRadiusOfTheWire;

  G4double truthPlaneHalfX;
  G4double truthPlaneHalfY;
  G4double truthPlaneHalfZ;

  G4double HalfYZLength;
  G4double HalfWorldLength;
  G4UserLimits* fUserLimits;
  BasicMagneticField* fMagField;   // pointer to the magnetic field 

  G4LogicalVolume**  logicTruthPlane;
  
  G4LogicalVolume**  logicMylarStrawU;
  G4LogicalVolume**  logicMylarStrawV;
  
  G4LogicalVolume**  logicGasStrawU;
  G4LogicalVolume**  logicGasStrawV;
  
  G4LogicalVolume**  logicWireStrawU;
  G4LogicalVolume**  logicWireStrawV;
  
  // BasicDetectorMessenger* fDetectorMessenger;  // pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

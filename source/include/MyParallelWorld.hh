#ifndef MyParallelWorld_hh
#define MyParallelWorld_hh 1

#include "globals.hh"
#include "G4VUserParallelWorld.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;

class MyParallelWorld : public G4VUserParallelWorld
{
public:
	MyParallelWorld(G4String worldName);
	virtual ~MyParallelWorld();

	virtual void Construct();
	virtual void ConstructSD();

	G4double GetTruthPlaneHalfX() const {return truthPlaneHalfX;};
  	G4double GetTruthPlaneHalfY() const {return truthPlaneHalfY;};
  	G4double GetTruthPlaneHalfZ() const {return truthPlaneHalfZ;};

private:
	G4LogicalVolume**  logicTruthPlane;

	G4double truthPlaneHalfX;
  	G4double truthPlaneHalfY;
  	G4double truthPlaneHalfZ;

};

#endif

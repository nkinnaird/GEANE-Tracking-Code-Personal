

#ifndef BasicTrackerSD_h
#define BasicTrackerSD_h 1

#include "G4VSensitiveDetector.hh"

#include "BasicTrackerHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;
// class BasicEventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// BasicTracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class BasicTrackerSD : public G4VSensitiveDetector
{
  public:
    BasicTrackerSD(const G4String& name, 
                const G4String& hitsCollectionName);
    virtual ~BasicTrackerSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

    // BasicTrackerHitsCollection* GetBTCollection() const {return fHitsCollection;};

  private:
    BasicTrackerHitsCollection* fHitsCollection;
    // BasicEventAction*  fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

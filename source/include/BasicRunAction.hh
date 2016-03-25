

#ifndef BasicRunAction_h
#define BasicRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class HistoManager;

/// Run action class

class BasicRunAction : public G4UserRunAction
{
  public:
    BasicRunAction(HistoManager*);
    virtual ~BasicRunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void   EndOfRunAction(const G4Run* run);

  private:
  	HistoManager* fHistoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

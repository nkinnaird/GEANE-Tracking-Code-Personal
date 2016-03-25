

#ifndef BasicEventAction_h
#define BasicEventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"

/// Event action class

class BasicRunAction;
class HistoManager;

class BasicEventAction : public G4UserEventAction
{
  public:
    BasicEventAction(BasicRunAction*, HistoManager*);
    virtual ~BasicEventAction();

    virtual void  BeginOfEventAction(const G4Event* );
    virtual void    EndOfEventAction(const G4Event* );


  private:
  	// BasicRunAction* fRunAct;
  	HistoManager* fHistoManager;

  	G4int fEventID, fCopyNo;
  	G4double fGlobalTime, fProperTime, fXpos, fYpos, fZpos, fXmom, fYmom, fZmom, fUpos, fVpos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

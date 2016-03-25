

#ifndef BasicActionInitialization_h
#define BasicActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class BasicActionInitialization : public G4VUserActionInitialization
{
  public:
    BasicActionInitialization();
    virtual ~BasicActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

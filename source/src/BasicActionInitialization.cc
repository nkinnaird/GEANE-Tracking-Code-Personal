

#include "BasicActionInitialization.hh"
#include "BasicPrimaryGeneratorAction.hh"
#include "BasicRunAction.hh"
#include "BasicEventAction.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicActionInitialization::BasicActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicActionInitialization::~BasicActionInitialization()
{
	// delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicActionInitialization::BuildForMaster() const
{
  // Histo manager
  HistoManager*  histo = new HistoManager();
  
  // Actions
  SetUserAction(new BasicRunAction(histo));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicActionInitialization::Build() const
{
  SetUserAction(new BasicPrimaryGeneratorAction);

  HistoManager*  histo = new HistoManager();

  BasicRunAction* runAction = new BasicRunAction(histo);  
  SetUserAction(runAction);
  
  BasicEventAction* eventAction = new BasicEventAction(runAction, histo);
  SetUserAction(eventAction);
  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#include "BasicRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicRunAction::BasicRunAction(HistoManager* histo)
 : G4UserRunAction(), fHistoManager(histo)
{ 
  // set printing event number per each 100 events
  // G4RunManager::GetRunManager()->SetPrintProgress(1000);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicRunAction::~BasicRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
   G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

   fHistoManager->book(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicRunAction::EndOfRunAction(const G4Run* /*aRun*/)
{

	fHistoManager->save();  
	// also need to consider deleting histo somewhere

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

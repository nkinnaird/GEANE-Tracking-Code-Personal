

#include "BasicEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4VHit.hh"
#include "G4THitsCollection.hh"

// #include "BasicTrackerHit.hh"
// #include "BasicTrackerSD.hh"

#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"

#include "G4UIcommand.hh"

#include "HistoManager.hh"
#include "BasicRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicEventAction::BasicEventAction(BasicRunAction* /*run*/, HistoManager* histo)
: G4UserEventAction(), /*fRunAct(run),*/ fHistoManager(histo),
  fEventID(0.), fCopyNo(0.),
  fGlobalTime(0.), fProperTime(0.), 
  fXpos(0.), fYpos(0.), fZpos(0.), 
  fXmom(0.), fYmom(0.), fZmom(0.),
  fUpos(0.), fVpos(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicEventAction::~BasicEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicEventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BasicEventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories <-- uneccessary and unwanted since this will include secondary trajectories

  // G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  // G4int n_trajectories = 0;
  // if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  fEventID = event->GetEventID();
  // G4cout << "Event ID: " << eventID << G4endl;

  // if ( eventID < 100 || eventID % 100 == 0) {
    // G4cout << ">>> Event: " << eventID  << G4endl;
    // if ( trajectoryContainer ) {
      // G4cout << "    " << n_trajectories
             // << " trajectories stored in this event." << G4endl;
    // }

    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4int size = hc->GetSize();

    G4cout << "    "  << size << " hits stored in this event" << G4endl;

    for ( G4int i=0; i<size; i++ ){
      G4cout << G4endl << "EventID: " << fEventID << std::setw(7);
      (hc->GetHit(i))->Print();
      G4cout  << G4endl;

      // This for accessing the map which refers to the ATTRIBUTES of the BasicTrackerHit variables, not the values themselves.
      
      // const std::map<G4String,G4AttDef>* store = (hc->GetHit(i))->GetAttDefs();
      // G4cout << "This is where I am printing some stuff. " << store->at(G4String("DetPos")).GetName() << G4endl;
      
      // This for accessing the values of the BasicTrackerHit variables.
      // Not sure if there is a better way to go about this, but values 0, 1, 2, 3 refer to Global Time, Proper Time, Position vector, and Momentum vector respectively.
      // Value 4 holds the CopyNo's with units of time to get things to work.
      // Value 5 and 6 hold U and V position doubles respectively.

      std::vector<G4AttValue>* values = (hc->GetHit(i))->CreateAttValues();
      
      /*G4cout << "Check to see if things are reprinted correctly. " << G4endl;
      for (G4int j = 0; j < 4; ++j)
      {
        G4cout << values->at(j).GetName() << " " << values->at(j).GetValue() << G4endl;
      }*/

      fGlobalTime = G4UIcommand::ConvertToDimensionedDouble(values->at(0).GetValue());
      G4cout << "Global time in ns: " << fGlobalTime << G4endl;

      fProperTime = G4UIcommand::ConvertToDimensionedDouble(values->at(1).GetValue());
      G4cout << "Proper time in ns: " << fProperTime << G4endl;

      G4ThreeVector XYZPos = G4UIcommand::ConvertToDimensioned3Vector(values->at(2).GetValue());
      G4cout << "Position in mm: " << XYZPos << G4endl;

      G4ThreeVector XYZMom = G4UIcommand::ConvertToDimensioned3Vector(values->at(3).GetValue());
      G4cout << "Momentum in MeV: " << XYZMom << G4endl;

      fXpos = XYZPos.getX();
      fYpos = XYZPos.getY();
      fZpos = XYZPos.getZ();

      fXmom = XYZMom.getX();
      fYmom = XYZMom.getY();
      fZmom = XYZMom.getZ();

      fCopyNo = G4UIcommand::ConvertToDimensionedDouble(values->at(4).GetValue());

      fUpos = G4UIcommand::ConvertToDimensionedDouble(values->at(5).GetValue());
      fVpos = G4UIcommand::ConvertToDimensionedDouble(values->at(6).GetValue());
      G4cout << "U Position: " << fUpos << " V Position: " << fVpos << G4endl;


      fHistoManager->FillNtuple(fEventID, fGlobalTime, fProperTime, fXpos, fYpos, fZpos, fXmom, fYmom, fZmom, fCopyNo, fUpos, fVpos);

      // Using this "Dimensioned" style here so that units stay the same from one hit to another - otherwise it takes the "Best" unit
      // changing from ps to ns, cm to m, etc as thresholds are passed.
      

    } 


  //}
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

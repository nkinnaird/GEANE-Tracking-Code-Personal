

#include "BasicTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
// #include "G4UIcommand.hh"
#include "G4AttCheck.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<BasicTrackerHit>* BasicTrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicTrackerHit::BasicTrackerHit()
 : G4VHit(),
   fTrackID(-1),
   fGlobalTime(0),
   fProperTime(0),
   fPos(G4ThreeVector()),
   fMom(G4ThreeVector()),
   fCopyNo(-1),
   fUpos(0), fVpos(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicTrackerHit::~BasicTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BasicTrackerHit::BasicTrackerHit(const BasicTrackerHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fGlobalTime = right.fGlobalTime;
  fProperTime = right.fProperTime;
  fPos       = right.fPos;
  fMom       = right.fMom;
  fCopyNo    = right.fCopyNo;
  fUpos      = right.fUpos;
  fVpos      = right.fVpos;
}

// BasicTrackerHit::BasicTrackerHit(const G4VHit& bhit){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const BasicTrackerHit& BasicTrackerHit::operator=(const BasicTrackerHit& right)
{
  fTrackID   = right.fTrackID;
  fGlobalTime = right.fGlobalTime;
  fProperTime = right.fProperTime;
  fPos       = right.fPos;
  fMom       = right.fMom;
  fCopyNo    = right.fCopyNo;
  fUpos      = right.fUpos;
  fVpos      = right.fVpos;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int BasicTrackerHit::operator==(const BasicTrackerHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void BasicTrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  const std::map<G4String,G4AttDef>* BasicTrackerHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("BasicTrackerHit",isNew);
  if (isNew) {

    G4String GTime("GTime");
    (*store)[GTime] = G4AttDef(GTime, "Global Time", "Physics", "G4BestUnit", "G4double");

    G4String PTime("PTime");
    (*store)[PTime] = G4AttDef(PTime, "Proper Time", "Physics", "G4BestUnit", "G4double");

    G4String DetPos("DetPos");
    (*store)[DetPos] = G4AttDef(DetPos, "Position on Detector", "Physics", "G4BestUnit", "G4ThreeVector");

    G4String DetMom("DetMom");
    (*store)[DetMom] = G4AttDef(DetMom, "Momentum on Detector", "Physics", "G4BestUnit", "G4ThreeVector");

    G4String CopyNo("CopyNo");
    (*store)[CopyNo] = G4AttDef(CopyNo, "Copy Number", "Physics", "G4BestUnit","G4double"); // Defining copy No in terms of time units since I'm not sure what the proper category/type key/value type to use...

    G4String UPosition("UPosition");
    (*store)[UPosition] = G4AttDef(UPosition, "U Position", "Physics", "G4BestUnit", "G4double");

    G4String VPosition("VPosition");
    (*store)[VPosition] = G4AttDef(VPosition, "V Position", "Physics", "G4BestUnit", "G4double");


  }
  return store;
}

std::vector<G4AttValue>* BasicTrackerHit::CreateAttValues() const
{ // Make sure these push_back statements align with the proper .at() values in BasicEventAction.cc
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("GTime",G4BestUnit(fGlobalTime,"Time"),""));

  values->push_back
    (G4AttValue("PTime",G4BestUnit(fProperTime,"Time"),""));

  values->push_back
    (G4AttValue("DetPos",G4BestUnit(fPos,"Length"),""));

  values->push_back
    (G4AttValue("DetMom",G4BestUnit(fMom,"Energy"),""));

  values->push_back
    (G4AttValue("CopyNo",G4BestUnit(fCopyNo,"Time"),""));

  values->push_back
    (G4AttValue("UPosition",G4BestUnit(fUpos,"Length"),""));  

  values->push_back
    (G4AttValue("VPosition",G4BestUnit(fVpos,"Length"),""));  


// #ifdef G4ATTDEBUG
  // G4cout << G4AttCheck(values,GetAttDefs());
// #endif

  return values;
}



void BasicTrackerHit::Print()
{
  // G4cout
  //    << "  trackID: " << fTrackID 
  //    // << " GlobalTime: " << fGlobalTime
  //    << " GlobalTime: " << std::setw(7) << G4BestUnit(fGlobalTime, "Time") 
  //    << " ProperTime: " << std::setw(7) << G4BestUnit(fProperTime, "Time")
  //    << " Position: " << std::setw(7) << G4BestUnit(fPos, "Length")
  //    << " Momentum: " << std::setw(7) << G4BestUnit(fMom, "Energy")
  //    // << " Postion: " << fPos
  //    // << " Momentum: " << fMom
  //    << " copyNo: " << fCopyNo
  //    << " UPos: " << fUpos
  //    << " VPos: " << fVpos
  //    << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

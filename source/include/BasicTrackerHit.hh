

#ifndef BasicTrackerHit_h
#define BasicTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///


class BasicTrackerHit : public G4VHit
{
  public:
    BasicTrackerHit();
    BasicTrackerHit(const BasicTrackerHit&);
    virtual ~BasicTrackerHit();

    // operators
    const BasicTrackerHit& operator=(const BasicTrackerHit&);
    G4int operator==(const BasicTrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    // virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID    (G4int track)           { fTrackID = track; };
    void SetGlobalTime (G4double gTime)        { fGlobalTime = gTime; };
    void SetProperTime (G4double pTime)        { fProperTime = pTime; };
    void SetPos        (G4ThreeVector xyz)     { fPos = xyz; };
    void SetMom        (G4ThreeVector pxpypz)  { fMom = pxpypz; };
    void SetCopyNo     (G4int copyNo)          { fCopyNo = copyNo; };
    void SetUPos       (G4double Upos)         { fUpos = Upos; };
    void SetVPos       (G4double Vpos)         { fVpos = Vpos; };

    // Get methods
    G4int     GetTrackID()    const   { return fTrackID; };
    G4double  GetGlobalTime() const   { return fGlobalTime; };
    G4double  GetProperTime() const   { return fProperTime; };
    G4ThreeVector GetPos()    const   { return fPos; };
    G4ThreeVector GetMom()    const   { return fMom; };
    G4int GetCopyNo()         const   { return fCopyNo; };
    G4double GetUPos()        const   { return fUpos; };
    G4double GetVPos()        const   { return fVpos; };

    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;

  private:

      G4int         fTrackID;
      G4double      fGlobalTime;
      G4double      fProperTime;
      G4ThreeVector fPos;
      G4ThreeVector fMom;
      G4int         fCopyNo;
      G4double      fUpos;
      G4double      fVpos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<BasicTrackerHit> BasicTrackerHitsCollection;

extern G4ThreadLocal G4Allocator<BasicTrackerHit>* BasicTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* BasicTrackerHit::operator new(size_t)
{
  if(!BasicTrackerHitAllocator)
      BasicTrackerHitAllocator = new G4Allocator<BasicTrackerHit>;
  return (void *) BasicTrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void BasicTrackerHit::operator delete(void *hit)
{
  BasicTrackerHitAllocator->FreeSingle((BasicTrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

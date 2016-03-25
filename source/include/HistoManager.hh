
#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 class TFile;
 class TTree;
 class TH1D;

  // const G4int MaxHisto = 5;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
  
    HistoManager();
   ~HistoManager();
   
    void book();
    void save();

    // void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
    // void Normalize(G4int id, G4double fac);    

    void FillNtuple(G4int IeventID, G4double IGlobalTime, G4double IProperTime,
                              G4double Ixpos, G4double Iypos, G4double Izpos, 
                              G4double Ixmom, G4double Iymom, G4double Izmom,
                              G4int IcopyNo, G4double Iupos, G4double Ivpos);
    
    // void PrintStatistic();
        
  private:
  
    TFile*   fRootFile;
    // TH1D*    fHisto[MaxHisto];            
    TTree*   fNtuple1;    
    // TTree*   fNtuple2;    

   /* G4double fEabs;
    G4double fEgap;
    G4double fLabs;
    G4double fLgap;*/

    G4int eventId;
    G4double GlobalTime;
    G4double ProperTime;
    G4double xpos;
    G4double ypos;
    G4double zpos;
    G4double xmom;
    G4double ymom;
    G4double zmom;
    G4int copyNo;
    G4double upos;
    G4double vpos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif




#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:fRootFile(0), 
     fNtuple1(0), 
     eventId(0),
     GlobalTime(0),
     ProperTime(0),
     xpos(0),
     ypos(0),
     zpos(0),
     xmom(0),
     ymom(0),
     zmom(0),
     copyNo(0),
     upos(0),
     vpos(0)
{
      
  // histograms
  // for (G4int k=0; k<MaxHisto; k++) fHisto[k] = 0;
    
  // ntuple
  fNtuple1 = 0;
  // fNtuple2 = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  if ( fRootFile ) delete fRootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
 G4String fileName = "basicntuple.root";
 fRootFile = new TFile(fileName,"RECREATE");
 if(!fRootFile) {
   G4cout << " HistoManager::book :" 
          << " problem creating the ROOT TFile "
          << G4endl;
   return;
 }
   
 /*fHisto[1] = new TH1D("1", "Edep in absorber", 100, 0., 800*CLHEP::MeV);
 if (!fHisto[1]) G4cout << "\n can't create histo 1" << G4endl;
 fHisto[2] = new TH1D("2", "Edep in gap", 100, 0., 100*CLHEP::MeV);
 if (!fHisto[2]) G4cout << "\n can't create histo 2" << G4endl;
 fHisto[3] = new TH1D("3", "trackL in absorber", 100, 0., 1*CLHEP::m);
 if (!fHisto[3]) G4cout << "\n can't create histo 3" << G4endl;
 fHisto[4] = new TH1D("4", "trackL in gap", 100, 0., 50*CLHEP::cm);
 if (!fHisto[4]) G4cout << "\n can't create histo 4" << G4endl;  */



 // create 1st ntuple in subdirectory "tuples"
 //
 fNtuple1 = new TTree("TrackInfo", "TrackingInformation");
 // fNtuple1->Branch("TrackingInformation", &fEabs, "Eabs/D"); <--- changes here
 fNtuple1->Branch("EventID", &eventId, "EventID/I"); // Last part seems to do with doubles/ints, etc.
 fNtuple1->Branch("GlobalTime", &GlobalTime, "GlobalTime/D"); 
 fNtuple1->Branch("ProperTime", &ProperTime, "ProperTime/D"); 
 fNtuple1->Branch("XPosition", &xpos, "XPosition/D"); 
 fNtuple1->Branch("YPosition", &ypos, "YPosition/D"); 
 fNtuple1->Branch("ZPosition", &zpos, "ZPosition/D"); 
 fNtuple1->Branch("XMomentum", &xmom, "XMomentum/D"); 
 fNtuple1->Branch("YMomentum", &ymom, "YMomentum/D"); 
 fNtuple1->Branch("ZMomentum", &zmom, "ZMomentum/D"); 
 fNtuple1->Branch("CopyNo", &copyNo, "CopyNo/I");
 fNtuple1->Branch("UPosition", &upos, "UPosition/D"); 
 fNtuple1->Branch("VPosition", &vpos, "VPosition/D"); 




 // fNtuple1->Branch("Egap", &fEgap, "Egap/D");

 // create 2nd ntuple in subdirectory "tuples"
 //
 /*fNtuple2 = new TTree("102", "TrackL");
 fNtuple2->Branch("Labs", &fLabs, "Labs/D");
 fNtuple2->Branch("Lgap", &fLgap, "Lgap/D");*/

 
 G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
  if (fRootFile) {
    fRootFile->Write();       // Writing the histograms to the file
    fRootFile->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " does not exist. (xbin=" << xbin << " weight=" << weight << ")"
           << G4endl;
    return;
  }
 if  (fHisto[ih]) { fHisto[ih]->Fill(xbin, weight); }
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }
  if (fHisto[ih]) fHisto[ih]->Scale(fac);
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              // G4double trackLAbs , G4double trackLGap )
void HistoManager::FillNtuple(G4int IeventID, G4double IGlobalTime, G4double IProperTime,
                              G4double Ixpos, G4double Iypos, G4double Izpos, 
                              G4double Ixmom, G4double Iymom, G4double Izmom,
                              G4int IcopyNo, G4double Iupos, G4double Ivpos)

{
 /*fEabs = energyAbs;
 fEgap = energyGap;
 fLabs = trackLAbs;
 fLgap = trackLGap;*/

    eventId = IeventID;
    GlobalTime = IGlobalTime;
    ProperTime = IProperTime;
    xpos = Ixpos;
    ypos = Iypos;
    zpos = Izpos;
    xmom = Ixmom;
    ymom = Iymom;
    zmom = Izmom;
    copyNo = IcopyNo;
    upos = Iupos;
    vpos = Ivpos;

  if (fNtuple1) fNtuple1->Fill(); // <-- pass parameters here
  // if (fNtuple2) fNtuple2->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void HistoManager::PrintStatistic()
{
  if(fHisto[1]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " EAbs : mean = " << G4BestUnit(fHisto[1]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(fHisto[1]->GetRMS(),  "Energy") << G4endl;
    G4cout                
    << " EGap : mean = " << G4BestUnit(fHisto[2]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(fHisto[2]->GetRMS(),  "Energy") << G4endl;
    G4cout 
    << " LAbs : mean = " << G4BestUnit(fHisto[3]->GetMean(), "Length") 
            << " rms = " << G4BestUnit(fHisto[3]->GetRMS(),  "Length") << G4endl;
    G4cout 
    << " LGap : mean = " << G4BestUnit(fHisto[4]->GetMean(), "Length") 
            << " rms = " << G4BestUnit(fHisto[4]->GetRMS(),  "Length") << G4endl;

  }
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



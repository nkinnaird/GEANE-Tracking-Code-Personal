//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jun 14 10:58:48 2015 by ROOT version 5.34/30
// from TTree TrackInfo/TrackingInformation
// found on file: basicntuple.root
//////////////////////////////////////////////////////////

#ifndef RootInput_h
#define RootInput_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class RootInput {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventID;
   Double_t        GlobalTime;
   Double_t        ProperTime;
   Double_t        XPosition;
   Double_t        YPosition;
   Double_t        ZPosition;
   Double_t        XMomentum;
   Double_t        YMomentum;
   Double_t        ZMomentum;
   Int_t           CopyNo;
   Double_t        UPosition;
   Double_t        VPosition;

   // Strings for map keys to pass back to main program. 
   std::string EventID_k   ;
   std::string GlobalTime_k;
   std::string ProperTime_k;
   std::string XPosition_k ;
   std::string YPosition_k ;
   std::string ZPosition_k ;
   std::string XMomentum_k ;
   std::string YMomentum_k ;
   std::string ZMomentum_k ;
   std::string CopyNo_k    ;
   std::string UPosition_k ;
   std::string VPosition_k ;

   // List of branches
   TBranch        *b_EventID;   //!
   TBranch        *b_GlobalTime;   //!
   TBranch        *b_ProperTime;   //!
   TBranch        *b_XPosition;   //!
   TBranch        *b_YPosition;   //!
   TBranch        *b_ZPosition;   //!
   TBranch        *b_XMomentum;   //!
   TBranch        *b_YMomentum;   //!
   TBranch        *b_ZMomentum;   //!
   TBranch        *b_CopyNo; // These branches I added myself later on. The others were auto generated.
   TBranch        *b_UPosition; 
   TBranch        *b_VPosition;

   RootInput(TTree *tree=0);
   virtual ~RootInput();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual std::map<std::string, std::vector<double> > LoopAndFill(std::map<std::string, std::vector<double> >);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootInput_cxx
RootInput::RootInput(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("basicntuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("basicntuple.root");
      }
      f->GetObject("TrackInfo",tree);

   }
   Init(tree);
}

RootInput::~RootInput()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootInput::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootInput::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RootInput::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("GlobalTime", &GlobalTime, &b_GlobalTime);
   fChain->SetBranchAddress("ProperTime", &ProperTime, &b_ProperTime);
   fChain->SetBranchAddress("XPosition", &XPosition, &b_XPosition);
   fChain->SetBranchAddress("YPosition", &YPosition, &b_YPosition);
   fChain->SetBranchAddress("ZPosition", &ZPosition, &b_ZPosition);
   fChain->SetBranchAddress("XMomentum", &XMomentum, &b_XMomentum);
   fChain->SetBranchAddress("YMomentum", &YMomentum, &b_YMomentum);
   fChain->SetBranchAddress("ZMomentum", &ZMomentum, &b_ZMomentum);
   fChain->SetBranchAddress("CopyNo", &CopyNo, &b_CopyNo);
   fChain->SetBranchAddress("UPosition", &UPosition, &b_UPosition);
   fChain->SetBranchAddress("VPosition", &VPosition, &b_VPosition);
   Notify();
}

Bool_t RootInput::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootInput::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootInput::Cut(Long64_t /*entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootInput_cxx

#define RootInput_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "TVector3.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TInterpreter.h"
#include <iostream>
#include <fstream>

#include "RootInput.hh"

std::map<std::string, std::vector<double> > RootInput::LoopAndFill(std::map<std::string, std::vector<double> > InputParameters)
{
//   In a ROOT session, you can do:
//      Root > .L RootInput.C
//      Root > RootInput t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   EventID_k = "EventID";
   GlobalTime_k = "GlobalTime";
   ProperTime_k = "ProperTime";
   XPosition_k = "XPosition";
   YPosition_k = "YPosition";
   ZPosition_k = "ZPosition";
   XMomentum_k = "XMomentum";
   YMomentum_k = "YMomentum";
   ZMomentum_k = "ZMomentum";
   CopyNo_k = "CopyNo";
   UPosition_k = "UPosition";
   VPosition_k = "VPosition";


   if (fChain == 0)
   {
   	std::cout << "Tree loading messed up somehow." << std::endl;
   	return InputParameters;
   } 

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // RootInput::Show(jentry);

      // Everytime GetEntry() is called, the variables linked by SetBranchAddress in the header file are filled accordingly.
      // Those are then put at the end of the corresponding vectors in the map that was passed in to this program.
      
      InputParameters[EventID_k].push_back(double(EventID));
      InputParameters[GlobalTime_k].push_back(double(GlobalTime));
      InputParameters[ProperTime_k].push_back(double(ProperTime));
      InputParameters[XPosition_k].push_back(double(XPosition));
      InputParameters[YPosition_k].push_back(double(YPosition));
      InputParameters[ZPosition_k].push_back(double(ZPosition));
      InputParameters[XMomentum_k].push_back(double(XMomentum));
      InputParameters[YMomentum_k].push_back(double(YMomentum));
      InputParameters[ZMomentum_k].push_back(double(ZMomentum));
      InputParameters[CopyNo_k].push_back(double(CopyNo));
      InputParameters[UPosition_k].push_back(double(UPosition));
      InputParameters[VPosition_k].push_back(double(VPosition));

   }

   return InputParameters;


}

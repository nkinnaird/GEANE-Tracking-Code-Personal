// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// ------------------------------------------------------------
//      GEANT 4 example main
// ------------------------------------------------------------
//
// History:
// - Created:   P. Arce   May 2007
//

// #include "ExErrorDetectorConstruction.hh"
#include "G4SteppingVerbose.hh"

#include "G4ErrorPropagator.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagatorManager.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorCylSurfaceTarget.hh"
#include "G4ErrorGeomVolumeTarget.hh"
#include "G4ErrorTrackLengthTarget.hh"
#include "G4ErrorFreeTrajState.hh"

#include "G4UImanager.hh"

#include "TRandom3.h"
#include "TMath.h"
#define NSTOP 20
#define READ 1 
#define WRITE 0



void Initialize();
G4ErrorTarget* BuildTarget( G4int iTarget );
void ProcessEvent( G4int iProp, size_t nEv, TFile* rootFile );
void Finalize();
void trycorr(int jev, G4ErrorMatrix* g4em, double xm[][20], double xp[][20],double dtvec[]);


G4ErrorTarget* theTarget;
G4ErrorTarget* theTarget1;
G4ErrorTarget* theTarget2;
G4ErrorMode theG4ErrorMode;
G4ErrorMode theG4ErrorMode1;
G4ErrorMode theG4ErrorMode2;

G4ErrorPropagatorManager* g4emgr;
G4ErrorPropagatorData* g4edata; 

#include <vector>
#include "TVector3.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TInterpreter.h"
#include <iostream>
#include <fstream>
using namespace std;


//-------------------------------------------------------------
int main() 
{

  //JPAdd
  //Here, we initialize the root objects to store the desired data.

  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#pragma link C++ class vector<vector<vector<float>>>+");
  gROOT->ProcessLine("#pragma link C++ class vector<TVector3>+");

  gInterpreter->GenerateDictionary("vector<TVector3>","TVector3.h");


  TFile *rootFile = new TFile("ErrPropLibrary.root", "RECREATE");
  rootFile->cd();

  TTree *PosTree = new TTree("PosTree","PosTree");
  TTree *MomTree = new TTree("MomTree","MomTree");
  TTree *ErrMatTree = new TTree("ErrMatTree","ErrMatTree");

  /////////////////////////////////////////////////////////

  Initialize();

  //----- Choose propagation mode
  // 0: propagate until target, all steps in one go
  // 1: propagate until target, returning control to the user at each step
  G4int iProp = 0;
  char* prop = getenv("G4ERROR_PROP");
  if( prop ) {
    if( G4String(prop) == G4String("UNTIL_TARGET") ){
      iProp = 0;
    } else if ( G4String(prop) == G4String("STEP_BY_STEP") ) {
      iProp = 1;
    } else {
      G4Exception("exG4eReco","Fatal error in Argument",FatalErrorInArgument,G4String("Variable G4ERROR_PROP = " + G4String(prop) + "   It must be: UNTIL_TARGET or STEP_BY_STEP").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",JustWarning,"Variable G4ERROR_PROP not defined, taking it = UNTIL_TARGET");
  } 

#if READ
  size_t nEvents = 1;
#endif
#if WRITE
  size_t nEvents = 1;
#endif
  for( size_t ii = 0; ii < nEvents; ii++ ){
    ProcessEvent( iProp, ii, rootFile );
  }
  printf("Call Finalize \n");
  Finalize();
  printf("Returned from Finalize \n");
  printf("End of main Routine \n");

  //JPAdd
  rootFile->cd();
  PosTree->Write();
  MomTree->Write();
  ErrMatTree->Write();

  rootFile->Close();
  
}


//-------------------------------------------------------------
void Initialize() 
{
  G4VSteppingVerbose::SetInstance(new G4SteppingVerbose);

  // Initialize the GEANT4e manager 
  //  G4ErrorPropagatorManager* g4emgr = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  //G4ErrorPropagatorData* g4edata = G4ErrorPropagatorData::GetErrorPropagatorData();
  g4emgr = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  g4edata = G4ErrorPropagatorData::GetErrorPropagatorData();

  g4emgr->SetUserInitialization(new ExErrorDetectorConstruction); 

  G4UImanager::GetUIpointer()->ApplyCommand("/exerror/setField -10. kilogauss");

  g4emgr->InitGeant4e();

  G4UImanager::GetUIpointer()->ApplyCommand("/control/verbose 1");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 100 mm");

  //----- Choose target type:
  // 1: PlaneSurfaceTarget
  // 2: CylSurfaceTarget
  // 3: GeomVolumeTarget
  // 4: TrackLengthTarget
  G4int iTarget = 1;
  char* target = getenv("G4ERROR_TARGET");
  if( target ) {
    if( G4String(target) == G4String("PLANE_SURFACE") ) {
      iTarget = 1;
    }else if( G4String(target) == G4String("CYL_SURFACE") ) {
      iTarget = 2;
    }else if( G4String(target) == G4String("VOLUME") ) {
      iTarget = 3;
    }else if( G4String(target) == G4String("TRKLEN") ) {
      iTarget = 4;
    }else {
      G4Exception("exG4eReco","Fatal error in Argument",FatalErrorInArgument,G4String("Variable G4ERROR_TARGET = " + G4String(target) + "   It must be:  PLANE_SURFACE, CYL_SURFACE, VOLUME, TRKLEN").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",JustWarning,"Variable G4ERROR_TARGET not defined, taking it = PLANE_SURFACE");
  } 
  iTarget=1;
  theTarget = BuildTarget( iTarget );
  iTarget=1;
  theTarget1 = BuildTarget( iTarget );
  iTarget=2;
  theTarget2 = BuildTarget( iTarget );
  g4edata->SetTarget( theTarget );

  theG4ErrorMode = G4ErrorMode_PropBackwards;
  char* mode = getenv("G4ERROR_MODE");
#if 0 
 if( mode ) {
    if( G4String(mode) == G4String("FORWARDS") ) {
      theG4ErrorMode = G4ErrorMode_PropForwards;
    } else if( G4String(mode) == G4String("BACKWARDS") ) {
      theG4ErrorMode = G4ErrorMode_PropBackwards;
    } else {
      G4Exception("exG4eReco","Fatal error in Argument",FatalErrorInArgument,G4String("Variable G4ERROR_MODE = " + G4String(mode) + "   It must be:  FORWARDS or BACKWARDS").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",JustWarning,"Variable G4ERROR_MODE not defined, taking it = BACKWARDS");
  } 
#endif
theG4ErrorMode =  G4ErrorMode_PropForwards;
      printf("Initialize: The G4 Error Mode is %d\n",theG4ErrorMode);      

}


void ProcessEvent( G4int iProp, size_t iev, TFile *rootFile )
{

  static TRandom3 r3;
  static double savex[20][3];
  static double savep[20][3];
  double t5tp[5][20];  //predicted track parameters
  static double t5tm[5][20];  //measured track parameters

  
  static double t5orig0[5];  
  static double t5final0[5];  
  static double t5orig1[5];  
  static double t5final1[5];  
  static double dt5o[5];
  static double dt5f[5];
  static double dt5fch[5];
  double rvec[20][3];
  double pvec[20][3];

  static double xposlast[3];
  static double xpos[3];

  double dl;
  double dtvec[5];

  G4ErrorMatrix G4EM[50];
  G4ErrorMatrix G4EMT[50];
  G4ErrorMatrix SigWireSqM2 = G4ErrorMatrix(5,5,0);
  SigWireSqM2(4,4) = 100.0; //inv mm squared (or 100 cm^-2 which goes to 1 mm^-2)
  SigWireSqM2(5,5) = 100.0;
// Set the starting trajectory.
  G4ThreeVector xv3,pv3;
  G4ThreeVector xv3save,pv3save;
  if(iev==0) {
    r3.SetSeed(3547459);
  }



#if  READ
  ifstream myfile;
#endif
#if  WRITE
  ofstream myfile;
#endif

  myfile.open("examplebacktrack.dat");
#if READ
  int dummy;
  if(iev==0) {

  for(int ipl=19;ipl>-1;ipl--) {
    myfile >> dummy >> rvec[ipl][0]  >> rvec[ipl][1] >> rvec[ipl][2] >>
	   pvec[ipl][0]  >> pvec[ipl][1] >> pvec[ipl][2];
    printf("plane %d position %f %f %f momentum %f %f %f \n", ipl, rvec[ipl][0], rvec[ipl][1], rvec[ipl][2],   pvec[ipl][0], pvec[ipl][1], pvec[ipl][2]);
    //for backwards tracking, reverse the momenta
    rvec[ipl][1] += 0.002*(r3.Rndm()-0.5); //+- 100 microns
    rvec[ipl][2] += 0.002*(r3.Rndm()-0.5); //+- 100 microns
    pvec[ipl][0] = -pvec[ipl][0];
    pvec[ipl][1] = -pvec[ipl][1];
    pvec[ipl][2] = -pvec[ipl][2];
    t5tm[0][ipl] = 0.0;
    t5tm[1][ipl] = 0.0;
    t5tm[2][ipl] = 0.0;
    t5tm[3][ipl] = rvec[ipl][1];
    t5tm[4][ipl] = rvec[ipl][2];
    savex[ipl][0] = rvec[ipl][0]; 
    savex[ipl][1] = rvec[ipl][1];
    savex[ipl][2] = rvec[ipl][2];
    savep[ipl][0] = pvec[ipl][0];
    savep[ipl][1] = pvec[ipl][1];
    savep[ipl][2] = pvec[ipl][2];

  }
  }
#if 0
  for(int ipl = 1;ipl<20;ipl++) {
    rvec[ipl][1] += 0.2*(r3.Rndm()-0.5); //+- 100 microns
    rvec[ipl][2] += 0.2*(r3.Rndm()-0.5); //+- 100 microns
  }
#endif
  //If this is the first event we set the 
  //  if(iev !=1) {
    xv3.set( rvec[0][0],rvec[0][1],rvec[0][2]); // starting trajectories set from 0th plane position values
#if READ
    //    pv3.set( 3.0*GeV, 0.0, 0.1-0.01*iev );
    //    pv3.set( pvec[0][0],pvec[0][1],pvec[0][2]);
        pv3.set( pvec[0][0]+300.0,pvec[0][1],pvec[0][2]); // starting momenta set from 0th plane values with additional x momenta
#endif
#if WRITE // WRITE set as 0 when Rob sent it to me
    pv3.set( 3.0*GeV, 0.0, 0.0 );
#endif


    //  theG4ErrorMode = G4ErrorMode_PropForwards;  
  theG4ErrorMode = G4ErrorMode_PropBackwards;  
  G4ErrorTrajErr error( 5, 0 );
  error(4,4) = 100.0; // Maybe from 100 cm^-2 which goes to 1 mm ^ -2
  error(5,5) = 100.0; // For inital unity y z (perp) errors

  //Initialize the state
  G4ErrorFreeTrajState* theG4ErrorTrajState = new G4ErrorFreeTrajState( "mu-", xv3, pv3, error );

  G4ErrorPropagatorManager* g4emgr1 = G4ErrorPropagatorManager::GetErrorPropagatorManager();

  int ierr = 0;

  //  G4Point3D surfPos(224.*cm,0.,0.);
  //Print out information on initial values for track

  G4Point3D posBeg = theG4ErrorTrajState->GetPosition();
  G4Normal3D momBeg = theG4ErrorTrajState->GetMomentum();
        G4ErrorFreeTrajParam* theG4TrackParam1 = new G4ErrorFreeTrajParam( posBeg,momBeg); // can also call getparameters method here I believe from G4errorfreetrajstate
        printf("Event number %d \n",iev);
        printf("Init Inv P value %f \n",theG4TrackParam1->GetInvP());
        printf("Init Lambda value %f \n",theG4TrackParam1->GetLambda());
        printf("Init Phi value %f \n",theG4TrackParam1->GetPhi());
        printf("Init YPerp value %f \n",theG4TrackParam1->GetYPerp());
        printf("Init Zperp value %f \n",theG4TrackParam1->GetZPerp());
      t5tp[0][0] = theG4TrackParam1->GetInvP();
      t5tp[1][0] = theG4TrackParam1->GetLambda();
      t5tp[2][0] = theG4TrackParam1->GetPhi();
      t5tp[3][0] = posBeg[1];
      t5tp[4][0] = posBeg[2];


	//load up t5orig
  xposlast[0] = posBeg[0];
  xposlast[1] = posBeg[1];
  xposlast[2] = posBeg[2];


  iProp=1;
  G4cout << "The track before" << theG4ErrorTrajState->GetParameters() << G4endl; 


#endif


  printf("initialize G4EM0 \n");  
  G4EM[0] = theG4ErrorTrajState->GetTheTransfMat(); // actually function name: GetTransfMatrix - maybe renamed in newer version of geant?
  for(int i=0;i<5;i++) {
    for(int j=0;j<5;j++) {
      G4EM[0][i][j] = 0.0;
      if(i==j) G4EM[0][i][j] = 1.0;
    }
  }

  for(int ipass=0;ipass<1;ipass++) {

    // NSTOP is 20 for 20 planes
  for(int ipl = 0;ipl<NSTOP-1;ipl++) {
   int itbstep = 0;


   //Forward tracking
   //    printf("Set the target number %d \n",ipl);
   //  G4Point3D surfPos(-474.0+10.0*(ipl+1),0.,0.);
   // G4Normal3D surfNorm(1.,0.,0.);
   /// G4ErrorTarget* theG4ErrorTarget = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
  

   //Backward tracking
   //    printf("Set the target number %d \n",ipl);
  G4Point3D surfPos(-274.0-10.0*(ipl+1),0.,0.);
  G4Normal3D surfNorm(-1.,0.,0.);
  G4ErrorTarget* theG4ErrorTarget = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
   // target set from Rob's geometry somehow, from plane to plane

  g4edata->SetTarget(theG4ErrorTarget);


  //  printf("beginning Plane %d \n",ipl);

  if( iProp == 0){
    // Propagate until G4ErrorTarget is found all in one go
     ierr = g4emgr1->Propagate( theG4ErrorTrajState, theTarget, theG4ErrorMode );
  } else if( iProp == 1){

    // Propagate until G4ErrorTarget is reached step by step
  
    g4emgr1->InitTrackPropagation();

    //G4Track* aTrack = G4EventManager::GetEventManager()->GetTrackingManager()->GetTrack();
    bool moreEvt = TRUE;
      while( moreEvt ){
       
      //      printf("The G4 Error Mode is %d\n",theG4ErrorMode);      
	//	      G4cout << "The track before One Step" << theG4ErrorTrajState->GetParameters() << G4endl; 
	      printf("On Plane %d \n",ipl);
      ierr = g4emgr1->PropagateOneStep( theG4ErrorTrajState, theG4ErrorMode );
      // G4cout << "The track after One Step" << theG4ErrorTrajState->GetParameters() << G4endl; 
      //      printf("Checking is target reached \n");
      //---- Check if target is reached

      // G4cout << "The Transport Matrix " << theG4ErrorTrajState->GetTheTransfMat() << G4endl; 

      //Add momentum and position at each step into the vector.
      G4Point3D tempPos = theG4ErrorTrajState->GetPosition();
      G4Normal3D tempMom = theG4ErrorTrajState->GetMomentum();

      xpos[0] = tempPos[0];
      xpos[1] = tempPos[1];
      xpos[2] = tempPos[2];
      dl = sqrt( (xpos[0]-xposlast[0])*(xpos[0]-xposlast[0]) +
                        (xpos[1]-xposlast[1])*(xpos[1]-xposlast[1]) +
			(xpos[2]-xposlast[2])*(xpos[2]-xposlast[2]) );

      printf("Load transport matrix 1 dl %f\n",dl);
       printf("itbstep %d \n",itbstep);
      if(dl>0.00001 && dl<100.0) {
      if(itbstep==0) {
	printf("0 case dl %e %d \n",dl,ipl);
        G4EM[ipl+1] = theG4ErrorTrajState->GetTheTransfMat();
      }
      else {
	printf("not 0 case %e %d \n",dl,ipl);
        G4EM[ipl+1] = theG4ErrorTrajState->GetTheTransfMat() * G4EM[ipl+1];
        // This is where I believe Rob multiplies transport matrices when stepping through, but not reaching the target in one go 
        // - not sure if I can just use the iProp = 0 propogate until target mode or not
      }
      }
      // G4cout << "The Transport Matrix " << theG4ErrorTrajState->GetTheTransfMat() << G4endl; 

      //       printf("Load transport matrix 2 \n");

      //       G4cout << "The track before check" << theG4ErrorTrajState->GetParameters() << G4endl; 
     
      

      xposlast[0] = xpos[0];
      xposlast[1] = xpos[1];
      xposlast[2] = xpos[2];
      //      printf("Set xposlast \n");
      itbstep++;
      if( g4emgr1->GetPropagator()->CheckIfLastStep( theG4ErrorTrajState->GetG4Track() )) {
	g4emgr1->GetPropagator()->InvokePostUserTrackingAction( theG4ErrorTrajState->GetG4Track() );  
	moreEvt = 0;
	//        G4cout << "The track at target" << theG4ErrorTrajState->GetParameters() << G4endl; 
	// G4cout << "The Transport Matrix " << theG4ErrorTrajState->GetTheTransfMat() << G4endl; 
	// G4cout << "The Global Transport Matrix " << G4EM[itbstep-1] << G4endl; 
        G4ErrorTrajErr errorMid = theG4ErrorTrajState->GetError();
	// G4cout << " Error: " << errorMid << G4endl; 


	//	G4cout << "STEP_BY_STEP propagation: Last Step " << G4endl;
      }
      }
      

  }
  //  G4cout << "The track before prediction" << theG4ErrorTrajState->GetParameters() << G4endl; 




      G4Point3D tempPos2 = theG4ErrorTrajState->GetPosition();
      G4Normal3D tempMom2 = theG4ErrorTrajState->GetMomentum();
      t5tp[0][ipl+1] = theG4TrackParam1->GetInvP();
      t5tp[1][ipl+1] = theG4TrackParam1->GetLambda();
      t5tp[2][ipl+1] = theG4TrackParam1->GetPhi();
      t5tp[3][ipl+1] = tempPos2[1];
      t5tp[4][ipl+1] = tempPos2[2];


#if WRITE // WRITE set as 0 when Rob sent it to me

      rvec[ipl][0] = tempPos2[0];
      rvec[ipl][1] = tempPos2[1];
      rvec[ipl][2] = tempPos2[2];
      pvec[ipl][0] = tempMom2[0];
      pvec[ipl][1] = tempMom2[1];
      pvec[ipl][2] = tempMom2[2];
      rvec[ipl][0] = tempPos2[0];
      rvec[ipl][1] = tempPos2[1];
      rvec[ipl][2] = tempPos2[2];
      pvec[ipl][0] = tempMom2[0];
      pvec[ipl][1] = tempMom2[1];
      pvec[ipl][2] = tempMom2[2];
      //      printf("Predicted Plane %d position %e %e %e momentum %e %e %e \n",    ipl,  rvec[ipl][0],
      //     rvec[ipl][1],
      //     rvec[ipl][2],
      //     pvec[ipl][0],
      //     pvec[ipl][1],
      //     pvec[ipl][2]);

      myfile <<  std::setw(6) << std::setprecision(3) << ipl << 
                       "  " << std::setw(17) << std::setprecision(9) <<
                       "  " << rvec[ipl][0] << " " << rvec[ipl][1]  << 
                       "  " << rvec[ipl][2] << " " << pvec[ipl][0]  << 
                       "  " << pvec[ipl][1] << " " << pvec[ipl][2]  << G4endl;
#endif  

    //  printf("Predicted Plane %d position %e %e \n", ipl,  
      //     t5tp[3][ipl],
      //     t5tp[4][ipl]);

      //update transport matrix

      //      printf("Update transport matrix for plane %d \n",ipl);
      //      G4cout << "The Transport Matrix for current plane" << G4EM[ipl] << G4endl; 
      //      G4cout << "The Transport Matrix for previous plane" << G4EM[ipl-1] << G4endl; 


      G4EM[ipl+1] = G4EM[ipl+1]*G4EM[ipl];
      //  printf("End of plane %d \n",ipl);
  }  //loop over planes
  //We need to create the measured track, predicted track and transport matrix

#if READ
  trycorr(iev,G4EM,t5tm,t5tp,dtvec);
#endif
  } //loop over passes

  // for(int i=0;i<5;i++) {
  //  printf(" parameter %d delta track vec %f \n",i,dtvec[i]);
  // }
  myfile.close();
  G4cout << " $$$ PROPAGATION ENDED " << G4endl;


}


//-------------------------------------------------------------
G4ErrorTarget* BuildTarget( G4int iTarget )
{

  G4ErrorTarget* target = 0;
  if( iTarget == 1 ) {
    G4Point3D surfPos(-2.66666667*cm,0.,0.);
    G4Normal3D surfNorm(1.,0.,0.);
    target = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
    printf("Going with target 1 \n");
  }else if( iTarget == 2 ) {
    G4Point3D surfPos(-47.4*cm,0.,0.);
    G4Normal3D surfNorm(-1.,0.,0.);
    target = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
    //    G4double radius = 222*cm;
    //target = new G4ErrorCylSurfaceTarget(radius);
  }else if( iTarget == 3 ) {
    target = new G4ErrorGeomVolumeTarget("MUON");
  }else if( iTarget == 4 ) {
    target = new G4ErrorTrackLengthTarget(223.*cm);
  }else {
    G4Exception("exG4eReco::BuildTarget. Target type has to be between 1 and 4");
  }
  return target;
}


//-------------------------------------------------------------
void Finalize()
{
  printf("Finalize \n");
  G4ErrorPropagatorManager::GetErrorPropagatorManager()->CloseGeometry();
  printf("Return from Finalize \n");

}


void trycorr(int iev, G4ErrorMatrix *G4EM, double trkparm[][20], double trkparp[][20], double *dtrack) {


//Pass track measured, track predicted, trp0i
//do matrix inversion

double res[5][NSTOP];
//double trkparp[5][NSTOP];
//double trkparm[5][NSTOP];
double tmp[5][5][NSTOP];
double cov[5][5][NSTOP];
double covtot[5][5];
double sol[5];
double vtmp1[5][NSTOP];
double vtmp2[5][NSTOP];
double vtmp3[5];
double vtmp4[5];
double vtmp5[5];
double covinv[5][5];
double dlocp[5][NSTOP];
double trp0i[5][5][NSTOP];
double trp0it[5][5][NSTOP];
double vv[5][5][NSTOP];
double G4EMMeVmm[50][5][5];
double G4temp[50][5][5];

 double scale34[5][5];
 printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n");

 // if(iev==0) {
   scale34[0][0] = 1.0e0;
   scale34[1][0] = 1.0e3;
   scale34[2][0] = 1.0e3;
   scale34[3][0] = 1.0e4;
   scale34[4][0] = 1.0e4;
   scale34[0][1] = 1.0e-3;
   scale34[1][1] = 1.0e0;
   scale34[2][1] = 1.0e0;
   scale34[3][1] = 1.0e1;
   scale34[4][1] = 1.0e1;
   scale34[0][2] = 1.0e-3;
   scale34[1][2] = 1.0e0;
   scale34[2][2] = 1.0e0;
   scale34[3][2] = 1.0e1;
   scale34[4][2] = 1.0e1;
   scale34[0][3] = 1.0e-4;
   scale34[1][3] = 1.0e-1;
   scale34[2][3] = 1.0e-1;
   scale34[3][3] = 1.0e0;
   scale34[4][3] = 1.0e0;
   scale34[0][4] = 1.0e-4;
   scale34[1][4] = 1.0e-1;
   scale34[2][4] = 1.0e-1;
   scale34[3][4] = 1.0e0;
   scale34[4][4] = 1.0e0;
   //} 

int nstop = NSTOP;
 
    for (int ipl=0;ipl<nstop;ipl++) {
     G4cout << "The Global Transport Matrix plane " << ipl << G4EM[ipl] << G4endl; 
    }



//Scale them

for(int itrp=0;itrp<nstop;itrp++) {
   for(int i=0;i<5;i++) {
     for(int j=0; j<5;j++) {
	   G4EMMeVmm[itrp][i][j] = G4EM[itrp][i][j]*scale34[i][j]; 
	   G4temp[itrp][i][j] = G4EM[itrp][i][j]*scale34[i][j]; 
    }
  }
 }

// for(int itrp=0;itrp<nstop;itrp++) {
//   for(int i=0;i<5;i++) {
//     for(int j=0; j<5;j++) {
//       G4EMMeVmm[itrp](i,j) = G4temp[itrp](i,j);
//      }
//   }
// }
	

 //Reorder
 // 1/p , y', z', y , z looks to be actually switched to 1/p, z', y', y, z (actually y' might be lambda, z' phi)
 // have to then switch matrix elements in certain ways
 for(int itrp=0;itrp<nstop;itrp++) {
   G4EMMeVmm[itrp][0][1] = G4temp[itrp][0][2];
   G4EMMeVmm[itrp][0][2] = G4temp[itrp][0][1];
   G4EMMeVmm[itrp][1][0] = G4temp[itrp][2][0];
   G4EMMeVmm[itrp][2][0] = G4temp[itrp][1][0];


   G4EMMeVmm[itrp][3][1] = G4temp[itrp][3][2];
   G4EMMeVmm[itrp][3][2] = G4temp[itrp][3][1];
   G4EMMeVmm[itrp][1][3] = G4temp[itrp][2][3];
   G4EMMeVmm[itrp][2][3] = G4temp[itrp][1][3];



   G4EMMeVmm[itrp][4][1] = G4temp[itrp][4][2];
   G4EMMeVmm[itrp][4][2] = G4temp[itrp][4][1];
   G4EMMeVmm[itrp][1][4] = G4temp[itrp][2][4];
   G4EMMeVmm[itrp][2][4] = G4temp[itrp][1][4];


   G4EMMeVmm[itrp][1][1] = G4temp[itrp][2][2];
   G4EMMeVmm[itrp][2][2] = G4temp[itrp][1][1];
   G4EMMeVmm[itrp][1][2] = G4temp[itrp][2][1];
   G4EMMeVmm[itrp][2][1] = G4temp[itrp][1][2];



 }


     for (int ipl=nstop-1;ipl>-1;ipl--) {
       printf("Plane %d \n",ipl);
       for(int i=0;i<5;i++) {
         printf(" %e %e %e %e %e \n",   G4EMMeVmm[ipl][i][0],G4EMMeVmm[ipl][i][1],G4EMMeVmm[ipl][i][2],G4EMMeVmm[ipl][i][3],G4EMMeVmm[ipl][i][4] );
       }
       printf("\n");
       printf("\n");
    }

//load up the errors
         for(int ipl=0;ipl<nstop;ipl++) {
            for(int j=0;j<5;j++) {
	      for(int k=0;k<5;k++) {
	        if( j==k && (k==3||k==4) ) {
		  vv[j][k][ipl] = 100.0; //0.01 cm resolution (is it cm?)
                }
                else
                  vv[j][k][ipl] = 0.0;
              }
            }
	  }


         for(int ipl=0;ipl<nstop;ipl++) {
            for(int j=0;j<5;j++) {
	      for(int k=0;k<5;k++) {
		trp0i[j][k][ipl] = G4EMMeVmm[ipl][j][k];
              }
            }
	  }

#if 0
         for (int ipl=0;ipl<nstop;ipl++) {
           printf("%d dy/ d invp %f \n",ipl,trp0i[3][0][ipl]);
         }
#endif
#if 0
         for (int ipl=0;ipl<nstop;ipl++) {
           printf("%d dz/ dz' %e \n",ipl,trp0i[4][2][ipl]);
         }
#endif
	 printf("Predict the delta tracks \n");

 
#if 1
         double pv3c = 3.0*GeV;
         double pv3p = (3.0-0.3)*GeV;
         printf("Correct momentum %f \n",pv3c);
         double dinvp = 1.0/pv3p - 1.0/pv3c; 
         double del5[5];
	 printf("Change in inverse p %f \n",dinvp);
#endif
 
#if 0
         double yawzc = 0.1/3.0;
	 double yawzp = (0.1-iev*0.01)/3.0;
         double del5[5];
         double dyz = yawzc-yawzp;
	 printf("Change in yawz %f \n",dyz);
#endif


#if 1
         for(int ipl=0;ipl<nstop;ipl++) {
            for(int j=0;j<5;j++) {
              del5[j] = 0.0;
	      for(int k=0;k<5;k++) {
                if(k==0) {
		  del5[j] += trp0i[j][k][ipl]*dinvp;
                }
              }
		//		trp0i[j][k][ipl] = G4EM[ipl][j][k];
              if(j==3) {
	      printf("plane %d coordinate %d predicted delta y %f \n",ipl,j,del5[j]);
              }
            }
	 }

#endif




#if 1
	 printf("Here are the delta-tracks \n");
   for(int k=0;k<nstop;k++) {
     for(int j=0;j<5;j++) {
      res[j][k] = trkparm[j][k]-trkparp[j][k];
     }
     printf("%d y:%f %f %f %e z: %f %f %f \n",k,trkparp[3][k],trkparm[3][k],res[3][k],
	    dinvp*trp0i[3][0][k],  trkparp[4][k],trkparm[4][k],res[4][k]);
   }
#endif

#if 0
	 printf("Here are the delta-tracks \n");
   for(int k=0;k<nstop;k++) {
     for(int j=0;j<5;j++) {
      res[j][k] = trkparm[j][k]-trkparp[j][k];
     }
     printf("%d y:%f %f %f  z: %f %f %f \n",k,trkparp[3][k],trkparm[3][k],res[3][k],
     trkparp[4][k],trkparm[4][k],res[4][k]);
   }
#endif

   printf("About to return \n");

//        DO IPL=0,NSTOP
//           DO J=1,5
//              DO K=1,5
//                 TRP0I(J,K,IPL) = TRP0IT(K,J,IPL)
//              ENDDO
//           ENDDO
//        ENDDO

//Calculate the transpose of the transport matrix

 for(int ipl=0;ipl<nstop;ipl++) {
            for(int j=0;j<5;j++) {
	      for(int k=0;k<5;k++) {
		trp0it[j][k][ipl] = trp0i[k][j][ipl];
              }
            }
	  }
//C
//C       Form the covariance matrix for the trial solution:
//C       ==================================================/
//        DO IPL=1,NSTOP
//           DO J=1,5
//              DO K=1,5
//                 TMP(J,K,IPL) = 0.0
//                 DO L =1,5
//                    TMP(J,K,IPL) = TMP(J,K,IPL) +
//     >                            VV(J,L,IPL)*TRP0I(L,K,IPL)
//                 ENDDO
//              ENDDO
//           ENDDO
//        ENDDO

 for(int ipl=0;ipl<nstop;ipl++) {
            for(int j=0;j<5;j++) {
	      for(int k=0;k<5;k++) {
                tmp[j][k][ipl] = 0.0;
                for(int l=0;l<5;l++) {
                  tmp[j][k][ipl] = tmp[j][k][ipl] +
		    vv[j][l][ipl]*trp0i[l][k][ipl];
		}
              }
            }
	  }




//        DO IPL=1,NSTOP
//           DO J=1,5
//              DO K=1,5
//                 COV(J,K,IPL) = 0.0
//                 DO L =1,5
//                    COV(J,K,IPL) = COV(J,K,IPL) +
//     >                        TRP0IT(J,L,IPL)*TMP(L,K,IPL)
//                 ENDDO
//              ENDDO
//           ENDDO
//        ENDDO

 for(int ipl=0;ipl<nstop;ipl++) {
            for(int j=0;j<5;j++) {
	      for(int k=0;k<5;k++) {
                cov[j][k][ipl] = 0.0;
                for(int l=0;l<5;l++) {
                  cov[j][k][ipl] = cov[j][k][ipl] +
		    trp0it[j][l][ipl]*tmp[l][k][ipl];
		}
           

              }
            }
	  }




//        DO J=1,5
//           DO K=1,5
//              COVTOT(J,K) = 0.0
//              DO IPL=1,NSTOP
//                 COVTOT(J,K) = COVTOT(J,K) + COV(J,K,IPL)
//              ENDDO
//           ENDDO
//        ENDDO

for(int j=0;j<5;j++) {
  for(int k=0;k<5;k++) {
    covtot[j][k] = 0.0;
    for(int ipl = 0; ipl<nstop; ipl++) {
      covtot[j][k] = covtot[j][k] + cov[j][k][ipl];
    }
  }
 }


//        DO J=1,5
//           DO K=1,5
//                 COVINV(J,K) = COVTOT(J,K)
//           ENDDO
//        ENDDO
//            for(int j=0;j<5;j++) {
//	      for(int k=0;k<5;k++) {
//                covinv[j][k] = covtot[j][k];/
//		}
//              }
    TMatrixD H_square(5,5);


    Double_t determ[1];
 #if 1
  G4ErrorMatrix covtotG4 = G4ErrorMatrix(5,5,0);
  for(int j=0;j<5;j++) {
    for(int k=0;k<5;k++) {
      H_square(j,k) = covtot[j][k];
    }
  }



  printf("Before inversion \n \n");
  for( int j=0;j<5;j++) {
    printf("%f %f %f %f %f \n",H_square(j,0),H_square(j,1),H_square(j,2),H_square(j,3),H_square(j,4));
  }

  H_square.InvertFast(determ);

  //  printf("Inversion return code %d \n",ierr);

  printf("After inversion \n \n");
  for( int j=0;j<5;j++) {
    printf("%f %f %f %f %f \n",H_square(j,0),H_square(j,1),H_square(j,2),H_square(j,3),H_square(j,4));
  }



//C
//C       Invert:
//C       =======
//        CALL DINV(5,COVTOT,5,RWORK,IFAIL,1,SOL)
//C
//C       Check the inversion:
//C       ====================
//        DO I=1,5
//           DO J=1,5
//             COVWRK(I,J) = 0.0
//                DO K=1,5
//                   COVWRK(I,J) = COVWRK(I,J) +
//     >               COVTOT(I,K)*COVINV(K,J)
//                ENDDO
//           ENDDO
//c           WRITE(6,23) I,(COVWRK(I,L),L=1,5)
//c23         FORMAT(1x,I4,5(1x,F10.3))
//        ENDDO
  for(int j=0;j<5;j++) {
    for(int k=0;k<5;k++) {
      covinv[j][k] = H_square(j,k);
    }
  }

  printf("Check the inversion \n");
  double unitmat[5][5];
  for(int j=0;j<5;j++) {
    for(int k=0;k<5;k++) {
      unitmat[j][k] = 0.0;
      for(int l=0;l<5;l++) {
        unitmat[j][k] = unitmat[j][k] + covinv[j][l]*covtot[l][k];
      }
    }
  }
  for( int j=0;j<5;j++) {
    printf("%f %f %f %f %f \n",unitmat[j][0],unitmat[j][1],unitmat[j][2],unitmat[j][3],unitmat[j][4]);
  }

  printf("This should be a unit matrix \n \n \n");



//C
//C       Now we can calculate the next guess:
//C       ====================================

//C
//C       Calculate the residuals at every step:
//C       ======================================
//        DO K=0,NSTOP
//           DO J=1,5
//              RES(J,K) = CSIZP(J,K)-CSIX(J,K)
//           ENDDO
//        ENDDO

   


//C
//C       Multiply the residuals by the measurement error matrix:
//C       =======================================================
//        DO IPL=1,NSTOP
//          DO J=1,5
//             VTMP1(J,IPL) = 0.0
//             DO K=1,5
//                VTMP1(J,IPL) = VTMP1(J,IPL) + VV(J,K,IPL)*RES(K,IPL)
//             ENDDO
//          ENDDO
//        ENDDO
	       for(int ipl=0;ipl<nstop;ipl++) {
		 for(int j=0;j<5;j++) {
		   vtmp1[j][ipl] = 0.0;
		   for(int k=0;k<5;k++) {
		      vtmp1[j][ipl] = vtmp1[j][ipl] +
			vv[j][k][ipl]*res[k][ipl];
		   }
		 }
	       }


//C
//C       Now multiply by the transpose of the transport matrix:
//C       ======================================================/
//        DO IPL=1,NSTOP
//          DO J=1,5
//             VTMP2(J,IPL) = 0.0
//             DO K=1,5
//                VTMP2(J,IPL) = VTMP2(J,IPL) +
//     >              TRP0IT(J,K,IPL)*VTMP1(K,IPL)
//             ENDDO
//          ENDDO
//        ENDDO


 for(int ipl=0;ipl<nstop;ipl++) {
   for(int j=0;j<5;j++) {
     vtmp2[j][ipl] = 0.0;
     for(int k=0;k<5;k++) {
       vtmp2[j][ipl] = vtmp2[j][ipl] + trp0it[j][k][ipl]*vtmp1[k][ipl];
     }
   }
 }

 //C
 //C       Sum over all planes:
 //C       ====================
 //       DO J=1,5
 //          VTMP3(J) = 0.0
 //          DO IPL=1,NSTOP
 //             VTMP3(J) = VTMP3(J) + VTMP2(J,IPL)
 //          ENDDO
 //       ENDDO
       
 for(int j=0;j<5;j++) {
   vtmp3[j] = 0.0;
   for(int ipl=0;ipl<nstop;ipl++) {
     vtmp3[j] = vtmp3[j] + vtmp2[j][ipl];
   }
 }


 
 //C
 //C       Calculate the solution:
 //C       =======================
 //       DO J=1,5
 //          SOL(J) = 0.0
 //          DO K=1,5
 //             SOL(J) = SOL(J) + COVTOT(J,K)*VTMP3(K)
 //          ENDDO
 //       ENDDO

 for(int j=0;j<5;j++) {
   sol[j] = 0.0;
   for(int k=0;k<5;k++) {
     sol[j] = sol[j] + covinv[j][k]*vtmp3[k];
   }
 }
 

 
 //c        DO I=1,5
 //c           WRITE(6,87) I,SOL(I)
 //c87      FORMAT(' Approximate: Parameter ',I6,' Delta ',F15.5)
 //c        ENDDO


 for(int i=0;i<5;i++) {
   printf("parameter %d delta %f \n",i,sol[i]);
   dtrack[i] = sol[i];
 }

 return;
 //C
 //C       Check this also:
 //C       ================
 //        DO I=1,5
 //          VTMP5(I) = 0.0
 //          DO K=1,5
 //             VTMP5(I) = VTMP5(I) + COVINV(I,K)*SOL(K)
 //          ENDDO
 //       ENDDO

 for(int i=0;i<5;i++) {
   vtmp5[i] = 0.0;
   for(int k=0;k<5;k++) {
     vtmp5[i] = vtmp5[i] + covinv[i][k]*sol[k];
   }
 }
 

 //C
 //C       Compare:
 //C       ========
 //C        PRINT *, ' Check solution'
 //        DO I=1,5
 //c           WRITE(6,45) I, VTMP3(I),VTMP5(I)
 // 45        FORMAT(' PARAMETER ',I5,' ANSWER ',F12.5,
 //    >              ' CHECK ',F12.5)
 //       ENDDO

 for(int i=0;i<5;i++) {
   printf("parameter %d answer %f check %f \n",i,vtmp3[i],vtmp5[i]);
 }

 // C
 //C     Calculate change at each plane:
 //C     ===============================
 //     DO I=0,NSTOP
 //       DO J=1,5
 //          DLOCP(J,I) = 0.0
 //          DO K=1,5
 //             DLOCP(J,I) = DLOCP(J,I) + TRP0I(J,K,I)*SOL(K)
 //           ENDDO
 //        ENDDO
 //     ENDDO

 for(int ipl=0;ipl<nstop;ipl++) {
   for(int j=0;j<5;j++) {
     dlocp[j][ipl] = 0.0;
     for(int k = 0;k<5;k++) {
       dlocp[j][ipl] = dlocp[j][ipl] + trp0i[j][k][ipl]*sol[k];
     }
   }
 }
#endif

 return;



 }





#if 0
  if(iev%4 ==0) {
    xv3.set( -474.0, 0, 0 );
    pv3.set( 3.0*GeV, 0.0, 0.0 );
    g4edata->SetTarget( theTarget1 );
    theG4ErrorMode = G4ErrorMode_PropForwards;  
  }
  else {
    xv3.set( -26.666667, 0, 0 );
    pv3.set( -2.9*GeV, 0.0, 0.0 );
    g4edata->SetTarget( theTarget2 );
    theG4ErrorMode = G4ErrorMode_PropBackwards;  
  }
#endif



  //JPAdd
  //There are a few ways to do this, but because of the file structure,
  //after the root file is written in the initialization, we'll pass it here, then
  //open the root file again for each event within this member function.
#if 0
  TTree *PosTree = (TTree*) rootFile->Get("PosTree");
  std::vector<TVector3> *PositionVect = new std::vector<TVector3>;
  TVector3 Position;
  TTree *MomTree = (TTree*) rootFile->Get("MomTree");
  std::vector<TVector3> *MomentumVect = new std::vector<TVector3>;
  TVector3 Momentum;
  TTree *ErrMatTree = (TTree*) rootFile->Get("ErrMatTree");
  std::vector< std::vector< std::vector<float> > > *ErrMatrixVect = new std::vector< std::vector< std::vector<float> > >;


  PosTree->Branch("PositionVect","vector<TVector3>",&PositionVect);
  MomTree->Branch("MomentumVect","vector<TVector3>",&MomentumVect);
  ErrMatTree->Branch("ErrMatrixVect","vector<vector<vector<float>>>",&ErrMatrixVect);
#endif
  ////////////////////////end JPAdd


#if 0


      Position.SetXYZ(tempPos.x(), tempPos.y(), tempPos.z());
      Momentum.SetXYZ(tempMom.x(), tempMom.y(), tempMom.z());

      PositionVect->push_back(Position);
      MomentumVect->push_back(Momentum);

      //Record the error matrix, but in vector form
      
      std::vector< std::vector< float > > ErrMatrix;
      
      
      ErrMatrix.resize(G4EM[ipl].num_row());
      for (int i=0; i!=ErrMatrix.size(); ++i) {
	ErrMatrix[i].resize(G4EM[ipl].num_col());
      }

      for (int i=0; i!=G4EM[ipl].num_row(); ++i) {      
	for (int j=0; j!=G4EM[ipl].num_col(); ++j) {      
	  ErrMatrix[i][j] = G4EM[ipl][i][j];
	}
      }
      ErrMatrixVect->push_back(ErrMatrix);
      ErrMatrix.clear();
#endif      
      /////////////////////////////////////////////////////endJPAdd


  //JPAdd
#if 0
  PosTree->Fill();
  MomTree->Fill();
  ErrMatTree->Fill();
#endif
  ////////////////////////////////////////////endJPAdd

#if 0

	if(iev==0) {
        t5orig0[0] = theG4TrackParam1->GetInvP();
        t5orig0[1] = pv3.y()/pv3.x();
        t5orig0[2] = pv3.z()/pv3.x();
        t5orig0[3] = xv3.y();
        t5orig0[4] = xv3.z();

        t5orig0[1] = theG4TrackParam1->GetLambda();
        t5orig0[2] = theG4TrackParam1->GetPhi();
        t5orig0[3] = theG4TrackParam1->GetYPerp();
        t5orig0[4] = theG4TrackParam1->GetZPerp();

        }
	if(iev==1) {
        t5orig1[0] = theG4TrackParam1->GetInvP();

        t5orig1[1] = pv3.y()/pv3.x();
        t5orig1[2] = pv3.z()/pv3.x();
        t5orig1[3] = xv3.y();
        t5orig1[4] = xv3.z();

        t5orig1[1] = theG4TrackParam1->GetLambda();
        t5orig1[2] = theG4TrackParam1->GetPhi();
        t5orig1[3] = theG4TrackParam1->GetYPerp();
        t5orig1[4] = theG4TrackParam1->GetZPerp();
        }



  G4cout << " Position: " << posEnd << G4endl
         << " Momentum: " << momEnd << G4endl
         << " Error: " << errorEnd << G4endl; 

  printf("Position %f %f %f \n",posEnd[0],posEnd[1],posEnd[2]);
	//Make vector out of final points and grab local coords
        G4ErrorFreeTrajParam* theG4TrackParam = new G4ErrorFreeTrajParam( posEnd,momEnd);
        printf("Inv P value %f \n",theG4TrackParam->GetInvP());
        printf("Lambda value %f \n",theG4TrackParam->GetLambda());
        printf("Phi value %f \n",theG4TrackParam->GetPhi());
        printf("YPerp value %f \n",theG4TrackParam->GetYPerp());
        printf("Zperp value %f \n",theG4TrackParam->GetZPerp());
        G4cout << "The final track" << theG4ErrorTrajState->GetParameters() << G4endl; 


	if(iev==0) {
	  printf("loading event 0 \n");
        t5final0[0] = theG4TrackParam1->GetInvP();
        t5final0[1] = theG4TrackParam1->GetLambda();
        t5final0[2] = theG4TrackParam1->GetPhi();
        t5final0[3] = theG4TrackParam1->GetYPerp();
        t5final0[4] = theG4TrackParam1->GetZPerp();

        t5final0[1] = momEnd[1]/momEnd[0];
        t5final0[2] = momEnd[2]/momEnd[0];
        t5final0[3] = posEnd[1];
        t5final0[4] = posEnd[2];

        }
	if(iev==1) {
	  //        G4cout << "The Global Transport Matrix " << G4EM[itbstep-1] << G4endl; 
	  //	printf("loading event 1 \n");
        t5final1[0] = theG4TrackParam1->GetInvP();
        t5final1[1] = theG4TrackParam1->GetLambda();
        t5final1[2] = theG4TrackParam1->GetPhi();
        t5final1[3] = theG4TrackParam1->GetYPerp();
        t5final1[4] = theG4TrackParam1->GetZPerp();

        t5final1[1] = momEnd[1]/momEnd[0];
        t5final1[2] = momEnd[2]/momEnd[0];
        t5final1[3] = posEnd[1];
        t5final1[4] = posEnd[2];

        }

        if(iev==1) {
	  for(int i=0;i<5;i++) {
	    dt5o[i] = t5orig1[i]-t5orig0[i];
          }
	  for(int i=0;i<5;i++) {
	    dt5f[i] = t5final1[i]-t5final0[i];
          }
          printf(" Plane 0 yprime %e yprime altered %e \n",t5orig0[1],t5orig1[1]);
          printf(" Plane 0 zprime %e zprime altered %e \n",t5orig0[2],t5orig1[2]);

          printf(" Plane N yprime %e yprime altered %e \n",t5final0[1],t5final1[1]);
          printf(" Plane N zprime %e zprime altered %e \n",t5final0[2],t5final1[2]);

          printf(" Plane 0 yperp %e yperp altered %e \n",t5orig0[3],t5orig1[3]);
          printf(" Plane 0 zperp %e zperp altered %e \n",t5orig0[4],t5orig1[4]);

          printf(" Plane N yperp %e yperp altered %e \n",t5final0[3],t5final1[3]);
          printf(" Plane N zperp %e zperp altered %e \n",t5final0[4],t5final1[4]);

     
          printf("\n \n");
          printf("P0 Change in inv mom %f \n",dt5o[0]);
          printf("P0 Change in y' %f \n",dt5o[1]);
          printf("P0 Change in z' %e \n",dt5o[2]);
          printf("P0 Change in y perp %f \n",dt5o[3]);
          printf("P0 Change in z perp %f \n",dt5o[4]);
          printf("\n" );

          printf("PN Change in inv mom %f \n",dt5f[0]);
          printf("PN Change in y' %f \n",dt5f[1]);
          printf("PN Change in z' %e \n",dt5f[2]);
          printf("PN Change in y perp %f \n",dt5f[3]);
          printf("PN Change in z perp %f \n",dt5f[4]);
          for(int i=0;i<5;i++) {
             dt5fch[i] = 0.0;
	    for(int j=0;j<5;j++) {
	    dt5fch[i] += G4EM[itbstep-1][i][j]*dt5o[j];
	    }
       
	  }
          printf("inv mom %f %f \n",dt5f[0],dt5fch[0]);
          printf("y prime %e %e \n",dt5f[1],dt5fch[1]);
          printf("z prime %e %e \n",dt5f[2],dt5fch[2]);
          printf("y %e %e \n",dt5f[3],dt5fch[3]);
          printf("z %e %e \n",dt5f[4],dt5fch[4]);

	}
#endif


#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>

#include <Eigen/Dense>

// Not all these needed probably.
#include "TVector3.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TInterpreter.h"
#include "TRandom3.h"
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TF2.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TAxis.h>

// #include "CLHEP/Matrix/Matrix.h" // Since I've switched to Eigen I shouldn't need CLHEP anymore. Leaving the library linking in the Cmake file regardless just in case.

#include "G4SteppingVerbose.hh"
#include "G4ErrorPropagator.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagatorManager.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorCylSurfaceTarget.hh"
#include "G4ErrorGeomVolumeTarget.hh"
#include "G4ErrorTrackLengthTarget.hh"
#include "G4ErrorFreeTrajState.hh"
#include "G4ErrorSurfaceTrajState.hh"
#include "G4ErrorMatrix.hh"

#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

// Visualization include
// #include "G4VisExecutive.hh"
// #include "G4UIExecutive.hh"

#include <RootInput.hh>
#include "BasicDetectorConstruction.hh"
#include "MyParallelWorld.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#define MATRIXDEBUG 0


static const int maxNumPlanes = 33;
static const int numPasses = 4; // This will create residuals on the detector planes based on numPasses-1 passes, and residuals on the starting point based on numPasses.
static const int numTrackParams = 5;

int numEventsSkipped = 0; // Number of events skipped due to bad momentum guesses.
std::vector<int> numEventsSkippedWithNumPlanes (maxNumPlanes, 0); // Not a method. Number of events skipped that hit a certain number of planes.
int numEventsSkippedBadStartingPosition = 0; // Number of events skipped due to bad starting position.

void Initialize();
G4ErrorTarget* BuildTarget( G4int iTarget );
void ProcessEvent( G4int iProp, int currentEventNum, int numPlanesHit, int arrayPosition );
void Finalize();
// void TrackCorrelation(G4ErrorMatrix *myTransferMatrices, G4ErrorMatrix* errorMatrices, double paramMeasured[][maxNumPlanes], double paramPredicted[][maxNumPlanes], int numPlanesHit, double deltaStartingTrack[5]);
void TrackCorrelation(std::vector<G4ErrorMatrix> myTransferMatrices, std::vector<G4ErrorMatrix> errorMatrices, std::vector <std::vector<G4double> > paramMeasured, std::vector <std::vector<G4double> > paramPredicted, int numPlanesHit, double deltaStartingTrack[5], double& eventChiSquared, int firstPlaneHit);

void ResidualPlots();

G4ErrorTarget* theTarget;
G4ErrorMode theG4ErrorMode;
G4ErrorPropagatorManager* g4emgr;
G4ErrorPropagatorData* g4edata; 


const G4ThreeVector vy(0,1,0); // Set the v and w axes. Here they are simply y and z.
const G4ThreeVector wz(0,0,1);

const double noHit = 900000.; // Unphysical double value to signify the lack of a hit within a tracker plane.

/////////////////////////////////////////////////////////////////////////////////////
// My global variables.
/////////////////////////////////////////////////////////////////////////////////////

// Build map of vectors of parameters which I will use later to compare to GEANE calculated values.
std::map<std::string, std::vector<double> > InputParameters;

// Blank position map to fill vector of maps.
std::map<std::string, double > MapOfDoubles;
std::vector<std::map<std::string, double > > PlanePositionMeasured; //This is a vector of maps that will be filled with the separate plane position measured maps.
std::vector<std::map<std::string, double > > PlaneParameterTruth; //This is a vector of maps that will be filled with the separate plane position truth maps.

/////////////////////////////////////////////////////////////////////////////////////
// Global histograms.

TH1F* TotalMomentumResiduals[maxNumPlanes];
TH1F* YMomentumResiduals[maxNumPlanes];
TH1F* ZMomentumResiduals[maxNumPlanes];
TH1F* YPositionResiduals[maxNumPlanes];
TH1F* ZPositionResiduals[maxNumPlanes];
TH2F* PositionMeasuredHist[maxNumPlanes];
TH2F* PositionTruthHist[maxNumPlanes];
TH1F* PositionMeasuredY[maxNumPlanes];
TH1F* MeasurementInaccuracyY[maxNumPlanes];
TH1F* MeasurementInaccuracyZ[maxNumPlanes];
TH1F* ChiSquaredHistogram;
  // TH1F* ChiSquaredProbHistogram;


/////////////////////////////////////////////////////////////////////////////////////

// std::vector<double> ChiSquaredValues; // .reserve() space for this
// std::vector<double> ChiSquaredProb;


// Blank residuals map to fill vector of maps.
std::map<std::string, double > TracebackResidualsMap;
std::vector<std::map<std::string, double > > PlaneTracebackResiduals; //This is a vector of maps that will be filled with the separate plane residual maps.


BasicDetectorConstruction* myPhysicalWorld = new BasicDetectorConstruction(); // To acquire geometry values.
MyParallelWorld* myGhostWorld = new MyParallelWorld("fillerString"); // It's possible making a parallel world object here could screw things up somehow later.

/////////////////////////////////////////////////////////////////////////////////////
// End global variables/
/////////////////////////////////////////////////////////////////////////////////////

int main(/*int argc,char** argv*/)
{
	RootInput *input = new RootInput(0); // Feed in the previously recorded truth values from the ROOT tree made using the basic.exe program.
	InputParameters = input->RootInput::LoopAndFill(InputParameters);

	std::cout << "Event parameters are: " << std::endl; 

	for (int i = 0; i < int(InputParameters["EventID"].size()); ++i)
	{
		std::cout << "Event ID: " << InputParameters["EventID"].at(i) << std::endl
	            << "GlobalTime: " << InputParameters["GlobalTime"].at(i) << std::endl
	            << "ProperTime: " << InputParameters["ProperTime"].at(i) << std::endl
	            << "XPosition: " << InputParameters["XPosition"].at(i) << std::endl
	            << "YPosition: " << InputParameters["YPosition"].at(i) << std::endl
	            << "ZPosition: " << InputParameters["ZPosition"].at(i) << std::endl
	            << "XMomentum: " << InputParameters["XMomentum"].at(i) << std::endl
	            << "YMomentum: " << InputParameters["YMomentum"].at(i) << std::endl
	            << "ZMomentum: " << InputParameters["ZMomentum"].at(i) << std::endl
              << "CopyNo: " << InputParameters["CopyNo"].at(i) << std::endl
              << "UPosition: " << InputParameters["UPosition"].at(i) << std::endl
              << "VPosition: " << InputParameters["VPosition"].at(i) << std::endl
	            << std::endl;
	}

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {
    PlanePositionMeasured.push_back(MapOfDoubles);
    PlaneParameterTruth.push_back(MapOfDoubles);

    PlaneTracebackResiduals.push_back(TracebackResidualsMap); 

  }

/////////////////////////////////////////////////////////////////////////////////////
  // Create array of histograms to be filled at each process event method.

  double zHistBound;
  double yHistBound;

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {

    if (planeNum==0) { zHistBound = 14; yHistBound = 3;}
    else{ zHistBound = 1; yHistBound = 1;}

  TotalMomentumResiduals[planeNum] = new TH1F(("Total Momentum Residuals Plane "+std::to_string(planeNum)).c_str(),"Total Momentum Residuals; Total Momentum Residual (MeV); Number of Events", 100, -300, 300); // Added a factor of 10 to all these numbers for the curved 0 plane.
  YMomentumResiduals[planeNum] = new TH1F(("Y Momentum Residuals Plane "+std::to_string(planeNum)).c_str(),"Y Momentum Residuals; Y Momentum Residual (MeV); Number of Events", 100, -50, 50);
  ZMomentumResiduals[planeNum] = new TH1F(("Z Momentum Residuals Plane "+std::to_string(planeNum)).c_str(),"Z Momentum Residuals; Z Momentum Residual (MeV); Number of Events", 100, -50, 50);
  YPositionResiduals[planeNum] = new TH1F(("Y Position Residuals Plane "+std::to_string(planeNum)).c_str(),"Y Position Residuals; Y Position Residual (mm); Number of Events",100,-yHistBound,yHistBound);
  ZPositionResiduals[planeNum] = new TH1F(("Z Position Residuals Plane "+std::to_string(planeNum)).c_str(),"Z Position Residuals; Z Position Residual (mm); Number of Events",100,-zHistBound,zHistBound);

  PositionMeasuredHist[planeNum] = new TH2F(("Position Measured Plane "+std::to_string(planeNum)).c_str(),"Position Measured; Z Position (mm); Y Position (mm); Number of Events",100,(-1.*myGhostWorld->GetTruthPlaneHalfZ()),(1.*myGhostWorld->GetTruthPlaneHalfZ()),100,(-1.*myGhostWorld->GetTruthPlaneHalfY()),(1.*myGhostWorld->GetTruthPlaneHalfY()));
  PositionTruthHist[planeNum] = new TH2F(("Position Truth Plane "+std::to_string(planeNum)).c_str(),"Position Truth; Z Position (mm); Y Position (mm); Number of Events",100,(-1.*myGhostWorld->GetTruthPlaneHalfZ()),(1.*myGhostWorld->GetTruthPlaneHalfZ()),100,(-1.*myGhostWorld->GetTruthPlaneHalfY()),(1.*myGhostWorld->GetTruthPlaneHalfY()));

  PositionMeasuredY[planeNum] = new TH1F(("Position Measured Y Plane "+std::to_string(planeNum)).c_str(),"Y Measurement; Y Measurement (mm); Number of Events",100,(-1.*myGhostWorld->GetTruthPlaneHalfY()),(1.*myGhostWorld->GetTruthPlaneHalfY()));

  MeasurementInaccuracyY[planeNum] = new TH1F(("Y Inaccuracy Plane "+std::to_string(planeNum)).c_str(),"Y Measurement Inaccuracy; Y Measurement Inaccuracy (mm); Number of Events",100,-.3,.3);
  MeasurementInaccuracyZ[planeNum] = new TH1F(("Z Inaccuracy Plane "+std::to_string(planeNum)).c_str(),"Z Measurement Inaccuracy; Z Measurement Inaccuracy (mm); Number of Events",100,-.3,.3);

  }

  ChiSquaredHistogram = new TH1F("Chi Squared Values for all Events", "Chi Squared Values; Chi Squared; Number of Events", 120, 0, 60);
  // TH1F* ChiSquaredProbHistogram = new TH1F("Chi Squared PDF", "Chi Squared PDF; Chi Squared; Number of Events", 100, 0, 50);
/////////////////////////////////////////////////////////////////////////////////////

  // Visualization disabled.
  // G4VisManager* visManager = new G4VisExecutive;
  // visManager->Initialize();

  // G4UIExecutive* ui = 0;
  // if ( argc == 1 ) {
  //   ui = new G4UIExecutive(argc, argv);
  // }
  //   G4UImanager* UImanager = G4UImanager::GetUIpointer();

/////////////////////////////////////////////////////////////////////////////////////

  Initialize();




  //----- Choose propagation mode // Leaving this here for now for future changes if necessary.
  // 0: propagate until target, all steps in one go
  // 1: propagate until target, returning control to the user at each step
  G4int iProp = 0; // Default to target.
  char* prop = getenv("G4ERROR_PROP");
  if( prop ) {
    if( G4String(prop) == G4String("UNTIL_TARGET") ){
      iProp = 0;
    } else if ( G4String(prop) == G4String("STEP_BY_STEP") ) {
      iProp = 1;
    } else {
      G4Exception("exG4eReco","Fatal error in Argument",
        FatalErrorInArgument,
        G4String("Variable G4ERROR_PROP = " + G4String(prop) + 
                 "   It must be: UNTIL_TARGET or STEP_BY_STEP").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",
      JustWarning,
      "Variable G4ERROR_PROP not defined, taking it = UNTIL_TARGET");
  } 


  // EventID array within the map goes like: 0, 0, 0, 1, 1, 2, 3, 3, 4, 4, 4, 4 ...
  // Need to get last value to see how many events I need to run, and then determine the amount of planes crossed for each event, and pass those numbers to ProcessEvent.
  int nEvents = InputParameters["EventID"].back(); // Access last element.
  // ChiSquaredValues.reserve(nEvents);

  int eventArraySize = InputParameters["EventID"].size();

  int numEventsSkippedMissingPlanes=0; // Number of events skipped that didn't hit enough planes.

  for( int currentEventNum = 0; currentEventNum <= nEvents; currentEventNum++ )
  {
    // Loop through number of events.
    int numPlanesHit = 0; // For each event start out with 0 planes hit.      
    int arrayPos = 0;  // Position within the InputParameters array.

      for (int arrayPosition = 0; arrayPosition < eventArraySize; ++arrayPosition)
      { 
        // Loop through EventID array.
        if (currentEventNum == InputParameters["EventID"].at(arrayPosition))
        {// If the EventID number matches the number in the array, increment the number of planes hit.
          numPlanesHit++;
        }
        
        if (currentEventNum < InputParameters["EventID"].at(arrayPosition))
        {
          // If I pass beyond the current event number in the array, break out of the for loop with the proper arrayPosition to pass.
          break;
        }

        arrayPos = arrayPosition; // This for the last element in the array, otherwise it goes one over.

        
        
      }

      // G4cout << "Num planes hit: " << numPlanesHit << G4endl;
      // G4cout << "Array postion: " << arrayPos << G4endl;

      if (numPlanesHit >= 6) // Particles always hit plane 0. 6 for the 0 plane hit, 1 full module hit, and one plane from another module.
      {
        ProcessEvent( iProp, currentEventNum, numPlanesHit, arrayPos ); // Process current event.
      }
      else{
        numEventsSkippedMissingPlanes++;
      }

  }



  Finalize();

  ResidualPlots();

  G4cout << "Number of events skipped due to not hitting enough planes: " << numEventsSkippedMissingPlanes << " out of: " << nEvents+1 << G4endl;
  G4cout << "Number of events skipped due to bad starting position: " << numEventsSkippedBadStartingPosition << " out of: " << nEvents+1-numEventsSkippedMissingPlanes << G4endl;
  G4cout << "Number of events skipped due to bad momentum guesses: " << numEventsSkipped << " out of: " << nEvents+1-numEventsSkippedMissingPlanes-numEventsSkippedBadStartingPosition << G4endl;

  for (int i = 0; i < maxNumPlanes; ++i)
  {
    G4cout << "Planes hit: " << i << " Events Skipped: " << numEventsSkippedWithNumPlanes.at(i) << G4endl; // 0 plane is excluded essentially.
  }

  // ui->SessionStart();
  // delete ui;

  // delete visManager;
  delete input;

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////


void Initialize() 
{
  G4VSteppingVerbose::SetInstance(new G4SteppingVerbose);

  // Initialize the GEANT4e manager 
  g4emgr = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  g4edata = G4ErrorPropagatorData::GetErrorPropagatorData();

  g4emgr->SetUserInitialization(myPhysicalWorld); 

  g4emgr->InitGeant4e();

  G4UImanager::GetUIpointer()->ApplyCommand("/control/verbose 1");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 100 mm");
  // G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/iverbose 4"); // GEANT4E specific verbose command.

  char* mode = getenv("G4ERROR_MODE");
  if( mode ) {
    if( G4String(mode) == G4String("FORWARDS") ) {
      theG4ErrorMode = G4ErrorMode_PropForwards;
    } else if( G4String(mode) == G4String("BACKWARDS") ) {
      theG4ErrorMode = G4ErrorMode_PropBackwards;
    } else {
      G4Exception("exG4eReco","Fatal error in Argument",
        FatalErrorInArgument,
        G4String("Variable G4ERROR_MODE = " + G4String(mode) + 
                 "   It must be:  FORWARDS or BACKWARDS").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",
      JustWarning,"Variable G4ERROR_MODE not defined, taking it = BACKWARDS");
  } 

 

}


void ProcessEvent( G4int iProp, int currentEventNum, int numPlanesHit, int arrayPosition)
{

  // Initialization of vectors, arrays, matrices, etc. 
  static TRandom3 r3;

//iplfix  // G4double trackParamMeasured[5][maxNumPlanes];  // What is given from input root file.
  // G4double trackParamPredicted[5][maxNumPlanes]; // What GEANT4E calculates as the average from propagation.

  std::vector <std::vector<G4double> > trackParamMeasured;  // What is given from input root file.
  std::vector <std::vector<G4double> > trackParamPredicted; // What GEANT4E calculates as the average from propagation.

  // Start with noHit arrays of maxNumPlanes size for track vectors, then later erase noHit entries before sending to Track Correlation method.
  trackParamMeasured.resize(numTrackParams, std::vector<G4double>(maxNumPlanes, noHit)); // Initialize 2D vector with 5 rows and 9 columns filled with noHit values.
  trackParamPredicted.resize(numTrackParams, std::vector<G4double>(maxNumPlanes, noHit)); // Initialize this only with the number of planes hit for the number of columns.

  // Fill the measured values array here, with noHit values for 1/p, py/p, and pz/p, and some small uncertainty in y z. Units are MeV mm.

  int increment = 0; // increment along copyNo array-piece
  int firstPlaneHit = 0; // That's not the 0 plane.

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum) // This planeNum stands for actual plane number.
  {

    trackParamMeasured[0][planeNum] = noHit;
    trackParamMeasured[1][planeNum] = noHit;
    trackParamMeasured[2][planeNum] = noHit;

    // G4cout << "Test increment: " << increment << " Copy no: " << InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1) << " First plane hit number: " << firstPlaneHit << G4endl;

      
    if (planeNum == InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1))
    {
        trackParamMeasured[3][planeNum] = InputParameters["UPosition"].at(arrayPosition-numPlanesHit+increment+1)+r3.Gaus(0,.1);//+0.2*(r3.Rndm()-0.5); // +- 100 RMS microns in the position measurements.
        trackParamMeasured[4][planeNum] = InputParameters["VPosition"].at(arrayPosition-numPlanesHit+increment+1)+r3.Gaus(0,.1);//+0.2*(r3.Rndm()-0.5); // 
        // Array position is the last position in the array for an events data. For the 0 plane, array position is numPlanesHit-1, so I have to add 1 to start at 0. Similarly for the rest of the planes.
      
        PlaneParameterTruth[planeNum]["XPos"] = (InputParameters["XPosition"].at(arrayPosition-numPlanesHit+increment+1)); // increment vs planeNum
        PlaneParameterTruth[planeNum]["YPos"] = (InputParameters["YPosition"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["ZPos"] = (InputParameters["ZPosition"].at(arrayPosition-numPlanesHit+increment+1));

        PlaneParameterTruth[planeNum]["UPos"] = (InputParameters["UPosition"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["VPos"] = (InputParameters["VPosition"].at(arrayPosition-numPlanesHit+increment+1));

        PlaneParameterTruth[planeNum]["XMom"] = (InputParameters["XMomentum"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["YMom"] = (InputParameters["YMomentum"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["ZMom"] = (InputParameters["ZMomentum"].at(arrayPosition-numPlanesHit+increment+1));

        PlaneParameterTruth[planeNum]["CopyNo"] = (InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1));

        if (firstPlaneHit == 0) // Make sure the previous hit plane was plane 0, and then don't change it after it's been changed once.
        { firstPlaneHit = int(InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1)); }

        if (increment < numPlanesHit-1) // -1 to make sure that it doesn't increment over the size of the vector for the above if statement check
        { increment++; } //only increase increment if a plane is hit
        
    }
    
    else
    {
        trackParamMeasured[3][planeNum] = noHit;
        trackParamMeasured[4][planeNum] = noHit;

        PlaneParameterTruth[planeNum]["XPos"] = (noHit);
        PlaneParameterTruth[planeNum]["YPos"] = (noHit);
        PlaneParameterTruth[planeNum]["ZPos"] = (noHit);

        PlaneParameterTruth[planeNum]["UPos"] = (noHit);
        PlaneParameterTruth[planeNum]["VPos"] = (noHit);

        PlaneParameterTruth[planeNum]["XMom"] = (noHit);
        PlaneParameterTruth[planeNum]["YMom"] = (noHit);
        PlaneParameterTruth[planeNum]["ZMom"] = (noHit);

        PlaneParameterTruth[planeNum]["CopyNo"] = (noHit);

    }

    PlanePositionMeasured[planeNum]["YPos"] = (trackParamMeasured[3][planeNum]);
    PlanePositionMeasured[planeNum]["ZPos"] = (trackParamMeasured[4][planeNum]);



  }


  // G4cout << " Filled into parameter truth arrays: " << G4endl;
  // for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  // {
  //   std::cout << "XPosition: " << PlaneParameterTruth[planeNum]["XPos"] << std::endl
  //             << "YPosition: " << PlaneParameterTruth[planeNum]["YPos"] << std::endl
  //             << "ZPosition: " << PlaneParameterTruth[planeNum]["ZPos"] << std::endl
  //             << "XMomentum: " << PlaneParameterTruth[planeNum]["XMom"] << std::endl
  //             << "YMomentum: " << PlaneParameterTruth[planeNum]["YMom"] << std::endl
  //             << "ZMomentum: " << PlaneParameterTruth[planeNum]["ZMom"] << std::endl
  //             << "Copy No: " << PlaneParameterTruth[planeNum]["CopyNo"] << std::endl
  //             << std::endl;
  // }

  // G4cout << " Filled into measured parameter arrays: " << G4endl;
  //     for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  //     {
  //           G4cout << " Measured Track Parameters Plane: " << planeNum << G4endl
  //           << " 1/P: " << trackParamMeasured[0][planeNum] << G4endl
  //           << " PV or Y momentum: " << trackParamMeasured[1][planeNum] << G4endl
  //           << " PW or Z momentum: " << trackParamMeasured[2][planeNum] << G4endl
  //           << " V or Y position: " << trackParamMeasured[3][planeNum] << G4endl
  //           << " W or Z position: " << trackParamMeasured[4][planeNum] << G4endl;
  //     }
  //     G4cout << G4endl << G4endl;





// .back() get the last element of the vector, which just so happens to be the current event number we are processing
  G4double totalMomentum = sqrt( (PlaneParameterTruth[firstPlaneHit]["XMom"] * PlaneParameterTruth[firstPlaneHit]["XMom"])
                                +(PlaneParameterTruth[firstPlaneHit]["YMom"] * PlaneParameterTruth[firstPlaneHit]["YMom"])
                                +(PlaneParameterTruth[firstPlaneHit]["ZMom"] * PlaneParameterTruth[firstPlaneHit]["ZMom"]) );
  
// Set the initial starting trajectory.
// These units should still be in MeV mm and not GeV cm. Only the error matrix should be in GeV cm.
// Take the starting plane values, and add changes to similulate a non-perfect starting trajectory from some initial weak trackfitting algorithm.

  G4double initialXPosChange = -0.01; // The starting plane is now right in front of the first plane hit and uses its recorded values as it's starting parameters.
  G4double initialYPosChange = 40.*(r3.Rndm()-0.5); 
  G4double initialZPosChange = 40.*(r3.Rndm()-0.5); 

  G4double startingXPos = PlaneParameterTruth[firstPlaneHit]["XPos"] + initialXPosChange;
  G4double startingYPos = PlaneParameterTruth[firstPlaneHit]["YPos"] + initialYPosChange;
  G4double startingZPos = PlaneParameterTruth[firstPlaneHit]["ZPos"] + initialZPosChange;


  G4double initialXMomChange = 200.*(r3.Rndm()-0.5);
  G4double initialYMomChange = 200.*(r3.Rndm()-0.5);
  G4double initialZMomChange = 200.*(r3.Rndm()-0.5);

  G4double startingXMom = 1.*PlaneParameterTruth[firstPlaneHit]["XMom"] + initialXMomChange;
  G4double startingYMom = 1.*PlaneParameterTruth[firstPlaneHit]["YMom"] + initialYMomChange;
  G4double startingZMom = 1.*PlaneParameterTruth[firstPlaneHit]["ZMom"] + initialZMomChange;


// Variables for track parameter residual determination.
double deltaStartingTrack[5];
G4double totalMomentumModified;
double totalMomentumResidual;
double YPositionResidual;
double ZPositionResidual;
double YMomentumResidual;
double ZMomentumResidual;

double eventChiSquared;


  totalMomentumModified = sqrt ((startingXMom)*(startingXMom) 
                              + (startingYMom)*(startingYMom)
                              + (startingZMom)*(startingZMom));

  // G4double initialTotalMomentumModified = totalMomentumModified;

        // G4cout << "Beginning starting XPos: " << startingXPos << G4endl
        //        << "Beginning starting YPos: " << startingYPos << G4endl
        //        << "Beginning starting ZPos: " << startingZPos << G4endl
        //        << "Beginning starting XMom: " << startingXMom << G4endl
        //        << "Beginning starting YMom: " << startingYMom << G4endl
        //        << "Beginning starting ZMom: " << startingZMom << G4endl;




printf("\n Event number %d \n",currentEventNum);

        // G4cout << "Plane nums hit: ";
        // for (int i = 0; i < maxNumPlanes; ++i)
        // {
        //   if (PlanePositionMeasured[i]["YPos"].back() != noHit)
        //   {
        //     G4cout << " " << i << " ";
        //   }
        // }
        // G4cout << G4endl;

for (int passNumber = 0; passNumber < numPasses; ++passNumber)
{

  std::vector<G4ErrorMatrix> transferMatrices; // set of transfer/transport matrices, one for each plane, describes transport from one plane to the next plane.
  std::vector<G4ErrorMatrix> transferMatricesFreeSystem;
  std::vector<G4ErrorMatrix> errorMatrices; // set of error matrices, one for each plane

  transferMatrices.clear(); // Clear vectors at the beginning of each pass so things can be push_backed correctly into the vectors.
  transferMatricesFreeSystem.clear();
  errorMatrices.clear();

  for (int i = 0; i < numTrackParams; ++i)
  {
    trackParamPredicted[i].clear();
    trackParamPredicted[i].resize(maxNumPlanes, noHit);
  }

  double reconstructedTotalMomentum;
  double reconstructedYMomentumChange;
  double reconstructedZMomentumChange;


  // if (startingYPos > .9*myGhostWorld->GetTruthPlaneHalfY()) {startingYPos = .9*myGhostWorld->GetTruthPlaneHalfY();} // These are to make sure the traced particle doesn't start outside of the detector or world volume. 
  // else if (startingYPos < -.9*myGhostWorld->GetTruthPlaneHalfY()) {startingYPos = -.9*myGhostWorld->GetTruthPlaneHalfY();} // .9 so that it's far enough away from the edge so it doesn't seg fault if no planes are hit during the traceback.
  // if (startingZPos > .9*myGhostWorld->GetTruthPlaneHalfZ()) {startingZPos = .9*myGhostWorld->GetTruthPlaneHalfZ();}
  // else if (startingZPos < -.9*myGhostWorld->GetTruthPlaneHalfZ()) {startingZPos = -.9*myGhostWorld->GetTruthPlaneHalfZ();}

  G4ThreeVector xv3( startingXPos, startingYPos, startingZPos );
  G4ThreeVector pv3( startingXMom, startingYMom, startingZMom );

  // std::cout << std::endl << "Pass number: " << passNumber+1 << std::endl;


     theG4ErrorMode = G4ErrorMode_PropForwards;  

  // Need to provide error matrix at the start of the trajectory state. Units for the error (and transport) matrices are in GeV cm.
  G4ErrorTrajErr error( 5, 0 ); 
  // error(4,4) = 0.0001; // Parentheses instead of brackets start at 1,1 instead of 0,0 for G4ErrorMatrix.
  // error(5,5) = 0.0001; 
  // I think I comment these out and start with a 0 error matrix, and add in the experimental part later. Either that or the other way around. Or do both, not sure, it's somewhat unclear but doesn't change much of anything.
  // This is NOT inverted first when plugging into the initialization of the trajectory state, these values should be cm^2 and not (cm^2)^-1.
  // So for .1mm resolution or .01cm I should have .0001 values for (4,4) and (5,5). 


/////////////////////////////////////////////////////////////////////////////////////
  // This for a free trajectory state with 1/p, lambda, phi, y perp, and z perp.

  G4ErrorFreeTrajState* myFreeTrajState 
    = new G4ErrorFreeTrajState( /*particleGeneratorAction->GetParticleName()*/"e+", xv3, pv3, error );

  std::vector<G4ErrorMatrix> SC2SDTransformationMatrix; // Matrices for transfmormation to and from free system and plane/momentum system.
  std::vector<G4ErrorMatrix> SC2SDTransformationMatrixInverse;

  SC2SDTransformationMatrix.clear();
  SC2SDTransformationMatrixInverse.clear();

  G4ErrorMatrix fillTransformationMatrix0;

  // This for a surface trajectory state with 1/p, v', w', v, and w. 
  // This creates the surface trajectory state from the free state.
  G4ErrorSurfaceTrajState* mySurfaceTrajState
    = new G4ErrorSurfaceTrajState(*myFreeTrajState, vy, wz, fillTransformationMatrix0); // Passing it the empty transformation matrix here fills it with that matrix for the starting plane point.
  
   SC2SDTransformationMatrix.push_back(fillTransformationMatrix0);

  int dog;
  SC2SDTransformationMatrixInverse.push_back(SC2SDTransformationMatrix[0].inverse(dog)); // In this case, the 0 here does not refer to the 0 plane, but the starting point plane. This is necessary for the transport matrix coordinate transformation later.


  g4emgr = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  
// Initial track values.
  G4Point3D posBeg = myFreeTrajState->GetPosition();
  G4Normal3D momBeg = myFreeTrajState->GetMomentum();
  // G4ErrorFreeTrajParam theG4TrackParam1 = myFreeTrajState->GetParameters();
        // printf("\n Event number %d \n",currentEventNum);
        // printf("Init Inv P value %f \n",theG4TrackParam1.GetInvP());
        // printf("Init Lambda value %f \n",theG4TrackParam1.GetLambda());
        // printf("Init Phi value %f \n",theG4TrackParam1.GetPhi());
        // printf("Init YPerp value %f \n",theG4TrackParam1.GetYPerp());
        // printf("Init Zperp value %f \n",theG4TrackParam1.GetZPerp());
      // trackParamPredicted[0][0] = theG4TrackParam1.GetInvP();
      // trackParamPredicted[1][0] = theG4TrackParam1.GetLambda();
      // trackParamPredicted[2][0] = theG4TrackParam1.GetPhi();
      // trackParamPredicted[3][0] = posBeg[1]; // These are different than YPerp and ZPerp. 
      // trackParamPredicted[4][0] = posBeg[2];
      // trackParamPredicted[3][0] = theG4TrackParam1.GetYPerp(); // Commented out for being the wrong system to work in.
      // trackParamPredicted[4][0] = theG4TrackParam1.GetZPerp();

/////////////////////////////////////////////////////////////////////////////////////

  // If I want to create a surface trajectory from scratch.
  // G4ErrorSurfaceTrajState* mySurfaceTrajState 
    // = new G4ErrorSurfaceTrajState( "mu-", xv3, pv3, vy, wz, error );


// Initial track values for the first stage GEANE tracing.
  G4Point3D posBegVW = mySurfaceTrajState->GetPosition();
  G4Normal3D momBegVW = mySurfaceTrajState->GetMomentum();
  G4ErrorSurfaceTrajParam theG4TrackParamVW = mySurfaceTrajState->GetParameters();
        // printf("Init Inv P value %f \n",theG4TrackParamVW.GetInvP());
        // printf("Init PV value %f \n",theG4TrackParamVW.GetPV());
        // printf("Init PW value %f \n",theG4TrackParamVW.GetPW());
        // printf("Init V value %f \n",theG4TrackParamVW.GetV());
        // printf("Init W value %f \n",theG4TrackParamVW.GetW());

/////////////////////////////////////////////////////////////////////////////////////

  G4ErrorMatrix zeroMatrix = myFreeTrajState->GetTransfMat(); // Transport matrix for the free state when no steps have been made is a zero matrix.
  // G4cout << "Initial zero matrix: " << zeroMatrix << G4endl;

  // printf("initialize matrices \n");  // Fill with 0 matrices to start with. 
  transferMatrices.push_back(zeroMatrix); // Have to fill with an object to start otherwise it seg faults. Surface transfer matrices can also just be filled with the 0 matrix.
  transferMatricesFreeSystem.push_back(zeroMatrix);
  // Free state system has this particular method while the surface system does not. On plane 0 it is a 0 matrix which will be excluded in the GEANE formulation later on.
  // (There is no transfer matrix to plane 0 because there is no plane before it.)
  errorMatrices.push_back(mySurfaceTrajState->GetError()); // Starts off with 0 matrix as well.

  iProp = 1; // Propagation mode, 0 for until target, 1 for step by step. 
  // Specifically set the propagation mode here to step by step for some reason - can't remember. Errors if I don't? Leaving for now.

  for (int planeNum = 1; planeNum < maxNumPlanes; ++planeNum) // This planeNum stands for plane number hit, noHit entries from above are simply ignored/erased.
      {
        // G4cout << G4endl << "Tracking to plane: " << planeNum << G4endl;
/////////////////////////////////////////////////////////////////////////////////////
/* // This one for skipping nonhit planes completely.
        if (PlaneParameterTruth[planeNum]["CopyNo"].back() == noHit)
        { 
          // G4cout << " A detector plane was not hit, skip planeNum: " << planeNum << G4endl;

          transferMatrices.push_back(zeroMatrix); // Fill matrix arrays with zero arrays for no-hit planes. 
          transferMatricesFreeSystem.push_back(zeroMatrix);
          errorMatrices.push_back(zeroMatrix); // Same for error matrices.

          continue; 
        } // Skip the rest of the steps within the loop if a plane is missed. 
*/

        // This one for skipping over non-hit planes in front of the first plane hit. (Since I start the tracking from just before the first plane hit.)
        if (planeNum < firstPlaneHit) 
        {
          // G4cout << "The first plane hit was: " << firstPlaneHit << " Skipping planenum " << planeNum << G4endl;
          transferMatrices.push_back(zeroMatrix); // Fill matrix arrays with zero arrays for skipped planes.
          transferMatricesFreeSystem.push_back(zeroMatrix);
          errorMatrices.push_back(zeroMatrix); // Same for error matrices.

          continue; 
        }

/////////////////////////////////////////////////////////////////////////////////////

          int stepNumBetweenPlanes = 0;

          // Re-Build and set my own target here without the build target method for plane to plane propagation.
          G4double truthPlaneTargetX=0.; // 50 mm is where I set things up in the geometry from - only hardcoded value that I should probably change.
          G4double truthPlaneTargetY=0.;
          G4double truthPlaneTargetZ=0.;

          // Have to do some arithmatic to determine how to change truthPlaneTargetX with a single for loop for 32 planes.

          truthPlaneTargetX += int((planeNum-1)/4.) * myPhysicalWorld->GetmoduleSeparation();

          bool layerSep = FALSE;
          if (planeNum % 2 == 0)
          {
            truthPlaneTargetX += myPhysicalWorld->GetlayerSeparation();
            layerSep = TRUE;
          }

          double intpart;
          double fractpartone = std::modf(planeNum/4. , &intpart);
          // double intparttwo;
          double fractparttwo = std::modf((planeNum+1)/4. , &intpart);
          bool UVSep = FALSE;

          if (fractpartone == 0 || fractparttwo == 0)
          {
            truthPlaneTargetX += myPhysicalWorld->GetUVSeparation();
            UVSep = TRUE;
          }

          // G4cout << "Plane num is: " << planeNum << " module separation is: " << int((planeNum-1)/4.) << " layer separation is: " << layerSep << " UV separation is: " << UVSep << G4endl;

          
          // truthPlaneTargetX = (planeNum-1) * myPhysicalWorld->GetmoduleSeparation() + (myPhysicalWorld->GetlayerSeparation() + myPhysicalWorld->GetUVSeparation())/2.;
          
          // G4cout << "Starting position of free state: " << myFreeTrajState->GetPosition() << G4endl;
          // G4cout << "Truth plane target X is: " << truthPlaneTargetX << G4endl;

          if (myFreeTrajState->GetPosition().x() >= truthPlaneTargetX)
          { 
            numEventsSkippedBadStartingPosition++;
            return; 
          } // If the Geant SD recorded a too large of a starting point (323.1 mm instead of 323 mm), skip the event.
          


          G4Normal3D surfNorm(-1.,0.,0.);

/////////////////////////////////////////////////////////////////////////////////////
          
          G4Point3D surfPos(truthPlaneTargetX,truthPlaneTargetY,truthPlaneTargetZ);

          theTarget = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
          g4edata->SetTarget( theTarget );
/////////////////////////////////////////////////////////////////////////////////////
          if( iProp == 0){ // This seems to fail as only the transfer matrix for the last step gets saved, and if it takes multiple steps to get to the target, then the full transfer matrix from the last plane to the current one does not get properly saved.
                                // See lines 566 - 573 in G4ErrorFreeTrajState.cc. Looks like the proper code is commented out for some reason - my guess is no one has really used the transport matrix for anything before.
            // Propagate until G4ErrorTarget is found all in one go
              
             g4emgr->Propagate( myFreeTrajState, theTarget, theG4ErrorMode ); // If I come back to use this, I will need to intregrate the surfacetrajstate stuff in properly and adapt to the changes I've made.
             transferMatrices[planeNum] = myFreeTrajState->GetTransfMat(); // Same for here.
          } 
/////////////////////////////////////////////////////////////////////////////////////
          else if( iProp == 1){
          // Propagate until G4ErrorTarget is reached step by step


          g4emgr->InitTrackPropagation();

          bool moreEvt = TRUE;
          while( moreEvt ){
            
            // G4cout << " Particle momentum is: " << myFreeTrajState->GetMomentum() << " with magnitude: " << myFreeTrajState->GetMomentum().mag() << G4endl;
            
            g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode ); // If I start on a plane and propagate to the same plane, the step length propagated is 0 as it should be.
            // g4emgr->PropagateOneStep( mySurfaceTrajState, theG4ErrorMode ); // Just use free state and transform later, works better.
            
            // Use step by step propagation to properly save the transfer matrices.

            // G4cout << "Propagating to Plane number: " << planeNum << " stepNumBetweenPlanes: " << stepNumBetweenPlanes << " Step length: " << myFreeTrajState->GetG4Track()->GetStepLength() << G4endl;
            G4ErrorMatrix fillTransformationMatrix;
            mySurfaceTrajState
              = new G4ErrorSurfaceTrajState(*myFreeTrajState, vy, wz, fillTransformationMatrix); // Create a new surface trajectory state each step of the way from the propagating free trajectory state. The passed empty matrix gets filled with the transformation matrix.

            SC2SDTransformationMatrix.push_back(fillTransformationMatrix);

            int cat;
            SC2SDTransformationMatrixInverse.push_back((SC2SDTransformationMatrix.back()).inverse(cat));

            // G4cout << " Transformation matrix " << SC2SDTransformationMatrix.back() << G4endl
            // << " Transformation matrix inverse " << SC2SDTransformationMatrixInverse.back() << G4endl;

            // printf("Check SC2SD inversion worked properly \n");
            // std::cout << SC2SDTransformationMatrix.back()*SC2SDTransformationMatrixInverse.back();

            // Transfer matrices uses push_back since they are only filled for hit planes.
                if(stepNumBetweenPlanes==0) {
                  // From matrix multiplication, I believe it should be: R10 = A1 * T10 * A0^-1, where R and T are the transport matrices in separate coordinate systems from plane 0 to 1, and A is the transformation between the two. 
                  transferMatricesFreeSystem.push_back(myFreeTrajState->GetTransfMat()); // Gets the transport or transfer matrix for the last step. If there is only a single step between planes this is fine.
                  transferMatrices.push_back(SC2SDTransformationMatrix.back()*(myFreeTrajState->GetTransfMat())*SC2SDTransformationMatrixInverse.at(SC2SDTransformationMatrixInverse.size()-1));
                }
                else {
                  G4ErrorMatrix fillTransportMatrix;
                  fillTransportMatrix = transferMatricesFreeSystem.back(); // Get the last transport matrix of previous steps.
                  transferMatricesFreeSystem.pop_back(); // Delete the matrix within the vector, and then fill with the accumulated version.
                  transferMatricesFreeSystem.push_back(myFreeTrajState->GetTransfMat() * fillTransportMatrix); // Accumulates transport matrices from step to step for steps from last plane to the next plane.
                  
                  transferMatrices.pop_back();
                  transferMatrices.push_back(SC2SDTransformationMatrix.back()*(transferMatricesFreeSystem.back())*SC2SDTransformationMatrixInverse.at(SC2SDTransformationMatrixInverse.size()-1)); // Have to be careful with plane to plane transport vs step by step transport. I think this is the right way to do it.
                    // Here I gather up the transport matrices for the free system for all steps between planes, then I multiply on the right by the previous plane transformation matrix, and on the left by the next step transportation matrix. It keeps looping until the next step transformation matrix is the same as the next plane transformation matrix.
                  
                  // Can replace these pop_back and push_back methods with .at() = if I want to now.
                }
              


            stepNumBetweenPlanes++; // Increment the step number between planes. (The 1st step is the same as stepNumBetweenPlanes==0.)


            //---- Check if target is reached
            if( g4emgr->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() )) 
            {
              g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );  
              moreEvt = 0;
              // G4cout << "STEP_BY_STEP propagation: Last Step " << G4endl;
            }

            delete mySurfaceTrajState; // Move through the steps up to the target, deleting and remaking the surfacetrajstate every step of the way.
          }


          }// end of else if
/////////////////////////////////////////////////////////////////////////////////////
            // After getting to the target, once again create the surface trajectory state on the target.
            G4ErrorMatrix unneededMatrix;
            mySurfaceTrajState
              = new G4ErrorSurfaceTrajState(*myFreeTrajState, vy, wz, unneededMatrix); // The idea is to propagate in steps with the free state, then at the end turn that into a surface state and grab associated parameters and the error.

            G4Point3D tempPos = mySurfaceTrajState->GetPosition();

            trackParamPredicted[0].at(planeNum) = mySurfaceTrajState->GetParameters().GetInvP();
            trackParamPredicted[1].at(planeNum) = (mySurfaceTrajState->GetParameters().GetPV())*trackParamPredicted[0].at(planeNum); // These are momentum projections in V and W directions, and not V' and W'. 
            trackParamPredicted[2].at(planeNum) = (mySurfaceTrajState->GetParameters().GetPW())*trackParamPredicted[0].at(planeNum); // This is a different parameterization than III in the manual, but the transport matrix takes this into account.
            trackParamPredicted[3].at(planeNum) = mySurfaceTrajState->GetParameters().GetV(); 
            trackParamPredicted[4].at(planeNum) = mySurfaceTrajState->GetParameters().GetW();

            // G4cout << " Constructed/Predicted Track Parameters Plane: " << planeNum << G4endl
            // << " 1/P: " << trackParamPredicted[0].at(planeNum) << G4endl
            // << " PV or Y momentum: " << trackParamPredicted[1].at(planeNum) << G4endl
            // << " PW or Z momentum: " << trackParamPredicted[2].at(planeNum) << G4endl
            // << " V or Y position: " << trackParamPredicted[3].at(planeNum) << G4endl
            // << " W or Z position: " << trackParamPredicted[4].at(planeNum) << G4endl << G4endl;
 

            G4ErrorTrajErr errorEnd;
            errorEnd = mySurfaceTrajState->GetError();
            errorMatrices.push_back(errorEnd);

      }  // End planeNum loop



/////////////////////////////////////////////////////////////////////////////////////



      // G4cout << " Transport matrices size: " << transferMatrices.size() << G4endl
      //        << " Transport matrices free system size: " << transferMatricesFreeSystem.size() << G4endl
      //        << " Error matrices size: " << errorMatrices.size() << G4endl << G4endl;

      // G4cout << G4endl << " Matrices and vectors before track correlation: " << G4endl;
      // for (int ipl = 0; ipl < numPlanesHit; ++ipl)
      // {
      //   G4cout << " Transport matrix plane hit " << ipl << G4endl
      //          << transferMatrices[ipl] << G4endl
      //          << " Error matrix plane hit " << ipl << G4endl
      //          << errorMatrices[ipl] << G4endl << G4endl;

      //       G4cout << " Predicted Track Parameters Plane Hit: " << ipl << G4endl
      //       << " 1/P: " << trackParamPredicted[0][ipl] << G4endl
      //       << " PV or Y momentum: " << trackParamPredicted[1][ipl] << G4endl
      //       << " PW or Z momentum: " << trackParamPredicted[2][ipl] << G4endl
      //       << " V or Y position: " << trackParamPredicted[3][ipl] << G4endl
      //       << " W or Z position: " << trackParamPredicted[4][ipl] << G4endl << G4endl;
      // }

      // G4cout << "Track param predicted size: " << trackParamPredicted[3].size() << G4endl;
      // G4cout << "Track param measured size: " << trackParamMeasured[3].size() << G4endl;



      TrackCorrelation(transferMatrices, errorMatrices, trackParamMeasured, trackParamPredicted, numPlanesHit, deltaStartingTrack, eventChiSquared, firstPlaneHit);
      delete mySurfaceTrajState;
      delete myFreeTrajState; // Might not need to have these here as I think they might die anyways once the loop cycles.


      // Pull out the new filled delta starting track to get a new starting track.

      reconstructedTotalMomentum = 1./(deltaStartingTrack[0]+(1./totalMomentumModified));
      // totalMomentumResidual = totalMomentum - reconstructedTotalMomentum; // For forward tracking from the 0 plane.

      reconstructedYMomentumChange = deltaStartingTrack[1]*reconstructedTotalMomentum;
      reconstructedZMomentumChange = deltaStartingTrack[2]*reconstructedTotalMomentum;

      // These are the new starting positions and momenta.
      startingXPos = startingXPos + 0.; // There are no returned variables that allow for changing the starting x position. Address this later when necessary.
      startingYPos = startingYPos + deltaStartingTrack[3];
      startingZPos = startingZPos + deltaStartingTrack[4];

      startingYMom = startingYMom + reconstructedYMomentumChange;
      startingZMom = startingZMom + reconstructedZMomentumChange;

      startingXMom = sqrt((reconstructedTotalMomentum*reconstructedTotalMomentum)
                        - (startingYMom*startingYMom)
                        - (startingZMom*startingZMom));

/////////////////////////////////////////////////////////////////////////////////////

        // G4cout << "Total momentum modified before change: " << totalMomentumModified << G4endl;
        // G4cout << "1 over tmm before change: " << 1./totalMomentumModified << " compared to delta 0: " << deltaStartingTrack[0] << G4endl;

        int numPlanesHitForThisEvent=0;
        if(std::abs(deltaStartingTrack[0]) > std::abs(1./totalMomentumModified))
        {
          G4cout << "Event skipped because a terrible fix for the momentum was made for some reason. " << G4endl;

          G4cout << "Total momentum modified before change: " << totalMomentumModified << G4endl;
          G4cout << "1 over tmm before change: " << 1./totalMomentumModified << " compared to delta 0: " << deltaStartingTrack[0] << G4endl;
          G4cout << "Plane nums hit: ";
          for (int i = 0; i < maxNumPlanes; ++i)
          {
            if (PlanePositionMeasured[i]["YPos"] != noHit)
            {
              G4cout << " " << i << " ";
              numPlanesHitForThisEvent++;
            }
          }

          G4cout << G4endl << G4endl;

          numEventsSkippedWithNumPlanes.at(numPlanesHitForThisEvent-1)++; // Subtract 1, so if only plane 0 was hit, then 0 planes were hit, etc.
          numEventsSkipped++;
          return;
        }

        // Need to refill this variable here for when the for loop gets back to the reconstructedTotalMomentum step. This just incorporates deltaStartingTrack[0].
        totalMomentumModified = sqrt ((startingXMom)*(startingXMom) 
                               + (startingYMom)*(startingYMom)
                               + (startingZMom)*(startingZMom));

        // G4cout << G4endl << " Delta starting track passed out: " << G4endl
        // << " 1/P: " << deltaStartingTrack[0] << G4endl
        // << " Py/P: " << deltaStartingTrack[1] << G4endl 
        // << " Pz/P: " << deltaStartingTrack[2] << G4endl
        // << " Y: " << deltaStartingTrack[3] << G4endl
        // << " Z: " << deltaStartingTrack[4] << G4endl << G4endl;

        // G4cout << "Total momentum modified with delta change: " << totalMomentumModified << G4endl;
        // G4cout << "Reconstructed total momentum: " << reconstructedTotalMomentum << G4endl;
        // G4cout << "Rtm squared: " << reconstructedTotalMomentum*reconstructedTotalMomentum << " compared to Y and Z mom squared: " << (startingYMom*startingYMom) << " " << (startingZMom*startingZMom) << G4endl;

        // G4cout << " New starting XPos: " << startingXPos << G4endl
        //        << " New starting YPos: " << startingYPos << G4endl
        //        << " New starting ZPos: " << startingZPos << G4endl
        //        << " New starting XMom: " << startingXMom << G4endl
        //        << " New starting YMom: " << startingYMom << G4endl
        //        << " New starting ZMom: " << startingZMom << G4endl;

/////////////////////////////////////////////////////////////////////////////////////        


/* This for the previous only forward style tracking with the 0 plane included in the tracking. Will have to be fixed if I go back to it.
        YPositionResidual = startingYPos - PlaneParameterTruth[0]["YPos"].back(); // These 2 for truth - best fit line residual.
        ZPositionResidual = startingZPos - PlaneParameterTruth[0]["ZPos"].back();

        YMomentumResidual = (startingYMom) - PlaneParameterTruth[0]["YMom"].back(); // Residuals formed from truth - best fit.
        ZMomentumResidual = (startingZMom) - PlaneParameterTruth[0]["ZMom"].back();
*/

 // G4cout << G4endl << "End of pass: " << passNumber+1  << G4endl << G4endl;
} // end of numPasses for loop

      // ChiSquaredValues.push_back(eventChiSquared); // Add the calculated chi squared to the vector after the number of passes is done.
      // ChiSquaredProb.push_back(TMath::Prob(eventChiSquared, 10));

    ChiSquaredHistogram->Fill(eventChiSquared);


      // G4cout << "Chi^2 value: " << eventChiSquared << " Chi^2 prob: " << TMath::Prob(eventChiSquared, 10) << G4endl;

/////////////////////////////////////////////////////////////////////////////////////
  // Need to do a second stage GEANE backward tracing here from the first tracker plane hit to plane 0 and form the residuals from that. Actually think it needs to be outside the numPasses loop in order to re-adjust the GEANE manager stuff.

  G4ThreeVector xv3Backwards( startingXPos, startingYPos, startingZPos );
  G4ThreeVector pv3Backwards( -startingXMom, -startingYMom, -startingZMom ); // Reverse the momenta for backwards tracing.

  G4ErrorTrajErr errorStorage( 5, 0 ); // Error matrix for backward tracing to storage region.

  G4ErrorFreeTrajState* backwardsFreeTrajState 
    = new G4ErrorFreeTrajState("e+", xv3Backwards, pv3Backwards, errorStorage );

  double angle = 25.; // Use the same hardcoded value for the 0 plane for now.
  G4double storagePlaneTargetX = (-460.*std::sin(angle*pi/180))*cm;
  G4double storagePlaneTargetZ = (460.-460.*std::cos(angle*pi/180))*cm;

  G4Normal3D storageNorm(1.,0.,0.);
            // surfNorm.setX(std::cos(angle*pi/180));
            // surfNorm.setZ(-std::sin(angle*pi/180));
  G4Point3D storagePos(storagePlaneTargetX,0.,storagePlaneTargetZ);
  
  theTarget = new G4ErrorPlaneSurfaceTarget(storageNorm, storagePos );
  g4edata->SetTarget( theTarget );
  theG4ErrorMode = G4ErrorMode_PropBackwards;

  g4emgr->Propagate( backwardsFreeTrajState, theTarget, theG4ErrorMode ); // This doesn't save the transport matrices properly, but that's not a problem here since I don't care about them here.

  G4ErrorMatrix matrixFill;
  G4ErrorSurfaceTrajState* backwardsSurfaceTrajState
            = new G4ErrorSurfaceTrajState(*backwardsFreeTrajState, vy, wz, matrixFill);

  trackParamPredicted[0][0] = backwardsSurfaceTrajState->GetParameters().GetInvP();
  trackParamPredicted[1][0] = (backwardsSurfaceTrajState->GetParameters().GetPV())*trackParamPredicted[0][0]; // Dividing by total momenta here to keep things consistent even though it's not necessary.
  trackParamPredicted[2][0] = (backwardsSurfaceTrajState->GetParameters().GetPW())*trackParamPredicted[0][0]; 
  trackParamPredicted[3][0] = backwardsSurfaceTrajState->GetParameters().GetV(); 
  trackParamPredicted[4][0] = backwardsSurfaceTrajState->GetParameters().GetW();

  totalMomentumResidual = sqrt ((PlaneParameterTruth[0]["XMom"])*(PlaneParameterTruth[0]["XMom"]) 
                             + (PlaneParameterTruth[0]["YMom"])*(PlaneParameterTruth[0]["YMom"])
                             + (PlaneParameterTruth[0]["ZMom"])*(PlaneParameterTruth[0]["ZMom"])) - 1./trackParamPredicted[0][0];


  YMomentumResidual = PlaneParameterTruth[0]["YMom"] - -1.*trackParamPredicted[1][0]*1./trackParamPredicted[0][0]; // -1 in front since I think it's momentum might be the negative of what I want.
  ZMomentumResidual = PlaneParameterTruth[0]["ZMom"] - -1.*trackParamPredicted[2][0]*1./trackParamPredicted[0][0];

  YPositionResidual = PlaneParameterTruth[0]["YPos"] - trackParamPredicted[3][0];
  ZPositionResidual = PlaneParameterTruth[0]["ZPos"] - trackParamPredicted[4][0];


/////////////////////////////////////////////////////////////////////////////////////

        // G4cout << G4endl << " Starting truth parameters: " << G4endl
        // << " 1/P: " << (1/totalMomentum) << G4endl
        // << " PV or Y momentum: " <<  (PlaneParamterTruth[0]["YMom"]) << G4endl
        // << " PW or Z momentum: " <<  (PlaneParamterTruth[0]["ZMom"])<< G4endl
        // << " V or Y position: " << (PlaneParamterTruth[0]["YPos"]) << G4endl
        // << " W or Z position: " << (PlaneParamterTruth[0]["ZPos"]) << G4endl << G4endl;


        PlaneTracebackResiduals[0]["P"] = (totalMomentumResidual);
        PlaneTracebackResiduals[0]["Py"] = (YMomentumResidual);
        PlaneTracebackResiduals[0]["Pz"] = (ZMomentumResidual);
        PlaneTracebackResiduals[0]["Y"] = (YPositionResidual);
        PlaneTracebackResiduals[0]["Z"] = (ZPositionResidual);

/////////////////////////////////////////////////////////////////////////////////////

  for (int planeNum = 1; planeNum < maxNumPlanes; ++planeNum) // Planes > 0 are treated separately from plane 0. I only want to pass residuals for planes that were hit.
  {

    if (PlaneParameterTruth[planeNum]["CopyNo"] == noHit)
    { continue; }

    G4double planeTotalMomentum = sqrt ((PlaneParameterTruth[planeNum]["XMom"])*(PlaneParameterTruth[planeNum]["XMom"]) 
                                      + (PlaneParameterTruth[planeNum]["YMom"])*(PlaneParameterTruth[planeNum]["YMom"])
                                      + (PlaneParameterTruth[planeNum]["ZMom"])*(PlaneParameterTruth[planeNum]["ZMom"]));  

    PlaneTracebackResiduals[planeNum]["P"] = (planeTotalMomentum - 1./trackParamPredicted[0].at(planeNum));
    PlaneTracebackResiduals[planeNum]["Py"] = ((PlaneParameterTruth[planeNum]["YMom"]) - trackParamPredicted[1].at(planeNum)*1./trackParamPredicted[0].at(planeNum));
    PlaneTracebackResiduals[planeNum]["Pz"] = ((PlaneParameterTruth[planeNum]["ZMom"]) - trackParamPredicted[2].at(planeNum)*1./trackParamPredicted[0].at(planeNum));
    // PlaneTracebackResiduals[planeNum]["Y"] = (PlaneParameterTruth[planeNum]["YPos"] - trackParamPredicted[3].at(planeNum);// Truth - best fit line residual.
    // PlaneTracebackResiduals[planeNum]["Z"] = (PlaneParameterTruth[planeNum]["ZPos"] - trackParamPredicted[4].at(planeNum);
    PlaneTracebackResiduals[planeNum]["Y"] = (PlanePositionMeasured[planeNum]["YPos"] - trackParamPredicted[3].at(planeNum)); // Measured - best fit line residual.
    PlaneTracebackResiduals[planeNum]["Z"] = (PlanePositionMeasured[planeNum]["ZPos"] - trackParamPredicted[4].at(planeNum));
    // Using .back() here because of how the maps are filled from event to event. back() gets the last element which we are currently processing.

  }
  // With how I've changed things, this upper for loop and lower for loop can probably be combined, and the intermediate vector-map done away with, but I'm leaving it in for now for future functionality.


  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {

    TotalMomentumResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["P"],1);
    YMomentumResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Py"],1);
    ZMomentumResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Pz"],1);
    YPositionResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Y"],1);
    ZPositionResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Z"],1);

    PositionMeasuredHist[planeNum]->Fill(PlanePositionMeasured[planeNum]["ZPos"],PlanePositionMeasured[planeNum]["YPos"]);
    PositionTruthHist[planeNum]->Fill(PlaneParameterTruth[planeNum]["ZPos"],PlaneParameterTruth[planeNum]["YPos"]);

    PositionMeasuredY[planeNum]->Fill(PlanePositionMeasured[planeNum]["YPos"]);

    MeasurementInaccuracyY[planeNum]->Fill(PlanePositionMeasured[planeNum]["YPos"]-PlaneParameterTruth[planeNum]["YPos"]);
    MeasurementInaccuracyZ[planeNum]->Fill(PlanePositionMeasured[planeNum]["ZPos"]-PlaneParameterTruth[planeNum]["ZPos"]);

  }







        // printf("\n End of Event \n ");
  return;
        
  }



//-------------------------------------------------------------
G4ErrorTarget* BuildTarget( G4int iTarget )
{
  // Not using this for now but leaving it in just in case. Might be useful for more complex target building in the future.
  G4ErrorTarget* target = 0;
  if( iTarget == 1 ) {
    G4Point3D surfPos(5.*cm,0.,0.); 
    G4Normal3D surfNorm(1.,0.,0.);
    target = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
  }else if( iTarget == 2 ) {
    G4double radius = 222*cm;
    target = new G4ErrorCylSurfaceTarget(radius);
  }else if( iTarget == 3 ) {
    target = new G4ErrorGeomVolumeTarget("DetectorPlane");
  }else if( iTarget == 4 ) {
    target = new G4ErrorTrackLengthTarget(223.*cm);
  }else {
    G4Exception("exG4eReco","Fatal error in Argument",
                FatalErrorInArgument,"Target type has to be between 1 and 4");
  }
  return target;
}


//-------------------------------------------------------------
void Finalize()
{
  G4ErrorPropagatorManager::GetErrorPropagatorManager()->CloseGeometry();
}

// I think for the final program where it is unclear which measured parameters will go with which predicted parameters, I would pass in
// multiple different sets of paramMeasured to then calculate chi^2s and find the smallest one.
// Here I am just passing in what I know to be the smallest chi^2 parameters.
// void TrackCorrelation(G4ErrorMatrix* myTransferMatrices, G4ErrorMatrix* errorMatrices, double paramMeasured[][maxNumPlanes], double paramPredicted[][maxNumPlanes], int numPlanesHit, double deltaStartingTrack[5]){
void TrackCorrelation(std::vector<G4ErrorMatrix> myTransferMatrices, std::vector<G4ErrorMatrix> errorMatrices, std::vector <std::vector<G4double> > paramMeasured, std::vector <std::vector<G4double> > paramPredicted, int numPlanesHit, double deltaStartingTrack[5], double& eventChiSquared, int firstPlaneHit){
//Pass in array of transformed transfer/transport matrices, error matrices, track measured, track predicted, and the number of planes hit for each correlation method called.

// All ipl values in this method monotonically increase over hit planes.
  // Ex: Hit planes are 0, 1, 2, 5, 6, 8
  // Ipl goes as:       0, 1, 2, 3, 4, 5

std::vector<Eigen::VectorXd> residuals;

std::vector<double> chiSquaredSinglePlane;
chiSquaredSinglePlane.resize(maxNumPlanes, 0.); // Not sure whether to initially fill with 0 or some high number to simulate a noHit entry. I'm thinking 0.
double chiSquaredTotal = 0.;

std::vector<Eigen::MatrixXd> cov;
Eigen::MatrixXd covarianceTotal(5,5);

Eigen::VectorXd deltaPsiNought;

std::vector<Eigen::VectorXd> trialTrajSinglePlane;
Eigen::VectorXd trialTrajectoryAllPlanes;

std::vector<Eigen::MatrixXd> transportMatrixBegToEnd;
std::vector<Eigen::MatrixXd> sigmaErrorMatrices;

std::vector<Eigen::MatrixXd> myTransferMatricesGeVcm;
std::vector<Eigen::MatrixXd> errorMatricesGeVcm;

/*
// Error matrix is in units of GeV cm, originally wanted to send that to MeV mm.
 Eigen::MatrixXd transfMatScaling(5,5);

   transfMatScaling(0,0) = 1.0e0;
   transfMatScaling(1,0) = 1.0e3;
   transfMatScaling(2,0) = 1.0e3;
   transfMatScaling(3,0) = 1.0e4;
   transfMatScaling(4,0) = 1.0e4;
   transfMatScaling(0,1) = 1.0e-3;
   transfMatScaling(1,1) = 1.0e0; 
   transfMatScaling(2,1) = 1.0e0;
   transfMatScaling(3,1) = 1.0e1;
   transfMatScaling(4,1) = 1.0e1;
   transfMatScaling(0,2) = 1.0e-3;
   transfMatScaling(1,2) = 1.0e0;
   transfMatScaling(2,2) = 1.0e0;
   transfMatScaling(3,2) = 1.0e1;
   transfMatScaling(4,2) = 1.0e1;
   transfMatScaling(0,3) = 1.0e-4;
   transfMatScaling(1,3) = 1.0e-1;
   transfMatScaling(2,3) = 1.0e-1;
   transfMatScaling(3,3) = 1.0e0;
   transfMatScaling(4,3) = 1.0e0;
   transfMatScaling(0,4) = 1.0e-4;
   transfMatScaling(1,4) = 1.0e-1;
   transfMatScaling(2,4) = 1.0e-1;
   transfMatScaling(3,4) = 1.0e0;
   transfMatScaling(4,4) = 1.0e0;
  
    // Debug printing
    // for (int ipl=0;ipl<numPlanesHit;ipl++) {
    //  G4cout << "The Global Transport Matrix plane (before scaling) " << ipl << myTransferMatrices[ipl] << G4endl; 
    // }

// Error matrix is in units of GeV cm, want to send that to MeV mm.
 Eigen::MatrixXd errMatScaling(5,5);

   errMatScaling(0,0) = 1.0e-6;
   errMatScaling(0,1) = 1.0e-3;
   errMatScaling(1,0) = 1.0e-3;
   errMatScaling(0,2) = 1.0e-3;
   errMatScaling(2,0) = 1.0e-3;
   errMatScaling(0,3) = 1.0e-2;
   errMatScaling(3,0) = 1.0e-2;
   errMatScaling(0,4) = 1.0e-2;
   errMatScaling(4,0) = 1.0e-2;
   errMatScaling(1,1) = 1.0e0;
   errMatScaling(1,2) = 1.0e0;
   errMatScaling(2,1) = 1.0e0;
   errMatScaling(1,3) = 1.0e1;
   errMatScaling(3,1) = 1.0e1;
   errMatScaling(1,4) = 1.0e1;
   errMatScaling(4,1) = 1.0e1;
   errMatScaling(2,2) = 1.0e0;
   errMatScaling(2,3) = 1.0e1;
   errMatScaling(3,2) = 1.0e1;
   errMatScaling(2,4) = 1.0e1;
   errMatScaling(4,2) = 1.0e1;
   errMatScaling(3,3) = 1.0e2;
   errMatScaling(3,4) = 1.0e2;
   errMatScaling(4,3) = 1.0e2;
   errMatScaling(4,4) = 1.0e2;
*/

// GeV cm are the better units for the matrix multiplication to work out.

for(int ipl=0;ipl<maxNumPlanes;ipl++) {
   for(int i=0;i<5;i++) {
     for(int j=0; j<5;j++) {
          myTransferMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
          myTransferMatricesGeVcm[ipl](i,j) = myTransferMatrices[ipl][i][j]*1.;//transfMatScaling(i,j); GEVCM (Don't scale them.)
          errorMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
          errorMatricesGeVcm[ipl](i,j) = errorMatrices[ipl][i][j]*1.;//errMatScaling(i,j); GEVCM
     }
   } 
 } 

/////////////////////////////////////////////////////////////////////////////////////
#if MATRIXDEBUG
     // Debug printing
     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("Transport matrices in GeV cm %d \n",ipl);
       std::cout << myTransferMatricesGeVcm[ipl];
       printf("\n \n");
     }
// #if MATRIXDEBUG
     // Debug printing
     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("Propagated error matrices in GeV cm %d \n",ipl);
       std::cout << errorMatricesGeVcm[ipl];
       printf("\n \n");
     }
#endif
///////////////////////////////////////////////////////////////////////////////////// 

      for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
      { 
            // errorMatricesMeVmm[ipl](3,3) += 0.01; // Add the experimental part. Might be adding it twice however depending on if I give it the initial error or not.
            // errorMatricesMeVmm[ipl](4,4) += 0.01;

            errorMatricesGeVcm[ipl](3,3) += 0.0001; // GEVCM
            errorMatricesGeVcm[ipl](4,4) += 0.0001; // GEVCM

            // For 2x2 case:
            // Fill sigma with zeroes, grab 2x2 part of errorMatricesGeVcm and invert it, then plug back into sigma.
            // The sigma error matrices needs to be filled with 0s in the rows and columns where there are no measurements.
            // In the case of a single U or single V measurement, I simply invert that single number in the error matrix, with everything else zeros.
            // If a plane is not hit, it is filled with zeros. This matrix then goes on to multiply the other matrix pieces (possibly non-zero) in the tracking code, serving to remove any effects from non-hit planes that might otherwise introduce errors into the tracking.

            sigmaErrorMatrices.push_back(Eigen::MatrixXd::Zero(5,5));

            if (paramMeasured[3][ipl] == noHit && paramMeasured[4][ipl] == noHit)
            {
              continue; // Leave the sigma error matrix as zeros.
            }

            else if (paramMeasured[3][ipl] == noHit && paramMeasured[4][ipl] != noHit)
            {
              sigmaErrorMatrices[ipl](4,4) = 1./errorMatricesGeVcm[ipl](4,4); // Only fill the corresponding measured error matrix element.
            }

            else if (paramMeasured[4][ipl] == noHit && paramMeasured[3][ipl] != noHit)
            {
              sigmaErrorMatrices[ipl](3,3) = 1./errorMatricesGeVcm[ipl](3,3);
            }

            else
            {
              sigmaErrorMatrices[ipl].bottomRightCorner<2,2>() = (errorMatricesGeVcm[ipl].bottomRightCorner<2,2>()).inverse(); // This selects the bottom right 2x2 block.
            }

      }

/////////////////////////////////////////////////////////////////////////////////////
#if MATRIXDEBUG
     // Debug printing
     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("Full inverted error matrices %d \n",ipl);
       std::cout << sigmaErrorMatrices[ipl];
       printf("\n \n");
     }
#endif
/////////////////////////////////////////////////////////////////////////////////////

      // Multiply transport matrices together. Be careful with Transport matrix = 0 matrix for the starting plane.

      for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
      { 
        transportMatrixBegToEnd.push_back(Eigen::MatrixXd::Zero(5,5)); // Initialize matrices.
      }

      // transportMatrixBegToEnd[1] = myTransferMatricesGeVcm[1]; // Fill the plane 1 matrix, ignore the plane 0 matrix since it's just 0's and won't be included in the sum later on. (Same for the first hit plane.)
      transportMatrixBegToEnd[firstPlaneHit] = myTransferMatricesGeVcm[firstPlaneHit]; // Fill the plane 1 matrix, ignore the plane 0 matrix since it's just 0's and won't be included in the sum later on. (Same for the first hit plane.)


      for (int ipl = (firstPlaneHit+1); ipl < maxNumPlanes; ++ipl)
      {
        transportMatrixBegToEnd[ipl] = myTransferMatricesGeVcm[ipl]*transportMatrixBegToEnd[ipl-1]; // Multiply previous full transport matrices on the left by the next single gap transport matrix.
      }

/////////////////////////////////////////////////////////////////////////////////////
#if MATRIXDEBUG
     // Debug printing
     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("transport matrix beg to end Plane %d \n",ipl);
       std::cout << transportMatrixBegToEnd[ipl];
       printf("\n \n");
     }
#endif
/////////////////////////////////////////////////////////////////////////////////////

//       Form the covariance matrix for the trial solution:
//       ==================================================

      for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
      {
        cov.push_back(Eigen::MatrixXd::Zero(5,5));
        cov[ipl] = (transportMatrixBegToEnd[ipl].transpose())*sigmaErrorMatrices[ipl]*transportMatrixBegToEnd[ipl];
        // cov is the equation within the sum in equation 27 in the geane manual paper.
      }

#if MATRIXDEBUG
     // Debug printing
     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("T transpose sigmainverse T Plane %d \n",ipl);
       std::cout << cov[ipl];
       printf("\n \n");
     }
#endif

      covarianceTotal = Eigen::MatrixXd::Zero(5,5);
      for(int ipl = 1; ipl<maxNumPlanes; ipl++) { // Sum starts at 1 for the sum from plane to plane.
               covarianceTotal = covarianceTotal + cov[ipl];
                // This whole thing will then need to be inverted or used in an Eigen solve method.
      }

/////////////////////////////////////////////////////////////////////////////////////

    //       Check the inversion:   
    //       ====================

      // Note: Since I switched to Eigen checking the inversion is less important or at the very least needs to be done differently than it was before since I don't even invert a matrix, and I just use a solve method.
      // Can revisit later if needed.

/////////////////////////////////////////////////////////////////////////////////////
    //
    //       Now we can calculate the next guess: - In all of this, the beginning plane (highest number plane) should have zeros for matrix (and maybe parameter) values, so should not contribute.
    //       ====================================
    //
    //       Calculate the residuals at every step: - Far right side of equation 23, 24, or 26 in the geane manual paper.
    //       ======================================
    //             
            for(int ipl=0;ipl<maxNumPlanes;ipl++) {
              residuals.push_back(Eigen::VectorXd::Zero(5));
               for(int j=0;j<5;j++) {
                  // j the number of parameters
                  residuals[ipl](j) = paramMeasured[j][ipl]-paramPredicted[j][ipl]; // There will be noHit (900000) entries that find their way into these residuals, but they will drop out when multiplied by the sigma error matrices.
                }
            }

            // GEVCM Convert residuals from MeV mm to GeV cm
            for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
            {
              residuals[ipl](0) = residuals[ipl](0) * 1.0e3; // It's possible raising this could reduce the accuracy of the tracking - I might want to change it later. Shouldn't affect results besides the computer matrix math calculations since the momenta is unmeasured, and the noHit contributions drop out when contracted with the 0 part of the sigma error matrix.
              residuals[ipl](1) = residuals[ipl](1) * 1.0e0;
              residuals[ipl](2) = residuals[ipl](2) * 1.0e0;
              residuals[ipl](3) = residuals[ipl](3) * 1.0e-1;
              residuals[ipl](4) = residuals[ipl](4) * 1.0e-1;
            } 


#if MATRIXDEBUG
             // Debug printing
             printf("\n");
             for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
               printf("params measured Plane %d \n",ipl);
               std::cout << paramMeasured[0][ipl] << G4endl
                         << paramMeasured[1][ipl] << G4endl
                         << paramMeasured[2][ipl] << G4endl
                         << paramMeasured[3][ipl] << G4endl
                         << paramMeasured[4][ipl] << G4endl;
               printf("\n \n");
             } 

             printf("\n");
             for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
               printf("params predicted Plane %d \n",ipl);
               std::cout << paramPredicted[0][ipl] << G4endl
                         << paramPredicted[1][ipl] << G4endl
                         << paramPredicted[2][ipl] << G4endl
                         << paramPredicted[3][ipl] << G4endl
                         << paramPredicted[4][ipl] << G4endl;
               printf("\n \n");
             } 

             printf("\n");
             for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
               printf("residuals Plane %d \n",ipl);
               std::cout << residuals[ipl];
               printf("\n \n");
             } 
#endif

    //       Now multiply out the right side of equation 26 (within the sum).
    //       ======================================================/

             for(int ipl=0;ipl<maxNumPlanes;ipl++) {
               trialTrajSinglePlane.push_back(Eigen::VectorXd::Zero(5));
               trialTrajSinglePlane[ipl] = (transportMatrixBegToEnd[ipl].transpose())*sigmaErrorMatrices[ipl]*residuals[ipl];
             }

     //       Sum over all planes: --- This is the right side of equation 26 or 32 from the geane manual paper.
     //       ====================

            trialTrajectoryAllPlanes = Eigen::VectorXd::Zero(5);
             for(int ipl = 1; ipl<maxNumPlanes; ipl++) { // Sum starts at 1.
               trialTrajectoryAllPlanes = trialTrajectoryAllPlanes + trialTrajSinglePlane[ipl];
             }

     //       Calculate the solution:   this is deltaPsiNought from equation 26 or 32 of the geane manual paper
     //       =======================

             // deltaPsiNought = Eigen::VectorXd::Zero(5);
             deltaPsiNought = covarianceTotal.colPivHouseholderQr().solve(trialTrajectoryAllPlanes); // Better way to go about solving for deltaPsiNought without having to invert the covariance matrix.
             // Solves the equation Ax = b for x, where A = covarianceTotal, b = trialTrajectoryAllPlanes, and x = deltaPsiNought.

#if MATRIXDEBUG
             // Debug printing
             printf("starting parameter deltas \n");
             for(int i=0;i<5;i++) {
               printf("parameter %d \n",i);
               std::cout << deltaPsiNought[i] << std::endl;
             }
#endif
             // The values filled into delta starting track will be available for the operations after the TrackCorrelation method is called up above.
             // for (int i = 0; i < 5; ++i){ 
               // deltaStartingTrack[i] = deltaPsiNought(i);
             // }

             // GEVCM
              deltaStartingTrack[0] = deltaPsiNought(0) * 1.0e-3;
              deltaStartingTrack[1] = deltaPsiNought(1) * 1.0e0;
              deltaStartingTrack[2] = deltaPsiNought(2) * 1.0e0;
              deltaStartingTrack[3] = deltaPsiNought(3) * 1.0e1;
              deltaStartingTrack[4] = deltaPsiNought(4) * 1.0e1;


            // Form the chi^2 for a single plane then add them all up to get a total chi^2 for the potential track - equation 23 in the GEANE manual.
             for(int ipl=0;ipl<maxNumPlanes;ipl++) {
               chiSquaredSinglePlane[ipl] = (residuals[ipl].transpose())*sigmaErrorMatrices[ipl]*residuals[ipl];
             }

             for(int ipl = 1; ipl<maxNumPlanes; ipl++) { // Sum starts at 1.
               chiSquaredTotal = chiSquaredTotal + chiSquaredSinglePlane[ipl];
             }

             eventChiSquared = chiSquaredTotal;

             // G4cout << "Chi^2 for the track is: " << chiSquaredTotal << G4endl;

 return;

}

void ResidualPlots(){

  TCanvas* myCanvas = new TCanvas("myCanvas","Canvas",200,10,1200,800);

  double YRMSarray[maxNumPlanes];
  double ZRMSarray[maxNumPlanes];

  double pointNo[maxNumPlanes];// = {0,1,2,3,4,5,6,7,8};

  for (int i = 0; i < maxNumPlanes; ++i)
  {
    pointNo[i] = i;
  }

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {



  YRMSarray[planeNum] = YPositionResiduals[planeNum]->GetRMS();
  ZRMSarray[planeNum] = ZPositionResiduals[planeNum]->GetRMS();

  MeasurementInaccuracyY[planeNum]->Fit("gaus");
  MeasurementInaccuracyZ[planeNum]->Fit("gaus");

  PositionMeasuredY[planeNum]->Fit("pol0","","",-.8*myGhostWorld->GetTruthPlaneHalfY(),.8*myGhostWorld->GetTruthPlaneHalfY());



  TotalMomentumResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/TotalMomentumResiduals.png").c_str());
  YMomentumResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/YMomentumResiduals.png").c_str());
  ZMomentumResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/ZMomentumResiduals.png").c_str());
  YPositionResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/YPositionResiduals.png").c_str());
  ZPositionResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/ZPositionResiduals.png").c_str());

  PositionMeasuredHist[planeNum]->Draw("colz");
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/PositionMeasured.png").c_str());
  PositionTruthHist[planeNum]->Draw("colz");
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/PositionTruth.png").c_str());

  PositionMeasuredY[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/PositionMeasuredY.png").c_str());

  MeasurementInaccuracyY[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/MeasurementInaccuracyY.png").c_str());
  MeasurementInaccuracyZ[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/MeasurementInaccuracyZ.png").c_str());


  delete TotalMomentumResiduals[planeNum];
  delete YMomentumResiduals[planeNum];
  delete ZMomentumResiduals[planeNum];
  delete YPositionResiduals[planeNum];
  delete ZPositionResiduals[planeNum];

  delete PositionMeasuredHist[planeNum];
  delete PositionTruthHist[planeNum];

  delete PositionMeasuredY[planeNum];

  delete MeasurementInaccuracyY[planeNum];
  delete MeasurementInaccuracyZ[planeNum];

  }

  // for (int i = 0; i < int(ChiSquaredValues.size()); ++i){
  //   ChiSquaredHistogram->Fill(ChiSquaredValues.at(i));
  //   // ChiSquaredProbHistogram->Fill(ChiSquaredProb.at(i));
  // }


  ChiSquaredHistogram->Draw();
  // ChiSquaredProbHistogram->SetMarkerColor(5);
  // ChiSquaredProbHistogram->Draw("SAME");
  myCanvas->SaveAs("../Plots/ChiSquareds.png");



  TGraph* RMSPerPlaneY = new TGraph(maxNumPlanes, pointNo, YRMSarray);
  RMSPerPlaneY->SetTitle("RMS; Plane Number; RMS (mm)");
  // RMSPerPlaneY->GetXaxis()->SetRangeUser(0,8);
  RMSPerPlaneY->GetYaxis()->SetRangeUser(0.,.5);
  RMSPerPlaneY->SetMarkerStyle(20);
  RMSPerPlaneY->SetMarkerColor(1);
  RMSPerPlaneY->Draw("AP"); //A for "all" or something, nothing plots without it. L for connected lines between points, and P for dot style points.

  TGraph* RMSPerPlaneZ = new TGraph(maxNumPlanes, pointNo, ZRMSarray);
  RMSPerPlaneZ->SetMarkerStyle(20);
  RMSPerPlaneZ->SetMarkerColor(2);
  RMSPerPlaneZ->Draw("PSAME");

  TLegend* leg = new TLegend(0.663,0.70,0.891,0.887);
  leg->AddEntry(RMSPerPlaneY,"Y position","p");
  leg->AddEntry(RMSPerPlaneZ,"Z position","p");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

  myCanvas->SaveAs("../Plots/PlaneRMS.png");

  delete RMSPerPlaneY;
  delete RMSPerPlaneZ;
  delete leg;

  delete ChiSquaredHistogram;

  delete myCanvas;

}




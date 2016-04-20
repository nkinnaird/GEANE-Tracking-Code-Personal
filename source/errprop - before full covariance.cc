// NOTE: Coordinate systems are weird, read carefully. In my own toy world, X is the particle direction, and Y and Z are the vertical and horizontal axes respectively. 
// The tracking group uses U and V to represent straw orientations, with U angled 7.5 degrees clockwise of the vertical, and V angled 7.5 degrees counterclockwise of the vertical. 
// Since the straws actually measured position in a plane perpendicular to their length, the U and V axes in the code correspond to axes angled 7.5 degrees from the horizontal (Z) clockwise and counterclockwise respectively, with positive U and V both being positive in Z. 
// GEANT4E uses a separate U, V, and W coordinate system (with correspondingly named varibles and methods) to correspond to arbitrarily oriented detector planes. 
// To connect to my system, U = X, V = V, and W = U, with the GEANE system on the left, and the tracker/my (how I've written things into the code here) system on the right. Hope you got that.
// Nick Kinnaird, 2/10/16


#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>

#include <Eigen/Dense>

// Not all these are needed probably.
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

#include "CLHEP/Matrix/Matrix.h" // I should replace the last CHLEP pieces I'm using with Eigen, as I've done with most things.

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
// const int iverbose=4;


static const int maxNumPlanes = 33;
static const int numPasses = 4; // This will create residuals on the detector planes based on numPasses-1 passes, and residuals on the starting point based on numPasses.
static const int numTrackParams = 5;

int numEventsSkipped = 0; // Number of events skipped due to bad momentum guesses.
std::vector<int> numEventsSkippedWithNumPlanes (maxNumPlanes, 0); // Not a method. Number of events skipped that hit a certain number of planes.
int numEventsSkippedBadStartingPosition = 0; // Number of events skipped due to bad starting position.

void Initialize();
G4ErrorTarget* BuildTarget( G4int iTarget );
void ProcessEvent(int currentEventNum, int numPlanesHit, int arrayPosition );
void Finalize();
// void TrackCorrelation(G4ErrorMatrix *myTransferMatrices, G4ErrorMatrix* errorMatrices, double paramMeasured[][maxNumPlanes], double paramPredicted[][maxNumPlanes], int numPlanesHit, double deltaStartingTrack[5]);
void TrackCorrelation(std::vector<G4ErrorMatrix> myTransferMatrices, std::vector<G4ErrorMatrix> errorMatrices, std::vector <std::vector<G4double> > paramMeasured, std::vector <std::vector<G4double> > paramPredicted, int numPlanesHit, double deltaStartingTrack[5], double& eventChiSquared, int firstPlaneHit);

void ResidualPlots();

G4ErrorTarget* theTarget;
G4ErrorMode theG4ErrorMode;
G4ErrorPropagatorManager* g4emgr;
G4ErrorPropagatorData* g4edata; 


G4ThreeVector vyV(0,1,0); // Set tracing vectors to be orthogonal in YZ, and transform to UV later. Necessary because of what's in the GEANT4 source code.
G4ThreeVector wzU(0,0,1);

double rotationAngle = 7.5*pi/180;
// G4ThreeVector vyV(0, std::sin(rotationAngle), std::cos(rotationAngle)); // Non-orthogonal vectors which GEANT4 doesn't like.
// G4ThreeVector wzU(0, -std::sin(rotationAngle), std::cos(rotationAngle));


Eigen::Matrix2d YZtoVUcoordinateTransformationMatrix(2,2);
Eigen::Matrix2d YZtoVUcoordinateTransformationMatrixInverse(2,2);

Eigen::MatrixXd YZtoVUcoordinateTransformFiveByFive(5,5);
Eigen::MatrixXd YZtoVUcoordinateTransformFiveByFiveInverse(5,5);

Eigen::MatrixXd JacobianToVU(5,5);


const double noHit = 900000.; // Unphysical double value to signify the lack of a hit within a tracker plane.

/////////////////////////////////////////////////////////////////////////////////////
// My global variables
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

TH1F* UPositionResiduals[maxNumPlanes];
TH1F* VPositionResiduals[maxNumPlanes];

TH2F* PositionMeasuredHist[maxNumPlanes];
TH2F* PositionTruthHist[maxNumPlanes];
// TH1F* PositionMeasuredY[maxNumPlanes];
// TH1F* MeasurementInaccuracyY[maxNumPlanes];
// TH1F* MeasurementInaccuracyZ[maxNumPlanes];
TH1F* ChiSquaredHistogram;
TH1F* ChiSquaredHistogramPlanesHit[maxNumPlanes];
  // TH1F* ChiSquaredProbHistogram;

int numEventsWith32Planes = 0;
double materialErrorPerPlaneSum[maxNumPlanes];
double materialErrorPerPlane[maxNumPlanes];

/////////////////////////////////////////////////////////////////////////////////////

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

	// std::cout << "Event parameters are: " << std::endl; 

	// for (int i = 0; i < int(InputParameters["EventID"].size()); ++i)
	// {
	// 	std::cout << "Event ID: " << InputParameters["EventID"].at(i) << std::endl
	//             << "GlobalTime: " << InputParameters["GlobalTime"].at(i) << std::endl
	//             << "ProperTime: " << InputParameters["ProperTime"].at(i) << std::endl
	//             << "XPosition: " << InputParameters["XPosition"].at(i) << std::endl
	//             << "YPosition: " << InputParameters["YPosition"].at(i) << std::endl
	//             << "ZPosition: " << InputParameters["ZPosition"].at(i) << std::endl
	//             << "XMomentum: " << InputParameters["XMomentum"].at(i) << std::endl
	//             << "YMomentum: " << InputParameters["YMomentum"].at(i) << std::endl
	//             << "ZMomentum: " << InputParameters["ZMomentum"].at(i) << std::endl
 //              << "CopyNo: "    << InputParameters["CopyNo"].at(i) << std::endl
 //              << "UPosition: " << InputParameters["UPosition"].at(i) << std::endl
 //              << "VPosition: " << InputParameters["VPosition"].at(i) << std::endl
	//             << std::endl;
	// }

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
  double uHistBound;
  double vHistBound;

      // zHistBound = yHistBound = uHistBound = vHistBound = -1; // Set the lower hist bound > upper hist bound for the histogram to auto set the x range based on the first 1000 (or maybe 100 or 10000) events.
      // -1 values takes care of that.

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {

    if (planeNum==0) 
    { 
      zHistBound = 15; 
      yHistBound = 15;
      uHistBound = 15;
      vHistBound = 15;
    }
    else
    { 
      zHistBound = 1; 
      yHistBound = 1;
      uHistBound = 1;
      vHistBound = 1;
    }



  TotalMomentumResiduals[planeNum] = new TH1F(("Total Momentum Residuals Plane "+std::to_string(planeNum)).c_str(),"Total Momentum Residuals; Total Momentum Residual (MeV); Number of Events", 100, -300, 300); // Added a factor of 10 to all these numbers for the curved 0 plane.
  YMomentumResiduals[planeNum] = new TH1F(("Y Momentum Residuals Plane "+std::to_string(planeNum)).c_str(),"Y Momentum Residuals; Y Momentum Residual (MeV); Number of Events", 100, -50, 50);
  ZMomentumResiduals[planeNum] = new TH1F(("Z Momentum Residuals Plane "+std::to_string(planeNum)).c_str(),"Z Momentum Residuals; Z Momentum Residual (MeV); Number of Events", 100, -50, 50);
  YPositionResiduals[planeNum] = new TH1F(("Y Position Residuals Plane "+std::to_string(planeNum)).c_str(),"Y Position Residuals; Y Position Residual (mm); Number of Events",100,-yHistBound,yHistBound);
  ZPositionResiduals[planeNum] = new TH1F(("Z Position Residuals Plane "+std::to_string(planeNum)).c_str(),"Z Position Residuals; Z Position Residual (mm); Number of Events",100,-zHistBound,zHistBound);

  UPositionResiduals[planeNum] = new TH1F(("U Position Residuals Plane "+std::to_string(planeNum)).c_str(),"U Position Residuals; U Position Residual (mm); Number of Events",100,-uHistBound,uHistBound);
  VPositionResiduals[planeNum] = new TH1F(("V Position Residuals Plane "+std::to_string(planeNum)).c_str(),"V Position Residuals; V Position Residual (mm); Number of Events",100,-vHistBound,vHistBound);


  PositionMeasuredHist[planeNum] = new TH2F(("Position Measured Plane "+std::to_string(planeNum)).c_str(),"Position Measured; Z Position (mm); Y Position (mm); Number of Events",100,(-1.*myGhostWorld->GetTruthPlaneHalfZ()),(1.*myGhostWorld->GetTruthPlaneHalfZ()),100,(-1.*myGhostWorld->GetTruthPlaneHalfY()),(1.*myGhostWorld->GetTruthPlaneHalfY()));
  PositionTruthHist[planeNum] = new TH2F(("Position Truth Plane "+std::to_string(planeNum)).c_str(),"Position Truth; Z Position (mm); Y Position (mm); Number of Events",100,(-1.*myGhostWorld->GetTruthPlaneHalfZ()),(1.*myGhostWorld->GetTruthPlaneHalfZ()),100,(-1.*myGhostWorld->GetTruthPlaneHalfY()),(1.*myGhostWorld->GetTruthPlaneHalfY()));

  // PositionMeasuredY[planeNum] = new TH1F(("Position Measured Y Plane "+std::to_string(planeNum)).c_str(),"Y Measurement; Y Measurement (mm); Number of Events",100,(-1.*myGhostWorld->GetTruthPlaneHalfY()),(1.*myGhostWorld->GetTruthPlaneHalfY()));

  // MeasurementInaccuracyY[planeNum] = new TH1F(("Y Inaccuracy Plane "+std::to_string(planeNum)).c_str(),"Y Measurement Inaccuracy; Y Measurement Inaccuracy (mm); Number of Events",100,-.3,.3);
  // MeasurementInaccuracyZ[planeNum] = new TH1F(("Z Inaccuracy Plane "+std::to_string(planeNum)).c_str(),"Z Measurement Inaccuracy; Z Measurement Inaccuracy (mm); Number of Events",100,-.3,.3);

  ChiSquaredHistogramPlanesHit[planeNum] = new TH1F(("Chi Squared Planes Hit "+std::to_string(planeNum)).c_str(), "Chi Squared Values; Chi Squared; Number of Events", 120, 0, 60);

  }

  ChiSquaredHistogram = new TH1F("Chi Squared Values for all Events", "Chi Squared Values; Chi Squared; Number of Events", 120, 0, 60);
  // TH1F* ChiSquaredProbHistogram = new TH1F("Chi Squared PDF", "Chi Squared PDF; Chi Squared; Number of Events", 100, 0, 50);


  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {
    materialErrorPerPlaneSum[planeNum] = 0.;
  }


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
/////////////////////////////////////////////////////////////////////////////////////

YZtoVUcoordinateTransformationMatrix(0,0)=std::sin(rotationAngle);
YZtoVUcoordinateTransformationMatrix(0,1)=std::cos(rotationAngle);
YZtoVUcoordinateTransformationMatrix(1,0)=-std::sin(rotationAngle);
YZtoVUcoordinateTransformationMatrix(1,1)=std::cos(rotationAngle);

YZtoVUcoordinateTransformationMatrixInverse = YZtoVUcoordinateTransformationMatrix.inverse();

YZtoVUcoordinateTransformFiveByFive = Eigen::MatrixXd::Zero(5,5);
YZtoVUcoordinateTransformFiveByFive.bottomRightCorner<2,2>() = YZtoVUcoordinateTransformationMatrix;
YZtoVUcoordinateTransformFiveByFiveInverse = Eigen::MatrixXd::Zero(5,5);
YZtoVUcoordinateTransformFiveByFiveInverse.bottomRightCorner<2,2>() = YZtoVUcoordinateTransformationMatrixInverse;

// G4cout << YZtoVUcoordinateTransformationMatrix << G4endl << YZtoVUcoordinateTransformationMatrixInverse << G4endl << YZtoVUcoordinateTransformFiveByFive << G4endl << YZtoVUcoordinateTransformFiveByFiveInverse << G4endl;

JacobianToVU = Eigen::MatrixXd::Zero(5,5);
JacobianToVU(0,0) = 1.;
JacobianToVU.block(1,1,2,2) = YZtoVUcoordinateTransformationMatrixInverse;
JacobianToVU.block(3,3,2,2) = YZtoVUcoordinateTransformationMatrixInverse;
// JacobianToVU = Eigen::MatrixXd::Identity(5,5);

G4cout << G4endl << "Jacobian between YZ and VU coordinate systems is: " << G4endl << JacobianToVU << G4endl;



/////////////////////////////////////////////////////////////////////////////////////  

  Initialize();

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

                  //Test statement for checking non-orthogonal Jacobian elements.

  G4ThreeVector xv3(-30., 100., -45.);
  G4ThreeVector pv3(900., 50., -200. );

  G4ErrorTrajErr fillerror( 5, 0 );
  G4ErrorFreeTrajState* testFreeState = new G4ErrorFreeTrajState("e+", xv3, pv3, fillerror );


  G4ErrorMatrix YZtransfMatrix;
  G4ErrorSurfaceTrajState* YZtestSurfaceState = new G4ErrorSurfaceTrajState(*testFreeState, vyV, wzU, YZtransfMatrix);


  G4ThreeVector vvv(0, std::sin(rotationAngle), std::cos(rotationAngle)); 
  G4ThreeVector uuu(0, -std::sin(rotationAngle), std::cos(rotationAngle));

  // G4cout << G4endl << "XXXXXXXXXXXXXXXX" << G4endl;
  G4ErrorMatrix VUtransfMatrix;
  G4ErrorSurfaceTrajState* VUtestSurfaceState = new G4ErrorSurfaceTrajState(*testFreeState, vvv, uuu, VUtransfMatrix);
  // G4cout << "Modified source code right above here ^^^ that I want to compare to below: " << G4endl;

  Eigen::MatrixXd YZtransfEigen(5,5);
  YZtransfEigen = Eigen::MatrixXd::Zero(5,5);

   for(int i=0;i<5;i++) {
     for(int j=0; j<5;j++) {
          YZtransfEigen(i,j) = YZtransfMatrix[i][j];
     }
   } 

  Eigen::MatrixXd YZtransfWithJacobian = YZtransfEigen.inverse() * JacobianToVU;
  // G4cout << "Above elements compared to: " << YZtransfWithJacobian(3,3) << " and " << YZtransfWithJacobian(3,4) << " and " << YZtransfWithJacobian(4,3) << " and " << YZtransfWithJacobian(4,4) << G4endl << G4endl;
  G4cout << "YZ transform matrix multiplied by VU Jacobian: " << G4endl << YZtransfWithJacobian << G4endl;



// std::exit(0);


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


  // EventID array within the map goes like: 0, 0, 0, 1, 1, 2, 3, 3, 4, 4, 4, 4 ...
  // Need to get last value to see how many events I need to run, and then determine the amount of planes crossed for each event, and pass those numbers to ProcessEvent.
  int nEvents = InputParameters["EventID"].back(); // Access last element.
  

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

      if (numPlanesHit >= 9) // Particles always hit plane 0. 9 for the 0 plane hit, and 2 full module hits (or an equivalent number of planes).
      {
        ProcessEvent(currentEventNum, numPlanesHit, arrayPos ); // Process current event.
      }
      else{
        numEventsSkippedMissingPlanes++;
      }

  }

  

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {
    materialErrorPerPlane[planeNum] = sqrt(materialErrorPerPlaneSum[planeNum]/numEventsWith32Planes) * 1e4; // Should go from cm^2 to um with this.
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

  // g4edata->SetVerbose(5);

  g4emgr->SetUserInitialization(myPhysicalWorld); 

  g4emgr->InitGeant4e();

  G4UImanager::GetUIpointer()->ApplyCommand("/control/verbose 1");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 100 mm");
  // G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/iverbose 4"); // GEANT4E specific verbose command, doesn't seem to work. Not sure how to set the iverbose level. Have to also set G4EVERBOSE in the env variables to print verbose statements.

  // char* mode = getenv("G4ERROR_MODE");
  // if( mode ) {
  //   if( G4String(mode) == G4String("FORWARDS") ) {
  //     theG4ErrorMode = G4ErrorMode_PropForwards;
  //   } else if( G4String(mode) == G4String("BACKWARDS") ) {
  //     theG4ErrorMode = G4ErrorMode_PropBackwards;
  //   } else {
  //     G4Exception("exG4eReco","Fatal error in Argument",
  //       FatalErrorInArgument,
  //       G4String("Variable G4ERROR_MODE = " + G4String(mode) + 
  //                "   It must be:  FORWARDS or BACKWARDS").c_str());
  //   }
  // } else {
  //   G4Exception("exG4eReco","Fatal error in Argument",
  //     JustWarning,"Variable G4ERROR_MODE not defined, taking it = BACKWARDS");
  // } 

 

}


void ProcessEvent(int currentEventNum, int numPlanesHit, int arrayPosition)
{

  G4int iProp = 1; // Set propogation to be step by step. iProp = 0 is until target.

  G4cout << G4endl << " Start Event number: " << currentEventNum << G4endl << G4endl;

  // Initialization of vectors, arrays, matrices, etc. 
  static TRandom3 r3;

  std::vector <std::vector<G4double> > trackParamMeasured;  // What is given from input root file.
  std::vector <std::vector<G4double> > trackParamPredicted; // What GEANT4E calculates as the average from propagation.

  std::vector <std::vector<G4double> > trackParamMeasuredExtra;
  std::vector <std::vector<G4double> > trackParamPredictedExtra; // Extra vector for non-measured track parameters that I might want to compare later.
  std::vector <std::vector<G4double> > trackParamPredictedMomentumExtra;

  // Start with noHit arrays of maxNumPlanes size for track vectors, then later erase noHit entries before sending to Track Correlation method.
  trackParamMeasured.resize(numTrackParams, std::vector<G4double>(maxNumPlanes, noHit)); // Initialize 2D vector with 5 rows and 9 columns filled with noHit values.
  trackParamPredicted.resize(numTrackParams, std::vector<G4double>(maxNumPlanes, noHit)); // Initialize this only with the number of planes hit for the number of columns.


  trackParamMeasuredExtra.resize(3, std::vector<G4double>(maxNumPlanes, noHit));
  trackParamPredictedExtra.resize(3, std::vector<G4double>(maxNumPlanes, noHit)); // For now just fill this with xyz recorded values.
  trackParamPredictedMomentumExtra.resize(3, std::vector<G4double>(maxNumPlanes, noHit));

  // Fill the measured values array here, with noHit values for 1/p, py/p, and pz/p, and some small uncertainty in y z. Units are MeV mm.

  int increment = 0; // increment along copyNo array-piece
  int firstPlaneHit = 0; // That's not the 0 plane.

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum) // This planeNum stands for actual plane number.
  {

    trackParamMeasured[0][planeNum] = noHit;
    trackParamMeasured[1][planeNum] = noHit;
    trackParamMeasured[2][planeNum] = noHit;

    // G4cout << "Test increment: " << increment << " Copy no: " << InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1) << " First plane hit number: " << firstPlaneHit << G4endl;

    // Array position is the last position in the array for an events data. For the 0 plane, array position is numPlanesHit-1, so I have to add 1 to start at 0. Similarly for the rest of the planes.



    if (planeNum == InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1))
    {
        PlaneParameterTruth[planeNum]["UPos"] = (InputParameters["UPosition"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["VPos"] = (InputParameters["VPosition"].at(arrayPosition-numPlanesHit+increment+1));

        if (PlaneParameterTruth[planeNum]["UPos"] == noHit)
        { trackParamMeasured[4][planeNum] = noHit; } // Don't smear the noHit value, as that messes things up later.
        else { trackParamMeasured[4][planeNum] = PlaneParameterTruth[planeNum]["UPos"]+r3.Gaus(0,.1); } //+0.2*(r3.Rndm()-0.5); // +- 100 RMS microns in the position measurements to simulate smearing. 

        if (PlaneParameterTruth[planeNum]["VPos"] == noHit)
        { trackParamMeasured[3][planeNum] = noHit; }
        else { trackParamMeasured[3][planeNum] = PlaneParameterTruth[planeNum]["VPos"]+r3.Gaus(0,.1); } //+0.2*(r3.Rndm()-0.5);

      
        PlaneParameterTruth[planeNum]["XPos"] = (InputParameters["XPosition"].at(arrayPosition-numPlanesHit+increment+1)); // increment vs planeNum
        PlaneParameterTruth[planeNum]["YPos"] = (InputParameters["YPosition"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["ZPos"] = (InputParameters["ZPosition"].at(arrayPosition-numPlanesHit+increment+1));

        trackParamMeasuredExtra[1][planeNum] = PlaneParameterTruth[planeNum]["YPos"];//+r3.Gaus(0,.1); // These are the simulated (but not actually measured) YZ values.
        trackParamMeasuredExtra[2][planeNum] = PlaneParameterTruth[planeNum]["ZPos"];//+r3.Gaus(0,.1);


        PlaneParameterTruth[planeNum]["XMom"] = (InputParameters["XMomentum"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["YMom"] = (InputParameters["YMomentum"].at(arrayPosition-numPlanesHit+increment+1));
        PlaneParameterTruth[planeNum]["ZMom"] = (InputParameters["ZMomentum"].at(arrayPosition-numPlanesHit+increment+1));

        PlaneParameterTruth[planeNum]["CopyNo"] = (InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1));

        if (firstPlaneHit == 0) // Make sure the previous hit plane was plane 0, and then don't change it after it's been changed once.
        { firstPlaneHit = int(InputParameters["CopyNo"].at(arrayPosition-numPlanesHit+increment+1)); }

        if (increment < numPlanesHit-1) // -1 to make sure that it doesn't increment over the size of the vector for the above if statement check
        { increment++; } //only increase increment if a plane is hit
        
    }
    
    else // Fill with noHit entries when planes/detectors are missed in the data generation. 
    {
        trackParamMeasured[3][planeNum] = noHit;
        trackParamMeasured[4][planeNum] = noHit;

        trackParamMeasuredExtra[1][planeNum] = noHit;
        trackParamMeasuredExtra[2][planeNum] = noHit;

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

    // Make these into U and VPos keys.
    PlanePositionMeasured[planeNum]["UPos"] = (trackParamMeasured[4][planeNum]);
    PlanePositionMeasured[planeNum]["VPos"] = (trackParamMeasured[3][planeNum]);

    PlanePositionMeasured[planeNum]["YPos"] = (trackParamMeasuredExtra[1][planeNum]); // Just fill with UV recorded values for now I guess.
    PlanePositionMeasured[planeNum]["ZPos"] = (trackParamMeasuredExtra[2][planeNum]);



  }


  // G4cout << " Filled into parameter truth arrays: " << G4endl;
  // for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  // {
  //   std::cout << "Copy No: " << PlaneParameterTruth[planeNum]["CopyNo"] << std::endl
  //             << "XPosition: " << PlaneParameterTruth[planeNum]["XPos"] << std::endl
  //             << "YPosition: " << PlaneParameterTruth[planeNum]["YPos"] << std::endl
  //             << "ZPosition: " << PlaneParameterTruth[planeNum]["ZPos"] << std::endl
  //             // << "XMomentum: " << PlaneParameterTruth[planeNum]["XMom"] << std::endl
  //             // << "YMomentum: " << PlaneParameterTruth[planeNum]["YMom"] << std::endl
  //             // << "ZMomentum: " << PlaneParameterTruth[planeNum]["ZMom"] << std::endl
  //             << "UPosition: " << PlaneParameterTruth[planeNum]["UPos"] << std::endl
  //             << "VPosition: " << PlaneParameterTruth[planeNum]["VPos"] << std::endl
  //             << std::endl;
  // }

  // G4cout << " Filled into parameter measured arrays: " << G4endl;
  // for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  // {
  //   std::cout << "Copy No: " << PlaneParameterTruth[planeNum]["CopyNo"] << std::endl
  //             << "YPosition: " << PlanePositionMeasured[planeNum]["YPos"] << std::endl
  //             << "ZPosition: " << PlanePositionMeasured[planeNum]["ZPos"] << std::endl
  //             << "UPosition: " << PlanePositionMeasured[planeNum]["UPos"] << std::endl
  //             << "VPosition: " << PlanePositionMeasured[planeNum]["VPos"] << std::endl
  //             << std::endl;
  // }





  // G4double totalMomentum = sqrt( (PlaneParameterTruth[firstPlaneHit]["XMom"] * PlaneParameterTruth[firstPlaneHit]["XMom"])
  //                               +(PlaneParameterTruth[firstPlaneHit]["YMom"] * PlaneParameterTruth[firstPlaneHit]["YMom"])
  //                               +(PlaneParameterTruth[firstPlaneHit]["ZMom"] * PlaneParameterTruth[firstPlaneHit]["ZMom"]) );
  


// Set the initial starting trajectory.
// These units should still be in MeV mm and not GeV cm. Only the error matrix should be in GeV cm.
// Take the starting plane values, and add changes to similulate a non-perfect starting trajectory from some initial weak trackfitting algorithm.

  G4double initialXPosChange = -.01; // The starting plane is now right in front of the first plane hit and uses its recorded values as it's starting parameters. Making this too low causes errors.
  // Can make it higher actually with all of the fixes I've made to the rest of the code, but I probably want to set as is since I set the starting guess to be close to the first hit measured offset by something, and I will do the same in the final code.

  // G4double initialUPosChange = 0.;//-5.;
  // G4double initialVPosChange = 0.;//-3.;

  // G4double initialYPosChange = initialVPosChange;
  // G4double initialZPosChange = initialUPosChange;

  // G4double initialYPosChange = YZtoVUcoordinateTransformationMatrixInverse(0,0)*initialVPosChange + YZtoVUcoordinateTransformationMatrixInverse(0,1)*initialUPosChange;
  // G4double initialZPosChange = YZtoVUcoordinateTransformationMatrixInverse(1,0)*initialVPosChange + YZtoVUcoordinateTransformationMatrixInverse(1,1)*initialUPosChange;

  G4double initialYPosChange = 40.*(r3.Rndm()-0.5); 
  G4double initialZPosChange = 40.*(r3.Rndm()-0.5); 

  // G4cout << "Initial y and z pos changes: " << initialYPosChange << " " << initialZPosChange << G4endl; 

  // G4cout << " U, V, Y, and Z starting offsets: " << G4endl
  // << " U: " << initialUPosChange << G4endl
  // << " V: " << initialVPosChange << G4endl
  // << " Y: " << initialYPosChange << G4endl
  // << " Z: " << initialZPosChange << G4endl << G4endl;

  G4double startingXPos = PlaneParameterTruth[firstPlaneHit]["XPos"] + initialXPosChange;
  G4double startingYPos = PlaneParameterTruth[firstPlaneHit]["YPos"] + initialYPosChange;
  G4double startingZPos = PlaneParameterTruth[firstPlaneHit]["ZPos"] + initialZPosChange;


  G4double initialXMomChange = 200.*(r3.Rndm()-0.5);
  G4double initialYMomChange = 200.*(r3.Rndm()-0.5);
  G4double initialZMomChange = 200.*(r3.Rndm()-0.5);

  G4double startingXMom = 1.*PlaneParameterTruth[firstPlaneHit]["XMom"] + initialXMomChange;
  G4double startingYMom = 1.*PlaneParameterTruth[firstPlaneHit]["YMom"] + initialYMomChange;
  G4double startingZMom = 1.*PlaneParameterTruth[firstPlaneHit]["ZMom"] + initialZMomChange;

  // Leave this in for if the source code is every fixed.
  // G4double startingVMom;
  // G4double startingUMom;


// Variables for track parameter residual determination.
double deltaStartingTrack[5];
G4double totalStartMomentumModified;
double totalMomentumResidual;

double UPositionResidual;
double VPositionResidual;
double YPositionResidual;
double ZPositionResidual;

double YMomentumResidual;
double ZMomentumResidual;

double eventChiSquared;


  totalStartMomentumModified = sqrt ((startingXMom)*(startingXMom) 
                              + (startingYMom)*(startingYMom)
                              + (startingZMom)*(startingZMom));


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    trackParamPredicted[i].clear(); // Clear at the beginning of each pass.
    trackParamPredicted[i].resize(maxNumPlanes, noHit);
  }

  double reconstructedTotalMomentum;


  // if (startingYPos > .9*myGhostWorld->GetTruthPlaneHalfY()) {startingYPos = .9*myGhostWorld->GetTruthPlaneHalfY();} // These are to make sure the traced particle doesn't start outside of the detector or world volume. 
  // else if (startingYPos < -.9*myGhostWorld->GetTruthPlaneHalfY()) {startingYPos = -.9*myGhostWorld->GetTruthPlaneHalfY();} // .9 so that it's far enough away from the edge so it doesn't seg fault if no planes are hit during the traceback.
  // if (startingZPos > .9*myGhostWorld->GetTruthPlaneHalfZ()) {startingZPos = .9*myGhostWorld->GetTruthPlaneHalfZ();}
  // else if (startingZPos < -.9*myGhostWorld->GetTruthPlaneHalfZ()) {startingZPos = -.9*myGhostWorld->GetTruthPlaneHalfZ();}

  G4ThreeVector xv3( startingXPos, startingYPos, startingZPos );
  G4ThreeVector pv3( startingXMom, startingYMom, startingZMom );

  // Leave this in for if the source code is ever fixed.
  // startingVMom = YZtoVUcoordinateTransformationMatrix(0,0)*startingYMom + YZtoVUcoordinateTransformationMatrix(0,1)*startingZMom; // Form the U and V momenta at the start of each pass from Y and Z.
  // startingUMom = YZtoVUcoordinateTransformationMatrix(1,0)*startingYMom + YZtoVUcoordinateTransformationMatrix(1,1)*startingZMom;

  // std::cout << std::endl << "Pass number: " << passNumber+1 << std::endl;


     theG4ErrorMode = G4ErrorMode_PropForwards; // Propagate forward from plane to plane.

  // Need to provide error matrix at the start of the trajectory state. Units for the error (and transport) matrices are in GeV cm.
  G4ErrorTrajErr error( 5, 0 ); 
  // error(4,4) = 0.0001; // Parentheses instead of brackets start at 1,1 instead of 0,0 for G4ErrorMatrix.
  // error(5,5) = 0.0001; 
  // I think I comment these out and start with a 0 error matrix, and add in the experimental part later. 
  // This is NOT inverted first when plugging into the initialization of the trajectory state, these values should be cm^2 and not (cm^2)^-1.
  // So for .1mm resolution or .01cm I should have .0001 values for (4,4) and (5,5). 


/////////////////////////////////////////////////////////////////////////////////////
  // This for a free trajectory state with 1/p, lambda, phi, y perp, and z perp.

  G4ErrorFreeTrajState* myFreeTrajState 
    = new G4ErrorFreeTrajState("e+", xv3, pv3, error );

  std::vector<G4ErrorMatrix> SC2SDTransformationMatrix; // Matrices for transfmormation to and from free system and plane/momentum system. These will be filled with per-step Jacobian matrices, not per-plane matrices.
  std::vector<G4ErrorMatrix> SC2SDTransformationMatrixInverse;

  SC2SDTransformationMatrix.clear();
  SC2SDTransformationMatrixInverse.clear();

  G4ErrorMatrix fillTransformationMatrix0; // 0th step transformation matrix.

  // This for a surface trajectory state with 1/p, v', w', v, and w. This needs to be orthogonal as currently written in the GEANT4 source code.
  // This creates the surface trajectory state from the free state.
  G4ErrorSurfaceTrajState* mySurfaceTrajState
    = new G4ErrorSurfaceTrajState(*myFreeTrajState, vyV, wzU, fillTransformationMatrix0); // Passing it the empty transformation matrix here fills it with that matrix for the starting plane point.

   SC2SDTransformationMatrix.push_back(fillTransformationMatrix0);

  int failure; // Int to pass and see if the inversion failed.
  SC2SDTransformationMatrixInverse.push_back(SC2SDTransformationMatrix[0].inverse(failure)); // In this case, the 0 here does not refer to the 0 plane, but the starting point plane. This is necessary for the transport matrix coordinate transformation later.


  g4emgr = G4ErrorPropagatorManager::GetErrorPropagatorManager();


/////////////////////////////////////////////////////////////////////////////////////

  G4ErrorMatrix zeroMatrix = myFreeTrajState->GetTransfMat(); // Transport matrix for the free state when no steps have been made is a zero matrix.
  // G4cout << "Check that is an initial zero matrix: " << zeroMatrix << G4endl;

  // printf("initialize matrices \n");  // Fill with 0 matrices to start with. 
  transferMatrices.push_back(zeroMatrix); // Have to fill with an object to start otherwise it seg faults. Surface transfer matrices can also just be filled with the 0 matrix.
  transferMatricesFreeSystem.push_back(zeroMatrix);
  // Free state system has this particular method while the surface system does not. On plane 0 it is a 0 matrix which will be excluded in the GEANE formulation later on.
  // There is no transfer matrix to plane 0 because there is no plane before it.
  errorMatrices.push_back(mySurfaceTrajState->GetError()); // Starts off with 0 matrix as well, but writing differently to have things consistent..

  for (int planeNum = 1; planeNum < maxNumPlanes; ++planeNum)
      {
        // G4cout << G4endl << "Tracking to plane: " << planeNum << G4endl;

/////////////////////////////////////////////////////////////////////////////////////

        // Skip over non-hit planes in front of the first plane hit, which is just before where I start the tracking. (Results seem to be a bit better that way.)
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
          G4double truthPlaneTargetX=0.;
          G4double truthPlaneTargetY=0.;
          G4double truthPlaneTargetZ=0.;

          // Have to do some arithmatic to determine how to change truthPlaneTargetX with a single for-loop for 32 planes.

          truthPlaneTargetX += int((planeNum-1)/4.) * myPhysicalWorld->GetmoduleSeparation(); // Add module separation.

          bool layerSep = FALSE;
          if (planeNum % 2 == 0)
          {
            truthPlaneTargetX += myPhysicalWorld->GetlayerSeparation(); // Add layer separation.
            layerSep = TRUE;
          }

          double intpart;
          double fractpartone = std::modf(planeNum/4. , &intpart);
          double fractparttwo = std::modf((planeNum+1)/4. , &intpart);
          bool UVSep = FALSE;
          if (fractpartone == 0 || fractparttwo == 0)
          {
            truthPlaneTargetX += myPhysicalWorld->GetUVSeparation(); // Add UV separation.
            UVSep = TRUE;
          }

          // G4cout << "Plane num is: " << planeNum << " module separation is: " << int((planeNum-1)/4.) << " layer separation is: " << layerSep << " UV separation is: " << UVSep << G4endl;

          
          // G4cout << "Starting position of free state: " << myFreeTrajState->GetPosition() << G4endl;
          // G4cout << "Truth plane target X is: " << truthPlaneTargetX << G4endl;

          if (myFreeTrajState->GetPosition().x() >= truthPlaneTargetX)
          { 
            numEventsSkippedBadStartingPosition++;
            return; 
          } // If the Geant SD recorded a too large of a starting point (323.1 mm instead of 323 mm), skip the event.
          


          G4Normal3D surfNorm(-1.,0.,0.);       
          G4Point3D surfPos(truthPlaneTargetX,truthPlaneTargetY,truthPlaneTargetZ);

          theTarget = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
          g4edata->SetTarget( theTarget );
/////////////////////////////////////////////////////////////////////////////////////
          if( iProp == 0){ // This fails as only the transfer matrix for the last step gets saved, and if it takes multiple steps to get to the target, then the full transfer matrix from the last plane to the current one does not get properly saved.
                                // See lines 566 - 573 in G4ErrorFreeTrajState.cc. Looks like the proper code (or part of it) is commented out.
            // Propagate until G4ErrorTarget is found all in one go.
              
             g4emgr->Propagate( myFreeTrajState, theTarget, theG4ErrorMode ); // If I come back to use this, I will need to intregrate the surfacetrajstate stuff in properly and adapt to the many changes I've made.
             transferMatrices[planeNum] = myFreeTrajState->GetTransfMat(); // Same for here.
          } 
/////////////////////////////////////////////////////////////////////////////////////
          else if( iProp == 1){
          // Propagate until G4ErrorTarget is reached step by step.


          g4emgr->InitTrackPropagation();

          bool moreEvt = TRUE;
          while( moreEvt )
          {
            
            // G4cout << " Particle momentum is: " << myFreeTrajState->GetMomentum() << " with magnitude: " << myFreeTrajState->GetMomentum().mag() << G4endl;
            
            g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode );
            
            // Use step by step propagation to properly save the transfer matrices.

            // G4cout << "Propagating to Plane number: " << planeNum << " stepNumBetweenPlanes: " << stepNumBetweenPlanes << " Step length: " << myFreeTrajState->GetG4Track()->GetStepLength() << G4endl;
            G4ErrorMatrix fillTransformationMatrix; // Remake this object everytime for each step before adding to the vector.
            mySurfaceTrajState = new G4ErrorSurfaceTrajState(*myFreeTrajState, vyV, wzU, fillTransformationMatrix); // Create a new surface trajectory state each step of the way from the propagating free trajectory state. The passed empty matrix gets filled with the Jacobian transformation matrix.

            SC2SDTransformationMatrix.push_back(fillTransformationMatrix);

            int failureInt;
            SC2SDTransformationMatrixInverse.push_back((SC2SDTransformationMatrix.back()).inverse(failureInt));

            // G4cout << " Transformation matrix " << SC2SDTransformationMatrix.back() << G4endl
            // << " Transformation matrix inverse " << SC2SDTransformationMatrixInverse.back() << G4endl;

            // printf("Check SC2SD inversion worked properly \n");
            // std::cout << SC2SDTransformationMatrix.back()*SC2SDTransformationMatrixInverse.back();

            // Transfer matrices uses push_back since they are only filled for hit planes.
                if(stepNumBetweenPlanes==0) {   // NOTE: The inverse of the Jacobian is on the right side in this formulation because of how the SC2SD matrix is generated in the first place. (Reverse from how one might do it on paper.)
                  // For matrix multiplication, it should be: R10 = A1 * T10 * A0^-1, where R and T are the transport matrices in separate coordinate systems from plane 0 to 1, and A is the transformation between the two. 
                  transferMatricesFreeSystem.push_back(myFreeTrajState->GetTransfMat()); // Gets the transport or transfer matrix for the last step. If there is only a single step between planes this is fine.
                  transferMatrices.push_back(SC2SDTransformationMatrix.back()*(myFreeTrajState->GetTransfMat())*SC2SDTransformationMatrixInverse.at(SC2SDTransformationMatrixInverse.size()-1));
                }
                else {
                  G4ErrorMatrix fillTransportMatrix;
                  fillTransportMatrix = transferMatricesFreeSystem.back(); // Get the last transport matrix of previous steps.
                  transferMatricesFreeSystem.pop_back(); // Delete the last matrix within the vector which equates to the last step, and then fill with the accumulated steps version.
                  transferMatricesFreeSystem.push_back(myFreeTrajState->GetTransfMat() * fillTransportMatrix); // Accumulates transport matrices from step to step for steps from last plane to the next plane.
                  
                  transferMatrices.pop_back();
                  transferMatrices.push_back(SC2SDTransformationMatrix.back()*(transferMatricesFreeSystem.back())*SC2SDTransformationMatrixInverse.at(SC2SDTransformationMatrixInverse.size()-1)); // Have to be careful with plane to plane transport vs step by step transport. This is the right way to do it.
                    // Here I gather up the transport matrices for the free system for all steps between planes, then I multiply on the right by the previous plane coordinate transformation matrix, and on the left by the next step coordinate transformation matrix. It keeps looping until the next step transformation matrix is the same as the next plane transformation matrix.
                  
                }


            stepNumBetweenPlanes++; // Increment the step number between planes. (The 1st step is the same as stepNumBetweenPlanes==0.)


            //---- Check if target is reached
            if( g4emgr->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() )) 
            {
              g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );  
              moreEvt = 0;
              // G4cout << "STEP_BY_STEP propagation: Last Step " << G4endl;
            }

            delete mySurfaceTrajState; // Move through the steps up to the target, deleting and remaking the surfacetrajstate every step of the way. Perhaps don't actually have to delete, might be automatic.
          } // end of while


          }// end of else if

/////////////////////////////////////////////////////////////////////////////////////
            // After getting to the target, once again create the surface trajectory state on the target.


            G4ErrorMatrix unneededMatrix;
            mySurfaceTrajState
              = new G4ErrorSurfaceTrajState(*myFreeTrajState, vyV, wzU, unneededMatrix); // The idea is to propagate in steps with the free state, then at the end on the plane turn that into a surface state and grab associated parameters and the error.

            trackParamPredicted[0].at(planeNum) = mySurfaceTrajState->GetParameters().GetInvP();

            G4double statePx = sqrt( 1./trackParamPredicted[0].at(planeNum)*1./trackParamPredicted[0].at(planeNum) - mySurfaceTrajState->GetParameters().GetPV()*mySurfaceTrajState->GetParameters().GetPV() - mySurfaceTrajState->GetParameters().GetPW()*mySurfaceTrajState->GetParameters().GetPW() );
            
            trackParamPredicted[1].at(planeNum) = (mySurfaceTrajState->GetParameters().GetPV()) / statePx;
            trackParamPredicted[2].at(planeNum) = (mySurfaceTrajState->GetParameters().GetPW()) / statePx;

            trackParamPredicted[3].at(planeNum) = mySurfaceTrajState->GetParameters().GetV(); 
            trackParamPredicted[4].at(planeNum) = mySurfaceTrajState->GetParameters().GetW();

            trackParamPredictedExtra[0].at(planeNum) = mySurfaceTrajState->GetPosition().x(); // Extra track parameters to grab the global xyz values for checks.
            trackParamPredictedExtra[1].at(planeNum) = mySurfaceTrajState->GetPosition().y(); // If the surface state is generated in a non-orthogonal frame, it's unlikely but possible that these positions will get messed up. In that case just grab these global values from the free state instead.
            trackParamPredictedExtra[2].at(planeNum) = mySurfaceTrajState->GetPosition().z();

            trackParamPredictedMomentumExtra[0].at(planeNum) = statePx;
            trackParamPredictedMomentumExtra[1].at(planeNum) = mySurfaceTrajState->GetParameters().GetPV(); // Extra track parameters to get the different momenta for easier use.
            trackParamPredictedMomentumExtra[2].at(planeNum) = mySurfaceTrajState->GetParameters().GetPW(); // Are currently XYZ momenta for XYZ tracing.

            // G4cout << " Predicted Track Parameters Plane: " << planeNum << G4endl
            // << " param 0: " << trackParamPredicted[0].at(planeNum) << G4endl
            // << " param 1: " << trackParamPredicted[1].at(planeNum) << G4endl
            // << " param 2: " << trackParamPredicted[2].at(planeNum) << G4endl
            // << " param 3: " << trackParamPredicted[3].at(planeNum) << G4endl
            // << " param 4: " << trackParamPredicted[4].at(planeNum) << G4endl
            // << " Y position: " << trackParamPredictedExtra[1].at(planeNum) << G4endl
            // << " Z position: " << trackParamPredictedExtra[2].at(planeNum) << G4endl 
            // << G4endl;
 

            G4ErrorTrajErr errorEnd;
            errorEnd = mySurfaceTrajState->GetError();
            errorMatrices.push_back(errorEnd);



      }  // End planeNum loop


/////////////////////////////////////////////////////////////////////////////////////



      if (passNumber == 2 && numPlanesHit == 33)
      {
          numEventsWith32Planes++;

          for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
          {
            materialErrorPerPlaneSum[planeNum] = materialErrorPerPlaneSum[planeNum] + errorMatrices[planeNum][4][4];
          }

      }


/////////////////////////////////////////////////////////////////////////////////////


      TrackCorrelation(transferMatrices, errorMatrices, trackParamMeasured, trackParamPredicted, numPlanesHit, deltaStartingTrack, eventChiSquared, firstPlaneHit);
      delete mySurfaceTrajState;
      delete myFreeTrajState; // Might not need to have these here as I think they might die anyways once the loop cycles.


      // Pull out the new filled delta starting track to get a new starting track. Could possibly transform deltaStartingTrack back to YZ before using it here.

      reconstructedTotalMomentum = 1./(deltaStartingTrack[0]+(1./totalStartMomentumModified));

/////////////////////////////////////////////////////////////////////////////////////

      double pvpOVERpup = deltaStartingTrack[1] + startingYMom/startingXMom;
      double pwpOVERpup = deltaStartingTrack[2] + startingZMom/startingXMom;
      startingXMom = sqrt( reconstructedTotalMomentum*reconstructedTotalMomentum / (1 + pvpOVERpup*pvpOVERpup + pwpOVERpup*pwpOVERpup) ); // This is only correct for an orthogonal system.

      // Leave in for if the source code is ever fixed.
      // startingVMom = pvpOVERpup * startingXMom; // This is the new startingXMom.
      // startingUMom = pwpOVERpup * startingXMom;
      // startingYMom = YZtoVUcoordinateTransformationMatrixInverse(0,0)*startingVMom + YZtoVUcoordinateTransformationMatrixInverse(0,1)*startingUMom;
      // startingZMom = YZtoVUcoordinateTransformationMatrixInverse(1,0)*startingVMom + YZtoVUcoordinateTransformationMatrixInverse(1,1)*startingUMom;

      startingYMom = pvpOVERpup * startingXMom; // This is the new startingXMom.
      startingZMom = pwpOVERpup * startingXMom;

      // G4cout << " Check momenta are equal: " << reconstructedTotalMomentum << " vs " << sqrt ((startingXMom)*(startingXMom) 
      //                                                                                       + (startingYMom)*(startingYMom)
      //                                                                                       + (startingZMom)*(startingZMom)) << G4endl;

      startingXPos = startingXPos + 0.; // There are no returned variables that allow for changing the starting x position. 
      startingYPos = startingYPos + deltaStartingTrack[3];
      startingZPos = startingZPos + deltaStartingTrack[4];


/////////////////////////////////////////////////////////////////////////////////////

        // NOTE: It's possible that the Geant4 thrown exceptions about particles not reaching the target come from bad momentum guesses that aren't bad enough to trigger the following if statement condition.
        // If so, I am incorrectly accounting for the number of skipped events.

        // G4cout << "Total momentum modified before change: " << totalStartMomentumModified << G4endl;
        // G4cout << "1 over tmm before change: " << 1./totalStartMomentumModified << " compared to delta 0: " << deltaStartingTrack[0] << G4endl;

        int numPlanesHitForThisEvent=0;
        if(std::abs(deltaStartingTrack[0]) > std::abs(1./totalStartMomentumModified) || reconstructedTotalMomentum*reconstructedTotalMomentum < ((startingYMom*startingYMom) + (startingZMom*startingZMom)))
        {
          G4cout << "Event skipped because a terrible fix for the momentum was made for some reason. " << G4endl;

          G4cout << "Total momentum modified before change: " << totalStartMomentumModified << G4endl;
          G4cout << "1 over tmm before change: " << 1./totalStartMomentumModified << " compared to delta 0: " << deltaStartingTrack[0] << G4endl;
          G4cout << "Reconstructed total momentum^2:  " << reconstructedTotalMomentum*reconstructedTotalMomentum << " compared to Ymom^2 + Zmom^2: " << (startingYMom*startingYMom) + (startingZMom*startingZMom) << G4endl;

          G4cout << "Plane nums hit: ";
          for (int i = 0; i < maxNumPlanes; ++i)
          {
            if (PlanePositionMeasured[i]["UPos"] != noHit || PlanePositionMeasured[i]["VPos"] != noHit)
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

        // G4cout << "total momentum before change: " << totalStartMomentumModified;
        // Need to refill totalStartMomentumModified for when the for loop gets back to the reconstructedTotalMomentum step. This properly incorporates deltaStartingTrack[0].
        totalStartMomentumModified = sqrt ((startingXMom)*(startingXMom) // Could also just set it equal to reconstructedTotalMomentum.
                               + (startingYMom)*(startingYMom)
                               + (startingZMom)*(startingZMom));
        // G4cout << " compared to total momentum after change : " << totalStartMomentumModified << G4endl;


        // G4cout << G4endl << " Delta starting track passed out: " << G4endl // This deltaStartingTrack is now actually back in YZ space.
        // << " 1/P: " << deltaStartingTrack[0] << G4endl
        // << " Pu/Px: " << deltaStartingTrack[2] << G4endl 
        // << " Pv/Px: " << deltaStartingTrack[1] << G4endl
        // << " U: " << deltaStartingTrack[4] << G4endl
        // << " V: " << deltaStartingTrack[3] << G4endl << G4endl;

        // G4cout << " New starting XPos: " << startingXPos << G4endl
        //        << " New starting YPos: " << startingYPos << G4endl
        //        << " New starting ZPos: " << startingZPos << G4endl
        //        << " New starting XMom: " << startingXMom << G4endl
        //        << " New starting YMom: " << startingYMom << G4endl
        //        << " New starting ZMom: " << startingZMom << G4endl;

/////////////////////////////////////////////////////////////////////////////////////        

 // G4cout << G4endl << "End of pass: " << passNumber+1  << G4endl << G4endl;
} // end of numPasses for loop

      
    // ChiSquaredProb.push_back(TMath::Prob(eventChiSquared, 10));
    ChiSquaredHistogram->Fill(eventChiSquared);
    ChiSquaredHistogramPlanesHit[numPlanesHit-1]->Fill(eventChiSquared,1);
    // G4cout << "Chi^2 value: " << eventChiSquared << " Chi^2 prob: " << TMath::Prob(eventChiSquared, 10) << G4endl;

/////////////////////////////////////////////////////////////////////////////////////
  // Need to do a second stage GEANE backward tracing here from the first tracker plane hit to plane 0 and form the residuals from that. It needs to be done once outside the numPasses loop.

  G4ThreeVector xv3Backwards( startingXPos, startingYPos, startingZPos );
  G4ThreeVector pv3Backwards( -startingXMom, -startingYMom, -startingZMom ); // Reverse the momenta for backwards tracing.

  G4ErrorTrajErr errorStorageRegion( 5, 0 ); // Error matrix for backward tracing to storage region.

  G4ErrorFreeTrajState* backwardsFreeTrajState 
    = new G4ErrorFreeTrajState("e+", xv3Backwards, pv3Backwards, errorStorageRegion );

  double angle = 25.; // Use the same hardcoded value for the 0 plane location for now. Same for 460, which equates to a radius of a circle of about 2 or 3 meters out.
  G4double storagePlaneTargetX = (-460.*std::sin(angle*pi/180))*cm;
  G4double storagePlaneTargetZ = (460.-460.*std::cos(angle*pi/180))*cm;

  G4Normal3D storageNorm(1.,0.,0.);
            // surfNorm.setX(std::cos(angle*pi/180));
            // surfNorm.setZ(-std::sin(angle*pi/180));
  G4Point3D storagePos(storagePlaneTargetX,0.,storagePlaneTargetZ);
  
  theTarget = new G4ErrorPlaneSurfaceTarget(storageNorm, storagePos );
  g4edata->SetTarget( theTarget );
  theG4ErrorMode = G4ErrorMode_PropBackwards;

  g4emgr->Propagate( backwardsFreeTrajState, theTarget, theG4ErrorMode ); // Propagate() doesn't save the transport matrices properly, but that's not a problem since I don't care about them here, and just want the residuals at the end.

  G4ErrorMatrix matrixFill;
  G4ErrorSurfaceTrajState* backwardsSurfaceTrajState
            = new G4ErrorSurfaceTrajState(*backwardsFreeTrajState, vyV, wzU, matrixFill);

  trackParamPredicted[0][0] = backwardsSurfaceTrajState->GetParameters().GetInvP();

  G4double zeroStatePx = sqrt( 1./trackParamPredicted[0].at(0)*1./trackParamPredicted[0].at(0) - backwardsSurfaceTrajState->GetParameters().GetPV()*backwardsSurfaceTrajState->GetParameters().GetPV() - backwardsSurfaceTrajState->GetParameters().GetPW()*backwardsSurfaceTrajState->GetParameters().GetPW() );

  trackParamPredicted[1][0] = (backwardsSurfaceTrajState->GetParameters().GetPV()) / zeroStatePx; 
  trackParamPredicted[2][0] = (backwardsSurfaceTrajState->GetParameters().GetPW()) / zeroStatePx; 
  trackParamPredicted[3][0] = backwardsSurfaceTrajState->GetParameters().GetV();
  trackParamPredicted[4][0] = backwardsSurfaceTrajState->GetParameters().GetW();

  trackParamPredictedExtra[0][0] = backwardsSurfaceTrajState->GetPosition().x();
  trackParamPredictedExtra[1][0] = backwardsSurfaceTrajState->GetPosition().y(); // trackParamPredicted 3,4 and trackParamPredictedExtra 1,2 are the same with XYZ tracing for the vyV and wzU vectors equal to Y and Z respectively. Leaving the redundancy in case the GEANT4 source code is fixed later.
  trackParamPredictedExtra[2][0] = backwardsSurfaceTrajState->GetPosition().z();

  trackParamPredictedMomentumExtra[0].at(0) = zeroStatePx;
  trackParamPredictedMomentumExtra[1].at(0) = backwardsSurfaceTrajState->GetParameters().GetPV(); // Current Y and Z momenta.
  trackParamPredictedMomentumExtra[2].at(0) = backwardsSurfaceTrajState->GetParameters().GetPW();

  totalMomentumResidual = sqrt ((PlaneParameterTruth[0]["XMom"])*(PlaneParameterTruth[0]["XMom"]) 
                              + (PlaneParameterTruth[0]["YMom"])*(PlaneParameterTruth[0]["YMom"])
                              + (PlaneParameterTruth[0]["ZMom"])*(PlaneParameterTruth[0]["ZMom"])) - 1./trackParamPredicted[0][0];


  YMomentumResidual = PlaneParameterTruth[0]["YMom"] - -1.*backwardsSurfaceTrajState->GetMomentum().y(); // -1 in front since it's momentum is the negative of what I want.
  ZMomentumResidual = PlaneParameterTruth[0]["ZMom"] - -1.*backwardsSurfaceTrajState->GetMomentum().z();

  // I don't currently produce UV momentum residuals for any planes.
  // I want truth - best fit for the storage region since that's what I ultimately care about.

  UPositionResidual = PlaneParameterTruth[0]["UPos"] - (YZtoVUcoordinateTransformationMatrix(1,0)*trackParamPredicted[3].at(0) + YZtoVUcoordinateTransformationMatrix(1,1)*trackParamPredicted[4].at(0)); 
  VPositionResidual = PlaneParameterTruth[0]["VPos"] - (YZtoVUcoordinateTransformationMatrix(0,0)*trackParamPredicted[3].at(0) + YZtoVUcoordinateTransformationMatrix(0,1)*trackParamPredicted[4].at(0)); 

  YPositionResidual = PlaneParameterTruth[0]["YPos"] - trackParamPredictedExtra[1][0]; // These are fine for the Y and Z position residuals since they come after all of the UV track fitting.
  ZPositionResidual = PlaneParameterTruth[0]["ZPos"] - trackParamPredictedExtra[2][0];


/////////////////////////////////////////////////////////////////////////////////////

        PlaneTracebackResiduals[0]["P"] = (totalMomentumResidual);
        PlaneTracebackResiduals[0]["Py"] = (YMomentumResidual);
        PlaneTracebackResiduals[0]["Pz"] = (ZMomentumResidual);

        PlaneTracebackResiduals[0]["Y"] = (YPositionResidual);
        PlaneTracebackResiduals[0]["Z"] = (ZPositionResidual);
        PlaneTracebackResiduals[0]["U"] = (UPositionResidual);
        PlaneTracebackResiduals[0]["V"] = (VPositionResidual);

/////////////////////////////////////////////////////////////////////////////////////



  for (int planeNum = 1; planeNum < maxNumPlanes; ++planeNum) // Planes > 0 are treated separately from plane 0. 
  {

    if (PlaneParameterTruth[planeNum]["CopyNo"] == noHit) // I only want to pass residuals for planes that were hit. For planes that don't record U or V values, those residuals end up being ~noHit = 900000 values, which end up way outside the histogram ranges.
    { continue; }


    G4double planeTotalMomentum = sqrt ((PlaneParameterTruth[planeNum]["XMom"])*(PlaneParameterTruth[planeNum]["XMom"]) 
                                      + (PlaneParameterTruth[planeNum]["YMom"])*(PlaneParameterTruth[planeNum]["YMom"])
                                      + (PlaneParameterTruth[planeNum]["ZMom"])*(PlaneParameterTruth[planeNum]["ZMom"]));  

    PlaneTracebackResiduals[planeNum]["P"] = (planeTotalMomentum - 1./trackParamPredicted[0].at(planeNum));

    PlaneTracebackResiduals[planeNum]["Py"] = ((PlaneParameterTruth[planeNum]["YMom"]) - trackParamPredictedMomentumExtra[1].at(planeNum)); // Correct for XYZ tracing where Pv = Py and Pw = Pz, filled into trackParamPredictedMomentumExtra.
    PlaneTracebackResiduals[planeNum]["Pz"] = ((PlaneParameterTruth[planeNum]["ZMom"]) - trackParamPredictedMomentumExtra[2].at(planeNum));

    PlaneTracebackResiduals[planeNum]["Y"] = (PlanePositionMeasured[planeNum]["YPos"] - trackParamPredictedExtra[1].at(planeNum)); // Measured - best fit line residual. Really is truth - best fit since I don't measure and smear YZ.
    PlaneTracebackResiduals[planeNum]["Z"] = (PlanePositionMeasured[planeNum]["ZPos"] - trackParamPredictedExtra[2].at(planeNum));  

    PlaneTracebackResiduals[planeNum]["U"] = (PlanePositionMeasured[planeNum]["UPos"] - (YZtoVUcoordinateTransformationMatrix(1,0)*trackParamPredicted[3].at(planeNum) + YZtoVUcoordinateTransformationMatrix(1,1)*trackParamPredicted[4].at(planeNum)));
    PlaneTracebackResiduals[planeNum]["V"] = (PlanePositionMeasured[planeNum]["VPos"] - (YZtoVUcoordinateTransformationMatrix(0,0)*trackParamPredicted[3].at(planeNum) + YZtoVUcoordinateTransformationMatrix(0,1)*trackParamPredicted[4].at(planeNum)));
    // PlaneTracebackResiduals[planeNum]["U"] = (PlanePositionMeasured[planeNum]["UPos"] - (trackParamPredicted[4].at(planeNum)));
    // PlaneTracebackResiduals[planeNum]["V"] = (PlanePositionMeasured[planeNum]["VPos"] - (trackParamPredicted[3].at(planeNum)));


  }
  // With how I've changed things, this upper for loop and lower for loop can probably be combined, and the intermediate vector-map done away with, but I'm leaving it in for now for future functionality.


  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {

    TotalMomentumResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["P"],1);
    YMomentumResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Py"],1);
    ZMomentumResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Pz"],1);

    YPositionResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Y"],1);
    ZPositionResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["Z"],1);
    UPositionResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["U"],1);
    VPositionResiduals[planeNum]->Fill(PlaneTracebackResiduals[planeNum]["V"],1);

    PositionMeasuredHist[planeNum]->Fill(PlanePositionMeasured[planeNum]["ZPos"],PlanePositionMeasured[planeNum]["YPos"]);
    PositionTruthHist[planeNum]->Fill(PlaneParameterTruth[planeNum]["ZPos"],PlaneParameterTruth[planeNum]["YPos"]);

    // PositionMeasuredY[planeNum]->Fill(PlanePositionMeasured[planeNum]["YPos"]);

    // MeasurementInaccuracyY[planeNum]->Fill(PlanePositionMeasured[planeNum]["YPos"]-PlaneParameterTruth[planeNum]["YPos"]);
    // MeasurementInaccuracyZ[planeNum]->Fill(PlanePositionMeasured[planeNum]["ZPos"]-PlaneParameterTruth[planeNum]["ZPos"]);

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

// For the final program where it is unclear which measured parameters will go with which predicted parameters, I would pass in
// multiple different sets of paramMeasured to then calculate chi^2s and find the smallest one.
// Here I am just passing in what I know to be the smallest chi^2 parameters, for a true track that's offset in position and momenta..
void TrackCorrelation(std::vector<G4ErrorMatrix> myTransferMatrices, std::vector<G4ErrorMatrix> errorMatrices, std::vector <std::vector<G4double> > paramMeasured, std::vector <std::vector<G4double> > paramPredicted, int numPlanesHit, double deltaStartingTrack[5], double& eventChiSquared, int firstPlaneHit){
//Pass in array of transformed transfer/transport matrices, error matrices, track measured, track predicted, the number of planes hit (not currently used), references for track fixes and a track chiSquared, and the first plane hit for each correlation method called.

std::vector<Eigen::VectorXd> paramMeasuredEigen;
std::vector<Eigen::VectorXd> paramPredictedEigen; // Already in YZ
std::vector<Eigen::VectorXd> paramPredictedInVU;
std::vector<Eigen::VectorXd> residuals;

std::vector<double> chiSquaredSinglePlane;
chiSquaredSinglePlane.resize(maxNumPlanes, 0.); // Not sure whether to initially fill with 0 or some high number to simulate a noHit entry. I'm thinking 0.
double chiSquaredTotal = 0.;

std::vector<Eigen::MatrixXd> cov;
Eigen::MatrixXd covarianceTotal(5,5);

Eigen::VectorXd deltaPsiNought;
Eigen::VectorXd deltaPsiNoughtYZ;

std::vector<Eigen::VectorXd> trialTrajSinglePlane;
Eigen::VectorXd trialTrajectoryAllPlanes;

std::vector<Eigen::MatrixXd> transportMatrixBegToEnd;
std::vector<Eigen::MatrixXd> INVsigmaErrorMatrices;

std::vector<Eigen::MatrixXd> myTransferMatricesGeVcm;
std::vector<Eigen::MatrixXd> errorMatricesGeVcm;

// GeV cm are the better units for the matrix multiplication to work out, otherwise numerical errors can start ocurring with large transport matrices.

for(int ipl=0;ipl<maxNumPlanes;ipl++) {
      myTransferMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
      errorMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
   for(int i=0;i<5;i++) {
     for(int j=0; j<5;j++) {
          myTransferMatricesGeVcm[ipl](i,j) = myTransferMatrices[ipl][i][j]*1.; // GEVCM (Don't scale them to MeV mm.)
          errorMatricesGeVcm[ipl](i,j) = errorMatrices[ipl][i][j]*1.; // GEVCM
     }
   } 
 } 

 for (int ipl = 0; ipl < maxNumPlanes; ++ipl) // Use this loop to convert YZ transport and error matrices to VU.
 {
   myTransferMatricesGeVcm[ipl] =  JacobianToVU.inverse() * myTransferMatricesGeVcm[ipl] * JacobianToVU;
   errorMatricesGeVcm[ipl] =  JacobianToVU.inverse() * errorMatricesGeVcm[ipl] * JacobianToVU;
   // NOTE: The inverse is on the left side here because of how we generate the Jacobian between the two bases.
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

     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("Propagated error matrices in GeV cm %d \n",ipl);
       std::cout << errorMatricesGeVcm[ipl];
       printf("\n \n");
     }
#endif


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// 

// NOTE: Use this part for material error matrices with experimental error properly added. This is a separate loop because I was doing things a little differently before, and might again someday.

     for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
     {
        errorMatricesGeVcm[ipl](3,3) += 0.01*0.01; // GEVCM // Add in experimental part properly, .01 cm ^ 2.
        errorMatricesGeVcm[ipl](4,4) += 0.01*0.01; // GEVCM

#if MATRIXDEBUG
       printf("Propagated error matrices in VU GeV cm ERROR ADDED %d \n",ipl);
       std::cout << errorMatricesGeVcm[ipl] << G4endl;

       printf("Propagated error matrices in VU GeV cm ERROR ADDED inverted %d \n",ipl);
       std::cout << errorMatricesGeVcm[ipl].inverse() << G4endl;
#endif

      }


/////////////////////////////////////////////////////////////////////////////////////

      for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
      { 

            // For 2x2 case:
            // Fill INVsigma with zeroes, grab 5x5 errorMatricesGeVcm and invert it, then plug 2x2 part back into INVsigma.
            // The sigma error matrices needs to be filled with 0s in the rows and columns where there are no measurements.
            // In the case of a single U or single V measurement, I simply grab the single element from the inverted 5x5, with everything else zeros.
            // If a plane is not hit, it is filled with zeros. This matrix then goes on to multiply the other matrix pieces (possibly non-zero) in the tracking code, serving to remove any effects from non-hit planes that might otherwise introduce errors into the tracking.

            INVsigmaErrorMatrices.push_back(Eigen::MatrixXd::Zero(5,5));

            // Eigen::Matrix2d tempMatrix = (errorMatricesGeVcm[ipl].bottomRightCorner<2,2>()).inverse(); // This grabs the 2x2 part and then inverts it, which is incorrect I believe, although approximately the same since material errors are so low.
            Eigen::Matrix2d tempMatrix = (errorMatricesGeVcm[ipl].inverse()).bottomRightCorner<2,2>(); // This is inverting a 5x5 which can probably be replaced with a solve method.

            // INVsigmaErrorMatrices[ipl].bottomRightCorner<2,2>() = tempMatrix;           

#if MATRIXDEBUG
 G4cout << "Plane " << ipl;
#endif

            if (paramMeasured[3][ipl] == noHit && paramMeasured[4][ipl] == noHit)
            { 
 #if MATRIXDEBUG            
  G4cout << " Hit neither. " << G4endl;
  #endif
              continue; // Leave the sigma error matrix as zeros.
            }

            else if (paramMeasured[3][ipl] == noHit && paramMeasured[4][ipl] != noHit)
            {
#if MATRIXDEBUG
 G4cout << " Hit U. " << G4endl;
#endif
              INVsigmaErrorMatrices[ipl](4,4) = tempMatrix(1,1); // 11 element for bottom right 2x2 matrix element.
              // INVsigmaErrorMatrices[ipl](4,4) = 10000.; // This for ignoring the material error and only using the experimental error.
            }

            else if (paramMeasured[4][ipl] == noHit && paramMeasured[3][ipl] != noHit)
            {
#if MATRIXDEBUG
  G4cout << " Hit V. " << G4endl;
#endif
              INVsigmaErrorMatrices[ipl](3,3) = tempMatrix(0,0); // 11 element for bottom right 2x2 matrix element.
              // INVsigmaErrorMatrices[ipl](3,3) = 10000.; // This for ignoring the material error and only using the experimental error.
            }

            else if (paramMeasured[3][ipl] != noHit && paramMeasured[4][ipl] != noHit)
            { 
#if MATRIXDEBUG
 G4cout << " Hit both U and V. " << G4endl;
#endif
                INVsigmaErrorMatrices[ipl].bottomRightCorner<2,2>() = tempMatrix; // This selects the bottom right 2x2 block.
                
                // INVsigmaErrorMatrices[ipl](4,4) = 10000.;// These for ignoring the material error and only using the experimental error.
                // INVsigmaErrorMatrices[ipl](3,3) = 10000.;
            }

      }



/////////////////////////////////////////////////////////////////////////////////////
#if MATRIXDEBUG
     // Debug printing
     printf("\n");
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       printf("Full inverted error matrices %d \n",ipl);
       std::cout << INVsigmaErrorMatrices[ipl];
       printf("\n \n");
     }
#endif
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////



      // Multiply transport matrices together. Be careful with Transport matrix = 0 matrix for the starting plane.

      for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
      { 
        transportMatrixBegToEnd.push_back(Eigen::MatrixXd::Zero(5,5)); // Initialize matrices.
      }

      transportMatrixBegToEnd[firstPlaneHit] = myTransferMatricesGeVcm[firstPlaneHit]; // Fill the first plane hit matrix, ignore the previous plane matrices since they are just 0's and won't be included in the sum later on.

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
        cov[ipl] = (transportMatrixBegToEnd[ipl].transpose())*INVsigmaErrorMatrices[ipl]*transportMatrixBegToEnd[ipl];
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

    //       Convert parameters arrays into Eigen objects and YZ measured to UV measured:   
    //       ====================
            for(int ipl=0;ipl<maxNumPlanes;ipl++) {
              paramMeasuredEigen.push_back(Eigen::VectorXd::Zero(5));
              paramPredictedEigen.push_back(Eigen::VectorXd::Zero(5));
              paramPredictedInVU.push_back(Eigen::VectorXd::Zero(5));

               for(int i=0;i<5;i++) {
                  paramMeasuredEigen[ipl](i) = paramMeasured[i][ipl];
                  paramPredictedEigen[ipl](i) = paramPredicted[i][ipl];

               } 
             }

            for(int ipl=0;ipl<maxNumPlanes;ipl++) {
                    paramPredictedInVU[ipl] = YZtoVUcoordinateTransformFiveByFive * paramPredictedEigen[ipl]; 
             } 

/////////////////////////////////////////////////////////////////////////////////////
    //
    //       Now we can calculate the next guess:
    //       ====================================
    //
    //       Calculate the residuals at every step: - Far right side of equation 23, 24, or 26 in the geane manual paper.
    //       ======================================
    //             

            for(int ipl=0;ipl<maxNumPlanes;ipl++) {
              residuals.push_back(Eigen::VectorXd::Zero(5));
              residuals[ipl] = paramMeasuredEigen[ipl]-paramPredictedInVU[ipl];
            }

            // GEVCM Convert residuals from MeV mm to GeV cm to multiply against GeV cm transport and error matrices.
            for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
            {
              residuals[ipl](0) = residuals[ipl](0) * 1.0e3;
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
               std::cout << "Eigen objects: " << G4endl << paramMeasuredEigen[ipl] << std::endl;
               printf("\n \n");
             } 

             printf("\n");
             for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
               printf("params predicted Plane %d \n",ipl);
               std::cout << "Eigen objects: " << G4endl << paramPredictedEigen[ipl] << G4endl;
               std::cout << " And from YZ to VU: " << G4endl << paramPredictedInVU[ipl]; // Also Eigen
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
               trialTrajSinglePlane[ipl] = (transportMatrixBegToEnd[ipl].transpose())*INVsigmaErrorMatrices[ipl]*residuals[ipl];
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

            // Form the chi^2 for a single plane then add them all up to get a total chi^2 for the potential track - equation 23 in the GEANE manual.
            for(int ipl=0;ipl<maxNumPlanes;ipl++) {
               chiSquaredSinglePlane[ipl] = (residuals[ipl].transpose())*INVsigmaErrorMatrices[ipl]*residuals[ipl];
             } 

             for(int ipl = 1; ipl<maxNumPlanes; ipl++) { // Sum starts at 1.
               chiSquaredTotal = chiSquaredTotal + chiSquaredSinglePlane[ipl];
             }

             eventChiSquared = chiSquaredTotal;

             // G4cout << "Chi^2 for the track is: " << chiSquaredTotal << G4endl;


              // I can multiply the Jacobian against deltaPsiNought here to immediately convert back to YZ for easire use in the main code. Careful with inverses.

              deltaPsiNoughtYZ = JacobianToVU * deltaPsiNought;

              // GEVCM
              deltaStartingTrack[0] = deltaPsiNoughtYZ(0) * 1.0e-3; // Convert back to MeV mm units for the main code to use.
              deltaStartingTrack[1] = deltaPsiNoughtYZ(1) * 1.0e0;
              deltaStartingTrack[2] = deltaPsiNoughtYZ(2) * 1.0e0;
              deltaStartingTrack[3] = deltaPsiNoughtYZ(3) * 1.0e1;
              deltaStartingTrack[4] = deltaPsiNoughtYZ(4) * 1.0e1;

 return;

}



void ResidualPlots(){

  TCanvas* myCanvas = new TCanvas("myCanvas","Canvas",200,10,1200,800);
  gStyle->SetOptStat(111111);

  double YRMSarray[maxNumPlanes];
  double ZRMSarray[maxNumPlanes];

  double URMSarray[maxNumPlanes];
  double VRMSarray[maxNumPlanes];

  double pointNo[maxNumPlanes];
  for (int i = 0; i < maxNumPlanes; ++i)
  {
    pointNo[i] = i;
  }

  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {


  YRMSarray[planeNum] = YPositionResiduals[planeNum]->GetRMS();
  ZRMSarray[planeNum] = ZPositionResiduals[planeNum]->GetRMS();

  URMSarray[planeNum] = UPositionResiduals[planeNum]->GetRMS();
  VRMSarray[planeNum] = VPositionResiduals[planeNum]->GetRMS();


  // MeasurementInaccuracyY[planeNum]->Fit("gaus");
  // MeasurementInaccuracyZ[planeNum]->Fit("gaus");

  // PositionMeasuredY[planeNum]->Fit("pol0","","",-.8*myGhostWorld->GetTruthPlaneHalfY(),.8*myGhostWorld->GetTruthPlaneHalfY());



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
  UPositionResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/UPositionResiduals.png").c_str());
  VPositionResiduals[planeNum]->Draw();
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/VPositionResiduals.png").c_str());


  PositionMeasuredHist[planeNum]->Draw("colz");
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/PositionMeasured.png").c_str());
  PositionTruthHist[planeNum]->Draw("colz");
  myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/PositionTruth.png").c_str());

  // PositionMeasuredY[planeNum]->Draw();
  // myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/PositionMeasuredY.png").c_str());

  // MeasurementInaccuracyY[planeNum]->Draw();
  // myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/MeasurementInaccuracyY.png").c_str());
  // MeasurementInaccuracyZ[planeNum]->Draw();
  // myCanvas->SaveAs(("../Plots/Plane"+std::to_string(planeNum)+"/MeasurementInaccuracyZ.png").c_str());

  ChiSquaredHistogramPlanesHit[planeNum]->Draw();
  TF1* chi2pdf = new TF1("chi2pdf","[2]*ROOT::Math::chisquared_pdf(x,[0],[1])",0,60);
  chi2pdf->SetParameters(planeNum-5, 0., ChiSquaredHistogramPlanesHit[planeNum]->Integral("WIDTH")); 

  // G4cout << "Num entries : " << ChiSquaredHistogramPlanesHit[planeNum]->GetEntries() << G4endl;
  // G4cout << "Integral : " << ChiSquaredHistogramPlanesHit[planeNum]->Integral() << G4endl;
  // G4cout << "Integral with width : " << ChiSquaredHistogramPlanesHit[planeNum]->Integral("WIDTH") << G4endl;
  // G4cout << "Integral of func: " << chi2pdf->Integral(0,60) << G4endl;

  chi2pdf->Draw("SAME");
  myCanvas->SaveAs(("../Plots/ChiSquaredsPerPlane/ChiSquaredsPlanesHit"+std::to_string(planeNum)+".png").c_str());

  // delete ChiSquaredHistogramPlanesHit[planeNum];

  delete TotalMomentumResiduals[planeNum];
  delete YMomentumResiduals[planeNum];
  delete ZMomentumResiduals[planeNum];

  delete YPositionResiduals[planeNum];
  delete ZPositionResiduals[planeNum];
  delete UPositionResiduals[planeNum];
  delete VPositionResiduals[planeNum];  

  delete PositionMeasuredHist[planeNum];
  delete PositionTruthHist[planeNum];

  // delete PositionMeasuredY[planeNum];

  // delete MeasurementInaccuracyY[planeNum];
  // delete MeasurementInaccuracyZ[planeNum];

  }


  ChiSquaredHistogram->Draw();
  // ChiSquaredProbHistogram->SetMarkerColor(5);
  // ChiSquaredProbHistogram->Draw("SAME");
  myCanvas->SaveAs("../Plots/ChiSquareds.png");

  delete ChiSquaredHistogram;



/////////////////////////////////////////////////////////////////////////////////////

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
  RMSPerPlaneZ->Draw("PSAME"); // SAME to draw onto the same canvas without overwriting the previous canvas.

  TLegend* legYZ = new TLegend(0.663,0.70,0.891,0.887);
  legYZ->AddEntry(RMSPerPlaneY,"Y position","p");
  legYZ->AddEntry(RMSPerPlaneZ,"Z position","p");
  legYZ->SetFillStyle(0);
  legYZ->SetBorderSize(0);
  legYZ->Draw();

  myCanvas->SaveAs("../Plots/YZPlaneRMS.png");

  delete RMSPerPlaneY;
  delete RMSPerPlaneZ;
  delete legYZ;

/////////////////////////////////////////////////////////////////////////////////////

  TGraph* materialErrorGraph = new TGraph(maxNumPlanes, pointNo, materialErrorPerPlane);
  materialErrorGraph->SetTitle("Material Error Per Plane; Plane Number; Error (um)");
  // materialErrorGraph->GetXaxis()->SetRangeUser(0,8);
  materialErrorGraph->GetYaxis()->SetRangeUser(0.,200);
  materialErrorGraph->SetMarkerStyle(20);
  materialErrorGraph->SetMarkerColor(1);
  materialErrorGraph->Draw("AP"); //A for "all" or something, nothing plots without it. L for connected lines between points, and P for dot style points.

  myCanvas->SaveAs("../Plots/MaterialErrorPerPlane.png");

  double planeXPositions[maxNumPlanes];
  planeXPositions[0] = 0.;

  for (int planeNum = 1; planeNum < maxNumPlanes; ++planeNum)
  {
          G4double truthPlaneTargetX=0.;
          // Have to do some arithmatic to determine how to change truthPlaneTargetX with a single for-loop for 32 planes.

          truthPlaneTargetX += int((planeNum-1)/4.) * myPhysicalWorld->GetmoduleSeparation(); // Add module separation.

          bool layerSep = FALSE;
          if (planeNum % 2 == 0)
          {
            truthPlaneTargetX += myPhysicalWorld->GetlayerSeparation(); // Add layer separation.
            layerSep = TRUE;
          }

          double intpart;
          double fractpartone = std::modf(planeNum/4. , &intpart);
          double fractparttwo = std::modf((planeNum+1)/4. , &intpart);
          bool UVSep = FALSE;
          if (fractpartone == 0 || fractparttwo == 0)
          {
            truthPlaneTargetX += myPhysicalWorld->GetUVSeparation(); // Add UV separation.
            UVSep = TRUE;
          }

          planeXPositions[planeNum] = truthPlaneTargetX;
  }

  TGraph* materialErrorGraphDistance = new TGraph(maxNumPlanes, planeXPositions, materialErrorPerPlane);
  materialErrorGraphDistance->SetTitle("Material Error Per Distance; Distance (cm); Error (um)");
  // materialErrorGraphDistance->GetXaxis()->SetRangeUser(0,8);
  materialErrorGraphDistance->GetYaxis()->SetRangeUser(0.,200);
  materialErrorGraphDistance->SetMarkerStyle(20);
  materialErrorGraphDistance->SetMarkerColor(1);
  materialErrorGraphDistance->Draw("AP"); //A for "all" or something, nothing plots without it. L for connected lines between points, and P for dot style points.

  myCanvas->SaveAs("../Plots/MaterialErrorDistance.png");


///////////////////////////////////////////////////////////////////////////////////// 

  TGraph* RMSPerPlaneU = new TGraph(maxNumPlanes, pointNo, URMSarray);
  RMSPerPlaneU->SetTitle("RMS; Plane Number; RMS (mm)");
  // RMSPerPlaneU->GetXaxis()->SetRangeUser(0,8);
  RMSPerPlaneU->GetYaxis()->SetRangeUser(0.,.15);
  RMSPerPlaneU->SetMarkerStyle(20);
  RMSPerPlaneU->SetMarkerColor(1);
  RMSPerPlaneU->Draw("AP"); //A for "all" or something, nothing plots without it. L for connected lines between points, and P for dot style points.

  TGraph* RMSPerPlaneV = new TGraph(maxNumPlanes, pointNo, VRMSarray);
  RMSPerPlaneV->SetMarkerStyle(20);
  RMSPerPlaneV->SetMarkerColor(2);
  RMSPerPlaneV->Draw("PSAME");

  TLegend* legUV = new TLegend(0.663,0.70,0.891,0.887);
  legUV->AddEntry(RMSPerPlaneU,"U position","p");
  legUV->AddEntry(RMSPerPlaneV,"V position","p");
  legUV->SetFillStyle(0);
  legUV->SetBorderSize(0);
  legUV->Draw();

  myCanvas->SaveAs("../Plots/UVPlaneRMS.png");

  delete RMSPerPlaneU;
  delete RMSPerPlaneV;
  delete legUV;

/////////////////////////////////////////////////////////////////////////////////////

  delete myCanvas;

}




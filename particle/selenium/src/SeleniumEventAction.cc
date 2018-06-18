// Selenium Event action code

#include "SeleniumEventAction.hh"
#include "SeleniumRunAction.hh"
#include "SeleniumAnalysis.hh"


#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

// Selenium Event Action constructor
SeleniumEventAction::SeleniumEventAction()
:G4UserEventAction(),
 fEnergyDepositSe(0.),
 fEnergyDepositAu(0.)
{

}

// Destructor
SeleniumEventAction::~SeleniumEventAction()
{

}

void SeleniumEventAction::BeginOfEventAction(const G4Event* event)
{
	// Initialize energy deposited in each material
	fEnergyDepositAu = 0.;
	fEnergyDepositSe = 0.;

}

void SeleniumEventAction::EndOfEventAction(const G4Event* event)
{
	// Get analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// Add deposited energy to histograms
	analysisManager->FillH1(0, fEnergyDepositSe);
	analysisManager->FillH1(1, fEnergyDepositAu);

}
// Selenium detector Run Action

#include "SeleniumRunAction.hh"
#include "SeleniumAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "G4SystemofUnits.hh"
#include "Randomize.hh"

#include "SimulationSettings.hh"

// Run Action Constructor
SeleniumRunAction::SeleniumRunAction(G4String name) : G4UserRunAction()
{
	// set printing event number per each event
	G4RunManager::GetRunManager()->SetPrintProgress(1);

	// Create analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetFirstHistoId(1);

	// Create the histograms
	analysisManager->CreateH1("EneDepSe", "Energy Deposited in Selenium", 100, 0., 140.*keV, "keV");
	analysisManager->CreateH1("EneDepAu", "Energy Deposited in Gold Electrode", 100, 0., 140.*keV, "keV");
	filename = OUTPUTDIR + name + "_keV";
	analysisManager->SetFileName(filename);

}

// Run Action Destructor
SeleniumRunAction::~SeleniumRunAction()
{
	delete G4AnalysisManager::Instance();
}

// BeginOfRunAction method. Executed before the first run
void SeleniumRunAction::BeginOfRunAction(const G4Run*)
{
	// Get the analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// Open file for output of data
	analysisManager->OpenFile();
}

// Method run at the end of the run
void SeleniumRunAction::EndOfRunAction(const G4Run*)
{
	// Save the histograms
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

}
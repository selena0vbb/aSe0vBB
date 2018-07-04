// Selenium Event action code

#include "SeleniumEventAction.hh"
#include "SeleniumRunAction.hh"
#include "SeleniumAnalysis.hh"
#include "SimulationSettings.hh"
#include "SeleniumHit.hh"


#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

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
	if ( fEnergyDepositSe > 0. )
	{
		analysisManager->FillH1(1, fEnergyDepositSe);	
	}
	if ( fEnergyDepositAu > 0. )
	{
		analysisManager->FillH1(2, fEnergyDepositAu);	
	}

	// Get the number of hits in the selenium. For each hit add 
	// it to the Energy vs position histogram
	G4HCofThisEvent* hce = event->GetHCofThisEvent();

	// Get the hits collection ID number
	G4int seleniumHCID = G4SDManager::GetSDMpointer()->GetCollectionID(seleniumHitsCollectionName);

	if (seleniumHCID >= 0)
	{
		SeleniumHitsCollection* hc = dynamic_cast<SeleniumHitsCollection*>(hce->GetHC(seleniumHCID));
		SeleniumHit* hit;
		G4int eventID = event->GetEventID();
		// Get number of hits. Iterate over all and add to Ntuple
		G4int numberOfHits = hc->GetSize();
		G4cout << numberOfHits << G4endl;
		for( int i = 0; i < numberOfHits; i++)
		{
			hit = (*hc)[i];
			G4cout << hit->GetCreatorProcessName() << G4endl;
			analysisManager->FillNtupleIColumn(0, eventID);
			analysisManager->FillNtupleIColumn(1, hit->GetTrackID());
			analysisManager->FillNtupleIColumn(2, hit->GetParentID());
			analysisManager->FillNtupleDColumn(3, hit->GetPosition().getX());
			analysisManager->FillNtupleDColumn(4, hit->GetPosition().getY());
			analysisManager->FillNtupleDColumn(5, hit->GetPosition().getZ());
			analysisManager->FillNtupleDColumn(6, hit->GetEdep());
			analysisManager->FillNtupleSColumn(7, hit->GetParticleDefinition()->GetParticleName());
			analysisManager->FillNtupleSColumn(8, hit->GetCreatorProcessName());
			analysisManager->AddNtupleRow();

		}

	}
	
	

}
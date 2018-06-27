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

		// Define the parameters we put in the histogram
		G4double zPosition = 0.;
		G4double edep = 0.;

		// Get number of hits. Iterate over all and add to histogram
		G4int numberOfHits = hc->GetSize();
		for( int i = 0; i < numberOfHits; i++)
		{
			SeleniumHit* hit = (*hc)[i];
			// Check to see if electron caused the hit
			if ( hit->GetParticleDefinition()->GetParticleName() == "e-")
			{
				zPosition = hit->GetPosition().getZ()*um;
				edep = hit->GetEdep()*eV;
				analysisManager->FillH1(3, zPosition, edep);
			}

		}

	}
	
	

}
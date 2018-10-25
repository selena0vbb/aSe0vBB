// Selenium Stepping Action code

#include "SeleniumSteppingAction.hh"
#include "SeleniumEventAction.hh"
#include "SeleniumDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"


// Constructor
SeleniumSteppingAction::SeleniumSteppingAction(
									SeleniumDetectorConstruction* detector,
									SeleniumEventAction* event)
	: G4UserSteppingAction(),
	  fDetectorConstruction(detector),
	  fEventAction(event)
{

}

// Destructor
SeleniumSteppingAction::~SeleniumSteppingAction()
{

}

// Stepping Action. Gets the energy and adds it to the current event
void SeleniumSteppingAction::UserSteppingAction(const G4Step* step)
{
	// Get the volume
	G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

	// Energy deposit
	G4double eDepStep = step->GetTotalEnergyDeposit();

	if( volume == fDetectorConstruction->GetSeleniumPV())
	{
		fEventAction->addEnergyDepSe(eDepStep);
	} else if( volume == fDetectorConstruction->GetGlassPV())
	{
		fEventAction->addEnergyDepSiO2(eDepStep);
	}

}

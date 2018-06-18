// Selenium Action Initialization code

#include "SeleniumActionInitialization.hh"
#include "SeleniumPrimaryGeneratorAction.hh"
#include "SeleniumRunAction.hh"
#include "SeleniumEventAction.hh"
#include "SeleniumSteppingAction.hh"
#include "SeleniumDetectorConstruction.hh"

// Constructor. Calls the parent class constructor and passes in the argument to the fDetectorConstruction member variable
SeleniumActionInitialization::SeleniumActionInitialization(G4String filePrefix,
	SeleniumDetectorConstruction* detector)
	: G4VUserActionInitialization(),
	  fDetectorConstruction(detector)
{
	name=filePrefix;
}

// Destructor
SeleniumActionInitialization::~SeleniumActionInitialization()
{

}

// Run on master. Executed on all threads
void SeleniumActionInitialization::BuildForMaster() const 
{
	SetUserAction(new SeleniumRunAction(name));
}

//
void SeleniumActionInitialization::Build() const
{
	// Creating and assigning the various action needed for the simulation
	// Run Action
	SetUserAction(new SeleniumRunAction(name));

	// Primary generator action
	SeleniumPrimaryGeneratorAction* spga = new SeleniumPrimaryGeneratorAction();
	SetUserAction(spga);

	// Event action
	SeleniumEventAction* event = new SeleniumEventAction;
	SetUserAction(event);

	// Stepping action
	SeleniumSteppingAction* stepping = new SeleniumSteppingAction(fDetectorConstruction, event);
	SetUserAction(stepping);

}


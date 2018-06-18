// Selenium Stepping Action Header

#ifndef SeleniumSteppingAction_H
#define SeleniumSteppingAction_H

#include "G4UserSteppingAction.hh"

class SeleniumDetectorConstruction;
class SeleniumEventAction;

// Define Stepping Action class. Collects energy deposited in different materials
class SeleniumSteppingAction : public G4UserSteppingAction
{

	public:
		SeleniumSteppingAction(SeleniumDetectorConstruction* detector, 
									SeleniumEventAction* event);
		virtual ~SeleniumSteppingAction();

		virtual void UserSteppingAction(const G4Step* step);

	private:
		SeleniumDetectorConstruction* fDetectorConstruction;
		SeleniumEventAction* fEventAction;
};

#endif
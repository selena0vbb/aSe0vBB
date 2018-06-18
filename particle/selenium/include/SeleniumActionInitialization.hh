// Selenium Action Initializer header

#ifndef SeleniumActionInitialization_H
#define SeleniumActionInitialization_H

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

// Define class for detector construction
class SeleniumDetectorConstruction;

// Define action initialization class
class SeleniumActionInitialization: public G4VUserActionInitialization
{

	public:
		SeleniumActionInitialization(G4String filePrefix, SeleniumDetectorConstruction* detector);
		virtual ~SeleniumActionInitialization();

		virtual void BuildForMaster() const;
		virtual void Build() const;
		G4String name;
		
	private:
		SeleniumDetectorConstruction* fDetectorConstruction;
};

#endif


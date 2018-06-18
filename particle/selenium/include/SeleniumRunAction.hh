// Selenium Detector Run Action

#ifndef SeleniumRunAction_H
#define SeleniumRunAction_H

#include "G4UserRunAction.hh"
#include "globals.hh"

// Define the G4Run class
class G4Run;

// Define the Run Action class
class SeleniumRunAction : public G4UserRunAction
{

	public:
		SeleniumRunAction(G4String name);
		virtual ~SeleniumRunAction();

		virtual void BeginOfRunAction(const G4Run*);
		virtual void EndOfRunAction(const G4Run*);
		
		// Defines a filename for saving output data
		G4String filename;

};

#endif
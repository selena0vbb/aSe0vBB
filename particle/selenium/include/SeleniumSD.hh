// Selenium sensitive detector header file

#ifndef SeleniumSD_H
#define SeleniumSD_H

#include "G4VSensitiveDetector.hh"
#include "SeleniumHit.hh"

#include <vector>

// Define the necessary classes for SeleniumSD class
class G4Step;
class G4HCofThisEvent;

// Class definition
class SeleniumSD : public G4VSensitiveDetector
{
	public:
		// Constructor and destructor
		SeleniumSD(const G4String& name, const G4String& hitsCollectionName);
		virtual ~SeleniumSD();

	public:
		// overloads the process hits base class
		virtual void Initialize(G4HCofThisEvent* hitCollection);
		virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory* history);
		virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

private:
	SeleniumHitsCollection* fHitsCollection;
	G4int fHitsCollectionID;

};


#endif
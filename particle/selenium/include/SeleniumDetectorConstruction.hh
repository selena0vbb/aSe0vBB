// Header file for Selenium Detector

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

// Detector Construction class to define materials and geometry of detector

class SeleniumDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		SeleniumDetectorConstruction();
		virtual ~SeleniumDetectorConstruction();


	public:
		virtual G4VPhysicalVolume* Construct();

		// Methods to get the physical volumes of the different parts of the detector
		const G4VPhysicalVolume* GetSeleniumPV() const;
		const G4VPhysicalVolume* GetElectrodePV() const;


	private:
		void defineMaterials();

		G4VPhysicalVolume* fSeleniumPV; // Selenium physical volume
		G4VPhysicalVolume* fElectrodePV; // Gold electrode physical volume
	
};


inline const G4VPhysicalVolume* SeleniumDetectorConstruction::GetSeleniumPV() const 
{
	return fSeleniumPV;
}

inline const G4VPhysicalVolume* SeleniumDetectorConstruction::GetElectrodePV() const 
{
	return fElectrodePV;
}
// SeleniumDetectorConstruction Implementation

#include "SeleniumDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>


// Define Constructor. Defines the materials used in the detector
SeleniumDetectorConstruction::SeleniumDetectorConstruction() : G4VUserDetectorConstruction()
{
	defineMaterials();
}

// Defines Destructor
SeleniumDetectorConstruction::~SeleniumDetectorConstruction()
{

}

// Define Materials. Creates Selenium for body of detector and gold for electrode
void SeleniumDetectorConstruction::defineMaterials()
{
	G4NistManager* matManager = G4NistManager::Instance();

	// Create the necessary elements. Selenium for body of detector, gold for top electrode
	G4Element* Se = matManager->FindOrBuildElement("G4_Se");
	G4Element* Au = matManager->FindOrBuildElement("G4_Au");

}


// Defines the geometry of the detecotr. Retuns a G4VPhyscalVolume class

G4VPhysicalVolume* SeleniumDetectorConstruction::Construct()
{

	// Define the geometric parameters that we know of
	G4double seleniumXY = 10 * mm;
	G4double seleniumDepth = 0.2 * mm;
	G4double goldElectrodeArea = 2 * mm2;
	G4double goldElectrodeDepth = 10 * nm;

	G4double goldElectrodeRadius = std::sqrt(goldElectrodeArea / pi) * mm;

	G4bool checkOverlaps = true;

	// Get materials for the detector
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* defaultMat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* Se = nist->FindOrBuildMaterial("G4_Se");
	G4Material* Au = nist->FindOrBuildMaterial("G4_Au");

	// Define the world size and create the world
	G4double worldSizeXY = 3*seleniumXY;
	G4double worldSizeZ = 10 * mm;

	G4Box* solidWorld = new G4Box(	"World",
									worldSizeXY/2, // world x dimension
									worldSizeXY/2, // world y dimension
									worldSizeZ/2); // world z dimension
	G4LogicalVolume* logicalWorld = 
		new G4LogicalVolume(solidWorld, 	// solid object associated with logical volume
							defaultMat, 	// material of the world
							"World");		//name

	G4VPhysicalVolume* physWorld = 
		new G4PVPlacement(	0, 					// rotation
							G4ThreeVector(), 	// position at 0,0,0
							logicalWorld, 		// logical volume
							"World",			// name
							0,
							false,
							0);


	// Define Selenium block for detector
	G4Box* seleniumS = new G4Box("Selenium",
								seleniumXY/2,
								seleniumXY/2,
								seleniumDepth/2);

	G4LogicalVolume* seleniumLV = 
		new G4LogicalVolume(seleniumS,
							Se,
							"Selenium");

	fSeleniumPV = 
		new G4PVPlacement(	0,
							G4ThreeVector(),
							seleniumLV,
							"Selenium",
							logicalWorld,
							false,
							0);

	// Defining and placing gold electrode on top of selenium
	// Use G4Tub object to define a cylinder
	G4Tubs* electrodeS = new G4Tubs(	"Electrode",
										0., 					// inner radius
										goldElectrodeRadius, 	// outer radius
										goldElectrodeDepth/2, 	// cylinder length
										0.,
										twopi); 				// Starting and ending angles (2pi to define whole cylinder)

	G4LogicalVolume* electrodeLV = new G4LogicalVolume(electrodeS, Au, "Electrode");

	G4ThreeVector electrodeOffset = G4ThreeVector(0, 0, seleniumDepth/2 + goldElectrodeDepth/2);
	fElectrodePV = 
		new G4PVPlacement( 	0,
							G4ThreeVector(0, 0, seleniumDepth/2 + goldElectrodeDepth/2),
							electrodeLV,
							"Electrode",
							logicalWorld,
							false,
							0);


	// Set up visualization attributes
		
	// returns the physical world
	return physWorld;


}

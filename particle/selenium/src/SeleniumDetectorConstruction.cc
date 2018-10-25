// SeleniumDetectorConstruction Implementation

#include "SeleniumDetectorConstruction.hh"
#include "SeleniumSD.hh"

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
#include "SimulationSettings.hh"



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
	G4Material* SiO2 = matManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	G4Element* Cu = matManager->FindOrBuildElement("G4_Cu");

}


// Defines the geometry of the detecotr. Retuns a G4VPhyscalVolume class

G4VPhysicalVolume* SeleniumDetectorConstruction::Construct()
{

	// Compute important geometric properties
	G4double goldElectrodeRadius = std::sqrt(goldElectrodeArea / pi) * mm;

	G4bool checkOverlaps = true;

	// Get materials for the detector
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* defaultMat = nist->FindOrBuildMaterial("G4_AIR");
	G4Element* Se = nist->FindOrBuildElement("G4_Se");
	G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	G4Element* Cu = nist->FindOrBuildElement("G4_Cu");

	// Define the world size and create the world
	G4double worldSizeXY = 3*copperXY;
	G4double worldSizeZ = 3*(seleniumDepth + glassDepth + copperDepth + airSpacingZ);

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



	// Define large copper top plate
	G4Box* copperS = new G4Box("CopperBlocking",
	            copperXY/2,
	            copperXY/2,
	            copperDepth/2);

	G4LogicalVolume* copperLV =
		new G4LogicalVolume(copperS,
		          Cu,
		          "CopperBlocking");

	G4VPhysicalVolume* copperPV =
		new G4PVPlacement(0,
		          G4ThreeVector(0, 0, -(seleniumDepth + copperDepth)/2 -glassDepth - airSpacingZ),
		          copperLV,
		          "CopperBlocking",
		          logicalWorld,
		          false,
		          0);


	// Create copper electrodes for pixel electrodes
	// Top Plate
	G4Box* elecTopS = new G4Box("TopElectrode",
	            seleniumXY/2,
	            seleniumXY/2,
	            electrodeDepth/2);

	G4LogicalVolume* elecTopLV =
		new G4LogicalVolume(elecTopS,
		          Cu,
		          "TopElectrode");

	G4VPhysicalVolume* elecTopPV =
		new G4PVPlacement(0,
		          G4ThreeVector(0, 0, (seleniumDepth + electrodeDepth)/2),
		          elecTopLV,
		          "TopElectrode",
		          logicalWorld,
		          false, 0);


	// Define Pixel electrode
	G4Solid* pixelS = new G4Tub("Pixel",
	            0,
	            pixelRadius,
	            electrodeDepth/2,
	            0, twopi);

	G4LogicalVolume* pixelLV =
		new G4LogicalVolume(pixelS,
		          Cu,
		          "Pixel");

	G4VPhyscalVolume* pixelPV =
		new G4PVPlacement(0,
		          G4ThreeVector(0, 0, -(seleniumDepth + electrodeDepth)/2),
		          pixelLV,
		          "Pixel",
		          logicalWorld,
		          false, 0);

	// Create Guard Ring
	G4Solid* guardS = new G4Tub("GuardRing",
	            guardInnerRadius,
	            guardOuterRadius,
	            electrodeDepth/2,
	            0, twopi);

	G4LogicalVolume* guardLV =
		new G4LogicalVolume(guardS,
		          Cu,
		          "GuardRing");

	G4VPhyscalVolume* guardPV =
		new G4PVPlacement(0,
		          G4ThreeVector(0, 0, -(seleniumDepth + electrodeDepth)/2),
		          guardLV,
		          "GuardRing",
		          logicalWorld,
		          false, 0);

	// Define glass substrate, need to use volume subtraction to not have it overlap with the

	// Define SiO2 Glass substrate
	G4Box* glassS = new G4Box("Glass",
	            seleniumXY/2,
	            seleniumXY/2,
	            glassDepth/2);

	// Remove the pixel
	G4Solid* glassMinusPixelS =
		new G4SubtractionSolid("Glass-Pixel",
		          glassS,
		          pixelS,
		          G4ThreeVector(0, 0, (glassDepth - electrodeDepth)/2));

	// Remove the guard ring
	G4Solid* glassMinusElecS =
		new G4SubtractionSolid("Glass-Elec",
		          glassMinusPixelS,
		          guardS,
		          G4ThreeVector(0, 0, (glassDepth - electrodeDepth)/2));

	G4LogicalVolume* glassLV =
		new G4LogicalVolume(glassMinusElecS,
		          SiO2,
		          "Glass");

	fGlassPV = new G4PVPlacement(0,
	            G4ThreeVector(0,0,-(seleniumDepth + glassDepth)/2),
	            glassLV,
	            "Glass",
	            logicalWorld,
	            false,
	            0);

	// Setting up sensitive detector for the selenium body. Only uses if SIMPLE is false
	if(!SIMPLE)
	{
		// Create a new selenium sensitive detector
		SeleniumSD* detector = new SeleniumSD("SeleniumSD", seleniumHitsCollectionName);

		// Get pointer to SD manager
		G4SDManager* sdManager = G4SDManager::GetSDMpointer();
		sdManager->AddNewDetector(detector);

		// Attach selenium logical volume with the new sensitive detector
		seleniumLV->SetSensitiveDetector(detector);
	}



	// returns the physical world
	return physWorld;


}

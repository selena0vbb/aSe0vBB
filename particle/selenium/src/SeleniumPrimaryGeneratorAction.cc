// Primary Generator Action code

#include "SeleniumPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"

SeleniumPrimaryGeneratorAction::SeleniumPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{

	// // Create a GPS and assign it to fParticleGun
	// fParticleGun = new G4GeneralParticleSource();

	// // Get the gamma particle from table
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleType = "gamma";
	G4ParticleDefinition* particle = particleTable->FindParticle(particleType);

	// // Print particle information
	// G4cout << particle << G4endl;

	// // Assign particle to source
	// fParticleGun->SetParticleDefinition(particle);
	// // Optional. Setting particle parameters. Can also do with macro
	// fParticleGun->SetParticleEnergy(122.*keV);
	// fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));

	// Using Particle gun for first pass simplicity
	G4int nOfParticles = 1;
	fParticleGun = new G4ParticleGun(nOfParticles);

	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	fParticleGun->SetParticleEnergy(122.*keV);


}

SeleniumPrimaryGeneratorAction::~SeleniumPrimaryGeneratorAction()
{
	delete fParticleGun;
}

void SeleniumPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
	G4double worldZHalf = 0.;
	// Get the world logical volume
	G4String logicalVolumeName = "World";
	G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume(logicalVolumeName);

	// Gets the box associated with the LV and finds the highest z point
	G4Box* worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
	worldZHalf = worldBox->GetZHalfLength();

	// Set particle source position
	fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., worldZHalf));
	fParticleGun->GeneratePrimaryVertex(event);


}
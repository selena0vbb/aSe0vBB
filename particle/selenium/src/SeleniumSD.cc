// Selenium SD source code

#include "SeleniumSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


SeleniumSD::SeleniumSD(const G4String& name, const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fHitsCollectionID(-1)
{

	collectionName.insert(hitsCollectionName);
}

SeleniumSD::~SeleniumSD()
{

}

void SeleniumSD::Initialize(G4HCofThisEvent* hce)
{
	// Create hits collection
	fHitsCollection = new SeleniumHitsCollection(SensitiveDetectorName, collectionName[0]);

	// Get the collection ID based on the actual collection. Set it to the class member
	if ( fHitsCollectionID < 0)
		fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

	// Add collection to hits collection of this event hce
	hce->AddHitsCollection(fHitsCollectionID, fHitsCollection);

}

G4bool SeleniumSD::ProcessHits(G4Step* step, G4TouchableHistory* history)
{

	// Check that energy was deposited in hit
	G4double edep = step->GetTotalEnergyDeposit();
	if ( edep == 0. ) return false;

	// Create a new hit. Add properties to it and add to the collection
	SeleniumHit* hit = new SeleniumHit();

	// Set track ID, energy deposited, position, and particle of the hit
	hit->SetTrackID(step->GetTrack()->GetTrackID());
	hit->SetParentID(step->GetTrack()->GetParentID());
	hit->SetEdep(edep);
	hit->SetPosition(step->GetPostStepPoint()->GetPosition());
	hit->SetParticleDefinition(step->GetTrack()->GetDefinition());

	// Try to get the creator process
	if(step->GetTrack()->GetCreatorProcess() == 0)
	{
		G4String message = "Primary Particle";
		hit->SetCreatorProcessName(message);
	}
	else
	{
		hit->SetCreatorProcessName(step->GetTrack()->GetCreatorProcess()->GetProcessName());
	}

	// Add hit to hit collection
	fHitsCollection->insert(hit);

	return true;

}

void SeleniumSD::EndOfEvent(G4HCofThisEvent* hce)
{
	if ( verboseLevel > 1 ) {
     	G4int nofHits = fHitsCollection->entries();
     	G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
	     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  	}

}

// Selenium hit source code

#include "SeleniumHit.hh"
#include "G4UnitsTable.hh"

#include <iomanip>


G4ThreadLocal G4Allocator<SeleniumHit>* SeleniumHitAllocator=0;

SeleniumHit::SeleniumHit()
 : G4VHit(),
   fTrackID(-1),
   fEdep(0.),
   fPos(G4ThreeVector()),
   fPosInit(G4ThreeVector()),
   fSecondaryTrack()
{

}

SeleniumHit::~SeleniumHit()
{

}

// Overloaded constructor. Accepts an instance of SeleniumHit. Copies member values over
SeleniumHit::SeleniumHit(const SeleniumHit& right) : G4VHit()
{
	fTrackID 			= right.fTrackID;
	fEdep 				= right.fEdep;
	fPos 				= right.fPos;
	fPosInit			= right.fPosInit;
	fParticleDefinition = right.fParticleDefinition;
	fSecondaryTrack   = right.fSecondaryTrack;
}

// Overload the equals operator. Sets all member variables to be the same
const SeleniumHit& SeleniumHit::operator=(const SeleniumHit& right)
{
	fTrackID  		    = right.fTrackID;
	fEdep     		    = right.fEdep;
  	fPos       			= right.fPos;
  	fPosInit			= right.fPosInit;
  	fParticleDefinition = right.fParticleDefinition;
  	fSecondaryTrack   =  right.fSecondaryTrack;

  	return *this;
}

// Overload the comparison operator
G4int SeleniumHit::operator==(const SeleniumHit& right) const
{
	return ( this == &right) ? 1 : 0;
}

// Print function
void SeleniumHit::Print()
{
	G4cout
		<< " trackID: " << fTrackID << "\n"
		<< " Particle type: " << fParticleDefinition->GetParticleName() << "\n"
		<< " Edep: " << std::setw(7) << G4BestUnit(fEdep, "Energy") << "\n"
		<< " Initial Position: " << std::setw(7) << G4BestUnit( fPosInit, "Length") << "\n"
		<< " Final Position: " << std::setw(7) << G4BestUnit( fPos, "Length") << "\n"
		<< " # of Secondaries: " << fSecondaryTrack.size() << "\n"
		<< G4endl;
}
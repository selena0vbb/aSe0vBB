// Selenium hit header

#ifndef SeleniumHit_H
#define SeleniumHit_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "tls.hh"

class SeleniumHit : public G4VHit
{
	public:
		SeleniumHit();
		SeleniumHit(const SeleniumHit&);
		virtual ~SeleniumHit();

		// Overloading certain operators
		const SeleniumHit& operator=(const SeleniumHit&);
		G4int operator==(const SeleniumHit&) const;

		inline void* operator new(size_t);
		inline void operator delete(void *);

		virtual void Print();

		// Set methods
		void SetTrackID(G4int track) { fTrackID = track; };
		void SetParentID(G4int parent) { fParentID = parent; };
		void SetEdep(G4double energy) { fEdep = energy; };
		void SetPosition(G4ThreeVector xyz) { fPos = xyz; };
		void SetInitialPosition(G4ThreeVector xyz) { fPosInit = xyz; };
		void SetParticleDefinition(G4ParticleDefinition* particle) { fParticleDefinition = particle; };
		void SetCreatorProcessName(G4String processName) { fCreatorProcessName = processName; };
		void SetSecondaryTrack(std::vector<const G4Track * > & secondaryTrack) { fSecondaryTrack = secondaryTrack; };

		// Get methods
		G4int GetTrackID() const { return fTrackID; };
		G4int GetParentID() const { return fParentID; };
		G4double GetEdep() const { return fEdep; };
		G4ThreeVector GetPosition() const { return fPos; };
		G4ThreeVector GetInitialPosition() const { return fPosInit; };
		G4ParticleDefinition* GetParticleDefinition() const{ return fParticleDefinition; };
		G4String GetCreatorProcessName() const { return fCreatorProcessName; };
		std::vector<const G4Track *> GetSecondaryTrack() { return fSecondaryTrack; };

	private:
	
		// Set up member of track ID, position and energy
		G4int fTrackID;
		G4int fParentID;
		G4double fEdep;
		G4ThreeVector fPos;
		G4ThreeVector fPosInit;
		G4ParticleDefinition* fParticleDefinition;
		G4String fCreatorProcessName;
		std::vector<const G4Track *> fSecondaryTrack;
};

typedef G4THitsCollection<SeleniumHit> SeleniumHitsCollection;
extern G4ThreadLocal G4Allocator<SeleniumHit>* SeleniumHitAllocator;

inline void* SeleniumHit::operator new(size_t)
{
	if (!SeleniumHitAllocator)
		SeleniumHitAllocator = new G4Allocator<SeleniumHit>;
	return (void *)SeleniumHitAllocator->MallocSingle();
}

inline void SeleniumHit::operator delete(void *hit)
{
	SeleniumHitAllocator->FreeSingle((SeleniumHit*) hit);
}
#endif
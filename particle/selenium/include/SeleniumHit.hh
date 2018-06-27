// Selenium hit header

#ifndef SeleniumHit_H
#define SeleniumHit_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
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
		void SetEdep(G4double energy) { fEdep = energy; };
		void SetPosition(G4ThreeVector xyz) { fPos = xyz; };
		void SetParticleDefinition(G4ParticleDefinition* particle) { fParticleDefinition = particle; }

		// Get methods
		G4int GetTrackID() const { return fTrackID; };
		G4double GetEdep() const { return fEdep; };
		G4ThreeVector GetPosition() const { return fPos; };
		G4ParticleDefinition* GetParticleDefinition() const{ return fParticleDefinition; };

	private:
	
		// Set up member of track ID, position and energy
		G4int fTrackID;
		G4double fEdep;
		G4ThreeVector fPos;
		G4ParticleDefinition* fParticleDefinition;
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
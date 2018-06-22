// Selenium Primary Generator Action header file

#ifndef SeleniumPrimaryGeneratorAction_H
#define SeleniumPrimaryGeneratorAction_H

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

//class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;

class SeleniumPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

	public:
		SeleniumPrimaryGeneratorAction();
		virtual ~SeleniumPrimaryGeneratorAction();

		virtual void GeneratePrimaries(G4Event* event);

		//const G4GeneralParticleSource* GetParticleGun() const {return fParticleGun;};
		const G4ParticleGun* GetParticleGun() const { return fParticleGun; };

	private:
		//G4GeneralParticleSource* fParticleGun; // G4 Particle gun pointer
		G4ParticleGun* fParticleGun;
};

#endif
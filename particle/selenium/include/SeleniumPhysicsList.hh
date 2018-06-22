// Selenium Detector Physics List

#ifndef SeleniumPhysicsList_H
#define SeleniumPhysicsList_H

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

class SeleniumPhysicsList: public G4VModularPhysicsList
{
	public:
		SeleniumPhysicsList();
		virtual ~SeleniumPhysicsList();

		virtual void ConstructParticle();
		virtual void ConstructProcess();

		virtual void SetCuts();

	private:
		G4VPhysicsConstructor* fEmPhysics;	
};


#endif
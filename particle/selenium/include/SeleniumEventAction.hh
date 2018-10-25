// Selenium Event Action definitions

#ifndef SeleniumEventAction_H
#define SeleniumEventAction_H

#include "G4UserEventAction.hh"
#include "globals.hh"

class SeleniumEventAction : public G4UserEventAction
{
	public:
		SeleniumEventAction();
		virtual ~SeleniumEventAction();

		virtual void BeginOfEventAction(const G4Event* event);
		virtual void EndOfEventAction(const G4Event* event);

		// Member functions to add energy to the materials
		void addEnergyDepSe(G4double energy){fEnergyDepositSe += energy;}
		void addEnergyDepCuPlate(G4double energy){fEnergyDepositCuPlate += energy;}
		void addEnergyDepSiO2(G4double energy){fEnergyDepositSiO2 += energy;}


	private:
		G4double fEnergyDepositSe;
		G4double fEnergyDepositCuPlate;
		G4double fEnergyDepositSiO2;

};

#endif

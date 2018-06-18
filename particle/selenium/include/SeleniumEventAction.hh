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

		void addEnergyDepAu(G4double energy){fEnergyDepositAu += energy;}
		void addEnergyDepSe(G4double energy){fEnergyDepositSe += energy;}


	private:
		G4double fEnergyDepositSe;
		G4double fEnergyDepositAu;

};

#endif
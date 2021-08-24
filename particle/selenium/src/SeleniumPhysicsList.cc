// Definition of the user defined physiscs list for the 

#include "SeleniumPhysicsList.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

SeleniumPhysicsList::SeleniumPhysicsList() 
 : G4VModularPhysicsList(), 
   fEmPhysics(nullptr)
{
	// Set the default cut value
	SetDefaultCutValue(2.*nm);

	SetVerboseLevel(1);

	// Creates the Penelope physics list
	fEmPhysics = new G4EmPenelopePhysics();
}

SeleniumPhysicsList::~SeleniumPhysicsList()
{

}

// Include headers of particles that we need 

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Baryons
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

// Nuclei
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

void SeleniumPhysicsList::ConstructParticle()
{

	// // Imaginary tracking particles
	// G4Geantino::GeantinoDefinition();
	// G4ChargedGeantino::ChargedGeantinoDefinition();

	// // Gamma
	// G4Gamma::GammaDefinition();

	// // Optical photon
	// G4OpticalPhoton::OpticalPhotonDefinition();

	// // Leptons
	// G4Electron::ElectronDefinition();
	// G4Positron::PositronDefinition();
	// G4MuonPlus::MuonPlusDefinition();
	// G4MuonMinus::MuonMinusDefinition();

	// G4NeutrinoE::NeutrinoEDefinition();
	// G4AntiNeutrinoE::AntiNeutrinoEDefinition();

	// // Barions
	// G4Proton::ProtonDefinition();
	// G4AntiProton::AntiProtonDefinition();
	// G4Neutron::NeutronDefinition();
	// G4AntiNeutron::AntiNeutronDefinition();

	// G4Deuteron::DeuteronDefinition();
	// G4Triton::TritonDefinition();
	// G4Alpha::AlphaDefinition();
	// G4GenericIon::GenericIonDefinition();
	fEmPhysics->ConstructParticle();

}

void SeleniumPhysicsList::ConstructProcess()
{

	// Transportation
    AddTransportation();

	// Register the penelope physics list
	fEmPhysics->ConstructProcess();
}

void SeleniumPhysicsList::SetCuts()
{
	// Set limits of cut
	G4double lowerLimit = 1*eV;
	G4double upperLimit = 100*GeV;
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowerLimit, upperLimit);

	// Set cuts by calling parent function. Sets to default value
	G4VUserPhysicsList::SetCuts();

	// Print the cut values if verbosity dictates
	if( verboseLevel>0 )
	{
		DumpCutValuesTable();
	}



}
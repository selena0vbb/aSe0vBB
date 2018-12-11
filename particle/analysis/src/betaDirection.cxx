#include "betaDirection.hh"

using namespace std;

TTree * createEventTree(TTree *tree){
/*
 * Turns the Geant4 Sensitive detector tree (every interaction is a separate leaf, no distiguishing of track except for 
 * trackID) into a tree where each leaf is a separate event and the positions are stored in TVector3s	
 */

	double xraw, yraw, zraw, eraw;
	int trackID;
	vector<TVector3> xyz;
	vector<double> x, y, z, e;
	int eventID;
	double totalE;
	TVector3 currentPos, initialPos;

	// Set the current tree branch addresses
	tree->SetBranchAddress("x", &xraw);
	tree->SetBranchAddress("y", &yraw);
	tree->SetBranchAddress("z", &zraw);
	tree->SetBranchAddress("energy", &eraw);
	tree->SetBranchAddress("TrackID", &trackID);

	int n = tree->GetEntries();

	// Set up the next tree
	TTree *eventTree = new TTree("aSeEvents", "aSeEvents");
	eventTree->Branch("eventID", &eventID);
	eventTree->Branch("x", &x);
	eventTree->Branch("y", &y);
	eventTree->Branch("z", &z);
	eventTree->Branch("energy", &e);
	eventTree->Branch("totalEnergy", &totalE);

	int currentTrackID = -1;
	int eventCounter = 0;

	// Iterate over the old tree and parse data into new structures
	for(int i=0; i<n; i++){
		tree->GetEntry(i);

		if(trackID != currentTrackID){
			if(currentTrackID == -1){
				currentTrackID = trackID;
			}else{
				eventID = eventCounter;
				tree->Fill();
				currentTrackID = trackID;
				eventCounter++;
				x.clear();
				y.clear();
				z.clear();
				e.clear();
				totalE = 0;
			}

		}

		totalE += eraw*keV;

		// Fill the vectors
		x.push_back(xraw);
		y.push_back(yraw);
		z.push_back(zraw);
		e.push_back(eraw*keV);		
	
	}

	return eventTree;

}


TH2D * plotBetaDirection(TTree *eventTree, double emin){
	
	
	TVector3  initialPosition, currentPosition, deltaPosition;
	vector<double> energy, x, y, z;
	double totalEnergy;
	int eventID;

	eventTree->SetBranchAddress("x", &x);
	eventTree->SetBranchAddress("y", &y);
	eventTree->SetBranchAddress("z", &z);
	eventTree->SetBranchAddress("energy", &energy);
	eventTree->SetBranchAddress("totalEnergy", &totalEnergy);
	eventTree->SetBranchAddress("eventID", &eventID);

	// Create histogram to store data in
	TH2D *hBeta = new TH2D("betaTrack", "Orientation of Beta Tracks xz", 200, -0.1, 0.1, 200, -0.1, 0.1);
	int nEvents = eventTree->GetEntries();

	// Iterate over events to create histogram of distribution of beta tracks
	for(int i=0; i<nEvents; i++){

		eventTree->GetEntry(i);
		initialPosition.SetXYZ(x[0], y[0], z[0]);
		if(totalEnergy*keV >= emin){
			for(int j=0; j<x.size(); j++){
				currentPosition.SetXYZ(x[j], y[j], z[j]);
				deltaPosition = currentPosition - initialPosition;
				hBeta->Fill(deltaPosition.X(), deltaPosition.Z(), energy[j]);
			}
		}
	}

	hBeta->Draw();
	return hBeta;

}

void savehist(TH2D *h, const char *filename, const char *histname){
/* saves a histogram to a file with name histname
 */

	TFile *f = new TFile(filename, "UPDATE");
	h->Write(histname);
	delete f;

	return;
}

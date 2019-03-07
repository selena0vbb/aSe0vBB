#include "betaDirection.hh"

using namespace std;

TTree * createEventTree(TTree *tree, const char *outfile){
/*
 * Turns the Geant4 Sensitive detector tree (every interaction is a separate leaf, no distiguishing of track except for 
 * trackID) into a tree where each leaf is a separate event and the positions are stored in TVector3s	
 */

	// Create file to save the tree to
	TFile *tOutfile = TFile::Open(outfile, "RECREATE");

	double xraw, yraw, zraw, eraw;
	int EventIDraw;
	vector<TVector3> xyz;
	vector<double> x, y, z, e, delx, dely, delz, delrho, delphi, deltheta;
	int eventID;
	double totalE;
	TVector3 currentPosition, initialPosition;

	// Set the current tree branch addresses
	tree->SetBranchAddress("x", &xraw);
	tree->SetBranchAddress("y", &yraw);
	tree->SetBranchAddress("z", &zraw);
	tree->SetBranchAddress("energy", &eraw);
	tree->SetBranchAddress("EventID", &EventIDraw);

	int n = tree->GetEntries();

	// Set up the next tree
	TTree *eventTree = new TTree("aSeEvents", "aSeEvents");
	eventTree->Branch("eventID", &eventID);
	eventTree->Branch("x", &x);
	eventTree->Branch("y", &y);
	eventTree->Branch("z", &z);
	eventTree->Branch("delx", &delx);
	eventTree->Branch("delx", &delx);
	eventTree->Branch("dely", &dely);
	eventTree->Branch("delz", &delz);
	eventTree->Branch("delrho", &delrho);
	eventTree->Branch("delphi", &delphi);
	eventTree->Branch("deltheta", &deltheta);
	eventTree->Branch("energy", &e);
	eventTree->Branch("totalEnergy", &totalE);

	int currentEventID = -1;
	int eventCounter = 0;
	cout << n << endl;
	// Iterate over the old tree and parse data into new structures
	for(int i=0; i<n; i++){
		tree->GetEntry(i);
		if(i%10000 == 0) cout << i << endl;
		if(EventIDraw != currentEventID){
			if(currentEventID > -1){
				eventID = eventCounter;
				eventTree->Fill();
				eventCounter++;
				x.clear(); y.clear(); z.clear(); e.clear();
				delx.clear(); dely.clear(); delz.clear(); delrho.clear(); delphi.clear(); deltheta.clear();				
				totalE = 0;
			}
			currentEventID = EventIDraw;
			initialPosition.SetXYZ(xraw, yraw, zraw);

		}
		currentPosition.SetXYZ(xraw, yraw, zraw);
		totalE += eraw*keV;

		// Fill the vectors
		x.push_back(xraw);
		y.push_back(yraw);
		z.push_back(zraw);
		delx.push_back((currentPosition-initialPosition).X());
		dely.push_back((currentPosition-initialPosition).Y());
		delz.push_back((currentPosition-initialPosition).Z());
		delrho.push_back((currentPosition-initialPosition).Perp());
		delphi.push_back((currentPosition-initialPosition).Phi());
		deltheta.push_back((currentPosition-initialPosition).Theta());
		e.push_back(eraw*keV);		
	
	}
	eventTree->Write();
	delete tOutfile;
	return eventTree;

}


TH2D * plotBetaDirection(TTree *eventTree, double emin){
	
	cout << 1 << endl;	
	TVector3  initialPosition, currentPosition, deltaPosition;
	vector<double> *energy = 0;
	vector<double> *x = 0;
	vector<double> *y = 0;
	vector<double> *z = 0;
	double totalEnergy;
	int eventID;

	eventTree->SetBranchAddress("x", &x);
	eventTree->SetBranchAddress("y", &y);
	eventTree->SetBranchAddress("z", &z);
	eventTree->SetBranchAddress("energy", &energy);
	eventTree->SetBranchAddress("totalEnergy", &totalEnergy);
	eventTree->SetBranchAddress("eventID", &eventID);
	cout << 2 << endl;
	// Create histogram to store data in
	TH2D *hBeta = new TH2D("betaTrack", "Orientation of Beta Tracks xz", 200, -0.1, 0.1, 200, -0.1, 0.1);
	int nEvents = eventTree->GetEntries();
	cout << 3 << endl;
	// Iterate over events to create histogram of distribution of beta tracks
	for(int i=0; i<nEvents; i++){

		eventTree->GetEntry(i);
		initialPosition.SetXYZ(x->at(0), y->at(0), z->at(0));
		if(totalEnergy*keV >= emin){
			for(int j=0; j<x->size(); j++){
				currentPosition.SetXYZ(x->at(j), y->at(j), z->at(j));
				deltaPosition = currentPosition - initialPosition;
				hBeta->Fill(deltaPosition.X(), deltaPosition.Z(), energy->at(j));
			}
		}
	}

	return hBeta;

}

TCanvas * plotBetaCanvas(TTree *eventTree, const char *selection){
// Plotting multiple different variable slices

	TCanvas *c = new TCanvas("c","Beta Track Direction", 800, 600);
	TStyle *gStyle = new TStyle();
	gStyle->SetPalette(53);
	gStyle->SetOptStat(0);
	c->Divide(2,2);

	// Create histograms
	TH2D *hxz = new TH2D("hxz", "XZ Projection", 500, -0.04, 0.04, 500, -0.04, 0.04);
	TH2D *hxy = new TH2D("hxy", "XY Projection", 500, -0.04, 0.04, 500, -0.04, 0.04);
	TH2D *hpz = new TH2D("hpz", "Phi Z Projection", 500, -TMath::Pi(), TMath::Pi(), 500, -0.04, 0.04);
	TH2D *hpt = new TH2D("hpt", "Phi Theta Proejction", 500, -TMath::Pi(), TMath::Pi(), 500, 0, TMath::Pi());

	// Set axis labels
	hxz->GetXaxis()->SetTitle("x (mm)");
	hxz->GetYaxis()->SetTitle("z (mm)");
	hpz->GetXaxis()->SetTitle("#phi");
	hpz->GetYaxis()->SetTitle("z (mm)");
	hpt->GetXaxis()->SetTitle("#phi");
	hpt->GetYaxis()->SetTitle("#theta");
	hxy->GetXaxis()->SetTitle("x (mm)");
	hxy->GetYaxis()->SetTitle("y (mm)");

	// Remove stat boxs
	hxz->SetStats(0);
	hpz->SetStats(0);
	hpt->SetStats(0);
	hxy->SetStats(0);

	// Plot histograms
	c->cd(1);
	eventTree->Draw("delz:delx>>hxz", selection, "colz");
	gPad->SetLogz(1);
	c->cd(2);
	eventTree->Draw("dely:delx>>hxy", selection, "colz");
	gPad->SetLogz(1);
	c->cd(3);
	eventTree->Draw("delz:delphi>>hpz", selection, "colz");
	gPad->SetLogz(1);
	c->cd(4);
	eventTree->Draw("deltheta:delphi>>hpt", selection, "colz");
	gPad->SetLogz(1);

	return c;
}

void savehist(TH2D *h, const char *filename, const char *histname){
/* saves a histogram to a file with name histname
 */

	TFile *f = new TFile(filename, "UPDATE");
	h->Write(histname);
	delete f;

	return;
}

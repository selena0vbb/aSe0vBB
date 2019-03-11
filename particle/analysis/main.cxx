#include "betaDirection.hh"

using namespace std;

int main(int argc, char *argv[]){

	TFile *f0 = new TFile("~/data/particle/pixel_sio2_122kev_0degangle_200k.root");
	TFile *f75 = new TFile("~/data/particle/pixel_sio2_122kev_75degangle_200k.root");

	TTree *t0 = (TTree*)f0->Get("aSeData");
	TTree *t75 = (TTree*)f75->Get("aSeData");
	
	
	// Make the event trees
	cout << "Get aSe Event tree." << endl;
	//TFile *f0 = new TFile("~/aSe0vBB/particle/analysis/output/0degEventTree.root");
	//TFile *f75 = new TFile("~/aSe0vBB/particle/analysis/output/75degEventTree.root");
	//TTree *t0event = (TTree*)f0->Get("aSeEvents");
	//TTree *t75event = (TTree*)f75->Get("aSeEvents"); 
	TTree *t0event = createEventTree(t0, "~/aSe0vBB/particle/analysis/output/0degEventTree.root");
	TTree *t75event = createEventTree(t75, "~/aSe0vBB/particle/analysis/output/75degEventTree.root");
/*	cout << "Plot directionality of beta particles." << endl;
	TH2D *h0 = plotBetaDirection(t0event, 100);
	cout << 2 << endl;
	TH2D *h75 = plotBetaDirection(t75event, 100);
	
	const char *outfilename = "~/aSe0vBB/particle/analysis/output/testTH2D.root";
	const char *h0Name = "0deg";
	const char *h75Name = "75deg";
	
	cout << "Save histograms." << endl;
	savehist(h0, outfilename, h0Name);
	savehist(h75, outfilename, h75Name);
*/
	cout << "Analysis completed." << endl;
	return 0;
}

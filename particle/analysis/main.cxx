#include "betaDirection.hh"

using namespace std;

int main(int argc, char *argv[]){
	
	// Open tree files
	TFile *f0 = new TFile("~/data/particle/pixel_sio2_122kev_0degangle_200k.root");
	TFile *f75 = new TFile("~/data/particle/pixel_sio2_122kev_75degangle_200k.root");

	// Make the trees
	TTree *t0 = (TTree*)f0->Get("aSeData");
	TTree *t75 = (TTree*)f75->Get("aSeData");
	
	// Make the event trees
	TTree *t0event = createEventTree(t0);
	TTree *t75event = createEventTree(t75); 

	TH2D *h0 = plotBetaDirection(t0event, 80);
	TH2D *h75 = plotBetaDirection(t75event, 80);
	
	const char *outfilename = "./output/testTH2D.root";
	const char *h0Name = "0deg";
	const char *h75Name = "75deg";

	savehist(h0, outfilename, h0Name);
	savehist(h75, outfilename, h75Name);

	cout << "Analysis completed." << endl;
	return 0;
}

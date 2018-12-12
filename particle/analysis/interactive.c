{
	gInterpreter->AddIncludePath("./include/");
	gROOT->ProcessLine(".L ./src/betaDirection.cxx");

	// Open tree files
	TFile *f0 = new TFile("~/data/particle/pixel_sio2_122kev_0degangle_200k.root");
	TFile *f75 = new TFile("~/data/particle/pixel_sio2_122kev_75degangle_200k.root");

	// Make the trees
	TTree *t0 = (TTree*)f0->Get("aSeData");
	TTree *t75 = (TTree*)f75->Get("aSeData");

	TTree *t0event = createEventTree(t0, "~/aSe0vBB/particle/analysis/output/0degEventTree");
	TTree *t75event = createEventTree(t75, "~/aSe0vBB/particle/analysis/output/75degEventTree"); 

}

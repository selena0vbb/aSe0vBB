{
	gInterpreter->AddIncludePath("./include/");
	gROOT->ProcessLine(".L ./src/betaDirection.cxx");

	// Open tree files
	TFile *f0 = new TFile("~/aSe0vBB/particle/analysis/output/0degEventTree.root");
	TFile *f75 = new TFile("~/aSe0vBB/particle/analysis/output/75degEventTree.root");

	// Make the trees
	TTree *t0 = (TTree*)f0->Get("aSeEvents");
	TTree *t75 = (TTree*)f75->Get("aSeEvents");

}

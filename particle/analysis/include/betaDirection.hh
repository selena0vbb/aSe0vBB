#ifndef betaDirection_H
#define betaDirection_H

#include <iostream>
#include <vector>
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"
#include "TPad.h"
#include "TStyle.h"
#include "seleniumConstants.hh"

TTree * createEventTree(TTree *tree, const char *outfile);

TH2D * plotBetaDirection(TTree *eventTree, double emin);

TCanvas * plotBetaCanvas(TTree *eventTree, const char *selection);

void savehist(TH2D *h, const char *filename, const char *histname);


#endif 

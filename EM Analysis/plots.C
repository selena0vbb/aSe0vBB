#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2.h"
#include "TPolyMarker3D.h"
using namespace std;

void plots(){
	SetMyStyle();
	gStyle->SetPalette();
	//Load files
	TFile *froot = new TFile("circularEMap.root","RECREATE");
	TTree *tCenter = new TTree("tcenter","tcenter");
	tCenter->ReadFile("centralaxis_prescan.txt");
	TTree *tVeff = new TTree("tVeff","tVeff");
	tVeff->ReadFile("Veff.txt");
	TTree *tEtrue = new TTree("tEtrue","tEtrue");
	tEtrue->ReadFile("Etrue.txt");
	tCenter->Write();	//r, z, map[6][8][4][2]
	tVeff->Write();		//r, z, map[6][8][2]
	tEtrue->Write();	//r, z, map[6][8][2][2]
	//map[6]:	r_g=10,20,50,100,50,100
	//			r_r=100,100,100,100,200,200
	//map[i][8]:r_c=15,40,65,90,115,140,165,190
	//map[i][j][4]:c1,c2,4000V+200V,4000V+500V; 0123, 01,    23
	//map[i][j][k][2]:							V/Ez, - , Ez/Er
	double rc[8] ={25,50,75,100,125,150,175,200 };
	double rg[6] ={10,20,50,100,50,100 };
	double rr[6] ={100,100,100,100,200,200 };
	
	if(false){
		double scale[6][8];	
		TCanvas *c0 = new TCanvas("c0","c0");
		for(int irg=0; irg<6; irg++){
			for(int irc=0; irc<8; irc++){
				tCenter->Draw(Form("map[%d][%d][0][0]:z",irg,irc),"r==0","l");
				double slopC1 = ((TGraph*)c0->GetPrimitive("Graph"))->Eval(20);
				tCenter->Draw(Form("map[%d][%d][1][0]:z",irg,irc),"r==0","l");
				double slopC2 = ((TGraph*)c0->GetPrimitive("Graph"))->Eval(20);
				scale[irg][irc] = slopC1/slopC2;
			}
		}
		TCanvas **c = new TCanvas*[6];//("c1","c1");
		for(int irg=0; irg<6; irg++){
			c[irg] = new TCanvas(Form("c%d",irg),Form("c%d",irg));
			tCenter->Draw(Form("map[%d][%d][0][0]-map[%d][%d][1][0]*%lf:z",irg,0,irg,0,scale[irg][0]),"r==0","l PLC PMC");
			tCenter->Draw(Form("map[%d][%d][3][1]*0.2e-7:z",irg,0),"r==0","lsame PLC PMC");
			for(int irc=1; irc<8; irc++){
				tCenter->Draw(Form("map[%d][%d][0][0]-map[%d][%d][1][0]*%lf:z",irg,irc,irg,irc,scale[irg][irc]),"r==0","lsame PLC PMC");
				tCenter->Draw(Form("map[%d][%d][3][1]*0.2e-7:z",irg,irc),"r==0","lsame PLC PMC");
			//	tCenter->Draw(Form("map[2][%d][1][0]:z",irc),"r==0","lsame PLC PMC");
			}
		}
	}else{
		double z, r, map[6][8][2][2], mapV[6][8][2];
		tEtrue->SetBranchAddress("z",&z);
		tEtrue->SetBranchAddress("r",&r);
		tEtrue->SetBranchAddress("map",map);
		tVeff->SetBranchAddress("map",mapV);
		long int N = tEtrue->GetEntries();

		double scale[6][8];	
		TCanvas *c0 = new TCanvas("c0","c0");
		for(int irg=0; irg<6; irg++){
			for(int irc=0; irc<8; irc++){
				tVeff->Draw(Form("map[%d][%d][0]:z",irg,irc),"r==0","l");
				double slopC1 = ((TGraph*)c0->GetPrimitive("Graph"))->Eval(20);
				tVeff->Draw(Form("map[%d][%d][1]:z",irg,irc),"r==0","l");
				double slopC2 = ((TGraph*)c0->GetPrimitive("Graph"))->Eval(20);
				scale[irg][irc] = slopC1/slopC2;
			}
		}
		double dtd=0.005;//us
		double mu =0.1;//um^2/us/V
		TCanvas *cE =new TCanvas("cE","cE");
		cE->Divide(2,2);
		TCanvas **c = new TCanvas*[6];//("c1","c1");
		TCanvas **cs = new TCanvas*[6];//("c1","c1");
		cE->Print("sig.pdf[");
		for(int irg=0; irg<6; irg++){
			c[irg] = new TCanvas(Form("c%d",irg),Form("c%d",irg));
			cs[irg] = new TCanvas(Form("cs%d",irg),Form("cs%d",irg));
			//tVeff->Draw(Form("map[%d][%d][3][1]*0.2e-7:z",irg,0),"r==0","lsame PLC PMC");
			for(int irc=0; irc<8; irc++){
				c[irg]->cd();
				tVeff->Draw(Form("map[%d][%d][0]-map[%d][%d][1]*%lf:z",irg,irc,irg,irc,scale[irg][irc]),Form("r==%d", 0),irc==0?"l PLC PMC":"lsame PLC PMC");
			//	tCenter->Draw(Form("map[%d][%d][3][1]*0.2e-7:z",irg,irc),"r==0","lsame PLC PMC");
			//	tCenter->Draw(Form("map[2][%d][1][0]:z",irc),"r==0","lsame PLC PMC");
				c[irg]->Update();
				TLine *lrc = new TLine(0,200,rc[irc],200);
				lrc->SetLineWidth(2);
				lrc->SetLineColor(kBlue);
				TLine *lrr = new TLine(rc[irc]+rg[irg],200,rc[irc]+rg[irg]+rr[irg],200);
				lrr->SetLineWidth(2);
				lrr->SetLineColor(kRed);
		TGraph2D *g2DEr = new TGraph2D(201*101);
		g2DEr->SetNpx(201);g2DEr->SetNpy(101);
		g2DEr->SetNameTitle("g2DEr", "Er [V/m]; r [#mum]; z [#mum]");
		TGraph2D *g2DEz = new TGraph2D(201*101);
		g2DEz->SetNpx(201);g2DEz->SetNpy(101);
		g2DEz->SetNameTitle("g2DEz", "Ez [V/m]; r [#mum]; z [#mum]");
		TGraph2D *g2DV1 = new TGraph2D(201*101);
		g2DV1->SetNpx(201);g2DV1->SetNpy(101);
		g2DV1->SetNameTitle("g2DV1", "V_{c1eff}; r [#mum]; z [#mum]");
		TGraph2D *g2DV2 = new TGraph2D(201*101);
		g2DV2->SetNpx(201);g2DV2->SetNpy(101);
		g2DV2->SetNameTitle("g2DV2", "V_{c2eff}; r [#mum]; z [#mum]");
				for(int iE=0; iE<N; iE++){
					tEtrue->GetEntry(iE);
					tVeff->GetEntry(iE);
					g2DEz->SetPoint(iE, r, z,  map[irg][irc][1][0]);
					g2DEr->SetPoint(iE, r, z,  map[irg][irc][1][1]);
					g2DV1->SetPoint(iE, r, z, mapV[irg][irc][0]);
					g2DV2->SetPoint(iE, r, z, mapV[irg][irc][1]);
				}
				cE->cd(1);
				g2DEz->Draw("colz");
				lrr->Draw("same");lrc->Draw("same");
				cE->cd(2);
				g2DEr->Draw("colz");
				lrr->Draw("same");lrc->Draw("same");
				cE->cd(3);
				g2DV1->Draw("colz");
				lrr->Draw("same");lrc->Draw("same");
				cE->cd(4);
				g2DV2->Draw("colz");
				lrr->Draw("same");lrc->Draw("same");
				cE->Update();
				cE->Print("sig.pdf");
				cout<<"start simulation..."<<endl;
				double zhole=0;
				double rhole=0;
				double rmax=0;
				vector<double> vhole;
				TH2D *hSig = new TH2D("hsig","hsig",150,0,150,170,-0.6,1.1);
				for(int ir=0; ir<(rc[irc]+rg[irg]+rr[irg]); ir++){
				//	cout<<ir<<endl;
					zhole=10; rhole = ir; 
					for(double td=0; td<200&&zhole<200; td+=dtd){
						vhole.push_back(g2DV1->Interpolate(rhole, zhole) - g2DV2->Interpolate(rhole, zhole)*scale[irg][irc]);
						zhole+=g2DEz->Interpolate(rhole, zhole)*1e-6*mu*dtd;
						rhole+=g2DEr->Interpolate(rhole, zhole)*1e-6*mu*dtd;
						rhole = rhole<0?-rhole:rhole;
					}
					//	cout<<"r"<<rhole<<endl;
					if(rhole < rc[irc] && zhole>=200){
						for(double td=0; vhole.size()>0; td+=dtd){
							hSig->Fill(td,vhole.at(0), ir);
							vhole.erase(vhole.begin());
						}
						rmax = rmax<ir? ir:rmax;
					}
				}
				cs[irg]->cd();
				hSig->SetNameTitle("hSig", Form("rc=%0.0lf #mum, rg=%0.0lf #mum, rr=%0.0lf #mum, max collection radius=%0.2lf;t[arbitary];V",rc[irc],rg[irg],rr[irg], rmax));
				hSig->Draw("COLZ");
				cs[irg]->Print("sig.pdf");
				cs[irg]->Update();
			}
		}
		cs[0]->Print("sig.pdf]");
	
	}



}

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2.h"
#include "TPolyMarker3D.h"
#include "TLine.h"
using namespace std;

void plotsStrip(){
	SetMyStyle();
	gStyle->SetPalette();
	//Load files
	TFile *froot = new TFile("circularEMapStrip.root","RECREATE");
	TTree *tCenter = new TTree("tcenter","tcenter");
	tCenter->ReadFile("strip_center.txt");
	TTree *tVeff = new TTree("tVeff","tVeff");
	tVeff->ReadFile("stripVeff.txt");
	TTree *tEtrue = new TTree("tEtrue","tEtrue");
	tEtrue->ReadFile("stripEtrue.txt");
	tCenter->Write();	//r, z, map[6][8][4][2]
	tVeff->Write();		//r, z, map[6][8][2]
	tEtrue->Write();	//r, z, map[6][8][2][2]
	//map[6]:	r_g=10,20,50,100,50,100
	//			r_r=20,20,20,20 ,50,50
	//map[i][8]:r_c=20,30,40,50,60,80,100,150
	//map[i][j][4]:c1,c2,4000V+200V,4000V+500V; 0123, 01,    23
	//map[i][j][k][2]:							V/Ez, - , Ez/Er
	double rc[8] ={20,30,40,50,60,80,100,150 };
	double rg[6] ={20,30,60,110,75,125};
	double rr[6] ={20,20,20,20,50,50};
	
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
		cE->Print("stripSig.pdf[");
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
				TGraph2D *g2DEr = new TGraph2D(101*101);
				g2DEr->SetNpx(101);g2DEr->SetNpy(101);
				g2DEr->SetNameTitle("g2DEr", "Ex [V/m]; r [#mum]; z [#mum]");
				TGraph2D *g2DEz = new TGraph2D(101*101);
				g2DEz->SetNpx(101);g2DEz->SetNpy(101);
				g2DEz->SetNameTitle("g2DEz", "Ez [V/m]; r [#mum]; z [#mum]");
				TGraph2D *g2DV1 = new TGraph2D(101*101);
				g2DV1->SetNpx(101);g2DV1->SetNpy(101);
				g2DV1->SetNameTitle("g2DV1", "V_{c1eff}; r [#mum]; z [#mum]");
				TGraph2D *g2DV2 = new TGraph2D(101*101);
				g2DV2->SetNpx(101);g2DV2->SetNpy(101);
				g2DV2->SetNameTitle("g2DV2", "V_{c2eff}; r [#mum]; z [#mum]");
				TLine *lrc = new TLine(0,200,rc[irc]*0.5,200);
				lrc->SetLineWidth(2);
				lrc->SetLineColor(kBlue);
				TLine *lrr = new TLine(rc[irc]*0.5+rg[irg]-rr[irg]*0.5,200,rg[irg]+rc[irc]*0.5,200);
				lrr->SetLineWidth(2);
				lrr->SetLineColor(kRed);
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
				lrc->Draw("same"); lrr->Draw("same");
				cE->cd(2);
				g2DEr->Draw("colz");
				lrc->Draw("same"); lrr->Draw("same");
				cE->cd(3);
				g2DV1->Draw("colz");
				lrc->Draw("same"); lrr->Draw("same");
				cE->cd(4);
				g2DV2->Draw("colz");
				lrc->Draw("same"); lrr->Draw("same");
				cE->Update();
				cE->Print("stripSig.pdf");
				cout<<"start simulation..."<<endl;
				double zhole=0;
				double rhole=0;
				double rrec = (rg[irg]+rc[irc]*0.5);
				vector<double> vhole;
				TH2D *hSig = new TH2D("hsig","hsig",150,0,150,170,-0.6,1.1);
				double Eff=0;
				for(int ir=0; ir<rrec*2; ir++){
				//	cout<<ir<<endl;
					zhole=10; rhole = ir*0.5; 
					//if(rhole>rc[irc]*0.5+50) continue;
					int irec = rhole/rrec/2;
					double rhole_rec = rhole - irec*rrec*2;
					rhole_rec = rhole_rec>rrec ? rrec*2-rhole_rec:rhole_rec;
					for(double td=0; td<200&&zhole<200; td+=dtd){
						vhole.push_back(g2DV1->Interpolate(rhole_rec, zhole) - g2DV2->Interpolate(rhole_rec, zhole)*scale[irg][irc]);
						zhole+=g2DEz->Interpolate(rhole_rec, zhole)*1e-6*mu*dtd;
						rhole_rec+=g2DEr->Interpolate(rhole_rec, zhole)*1e-6*mu*dtd;
						rhole_rec = rhole_rec<0?-rhole_rec:rhole_rec;
						rhole_rec = rhole_rec>rrec?2*rrec-rhole_rec:rhole_rec;
					}
					//	cout<<"r"<<rhole<<endl;
					if(rhole_rec < rc[irc]*0.5 && zhole>=200){
						for(double td=0; vhole.size()>0; td+=dtd){
							hSig->Fill(td,vhole.at(0));
							vhole.erase(vhole.begin());
						}
						Eff+=100./rrec/2;
					}
				}
				cs[irg]->cd();
				hSig->SetNameTitle("hSig", Form("d1=%0.0lf #mum, d2=%0.0lf #mum, dg=%0.0lf, with collection of%0.2lf%%;t[arbitary unit];V",rc[irc],rr[irg],rg[irg]-rr[irg]*0.5,Eff));
				hSig->Draw("COLZ");
				cs[irg]->Print("stripSig.pdf");
				cs[irg]->Update();
			}
		}
		cs[0]->Print("stripSig.pdf]");
	
	}



}

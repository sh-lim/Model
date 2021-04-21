#include "/phenix/u/shlim/Style.h"

void Draw_init_pPb_PbPb(){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);

	TCanvas *c0 = new TCanvas("c0","c0",1.1*2*600,600);
	c0->SetFillColor(1);
	c0->Divide(2,1);

	c0->cd(2);
	SetPadStyle();
	gPad->SetTopMargin(0.0);
	gPad->SetRightMargin(0.0);
	gPad->SetBottomMargin(0.09);
	gPad->SetLeftMargin(0.07);
	
	htmp = (TH1D*)gPad->DrawFrame(-10.5,-10.5,10.5,10.5);
	SetHistoStyle("x [fm]","y [fm]");
	htmp->GetXaxis()->CenterTitle();
	htmp->GetYaxis()->CenterTitle();
	htmp->GetXaxis()->SetTitleColor(0);
	htmp->GetXaxis()->SetTitleFont(63);
	htmp->GetXaxis()->SetTitleSize(24);
	htmp->GetXaxis()->SetLabelColor(0);
	htmp->GetXaxis()->SetLabelFont(63);
	htmp->GetXaxis()->SetLabelSize(24);
	htmp->GetYaxis()->SetTitleColor(0);
	htmp->GetYaxis()->SetTitleFont(63);
	htmp->GetYaxis()->SetTitleSize(24);
	htmp->GetYaxis()->SetLabelColor(0);
	htmp->GetYaxis()->SetLabelFont(63);
	htmp->GetYaxis()->SetLabelSize(24);
	htmp->GetYaxis()->SetTitleOffset(0.8);

	TFile *infile0 = new TFile("mcglauber_outfile_pPb_8160GeV_2_1fm_00.root","read");
	TH2D *hhe3 = (TH2D*)infile0->Get("inited_event1_translated"); //24, 51, 68, 74, 79, 84

	for (int ix=0; ix<hhe3->GetNbinsX(); ix++){
		for (int iy=0; iy<hhe3->GetNbinsY(); iy++){
			if ( hhe3->GetBinContent(ix+1,iy+1)<0.1 ){
				hhe3->SetBinContent(ix+1,iy+1,0);
			}
		}
	}

	hhe3->Draw("col same");

	TLatex *tex = new TLatex(-8.5,8.5,"#color[0]{p+Pb b=0.9 fm}");
	tex->SetTextFont(63);
	tex->SetTextSize(32);
	tex->Draw();

	//return;


	c0->cd(1);
	SetPadStyle();
	gPad->SetTopMargin(0.0);
	gPad->SetRightMargin(0.0);
	gPad->SetBottomMargin(0.09);
	gPad->SetLeftMargin(0.07);
	
	htmp = (TH1D*)gPad->DrawFrame(-10.5,-10.5,10.5,10.5);
	SetHistoStyle("x [fm]","y [fm]");
	htmp->GetXaxis()->CenterTitle();
	htmp->GetYaxis()->CenterTitle();
	htmp->GetXaxis()->SetTitleColor(0);
	htmp->GetXaxis()->SetTitleFont(63);
	htmp->GetXaxis()->SetTitleSize(24);
	htmp->GetXaxis()->SetLabelColor(0);
	htmp->GetXaxis()->SetLabelFont(63);
	htmp->GetXaxis()->SetLabelSize(24);
	htmp->GetYaxis()->SetTitleColor(0);
	htmp->GetYaxis()->SetTitleFont(63);
	htmp->GetYaxis()->SetTitleSize(24);
	htmp->GetYaxis()->SetLabelColor(0);
	htmp->GetYaxis()->SetLabelFont(63);
	htmp->GetYaxis()->SetLabelSize(24);
	htmp->GetYaxis()->SetTitleOffset(0.8);

	TFile *infile1 = new TFile("mcglauber_outfile_PbPb_5020GeV_0_1fm_00.root","read");
	TH2D *hhe4 = (TH2D*)infile1->Get("inited_event7_translated"); //24, 51, 68, 74, 79, 84

	for (int ix=0; ix<hhe4->GetNbinsX(); ix++){
		for (int iy=0; iy<hhe4->GetNbinsY(); iy++){
			if ( hhe4->GetBinContent(ix+1,iy+1)<0.1 ){
				hhe4->SetBinContent(ix+1,iy+1,0);
			}
		}
	}

	hhe4->Draw("col same");

	TLatex *tex = new TLatex(-8.5,8.5,"#color[0]{Pb+Pb b=0.9 fm}");
	tex->SetTextFont(63);
	tex->SetTextSize(32);
	tex->Draw();

}

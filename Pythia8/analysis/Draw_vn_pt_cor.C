#include "Style.h"

void Draw_vn_pt_cor(){

	gStyle->SetOptStat(0);

	TFile *infile = new TFile("/alice/home/shlim/work/Pythia/jobs/outfile_hist_pPb8TeV_set00_grp000_try103.root","read");

	TProfile *hprof_ck = (TProfile*)infile->Get("hprof_ck");
	TProfile *hprof_corr2 = (TProfile*)infile->Get("hprof_corr2");
	TProfile *hprof_corr4 = (TProfile*)infile->Get("hprof_corr4");
	TProfile *hprof_ptcorr2 = (TProfile*)infile->Get("hprof_ptcorr2");

	TH1D *hprof_rho = new TH1D("hprof_rho","",hprof_ck->GetNbinsX(),0,hprof_ck->GetBinWidth(1)*hprof_ck->GetNbinsX());

	for (int ii=0; ii<hprof_ck->GetNbinsX(); ii++){
		float A = hprof_ptcorr2->GetBinContent(ii+1);
		float B2 = hprof_corr2->GetBinContent(ii+1);
		float B4 = hprof_corr4->GetBinContent(ii+1);
		float B = B4 - B2*B2;
		float C = hprof_ck->GetBinContent(ii+1);

		float rho = A / sqrt(B) / sqrt(C);

		cout << A << " " << B << " " << C << " " << rho << endl;
		
		if ( B>0 && C>0 ){
			hprof_rho->SetBinContent(ii+1, rho);
		}
	}


	TCanvas *c1 = new TCanvas("c1","c1",1.1*3*400,2*400);
	c1->Divide(3,2);

	c1->cd(1);
	SetPadStyle();
	gPad->SetTopMargin(0.05);
	hprof_ck->Draw();

	c1->cd(2);
	SetPadStyle();
	gPad->SetTopMargin(0.05);
	hprof_ptcorr2->Draw();

	c1->cd(4);
	SetPadStyle();
	gPad->SetTopMargin(0.05);
	hprof_corr2->Draw();

	c1->cd(5);
	SetPadStyle();
	gPad->SetTopMargin(0.05);
	hprof_corr4->Draw();

	c1->cd(6);
	SetPadStyle();
	gPad->SetTopMargin(0.05);
	hprof_rho->Draw();



}

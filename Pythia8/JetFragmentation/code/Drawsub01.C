#include "Style.h"

void Drawsub01(){

	const int npt = 3;
	const float ptbin[npt+1] = {25, 40, 60, 80};

	TFile *infile = new TFile("../outfile_hist_pp5TeV_set25_grp000_try200.root","read");

	TH1D *jeteta[npt];
	TH1D *perpeta[npt];

	TH1D *incljt[npt];
	TH1D *inclz[npt];
	TH1D *perp_incljt[npt];
	TH1D *perp_inclz[npt];

	TH1D *sub_incljt[npt];
	TH1D *sub_inclz[npt];

	float njet[npt] = {0.};
	float nperpjet[npt] = {0.};

	for (int ii=0; ii<npt; ii++){
		jeteta[ii] = (TH1D*)infile->Get(Form("jeteta_pt%d",ii));
		perpeta[ii] = (TH1D*)infile->Get(Form("perpeta_pt%d",ii));

		int etabin_lo = jeteta[ii]->FindBin(-0.5+0.001);
		int etabin_hi = jeteta[ii]->FindBin(+0.5-0.001);

		njet[ii] = jeteta[ii]->Integral(etabin_lo, etabin_hi);
		nperpjet[ii] = perpeta[ii]->Integral(etabin_lo,etabin_hi);

		incljt[ii] = (TH1D*)infile->Get(Form("incljt_pt%d",ii));
		inclz[ii] = (TH1D*)infile->Get(Form("inclz_pt%d",ii));

		incljt[ii]->Sumw2();
		inclz[ii]->Sumw2();

		perp_incljt[ii] = (TH1D*)infile->Get(Form("perp_incljt_pt%d",ii));
		perp_inclz[ii] = (TH1D*)infile->Get(Form("perp_inclz_pt%d",ii));

		perp_incljt[ii]->Sumw2();
		perp_inclz[ii]->Sumw2();
	}

	for (int ii=0; ii<npt; ii++){
		for (int ijt=0; ijt<incljt[ii]->GetNbinsX(); ijt++){
			float djt = incljt[ii]->GetBinWidth(ijt+1);
			float xx = incljt[ii]->GetBinContent(ijt+1);
			float xx_err = incljt[ii]->GetBinError(ijt+1);

			incljt[ii]->SetBinContent(ijt+1, xx/djt/njet[ii]);
			incljt[ii]->SetBinError(ijt+1, xx_err/djt/njet[ii]);

			xx = perp_incljt[ii]->GetBinContent(ijt+1);
			xx_err = perp_incljt[ii]->GetBinError(ijt+1);

			perp_incljt[ii]->SetBinContent(ijt+1, xx/djt/nperpjet[ii]);
			perp_incljt[ii]->SetBinError(ijt+1, xx_err/djt/nperpjet[ii]);
		}//ijt

		for (int iz=0; iz<inclz[ii]->GetNbinsX(); iz++){
			float dz = inclz[ii]->GetBinWidth(iz+1);
			float xx = inclz[ii]->GetBinContent(iz+1);
			float xx_err = inclz[ii]->GetBinError(iz+1);

			inclz[ii]->SetBinContent(iz+1, xx/dz/njet[ii]);
			inclz[ii]->SetBinError(iz+1, xx_err/dz/njet[ii]);

			xx = perp_inclz[ii]->GetBinContent(iz+1);
			xx_err = perp_inclz[ii]->GetBinError(iz+1);

			perp_inclz[ii]->SetBinContent(iz+1, xx/dz/nperpjet[ii]);
			perp_inclz[ii]->SetBinError(iz+1, xx_err/dz/nperpjet[ii]);
		}//iz

		sub_incljt[ii] = (TH1D*)incljt[ii]->Clone(Form("sub_incljt_pt%d",ii));
		sub_incljt[ii]->Add(perp_incljt[ii], -1);

		sub_inclz[ii] = (TH1D*)inclz[ii]->Clone(Form("sub_inclz_pt%d",ii));
		sub_inclz[ii]->Add(perp_inclz[ii], -1);

	}//ii

	TCanvas *c1 = new TCanvas("c1","c1",600,700);
	TPad *p1_1 = new TPad("p1_1","p1_1",0,0.4,1,1);
	p1_1->Draw();
	p1_1->cd();
	SetPadStyle();
	gPad->SetBottomMargin(0);
	gPad->SetLeftMargin(0.13);
	p1_1->SetLogy(1);
	p1_1->SetLogx(1);

	htmp = (TH1D*)gPad->DrawFrame(0.09, 8e-4, 2.5, 5e1);
	SetHistoStyle("","1/N_{jet} dN/d#it{j}_{T} (GeV/c)^{-1}","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.5);

	sub_incljt[0]->SetMarkerStyle(24);
	sub_incljt[0]->SetMarkerColor(2);
	sub_incljt[0]->SetLineColor(2);
	sub_incljt[0]->Draw("p same");

	/*
	c1->cd();
	TPad *p1_2 = new TPad("p1_2","p1_2",0,0.0,1,0.4);
	p1_2->Draw();
	p1_2->cd();
	SetPadStyle();
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.13);
	p1_2->SetLogx(1);

	htmp = (TH1D*)gPad->DrawFrame(0.09, 0.2, 2.5, 1.8);
	SetHistoStyle("#it{j}_{T} (GeV/c)","Ratio","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.5);
	htmp->GetYaxis()->SetNdivisions(7,5,0);
	htmp->GetYaxis()->CenterTitle();
	htmp->GetXaxis()->SetTitleOffset(3.0);

	TH1D *ratiojt = (TH1D*)sub_incljt->Clone("ratiojt");
	ratiojt->Divide(sub_sigjt);
	ratiojt->Draw("p same");

	TCanvas *c2 = new TCanvas("c2","c2",600,700);
	TPad *p2_1 = new TPad("p2_1","p2_1",0,0.4,1,1);
	p2_1->Draw();
	p2_1->cd();
	SetPadStyle();
	gPad->SetBottomMargin(0);
	gPad->SetLeftMargin(0.13);
	p2_1->SetLogy(1);
	p2_1->SetLogx(1);

	htmp = (TH1D*)gPad->DrawFrame(0.02, 5e-3, 1, 1e2);
	SetHistoStyle("","1/N_{jet} dN/dz","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.5);

	sub_sigz->SetMarkerStyle(20);
	sub_sigz->SetMarkerColor(1);
	sub_sigz->SetLineColor(1);
	sub_sigz->Draw("p same");

	sub_inclz->SetMarkerStyle(24);
	sub_inclz->SetMarkerColor(2);
	sub_inclz->SetLineColor(2);
	sub_inclz->Draw("p same");

	c2->cd();
	TPad *p2_2 = new TPad("p2_2","p2_2",0,0.0,1,0.4);
	p2_2->Draw();
	p2_2->cd();
	SetPadStyle();
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.13);
	p2_2->SetLogx(1);

	htmp = (TH1D*)gPad->DrawFrame(0.02, 0.2, 1, 1.8);
	SetHistoStyle("z","Ratio","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.5);
	htmp->GetYaxis()->SetNdivisions(7,5,0);
	htmp->GetYaxis()->CenterTitle();
	htmp->GetXaxis()->SetTitleOffset(3.0);

	TH1D *ratioz = (TH1D*)sub_inclz->Clone("ratioz");
	ratioz->Divide(sub_sigz);
	ratioz->Draw("p same");
	*/

	return;

	/*
	//delete infile;
	TFile *outfile = new TFile("subjt.root", "RECREATE");
	sigjt->Write();
	incljt->Write();
	subjt->Write();
	outfile-> Close();
	*/
}

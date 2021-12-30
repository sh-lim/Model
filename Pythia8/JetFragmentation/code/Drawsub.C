#include "Style.h"

void Drawsub(){

	TFile *infile = new TFile("jtresult.root","read");
	TH1D *jeteta = (TH1D*)infile->Get("jeteta");
	TH1D *perpeta = (TH1D*)infile->Get("perpeta");

	int etabin_lo = jeteta->FindBin(-0.5+0.001);
	int etabin_hi = jeteta->FindBin(+0.5-0.001);

	float njet = jeteta->Integral(etabin_lo, etabin_hi);
	float nperpjet = perpeta->Integral(etabin_lo,etabin_hi);

	TH1D *sigjt = (TH1D*)infile->Get("sigjt");
	TH1D *incljt = (TH1D*)infile->Get("incljt");

	TH1D *sigz = (TH1D*)infile->Get("sigz");
	TH1D *inclz = (TH1D*)infile->Get("inclz");

	TH1D *perp_sigjt = (TH1D*)infile->Get("perp_sigjt");
	TH1D *perp_incljt = (TH1D*)infile->Get("perp_incljt");
	perp_sigjt->Sumw2();
	perp_incljt->Sumw2();

	TH1D *perp_sigz = (TH1D*)infile->Get("perp_sigz");
	TH1D *perp_inclz = (TH1D*)infile->Get("perp_inclz");
	perp_sigz->Sumw2();
	perp_inclz->Sumw2();

	for (int ijt=0; ijt<sigjt->GetNbinsX(); ijt++){
		float djt = sigjt->GetBinWidth(ijt+1);
		float xx = sigjt->GetBinContent(ijt+1);
		float xx_err = sigjt->GetBinError(ijt+1);

		sigjt->SetBinContent(ijt+1, xx/djt/njet);
		sigjt->SetBinError(ijt+1, xx_err/djt/njet);

		xx = incljt->GetBinContent(ijt+1);
		xx_err = incljt->GetBinError(ijt+1);

		incljt->SetBinContent(ijt+1, xx/djt/njet);
		incljt->SetBinError(ijt+1, xx_err/djt/njet);

		xx = perp_sigjt->GetBinContent(ijt+1);
		xx_err = perp_sigjt->GetBinError(ijt+1);

		perp_sigjt->SetBinContent(ijt+1, xx/djt/nperpjet);
		perp_sigjt->SetBinError(ijt+1, xx_err/djt/nperpjet);

		xx = perp_incljt->GetBinContent(ijt+1);
		xx_err = perp_incljt->GetBinError(ijt+1);

		perp_incljt->SetBinContent(ijt+1, xx/djt/nperpjet);
		perp_incljt->SetBinError(ijt+1, xx_err/djt/nperpjet);
	}

	for (int iz=0; iz<sigz->GetNbinsX(); iz++){
		float dz = sigz->GetBinWidth(iz+1);
		float xx = sigz->GetBinContent(iz+1);
		float xx_err = sigz->GetBinError(iz+1);

		sigz->SetBinContent(iz+1, xx/dz/njet);
		sigz->SetBinError(iz+1, xx_err/dz/njet);

		xx = inclz->GetBinContent(iz+1);
		xx_err = inclz->GetBinError(iz+1);

		inclz->SetBinContent(iz+1, xx/dz/njet);
		inclz->SetBinError(iz+1, xx_err/dz/njet);

		xx = perp_sigz->GetBinContent(iz+1);
		xx_err = perp_sigz->GetBinError(iz+1);

		perp_sigz->SetBinContent(iz+1, xx/dz/nperpjet);
		perp_sigz->SetBinError(iz+1, xx_err/dz/nperpjet);

		xx = perp_inclz->GetBinContent(iz+1);
		xx_err = perp_inclz->GetBinError(iz+1);

		perp_inclz->SetBinContent(iz+1, xx/dz/nperpjet);
		perp_inclz->SetBinError(iz+1, xx_err/dz/nperpjet);
	}

	//return;

	auto sub_incljt = (TH1D*)incljt->Clone("sub_incljt");
	sub_incljt->Add(perp_incljt, -1);

	auto sub_sigjt = (TH1D*)sigjt->Clone("sub_sigjt");
	sub_sigjt->Add(perp_sigjt, -1);

	auto sub_inclz = (TH1D*)inclz->Clone("sub_inclz");
	sub_inclz->Add(perp_inclz, -1);

	auto sub_sigz = (TH1D*)sigz->Clone("sub_sigz");
	sub_sigz->Add(perp_sigz, -1);

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

	sub_sigjt->SetMarkerStyle(20);
	sub_sigjt->SetMarkerColor(1);
	sub_sigjt->SetLineColor(1);
	sub_sigjt->Draw("p same");

	sub_incljt->SetMarkerStyle(24);
	sub_incljt->SetMarkerColor(2);
	sub_incljt->SetLineColor(2);
	sub_incljt->Draw("p same");

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

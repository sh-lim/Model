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

	TH1D *perp_sigjt = (TH1D*)infile->Get("perp_sigjt");
	TH1D *perp_incljt = (TH1D*)infile->Get("perp_incljt");
	perp_sigjt->Sumw2();
	perp_incljt->Sumw2();

	for (int sjt=0; sjt<sigjt->GetNbinsX(); sjt++){
		float djt = sigjt->GetBinWidth(sjt+1);
		float xx = sigjt->GetBinContent(sjt+1);
		float xx_err = sigjt->GetBinError(sjt+1);

		sigjt->SetBinContent(sjt+1, xx/djt/njet);
		sigjt->SetBinError(sjt+1, xx_err/djt/njet);
	}


	for (int ijt=0; ijt<incljt->GetNbinsX(); ijt++){
		float djt = incljt->GetBinWidth(ijt+1);
		float xx = incljt->GetBinContent(ijt+1);
		float xx_err = incljt->GetBinError(ijt+1);

		incljt->SetBinContent(ijt+1, xx/djt/njet);
		incljt->SetBinError(ijt+1, xx_err/djt/njet);
	}

	for (int pjt=0; pjt<perp_sigjt->GetNbinsX(); pjt++){
		float dpjt = perp_sigjt->GetBinWidth(pjt+1);
		float pxx = perp_sigjt->GetBinContent(pjt+1);
		float pxx_err = perp_sigjt->GetBinError(pjt+1);

		perp_sigjt->SetBinContent(pjt+1, pxx/dpjt/nperpjet);
		perp_sigjt->SetBinError(pjt+1, pxx_err/dpjt/nperpjet);
	}

	for (int pjt=0; pjt<perp_incljt->GetNbinsX(); pjt++){
		float dpjt = perp_incljt->GetBinWidth(pjt+1);
		float pxx = perp_incljt->GetBinContent(pjt+1);
		float pxx_err = perp_incljt->GetBinError(pjt+1);

		perp_incljt->SetBinContent(pjt+1, pxx/dpjt/nperpjet);
		perp_incljt->SetBinError(pjt+1, pxx_err/dpjt/nperpjet);
	}

	auto sub_incljt = (TH1D*)incljt->Clone("sub_incljt");
	sub_incljt->Add(perp_incljt, -1);

	auto sub_sigjt = (TH1D*)sigjt->Clone("sub_sigjt");
	sub_sigjt->Add(perp_sigjt, -1);

	TCanvas *c1 = new TCanvas("c1","c1",600,700);
	TPad *p1_1 = new TPad("p1_1","p1_1",0,0.4,1,1);
	p1_1->Draw();
	p1_1->cd();
	SetPadStyle();
	gPad->SetBottomMargin(0);
	gPad->SetLeftMargin(0.13);
	p1_1->SetLogy(1);
	p1_1->SetLogx(1);

	htmp = (TH1D*)gPad->DrawFrame(0.07, 5e-3, 2, 1e2);
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

	TH1D *ratiojt = (TH1D*)sub_sigjt->Clone("ratiojt");
	ratiojt->Divide(sub_incljt);

	c1->cd();
	TPad *p1_2 = new TPad("p1_2","p1_2",0,0.0,1,0.4);
	p1_2->Draw();
	p1_2->cd();
	SetPadStyle();
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.13);
	p1_2->SetLogx(1);

	htmp = (TH1D*)gPad->DrawFrame(0.07, 0, 2, 2);
	SetHistoStyle("#it{j}_{T} (GeV/c)","Ratio","",20,18);
	htmp->GetYaxis()->SetTitleOffset(1.5);
	htmp->GetXaxis()->SetTitleOffset(3.0);

	ratiojt->Draw("p same");



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

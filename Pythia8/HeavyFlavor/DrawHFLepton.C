#include "Style.h"

void DrawHFLepton(){

	const int nMarker[5] = {24, 25, 26, 20, 21};
	const int nColor[5] = {1, 2, 4, 6, 8};

	const int ppXection = 42.0;
	const float nEvent = 100*10000000;
	const float deta = 0.7;

	const float twopi = TMath::TwoPi();

	const int nset = 1;

	TFile *infile0 = new TFile("ppg223_refoldelectrons.root","read");
	TGraphAsymmErrors *gPHENIXcc = (TGraphAsymmErrors*)infile0->Get("cross_charm");
	TGraphAsymmErrors *gPHENIXbb = (TGraphAsymmErrors*)infile0->Get("cross_bottom");

	TFile *infile[nset];
	infile[0] = new TFile("outfile.root","read");

	TH1D *hHFemu_cc[nset];
	TH1D *hHFemu_bb[nset];

	for (int iset=0; iset<nset; iset++){

		hHFemu_cc[iset] = (TH1D*)infile[iset]->Get("hHFemu_cc");
		hHFemu_bb[iset] = (TH1D*)infile[iset]->Get("hHFemu_bb");

		hHFemu_cc[iset]->Rebin();
		hHFemu_bb[iset]->Rebin();

		for (int ii=0; ii<hHFemu_cc[iset]->GetNbinsX(); ii++){

			float pt = hHFemu_cc[iset]->GetBinCenter(ii+1);
			float dpt = hHFemu_cc[iset]->GetBinWidth(ii+1);

			float yy = hHFemu_cc[iset]->GetBinContent(ii+1);
			float yy_err = hHFemu_cc[iset]->GetBinError(ii+1);

			float yield = yy*ppXection/twopi/pt/dpt/deta/nEvent/4.0;
			float yield_err = yy_err*ppXection/twopi/pt/dpt/deta/nEvent/4.0;

			hHFemu_cc[iset]->SetBinContent(ii+1, yield);
			hHFemu_cc[iset]->SetBinError(ii+1, yield_err);

			yy = hHFemu_bb[iset]->GetBinContent(ii+1);
			yy_err = hHFemu_bb[iset]->GetBinError(ii+1);

			yield = yy*ppXection/twopi/pt/dpt/deta/nEvent/4.0;
			yield_err = yy_err*ppXection/twopi/pt/dpt/deta/nEvent/4.0;

			hHFemu_bb[iset]->SetBinContent(ii+1, yield);
			hHFemu_bb[iset]->SetBinError(ii+1, yield_err);

		}//ii
	}//iset

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*500,500);
	c1->Divide(2,1);

	{
		c1->cd(1);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetLeftMargin(0.18);

		htmp = (TH1D*)gPad->DrawFrame(0,1e-10,10,1e-2);
		SetHistoStyle("p_{T} [GeV/c]","#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} [mb (GeV/c)^{-2}]","",22,20);
		htmp->GetYaxis()->SetTitleOffset(1.7);

		for (int iset=0; iset<nset; iset++){
			hHFemu_cc[iset]->SetLineColor(nColor[iset+1]);
			hHFemu_cc[iset]->SetMarkerColor(nColor[iset+1]);
			hHFemu_cc[iset]->SetMarkerStyle(20);
			hHFemu_cc[iset]->SetMarkerSize(0.9);
			hHFemu_cc[iset]->Draw("p same");
		}

		gPHENIXcc->SetMarkerStyle(24);
		gPHENIXcc->SetLineColor(1);
		gPHENIXcc->SetMarkerColor(1);
		gPHENIXcc->Draw("p");
	}

	{
		c1->cd(2);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetLeftMargin(0.18);

		htmp = (TH1D*)gPad->DrawFrame(0,1e-10,10,1e-2);
		SetHistoStyle("p_{T} [GeV/c]","#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} [mb (GeV/c)^{-2}]","",22,20);
		htmp->GetYaxis()->SetTitleOffset(1.7);

		for (int iset=0; iset<nset; iset++){
			hHFemu_bb[iset]->SetLineColor(nColor[iset+1]);
			hHFemu_bb[iset]->SetMarkerColor(nColor[iset+1]);
			hHFemu_bb[iset]->SetMarkerStyle(20);
			hHFemu_bb[iset]->SetMarkerSize(0.9);
			hHFemu_bb[iset]->Draw("p same");
		}

		gPHENIXbb->SetMarkerStyle(24);
		gPHENIXbb->SetLineColor(1);
		gPHENIXbb->SetMarkerColor(1);
		gPHENIXbb->Draw("p");
	}


}

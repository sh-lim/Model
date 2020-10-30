#include "Style.h"

void Draw_frag2(){

	TFile *infile = new TFile("outfile_hist.root","read");

	TH2D *hjetR04_eta_pt = (TH2D*)infile->Get("hjetR04_eta_pt");

	TH1D *hjetR04_eta = (TH1D*)hjetR04_eta_pt->ProjectionX("hjetR04_eta");
	int etabin_lo = hjetR04_eta->FindBin(-0.5+0.001);
	int etabin_hi = hjetR04_eta->FindBin(+0.5-0.001);
	TH1D *hjetR04_pt = (TH1D*)hjetR04_eta_pt->ProjectionY("hjetR04_pt",etabin_lo,etabin_hi);

	{
		TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,1*400);
		c1->Divide(2,1);

		c1->cd(1);
		SetPadStyle();
		gPad->SetRightMargin(0.15);

		htmp = (TH1D*)gPad->DrawFrame(-1.5,0,1.5,1.2*hjetR04_eta->GetMaximum());
		SetHistoStyle("Jet #eta","N","",20,16);
		hjetR04_eta->SetLineColor(1);
		hjetR04_eta->SetLineWidth(2);
		hjetR04_eta->Draw("same");

		c1->cd(2);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(0,1,200,1.5*hjetR04_pt->GetMaximum());
		SetHistoStyle("Jet p_{T} (GeV/c)","N","",20,16);
		hjetR04_pt->SetLineColor(1);
		hjetR04_pt->SetLineWidth(2);
		hjetR04_pt->Draw("same");
	}

	const int nptbin = 4;
	const int jetpt[nptbin+1] = {25, 40, 60, 80, 110};

	const int nColor[nptbin] = {1, 2, 4, 6};
	const int nMarker[nptbin] = {20, 21, 24, 25};

	TFile *infileZ[nptbin];
	TFile *infileJT[nptbin];

	TGraphAsymmErrors *gZ[nptbin];
	TGraphAsymmErrors *gJT[nptbin];

	for (int iset=0; iset<nptbin; iset++){
		infileZ[iset] = new TFile(Form("/Users/shlim/Work/PNU-NPL/Research/Pythia8/Jet/HEPData-ins929691-v1-Table_%d.root",iset+1),"read");
		infileJT[iset] = new TFile(Form("/Users/shlim/Work/PNU-NPL/Research/Pythia8/Jet/HEPData-ins929691-v1-Table_%d.root",iset+21),"read");
		TDirectoryFile *tdf = (TDirectoryFile*)infileZ[iset]->Get(Form("Table %d",iset+1));
		gZ[iset] = (TGraphAsymmErrors*)tdf->Get("Graph1D_y1");

		tdf = (TDirectoryFile*)infileJT[iset]->Get(Form("Table %d",iset+21));
		gJT[iset] = (TGraphAsymmErrors*)tdf->Get("Graph1D_y1");
	}

	TH2D *hjetR04_jt_pt = (TH2D*)infile->Get("hjetR04_jt_pt");
	hjetR04_jt_pt->Sumw2();

	TH2D *hjetR04_zh_pt = (TH2D*)infile->Get("hjetR04_zh_pt");
	hjetR04_zh_pt->Sumw2();

	TH1D *hjetR04_jt[nptbin];
	TH1D *hjetR04_zh[nptbin];

	for (int ipt=0; ipt<nptbin; ipt++){

		int ptbin_lo = hjetR04_jt_pt->GetYaxis()->FindBin(jetpt[ipt]+0.1);
		int ptbin_hi = hjetR04_jt_pt->GetYaxis()->FindBin(jetpt[ipt+1]-0.1);
		hjetR04_jt[ipt] = (TH1D*)hjetR04_jt_pt->ProjectionX(Form("hjetR04_jt_pt%d",ipt),ptbin_lo,ptbin_hi);
		hjetR04_zh[ipt] = (TH1D*)hjetR04_zh_pt->ProjectionX(Form("hjetR04_zh_pt%d",ipt),ptbin_lo,ptbin_hi);

		ptbin_lo = hjetR04_pt->FindBin(jetpt[ipt]+0.1);
		ptbin_hi = hjetR04_pt->FindBin(jetpt[ipt+1]-0.1);
		float njet_eta05 = hjetR04_pt->Integral(ptbin_lo, ptbin_hi);

		for (int ijt=0; ijt<hjetR04_jt[ipt]->GetNbinsX(); ijt++){
			float xx = hjetR04_jt[ipt]->GetBinContent(ijt+1);
			float xx_err = hjetR04_jt[ipt]->GetBinError(ijt+1);
			float djt = hjetR04_jt[ipt]->GetBinWidth(ijt+1);


			hjetR04_jt[ipt]->SetBinContent(ijt+1, xx/djt/njet_eta05);
			hjetR04_jt[ipt]->SetBinError(ijt+1, xx_err/djt/njet_eta05);
		}
		hjetR04_jt[ipt]->SetLineWidth(2);
		hjetR04_jt[ipt]->SetLineColor(nColor[ipt]);

		for (int izh=0; izh<hjetR04_zh[ipt]->GetNbinsX(); izh++){
			float xx = hjetR04_zh[ipt]->GetBinContent(izh+1);
			float xx_err = hjetR04_zh[ipt]->GetBinError(izh+1);
			float dzh = hjetR04_zh[ipt]->GetBinWidth(izh+1);

			hjetR04_zh[ipt]->SetBinContent(izh+1, xx/dzh/njet_eta05);
			hjetR04_zh[ipt]->SetBinError(izh+1, xx_err/dzh/njet_eta05);

		}

		hjetR04_zh[ipt]->SetLineWidth(2);
		hjetR04_zh[ipt]->SetLineColor(nColor[ipt]);
	}

	{

		TCanvas *c2 = new TCanvas("c2","c2",1.1*2*400,1*400);
		c2->Divide(2,1);

		c2->cd(1);
		SetPadStyle();
		//gPad->SetLogy();
		//gPad->SetLogx();

		htmp = (TH1D*)gPad->DrawFrame(0.0,0.0,3,20);
		SetHistoStyle("p_{T}^{rel} (GeV/c)","f(p_{T}^{rel})","",20,16);

		hjetR04_jt[0]->Draw("same");
		hjetR04_jt[1]->Draw("same");
		hjetR04_jt[2]->Draw("same");
		hjetR04_jt[3]->Draw("same");

		for (int iset=0; iset<nptbin; iset++){
			gJT[iset]->SetMarkerStyle(nMarker[iset]);
			gJT[iset]->SetMarkerColor(nColor[iset]);
			gJT[iset]->SetLineColor(nColor[iset]);
			gJT[iset]->Draw("P");
		}

		c2->cd(2);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetLogx();

		htmp = (TH1D*)gPad->DrawFrame(0.01,0.02,1,1e3);
		SetHistoStyle("z","F(z)","",20,16);

		hjetR04_zh[0]->Draw("same");
		hjetR04_zh[1]->Draw("same");
		hjetR04_zh[2]->Draw("same");
		hjetR04_zh[3]->Draw("same");

		for (int iset=0; iset<nptbin; iset++){
			gZ[iset]->SetMarkerStyle(nMarker[iset]);
			gZ[iset]->SetMarkerColor(nColor[iset]);
			gZ[iset]->SetLineColor(nColor[iset]);
			gZ[iset]->Draw("P");
		}

	}


}

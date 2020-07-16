#include "Style.h"

void Draw_frag(){

	TFile *infile = new TFile("outfile_hist.root","read");

	TH2D *hjetR04_eta_pt = (TH2D*)infile->Get("hjetR04_eta_pt");
	//hjetR04_eta_pt->Sumw2();
	TH1D *hjetR04_eta = (TH1D*)hjetR04_eta_pt->ProjectionX("hjetR04_eta");
	TH1D *hjetR04_pt = (TH1D*)hjetR04_eta_pt->ProjectionY("hjetR04_pt");

	int etabin_lo = hjetR04_eta->FindBin(-0.5+0.001);
	int etabin_hi = hjetR04_eta->FindBin(+0.5-0.001);
	float njet_eta05 = hjetR04_eta->Integral(etabin_lo, etabin_hi);

	TH2D *hjetR04_jt_pt = (TH2D*)infile->Get("hjetR04_jt_pt");
	hjetR04_jt_pt->Sumw2();
	TH1D *hjetR04_jt = (TH1D*)hjetR04_jt_pt->ProjectionX("hjetR04_jt");
	for (int ijt=0; ijt<hjetR04_jt->GetNbinsX(); ijt++){
		float djt = hjetR04_jt->GetBinWidth(ijt+1);
		float jt = hjetR04_jt->GetBinCenter(ijt+1);

		float xx = hjetR04_jt->GetBinContent(ijt+1);
		float xx_err = hjetR04_jt->GetBinError(ijt+1);


		hjetR04_jt->SetBinContent(ijt+1, xx/djt/njet_eta05);
		hjetR04_jt->SetBinError(ijt+1, xx_err/djt/njet_eta05);
	}

	TH2D *hjetR04_zh_pt = (TH2D*)infile->Get("hjetR04_zh_pt");
	hjetR04_zh_pt->Sumw2();
	TH1D *hjetR04_zh = (TH1D*)hjetR04_zh_pt->ProjectionX("hjetR04_zh");
	for (int izh=0; izh<hjetR04_zh->GetNbinsX(); izh++){
		float dzh = hjetR04_zh->GetBinWidth(izh+1);
		float zh = hjetR04_zh->GetBinCenter(izh+1);

		float xx = hjetR04_zh->GetBinContent(izh+1);
		float xx_err = hjetR04_zh->GetBinError(izh+1);


		hjetR04_zh->SetBinContent(izh+1, xx/dzh/njet_eta05);
		hjetR04_zh->SetBinError(izh+1, xx_err/dzh/njet_eta05);
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.2*3*300,2*300);
	c1->Divide(3,2);

	c1->cd(1);
	SetPadStyle();
	gPad->SetRightMargin(0.15);

	htmp = (TH1D*)gPad->DrawFrame(-1.5,0,1.5,200);
	SetHistoStyle("Jet #eta","Jet p_{T} (GeV/c)","",16,16);

	hjetR04_eta_pt->Draw("colz same");

	c1->cd(2);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-1.5,0,1.5,1.2*hjetR04_eta->GetMaximum());
	SetHistoStyle("Jet #eta","N","",16,16);
	hjetR04_eta->SetLineColor(1);
	hjetR04_eta->SetLineWidth(2);
	hjetR04_eta->Draw("same");

	c1->cd(3);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1D*)gPad->DrawFrame(0,1,200,1.5*hjetR04_pt->GetMaximum());
	SetHistoStyle("Jet p_{T} (GeV/c)","N","",16,16);
	hjetR04_pt->SetLineColor(1);
	hjetR04_pt->SetLineWidth(2);
	hjetR04_pt->Draw("same");

	c1->cd(5);
	SetPadStyle();
	//gPad->SetLogy();
	//gPad->SetLogx();

	//htmp = (TH1D*)gPad->DrawFrame(0.1,1e-3,3,1.5*hjetR04_jt->GetMaximum());
	htmp = (TH1D*)gPad->DrawFrame(0.0,0.0,3,1.2*hjetR04_jt->GetMaximum());
	SetHistoStyle("p_{T}^{rel} (GeV/c)","f(p_{T}^{rel})","",16,16);

	hjetR04_jt->SetLineColor(1);
	hjetR04_jt->SetLineWidth(2);
	hjetR04_jt->Draw("same");

	c1->cd(6);
	SetPadStyle();
	gPad->SetLogy();
	gPad->SetLogx();

	htmp = (TH1D*)gPad->DrawFrame(0.005,0.02,0.8,1.5*hjetR04_zh->GetMaximum());
	SetHistoStyle("z","F(z)","",16,16);

	hjetR04_zh->SetLineColor(1);
	hjetR04_zh->SetLineWidth(2);
	hjetR04_zh->Draw("same");
}

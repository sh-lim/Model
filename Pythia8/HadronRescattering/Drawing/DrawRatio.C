#include "Style.h"

void DrawRatio(){

	gStyle -> SetOptStat(0);

	const int had = 2;
	const int nfile = 2;
	const int nmult = 6;
	const int npt = 15;

	TFile * infile[had][nfile];

	infile[0][0] = new TFile("outfile_kaon_on_0_changeptmult.root", "read");
	infile[0][1] = new TFile("outfile_kaon_on_1_changeptmult.root", "read");
	infile[1][0] = new TFile("outfile_kaon_off_0_changeptmult.root", "read");
	infile[1][1] = new TFile("outfile_kaon_off_1_changeptmult.root", "read");

	TH1D * hevent[had][nfile];
	TH1D * hkaon_pt[had][nfile][nmult];
	TH1D * hpion_pt[had][nfile][nmult];

	TH1D * hkaon_ratio[nmult];
	TH1D * hpion_ratio[nmult];

	for( int j = 0 ; j < had ; j++ ){
		for( int k = 0 ; k < nfile ; k++ ){
			for( int i = 0 ; i < nmult ; i++ ){

				hpion_pt[j][k][i] = ( TH1D * ) infile[j][k] -> Get(Form("histo_pion_pT_m%d", i));
				hkaon_pt[j][k][i] = ( TH1D * ) infile[j][k] -> Get(Form("histo_kaon_pT_m%d", i));

				hpion_pt[j][k][i]->Sumw2();
				hkaon_pt[j][k][i]->Sumw2();

				if ( k>0 ){
					hpion_pt[j][0][i] -> Add(hpion_pt[j][k][i], 1);
					hkaon_pt[j][0][i] -> Add(hkaon_pt[j][k][i], 1);
				}
			}//i

			hevent[j][k] = (TH1D*)infile[j][k]->Get("hist_event");

			if ( k>0 ){
				hevent[j][0]->Add(hevent[j][k]);
			}
		}//k
	}//j

	TGraphErrors *gkaon_ratio[nmult];
	TGraphErrors *gpion_ratio[nmult];

	for( int i = 0 ; i < nmult ; i++ ){

		hkaon_ratio[i] = (TH1D*)hkaon_pt[0][0][i]->Clone(Form("hkaon_ratio_%d",i));
		hpion_ratio[i] = (TH1D*)hpion_pt[0][0][i]->Clone(Form("hpion_ratio_%d",i));

		hkaon_ratio[i]->Divide(hkaon_pt[1][0][i]);
		hpion_ratio[i]->Divide(hpion_pt[1][0][i]);

		float norm = hevent[1][0]->GetBinContent(i+1)/hevent[0][0]->GetBinContent(i+1);

		hkaon_ratio[i]->Scale(norm);
		hpion_ratio[i]->Scale(norm);

		gkaon_ratio[i] = new TGraphErrors;
		gpion_ratio[i] = new TGraphErrors;

		for (int ipt=0; ipt<npt; ipt++){

			float xx = hkaon_ratio[i]->GetBinCenter(ipt+1);
			float yy = hkaon_ratio[i]->GetBinContent(ipt+1);
			float yy_err = hkaon_ratio[i]->GetBinError(ipt+1);

			gkaon_ratio[i]->SetPoint(ipt, xx, yy);
			gkaon_ratio[i]->SetPointError(ipt, xx, yy_err);

			yy = hpion_ratio[i]->GetBinContent(ipt+1);
			yy_err = hpion_ratio[i]->GetBinError(ipt+1);
			gpion_ratio[i]->SetPoint(ipt, xx, yy);
			gpion_ratio[i]->SetPointError(ipt, xx, yy_err);

		}//ipt

	}//i

	TCanvas *c0 = new TCanvas("c0","c0",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0.7,6,1.3);
	SetHistoStyle("p_{T} (GeV/c)","Ratio","",24,22);

	gpion_ratio[0]->SetLineColorAlpha(1,0.3);
	gpion_ratio[0]->SetLineWidth(4);
	gpion_ratio[0]->SetFillColorAlpha(1,0.3);
	gpion_ratio[0]->Draw("L3");

	gpion_ratio[5]->SetLineColorAlpha(2,0.3);
	gpion_ratio[5]->SetLineWidth(4);
	gpion_ratio[5]->SetFillColorAlpha(2,0.3);
	gpion_ratio[5]->Draw("L3");

	{
		TLegend *leg = new TLegend(0.4,0.9-0.05*5,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","Charged pion, |y|<1","h");
		leg->AddEntry("","Hadron rescattering on/off","h");
		leg->AddEntry(gpion_ratio[5],"0-5%","LF");
		leg->AddEntry(gpion_ratio[0],"60-100%","LF");
		leg->Draw();
	}


	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0.7,6,1.3);
	SetHistoStyle("p_{T} (GeV/c)","Ratio","",24,22);

	gkaon_ratio[0]->SetLineColorAlpha(1,0.3);
	gkaon_ratio[0]->SetLineWidth(4);
	gkaon_ratio[0]->SetFillColorAlpha(1,0.3);
	gkaon_ratio[0]->Draw("L3");

	gkaon_ratio[5]->SetLineColorAlpha(2,0.3);
	gkaon_ratio[5]->SetLineWidth(4);
	gkaon_ratio[5]->SetFillColorAlpha(2,0.3);
	gkaon_ratio[5]->Draw("L3");

	{
		TLegend *leg = new TLegend(0.4,0.9-0.05*5,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","Charged kaon, |y|<1","h");
		leg->AddEntry("","Hadron rescattering on/off","h");
		leg->AddEntry(gkaon_ratio[5],"0-5%","LF");
		leg->AddEntry(gkaon_ratio[0],"60-100%","LF");
		leg->Draw();
	}

}

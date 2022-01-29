#include "Style.h"

void ScanV0M_70(){

	const int ngrp = 2;

	TFile *infile[ngrp];

	TH1D *href[ngrp];
	TH2D *hLc[ngrp];
	TH2D *hD0[ngrp];

	for (int ig=0; ig<ngrp; ig++){

		infile[ig] = new TFile(Form("outfile_hist_set70_grp%03d.root",ig),"read");

		href[ig] = (TH1D*)infile[ig]->Get("hmult_v0m");
		hLc[ig] = (TH2D*)infile[ig]->Get("hLc_pt_v0m");
		hD0[ig] = (TH2D*)infile[ig]->Get("hD0_pt_v0m");

	}//ig

	const int nmult = 4;
	const int mult_dn[nmult] = {0, 0, 92, 242}; 
	const int mult_up[nmult] = {500, 92, 242, 500}; 

	const int nMarker[nmult] = {24, 25, 26, 27}; 
	const int nColor[nmult] = {1, 2, 4, 6}; 

	TH1D *hLc_pt[ngrp][nmult];
	TH1D *hD0_pt[ngrp][nmult];

	TH1D *hLc_D0_ratio[ngrp][nmult];

	float nevent[ngrp][nmult] = {0};

	for (int ig=0; ig<ngrp; ig++){
		for (int ii=0; ii<nmult; ii++){

			nevent[ig][ii] = href[ig]->Integral(mult_dn[ii]+1, mult_up[ii]);

			cout << ig << " " << ii << " " << nevent[ig][ii] << endl;

			hLc_pt[ig][ii] = hLc[ig]->ProjectionX(Form("hLc_pt_grp%d_mult%d",ig,ii),mult_dn[ii]+1,mult_up[ii]); 
			hLc_pt[ig][ii]->Rebin(4);
			hLc_pt[ig][ii]->Sumw2();

			float norm = hLc_pt[ig][ii]->GetBinWidth(1)*nevent[ig][ii];
			hLc_pt[ig][ii]->Scale(1./norm);

			hD0_pt[ig][ii] = hD0[ig]->ProjectionX(Form("hD0_pt_grp%d_mult%d",ig,ii),mult_dn[ii]+1,mult_up[ii]); 
			hD0_pt[ig][ii]->Rebin(4);
			hD0_pt[ig][ii]->Sumw2();

			norm = hD0_pt[ig][ii]->GetBinWidth(1)*nevent[ig][ii];
			hD0_pt[ig][ii]->Scale(1./norm);

			hLc_D0_ratio[ig][ii] = (TH1D*)hLc_pt[ig][ii]->Clone(Form("hLc_D0_ratio_grp%d_mult%d",ig,ii));
			hLc_D0_ratio[ig][ii]->Divide(hD0_pt[ig][ii]);

			hLc_D0_ratio[ig][ii]->SetMarkerStyle(nMarker[ii]);
			hLc_D0_ratio[ig][ii]->SetMarkerColor(nColor[ii]);
			hLc_D0_ratio[ig][ii]->SetLineColor(nColor[ii]);

		}//ii
	}//ig

	//return;

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*500,500);
	c1->Divide(2,1);

	{
		c1->cd(1);
		SetPadStyle();
		gPad->SetTicks();
		htmp = (TH1D*)gPad->DrawFrame(0, 0, 20, 1.0);
		SetHistoStyle("p_{T} (GeV/c)","#Lambda_{c}^{+}/D^{0}","");

		hLc_D0_ratio[0][0]->Draw("p same");
		hLc_D0_ratio[0][1]->Draw("p same");
		hLc_D0_ratio[0][2]->Draw("p same");
		hLc_D0_ratio[0][3]->Draw("p same");

		TLegend *leg = new TLegend(0.6,0.90-0.05*6,0.9,0.90);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","Monash","h");
		leg->AddEntry(hLc_D0_ratio[0][0],"0-100%","P");
		leg->AddEntry(hLc_D0_ratio[0][1],"30-100%","P");
		leg->AddEntry(hLc_D0_ratio[0][2],"0.1-30%","P");
		leg->AddEntry(hLc_D0_ratio[0][3],"0-0.1%","P");
		leg->Draw();
	}

	{
		c1->cd(2);
		SetPadStyle();
		gPad->SetTicks();
		htmp = (TH1D*)gPad->DrawFrame(0, 0, 20, 1.0);
		SetHistoStyle("p_{T} (GeV/c)","#Lambda_{c}^{+}/D^{0}","");

		hLc_D0_ratio[1][0]->Draw("p same");
		hLc_D0_ratio[1][1]->Draw("p same");
		hLc_D0_ratio[1][2]->Draw("p same");
		hLc_D0_ratio[1][3]->Draw("p same");

		TLegend *leg = new TLegend(0.6,0.90-0.05*6,0.9,0.90);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","CR mode 2","h");
		leg->AddEntry(hLc_D0_ratio[0][0],"0-100%","P");
		leg->AddEntry(hLc_D0_ratio[0][1],"30-100%","P");
		leg->AddEntry(hLc_D0_ratio[0][2],"0.1-30%","P");
		leg->AddEntry(hLc_D0_ratio[0][3],"0-0.1%","P");
		leg->Draw();
	}


	/*
	for (int ii=0; ii<href->GetNbinsX(); ii++){
		cout << "0 - " << ii+1 << " " << href->Integral(1,ii+1)/href->Integral() << endl;
	}
	*/

}

void ScanV0M(){

	const int ngrp = 2;

	TFile *infile[ngrp];

	TH1D *href[ngrp];
	TH2D *hXib[ngrp];

	for (int ig=0; ig<ngrp; ig++){

		infile[ig] = new TFile(Form("outfile_hist_set60_grp%03d.root",ig),"read");

		href[ig] = (TH1D*)infile[ig]->Get("hmult_v0m");
		hXib[ig] = (TH2D*)infile[ig]->Get("hXib_pt_v0m");

		if ( ig>0 ){
			href[0]->Add(href[ig]);
			hXib[0]->Add(hXib[ig]);
		}//igrp
	}

	const int nmult = 4;
	const int mult_dn[nmult] = {0, 0, 92, 242}; 
	const int mult_up[nmult] = {500, 92, 242, 500}; 

	const int nMarker[nmult] = {24, 25, 26, 27}; 
	const int nColor[nmult] = {1, 2, 4, 6}; 

	TH1D *hXib_pt[nmult];
	TH1D *hXib_pt_ratio[nmult];

	float nevent[nmult] = {0};

	for (int ii=0; ii<nmult; ii++){

		nevent[ii] = href[0]->Integral(mult_dn[ii]+1, mult_up[ii]);

		cout << ii << " " << nevent[ii] << endl;

		hXib_pt[ii] = hXib[0]->ProjectionX(Form("hXib_pt_mult%d",ii),mult_dn[ii]+1,mult_up[ii]); 
		hXib_pt[ii]->Rebin(4);
		hXib_pt[ii]->Sumw2();

		/*
		float norm = hXib_pt[ii]->GetBinWidth(1)*nevent[ii];
		hXib_pt[ii]->Scale(1./norm);
		*/

		hXib_pt[ii]->SetMarkerStyle(nMarker[ii]);
		hXib_pt[ii]->SetMarkerColor(nColor[ii]);
		hXib_pt[ii]->SetLineColor(nColor[ii]);

		hXib_pt_ratio[ii] = (TH1D*)hXib_pt[ii]->Clone(Form("hXib_pt_ratio_mult%d",ii));
		hXib_pt_ratio[ii]->Divide(hXib_pt[0]);

	}//ii

	//return;

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*500,500);
	c1->Divide(2,1);

	{
		c1->cd(1);
		gPad->SetMargin(0.12,0.03,0.12,0.03);
		gPad->SetLogy();
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0, 1, 40, 2*hXib_pt[0]->GetMaximum());
		htmp->GetXaxis()->SetTitle("#Xi_{b} p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.045);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.045);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hXib_pt[0]->Draw("p same");
		hXib_pt[1]->Draw("p same");
		hXib_pt[2]->Draw("p same");
		hXib_pt[3]->Draw("p same");

		TLegend *leg = new TLegend(0.15,0.15,0.6,0.15+0.05*6);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","Monash","h");
		leg->AddEntry(hXib_pt[0],"0-100%","P");
		leg->AddEntry(hXib_pt[1],"30-100%","P");
		leg->AddEntry(hXib_pt[2],"0.1-30%","P");
		leg->AddEntry(hXib_pt[3],"0-0.1%","P");
		leg->Draw();
	}

	{

		c1->cd(2);
		gPad->SetMargin(0.12,0.03,0.12,0.03);
		gPad->SetLogy();
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0, 0.01, 40, 2);
		htmp->GetXaxis()->SetTitle("#Xi_{b} p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.045);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.045);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hXib_pt_ratio[0]->Draw("p same");
		hXib_pt_ratio[1]->Draw("p same");
		hXib_pt_ratio[2]->Draw("p same");
		hXib_pt_ratio[3]->Draw("p same");
	}


	/*
	for (int ii=0; ii<href->GetNbinsX(); ii++){
		cout << "0 - " << ii+1 << " " << href->Integral(1,ii+1)/href->Integral() << endl;
	}
	*/

}

void ScanV0M_01(){

	const int ngrp = 1;

	TFile *infile[ngrp];

	TH1D *href[ngrp];
	TH2D *hXic[ngrp][2];
	TH2D *hLc[ngrp][2];

	for (int ig=0; ig<ngrp; ig++){

		infile[ig] = new TFile(Form("outfile_hist_set62_grp%03d.root",ig),"read");

		href[ig] = (TH1D*)infile[ig]->Get("hmult_v0m");

		for (int ii=0; ii<2; ii++){
			hXic[ig][ii] = (TH2D*)infile[ig]->Get(Form("hXic_pt_v0m_%d",ii));
			hLc[ig][ii] = (TH2D*)infile[ig]->Get(Form("hLc_pt_v0m_%d",ii));
		}


		if ( ig>0 ){
			href[0]->Add(href[ig]);

			for (int ii=0; ii<2; ii++){
				hXic[0][ii]->Add(hXic[ig][ii]);
				hLc[0][ii]->Add(hLc[ig][ii]);
			}
		}//igrp
	}

	const int nmult = 4;
	const int mult_dn[nmult] = {0, 0, 92, 242}; 
	const int mult_up[nmult] = {500, 92, 242, 500}; 

	const int nMarker[nmult] = {24, 25, 26, 27}; 
	const int nColor[nmult] = {1, 2, 4, 6}; 

	TH1D *hXic_pt[nmult][2];
	TH1D *hLc_pt[nmult][2];

	TH1D *hXic_pt_ratio[nmult][2];
	TH1D *hLc_pt_ratio[nmult][2];

	float nevent[nmult] = {0};

	for (int ii=0; ii<nmult; ii++){

		nevent[ii] = href[0]->Integral(mult_dn[ii]+1, mult_up[ii]);

		cout << ii << " " << nevent[ii] << endl;

		for (int jj=0; jj<2; jj++){
			hXic_pt[ii][jj] = hXic[0][jj]->ProjectionX(Form("hXic_pt_mult%d_%d",ii,jj),mult_dn[ii]+1,mult_up[ii]); 
			hXic_pt[ii][jj]->Rebin(4);
			hXic_pt[ii][jj]->Sumw2();

			hLc_pt[ii][jj] = hLc[0][jj]->ProjectionX(Form("hLc_pt_mult%d_%d",ii,jj),mult_dn[ii]+1,mult_up[ii]); 
			hLc_pt[ii][jj]->Rebin(4);
			hLc_pt[ii][jj]->Sumw2();

			hXic_pt[ii][jj]->SetMarkerStyle(nMarker[ii]);
			hXic_pt[ii][jj]->SetMarkerColor(nColor[ii]);
			hXic_pt[ii][jj]->SetLineColor(nColor[ii]);

			hLc_pt[ii][jj]->SetMarkerStyle(nMarker[ii]);
			hLc_pt[ii][jj]->SetMarkerColor(nColor[ii]);
			hLc_pt[ii][jj]->SetLineColor(nColor[ii]);

			hXic_pt_ratio[ii][jj] = (TH1D*)hXic_pt[ii][jj]->Clone(Form("hXic_pt_ratio_mult%d_%d",ii,jj));
			hLc_pt_ratio[ii][jj] = (TH1D*)hLc_pt[ii][jj]->Clone(Form("hLc_pt_ratio_mult%d_%d",ii,jj));
		}

		hXic_pt_ratio[ii][0]->Add(hXic_pt_ratio[ii][1]);
		hLc_pt_ratio[ii][0]->Add(hLc_pt_ratio[ii][1]);

		hXic_pt_ratio[ii][1]->Divide(hXic_pt_ratio[ii][0]);
		hLc_pt_ratio[ii][1]->Divide(hLc_pt_ratio[ii][0]);

	}//ii

	//return;

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*500,500);
	c1->Divide(2,1);

	{
		c1->cd(1);
		gPad->SetMargin(0.12,0.03,0.12,0.03);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0, 0.00, 40, 0.5);
		htmp->GetXaxis()->SetTitle("#Lambda_{c}^{+} p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.045);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.045);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hLc_pt_ratio[0][1]->Draw("p same");
		hLc_pt_ratio[1][1]->Draw("p same");
		hLc_pt_ratio[2][1]->Draw("p same");
		hLc_pt_ratio[3][1]->Draw("p same");
	}

	{
		c1->cd(2);
		gPad->SetMargin(0.12,0.03,0.12,0.03);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0, 0.00, 40, 0.5);
		htmp->GetXaxis()->SetTitle("#Xi_{c}^{0} p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.045);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.045);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hXic_pt_ratio[0][1]->Draw("p same");
		hXic_pt_ratio[1][1]->Draw("p same");
		hXic_pt_ratio[2][1]->Draw("p same");
		hXic_pt_ratio[3][1]->Draw("p same");
	}

	/*
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
	*/


	/*
	for (int ii=0; ii<href->GetNbinsX(); ii++){
		cout << "0 - " << ii+1 << " " << href->Integral(1,ii+1)/href->Integral() << endl;
	}
	*/

}

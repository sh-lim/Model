void DrawAMPTratio(){

	const int nset = 3;
	const int nColor[3] = {1, 2, 4};

	string setname[nset] = {"p+Au", "d+Au", "^{3}He+Au"};

	TFile *infile = new TFile("outfile_AMPT_ratio.root","read");

	TGraphErrors *gon_ratio[nset];
	TGraphErrors *goff_ratio[nset];

	for (int iset=0; iset<nset; iset++){

		gon_ratio[iset] = (TGraphErrors*)infile->Get(Form("gon_ratio_set%d",iset));
		goff_ratio[iset] = (TGraphErrors*)infile->Get(Form("goff_ratio_set%d",iset));

		gon_ratio[iset]->SetFillColorAlpha(nColor[iset],0.5);
		goff_ratio[iset]->SetFillColorAlpha(nColor[iset],0.5);

	}

	TCanvas *c3 = new TCanvas("c3","c3",1.1*400,400);
	{
		gPad->SetMargin(0.14,0.03,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.0,3.0,1.5);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("Ratio (0-20% / 60-100%)");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.2,0.15,0.5,0.4);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","AMPT 200 GeV","h");
		leg->AddEntry("","parton scattering on","h");


		for (int iset=0; iset<nset; iset++){
			gon_ratio[iset]->Draw("e3");
			leg->AddEntry(gon_ratio[iset],setname[iset].c_str(),"F");
		}//iset
		leg->Draw();
	}//jj


	TCanvas *c4 = new TCanvas("c4","c4",1.1*400,400);
	{
		gPad->SetMargin(0.14,0.03,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.0,3.0,1.5);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("Ratio (0-20% / 60-100%)");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.2,0.15,0.5,0.4);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","AMPT 200 GeV","h");
		leg->AddEntry("","parton scattering off","h");


		for (int iset=0; iset<nset; iset++){
			goff_ratio[iset]->Draw("e3");
			leg->AddEntry(goff_ratio[iset],setname[iset].c_str(),"F");
		}//iset
		leg->Draw();
	}//jj

}

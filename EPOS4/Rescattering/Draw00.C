void Draw00(){

	const int nset = 4;

	int nMarker[4] = {20, 21, 24, 25};
	int nColor[4] = {1, 2, 4, 6};

	TFile *infile[nset];
	infile[0] = new TFile("Ana00_pp13TeV_opt1_grp110.root","read");
	infile[1] = new TFile("Ana00_pp13TeV_opt2_grp110.root","read");
	infile[2] = new TFile("Ana00_pp13TeV_opt3_grp110.root","read");
	infile[3] = new TFile("Ana00_pp13TeV_opt4_grp110.root","read");

	TH1D *hmult[nset];
	TH1D *heta[nset];
	TH1D *hmult_211[nset];
	TH1D *hmult_321[nset];
	TH1D *hmult_113[nset];
	TH1D *hmult_313[nset];
	TH1D *hmult_333[nset];

	TH1D *hratio_333_321[nset];
	TH1D *hratio_313_321[nset];
	TH1D *hratio_113_211[nset];
	TH1D *hratio_321_211[nset];

	for (int iset=0; iset<nset; iset++){

		hmult[iset] = (TH1D*)infile[iset]->Get("hmult");
		hmult[iset]->SetLineColor(nColor[iset]);

		cout << iset << " " << hmult[iset]->Integral() << " " << hmult[iset]->Integral(2,-1) << endl;

		heta[iset] = (TH1D*)infile[iset]->Get("heta");
		heta[iset]->Scale(1./(hmult[iset]->Integral(2,-1)*heta[iset]->GetBinWidth(1)));
		heta[iset]->SetLineColor(nColor[iset]);
		heta[iset]->SetMarkerColor(nColor[iset]);
		heta[iset]->SetMarkerStyle(nMarker[iset]);
		heta[iset]->SetMarkerSize(0.6);

		hmult_211[iset] = (TH1D*)infile[iset]->Get("hmult_211");
		hmult_321[iset] = (TH1D*)infile[iset]->Get("hmult_321");
		hmult_113[iset] = (TH1D*)infile[iset]->Get("hmult_113");
		hmult_313[iset] = (TH1D*)infile[iset]->Get("hmult_313");
		hmult_333[iset] = (TH1D*)infile[iset]->Get("hmult_333");

		hmult_211[iset]->Rebin(10);
		hmult_321[iset]->Rebin(10);
		hmult_113[iset]->Rebin(10);
		hmult_313[iset]->Rebin(10);
		hmult_333[iset]->Rebin(10);

		hmult_211[iset]->Sumw2();
		hmult_321[iset]->Sumw2();
		hmult_113[iset]->Sumw2();
		hmult_313[iset]->Sumw2();
		hmult_333[iset]->Sumw2();

		if ( iset<2 ){
			hmult_333[iset]->Scale(1./0.49);
			hmult_313[iset]->Scale(1./0.67);
		}
		hmult_211[iset]->Scale(0.5);
		hmult_321[iset]->Scale(0.5);
		hmult_313[iset]->Scale(0.5);

		hratio_313_321[iset] = (TH1D*)hmult_313[iset]->Clone(Form("hratio_313_321_%d",iset));
		hratio_313_321[iset]->Divide(hmult_321[iset]);
		hratio_313_321[iset]->SetMarkerStyle(nMarker[iset]);
		hratio_313_321[iset]->SetMarkerColor(nColor[iset]);
		hratio_313_321[iset]->SetLineColor(nColor[iset]);

		hratio_333_321[iset] = (TH1D*)hmult_333[iset]->Clone(Form("hratio_333_321_%d",iset));
		hratio_333_321[iset]->Divide(hmult_321[iset]);
		hratio_333_321[iset]->SetMarkerStyle(nMarker[iset]);
		hratio_333_321[iset]->SetMarkerColor(nColor[iset]);
		hratio_333_321[iset]->SetLineColor(nColor[iset]);

		hratio_321_211[iset] = (TH1D*)hmult_321[iset]->Clone(Form("hratio_321_211_%d",iset));
		hratio_321_211[iset]->Divide(hmult_211[iset]);
		hratio_321_211[iset]->SetMarkerStyle(nMarker[iset]);
		hratio_321_211[iset]->SetMarkerColor(nColor[iset]);
		hratio_321_211[iset]->SetLineColor(nColor[iset]);

		hratio_113_211[iset] = (TH1D*)hmult_113[iset]->Clone(Form("hratio_113_211_%d",iset));
		hratio_113_211[iset]->Divide(hmult_211[iset]);
		hratio_113_211[iset]->SetMarkerStyle(nMarker[iset]);
		hratio_113_211[iset]->SetMarkerColor(nColor[iset]);
		hratio_113_211[iset]->SetLineColor(nColor[iset]);

	}//iset

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	c1->SetMargin(0.14,0.03,0.12,0.03);

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1,200,2*hmult[0]->GetMaximum());
		gPad->SetLogy();

		hmult[0]->Draw("same");
		hmult[1]->Draw("same");
		//hmult[2]->Draw("same");
		//hmult[3]->Draw("same");

	}

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	c2->SetMargin(0.14,0.03,0.12,0.03);

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(-5,0,5,10);
		htmp->GetXaxis()->SetTitle("#eta");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("dN_{ch}/d#eta");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.75,0.25,0.9,0.45);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		for (int ii=0; ii<nset; ii++){
			heta[ii]->Draw("p same");
			leg->AddEntry(heta[ii],Form("opt%d",ii+1),"p");
		}
		leg->Draw();
	}

	return;

	TCanvas *c6 = new TCanvas("c6","c6",1.2*500,500);
	c6->SetMargin(0.14,0.03,0.12,0.03);

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.3);
		htmp->GetXaxis()->SetTitle("N_{ch, |eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("#rho/#pi");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.75,0.25,0.9,0.45);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		for (int ii=0; ii<nset; ii++){
			hratio_113_211[ii]->Draw("p same");
			leg->AddEntry(hratio_313_321[ii],Form("opt%d",ii+1),"p");
		}
		leg->Draw();
	}

	//return;

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	c3->SetMargin(0.14,0.03,0.12,0.03);

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.2);
		htmp->GetXaxis()->SetTitle("N_{ch, |eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("K/#pi");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.75,0.25,0.9,0.45);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		for (int ii=0; ii<nset; ii++){
			hratio_321_211[ii]->Draw("p same");
			leg->AddEntry(hratio_321_211[ii],Form("opt%d",ii+1),"p");
		}
		leg->Draw();
	}

	TCanvas *c4 = new TCanvas("c4","c4",1.2*500,500);
	c4->SetMargin(0.14,0.03,0.12,0.03);

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.25);
		htmp->GetXaxis()->SetTitle("N_{ch, |eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("#phi/K");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.75,0.25,0.9,0.45);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		for (int ii=0; ii<nset; ii++){
			hratio_333_321[ii]->Draw("p same");
			leg->AddEntry(hratio_333_321[ii],Form("opt%d",ii+1),"p");
		}
		leg->Draw();
	}

	TCanvas *c5 = new TCanvas("c5","c5",1.2*500,500);
	c5->SetMargin(0.14,0.03,0.12,0.03);

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.6);
		htmp->GetXaxis()->SetTitle("N_{ch, |eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("K^{*0}/K");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.75,0.25,0.9,0.45);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		for (int ii=0; ii<nset; ii++){
			hratio_313_321[ii]->Draw("p same");
			leg->AddEntry(hratio_313_321[ii],Form("opt%d",ii+1),"p");
		}
		leg->Draw();
	}


	return;

	/*
	*/

	/*
	*/



}

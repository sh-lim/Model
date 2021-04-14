void DrawPt(){

	const int nset = 3;
	const int neta = 3;
	const int norder = 2;

	const int nColor[3] = {1, 2, 4};

	TFile *infile[nset];
	TFile *infile_off[nset];
	infile[0] = new TFile("outfile_pt_pAu200GeV_grp2.root","read");
	infile[1] = new TFile("outfile_pt_dAu200GeV_grp2.root","read");
	infile[2] = new TFile("outfile_pt_He3Au200GeV_grp2.root","read");

	infile_off[0] = new TFile("outfile_pt_pAu200GeV_grp3.root","read");
	infile_off[1] = new TFile("outfile_pt_dAu200GeV_grp3.root","read");
	infile_off[2] = new TFile("outfile_pt_He3Au200GeV_grp3.root","read");

	string setname[nset] = {"p+Au", "d+Au", "^{3}He+Au"};

	TProfile *hprof_meanpt[nset];
	TProfile *hprof_meanpt_sT[nset];

	TProfile *hprof_meanpt_off[nset];
	TProfile *hprof_meanpt_sT_off[nset];

	TGraphErrors *gmeanpt[nset];

	for (int iset=0; iset<nset; iset++){
		hprof_meanpt[iset] = (TProfile*)infile[iset]->Get(Form("hprof_meanpt"));
		hprof_meanpt_sT[iset] = (TProfile*)infile[iset]->Get(Form("hprof_meanpt_sT"));

		hprof_meanpt_off[iset] = (TProfile*)infile_off[iset]->Get(Form("hprof_meanpt"));
		hprof_meanpt_sT_off[iset] = (TProfile*)infile_off[iset]->Get(Form("hprof_meanpt_sT"));
	}


	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,400);
	c1->Divide(2,1);

	{
		c1->cd(1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.3,100,0.8);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta|_{#eta=0}");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("#LTp_{T}#GT");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.6,0.15,0.9,0.4);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","AMPT 200 GeV","h");
		leg->AddEntry("","parton scattering on","h");

		for (int iset=0; iset<nset; iset++){
			hprof_meanpt[iset]->SetMarkerStyle(24);
			hprof_meanpt[iset]->SetLineColor(nColor[iset]);
			hprof_meanpt[iset]->SetMarkerColor(nColor[iset]);
			hprof_meanpt[iset]->SetLineWidth(2);
			hprof_meanpt[iset]->Draw("p same");
			leg->AddEntry(hprof_meanpt[iset],setname[iset].c_str(),"P");
		}//iset
		leg->Draw();
	}//jj


	{
		c1->cd(2);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.3,100,0.8);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta|_{#eta=0}/#LTS_{T}#GT");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("#LTp_{T}#GT");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.6,0.15,0.9,0.4);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","AMPT 200 GeV","h");
		leg->AddEntry("","parton scattering on","h");

		for (int iset=0; iset<nset; iset++){
			hprof_meanpt_sT[iset]->SetMarkerStyle(24);
			hprof_meanpt_sT[iset]->SetLineColor(nColor[iset]);
			hprof_meanpt_sT[iset]->SetMarkerColor(nColor[iset]);
			hprof_meanpt_sT[iset]->SetLineWidth(2);
			hprof_meanpt_sT[iset]->Draw("p same");
			leg->AddEntry(hprof_meanpt_sT[iset],setname[iset].c_str(),"P");
		}//iset
		leg->Draw();
	}//jj

	TCanvas *c2 = new TCanvas("c2","c2",1.1*2*400,400);
	c2->Divide(2,1);

	{
		c2->cd(1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.3,100,0.8);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta|_{#eta=0}");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("#LTp_{T}#GT");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.6,0.15,0.9,0.4);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","AMPT 200 GeV","h");
		leg->AddEntry("","parton scattering off","h");

		for (int iset=0; iset<nset; iset++){
			hprof_meanpt_off[iset]->SetMarkerStyle(24);
			hprof_meanpt_off[iset]->SetLineColor(nColor[iset]);
			hprof_meanpt_off[iset]->SetMarkerColor(nColor[iset]);
			hprof_meanpt_off[iset]->SetLineWidth(2);
			hprof_meanpt_off[iset]->Draw("p same");
			leg->AddEntry(hprof_meanpt_off[iset],setname[iset].c_str(),"P");
		}//iset
		leg->Draw();
	}//jj

	{
		c2->cd(2);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.3,100,0.8);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta|_{#eta=0}/#LTS_{T}#GT");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("#LTp_{T}#GT");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.6,0.15,0.9,0.4);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","AMPT 200 GeV","h");
		leg->AddEntry("","parton scattering off","h");

		for (int iset=0; iset<nset; iset++){
			hprof_meanpt_sT_off[iset]->SetMarkerStyle(24);
			hprof_meanpt_sT_off[iset]->SetLineColor(nColor[iset]);
			hprof_meanpt_sT_off[iset]->SetMarkerColor(nColor[iset]);
			hprof_meanpt_sT_off[iset]->SetLineWidth(2);
			hprof_meanpt_sT_off[iset]->Draw("p same");
			leg->AddEntry(hprof_meanpt_sT_off[iset],setname[iset].c_str(),"P");
		}//iset
		leg->Draw();
	}//jj

}

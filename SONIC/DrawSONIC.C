void DrawSONIC(){

	const int nset = 3;
	const int nColor[nset] = {1, 2, 4};

	const char *setname[nset] = {"pO","OO","pPb"};

	TFile *infile_pre[nset];
	TFile *infile_vn[nset];

	TH1D *hdNchdeta[nset];
	TH1D *hvn[nset][3];

	for (int iset=0; iset<nset; iset++){

		infile_pre[iset] = new TFile(Form("pre_outfile_SONIC_%s8160GeV.root",setname[iset]),"read");
		infile_vn[iset] = new TFile(Form("outfile_vn_SONIC_%s8160GeV.root",setname[iset]),"read");

		hdNchdeta[iset] = (TH1D*)infile_pre[iset]->Get("hdNchdeta1");
		hdNchdeta[iset]->SetLineColor(nColor[iset]);

		for (int ii=0; ii<3; ii++){
			hvn[iset][ii] = (TH1D*)infile_vn[iset]->Get(Form("hv%d",ii+2));
			hvn[iset][ii]->SetLineColor(nColor[iset]);
			hvn[iset][ii]->SetMarkerColor(nColor[iset]);
			hvn[iset][ii]->SetMarkerStyle(24);
		}



	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,2*400);
	c1->Divide(2,2);

	{
		c1->cd(1);
		gPad->SetMargin(0.14,0.03,0.12,0.03);
		gPad->SetLogy();
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1,150,2*hdNchdeta[0]->GetMaximum());
		htmp->GetYaxis()->SetTitle("N_{event}");
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta |_{#eta=0}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		hdNchdeta[0]->Draw("same");
		hdNchdeta[1]->Draw("same");
		hdNchdeta[2]->Draw("same");

		TLegend *leg = new TLegend(0.55,0.75,0.9,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.05);
		leg->AddEntry("","SONIC 8.16 TeV","");
		leg->AddEntry(hdNchdeta[0],"p+O","L");
		leg->AddEntry(hdNchdeta[1],"O+O","L");
		leg->AddEntry(hdNchdeta[2],"p+Pb","L");
		leg->Draw();

	}

	for (int ii=0; ii<3; ii++){

		c1->cd(ii+2);

		gPad->SetMargin(0.14,0.03,0.12,0.03);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.15);
		htmp->GetYaxis()->SetTitle(Form("v_{%d}",ii+2));
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		hvn[0][ii]->Draw("p same");
		hvn[1][ii]->Draw("p same");
		hvn[2][ii]->Draw("p same");

	}



}

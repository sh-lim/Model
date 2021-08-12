void DrawSONICV5(){

	const int nset = 4;
	const int nColor[nset] = {1, 2, 1, 2};
	const int nMarker[nset] = {20, 20, 24, 24};

	const char *setname[nset] = {"pAu","dAu","He3Au",""};

	TFile *infile_pre[nset];
	TFile *infile_vn[nset];

	TH1D *hdNchdeta[nset];
	TH1D *hvn[nset][3];

	infile_vn[0] = new TFile(Form("outfile_vn_SONIC_pAu200GeV_IPGlasma_cent05.root"),"read");
	infile_vn[1] = new TFile(Form("outfile_vn_SONIC_pAu200GeV_IPGlasma2_cent05.root"),"read");
	infile_vn[2] = new TFile(Form("outfile_vn_SONIC_pAu200GeV_IPGlasma_rmax3_cent05.root"),"read");
	infile_vn[3] = new TFile(Form("outfile_vn_SONIC_pAu200GeV_IPGlasma2_rmax3_cent05.root"),"read");

	TGraphErrors *gvn[nset][3];

	for (int iset=0; iset<nset; iset++){
		for (int ii=0; ii<3; ii++){
			hvn[iset][ii] = (TH1D*)infile_vn[iset]->Get(Form("hv%d_0",ii+2));
			hvn[iset][ii]->SetLineColor(nColor[iset]);
			hvn[iset][ii]->SetMarkerColor(nColor[iset]);
			hvn[iset][ii]->SetMarkerStyle(24);

			gvn[iset][ii] = new TGraphErrors();
			gvn[iset][ii]->SetLineColorAlpha(nColor[iset], 0.5);
			gvn[iset][ii]->SetFillColorAlpha(nColor[iset], 0.5);

			for (int ipt=0; ipt<hvn[iset][ii]->GetNbinsX(); ipt++){
				float xx = hvn[iset][ii]->GetBinCenter(ipt+1);
				float yy = hvn[iset][ii]->GetBinContent(ipt+1);
				float yy_err = hvn[iset][ii]->GetBinError(ipt+1);

				gvn[iset][ii]->SetPoint(ipt, xx, yy);
				gvn[iset][ii]->SetPointError(ipt, xx, yy_err);
			}//
		}

	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,1*400);
	c1->Divide(2,1);

	for (int ii=0; ii<2; ii++){

		c1->cd(ii+1);

		gPad->SetMargin(0.14,0.03,0.12,0.03);
		float ymax = 0.20;
		if ( ii==1 || ii==2 ){
			ymax = 0.10;
		}
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,ymax);
		htmp->GetYaxis()->SetTitle(Form("v_{%d}",ii+2));
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		gvn[0][ii]->Draw("3");
		gvn[1][ii]->Draw("3");
		gvn[2][ii]->Draw("3");
		gvn[3][ii]->Draw("3");

	}



}

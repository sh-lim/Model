void CalPsiRes(){

	const int nset = 3;
	const int neta = 3;
	const int norder = 2;

	const int nColor[3] = {1, 2, 4};

	TFile *infile[nset];
	infile[0] = new TFile("outfile_pAu200GeV_grp0.root","read");
	infile[1] = new TFile("outfile_dAu200GeV_grp2.root","read");
	infile[2] = new TFile("outfile_He3Au200GeV_grp2.root","read");

	string setname[nset] = {"p+Au", "d+Au", "^{3}He+Au"};

	TH1D *hcosnPsi[nset][neta][norder];
	TH1D *hcosndphi[nset][neta][norder];

	TGraphErrors *gvn[nset][neta][norder];

	for (int iset=0; iset<nset; iset++){
		for (int ieta=0; ieta<neta; ieta++){
			for (int io=0; io<norder; io++){
				hcosnPsi[iset][ieta][io] = (TH1D*)infile[iset]->Get(Form("hcosnPsi_eta%d_order%d",ieta,io+2));
				hcosndphi[iset][ieta][io] = (TH1D*)infile[iset]->Get(Form("hcosndphi_eta%d_order%d",ieta,io+2));
			}
		}
	}

	TH1D *hppEcc[nset][norder];
	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){
			hppEcc[iset][io] = (TH1D*)infile[iset]->Get(Form("hppEcc_order%d",io+2));
			hppEcc[iset][io]->SetLineColor(nColor[iset]);
			hppEcc[iset][io]->Scale(1./hppEcc[iset][io]->Integral());
		}
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,400);
	c1->Divide(2,1);

	for (int jj=0; jj<2; jj++){
		c1->cd(jj+1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,1,0.03);
		htmp->GetXaxis()->SetTitle(Form("#varepsilon_{%d}",jj+2));
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle();
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
		leg->SetTextSize(0.045);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);

		for (int iset=0; iset<nset; iset++){
			hppEcc[iset][jj]->Draw("same");
			leg->AddEntry(hppEcc[iset][jj],setname[iset].c_str(),"L");
		}//iset
		leg->Draw();
	}//jj


	//return;
	float EPres[nset][neta][norder];


	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){

			//CNT-FVTXS
			float mean0 = hcosnPsi[iset][0][io]->GetMean();
			//FVTXS-BBCS
			float mean1 = hcosnPsi[iset][1][io]->GetMean();
			//BBCS-CNT
			float mean2 = hcosnPsi[iset][2][io]->GetMean();

			//CNT
			EPres[iset][0][io] = sqrt(mean0*mean2/mean1);

			//FVTX-S
			EPres[iset][1][io] = sqrt(mean0*mean1/mean2);

			//BBC-S
			EPres[iset][2][io] = sqrt(mean1*mean2/mean0);

			cout << "sys: " << iset << ", order: " << io+2 << endl;
			cout << "CNT: " << EPres[iset][0][io] << ", FVTXS: " << EPres[iset][1][io] << ", BBCS: " << EPres[iset][2][io] << endl;

		}
	}

	for (int iset=0; iset<nset; iset++){
		for (int ieta=0; ieta<neta; ieta++){
			for (int io=0; io<norder; io++){
				hcosndphi[iset][ieta][io]->Scale(1./EPres[iset][ieta][io]);
				hcosndphi[iset][ieta][io]->SetMarkerStyle(23+ieta);
				hcosndphi[iset][ieta][io]->SetMarkerColor(1+io);
				hcosndphi[iset][ieta][io]->SetLineColor(1+io);

				gvn[iset][ieta][io] = new TGraphErrors;
				gvn[iset][ieta][io]->SetLineColor(1+io);
				gvn[iset][ieta][io]->SetLineWidth(2);
				gvn[iset][ieta][io]->SetFillColorAlpha(1+io,0.2);

				for (int ib=0; ib<hcosndphi[iset][ieta][io]->GetNbinsX(); ib++){
					float xx = hcosndphi[iset][ieta][io]->GetBinCenter(ib+1);
					float xx_err = 0.5*hcosndphi[iset][ieta][io]->GetBinWidth(ib+1);
					float yy = hcosndphi[iset][ieta][io]->GetBinContent(ib+1);
					float yy_err = hcosndphi[iset][ieta][io]->GetBinError(ib+1);

					gvn[iset][ieta][io]->SetPoint(ib, xx, yy);
					gvn[iset][ieta][io]->SetPointError(ib, xx_err, yy_err);
				}

			}//io
		}//ieta
	}//iset

	TCanvas *c2 = new TCanvas("c2","c2",1.1*3*400,400);
	c2->Divide(3,1);

	for (int iset=0; iset<nset; iset++){
		c2->cd(iset+1);

		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.25);
		htmp->GetXaxis()->SetTitle();
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle();
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		hcosndphi[iset][1][0]->Draw("p same");
		hcosndphi[iset][1][1]->Draw("p same");

		TLegend *leg = new TLegend(0.6,0.75,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",setname[iset].c_str(),"h");
		leg->AddEntry(hcosndphi[iset][1][0],"v_{2}","P");
		leg->AddEntry(hcosndphi[iset][1][1],"v_{3}","P");
		leg->Draw();

	}

	TCanvas *c3 = new TCanvas("c3","c3",1.1*3*400,400);
	c3->Divide(3,1);

	for (int iset=0; iset<nset; iset++){
		c3->cd(iset+1);

		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.25);
		htmp->GetXaxis()->SetTitle();
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle();
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		gvn[iset][1][0]->Draw("e3");
		gvn[iset][1][1]->Draw("e3");

		TLegend *leg = new TLegend(0.6,0.75,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",setname[iset].c_str(),"h");
		leg->AddEntry(gvn[iset][1][0],"v_{2}","L");
		leg->AddEntry(gvn[iset][1][1],"v_{3}","L");
		leg->Draw();

	}

}

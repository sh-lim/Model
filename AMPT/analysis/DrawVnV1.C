void DrawVnV1(){

	const int nset = 2;
	const int neta = 3;
	const int norder = 1;
	const int npid = 4;

	TFile *infile[nset];

	infile[0] = new TFile("outfileVnV1_pPb5TeV_grp0.root","read");
	infile[1] = new TFile("outfileVnV1_pPb5TeV_grp10.root","read");

	const int nColor[3] = {1, 2, 4};

	TH1D *hcosndPsi[nset][neta][norder];

	TProfile *hcosndphi[nset][neta][norder][npid];
	TProfile *hcosndphiP[nset][norder][npid];

	TGraphErrors *gvn[nset][neta][norder][npid];
	TGraphErrors *gvnP[nset][norder][npid];

	for (int iset=0; iset<nset; iset++){
		for (int ieta=0; ieta<neta; ieta++){
			for (int io=0; io<norder; io++){
				hcosndPsi[iset][ieta][io] = (TH1D*)infile[iset]->Get(Form("hcosndPsi_eta%d_order%d",ieta,io+2));

				for (int ipid=0; ipid<npid; ipid++){
					hcosndphi[iset][ieta][io][ipid] = (TProfile*)infile[iset]->Get(Form("hcosndphi_eta%d_order%d_pid%d",ieta,io+2,ipid));
				}
			}
		}
	}

	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){
			for (int ipid=0; ipid<npid; ipid++){
				hcosndphiP[iset][io][ipid] = (TProfile*)infile[iset]->Get(Form("hcosndphiP_order%d_pid%d",io+2,ipid));
			}
		}
	}


	//return;
	float EPres[nset][neta][norder];

	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){

			//CNT-FVTXS
			float mean0 = hcosndPsi[iset][0][io]->GetMean();
			//FVTXS-BBCS
			float mean1 = hcosndPsi[iset][1][io]->GetMean();
			//BBCS-CNT
			float mean2 = hcosndPsi[iset][2][io]->GetMean();

			//CNT
			EPres[iset][0][io] = sqrt(mean0*mean2/mean1);

			//FVTX-S
			EPres[iset][1][io] = sqrt(mean0*mean1/mean2);

			//BBC-S
			EPres[iset][2][io] = sqrt(mean1*mean2/mean0);

			cout << "set: " << iset << ", order: " << io+2 << endl;
			cout << "CNT: " << EPres[iset][0][io] << ", FVTXS: " << EPres[iset][1][io] << ", BBCS: " << EPres[iset][2][io] << endl;

		}
	}


	for (int iset=0; iset<nset; iset++){
		for (int ieta=0; ieta<neta; ieta++){
			for (int io=0; io<norder; io++){
				for (int ipid=0; ipid<npid; ipid++){
					hcosndphi[iset][ieta][io][ipid]->Scale(1./EPres[iset][ieta][io]);
					hcosndphi[iset][ieta][io][ipid]->SetMarkerStyle(23+ieta);
					hcosndphi[iset][ieta][io][ipid]->SetMarkerColor(1+io);
					hcosndphi[iset][ieta][io][ipid]->SetLineColor(1+io);

					gvn[iset][ieta][io][ipid] = new TGraphErrors;
					gvn[iset][ieta][io][ipid]->SetLineColor(nColor[ieta]);
					gvn[iset][ieta][io][ipid]->SetLineWidth(2);
					gvn[iset][ieta][io][ipid]->SetFillColorAlpha(nColor[ieta],0.2);

					for (int ib=0; ib<hcosndphi[iset][ieta][io][ipid]->GetNbinsX(); ib++){
						float xx = hcosndphi[iset][ieta][io][ipid]->GetBinCenter(ib+1);
						float xx_err = 0.5*hcosndphi[iset][ieta][io][ipid]->GetBinWidth(ib+1);
						float yy = hcosndphi[iset][ieta][io][ipid]->GetBinContent(ib+1);
						float yy_err = hcosndphi[iset][ieta][io][ipid]->GetBinError(ib+1);

						if ( iset==1 && ieta==0 && io==0 ){
							cout << xx << " " << yy << endl;
						}

						gvn[iset][ieta][io][ipid]->SetPoint(ib, xx, yy);
						gvn[iset][ieta][io][ipid]->SetPointError(ib, xx_err, 0.005+yy_err);
					}

				}//ipid
			}//io
		}//ieta
	}

	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){
			for (int ipid=0; ipid<npid; ipid++){
				hcosndphiP[iset][io][ipid]->SetMarkerStyle(kGreen+2);
				hcosndphiP[iset][io][ipid]->SetMarkerColor(kGreen+2);
				hcosndphiP[iset][io][ipid]->SetLineColor(kGreen+2);

				gvnP[iset][io][ipid] = new TGraphErrors;
				gvnP[iset][io][ipid]->SetLineColor(kGreen+2);
				gvnP[iset][io][ipid]->SetLineWidth(2);
				gvnP[iset][io][ipid]->SetFillColorAlpha(kGreen+2,0.2);

				for (int ib=0; ib<hcosndphiP[iset][io][ipid]->GetNbinsX(); ib++){
					float xx = hcosndphiP[iset][io][ipid]->GetBinCenter(ib+1);
					float xx_err = 0.5*hcosndphiP[iset][io][ipid]->GetBinWidth(ib+1);
					float yy = hcosndphiP[iset][io][ipid]->GetBinContent(ib+1);
					float yy_err = hcosndphiP[iset][io][ipid]->GetBinError(ib+1);

					gvnP[iset][io][ipid]->SetPoint(ib, xx, yy);
					gvnP[iset][io][ipid]->SetPointError(ib, xx_err, 0.005+yy_err);
				}

			}//ipid
		}//io
	}


	TCanvas *c3 = new TCanvas("c3","c3",1.1*2*500,500);
	c3->Divide(2,1);

	for (int iset=0; iset<nset; iset++){

		c3->cd(iset+1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,-0.05,10,0.8);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("v_{n}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		gvn[iset][0][0][0]->Draw("e3");
		gvn[iset][1][0][0]->Draw("e3");
		gvn[iset][2][0][0]->Draw("e3");
		gvnP[iset][0][0]->Draw("e3");

		TLegend *leg = new TLegend(0.6,0.5,0.9,0.85);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","AMPT pPb 5.02 TeV","h");
		if ( iset==0 ){
			leg->AddEntry("","Parton scattering on","h");
		}else{
			leg->AddEntry("","Parton scattering off","h");
		}
		leg->AddEntry("","Charged hadron, |#eta|<1","h");
		leg->AddEntry(gvn[iset][0][0][0],"EP, -1<#eta<1","F");
		leg->AddEntry(gvn[iset][1][0][0],"EP, -3<#eta<-1","F");
		leg->AddEntry(gvn[iset][2][0][0],"EP, -5<#eta<-3","F");
		leg->AddEntry(gvnP[iset][0][0],"PP","F");
		leg->Draw();
	}

	/*
	*/

}

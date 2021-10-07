void DrawVnV2(){

	const int nset = 4;
	const int neta = 3;
	const int norder = 1;
	const int npid = 4;

	TFile *infile[nset];

	infile[0] = new TFile("outfileVnV2_AMPT_pPb_5020GeV_grp005_try101.root","read");
	infile[1] = new TFile("outfileVnV2_AMPT_pPb_5020GeV_grp015_try101.root","read");
	infile[2] = new TFile("outfileVnV2_AMPT_pPb_5020GeV_grp025_try101.root","read");
	infile[3] = new TFile("outfileVnV2_AMPT_pPb_5020GeV_grp035_try101.root","read");

	const int nColor[4] = {1, 2, 4, 6};

	TH1D *hcosndPsi[nset][neta][norder];

	TProfile *hcosndphi[nset][neta][norder][npid];
	TProfile *hcosndphiP[nset][norder][npid];
	TProfile *hcosndphiInitP[nset][norder][2];
	TProfile *hcosndphiFinalP[nset][norder][2];
	TProfile *hcosndphiFinalPH[nset][norder][2];

	TGraphErrors *gvn[nset][neta][norder][npid];
	TGraphErrors *gvnP[nset][norder][npid];
	TGraphErrors *gvnInitP[nset][norder][npid];
	TGraphErrors *gvnFinalP[nset][norder][npid];
	TGraphErrors *gvnFinalPH[nset][norder][npid];

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

	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){
			for (int ipid=0; ipid<2; ipid++){
				hcosndphiInitP[iset][io][ipid] = (TProfile*)infile[iset]->Get(Form("hcosndphiInitParton_order%d_pid%d",io+2,ipid));
				hcosndphiFinalP[iset][io][ipid] = (TProfile*)infile[iset]->Get(Form("hcosndphiFinalParton_order%d_pid%d",io+2,ipid));
				hcosndphiFinalPH[iset][io][ipid] = (TProfile*)infile[iset]->Get(Form("hcosndphiFinalPartonH_order%d_pid%d",io+2,ipid));
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
					gvn[iset][ieta][io][ipid]->SetLineWidth(0);
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
						gvn[iset][ieta][io][ipid]->SetPointError(ib, xx_err, 0.001+yy_err);
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
				gvnP[iset][io][ipid]->SetLineWidth(0);
				gvnP[iset][io][ipid]->SetFillColorAlpha(kGreen+2,0.2);

				for (int ib=0; ib<hcosndphiP[iset][io][ipid]->GetNbinsX(); ib++){
					float xx = hcosndphiP[iset][io][ipid]->GetBinCenter(ib+1);
					float xx_err = 0.5*hcosndphiP[iset][io][ipid]->GetBinWidth(ib+1);
					float yy = hcosndphiP[iset][io][ipid]->GetBinContent(ib+1);
					float yy_err = hcosndphiP[iset][io][ipid]->GetBinError(ib+1);

					gvnP[iset][io][ipid]->SetPoint(ib, xx, yy);
					gvnP[iset][io][ipid]->SetPointError(ib, xx_err, 0.001+yy_err);
				}

			}//ipid
		}//io
	}

	for (int iset=0; iset<nset; iset++){
		for (int io=0; io<norder; io++){
			for (int ipid=0; ipid<2; ipid++){
				hcosndphiInitP[iset][io][ipid]->SetMarkerStyle(kGreen+2);
				hcosndphiInitP[iset][io][ipid]->SetMarkerColor(kGreen+2);
				hcosndphiInitP[iset][io][ipid]->SetLineColor(kGreen+2);

				hcosndphiFinalP[iset][io][ipid]->SetMarkerStyle(kGreen+2);
				hcosndphiFinalP[iset][io][ipid]->SetMarkerColor(kGreen+2);
				hcosndphiFinalP[iset][io][ipid]->SetLineColor(kGreen+2);

				gvnInitP[iset][io][ipid] = new TGraphErrors;
				gvnInitP[iset][io][ipid]->SetLineColor(nColor[iset]);
				gvnInitP[iset][io][ipid]->SetLineWidth(0);
				gvnInitP[iset][io][ipid]->SetFillColorAlpha(nColor[iset],0.2);

				gvnFinalP[iset][io][ipid] = new TGraphErrors;
				gvnFinalP[iset][io][ipid]->SetLineColor(nColor[iset]);
				gvnFinalP[iset][io][ipid]->SetLineWidth(0);
				gvnFinalP[iset][io][ipid]->SetFillColorAlpha(nColor[iset],0.2);

				gvnFinalPH[iset][io][ipid] = new TGraphErrors;
				gvnFinalPH[iset][io][ipid]->SetLineColor(nColor[iset]);
				gvnFinalPH[iset][io][ipid]->SetLineWidth(0);
				gvnFinalPH[iset][io][ipid]->SetFillColorAlpha(nColor[iset],0.2);

				for (int ib=0; ib<hcosndphiInitP[iset][io][ipid]->GetNbinsX(); ib++){
					float xx = hcosndphiInitP[iset][io][ipid]->GetBinCenter(ib+1);
					float xx_err = 0.5*hcosndphiInitP[iset][io][ipid]->GetBinWidth(ib+1);
					float yy = hcosndphiInitP[iset][io][ipid]->GetBinContent(ib+1);
					float yy_err = hcosndphiInitP[iset][io][ipid]->GetBinError(ib+1);

					gvnInitP[iset][io][ipid]->SetPoint(ib, xx, yy);
					gvnInitP[iset][io][ipid]->SetPointError(ib, xx_err, 0.001+yy_err);
				}

				for (int ib=0; ib<hcosndphiFinalP[iset][io][ipid]->GetNbinsX(); ib++){
					float xx = hcosndphiFinalP[iset][io][ipid]->GetBinCenter(ib+1);
					float xx_err = 0.5*hcosndphiFinalP[iset][io][ipid]->GetBinWidth(ib+1);
					float yy = hcosndphiFinalP[iset][io][ipid]->GetBinContent(ib+1);
					float yy_err = hcosndphiFinalP[iset][io][ipid]->GetBinError(ib+1);

					gvnFinalP[iset][io][ipid]->SetPoint(ib, xx, yy);
					gvnFinalP[iset][io][ipid]->SetPointError(ib, xx_err, 0.001+yy_err);
				}

				for (int ib=0; ib<hcosndphiFinalPH[iset][io][ipid]->GetNbinsX(); ib++){
					float xx = hcosndphiFinalPH[iset][io][ipid]->GetBinCenter(ib+1);
					float xx_err = 0.5*hcosndphiFinalPH[iset][io][ipid]->GetBinWidth(ib+1);
					float yy = hcosndphiFinalPH[iset][io][ipid]->GetBinContent(ib+1);
					float yy_err = hcosndphiFinalPH[iset][io][ipid]->GetBinError(ib+1);

					gvnFinalPH[iset][io][ipid]->SetPoint(ib, xx, yy);
					gvnFinalPH[iset][io][ipid]->SetPointError(ib, xx_err, 0.001+yy_err);
				}

			}//ipid
		}//io
	}

	TCanvas *c0 = new TCanvas("c0","c0",1.1*3*500,500);
	c0->Divide(3,1);

	{
		c0->cd(1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,-0.05,10,0.2);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("v_{n}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitleOffset(1.15);
		htmp->GetYaxis()->SetLabelSize(0.04);

		for (int iset=0; iset<nset; iset++){
			gvnP[iset][0][0]->SetLineColorAlpha(nColor[iset],0.3);
			gvnP[iset][0][0]->SetFillColorAlpha(nColor[iset],0.3);
			gvnP[iset][0][0]->Draw("e3");
		}

		TLegend *leg = new TLegend(0.5,0.55,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","AMPT pPb 5.02 TeV","h");
		leg->AddEntry("","0-10% (-5<#eta<-3)","h");
		leg->AddEntry("","Charged hadron, |#eta|<1","h");
		leg->AddEntry("","Participant plane","h");
		leg->AddEntry(gvnP[0][0][0],"Parton on, hadron on","F");
		leg->AddEntry(gvnP[1][0][0],"Parton off, hadron on","F");
		leg->AddEntry(gvnP[2][0][0],"Parton on, hadron off","F");
		leg->AddEntry(gvnP[3][0][0],"Parton off, hadron off","F");
		leg->Draw();
	}

	{
		c0->cd(2);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,-0.05,10,0.2);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("v_{n}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitleOffset(1.15);
		htmp->GetYaxis()->SetLabelSize(0.04);

		for (int iset=0; iset<nset; iset++){
			gvnInitP[iset][0][0]->SetLineColorAlpha(nColor[iset],0.3);
			gvnInitP[iset][0][0]->SetFillColorAlpha(nColor[iset],0.3);
			gvnInitP[iset][0][0]->Draw("e3");
		}

		TLegend *leg = new TLegend(0.5,0.55,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","AMPT pPb 5.02 TeV","h");
		leg->AddEntry("","0-10% (-5<#eta<-3)","h");
		leg->AddEntry("","Initial parton (u, d), |#eta|<1","h");
		leg->AddEntry("","Participant plane","h");
		leg->AddEntry(gvnP[0][0][0],"Parton on, hadron on","F");
		leg->AddEntry(gvnP[1][0][0],"Parton off, hadron on","F");
		leg->AddEntry(gvnP[2][0][0],"Parton on, hadron off","F");
		leg->AddEntry(gvnP[3][0][0],"Parton off, hadron off","F");
		leg->Draw();
	}

	{
		c0->cd(3);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,-0.05,10,0.2);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("v_{n}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitleOffset(1.15);
		htmp->GetYaxis()->SetLabelSize(0.04);

		for (int iset=0; iset<nset; iset++){
			gvnFinalP[iset][0][0]->SetLineColorAlpha(nColor[iset],0.3);
			gvnFinalP[iset][0][0]->SetFillColorAlpha(nColor[iset],0.3);
			gvnFinalP[iset][0][0]->Draw("e3");
		}

		TLegend *leg = new TLegend(0.5,0.55,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.035);
		leg->AddEntry("","AMPT pPb 5.02 TeV","h");
		leg->AddEntry("","0-10% (-5<#eta<-3)","h");
		leg->AddEntry("","Final parton (u, d), |#eta|<1","h");
		leg->AddEntry("","Participant plane","h");
		leg->AddEntry(gvnP[0][0][0],"Parton on, hadron on","F");
		leg->AddEntry(gvnP[1][0][0],"Parton off, hadron on","F");
		leg->AddEntry(gvnP[2][0][0],"Parton on, hadron off","F");
		leg->AddEntry(gvnP[3][0][0],"Parton off, hadron off","F");
		leg->Draw();
	}


	TCanvas *c3 = new TCanvas("c3","c3",1.1*2*500,2*500);
	c3->Divide(2,2);

	for (int iset=0; iset<nset; iset++){

		c3->cd(iset+1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,-0.05,10,0.35);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("v_{n}");
		htmp->GetYaxis()->SetTitleOffset(1.15);
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		gvn[iset][1][0][1]->SetFillColorAlpha(1,0.3);
		gvn[iset][1][0][3]->SetFillColorAlpha(2,0.3);
		gvn[iset][1][0][1]->Draw("e3");
		gvn[iset][1][0][3]->Draw("e3");

		gvnP[iset][0][1]->SetFillColorAlpha(1,0.3);
		gvnP[iset][0][3]->SetFillColorAlpha(2,0.3);
		gvnP[iset][0][1]->SetFillStyle(3001);
		gvnP[iset][0][3]->SetFillStyle(3001);
		gvnP[iset][0][1]->Draw("e3");
		gvnP[iset][0][3]->Draw("e3");

		TLegend *leg = new TLegend(0.2,0.55,0.5,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.035);
		leg->AddEntry("","AMPT pPb 5.02 TeV","h");
		leg->AddEntry("","0-10% (-5<#eta<-3)","h");
		if ( iset==0 ){
			leg->AddEntry("","Parton on, hadron on","h");
		}else if ( iset==1 ){
			leg->AddEntry("","Parton off, hadron on","h");
		}else if ( iset==2 ){
			leg->AddEntry("","Parton on, hadron off","h");
		}else if ( iset==3 ){
			leg->AddEntry("","Parton off, hadron off","h");
		}
		leg->AddEntry("","Charged hadron, |#eta|<1","h");
		leg->AddEntry(gvn[iset][1][0][1],"#pi (EP, -3<#eta<-1)","F");
		leg->AddEntry(gvn[iset][1][0][3],"p (EP, -3<#eta<-1)","F");
		leg->AddEntry(gvnP[iset][0][1],"#pi (PP)","F");
		leg->AddEntry(gvnP[iset][0][3],"p (PP)","F");
		leg->Draw();
	}


	TCanvas *c4 = new TCanvas("c4","c4",1.1*2*500,2*500);
	c4->Divide(2,2);

	for (int iset=0; iset<nset; iset++){

		c4->cd(iset+1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,-0.01,10,0.15);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		htmp->GetYaxis()->SetTitle("v_{n}");
		htmp->GetYaxis()->SetTitleOffset(1.15);
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		gvnFinalPH[iset][0][0]->SetFillColorAlpha(1,0.3);
		gvnFinalPH[iset][0][1]->SetFillColorAlpha(2,0.3);
		gvnFinalPH[iset][0][0]->SetFillStyle(3001);
		gvnFinalPH[iset][0][1]->SetFillStyle(3001);
		gvnFinalPH[iset][0][0]->Draw("e3");
		gvnFinalPH[iset][0][1]->Draw("e3");

		TLegend *leg = new TLegend(0.2,0.55,0.5,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.035);
		leg->AddEntry("","AMPT pPb 5.02 TeV","h");
		leg->AddEntry("","0-10% (-5<#eta<-3)","h");
		leg->AddEntry("","Participant plane","h");
		if ( iset==0 ){
			leg->AddEntry("","Parton on, hadron on","h");
		}else if ( iset==1 ){
			leg->AddEntry("","Parton off, hadron on","h");
		}else if ( iset==2 ){
			leg->AddEntry("","Parton on, hadron off","h");
		}else if ( iset==3 ){
			leg->AddEntry("","Parton off, hadron off","h");
		}
		leg->AddEntry("","Final parton (u, d), |#eta|<1","h");
		leg->AddEntry(gvnFinalPH[iset][0][0],"u, d for #pi","F");
		leg->AddEntry(gvnFinalPH[iset][0][1],"u, d for p","F");
		leg->AddEntry("","","");
		leg->AddEntry("","","");
		leg->Draw();
	}

}

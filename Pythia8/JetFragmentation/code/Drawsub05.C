#include "Style.h"

void Drawsub05(const int opt = 0){

	const int nset = 3;
	const int nmult = 1;
	const int npt = 4;
	const int nz = 3;
	const int nColor[npt] = {1, 4, 8, 2};
	const int ptbin[npt+1] = {20, 40, 60, 80, 100};
	//const float ptbin[npt+1] = {10, 20, 40, 60, 100};
	const float scaler[npt] = {1, 10, 100, 1000};

	const float eta_edge = 0.25;

	//char *setname[nset] = {"New", "Old"};
	//char *setname[nset] = {"Default", "#tau_{0}^{max}=10 mm/c"};
	//char *setname[nset] = {"SoftQCD", "HardQCD"};
	//char *setname[nset] = {"Tune4C"};
	const char *setname[nset] = {"Set0", "Set1", "Set2"};

	TFile *infile[nset];
	//infile[0] = new TFile("../outfile_hist_pp5TeV_set34_grp000_try400.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[1] = new TFile("../outfile_hist_pp5TeV_set33_grp000_try400.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[2] = new TFile("../outfile_hist_pp5TeV_set30_grp000_try400.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	//infile[0] = new TFile("outfile_hist_pp5TeV_set34_grp000_try401.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[1] = new TFile("outfile_hist_pp5TeV_set33_grp000_try401.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[2] = new TFile("outfile_hist_pp5TeV_set30_grp000_try401.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	infile[0] = new TFile("../outfile_hist_pp5TeV_set34_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	infile[1] = new TFile("../outfile_hist_pp5TeV_set33_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	infile[2] = new TFile("../outfile_hist_pp5TeV_set30_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	TH1D *jeteta[nset][nmult][npt];
	TH1D *perpeta[nset][nmult][npt];

	TH1D *incljt[nset][nmult][npt];
	TH1D *perp_incljt[nset][nmult][npt];

	TH1D *zdepjt[nset][nmult][npt][3];
	TH1D *perp_zdepjt[nset][nmult][npt][3];

	TH1D *sub_incljt[nset][nmult][npt];
	TH1D *sub_zdepjt[nset][nmult][npt][3];

	TH1D *fbkg_zdepjt[nset][nmult][npt][3];

	TH1D *ratiojt[nset][nmult][npt];
	TH1D *ratiozdepjt[nset][nmult][npt][3];

	float njet[nset][nmult][npt] = {0.};
	float nperpjet[nset][nmult][npt] = {0.};

	for (int iset=0; iset<nset; iset++){
		for (int im=0; im<nmult; im++){
			for (int ii=0; ii<npt; ii++){
				jeteta[iset][im][ii] = (TH1D*)infile[iset]->Get(Form("jeteta_m%d_pt%d",im,ii));
				perpeta[iset][im][ii] = (TH1D*)infile[iset]->Get(Form("perpeta_m%d_pt%d",im,ii));

				int etabin_lo = jeteta[iset][im][ii]->FindBin(-eta_edge+0.001);
				int etabin_hi = jeteta[iset][im][ii]->FindBin(+eta_edge-0.001);

				njet[iset][im][ii] = jeteta[iset][im][ii]->Integral(etabin_lo, etabin_hi);
				nperpjet[iset][im][ii] = perpeta[iset][im][ii]->Integral(etabin_lo,etabin_hi);

				incljt[iset][im][ii] = (TH1D*)infile[iset]->Get(Form("incljt_m%d_pt%d",im,ii));
				incljt[iset][im][ii]->Sumw2();

				perp_incljt[iset][im][ii] = (TH1D*)infile[iset]->Get(Form("perp_incljt_m%d_pt%d",im,ii));
				perp_incljt[iset][im][ii]->Sumw2();

				for (int iz=0; iz<nz; iz++){
					zdepjt[iset][im][ii][iz] = (TH1D*)infile[iset]->Get(Form("zdepjt_m%d_pt%d_z%d",im,ii,iz));
					zdepjt[iset][im][ii][iz]->Sumw2();

					perp_zdepjt[iset][im][ii][iz] = (TH1D*)infile[iset]->Get(Form("perp_zdepjt_m%d_pt%d_z%d",im,ii,iz));
					perp_zdepjt[iset][im][ii][iz]->Sumw2();
				}//iz

			}//ii
		}//im
	}//iset

	//return;

	for (int iset=0; iset<nset; iset++){
		for (int im=0; im<nmult; im++){
			for (int ii=0; ii<npt; ii++){
				for (int ijt=0; ijt<incljt[iset][im][ii]->GetNbinsX(); ijt++){
					float jt = incljt[iset][im][ii]->GetBinCenter(ijt+1);
					float djt = incljt[iset][im][ii]->GetBinWidth(ijt+1);
					float xx = incljt[iset][im][ii]->GetBinContent(ijt+1);
					float xx_err = incljt[iset][im][ii]->GetBinError(ijt+1);

					incljt[iset][im][ii]->SetBinContent(ijt+1, xx/djt/jt/njet[iset][im][ii]);
					incljt[iset][im][ii]->SetBinError(ijt+1, xx_err/djt/jt/njet[iset][im][ii]);

					xx = perp_incljt[iset][im][ii]->GetBinContent(ijt+1);
					xx_err = perp_incljt[iset][im][ii]->GetBinError(ijt+1);

					perp_incljt[iset][im][ii]->SetBinContent(ijt+1, xx/djt/jt/nperpjet[iset][im][ii]);
					perp_incljt[iset][im][ii]->SetBinError(ijt+1, xx_err/djt/jt/nperpjet[iset][im][ii]);

					if ( jt<0.05 ){
						incljt[iset][im][ii]->SetBinContent(ijt+1, 0);
						incljt[iset][im][ii]->SetBinError(ijt+1, 0);

						perp_incljt[iset][im][ii]->SetBinContent(ijt+1, 0);
						perp_incljt[iset][im][ii]->SetBinError(ijt+1, 0);
					}
				}//ijt

				sub_incljt[iset][im][ii] = (TH1D*)incljt[iset][im][ii]->Clone(Form("sub_incljt_set%d_m%d_pt%d",iset,im,ii));
				sub_incljt[iset][im][ii]->Add(perp_incljt[iset][im][ii], -1);

				sub_incljt[iset][im][ii]->SetMarkerStyle(20+4*iset);
				sub_incljt[iset][im][ii]->SetMarkerColor(nColor[ii]);
				sub_incljt[iset][im][ii]->SetLineColor(nColor[ii]);
				//sub_incljt[iset][im][ii]->Scale(scaler[ii]);

				for (int iz=0; iz<nz; iz++){
					for (int ijt=0; ijt<zdepjt[iset][im][ii][iz]->GetNbinsX(); ijt++){
						float jt = zdepjt[iset][im][ii][iz]->GetBinCenter(ijt+1);
						float djt = zdepjt[iset][im][ii][iz]->GetBinWidth(ijt+1);
						float xx = zdepjt[iset][im][ii][iz]->GetBinContent(ijt+1);
						float xx_err = zdepjt[iset][im][ii][iz]->GetBinError(ijt+1);

						zdepjt[iset][im][ii][iz]->SetBinContent(ijt+1, xx/djt/jt/njet[iset][im][ii]);
						zdepjt[iset][im][ii][iz]->SetBinError(ijt+1, xx_err/djt/jt/njet[iset][im][ii]);

						xx = perp_zdepjt[iset][im][ii][iz]->GetBinContent(ijt+1);
						xx_err = perp_zdepjt[iset][im][ii][iz]->GetBinError(ijt+1);

						perp_zdepjt[iset][im][ii][iz]->SetBinContent(ijt+1, xx/djt/jt/nperpjet[iset][im][ii]);
						perp_zdepjt[iset][im][ii][iz]->SetBinError(ijt+1, xx_err/djt/jt/nperpjet[iset][im][ii]);

						if ( jt<0.05 ){
							zdepjt[iset][im][ii][iz]->SetBinContent(ijt+1, 0);
							zdepjt[iset][im][ii][iz]->SetBinError(ijt+1, 0);

							perp_zdepjt[iset][im][ii][iz]->SetBinContent(ijt+1, 0);
							perp_zdepjt[iset][im][ii][iz]->SetBinError(ijt+1, 0);
						}
					}//ijt

					fbkg_zdepjt[iset][im][ii][iz] = (TH1D*)perp_zdepjt[iset][im][ii][iz]->Clone(Form("fbkg_zdepjt_set%d_m%d_pt%d_z%d",iset,im,ii,iz));
					fbkg_zdepjt[iset][im][ii][iz]->Divide(zdepjt[iset][im][ii][iz]);

					fbkg_zdepjt[iset][im][ii][iz]->SetMarkerStyle(20+4*iset);
					fbkg_zdepjt[iset][im][ii][iz]->SetMarkerColor(nColor[ii]);
					fbkg_zdepjt[iset][im][ii][iz]->SetLineColor(nColor[ii]);

					sub_zdepjt[iset][im][ii][iz] = (TH1D*)zdepjt[iset][im][ii][iz]->Clone(Form("sub_zdepjt_set%d_m%d_pt%d_z%d",iset,im,ii,iz));
					sub_zdepjt[iset][im][ii][iz]->Add(perp_zdepjt[iset][im][ii][iz], -1);

					sub_zdepjt[iset][im][ii][iz]->SetMarkerStyle(20+4*iset);
					sub_zdepjt[iset][im][ii][iz]->SetMarkerColor(nColor[ii]);
					sub_zdepjt[iset][im][ii][iz]->SetLineColor(nColor[ii]);
					sub_zdepjt[iset][im][ii][iz]->Scale(scaler[ii]);

				}//iz

			}//ii
		}//im
	}//iset

	//return;

	TCanvas *c1[nmult];
	TCanvas *c2[nmult][nz];
	TCanvas *c3[nmult];

	for (int im=0; im<nmult; im++){

		c1[im] = new TCanvas(Form("c1_m%d",im),Form("c1_m%d",im),600,700);
		TPad *p1_1 = new TPad("p1_1","p1_1",0,0.4,1,1);
		p1_1->Draw();
		p1_1->cd();
		SetPadStyle();
		gPad->SetBottomMargin(0);
		gPad->SetLeftMargin(0.13);
		p1_1->SetLogy(1);
		p1_1->SetLogx(1);

		htmp = (TH1D*)gPad->DrawFrame(0.052, 5e-5, 3.02, 5e6);
		//htmp = (TH1D*)gPad->DrawFrame(0.052, 5e-5, 3.02, 150);
		SetHistoStyle("","1/#it{N}_{jet} 1/#it{j_{T}} d#it{N}_{ch}^{ jet}/d#it{j_{T}} (GeV/#it{c})^{-1}","",20,18);
		htmp->GetYaxis()->SetTitleOffset(1.6);

		for (int iset=0; iset<nset; iset++){
			for (int ii=0; ii<npt; ii++){
				sub_incljt[iset][im][ii]->Draw("p same");
			}
		}

		{
			TLegend *leg = new TLegend(0.165,0.05,0.71,0.05+0.35); 
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.04);
			leg->SetNColumns(4);
			leg->AddEntry("","PYTHIA8 pp 5 TeV","h");
			leg->AddEntry("","Inclusive z","h");
			leg->AddEntry("",Form("%d<#it{p_{T}}<%d GeV/#it{c} #times10^{0}",ptbin[0],ptbin[1]),"");
			leg->AddEntry(sub_incljt[0][im][0],setname[0],"P");
			leg->AddEntry(sub_incljt[1][im][0],setname[1],"P");
			leg->AddEntry(sub_incljt[2][im][0],setname[2],"P");
			leg->AddEntry("",Form("%d<#it{p_{T}}<%d GeV/#it{c} #times10^{0}",ptbin[1],ptbin[2]),"");
			leg->AddEntry(sub_incljt[0][im][1],setname[0],"P");
			leg->AddEntry(sub_incljt[1][im][1],setname[1],"P");
			leg->AddEntry(sub_incljt[2][im][1],setname[2],"P");
			leg->AddEntry("",Form("%d<#it{p_{T}}<%d GeV/#it{c} #times10^{0}",ptbin[2],ptbin[3]),"");
			leg->AddEntry(sub_incljt[0][im][2],setname[0],"P");
			leg->AddEntry(sub_incljt[1][im][2],setname[1],"P");
			leg->AddEntry(sub_incljt[2][im][2],setname[2],"P");
			leg->AddEntry("",Form("%d<#it{p_{T}}<%d GeV/#it{c} #times10^{0}",ptbin[3],ptbin[4]),"");
			leg->AddEntry(sub_incljt[0][im][3],setname[0],"P");
			leg->AddEntry(sub_incljt[1][im][3],setname[1],"P");
			leg->AddEntry(sub_incljt[2][im][3],setname[2],"P");
			leg->Draw();
		}

		c1[im]->cd();
		TPad *p1_2 = new TPad("p1_2","p1_2",0,0.0,1,0.4);
		p1_2->Draw();
		p1_2->cd();
		SetPadStyle();
		gPad->SetTopMargin(0);
		gPad->SetBottomMargin(0.2);
		gPad->SetLeftMargin(0.13);
		p1_2->SetLogx(1);

		htmp = (TH1D*)gPad->DrawFrame(0.052, 0.2, 3.02, 1.8);
		SetHistoStyle("#it{j_{T}} (GeV/#it{c})","Ratio","",20,18);
		htmp->GetYaxis()->SetTitleOffset(1.6);
		htmp->GetYaxis()->SetNdivisions(7,5,0);
		htmp->GetYaxis()->CenterTitle();
		htmp->GetXaxis()->SetTitleOffset(1.3);

		for (int iset=1; iset<nset; iset++){
			for (int ii=0; ii<npt; ii++){
				ratiojt[iset][im][ii] = (TH1D*)sub_incljt[iset][im][ii]->Clone(Form("ratiojt_%d_m%d_pt%d",iset,im,ii));
				ratiojt[iset][im][ii]->Divide(sub_incljt[0][im][ii]);
				ratiojt[iset][im][ii]->Draw("p same");
			}//ii
		}//iset

		{
			TLegend *leg = new TLegend(0.6,0.7,0.95,0.94); 
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.07);
			//leg->AddEntry("",Form("#frac{%s}{%s}",setname[1],setname[0]),"h");
			leg->AddEntry("","Ratio to Set0","h");
			leg->Draw();
		}

		//continue;

		//for (int iz=0; iz<nz; iz++){
		for (int iz=0; iz<0; iz++){

			if ( iz==1 ) continue;

			c2[im][iz] = new TCanvas(Form("c2_m%d_z%d",im,iz),Form("c2_m%d_z%d",im,iz),600,700);
			TPad *p1_1 = new TPad("p1_1","p1_1",0,0.4,1,1);
			p1_1->Draw();
			p1_1->cd();
			SetPadStyle();
			gPad->SetBottomMargin(0);
			gPad->SetLeftMargin(0.13);
			p1_1->SetLogy(1);
			p1_1->SetLogx(1);

			htmp = (TH1D*)gPad->DrawFrame(0.052, 5e-5, 3.02, 5e6);
			SetHistoStyle("","1/#it{N}_{jet} 1/#it{j_{T}} d#it{N}_{ch}^{ jet}/d#it{j_{T}} (GeV/#it{c})^{-1}","",20,18);
			htmp->GetYaxis()->SetTitleOffset(1.6);

			for (int iset=0; iset<nset; iset++){
				for (int ii=0; ii<npt; ii++){
					sub_zdepjt[iset][im][ii][iz]->Draw("p same");
				}
			}

			{
				TLegend *leg = new TLegend(0.165,0.05,0.71,0.05+0.35); 
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.04);
				leg->SetNColumns(3);
				leg->AddEntry("","PYTHIA8 pp 5 TeV","h");
				if ( iz==0 ){
					leg->AddEntry("","0<z<0.2","h");
				}else if ( iz==1 ){
				}else{
					leg->AddEntry("","0.4<z#leq1","h");
				}
				leg->AddEntry("","10<#it{p_{T}}<20 GeV/#it{c} #times10^{0}","");
				leg->AddEntry(sub_incljt[0][im][0],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][0],setname[1],"P");
				leg->AddEntry("","20<#it{p_{T}}<40 GeV/#it{c} #times10^{1}","");
				leg->AddEntry(sub_incljt[0][im][1],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][1],setname[1],"P");
				leg->AddEntry("","40<#it{p_{T}}<60 GeV/#it{c} #times10^{2}","");
				leg->AddEntry(sub_incljt[0][im][2],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][2],setname[1],"P");
				leg->AddEntry("","60<#it{p_{T}}<100 GeV/#it{c} #times10^{3}","");
				leg->AddEntry(sub_incljt[0][im][3],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][3],setname[1],"P");
				leg->Draw();
			}

			c2[im][iz]->cd();
			TPad *p1_2 = new TPad("p1_2","p1_2",0,0.0,1,0.4);
			p1_2->Draw();
			p1_2->cd();
			SetPadStyle();
			gPad->SetTopMargin(0);
			gPad->SetBottomMargin(0.2);
			gPad->SetLeftMargin(0.13);
			p1_2->SetLogx(1);

			htmp = (TH1D*)gPad->DrawFrame(0.052, 0.0, 3.02, 1.3);
			SetHistoStyle("#it{j_{T}} (GeV/#it{c})","Ratio to inclusive z","",20,18);
			htmp->GetYaxis()->SetTitleOffset(1.6);
			htmp->GetYaxis()->SetNdivisions(7,5,0);
			htmp->GetYaxis()->CenterTitle();
			htmp->GetXaxis()->SetTitleOffset(3.0);

			for (int iset=0; iset<nset; iset++){
				for (int ii=0; ii<npt; ii++){
					ratiozdepjt[iset][im][ii][iz] = (TH1D*)sub_zdepjt[iset][im][ii][iz]->Clone(Form("ratiozdepjt_set%d_m%d_pt%d_z%d",iset,im,ii,iz));
					ratiozdepjt[iset][im][ii][iz]->Divide(sub_incljt[iset][im][ii]);
					ratiozdepjt[iset][im][ii][iz]->SetMarkerStyle(20+4*iset);
					ratiozdepjt[iset][im][ii][iz]->Draw("p same");
				}//ii
			}

		}//iz

		c3[im] = new TCanvas(Form("c3_m%d",im),Form("c3_m%d",im),3*1.1*500,500);
		c3[im]->Divide(3,1);

		for (int iz=0; iz<nz; iz++){

			c3[im]->cd(iz+1);
			SetPadStyle();
			gPad->SetLeftMargin(0.12);
			gPad->SetLogx();
			gPad->SetLogy();

			htmp = (TH1D*)gPad->DrawFrame(0.052,0.001,3.02,1);
			SetHistoStyle("#it{j_{T}} (GeV/#it{c})","Background fraction","",20,18);
			htmp->GetYaxis()->SetTitleOffset(1.4);
			htmp->GetYaxis()->SetNdivisions(7,5,0);
			htmp->GetYaxis()->CenterTitle();
			htmp->GetXaxis()->SetTitleOffset(1.2);

			for (int iset=0; iset<nset; iset++){
				for (int ii=0; ii<npt; ii++){
					fbkg_zdepjt[iset][0][ii][iz]->Draw("p same");
				}
			}

			if ( iz==0 )
			{
				TLegend *leg = new TLegend(0.165,0.15,0.91,0.15+0.3); 
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.04);
				leg->SetNColumns(3);
				leg->AddEntry("","PYTHIA8 pp 5 TeV","h");
				leg->AddEntry("","10<#it{p_{T}}<20 GeV/#it{c} #times10^{0}","");
				leg->AddEntry(sub_incljt[0][im][0],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][0],setname[1],"P");
				leg->AddEntry("","20<#it{p_{T}}<40 GeV/#it{c} #times10^{1}","");
				leg->AddEntry(sub_incljt[0][im][1],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][1],setname[1],"P");
				leg->AddEntry("","40<#it{p_{T}}<60 GeV/#it{c} #times10^{2}","");
				leg->AddEntry(sub_incljt[0][im][2],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][2],setname[1],"P");
				leg->AddEntry("","60<#it{p_{T}}<100 GeV/#it{c} #times10^{3}","");
				leg->AddEntry(sub_incljt[0][im][3],setname[0],"P");
				leg->AddEntry(sub_incljt[1][im][3],setname[1],"P");
				leg->Draw();
			}

			{
				TLegend *leg = new TLegend(0.30,0.85,0.6,0.90); 
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.04);
				if ( iz==0 ){
					leg->AddEntry("","0<z<0.2","h");
				}else if ( iz==1 ){
					leg->AddEntry("","0.2<z<0.4","h");
				}else{
					leg->AddEntry("","0.4<z#leq1","h");
				}
				leg->Draw();
			}

		}//iz

	}//im

	/*
	*/

	//return;

	TFile *outfile = new TFile("PYTHIA8_pp5TeV_inclz_jt_fulljet_jtbin2.root", "RECREATE");
	for (int iset=0; iset<nset; iset++){
		for (int ipt=0; ipt<npt; ipt++){
			sub_incljt[iset][0][ipt]->Write();
		}

		/*
		for (int ipt=0; ipt<npt; ipt++){
			for (int iz=0; iz<nz; iz++){
				sub_zdepjt[iset][0][ipt][iz]->Write();
			}
		}
		*/
	}
	outfile-> Close();

	/*
	TH1D *ratioz = (TH1D*)sub_inclz->Clone("ratioz");
	ratioz->Divide(sub_sigz);
	ratioz->Draw("p same");
	*/

	return;

	//delete infile;
}

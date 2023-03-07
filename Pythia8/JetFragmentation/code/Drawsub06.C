#include "Style.h"

void Drawsub06(const int opt = 0){

	const int nset = 2;
	const int nmult = 1;
	const int npt = 2;
	const int nz = 3;
	const int nColor[4] = {1, 4, 8, 2};
	//const int ptbin[npt+1] = {20, 40, 60, 80, 100};
	//const float ptbin[npt+1] = {10, 20, 40, 60, 100};
	//const float scaler[npt] = {1, 10, 100, 1000};

	const int ptbin[npt+1] = {40, 60, 100};
	const float scaler[npt] = {1, 10};

	const float eta_edge = 0.50;

	//char *setname[nset] = {"New", "Old"};
	//char *setname[nset] = {"Default", "#tau_{0}^{max}=10 mm/c"};
	//char *setname[nset] = {"SoftQCD", "HardQCD"};
	//char *setname[nset] = {"Tune4C"};
	//const char *setname[nset] = {"Set0", "Set1", "Set2"};
	const char *setname[nset] = {"5 GeV","40 GeV"};

	TFile *infile[nset];
	//infile[0] = new TFile("../outfile_hist_pp5TeV_set34_grp000_try400.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[1] = new TFile("../outfile_hist_pp5TeV_set33_grp000_try400.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[2] = new TFile("../outfile_hist_pp5TeV_set30_grp000_try400.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	//infile[0] = new TFile("outfile_hist_pp5TeV_set34_grp000_try401.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[1] = new TFile("outfile_hist_pp5TeV_set33_grp000_try401.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[2] = new TFile("outfile_hist_pp5TeV_set30_grp000_try401.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	//infile[0] = new TFile("../outfile_hist_pp5TeV_set34_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[1] = new TFile("../outfile_hist_pp5TeV_set33_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[2] = new TFile("../outfile_hist_pp5TeV_set30_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	infile[0] = new TFile("../outfile_hist_set33_grp000_v1.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	infile[1] = new TFile("../outfile_hist_set34_grp000_v1.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 
	//infile[2] = new TFile("../outfile_hist_pp5TeV_set30_grp000_try402.root","read"); //Monash, SoftQCD, pt cut for jet reco, tau0max cut 

	TH1D *jeteta[nset][npt];
	TH1D *perpeta[nset][npt];

	TH1D *incljt[nset][nmult][npt];
	TH1D *perp_incljt[nset][nmult][npt];

	TH1D *zdepjt[nset][nmult][npt][3];
	TH1D *perp_zdepjt[nset][nmult][npt][3];

	TH1D *sub_incljt[nset][nmult][npt];
	TH1D *sub_zdepjt[nset][nmult][npt][3];

	TH1D *fbkg_zdepjt[nset][nmult][npt][3];

	TH1D *ratiozdepjt[nset][nmult][npt][3];

	TH2D *hjet_eta_pt[nset];
	TH2D *hjet_jt_pt[nset];

	TH1D *ratiojt[nset][npt];
	TH1D *hjet_jt[nset][npt];

	float njet[nset][npt] = {0.};

	for (int iset=0; iset<nset; iset++){

		hjet_eta_pt[iset] = (TH2D*)infile[iset]->Get("hjetR04_eta_pt");
		hjet_jt_pt[iset] = (TH2D*)infile[iset]->Get("hjetR04_jt_pt");

		TH1D *hjet_eta = (TH1D*)hjet_eta_pt[iset]->ProjectionX(Form("hjet_eta_%d",iset));

		int etabin_lo = hjet_eta->FindBin(-eta_edge+0.001);
		int etabin_hi = hjet_eta->FindBin(+eta_edge-0.001);

		TH1D *hjet_pt = (TH1D*)hjet_eta_pt[iset]->ProjectionY(Form("hjet_pt_%d",iset),etabin_lo,etabin_hi);
		TH1D *hjet_pt2 = (TH1D*)hjet_jt_pt[iset]->ProjectionY(Form("hjet_pt2_%d",iset));

		for (int ii=0; ii<npt; ii++){

			int ptbin_lo = hjet_pt->FindBin(ptbin[ii]+0.001);
			int ptbin_hi = hjet_pt->FindBin(ptbin[ii+1]-0.001);

			njet[iset][ii] = hjet_pt->Integral(ptbin_lo, ptbin_hi);
			cout << "set: " << iset << ", njet: " << njet[iset][ii] << endl;

			ptbin_lo = hjet_pt2->FindBin(ptbin[ii]+0.001);
			ptbin_hi = hjet_pt2->FindBin(ptbin[ii+1]-0.001);

			hjet_jt[iset][ii] = (TH1D*)hjet_jt_pt[iset]->ProjectionX(Form("hjet_jt_%d_%d",iset,ii),ptbin_lo, ptbin_hi);
			hjet_jt[iset][ii]->Sumw2();

			hjet_jt[iset][ii]->SetMarkerStyle(20+4*iset);
			hjet_jt[iset][ii]->SetMarkerColor(nColor[ii]);
			hjet_jt[iset][ii]->SetLineColor(nColor[ii]);

		}//ii
		
	}//iset

	//return;
	for (int iset=0; iset<nset; iset++){
		for (int ii=0; ii<npt; ii++){
			for (int ijt=0; ijt<hjet_jt[iset][ii]->GetNbinsX(); ijt++){
				float jt = hjet_jt[iset][ii]->GetBinCenter(ijt+1);
				float djt = hjet_jt[iset][ii]->GetBinWidth(ijt+1);
				float xx = hjet_jt[iset][ii]->GetBinContent(ijt+1);
				float xx_err = hjet_jt[iset][ii]->GetBinError(ijt+1);

				hjet_jt[iset][ii]->SetBinContent(ijt+1, xx/djt/jt/njet[iset][ii]);
				hjet_jt[iset][ii]->SetBinError(ijt+1, xx_err/djt/jt/njet[iset][ii]);

				if ( jt<0.05 ){
					hjet_jt[iset][ii]->SetBinContent(ijt+1, 0);
					hjet_jt[iset][ii]->SetBinError(ijt+1, 0);
				}
			}//ijt

		}//ii
	}//iset

	TCanvas *c1 = new TCanvas(Form("c1"),Form("c1"),600,700);

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
			hjet_jt[iset][ii]->Draw("p same");
		}
	}

	{
		TLegend *leg = new TLegend(0.165,0.05,0.71,0.05+0.2); 
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->SetNColumns(3);
		leg->AddEntry("","PYTHIA8 pp 5 TeV","h");
		leg->AddEntry("","Inclusive z","h");
		leg->AddEntry("",Form("%d<#it{p_{T}}<%d GeV/#it{c} #times10^{0}",ptbin[0],ptbin[1]),"");
		leg->AddEntry(hjet_jt[0][0],setname[0],"P");
		leg->AddEntry(hjet_jt[1][0],setname[1],"P");
		leg->AddEntry("",Form("%d<#it{p_{T}}<%d GeV/#it{c} #times10^{1}",ptbin[1],ptbin[2]),"");
		leg->AddEntry(hjet_jt[0][1],setname[0],"P");
		leg->AddEntry(hjet_jt[1][1],setname[1],"P");
		leg->Draw();
	}

	c1->cd();

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
			ratiojt[iset][ii] = (TH1D*)hjet_jt[iset][ii]->Clone(Form("ratiojt_%d_pt%d",iset,ii));
			ratiojt[iset][ii]->Divide(hjet_jt[0][ii]);
			ratiojt[iset][ii]->Draw("p same");
		}//ii
	}//iset

	return;

	TCanvas *c2[nmult][nz];
	TCanvas *c3[nmult];

	for (int im=0; im<nmult; im++){



		{
			TLegend *leg = new TLegend(0.6,0.7,0.95,0.94); 
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.07);
			//leg->AddEntry("",Form("#frac{%s}{%s}",setname[1],setname[0]),"h");
			leg->AddEntry("","Ratio to Set0","h");
			leg->Draw();
		}

		continue;

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

	return;

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

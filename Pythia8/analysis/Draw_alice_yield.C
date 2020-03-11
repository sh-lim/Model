void Draw_alice_yield(){


	gStyle->SetOptStat(0);

	//TFile *infile = new TFile("../ALICE-kinematics/outfile_hist.root","read");
	//TFile *infile = new TFile("../ALICE-kinematics/outfile_hist_set01_grp000.root","read");
	//TFile *infile = new TFile("../ALICE-kinematics/outfile_hist_set01_grp003.root","read");
	TFile *infile = new TFile("../ALICE-kinematics/outfile_hist_set01_grp003-0-10.root","read");

	const int npt = 6;
	const float pt_cut[npt] = {10, 20, 30, 40, 50, 60};

	TH2D *h2d_same[npt];
	TH2D *h2d_mixed[npt];

	TH1D *hntrig_same = (TH1D*)infile->Get("hntrig_same_mult00");
	TH1D *hntrig_mixed = (TH1D*)infile->Get("hntrig_mixed_mult00");

	TH1D *h1d_deta_same[npt];
	TH1D *h1d_deta_mixed[npt];

	TH1D *h1d_dphi_same[npt];
	TH1D *h1d_dphi_mixed[npt];

	for (int ipt=0; ipt<npt; ipt++){
		h2d_same[ipt] = (TH2D*)infile->Get(Form("h2d_same_dphi_deta_mult00_ptlead%02d",ipt));
		h2d_mixed[ipt] = (TH2D*)infile->Get(Form("h2d_mixed_dphi_deta_mult00_ptlead%02d",ipt));

		h2d_same[ipt]->RebinX(4);
		h2d_mixed[ipt]->RebinX(4);

		float ntrig_same = hntrig_same->Integral(hntrig_same->FindBin(pt_cut[ipt]+0.1), hntrig_same->GetNbinsX());
		float ntrig_mixed = hntrig_mixed->Integral(hntrig_mixed->FindBin(pt_cut[ipt]+0.1), hntrig_mixed->GetNbinsX());
		float nnorm_mixed = h2d_mixed[ipt]->GetBinContent(h2d_mixed[ipt]->FindBin(0,0));
		nnorm_mixed /= ntrig_mixed;
		nnorm_mixed /= h2d_mixed[ipt]->GetXaxis()->GetBinWidth(1);
		nnorm_mixed /= h2d_mixed[ipt]->GetYaxis()->GetBinWidth(1);

		//cout << nnorm_mixed << endl;
		//cout << ntrig_same << " " << ntrig_mixed << endl;

		//1D projection to near-side 
		int phibin_min = h2d_same[ipt]->GetXaxis()->FindBin(-TMath::Pi()/6.0);
		int phibin_max = h2d_same[ipt]->GetXaxis()->FindBin(+TMath::Pi()/6.0);

		h1d_deta_same[ipt] = (TH1D*)h2d_same[ipt]->ProjectionY(Form("h1d_deta_same_pt%02d",ipt),phibin_min,phibin_max);
		h1d_deta_mixed[ipt] = (TH1D*)h2d_mixed[ipt]->ProjectionY(Form("h1d_deta_mixed_pt%02d",ipt),phibin_min,phibin_max);

		//1D projection to long-side 
		int etabin_min = h2d_same[ipt]->GetYaxis()->FindBin(-1.8+0.01);
		int etabin_max = h2d_same[ipt]->GetYaxis()->FindBin(-1.5-0.01);

		h1d_dphi_same[ipt] = (TH1D*)h2d_same[ipt]->ProjectionX(Form("h1d_dphi_same_pt%02d",ipt),etabin_min,etabin_max);
		h1d_dphi_mixed[ipt] = (TH1D*)h2d_mixed[ipt]->ProjectionX(Form("h1d_dphi_mixed_pt%02d",ipt),etabin_min,etabin_max);

		etabin_min = h2d_same[ipt]->GetYaxis()->FindBin(1.5+0.01);
		etabin_max = h2d_same[ipt]->GetYaxis()->FindBin(1.8-0.01);

		TH1D *htmp_same = (TH1D*)h2d_same[ipt]->ProjectionX(Form("h1d_dphi_same_pt%02d",ipt),etabin_min,etabin_max);
		TH1D *htmp_mixed = (TH1D*)h2d_mixed[ipt]->ProjectionX(Form("h1d_dphi_mixed_pt%02d",ipt),etabin_min,etabin_max);

		h1d_dphi_same[ipt]->Add(htmp_same);
		h1d_dphi_mixed[ipt]->Add(htmp_mixed);


		//normalization
		//h2d_same[ipt]->RebinY(2);
		//h2d_mixed[ipt]->RebinY(2);

		h2d_same[ipt]->Scale(1./ntrig_same);
		h2d_mixed[ipt]->Scale(1./ntrig_mixed);
		h2d_same[ipt]->Divide(h2d_mixed[ipt]);
		h2d_same[ipt]->Scale(nnorm_mixed);

		h1d_deta_same[ipt]->Scale(1./ntrig_same);
		h1d_deta_mixed[ipt]->Scale(1./ntrig_mixed);
		h1d_deta_same[ipt]->Divide(h1d_deta_mixed[ipt]);

		h1d_dphi_same[ipt]->Scale(1./ntrig_same);
		h1d_dphi_mixed[ipt]->Scale(1./ntrig_mixed);
		h1d_dphi_same[ipt]->Divide(h1d_dphi_mixed[ipt]);
		h1d_dphi_same[ipt]->Scale(nnorm_mixed);

		//ZYAM subtraction

		float associated_yield = 0.0;
		float zyam = h1d_dphi_same[ipt]->GetBinContent(h1d_dphi_same[ipt]->GetMinimumBin());
		for (int iphi=0; iphi<h1d_dphi_same[ipt]->GetNbinsX(); iphi++){
			h1d_dphi_same[ipt]->SetBinContent(iphi+1, h1d_dphi_same[ipt]->GetBinContent(iphi+1)-zyam);

			if ( fabs(h1d_dphi_same[ipt]->GetBinCenter(iphi+1))<TMath::Pi()/6 ){
				associated_yield += h1d_dphi_same[ipt]->GetBinContent(iphi+1);
			}
		}

		cout << "AY: " << associated_yield << endl;

	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*300*3,300*2);
	c1->Divide(3,2);

	for (int ipt=0; ipt<npt; ipt++){
		c1->cd(ipt+1);
		gPad->SetPhi(135);

		h2d_same[ipt]->SetAxisRange(-1.8+0.01,1.8-0.01,"Y");
		h2d_same[ipt]->GetYaxis()->SetTitle("|#Delta#eta|");
		h2d_same[ipt]->GetYaxis()->CenterTitle();
		h2d_same[ipt]->GetYaxis()->SetTitleOffset(1.8);
		h2d_same[ipt]->GetXaxis()->SetTitle("|#Delta#phi|");
		h2d_same[ipt]->GetXaxis()->CenterTitle();
		h2d_same[ipt]->GetXaxis()->SetTitleOffset(1.8);
		//h2d_same[ipt]->SetMaximum(0.6*h2d_same[ipt]->GetMaximum());
		h2d_same[ipt]->Draw("surf1");

		TLegend *leg = new TLegend(0.05,0.75,0.5,0.99);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");
		leg->AddEntry("","HardQCD:all","h");
		leg->AddEntry("","MPI on","h");
		leg->AddEntry("",Form("p_{T, jet}>%g GeV/c",pt_cut[ipt]),"h");
		leg->AddEntry("","1<p_{T, trig}<2 GeV/c","h");
		leg->AddEntry("","1<p_{T, assoc}<2 GeV/c","h");
		leg->Draw();
	}

	TCanvas *c2 = new TCanvas("c2","c2",1.1*300*3,300*2);
	c2->Divide(3,2);

	for (int ipt=0; ipt<npt; ipt++){
		c2->cd(ipt+1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.05);
		gPad->SetTopMargin(0.05);

		//TF1 *f1 = new TF1("f1","[0]*( 1 + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x))",-TMath::Pi()/2,3*TMath::Pi()/2);
		//h1d_dphi_same[ipt]->Fit(f1,"R0Q");

		h1d_dphi_same[ipt]->GetXaxis()->SetTitle("|#Delta#phi|");
		h1d_dphi_same[ipt]->GetYaxis()->SetTitle("(1/N_{trig})(d^{2}N/d#Delta#phid#Delta#eta) - C_{ZYAM}");
		h1d_dphi_same[ipt]->Draw();

		TLegend *leg = new TLegend(0.15,0.63,0.6,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");
		leg->AddEntry("","HardQCD:all","h");
		leg->AddEntry("","MPI on, 0%-10%","h");
		leg->AddEntry("",Form("p_{T, jet}>%g GeV/c",pt_cut[ipt]),"h");
		leg->AddEntry("","1<p_{T, track}<2 GeV/c","h");
		leg->AddEntry("","1.5<|#Delta#eta|<1.8","h");
		leg->Draw();

		TBox *box = new TBox(-TMath::Pi()/6, h1d_dphi_same[ipt]->GetMinimum(), TMath::Pi()/6, h1d_dphi_same[ipt]->GetMaximum());
		box->SetFillColorAlpha(1,0.3);
		//box->Draw();

		//f1->Draw("same");
	}


}

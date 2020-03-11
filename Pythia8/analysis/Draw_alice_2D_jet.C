void Draw_alice_2D_jet(){


	gStyle->SetOptStat(0);

	TFile *infile = new TFile("../ALICE-kinematics/outfile_hist_set01_grp000.root","read");

	TH2D *h2d_same = (TH2D*)infile->Get("h2d_same_dphi_deta_mult00_ptlead00");
	TH2D *h2d_mixed = (TH2D*)infile->Get("h2d_mixed_dphi_deta_mult00_ptlead00");

	h2d_same->RebinX(2);
	h2d_mixed->RebinX(2);

	h2d_same->RebinY(4);
	h2d_mixed->RebinY(4);

	h2d_same->Scale(1./h2d_same->Integral());
	h2d_mixed->Scale(1./h2d_mixed->Integral());

	h2d_same->Divide(h2d_mixed);
	h2d_same->SetAxisRange(-1.9,1.9,"Y");
	h2d_same->GetYaxis()->SetTitle("|#Delta#eta|");
	h2d_same->GetYaxis()->CenterTitle();
	h2d_same->GetYaxis()->SetTitleOffset(1.8);
	h2d_same->GetXaxis()->SetTitle("|#Delta#phi|");
	h2d_same->GetXaxis()->CenterTitle();
	h2d_same->GetXaxis()->SetTitleOffset(1.8);
	h2d_same->SetMaximum(0.6*h2d_same->GetMaximum());

	TCanvas *c1 = new TCanvas("c1","c1",1.1*500,500);
	gPad->SetPhi(135);
	h2d_same->Draw("surf1");

	TLegend *leg = new TLegend(0.05,0.75,0.5,0.99);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("","Pythia8 pp 13 TeV","h");
	leg->AddEntry("","HardQCD:all","h");
	leg->AddEntry("","MPI,ISR,FSR off","h");
	leg->AddEntry("","p_{T}^{jet}>10 GeV/c, |#eta^{jet}|<0.4","h");
	leg->AddEntry("","1<p_{T}^{Trig, Assoc}<2 GeV/c","h");
	leg->Draw();

}

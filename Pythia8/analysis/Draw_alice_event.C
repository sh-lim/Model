#include "Style.h"

void Draw_alice_event(){

	const bool bSAVE = false;

	gStyle->SetOptStat(0);
	//gStyle->SetPalette(55);

	const int nmult = 5;
	const int cut_mult[nmult+1] = {100, 60, 20, 5, 1, 0};

	TFile *infile = new TFile("../ALICE-kinematics/outfile_hist_pp13TeV_set00_grp000_try001.root","read");
	//TFile *infile = new TFile("../ALICE-kinematics/outfile_hist.root","read");

	TH2D *hevent_Q[nmult];
	TH2D *hevent_bMPI[nmult];
	TH2D *hevent_nMPI[nmult];

	TH2D *hevent_Q_bMPI[nmult];
	TH2D *hevent_Q_nMPI[nmult];

	for (int im=0; im<nmult; im++){

		hevent_Q[im] = (TH2D*)infile->Get(Form("hevent_Q_mult%02d",im));
		hevent_bMPI[im] = (TH2D*)infile->Get(Form("hevent_bMPI_mult%02d",im));
		hevent_nMPI[im] = (TH2D*)infile->Get(Form("hevent_nMPI_mult%02d",im));

		hevent_Q_bMPI[im] = (TH2D*)infile->Get(Form("hevent_Q_bMPI_mult%02d",im));
		hevent_Q_nMPI[im] = (TH2D*)infile->Get(Form("hevent_Q_nMPI_mult%02d",im));
	}

	TH2D *hevent_mult_mid_fwd = (TH2D*)infile->Get("hevent_mult_mid_fwd");

	TCanvas *c100 = new TCanvas("c100","c100",1.1*400,400);
	gPad->SetLogz();
	SetPadStyle(1);

	htmp = (TH1D*)gPad->DrawFrame(0,0,50,50);
	SetHistoStyle("N_{ch} in |#eta|<0.9, p_{T}>0.2 GeV/c","N_{ch} in V0","",18,18);
	htmp->GetYaxis()->SetTitleOffset(1.4);
	htmp->GetXaxis()->SetTitleOffset(1.1);
	hevent_mult_mid_fwd->Draw("colz same");

	//return;

	TCanvas *c0 = new TCanvas("c0","c0",1.1*300*3,300*2);
	c0->Divide(3,2);
	
	for (int im=0; im<nmult; im++){
		c0->cd(im+1);
		gPad->SetLogz();
		SetPadStyle(1);

		htmp = (TH1D*)gPad->DrawFrame(0,0,30,50);
		SetHistoStyle("p_{T, lead} (GeV/c)","Q_{hard} (GeV)","",18,18);
		htmp->GetYaxis()->SetTitleOffset(2.0);
		htmp->GetXaxis()->SetTitleOffset(2.0);

		hevent_Q[im]->Draw("colz same");

		TProfile *hprofx = (TProfile*)hevent_Q[im]->ProfileX(Form("hprofx_mult%02d",im));
		hprofx->SetMarkerStyle(20);
		hprofx->SetLineWidth(2);
		hprofx->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.7,0.7,0.9);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");
		leg->AddEntry("","SoftQCD:nonDiffractive","h");
		leg->AddEntry("",Form("V0M %d%c-%d%c",cut_mult[im+1],'%',cut_mult[im],'%'),"h");
		leg->Draw();
	}


	TCanvas *c1 = new TCanvas("c1","c1",1.1*300*3,300*2);
	c1->Divide(3,2);
	
	for (int im=0; im<nmult; im++){
		c1->cd(im+1);
		gPad->SetLogz();
		SetPadStyle(1);

		htmp = (TH1D*)gPad->DrawFrame(0,0,30,30);
		SetHistoStyle("p_{T, lead} (GeV/c)","# of MPI","",18,18);
		htmp->GetYaxis()->SetTitleOffset(2.0);
		htmp->GetXaxis()->SetTitleOffset(2.0);

		hevent_nMPI[im]->Draw("colz same");

		TProfile *hprofx = (TProfile*)hevent_nMPI[im]->ProfileX(Form("hprofx2_mult%02d",im));
		hprofx->SetMarkerStyle(20);
		hprofx->SetLineWidth(2);
		hprofx->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.7,0.7,0.9);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");
		leg->AddEntry("","SoftQCD:nonDiffractive","h");
		leg->AddEntry("",Form("V0M %d%c-%d%c",cut_mult[im+1],'%',cut_mult[im],'%'),"h");
		leg->Draw();
	}

	TCanvas *c2 = new TCanvas("c2","c2",1.1*300*3,300*2);
	c2->Divide(3,2);
	
	for (int im=0; im<nmult; im++){
		c2->cd(im+1);
		gPad->SetLogz();
		SetPadStyle(1);

		htmp = (TH1D*)gPad->DrawFrame(0,0,30,2.5);
		SetHistoStyle("p_{T, lead} (GeV/c)","Impact paramter for MPI","",18,18);
		htmp->GetYaxis()->SetTitleOffset(2.0);
		htmp->GetXaxis()->SetTitleOffset(2.0);

		hevent_bMPI[im]->Draw("colz same");

		TProfile *hprofx = (TProfile*)hevent_bMPI[im]->ProfileX(Form("hprofx3_mult%02d",im));
		hprofx->SetMarkerStyle(20);
		hprofx->SetLineWidth(2);
		hprofx->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.7,0.7,0.9);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");
		leg->AddEntry("","SoftQCD:nonDiffractive","h");
		leg->AddEntry("",Form("V0M %d%c-%d%c",cut_mult[im+1],'%',cut_mult[im],'%'),"h");
		leg->Draw();
	}

	TCanvas *c2_ = new TCanvas("c2_","c2_",1.1*300*3,300*2);
	c2_->Divide(3,2);
	
	for (int im=0; im<nmult; im++){
		c2_->cd(im+1);
		gPad->SetLogz();
		SetPadStyle(1);

		htmp = (TH1D*)gPad->DrawFrame(0,0,30,2.5);
		SetHistoStyle("Q_{hard} (GeV/c)","Impact paramter for MPI","",18,18);
		htmp->GetYaxis()->SetTitleOffset(2.0);
		htmp->GetXaxis()->SetTitleOffset(2.0);

		hevent_Q_bMPI[im]->Draw("colz same");

		TProfile *hprofx = (TProfile*)hevent_Q_bMPI[im]->ProfileX(Form("hprofx4_mult%02d",im));
		hprofx->SetMarkerStyle(20);
		hprofx->SetLineWidth(2);
		hprofx->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.7,0.7,0.9);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");
		leg->AddEntry("","SoftQCD:nonDiffractive","h");
		leg->AddEntry("",Form("V0M %d%c-%d%c",cut_mult[im+1],'%',cut_mult[im],'%'),"h");
		leg->Draw();
	}

	if ( bSAVE ){
		c0->cd();
		c0->SaveAs("plots/Pythia8_pp13TeV_set00_grp000_event_info_00.pdf");

		c1->cd();
		c1->SaveAs("plots/Pythia8_pp13TeV_set00_grp000_event_info_01.pdf");

		c2->cd();
		c2->SaveAs("plots/Pythia8_pp13TeV_set00_grp000_event_info_02.pdf");

		c2_->cd();
		c2_->SaveAs("plots/Pythia8_pp13TeV_set00_grp000_event_info_03.pdf");
	}

}

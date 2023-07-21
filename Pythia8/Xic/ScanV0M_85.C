#include "Style.h"

using namespace RooFit;

TH1D *h3body;
TH1D *h4body;
TH1D *hbkg;

double fitfunc(Double_t *x, Double_t *par) {
	double xx = x[0];
	int bin = h3body->GetXaxis()->FindBin(xx);

	double y3body = h3body->GetBinContent(bin);
	double y4body = h4body->GetBinContent(bin);
	double ybkg = hbkg->GetBinContent(bin);

	return par[0]*(par[1]*y3body + par[2]*y4body + (1-par[1]-par[2])*ybkg);
}

void ScanV0M_85(int opt=0){

	const float legtextsize = 0.04;

	const int npt = 6;
	const float ptbin[npt+1] = {1, 2, 4, 6, 8, 12, 24};

	//0-1% setting
	//const float bSFBKG = 15.0;
	//const float bSF4BODY = 10.0;

	const float bSFBKG = 7.5;
	const float bSF4BODY = 10.0;

	bool bSAVE = false;

	TFile *infile;
	string multname;

	if ( opt==0 ){
		infile = new TFile("outfile_hist_pp13TeV_set85_grp000_try100.root","read"); //MB
		multname = "0-100% V0M";
	}else if ( opt==1 ){
		infile = new TFile("outfile_hist_pp13TeV_set85_grp000_try101.root","read"); //0-5%
		multname = "0-5% V0M";
	}else if ( opt==2 ){
		infile = new TFile("outfile_hist_pp13TeV_set85_grp000_try102.root","read"); //0-1%
		multname = "0-1% V0M";
	}else if ( opt==3 ){
		infile = new TFile("outfile_hist_pp13TeV_set85_grp000_try103.root","read"); //30-100%
		multname = "30-100% V0M";
	}

	TH2D *hmass_pt_sig3Xic0 = (TH2D*)infile->Get("hmass_pt_sig3Xic0");
	TH2D *hmass_pt_sig4Xic0 = (TH2D*)infile->Get("hmass_pt_sig4Xic0");
	TH2D *hmass_pt_sig4Xicp = (TH2D*)infile->Get("hmass_pt_sig4Xicp");

	TH2D *hmass_pt_cut_sig3Xic0;
	TH2D *hmass_pt_cut_sig4Xic0;
	TH2D *hmass_pt_cut_sig4Xicp;

	hmass_pt_cut_sig3Xic0 = (TH2D*)infile->Get("hmass_pt_cut_sig3Xic0");
	hmass_pt_cut_sig4Xic0 = (TH2D*)infile->Get("hmass_pt_cut_sig4Xic0");
	hmass_pt_cut_sig4Xicp = (TH2D*)infile->Get("hmass_pt_cut_sig4Xicp");

	TH1D *hmass_sig3Xic0 = (TH1D*)hmass_pt_sig3Xic0->ProjectionX("hmass_sig3Xic0");
	TH1D *hmass_sig4Xic0 = (TH1D*)hmass_pt_sig4Xic0->ProjectionX("hmass_sig4Xic0");
	TH1D *hmass_sig4Xicp = (TH1D*)hmass_pt_sig4Xicp->ProjectionX("hmass_sig4Xicp");

	TH1D *hmass_cut_sig3Xic0 = (TH1D*)hmass_pt_cut_sig3Xic0->ProjectionX("hmass_cut_sig3Xic0");
	TH1D *hmass_cut_sig4Xic0 = (TH1D*)hmass_pt_cut_sig4Xic0->ProjectionX("hmass_cut_sig4Xic0");
	TH1D *hmass_cut_sig4Xicp = (TH1D*)hmass_pt_cut_sig4Xicp->ProjectionX("hmass_cut_sig4Xicp");

	TH1D *hpt_cut_sig3Xic0 = (TH1D*)hmass_pt_cut_sig3Xic0->ProjectionY("hpt_cut_sig3Xic0");
	TH1D *hpt_cut_sig4Xic0 = (TH1D*)hmass_pt_cut_sig4Xic0->ProjectionY("hpt_cut_sig4Xic0");
	TH1D *hpt_cut_sig4Xicp = (TH1D*)hmass_pt_cut_sig4Xicp->ProjectionY("hpt_cut_sig4Xicp");

	hmass_sig3Xic0->SetLineColor(1);
	hmass_sig4Xic0->SetLineColor(2);
	hmass_sig4Xicp->SetLineColor(4);

	hmass_sig3Xic0->SetLineWidth(2);
	hmass_sig4Xic0->SetLineWidth(2);
	hmass_sig4Xicp->SetLineWidth(2);

	hmass_cut_sig3Xic0->SetLineColor(1);
	hmass_cut_sig4Xic0->SetLineColor(2);
	hmass_cut_sig4Xicp->SetLineColor(4);

	hmass_cut_sig3Xic0->SetLineWidth(2);
	hmass_cut_sig4Xic0->SetLineWidth(2);
	hmass_cut_sig4Xicp->SetLineWidth(2);

	TH2D *hmass_pt_comb[3];
	TH2D *hmass_pt_mixed1[3];
	TH2D *hmass_pt_mixed2[3];
	TH2D *hmass_pt_corrJet[3];
	TH2D *hmass_pt_corrXibm[3];
	TH2D *hmass_pt_corrDiquark[3];

	TH2D *hmass_pt_cut_comb[3];
	TH2D *hmass_pt_cut_mixed1[3];
	TH2D *hmass_pt_cut_mixed2[3];
	TH2D *hmass_pt_cut_corrJet[3];
	TH2D *hmass_pt_cut_corrXibm[3];
	TH2D *hmass_pt_cut_corrDiquark[3];

	TH2D *hoa_pt_comb[3];
	TH2D *hoa_pt_mixed1[3];
	TH2D *hoa_pt_mixed2[3];
	TH2D *hoa_pt_sig3Xic0[3];
	TH2D *hoa_pt_sig4Xic0[3];
	TH2D *hoa_pt_sig4Xicp[3];
	TH2D *hoa_pt_corrJet[3];
	TH2D *hoa_pt_corrXibm[3];
	TH2D *hoa_pt_corrDiquark[3];

	TH1D *hmass_comb[3];
	TH1D *hmass_mixed1[3];
	TH1D *hmass_mixed2[3];
	TH1D *hmass_corrJet[3];
	TH1D *hmass_corrXibm[3];
	TH1D *hmass_corrDiquark[3];

	TH1D *hmass_cut_comb[3];
	TH1D *hmass_cut_mixed1[3];
	TH1D *hmass_cut_mixed2[3];
	TH1D *hmass_cut_corrJet[3];
	TH1D *hmass_cut_corrXibm[3];
	TH1D *hmass_cut_corrDiquark[3];

	TH2D *hmass_pt_rs_cut = (TH2D*)hmass_pt_cut_sig3Xic0->Clone("hmass_pt_rs_cut");
	hmass_pt_rs_cut->Add(hmass_pt_cut_sig4Xic0, bSF4BODY);
	hmass_pt_rs_cut->Add(hmass_pt_cut_sig4Xicp, bSF4BODY);

	TH2D *hmass_pt_rs_cut_sig = (TH2D*)hmass_pt_cut_sig3Xic0->Clone("hmass_pt_rs_cut_sig");
	hmass_pt_rs_cut_sig->Add(hmass_pt_cut_sig4Xic0, bSF4BODY);
	hmass_pt_rs_cut_sig->Add(hmass_pt_cut_sig4Xicp, bSF4BODY);

	TH2D *hmass_pt_ws_cut = (TH2D*)hmass_pt_cut_sig3Xic0->Clone("hmass_pt_ws_cut");
	hmass_pt_ws_cut->Reset();
	cout << hmass_pt_ws_cut->Integral() << endl;

	TH2D *hmass_pt_rs_cut_bkg = (TH2D*)hmass_pt_cut_sig3Xic0->Clone("hmass_pt_rs_cut_bkg");
	hmass_pt_rs_cut_bkg->Reset();
	cout << hmass_pt_rs_cut_bkg->Integral() << endl;

	TH2D *hmass_pt_cs_cut_mixed = (TH2D*)hmass_pt_cut_sig3Xic0->Clone("hmass_pt_cs_cut_mixed");
	hmass_pt_cs_cut_mixed->Reset();
	cout << hmass_pt_cs_cut_mixed->Integral() << endl;

	for (int ii=0; ii<3; ii++){

		hmass_pt_comb[ii] = (TH2D*)infile->Get(Form("hmass_pt_comb_chg%d",ii));
		hmass_pt_mixed1[ii] = (TH2D*)infile->Get(Form("hmass_pt_mixed1_chg%d",ii));
		hmass_pt_mixed2[ii] = (TH2D*)infile->Get(Form("hmass_pt_mixed2_chg%d",ii));
		hmass_pt_corrJet[ii] = (TH2D*)infile->Get(Form("hmass_pt_corrJet_chg%d",ii));
		hmass_pt_corrXibm[ii] = (TH2D*)infile->Get(Form("hmass_pt_corrXibm_chg%d",ii));
		hmass_pt_corrDiquark[ii] = (TH2D*)infile->Get(Form("hmass_pt_corrDiquark_chg%d",ii));

		hmass_pt_cut_comb[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut_comb_chg%d",ii));
		hmass_pt_cut_mixed1[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut_mixed1_chg%d",ii));
		hmass_pt_cut_mixed2[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut_mixed2_chg%d",ii));
		hmass_pt_cut_corrJet[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut_corrJet_chg%d",ii));
		hmass_pt_cut_corrXibm[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut_corrXibm_chg%d",ii));
		hmass_pt_cut_corrDiquark[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut_corrDiquark_chg%d",ii));

		hoa_pt_comb[ii] = (TH2D*)infile->Get(Form("hoa_pt_comb_chg%d",ii));
		hoa_pt_mixed1[ii] = (TH2D*)infile->Get(Form("hoa_pt_mixed1_chg%d",ii));
		hoa_pt_mixed2[ii] = (TH2D*)infile->Get(Form("hoa_pt_mixed2_chg%d",ii));
		hoa_pt_sig3Xic0[ii] = (TH2D*)infile->Get(Form("hoa_pt_sig3Xic0_chg%d",ii));
		hoa_pt_sig4Xic0[ii] = (TH2D*)infile->Get(Form("hoa_pt_sig4Xic0_chg%d",ii));
		hoa_pt_sig4Xicp[ii] = (TH2D*)infile->Get(Form("hoa_pt_sig4Xicp_chg%d",ii));
		hoa_pt_corrJet[ii] = (TH2D*)infile->Get(Form("hoa_pt_corrJet_chg%d",ii));
		hoa_pt_corrXibm[ii] = (TH2D*)infile->Get(Form("hoa_pt_corrXibm_chg%d",ii));
		hoa_pt_corrDiquark[ii] = (TH2D*)infile->Get(Form("hoa_pt_corrDiquark_chg%d",ii));

		if ( ii==0 ){
			hmass_pt_rs_cut->Add(hmass_pt_cut_comb[ii], bSFBKG);
			hmass_pt_rs_cut->Add(hmass_pt_cut_corrJet[ii], bSFBKG);
			hmass_pt_rs_cut->Add(hmass_pt_cut_corrDiquark[ii], bSFBKG);
			hmass_pt_rs_cut->Add(hmass_pt_cut_corrXibm[ii]);

			hmass_pt_rs_cut_bkg->Add(hmass_pt_cut_comb[ii], bSFBKG);
			hmass_pt_rs_cut_bkg->Add(hmass_pt_cut_corrJet[ii], bSFBKG);
			hmass_pt_rs_cut_bkg->Add(hmass_pt_cut_corrDiquark[ii], bSFBKG);
			hmass_pt_rs_cut_bkg->Add(hmass_pt_cut_corrXibm[ii]);

			hmass_pt_cs_cut_mixed->Add(hmass_pt_cut_mixed1[ii]);
			hmass_pt_cs_cut_mixed->Add(hmass_pt_cut_mixed2[ii]);
		}else{
			hmass_pt_ws_cut->Add(hmass_pt_cut_comb[ii], bSFBKG);
			hmass_pt_ws_cut->Add(hmass_pt_cut_corrJet[ii], bSFBKG);
			hmass_pt_ws_cut->Add(hmass_pt_cut_corrDiquark[ii], bSFBKG);
			hmass_pt_ws_cut->Add(hmass_pt_cut_corrXibm[ii]);

			hmass_pt_cs_cut_mixed->Add(hmass_pt_cut_mixed1[ii]);
			hmass_pt_cs_cut_mixed->Add(hmass_pt_cut_mixed2[ii]);
		}

		hmass_comb[ii] = (TH1D*)hmass_pt_comb[ii]->ProjectionX(Form("hmass_comb_chg%d",ii));
		hmass_comb[ii]->SetLineColor(ii+1);
		hmass_comb[ii]->SetLineWidth(2);

		hmass_mixed1[ii] = (TH1D*)hmass_pt_mixed1[ii]->ProjectionX(Form("hmass_mixed1_chg%d",ii));
		hmass_mixed1[ii]->SetLineColor(ii+1);
		hmass_mixed1[ii]->SetMarkerStyle(24);
		hmass_mixed1[ii]->SetMarkerColor(ii+1);
		hmass_mixed1[ii]->SetMarkerSize(0.8);
		hmass_mixed1[ii]->SetLineWidth(2);

		hmass_mixed2[ii] = (TH1D*)hmass_pt_mixed2[ii]->ProjectionX(Form("hmass_mixed2_chg%d",ii));
		hmass_mixed2[ii]->SetLineColor(ii+1);
		hmass_mixed2[ii]->SetMarkerStyle(25);
		hmass_mixed2[ii]->SetMarkerColor(ii+1);
		hmass_mixed2[ii]->SetMarkerSize(0.8);
		hmass_mixed2[ii]->SetLineWidth(2);

		hmass_corrJet[ii] = (TH1D*)hmass_pt_corrJet[ii]->ProjectionX(Form("hmass_corrJet_chg%d",ii));
		hmass_corrJet[ii]->SetLineColor(ii+1);
		hmass_corrJet[ii]->SetLineWidth(2);

		hmass_corrXibm[ii] = (TH1D*)hmass_pt_corrXibm[ii]->ProjectionX(Form("hmass_corrXibm_chg%d",ii));
		hmass_corrXibm[ii]->SetLineColor(ii+1);
		hmass_corrXibm[ii]->SetLineWidth(2);

		hmass_corrDiquark[ii] = (TH1D*)hmass_pt_corrDiquark[ii]->ProjectionX(Form("hmass_corrDiquark_chg%d",ii));
		hmass_corrDiquark[ii]->SetLineColor(ii+1);
		hmass_corrDiquark[ii]->SetLineWidth(2);

		hmass_cut_comb[ii] = (TH1D*)hmass_pt_cut_comb[ii]->ProjectionX(Form("hmass_cut_comb_chg%d",ii));
		hmass_cut_comb[ii]->SetLineColor(ii+1);
		hmass_cut_comb[ii]->SetLineWidth(2);

		hmass_cut_mixed1[ii] = (TH1D*)hmass_pt_cut_mixed1[ii]->ProjectionX(Form("hmass_cut_mixed1_chg%d",ii));
		hmass_cut_mixed1[ii]->SetLineColor(ii+1);
		hmass_cut_mixed1[ii]->SetMarkerStyle(24);
		hmass_cut_mixed1[ii]->SetMarkerColor(ii+1);
		hmass_cut_mixed1[ii]->SetMarkerSize(0.8);
		hmass_cut_mixed1[ii]->SetLineWidth(2);

		hmass_cut_mixed2[ii] = (TH1D*)hmass_pt_cut_mixed2[ii]->ProjectionX(Form("hmass_cut_mixed2_chg%d",ii));
		hmass_cut_mixed2[ii]->SetLineColor(ii+1);
		hmass_cut_mixed2[ii]->SetMarkerStyle(25);
		hmass_cut_mixed2[ii]->SetMarkerColor(ii+1);
		hmass_cut_mixed2[ii]->SetMarkerSize(0.8);
		hmass_cut_mixed2[ii]->SetLineWidth(2);

		hmass_cut_corrJet[ii] = (TH1D*)hmass_pt_cut_corrJet[ii]->ProjectionX(Form("hmass_cut_corrJet_chg%d",ii));
		hmass_cut_corrJet[ii]->SetLineColor(ii+1);
		hmass_cut_corrJet[ii]->SetLineWidth(2);

		hmass_cut_corrXibm[ii] = (TH1D*)hmass_pt_cut_corrXibm[ii]->ProjectionX(Form("hmass_cut_corrXibm_chg%d",ii));
		hmass_cut_corrXibm[ii]->SetLineColor(ii+1);
		hmass_cut_corrXibm[ii]->SetLineWidth(2);

		hmass_cut_corrDiquark[ii] = (TH1D*)hmass_pt_cut_corrDiquark[ii]->ProjectionX(Form("hmass_cut_corrDiquark_chg%d",ii));
		hmass_cut_corrDiquark[ii]->SetLineColor(ii+1);
		hmass_cut_corrDiquark[ii]->SetLineWidth(2);
	}

	if ( 0 )
	{
		TCanvas *c0 = new TCanvas("c0","c0",1.2*2*400,400);
		c0->Divide(2,1);

		c0->cd(1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(0,1,15,2*hpt_cut_sig3Xic0->GetMaximum());
		SetHistoStyle("e#Xi p_{T} (GeV/c)","","",20,18);

		hpt_cut_sig4Xicp->Scale(3);

		hpt_cut_sig3Xic0->Sumw2();
		hpt_cut_sig3Xic0->SetLineColor(1);
		hpt_cut_sig3Xic0->SetMarkerSize(0.8);
		hpt_cut_sig3Xic0->SetMarkerStyle(20);
		hpt_cut_sig3Xic0->SetMarkerColor(1);
		hpt_cut_sig3Xic0->Draw("p same");

		hpt_cut_sig4Xic0->Sumw2();
		hpt_cut_sig4Xic0->SetLineColor(2);
		hpt_cut_sig4Xic0->SetMarkerSize(0.8);
		hpt_cut_sig4Xic0->SetMarkerStyle(21);
		hpt_cut_sig4Xic0->SetMarkerColor(2);
		hpt_cut_sig4Xic0->Draw("p same");

		hpt_cut_sig4Xicp->Sumw2();
		hpt_cut_sig4Xicp->SetLineColor(4);
		hpt_cut_sig4Xicp->SetMarkerSize(0.8);
		hpt_cut_sig4Xicp->SetMarkerStyle(21);
		hpt_cut_sig4Xicp->SetMarkerColor(4);
		hpt_cut_sig4Xicp->Draw("p same");

		TH1D *hsum = (TH1D*)hpt_cut_sig3Xic0->Clone("hsum");
		hsum->Add(hpt_cut_sig4Xic0);
		hsum->Add(hpt_cut_sig4Xicp);

		TH1D *hpt_R_sig3Xic0 = (TH1D*)hpt_cut_sig3Xic0->Clone("hpt_R_sig3Xic0");
		TH1D *hpt_R_sig4Xic0 = (TH1D*)hpt_cut_sig4Xic0->Clone("hpt_R_sig4Xic0");
		TH1D *hpt_R_sig4Xicp = (TH1D*)hpt_cut_sig4Xicp->Clone("hpt_R_sig4Xicp");

		hpt_R_sig3Xic0->Divide(hsum);
		hpt_R_sig4Xic0->Divide(hsum);
		hpt_R_sig4Xicp->Divide(hsum);

		TLegend *leg = new TLegend(0.2,0.2,0.6,0.4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("",multname.c_str(),"");
		leg->AddEntry(hpt_cut_sig3Xic0,Form("#Xi_{c}^{0} 3body"),"PL");
		leg->AddEntry(hpt_cut_sig4Xic0,Form("#Xi_{c}^{0} 4body"),"PL");
		leg->AddEntry(hpt_cut_sig4Xicp,Form("#Xi_{c}^{+} 4body"),"PL");
		leg->Draw();

		c0->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,15,1);
		SetHistoStyle("e#Xi p_{T} (GeV/c)","Fraction","",20,18);
		hpt_R_sig3Xic0->Draw("p same");
		hpt_R_sig4Xic0->Draw("p same");
		hpt_R_sig4Xicp->Draw("p same");

		TF1 *f1 = new TF1("f1","[0]*([1] + x)",0,15);
		hpt_R_sig3Xic0->Fit(f1,"R0");
		f1->Draw("same");

	}

	//return;


	TCanvas *c1 = new TCanvas("c1","c1",1.2*2*400,2*400);
	c1->Divide(2,2);
	{
		c1->cd(1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_sig3Xic0->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(1.0);
		hmass_sig3Xic0->Draw("same");
		hmass_sig4Xic0->Draw("same");
		hmass_sig4Xicp->Draw("same");

		TLegend *leg = new TLegend(0.55,0.93-0.045*3,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry(hmass_sig3Xic0,Form("#Xi_{c}^{0} 3body, %d",int(hmass_sig3Xic0->Integral())),"L");
		leg->AddEntry(hmass_sig4Xic0,Form("#Xi_{c}^{0} 4body, %d",int(hmass_sig4Xic0->Integral())),"L");
		leg->AddEntry(hmass_sig4Xicp,Form("#Xi_{c}^{+} 4body, %d",int(hmass_sig4Xicp->Integral())),"L");
		leg->Draw();
	}

	{
		c1->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_comb[0]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(1.0);
		hmass_comb[0]->Draw("same");
		hmass_comb[1]->Draw("same");
		hmass_comb[2]->Draw("same");

		hmass_mixed1[0]->Scale(hmass_comb[0]->Integral()/hmass_mixed1[0]->Integral());
		hmass_mixed1[1]->Scale(hmass_comb[0]->Integral()/hmass_mixed1[1]->Integral());

		hmass_mixed1[0]->Draw("p same");
		hmass_mixed1[1]->Draw("p same");

		hmass_mixed2[0]->Scale(hmass_comb[0]->Integral()/hmass_mixed2[0]->Integral());
		hmass_mixed2[1]->Scale(hmass_comb[0]->Integral()/hmass_mixed2[1]->Integral());

		hmass_mixed2[0]->Draw("p same");
		hmass_mixed2[1]->Draw("p same");

		TLegend *leg = new TLegend(0.45,0.93-0.045*8,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","Combinatorial e#Xi pair","");
		leg->AddEntry(hmass_comb[0],Form("Unlike-sign, %d",int(hmass_comb[0]->Integral())),"L");
		leg->AddEntry(hmass_comb[1],Form("Like-sign --, %d",int(hmass_comb[1]->Integral())),"L");
		leg->AddEntry(hmass_comb[2],Form("Like-sign ++, %d",int(hmass_comb[2]->Integral())),"L");
		leg->AddEntry(hmass_mixed1[0],"Mixed-event 1, Unlike-sign","P");
		leg->AddEntry(hmass_mixed1[1],"Mixed-event 1, Like-sign","P");
		leg->AddEntry(hmass_mixed2[0],"Mixed-event 2, Unlike-sign","P");
		leg->AddEntry(hmass_mixed2[1],"Mixed-event 2, Like-sign","P");
		leg->Draw();
	}

	{
		c1->cd(3);
		SetPadStyle();

		hmass_corrJet[0]->Add(hmass_corrDiquark[0]);
		hmass_corrJet[1]->Add(hmass_corrDiquark[1]);
		hmass_corrJet[2]->Add(hmass_corrDiquark[2]);

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_corrJet[0]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);

		hmass_corrJet[0]->Draw("same");
		hmass_corrJet[1]->Draw("same");
		hmass_corrJet[2]->Draw("same");

		TLegend *leg = new TLegend(0.45,0.93-0.045*4,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","Correlated e#Xi pair from jet","");
		leg->AddEntry(hmass_corrJet[0],Form("Unlike-sign, %d",int(hmass_corrJet[0]->Integral())),"L");
		leg->AddEntry(hmass_corrJet[1],Form("Like-sign --, %d",int(hmass_corrJet[1]->Integral())),"L");
		leg->AddEntry(hmass_corrJet[2],Form("Like-sign ++, %d",int(hmass_corrJet[2]->Integral())),"L");
		leg->Draw();
	}

	{
		c1->cd(4);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.3*hmass_corrXibm[1]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		hmass_corrXibm[0]->Draw("same");
		hmass_corrXibm[1]->Draw("same");

		TLegend *leg = new TLegend(0.45,0.93-0.05*3,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","Correlated e#Xi pair from #Xi_{b}^{-}","");
		leg->AddEntry(hmass_corrXibm[0],Form("Unlike-sign, %d",int(hmass_corrXibm[0]->Integral())),"L");
		leg->AddEntry(hmass_corrXibm[1],Form("Like-sign, %d",int(hmass_corrXibm[1]->Integral())),"L");
		leg->Draw();
	}

	//return;

	TCanvas *c2 = new TCanvas("c2","c2",1.2*2*400,2*400);
	c2->Divide(2,2);
	{
		c2->cd(1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_cut_sig3Xic0->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		hmass_cut_sig3Xic0->Draw("same");
		hmass_cut_sig4Xic0->Draw("same");
		hmass_cut_sig4Xicp->Draw("same");

		TLegend *leg = new TLegend(0.55,0.93-0.05*3,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry(hmass_cut_sig3Xic0,Form("#Xi_{c}^{0} 3body, %d",int(hmass_cut_sig3Xic0->Integral())),"L");
		leg->AddEntry(hmass_cut_sig4Xic0,Form("#Xi_{c}^{0} 4body, %d",int(hmass_cut_sig4Xic0->Integral())),"L");
		leg->AddEntry(hmass_cut_sig4Xicp,Form("#Xi_{c}^{+} 4body, %d",int(hmass_cut_sig4Xicp->Integral())),"L");
		leg->Draw();
	}

	{
		c2->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_cut_comb[0]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		hmass_cut_comb[0]->Draw("same");
		hmass_cut_comb[1]->Draw("same");
		hmass_cut_comb[2]->Draw("same");

		int massmin = hmass_cut_comb[0]->FindBin(1.0+0.001);
		int massmax = hmass_cut_comb[0]->FindBin(5.0-0.001);

		hmass_cut_mixed1[0]->Scale(hmass_cut_comb[0]->Integral(massmin,massmax)/hmass_cut_mixed1[0]->Integral(massmin,massmax));
		hmass_cut_mixed1[1]->Scale(hmass_cut_comb[0]->Integral(massmin,massmax)/hmass_cut_mixed1[1]->Integral(massmin,massmax));
		hmass_cut_mixed1[0]->Draw("p same");
		hmass_cut_mixed1[1]->Draw("p same");

		hmass_cut_mixed2[0]->Scale(hmass_cut_comb[0]->Integral(massmin,massmax)/hmass_cut_mixed2[0]->Integral(massmin,massmax));
		hmass_cut_mixed2[1]->Scale(hmass_cut_comb[0]->Integral(massmin,massmax)/hmass_cut_mixed2[1]->Integral(massmin,massmax));
		hmass_cut_mixed2[0]->Draw("p same");
		hmass_cut_mixed2[1]->Draw("p same");

		TLegend *leg = new TLegend(0.45,0.93-0.05*8,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","Combinatorial e#Xi pair","");
		leg->AddEntry(hmass_cut_comb[0],Form("Unlike-sign, %d",int(hmass_cut_comb[0]->Integral())),"L");
		leg->AddEntry(hmass_cut_comb[1],Form("Like-sign --, %d",int(hmass_cut_comb[1]->Integral())),"L");
		leg->AddEntry(hmass_cut_comb[2],Form("Like-sign ++, %d",int(hmass_cut_comb[2]->Integral())),"L");
		leg->AddEntry(hmass_cut_mixed1[0],"Mixed-event 1, Unlike-sign","P");
		leg->AddEntry(hmass_cut_mixed1[1],"Mixed-event 1, Like-sign","P");
		leg->AddEntry(hmass_cut_mixed2[0],"Mixed-event 2, Unlike-sign","P");
		leg->AddEntry(hmass_cut_mixed2[1],"Mixed-event 2, Like-sign","P");
		leg->Draw();
	}

	{
		c2->cd(3);
		SetPadStyle();
	
		hmass_cut_corrJet[0]->Add(hmass_cut_corrDiquark[0]);
		hmass_cut_corrJet[1]->Add(hmass_cut_corrDiquark[1]);

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_cut_corrJet[0]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);

		TH1D *hmass_cut_mixed1_cj[2];
		for(int ii=0; ii<2; ii++){
			hmass_cut_mixed1_cj[ii] = (TH1D*)hmass_cut_mixed1[ii]->Clone(Form("hmass_cut_mixed1_cj_%d",ii));

			int massmin = hmass_cut_mixed1_cj[ii]->FindBin(3.0+0.001);
			int massmax = hmass_cut_mixed1_cj[ii]->FindBin(5.0-0.001);

			hmass_cut_mixed1_cj[ii]->Scale(hmass_cut_corrJet[0]->Integral(massmin,massmax)/hmass_cut_mixed1_cj[ii]->Integral(massmin,massmax));
		}

		hmass_cut_corrJet[0]->Draw("same");
		hmass_cut_corrJet[1]->Draw("same");

		hmass_cut_mixed1_cj[0]->Draw("p same");

		TLegend *leg = new TLegend(0.45,0.93-0.05*4,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","Correlated e#Xi pair from jet","");
		leg->AddEntry(hmass_cut_corrJet[0],Form("Unlike-sign, %d",int(hmass_cut_corrJet[0]->Integral())),"L");
		leg->AddEntry(hmass_cut_corrJet[1],Form("Like-sign, %d",int(hmass_cut_corrJet[1]->Integral())),"L");
		leg->AddEntry(hmass_cut_mixed1_cj[0],Form("Mixed-event 1, Unlike-sign"),"P");
		leg->Draw();
	}

	{
		c2->cd(4);
		SetPadStyle();
	
		TH1D *hmass_cut_bkg_all[2];
		TH1D *hmass_cut_mixed1_all[2];
		for(int ii=0; ii<2; ii++){

			hmass_cut_bkg_all[ii] = (TH1D*)hmass_cut_corrJet[ii]->Clone(Form("hmass_cut_bkg_all_%d",ii));
			hmass_cut_bkg_all[ii]->Add(hmass_cut_comb[ii]);
			hmass_cut_mixed1_all[ii] = (TH1D*)hmass_cut_mixed1[ii]->Clone(Form("hmass_cut_mixed1_all_%d",ii));

			int massmin = hmass_cut_mixed1_all[ii]->FindBin(3.0+0.001);
			int massmax = hmass_cut_mixed1_all[ii]->FindBin(5.0-0.001);

			hmass_cut_mixed1_all[ii]->Scale(hmass_cut_bkg_all[0]->Integral(massmin,massmax)/hmass_cut_mixed1_all[ii]->Integral(massmin,massmax));
		}

		htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_cut_bkg_all[0]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);

		hmass_cut_bkg_all[0]->Draw("same");
		hmass_cut_bkg_all[1]->Draw("same");

		hmass_cut_mixed1_all[0]->Draw("p same");

		TLegend *leg = new TLegend(0.4,0.93-0.05*4,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","Combinatorial + Correlated bkg.","");
		leg->AddEntry(hmass_cut_corrJet[0],Form("Unlike-sign, %d",int(hmass_cut_bkg_all[0]->Integral())),"L");
		leg->AddEntry(hmass_cut_corrJet[1],Form("Like-sign, %d",int(hmass_cut_bkg_all[1]->Integral())),"L");
		leg->AddEntry(hmass_cut_mixed1_all[0],Form("Mixed-event 1, Unlike-sign"),"P");
		leg->Draw();
	}

	//return;

	TH1D *hmass_sig3Xic0_cut[npt];
	TH1D *hmass_sig4Xic0_cut[npt];
	TH1D *hmass_sig4Xicp_cut[npt];
	TH1D *hmass_rs_cut[npt];
	TH1D *hmass_rs_cut_sig[npt];
	TH1D *hmass_rs_cut_bkg[npt];
	TH1D *hmass_ws_cut[npt];
	TH1D *hmass_cs_cut_mixed[npt];

	TH1D *hmass_sub1_cut[npt];
	TH1D *hmass_sub2_cut[npt];

	TH1D *hoa_comb[2][npt];
	TH1D *hoa_mixed1[2][npt];
	TH1D *hoa_mixed2[2][npt];
	TH1D *hoa_sig3Xic0[2][npt];
	TH1D *hoa_sig4Xic0[2][npt];
	TH1D *hoa_sig4Xicp[2][npt];
	TH1D *hoa_corrJet[2][npt];
	TH1D *hoa_corrXibm[2][npt];
	TH1D *hoa_corrDiquark[2][npt];

	TH1D *hpt_sig_truth = new TH1D("hpt_sig_truth","",npt,ptbin);
	TH1D *hpt_3body_truth = new TH1D("hpt_3body_truth","",npt,ptbin);
	TH1D *hpt_3body_fit1 = new TH1D("hpt_3body_fit1","",npt,ptbin);
	TH1D *hpt_3body_fit2 = new TH1D("hpt_3body_fit2","",npt,ptbin);
	TH1D *hpt_3body_fit3 = new TH1D("hpt_3body_fit3","",npt,ptbin);
	TH1D *hpt_sub1 = new TH1D("hpt_sub1","",npt,ptbin);
	TH1D *hpt_sub2 = new TH1D("hpt_sub2","",npt,ptbin);

	TH1D *hpt_rs_cut = (TH1D*)hmass_pt_rs_cut->ProjectionY("hpt_rs_cut");

	TCanvas *c3;
	TCanvas *c4;
	TCanvas *c5;
	TCanvas *c7;
	TCanvas *c8;
	TCanvas *c9;

	c3 = new TCanvas("c3","c3",1.2*3*400,2*400);
	c3->Divide(3,2);

	c4 = new TCanvas("c4","c4",1.2*3*400,2*400);
	c4->Divide(3,2);

	//c5 = new TCanvas("c5","c5",1.2*3*400,2*400);
	//c5->Divide(3,2);

	c7 = new TCanvas("c7","c7",1.2*3*400,2*400);
	c7->Divide(3,2);

	c8 = new TCanvas("c8","c8",1.2*3*400,2*400);
	c8->Divide(3,2);

	c9 = new TCanvas("c9","c9",1.2*3*400,2*400);
	c9->Divide(3,2);

	for (int ii=0; ii<2; ii++){
		for (int ipt=0; ipt<npt; ipt++){

			int ptmin = hpt_rs_cut->FindBin(ptbin[ipt]);
			int ptmax = hpt_rs_cut->FindBin(ptbin[ipt+1]);

			hoa_comb[ii][ipt] = (TH1D*)hoa_pt_comb[ii]->ProjectionY(Form("hoa_comb_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_mixed1[ii][ipt] = (TH1D*)hoa_pt_mixed1[ii]->ProjectionY(Form("hoa_mixed1_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_mixed2[ii][ipt] = (TH1D*)hoa_pt_mixed2[ii]->ProjectionY(Form("hoa_mixed2_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_sig3Xic0[ii][ipt] = (TH1D*)hoa_pt_sig3Xic0[ii]->ProjectionY(Form("hoa_sig3Xic0_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_sig4Xic0[ii][ipt] = (TH1D*)hoa_pt_sig4Xic0[ii]->ProjectionY(Form("hoa_sig4Xic0_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_sig4Xicp[ii][ipt] = (TH1D*)hoa_pt_sig4Xicp[ii]->ProjectionY(Form("hoa_sig4Xicp_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_corrXibm[ii][ipt] = (TH1D*)hoa_pt_corrXibm[ii]->ProjectionY(Form("hoa_corrXibm_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_corrJet[ii][ipt] = (TH1D*)hoa_pt_corrJet[ii]->ProjectionY(Form("hoa_corrJet_%d_%d",ii,ipt),ptmin,ptmax);
			hoa_corrDiquark[ii][ipt] = (TH1D*)hoa_pt_corrDiquark[ii]->ProjectionY(Form("hoa_corrDiquark_%d_%d",ii,ipt),ptmin,ptmax);

			hoa_sig4Xicp[ii][ipt]->Add(hoa_sig4Xic0[ii][ipt]);
			hoa_corrJet[ii][ipt]->Add(hoa_corrDiquark[ii][ipt]);

			hoa_sig3Xic0[ii][ipt]->Rebin(3);
			hoa_sig4Xicp[ii][ipt]->Rebin(3);
			hoa_comb[ii][ipt]->Rebin(3);
			hoa_mixed1[ii][ipt]->Rebin(3);
			hoa_mixed2[ii][ipt]->Rebin(3);
			hoa_corrJet[ii][ipt]->Rebin(3);
			hoa_corrXibm[ii][ipt]->Rebin(3);

			hoa_comb[ii][ipt]->Scale(1./(hoa_comb[ii][ipt]->Integral()+0.01));
			hoa_mixed1[ii][ipt]->Scale(1./(hoa_mixed1[ii][ipt]->Integral()+0.01));
			hoa_mixed2[ii][ipt]->Scale(1./(hoa_mixed2[ii][ipt]->Integral()+0.01));
			hoa_sig3Xic0[ii][ipt]->Scale(1./(hoa_sig3Xic0[ii][ipt]->Integral()+0.01));
			hoa_sig4Xicp[ii][ipt]->Scale(1./(hoa_sig4Xicp[ii][ipt]->Integral()+0.01));
			hoa_corrJet[ii][ipt]->Scale(1./(hoa_corrJet[ii][ipt]->Integral()+0.01));

			hoa_sig3Xic0[ii][ipt]->SetMarkerStyle(20);
			hoa_sig3Xic0[ii][ipt]->SetMarkerColor(1);
			hoa_sig3Xic0[ii][ipt]->SetMarkerSize(0.8);
			hoa_sig3Xic0[ii][ipt]->SetLineColor(1);

			hoa_sig4Xicp[ii][ipt]->SetMarkerStyle(20);
			hoa_sig4Xicp[ii][ipt]->SetMarkerColor(2);
			hoa_sig4Xicp[ii][ipt]->SetMarkerSize(0.8);
			hoa_sig4Xicp[ii][ipt]->SetLineColor(2);

			hoa_comb[ii][ipt]->SetMarkerStyle(24+ii);
			hoa_comb[ii][ipt]->SetMarkerColor(4);
			hoa_comb[ii][ipt]->SetMarkerSize(0.8);
			hoa_comb[ii][ipt]->SetLineColor(4);

			hoa_mixed1[ii][ipt]->SetMarkerStyle(24+ii);
			hoa_mixed1[ii][ipt]->SetMarkerColor(6);
			hoa_mixed1[ii][ipt]->SetMarkerSize(0.8);
			hoa_mixed1[ii][ipt]->SetLineColor(6);

			hoa_mixed2[ii][ipt]->SetMarkerStyle(24+ii);
			hoa_mixed2[ii][ipt]->SetMarkerColor(kGreen+2);
			hoa_mixed2[ii][ipt]->SetMarkerSize(0.8);
			hoa_mixed2[ii][ipt]->SetLineColor(kGreen+2);

			hoa_corrJet[ii][ipt]->SetMarkerStyle(27+6*ii);
			hoa_corrJet[ii][ipt]->SetMarkerColor(1);
			hoa_corrJet[ii][ipt]->SetMarkerSize(1.2);
			hoa_corrJet[ii][ipt]->SetLineColor(1);

		}//ipt
	}//ii

	//return;

	for (int ipt=0; ipt<npt; ipt++){
	//for (int ipt=0; ipt<1; ipt++){

		int ptmin = hpt_rs_cut->FindBin(ptbin[ipt]);
		int ptmax = hpt_rs_cut->FindBin(ptbin[ipt+1]);

		hmass_rs_cut[ipt] = (TH1D*)hmass_pt_rs_cut->ProjectionX(Form("hmass_rs_cut_pt%d",ipt),ptmin,ptmax);
		hmass_rs_cut_bkg[ipt] = (TH1D*)hmass_pt_rs_cut_bkg->ProjectionX(Form("hmass_rs_cut_bkg_pt%d",ipt),ptmin,ptmax);
		hmass_rs_cut_sig[ipt] = (TH1D*)hmass_pt_rs_cut_sig->ProjectionX(Form("hmass_rs_cut_sig_pt%d",ipt),ptmin,ptmax);
		hmass_ws_cut[ipt] = (TH1D*)hmass_pt_ws_cut->ProjectionX(Form("hmass_ws_cut_pt%d",ipt),ptmin,ptmax);
		hmass_cs_cut_mixed[ipt] = (TH1D*)hmass_pt_cs_cut_mixed->ProjectionX(Form("hmass_cs_cut_mixed_pt%d",ipt),ptmin,ptmax);
		hmass_sig3Xic0_cut[ipt] = (TH1D*)hmass_pt_cut_sig3Xic0->ProjectionX(Form("hmass_sig3Xic0_cut_pt%d",ipt),ptmin,ptmax);
		hmass_sig4Xic0_cut[ipt] = (TH1D*)hmass_pt_cut_sig4Xic0->ProjectionX(Form("hmass_sig4Xic0_cut_pt%d",ipt),ptmin,ptmax);
		hmass_sig4Xicp_cut[ipt] = (TH1D*)hmass_pt_cut_sig4Xicp->ProjectionX(Form("hmass_sig4Xicp_cut_pt%d",ipt),ptmin,ptmax);

		hmass_sig4Xicp_cut[ipt]->Add(hmass_sig4Xic0_cut[ipt]);

		hmass_rs_cut[ipt]->SetMarkerStyle(20);
		hmass_rs_cut[ipt]->SetMarkerColor(1);
		hmass_rs_cut[ipt]->SetMarkerSize(0.8);
		hmass_rs_cut[ipt]->SetLineColor(1);
		hmass_rs_cut[ipt]->Sumw2();

		hmass_rs_cut_sig[ipt]->SetMarkerStyle(20);
		hmass_rs_cut_sig[ipt]->SetMarkerColor(1);
		hmass_rs_cut_sig[ipt]->SetMarkerSize(0.8);
		hmass_rs_cut_sig[ipt]->SetLineColor(1);
		hmass_rs_cut_sig[ipt]->Sumw2();

		hmass_rs_cut_bkg[ipt]->SetMarkerStyle(24);
		hmass_rs_cut_bkg[ipt]->SetMarkerColor(4);
		hmass_rs_cut_bkg[ipt]->SetMarkerSize(0.8);
		hmass_rs_cut_bkg[ipt]->SetLineColor(4);
		hmass_rs_cut_bkg[ipt]->Sumw2();

		hmass_ws_cut[ipt]->SetMarkerStyle(24);
		hmass_ws_cut[ipt]->SetMarkerColor(2);
		hmass_ws_cut[ipt]->SetMarkerSize(0.8);
		hmass_ws_cut[ipt]->SetLineColor(2);
		hmass_ws_cut[ipt]->Sumw2();

		hmass_cs_cut_mixed[ipt]->SetMarkerStyle(27);
		hmass_cs_cut_mixed[ipt]->SetMarkerColor(kGreen+3);
		hmass_cs_cut_mixed[ipt]->SetLineColor(kGreen+3);
		hmass_cs_cut_mixed[ipt]->Sumw2();

		hmass_sig3Xic0_cut[ipt]->SetLineWidth(2);
		hmass_sig3Xic0_cut[ipt]->SetLineColor(4);

		hmass_sig4Xicp_cut[ipt]->SetLineWidth(2);
		hmass_sig4Xicp_cut[ipt]->SetLineColor(6);

		int massmin = hmass_rs_cut_bkg[ipt]->FindBin(3.0+0.001);
		int massmax = hmass_rs_cut_bkg[ipt]->FindBin(5.0-0.001);

		hmass_cs_cut_mixed[ipt]->Scale(hmass_rs_cut_bkg[ipt]->Integral(massmin,massmax)/hmass_cs_cut_mixed[ipt]->Integral(massmin,massmax));

		hmass_sub1_cut[ipt] = (TH1D*)hmass_rs_cut[ipt]->Clone(Form("hmass_sub1_cut_pt%d",ipt));
		hmass_sub1_cut[ipt]->Add(hmass_ws_cut[ipt], -1);
		hmass_sub1_cut[ipt]->SetMarkerStyle(24);
		hmass_sub1_cut[ipt]->SetMarkerSize(0.8);
		hmass_sub1_cut[ipt]->SetMarkerColor(2);
		hmass_sub1_cut[ipt]->SetLineColor(2);

		hmass_sub2_cut[ipt] = (TH1D*)hmass_rs_cut[ipt]->Clone(Form("hmass_sub2_cut_pt%d",ipt));
		hmass_sub2_cut[ipt]->Add(hmass_cs_cut_mixed[ipt], -1);
		hmass_sub2_cut[ipt]->SetMarkerStyle(27);
		hmass_sub2_cut[ipt]->SetMarkerSize(1.2);
		hmass_sub2_cut[ipt]->SetMarkerColor(kGreen+3);
		hmass_sub2_cut[ipt]->SetLineColor(kGreen+3);

		if ( c3 ){
			c3->cd(ipt+1);
			SetPadStyle();

			htmp = (TH1D*)gPad->DrawFrame(1,0,5,1.2*hmass_rs_cut[ipt]->GetMaximum());
			SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
			htmp->GetXaxis()->SetTitleOffset(1.0);
			hmass_rs_cut[ipt]->Draw("p same");
			hmass_rs_cut_bkg[ipt]->Draw("p same");
			hmass_ws_cut[ipt]->Draw("p same");
			hmass_cs_cut_mixed[ipt]->Draw("p same");

			TLegend *leg = new TLegend(0.45,0.93-0.05*8,0.93,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(legtextsize);
			leg->AddEntry("","PYTHIA8 pp 13 TeV","");
			leg->AddEntry("",multname.c_str(),"");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
			leg->AddEntry(hmass_rs_cut[ipt],"Unlike-sign all","p");
			leg->AddEntry(hmass_rs_cut_bkg[ipt],"Unlike-sign bkg","p");
			leg->AddEntry(hmass_ws_cut[ipt],"Like-sign all","p");
			leg->AddEntry(hmass_cs_cut_mixed[ipt],"Mixed-event","p");
			leg->AddEntry("","Scaled in 3<M<5 GeV/c^{2}","");
			leg->Draw();
		}

		TLine *line = new TLine(1,0,5,0);
		line->SetLineWidth(2);
		line->SetLineStyle(2);
		line->Draw();

		if ( c4 ){
			c4->cd(ipt+1);
			SetPadStyle();

			htmp = (TH1D*)gPad->DrawFrame(1,-0.3*hmass_rs_cut_sig[ipt]->GetMaximum(),5,1.5*hmass_rs_cut_sig[ipt]->GetMaximum());
			SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
			htmp->GetXaxis()->SetTitleOffset(1.0);

			hmass_rs_cut_sig[ipt]->Draw("p same");
			hmass_sub1_cut[ipt]->Draw("p same");
			hmass_sub2_cut[ipt]->Draw("p same");

			massmin = hmass_rs_cut_bkg[ipt]->FindBin(1.3+0.001);
			massmax = hmass_rs_cut_bkg[ipt]->FindBin(2.5-0.001);

			float massmin2 = hmass_rs_cut_bkg[ipt]->FindBin(1.8+0.001);
			float massmax2 = hmass_rs_cut_bkg[ipt]->FindBin(2.5-0.001);

			TLegend *leg = new TLegend(0.45,0.93-0.05*6,0.93,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(legtextsize);
			leg->AddEntry("","PYTHIA8 pp 13 TeV","");
			leg->AddEntry("",multname.c_str(),"");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
			//leg->AddEntry(hmass_rs_cut_sig[ipt],Form("Truth signal, %d (%d)",int(hmass_rs_cut_sig[ipt]->Integral(massmin,massmax)),int(hmass_rs_cut_sig[ipt]->Integral(massmin2,massmax2))),"p");
			//leg->AddEntry(hmass_sub1_cut[ipt],Form("LS subtracted, %d (%d)",int(hmass_sub1_cut[ipt]->Integral(massmin,massmax)),int(hmass_sub1_cut[ipt]->Integral(massmin2,massmax2))),"p");
			//leg->AddEntry(hmass_sub2_cut[ipt],Form("ME subtracted, %d (%d)",int(hmass_sub2_cut[ipt]->Integral(massmin,massmax)),int(hmass_sub2_cut[ipt]->Integral(massmin2,massmax2))),"p");
			leg->AddEntry(hmass_rs_cut_sig[ipt],Form("Truth signal, %d",int(hmass_rs_cut_sig[ipt]->Integral(massmin,massmax))),"p");
			leg->AddEntry(hmass_sub1_cut[ipt],Form("LS subtracted, %d",int(hmass_sub1_cut[ipt]->Integral(massmin,massmax))),"p");
			leg->AddEntry(hmass_sub2_cut[ipt],Form("ME subtracted, %d",int(hmass_sub2_cut[ipt]->Integral(massmin,massmax))),"p");
			leg->Draw();
		}

		if ( c5 ){
			c5->cd(ipt+1);
			SetPadStyle();

			htmp = (TH1D*)gPad->DrawFrame(1,-0.3*hmass_rs_cut_sig[ipt]->GetMaximum(),5,1.5*hmass_rs_cut_sig[ipt]->GetMaximum());
			SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
			htmp->GetXaxis()->SetTitleOffset(2.3);

			line->Draw();

			hmass_rs_cut_sig[ipt]->Draw("p same");
			hmass_sub1_cut[ipt]->Draw("p same");
			hmass_sub2_cut[ipt]->Draw("p same");

			hmass_sig3Xic0_cut[ipt]->SetLineWidth(1);
			hmass_sig4Xicp_cut[ipt]->SetLineWidth(1);

			hmass_sig3Xic0_cut[ipt]->Draw("same");
			hmass_sig4Xicp_cut[ipt]->Draw("same");

			massmin = hmass_rs_cut_bkg[ipt]->FindBin(1.3+0.001);
			massmax = hmass_rs_cut_bkg[ipt]->FindBin(2.5-0.001);

			TLegend *leg = new TLegend(0.45,0.93-0.05*8,0.93,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(legtextsize);
			leg->AddEntry("","PYTHIA8 pp 13 TeV","");
			leg->AddEntry("",multname.c_str(),"");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
			leg->AddEntry(hmass_rs_cut_sig[ipt],Form("Truth signal, %d",int(hmass_rs_cut_sig[ipt]->Integral(massmin,massmax))),"p");
			leg->AddEntry(hmass_sub1_cut[ipt],Form("LS subtracted, %d",int(hmass_sub1_cut[ipt]->Integral(massmin,massmax))),"p");
			leg->AddEntry(hmass_sub2_cut[ipt],Form("ME subtracted, %d",int(hmass_sub2_cut[ipt]->Integral(massmin,massmax))),"p");
			leg->AddEntry(hmass_sig3Xic0_cut[ipt],Form("Truth signal 3body, %d",int(hmass_sig3Xic0_cut[ipt]->Integral(massmin,massmax))),"l");
			leg->AddEntry(hmass_sig4Xicp_cut[ipt],Form("Truth signal 4body, %d",int(hmass_sig4Xicp_cut[ipt]->Integral(massmin,massmax))),"l");
			leg->Draw();
		}

		if ( c7 ){
			c7->cd(ipt+1);
			SetPadStyle();

			//Fit
			TH1D *hdata = (TH1D*)hmass_rs_cut[ipt]->Clone();

			htmp = (TH1D*)gPad->DrawFrame(1,-0.05*hdata->GetMaximum(),5,1.5*hdata->GetMaximum());
			SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
			htmp->GetXaxis()->SetTitleOffset(1.0);

			h3body = (TH1D*)hmass_sig3Xic0_cut[ipt]->Clone();
			h4body = (TH1D*)hmass_sig4Xicp_cut[ipt]->Clone();
			hbkg = (TH1D*)hmass_rs_cut_bkg[ipt]->Clone();
			h3body->Scale(1./h3body->Integral());
			h4body->Scale(1./h4body->Integral());
			hbkg->Scale(1./hbkg->Integral());

			TF1 *func = new TF1("func",fitfunc,1.3,5.0,3);
			func->SetParameters(hdata->Integral(), 0.2, 0.1);
			hdata->Fit("func","R0QL");

			TH1D *hfunc = (TH1D*)hdata->Clone();
			hfunc->Reset();
			for (int ii=0; ii<hfunc->GetNbinsX(); ii++){
				float xx = hfunc->GetBinCenter(ii+1);
				float yy = func->Eval(xx);

				hfunc->SetBinContent(ii+1, yy);
			}
			hfunc->SetLineColor(2);
			hfunc->SetLineWidth(2);

			float n3body = func->GetParameter(0) * func->GetParameter(1);
			float n4body = func->GetParameter(0) * func->GetParameter(2);
			float nbkg = func->GetParameter(0) * (1 - func->GetParameter(1) - func->GetParameter(2));

			TH1D *h3body_fit = (TH1D*)hmass_sig3Xic0_cut[ipt]->Clone();
			TH1D *h4body_fit = (TH1D*)hmass_sig4Xicp_cut[ipt]->Clone();
			TH1D *hbkg_fit = (TH1D*)hmass_rs_cut_bkg[ipt]->Clone();

			h3body_fit->Scale(n3body/h3body_fit->Integral());
			h4body_fit->Scale(n4body/h4body_fit->Integral());
			hbkg_fit->Scale(nbkg/hbkg_fit->Integral());

			h3body_fit->SetLineColor(4);
			h4body_fit->SetLineColor(6);
			hbkg_fit->SetLineColor(8);
			hbkg_fit->SetLineWidth(2);

			hdata->Draw("same");
			h3body_fit->Draw("hist same");
			h4body_fit->Draw("hist same");
			hbkg_fit->Draw("hist same");
			hfunc->Draw("same");

			hmass_sig3Xic0_cut[ipt]->SetMarkerStyle(24);
			hmass_sig3Xic0_cut[ipt]->SetMarkerSize(0.8);
			hmass_sig3Xic0_cut[ipt]->Draw("p same");

			hmass_sig4Xicp_cut[ipt]->SetMarkerStyle(25);
			hmass_sig4Xicp_cut[ipt]->SetMarkerSize(0.8);
			hmass_sig4Xicp_cut[ipt]->Draw("p same");

			hmass_rs_cut_bkg[ipt]->SetMarkerStyle(26);
			hmass_rs_cut_bkg[ipt]->SetMarkerColor(8);
			hmass_rs_cut_bkg[ipt]->SetMarkerSize(0.8);
			hmass_rs_cut_bkg[ipt]->Draw("p same");

			TLegend *leg = new TLegend(0.5,0.93-0.05*10,0.93,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(legtextsize);
			leg->AddEntry("","PYTHIA8 pp 13 TeV","");
			leg->AddEntry("",multname.c_str(),"");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
			leg->AddEntry(hmass_sig3Xic0_cut[ipt],Form("Truth signal 3body, %d",int(hmass_sig3Xic0_cut[ipt]->Integral())),"p");
			leg->AddEntry(hmass_sig4Xicp_cut[ipt],Form("Truth signal 4body, %d",int(hmass_sig4Xicp_cut[ipt]->Integral())),"p");
			leg->AddEntry(hmass_rs_cut_bkg[ipt],Form("Truth background"),"p");
			leg->AddEntry(hfunc,Form("Fit total, %d",int(n3body+n4body)),"l");
			leg->AddEntry(h3body_fit,Form("Fit 3body, %d",int(n3body)),"l");
			leg->AddEntry(h4body_fit,Form("Fit 4body, %d",int(n4body)),"l");
			leg->AddEntry(hbkg_fit,"Fit background","l");
			leg->Draw();

			//summary
			hpt_sig_truth->SetBinContent(ipt+1, hmass_rs_cut_sig[ipt]->Integral()); 
			hpt_3body_truth->SetBinContent(ipt+1, hmass_sig3Xic0_cut[ipt]->Integral()); 
			hpt_3body_fit1->SetBinContent(ipt+1, n3body);
			hpt_sub1->SetBinContent(ipt+1, hmass_sub1_cut[ipt]->Integral(massmin,massmax));
			hpt_sub2->SetBinContent(ipt+1, hmass_sub2_cut[ipt]->Integral(massmin,massmax));
		}

		if ( c8 ){
			c8->cd(ipt+1);
			SetPadStyle();

			//Fit
			TH1D *hdata = (TH1D*)hmass_rs_cut[ipt]->Clone();

			htmp = (TH1D*)gPad->DrawFrame(1,-0.05*hdata->GetMaximum(),5,1.5*hdata->GetMaximum());
			SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
			htmp->GetXaxis()->SetTitleOffset(1.0);

			h3body = (TH1D*)hmass_sig3Xic0_cut[ipt]->Clone();
			h4body = (TH1D*)hmass_sig4Xicp_cut[ipt]->Clone();
			hbkg = (TH1D*)hmass_ws_cut[ipt]->Clone();
			h3body->Scale(1./h3body->Integral());
			h4body->Scale(1./h4body->Integral());
			hbkg->Scale(1./hbkg->Integral());

			TF1 *func = new TF1("func",fitfunc,1.3,5.0,3);
			func->SetParameters(hdata->Integral(), 0.2, 0.1);
			hdata->Fit("func","R0QL");

			TH1D *hfunc = (TH1D*)hdata->Clone();
			hfunc->Reset();
			for (int ii=0; ii<hfunc->GetNbinsX(); ii++){
				float xx = hfunc->GetBinCenter(ii+1);
				float yy = func->Eval(xx);

				hfunc->SetBinContent(ii+1, yy);
			}
			hfunc->SetLineColor(2);
			hfunc->SetLineWidth(2);

			float n3body = func->GetParameter(0) * func->GetParameter(1);
			float n4body = func->GetParameter(0) * func->GetParameter(2);
			float nbkg = func->GetParameter(0) * (1 - func->GetParameter(1) - func->GetParameter(2));

			TH1D *h3body_fit = (TH1D*)hmass_sig3Xic0_cut[ipt]->Clone();
			TH1D *h4body_fit = (TH1D*)hmass_sig4Xicp_cut[ipt]->Clone();
			TH1D *hbkg_fit = (TH1D*)hmass_ws_cut[ipt]->Clone();

			h3body_fit->Scale(n3body/h3body_fit->Integral());
			h4body_fit->Scale(n4body/h4body_fit->Integral());
			hbkg_fit->Scale(nbkg/hbkg_fit->Integral());

			h3body_fit->SetLineColor(4);
			h4body_fit->SetLineColor(6);
			hbkg_fit->SetLineColor(8);
			hbkg_fit->SetLineWidth(2);

			hdata->Draw("same");
			h3body_fit->Draw("hist same");
			h4body_fit->Draw("hist same");
			hbkg_fit->Draw("hist same");
			hfunc->Draw("same");

			hmass_sig3Xic0_cut[ipt]->SetMarkerStyle(24);
			hmass_sig3Xic0_cut[ipt]->SetMarkerSize(0.8);
			hmass_sig3Xic0_cut[ipt]->Draw("p same");

			hmass_sig4Xicp_cut[ipt]->SetMarkerStyle(25);
			hmass_sig4Xicp_cut[ipt]->SetMarkerSize(0.8);
			hmass_sig4Xicp_cut[ipt]->Scale(bSF4BODY);
			hmass_sig4Xicp_cut[ipt]->Draw("p same");

			hmass_rs_cut_bkg[ipt]->SetMarkerStyle(26);
			hmass_rs_cut_bkg[ipt]->SetMarkerColor(1);
			hmass_rs_cut_bkg[ipt]->SetMarkerSize(0.8);
			hmass_rs_cut_bkg[ipt]->Draw("p same");

			TLegend *leg = new TLegend(0.5,0.93-0.05*10,0.93,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(legtextsize);
			leg->AddEntry("","PYTHIA8 pp 13 TeV","");
			leg->AddEntry("",multname.c_str(),"");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
			leg->AddEntry(hmass_sig3Xic0_cut[ipt],Form("Truth signal 3body, %d",int(hmass_sig3Xic0_cut[ipt]->Integral())),"p");
			leg->AddEntry(hmass_sig4Xicp_cut[ipt],Form("Truth signal 4body, %d",int(hmass_sig4Xicp_cut[ipt]->Integral())),"p");
			leg->AddEntry(hmass_rs_cut_bkg[ipt],Form("Truth background"),"p");
			leg->AddEntry(hfunc,Form("Fit total, %d",int(n3body+n4body)),"l");
			leg->AddEntry(h3body_fit,Form("Fit 3body, %d",int(n3body)),"l");
			leg->AddEntry(h4body_fit,Form("Fit 4body, %d",int(n4body)),"l");
			leg->AddEntry(hbkg_fit,"Fit background (LS)","l");
			leg->Draw();

			//summary
			hpt_3body_fit2->SetBinContent(ipt+1, n3body);
		}

		if ( c9 ){
			c9->cd(ipt+1);
			SetPadStyle();

			//Fit
			TH1D *hdata = (TH1D*)hmass_rs_cut[ipt]->Clone();

			htmp = (TH1D*)gPad->DrawFrame(1,-0.05*hdata->GetMaximum(),4,1.5*hdata->GetMaximum());
			SetHistoStyle("Mass (GeV/c^{2})","","",22,20);
			htmp->GetXaxis()->SetTitleOffset(1.0);

			h3body = (TH1D*)hmass_sig3Xic0_cut[ipt]->Clone();
			h4body = (TH1D*)hmass_sig4Xicp_cut[ipt]->Clone();
			hbkg = (TH1D*)hmass_cs_cut_mixed[ipt]->Clone();
			h3body->Scale(1./h3body->Integral());
			h4body->Scale(1./h4body->Integral());
			hbkg->Scale(1./hbkg->Integral());

			TF1 *func = new TF1("func",fitfunc,1.3,4.0,3);
			func->SetParameters(hdata->Integral(), 0.2, 0.1);
			hdata->Fit("func","R0QL");

			TH1D *hfunc = (TH1D*)hdata->Clone();
			hfunc->Reset();
			for (int ii=0; ii<hfunc->GetNbinsX(); ii++){
				float xx = hfunc->GetBinCenter(ii+1);
				float yy = func->Eval(xx);

				hfunc->SetBinContent(ii+1, yy);
			}
			hfunc->SetLineColor(2);
			hfunc->SetLineWidth(2);

			float n3body = func->GetParameter(0) * func->GetParameter(1);
			float n4body = func->GetParameter(0) * func->GetParameter(2);
			float nbkg = func->GetParameter(0) * (1 - func->GetParameter(1) - func->GetParameter(2));

			TH1D *h3body_fit = (TH1D*)hmass_sig3Xic0_cut[ipt]->Clone();
			TH1D *h4body_fit = (TH1D*)hmass_sig4Xicp_cut[ipt]->Clone();
			TH1D *hbkg_fit = (TH1D*)hmass_cs_cut_mixed[ipt]->Clone();

			h3body_fit->Scale(n3body/h3body_fit->Integral());
			h4body_fit->Scale(n4body/h4body_fit->Integral());
			hbkg_fit->Scale(nbkg/hbkg_fit->Integral());

			h3body_fit->SetLineColor(4);
			h4body_fit->SetLineColor(6);
			hbkg_fit->SetLineColor(8);
			hbkg_fit->SetLineWidth(2);

			hdata->Draw("same");
			h3body_fit->Draw("hist same");
			h4body_fit->Draw("hist same");
			hbkg_fit->Draw("hist same");
			hfunc->Draw("same");

			hmass_sig3Xic0_cut[ipt]->SetMarkerStyle(24);
			hmass_sig3Xic0_cut[ipt]->SetMarkerSize(0.8);
			hmass_sig3Xic0_cut[ipt]->Draw("p same");

			hmass_sig4Xicp_cut[ipt]->SetMarkerStyle(25);
			hmass_sig4Xicp_cut[ipt]->SetMarkerSize(0.8);
			hmass_sig4Xicp_cut[ipt]->Draw("p same");

			hmass_rs_cut_bkg[ipt]->SetMarkerStyle(26);
			hmass_rs_cut_bkg[ipt]->SetMarkerColor(1);
			hmass_rs_cut_bkg[ipt]->SetMarkerSize(0.8);
			hmass_rs_cut_bkg[ipt]->Draw("p same");

			TLegend *leg = new TLegend(0.5,0.93-0.05*10,0.93,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(legtextsize);
			leg->AddEntry("","PYTHIA8 pp 13 TeV","");
			leg->AddEntry("",multname.c_str(),"");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
			leg->AddEntry(hmass_sig3Xic0_cut[ipt],Form("Truth signal 3body, %d",int(hmass_sig3Xic0_cut[ipt]->Integral())),"p");
			leg->AddEntry(hmass_sig4Xicp_cut[ipt],Form("Truth signal 4body, %d",int(hmass_sig4Xicp_cut[ipt]->Integral())),"p");
			leg->AddEntry(hmass_rs_cut_bkg[ipt],Form("Truth background"),"p");
			leg->AddEntry(hfunc,Form("Fit total, %d",int(n3body+n4body)),"l");
			leg->AddEntry(h3body_fit,Form("Fit 3body, %d",int(n3body)),"l");
			leg->AddEntry(h4body_fit,Form("Fit 4body, %d",int(n4body)),"l");
			leg->AddEntry(hbkg_fit,"Fit background (ME)","l");
			leg->Draw();

			//summary
			hpt_3body_fit3->SetBinContent(ipt+1, n3body);
		}

	}//ipt

	//return;
	hpt_3body_truth->SetMarkerStyle(24);
	hpt_3body_truth->SetMarkerSize(0.8);
	hpt_3body_truth->SetLineColor(1);
	hpt_3body_truth->Sumw2();

	hpt_3body_fit1->SetMarkerStyle(25);
	hpt_3body_fit1->SetMarkerSize(0.8);
	hpt_3body_fit1->SetLineColor(2);

	hpt_3body_fit2->SetMarkerStyle(26);
	hpt_3body_fit2->SetMarkerColor(4);
	hpt_3body_fit2->SetMarkerSize(0.8);
	hpt_3body_fit2->SetLineColor(4);

	hpt_3body_fit3->SetMarkerStyle(32);
	hpt_3body_fit3->SetMarkerColor(6);
	hpt_3body_fit3->SetMarkerSize(0.8);
	hpt_3body_fit3->SetLineColor(6);

	hpt_sub1->SetMarkerStyle(26);
	hpt_sub1->SetLineColor(4);

	hpt_sub2->SetMarkerStyle(32);
	hpt_sub2->SetLineColor(6);

	TH1D *hpt_ratio_3body_fit1 = (TH1D*)hpt_3body_fit1->Clone();
	TH1D *hpt_ratio_3body_fit2 = (TH1D*)hpt_3body_fit2->Clone();
	TH1D *hpt_ratio_3body_fit3 = (TH1D*)hpt_3body_fit3->Clone();
	TH1D *hpt_ratio_sub1 = (TH1D*)hpt_sub1->Clone();
	TH1D *hpt_ratio_sub2 = (TH1D*)hpt_sub2->Clone();

	hpt_ratio_3body_fit1->Divide(hpt_3body_truth);
	hpt_ratio_3body_fit2->Divide(hpt_3body_truth);
	hpt_ratio_3body_fit3->Divide(hpt_3body_truth);
	hpt_ratio_sub1->Divide(hpt_3body_truth);
	hpt_ratio_sub2->Divide(hpt_3body_truth);


	TCanvas *c10 = new TCanvas("c10","c10",1.2*2*500,500);
	c10->Divide(2,1);

	{
		c10->cd(1);
		SetPadStyle();
		gPad->SetLogy();

		hpt_sig_truth->Scale(1,"width");
		hpt_3body_truth->Scale(1,"width");
		hpt_3body_fit1->Scale(1,"width");
		hpt_3body_fit2->Scale(1,"width");
		hpt_3body_fit3->Scale(1,"width");
		hpt_sub1->Scale(1,"width");
		hpt_sub2->Scale(1,"width");

		htmp = (TH1D*)gPad->DrawFrame(0,5,24,10*hpt_3body_truth->GetMaximum());
		SetHistoStyle("p_{T} (GeV/c)","dN/dp_{T} (GeV/c)^{-1}","",22,20);
		htmp->GetXaxis()->SetTitleOffset(1.3);

		//hpt_sig_truth->Draw("same");
		hpt_3body_truth->Draw("same");
		hpt_3body_fit1->Draw("same");
		hpt_3body_fit2->Draw("same");
		//hpt_3body_fit3->Draw("same");
		//hpt_3body_fit2->Draw("same");
		//hpt_3body_fit3->Draw("same");
		//hpt_sub1->Draw("same");
		//hpt_sub2->Draw("same");

		TLegend *leg = new TLegend(0.5,0.93-0.05*7,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("",multname.c_str(),"h");
		leg->AddEntry("","1.3<e#Xi mass<2.5 GeV/c^{2}","h");
		leg->AddEntry(hpt_3body_truth,"Truth signal 3body","PL");
		leg->AddEntry(hpt_3body_fit1,"Fit 3body (US template)","PL");
		leg->AddEntry(hpt_3body_fit2,"Fit 3body (LS template)","PL");
		//leg->AddEntry(hpt_3body_fit3,"Fit 3body (ME template)","PL");
		//leg->AddEntry(hpt_sub1,"LS subtraction only","PL");
		//leg->AddEntry(hpt_sub2,"ME subtraction only","PL");
		leg->Draw();
	}

	{
		c10->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,24,2);
		SetHistoStyle("p_{T} (GeV/c)","Ratio","",22,20);
		htmp->GetXaxis()->SetTitleOffset(1.3);

		hpt_ratio_3body_fit1->Draw("p same");
		hpt_ratio_3body_fit2->Draw("p same");
		//hpt_ratio_3body_fit3->Draw("p same");
		//hpt_ratio_sub1->Draw("p same");
		//hpt_ratio_sub2->Draw("p same");
	}


	/*
	TCanvas *c5 = new TCanvas("c5","c5",1.2*3*500,2*500);
	c5->Divide(3,2);

	{
		c5->cd(1);
		SetPadStyle(1);
	
		htmp = (TH1D*)gPad->DrawFrame(0,0,10,180);
		SetHistoStyle("p_{T} (GeV/c)","Opening angle (deg)","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		htmp->GetYaxis()->SetTitleOffset(2.3);

		hoa_pt_sig3Xic0[0]->Draw("colz same");
	}

	{
		c5->cd(4);
		SetPadStyle(1);
	
		htmp = (TH1D*)gPad->DrawFrame(0,0,10,180);
		SetHistoStyle("p_{T} (GeV/c)","Opening angle (deg)","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		htmp->GetYaxis()->SetTitleOffset(2.3);

		hoa_pt_sig4Xicp[0]->Draw("colz same");
	}

	{
		c5->cd(2);
		SetPadStyle(1);
	
		htmp = (TH1D*)gPad->DrawFrame(0,0,10,180);
		SetHistoStyle("p_{T} (GeV/c)","Opening angle (deg)","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		htmp->GetYaxis()->SetTitleOffset(2.3);

		hoa_pt_comb[0]->Draw("colz same");
	}

	{
		c5->cd(5);
		SetPadStyle(1);
	
		htmp = (TH1D*)gPad->DrawFrame(0,0,10,180);
		SetHistoStyle("p_{T} (GeV/c)","Opening angle (deg)","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		htmp->GetYaxis()->SetTitleOffset(2.3);

		hoa_pt_comb[1]->Draw("colz same");
	}

	{
		c5->cd(3);
		SetPadStyle(1);
	
		htmp = (TH1D*)gPad->DrawFrame(0,0,10,180);
		SetHistoStyle("p_{T} (GeV/c)","Opening angle (deg)","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		htmp->GetYaxis()->SetTitleOffset(2.3);

		hoa_pt_corrJet[0]->Draw("colz same");
	}

	{
		c5->cd(6);
		SetPadStyle(1);
	
		htmp = (TH1D*)gPad->DrawFrame(0,0,10,180);
		SetHistoStyle("p_{T} (GeV/c)","Opening angle (deg)","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);
		htmp->GetYaxis()->SetTitleOffset(2.3);

		hoa_pt_corrJet[1]->Draw("colz same");
	}
	*/

	//return;

	/*
	TCanvas *c6 = new TCanvas("c6","c6",1.2*3*500,2*500);
	c6->Divide(3,2);

	for (int ipt=0; ipt<npt; ipt++){
		c6->cd(ipt+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,180,2*hoa_sig3Xic0[0][ipt]->GetMaximum());
		SetHistoStyle("Opening angle (deg)","","",22,20);
		htmp->GetXaxis()->SetTitleOffset(2.3);

		hoa_sig3Xic0[0][ipt]->Draw("p same");
		hoa_sig4Xicp[0][ipt]->Draw("p same");
		hoa_comb[0][ipt]->Draw("p same");
		hoa_comb[1][ipt]->Draw("p same");
		hoa_mixed1[0][ipt]->Draw("p same");
		hoa_mixed1[1][ipt]->Draw("p same");
		//hoa_mixed2[0][ipt]->Draw("p same");
		//hoa_mixed2[1][ipt]->Draw("p same");
		hoa_corrJet[0][ipt]->Draw("p same");
		hoa_corrJet[1][ipt]->Draw("p same");

		TLegend *leg = new TLegend(0.4,0.93-0.05*9,0.93,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(legtextsize);
		leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
		leg->AddEntry(hoa_sig3Xic0[0][ipt],"#Xi_{c}^{0} 3body","p");
		leg->AddEntry(hoa_sig4Xicp[0][ipt],"#Xi_{c}^{+} 4body","p");
		leg->AddEntry(hoa_comb[0][ipt],"Combinatorial bkg, Unlike-sign","p");
		leg->AddEntry(hoa_comb[1][ipt],"Combinatorial bkg, Like-sign","p");
		leg->AddEntry(hoa_mixed1[0][ipt],"Mixed events, Unlike-sign","p");
		leg->AddEntry(hoa_mixed1[1][ipt],"Mixed events, Like-sign","p");
		leg->AddEntry(hoa_corrJet[0][ipt],"Correlated bkg, Unlike-sign","p");
		leg->AddEntry(hoa_corrJet[1][ipt],"Correlated bkg, Like-sign","p");
		leg->Draw();
	}
	*/


	if ( bSAVE ){

		//c1->cd();
		//c1->SaveAs("PYTHIA8pp13TeV_Scan80_c1.eps");

		//c2->cd();
		//c2->SaveAs("PYTHIA8pp13TeV_Scan80_c2.eps");

		c3->cd();
		c3->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c3.png",opt));
		c3->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c3.pdf",opt));

		c4->cd();
		c4->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c4.png",opt));
		c4->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c4.pdf",opt));

		//c5->cd();
		//c5->SaveAs("PYTHIA8pp13TeV_Scan80_c5.eps");

		c7->cd();
		c7->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c7.png",opt));
		c7->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c7.pdf",opt));

		c8->cd();
		c8->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c8.png",opt));
		c8->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c8.pdf",opt));

		c9->cd();
		c9->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c9.png",opt));
		c9->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c9.pdf",opt));

		c10->cd();
		c10->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c10.png",opt));
		c10->SaveAs(Form("PYTHIA8pp13TeV_Scan85_opt%d_c10.pdf",opt));
	}

}

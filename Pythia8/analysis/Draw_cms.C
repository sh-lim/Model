#include "Style.h"

void Draw_cms(const char *dataset="pp13TeV_nondiffractive_grp0_try3", const bool bSAVE=false){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TLegendEntry *le;

  const char *system = "p+p";

	TFile *infile = new TFile(Form("outfile_hist_%s.root",dataset),"read");

	const int nmult = 15;
	const float const_pi = TMath::Pi();

	const int nref_1 = 0;
	const int nref_2 = 1;
	const int nref_3 = 2;

	TH2D *h2d_same_dphi_deta[nmult];
	TH2D *h2d_mixed_dphi_deta[nmult];

	TH1D *h1d_same_dphi_long[nmult][2];
	TH1D *h1d_mixed_dphi_long[nmult][2];

	TH1D *h1d_same_dphi_short[nmult][2];
	TH1D *h1d_mixed_dphi_short[nmult][2];

	TH1D *h1d_same_dphi_sub[nmult];
	TH1D *h1d_same_dphi_sub_zyam[nmult];

	TH1D *hntrig = (TH1D*)infile->Get("hntrig_mid");
	TH1D *hntrig_mixed = (TH1D*)infile->Get("hntrig_mixed_mid");

	for (int im=0; im<nmult; im++){
		h2d_same_dphi_deta[im] = (TH2D*)infile->Get(Form("h2d_same_dphi_deta_m%02d",im));
		h2d_mixed_dphi_deta[im] = (TH2D*)infile->Get(Form("h2d_mixed_dphi_deta_m%02d",im));

    h2d_same_dphi_deta[im]->RebinX(3);
    h2d_mixed_dphi_deta[im]->RebinX(3);
    h2d_same_dphi_deta[im]->RebinY(2);
    h2d_mixed_dphi_deta[im]->RebinY(2);
    h2d_same_dphi_deta[im]->Sumw2();
    h2d_mixed_dphi_deta[im]->Sumw2();
		h2d_same_dphi_deta[im]->Divide(h2d_mixed_dphi_deta[im]);

		float ntrig = hntrig->Integral(10*im+1, 10*im+10);
		float ntrig_mixed = hntrig_mixed->Integral(10*im+1, 10*im+10);
    float norm = h2d_mixed_dphi_deta[im]->GetBinContent(h2d_mixed_dphi_deta[im]->FindBin(0,0));
    h2d_same_dphi_deta[im]->Scale(norm/ntrig*ntrig_mixed/ntrig);

		int etamin = h2d_same_dphi_deta[im]->GetYaxis()->FindBin(-5.0+0.0001);
		int etamax = h2d_same_dphi_deta[im]->GetYaxis()->FindBin(-2.0-0.0001);

		h1d_same_dphi_long[im][0] = (TH1D*)h2d_same_dphi_deta[im]->ProjectionX(Form("h1d_same_dphi_long_%d_0",im),etamin,etamax);

		etamin = h2d_same_dphi_deta[im]->GetYaxis()->FindBin(+2.0+0.0001);
		etamax = h2d_same_dphi_deta[im]->GetYaxis()->FindBin(+5.0-0.0001);
		h1d_same_dphi_long[im][1] = (TH1D*)h2d_same_dphi_deta[im]->ProjectionX(Form("h1d_same_dphi_long_%d_1",im),etamin,etamax);

		h1d_same_dphi_long[im][0]->Add(h1d_same_dphi_long[im][1]);
		h1d_same_dphi_long[im][0]->Scale(1./6.0/h1d_same_dphi_long[im][0]->GetBinWidth(1));

		etamin = h2d_same_dphi_deta[im]->GetYaxis()->FindBin(-1.0+0.0001);
		etamax = h2d_same_dphi_deta[im]->GetYaxis()->FindBin(+1.0-0.0001);

		h1d_same_dphi_short[im][0] = (TH1D*)h2d_same_dphi_deta[im]->ProjectionX(Form("h1d_same_dphi_short_%d_0",im),etamin,etamax);
		h1d_same_dphi_short[im][0]->Scale(1./2.0/h1d_same_dphi_short[im][0]->GetBinWidth(1));

		h1d_same_dphi_sub[im] = (TH1D*)h1d_same_dphi_short[im][0]->Clone(Form("h1d_same_dphi_sub_%d",im));
		h1d_same_dphi_sub[im]->Add(h1d_same_dphi_long[im][0],-1);

	}//

  //return;

	float Nassoc[nmult];

	TH1D *hv22_dir = new TH1D("hv22_dir","",20,0,200);

	TGraphErrors *gv22_dir = new TGraphErrors;
	TGraphErrors *gv22_sub_nref_1 = new TGraphErrors;
	TGraphErrors *gv22_sub_nref_2 = new TGraphErrors;
	TGraphErrors *gv22_sub_nref_3 = new TGraphErrors;

	TH1D *hcoeff_dir[5];
	for (int io=0; io<5; io++){
		hcoeff_dir[io] = new TH1D(Form("hcoeff_dir_%d",io),"",20,0,200);
	}

	TF1 *fdir[nmult];

	TCanvas *c1 = new TCanvas("c1","c1",1.1*300*4,300*3);
	c1->Divide(4,3);

	TCanvas *c2 = new TCanvas("c2","c2",1.1*250*4,250*3);
	c2->Divide(4,3);

	for (int im=0; im<nmult; im++){

		fdir[im] = new TF1(Form("fdir_%d",im),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x) + 2*[5]*cos(5*x))",-const_pi/2,3*const_pi/2);
		fdir[im]->SetLineWidth(3);
		fdir[im]->SetLineStyle(2);
		fdir[im]->SetLineColor(1);

		h1d_same_dphi_long[im][0]->Fit(fdir[im], "R0Q");

		cout << "Mult: " << im
			<< ", cos(x): " << fdir[im]->GetParameter(1)
			<< ", cos(2x): " << fdir[im]->GetParameter(2)
			<< ", cos(3x): " << fdir[im]->GetParameter(3)
			<< ", cos(4x): " << fdir[im]->GetParameter(4)
			<< ", cos(5x): " << fdir[im]->GetParameter(5)
			<< endl;

		for (int io=0; io<5; io++){
			hcoeff_dir[io]->SetBinContent(im+1, fdir[im]->GetParameter(io+1));
			hcoeff_dir[io]->SetBinError(im+1, fdir[im]->GetParError(io+1));
		}

		hv22_dir->SetBinContent(im+1, fdir[im]->GetParameter(2));
		hv22_dir->SetBinError(im+1, fdir[im]->GetParError(2));
		Nassoc[im] = fdir[im]->GetParameter(0)*2*const_pi;

		if ( im>=12 ) continue;

		c1->cd(im+1);
		SetPadStyle(0.20,0.01,0.14,0.01);

		float max = TMath::Max(h1d_same_dphi_long[im][0]->GetMaximum(),h1d_same_dphi_short[im][0]->GetMaximum());
		float min = TMath::Min(h1d_same_dphi_long[im][0]->GetMinimum(),h1d_same_dphi_short[im][0]->GetMinimum());

		//float max = TMath::Max(h1d_same_dphi_long[im][0]->GetMaximum(),h1d_same_dphi_long[im][0]->GetMaximum());
		//float min = TMath::Min(h1d_same_dphi_long[im][0]->GetMinimum(),h1d_same_dphi_long[im][0]->GetMinimum());

    //htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, 0.98*min, 3*const_pi/2, 1.05*max);
    htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.1*(max-min), 3*const_pi/2, max+0.6*(max-min));
		SetHistoStyle("#Delta#phi (rad)","(1/N^{Trig}) dN^{Pair}/(d#Delta#phi)","",0.08);
		htmp->GetXaxis()->CenterTitle();
		htmp->GetXaxis()->SetTitleOffset(0.9);
		htmp->GetYaxis()->CenterTitle();
		htmp->GetYaxis()->SetNdivisions(5,3,0);
		htmp->GetYaxis()->SetTitleOffset(1.15);
		h1d_same_dphi_long[im][0]->SetMarkerStyle(24);
		h1d_same_dphi_long[im][0]->Draw("p same");

		h1d_same_dphi_short[im][0]->SetMarkerStyle(24);
		h1d_same_dphi_short[im][0]->SetMarkerColor(2);
		h1d_same_dphi_short[im][0]->Draw("p same");

		fdir[im]->SetLineWidth(3);
		fdir[im]->Draw("same");

		if ( im==0 ){
			TLegend *leg = new TLegend(0.10,0.80,0.8,0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("",Form("h^{#pm}-h^{#pm}"),"");
			le->SetTextFont(62);
			le->SetTextSize(0.07);
			le = leg->AddEntry("",Form("PYTHIA8 %s 13 TeV",system),"");
			le->SetTextFont(62);
			le->SetTextSize(0.07);
			leg->Draw();
    }else if ( im==1 ){
			TLegend *leg = new TLegend(0.05,0.75,0.8,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry("",Form("0.5<p_{T}^{Trig}, p_{T}^{Assoc}<5.0 GeV/c"),"");
			le->SetTextSize(0.065);
			le = leg->AddEntry("",Form("0<|#eta|^{Trig}, |#eta|^{Assoc}<2.5"),"");
			le->SetTextSize(0.065);
			leg->Draw();
    }else if ( im==2 ){
			TLegend *leg = new TLegend(0.18,0.70,0.9,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			le = leg->AddEntry(h1d_same_dphi_long[im][0],"Long-range 2<|#Delta#eta|<5","P");
			le->SetTextSize(0.065);
			le = leg->AddEntry(h1d_same_dphi_short[im][0],"Short-range 0<|#Delta#eta|<1","P");
			le->SetTextSize(0.065);
			le = leg->AddEntry(fdir[im],"Fourier fits","L");
			le->SetTextSize(0.065);
			leg->Draw();
		}

		{
			TLatex *tex = new TLatex(4.5,max,Form("%d#leqN^{ch}<%d",10*im,10*(im+1)));
			tex->SetTextFont(62);
			tex->SetTextAlign(32);
			tex->SetTextSize(0.07);
			tex->Draw();
		}

		{
			TLatex *tex = new TLatex(4.5,max+0.15*(max-min),Form("(%c)",97+im));
			tex->SetTextFont(62);
			tex->SetTextAlign(32);
			tex->SetTextSize(0.07);
			tex->Draw();
		}

		c2->cd(im+1);
		SetPadStyle();

		max = h1d_same_dphi_sub[im]->GetMaximum();
		min = h1d_same_dphi_sub[im]->GetMinimum();

		htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, -1.1*fabs(min), 3*const_pi/2, 1.1*max);
		SetHistoStyle("","");

		h1d_same_dphi_sub[im]->SetMarkerStyle(24);
		h1d_same_dphi_sub[im]->Draw("p same");
	}//

	//return;

  int g_count = 0;
	for (int im=nref_1+2; im<nmult; im++){
		float R_assoc = Nassoc[nref_1]/Nassoc[im];

		float Y_peri = h1d_same_dphi_sub[nref_1]->Integral(h1d_same_dphi_sub[nref_1]->FindBin(-1.2),h1d_same_dphi_sub[nref_1]->FindBin(+1.2)); 
		float Y_cent = h1d_same_dphi_sub[im]->Integral(h1d_same_dphi_sub[im]->FindBin(-1.2),h1d_same_dphi_sub[im]->FindBin(+1.2)); 

		float v22_sub = hv22_dir->GetBinContent(im+1) - hv22_dir->GetBinContent(nref_1+1)*R_assoc*Y_cent/Y_peri;

		gv22_sub_nref_1->SetPoint(g_count, 10*im+5, v22_sub);
		gv22_sub_nref_1->SetPointError(g_count, 0, 0);

		cout << "mult: " << im 
			<< ", R_assoc: " << R_assoc 
			<< ", Y_cent: " << Y_cent 
			<< ", Y_peri: " << Y_peri 
			<< ", R_jet: " << Y_cent/Y_peri 
			<< ", v22 cent: " << hv22_dir->GetBinContent(im+1) 
			<< ", v22 peri: " << hv22_dir->GetBinContent(nref_1+1) 
			<< ", v22_sub: " << v22_sub
			<< endl;
    g_count++;
	}

  g_count = 0;
	for (int im=nref_2+1; im<nmult; im++){
		float R_assoc = Nassoc[nref_2]/Nassoc[im];

		float Y_peri = h1d_same_dphi_sub[nref_2]->Integral(h1d_same_dphi_sub[nref_2]->FindBin(-1.2),h1d_same_dphi_sub[nref_2]->FindBin(+1.2)); 
		float Y_cent = h1d_same_dphi_sub[im]->Integral(h1d_same_dphi_sub[im]->FindBin(-1.2),h1d_same_dphi_sub[im]->FindBin(+1.2)); 

		float v22_sub = hv22_dir->GetBinContent(im+1) - hv22_dir->GetBinContent(nref_2+1)*R_assoc*Y_cent/Y_peri;

		gv22_sub_nref_2->SetPoint(g_count, 10*im+5, v22_sub);
		gv22_sub_nref_2->SetPointError(g_count, 0, 0);

		cout << "mult: " << im 
			<< ", R_assoc: " << R_assoc 
			<< ", Y_cent: " << Y_cent 
			<< ", Y_peri: " << Y_peri 
			<< ", R_jet: " << Y_cent/Y_peri 
			<< ", v22 cent: " << hv22_dir->GetBinContent(im+1) 
			<< ", v22 peri: " << hv22_dir->GetBinContent(nref_2+1) 
			<< ", v22_sub: " << v22_sub
			<< endl;
    g_count++;
	}

  g_count = 0;
	for (int im=nref_3+1; im<nmult; im++){
		float R_assoc = Nassoc[nref_3]/Nassoc[im];

		float Y_peri = h1d_same_dphi_sub[nref_3]->Integral(h1d_same_dphi_sub[nref_3]->FindBin(-1.2),h1d_same_dphi_sub[nref_3]->FindBin(+1.2)); 
		float Y_cent = h1d_same_dphi_sub[im]->Integral(h1d_same_dphi_sub[im]->FindBin(-1.2),h1d_same_dphi_sub[im]->FindBin(+1.2)); 

		float v22_sub = hv22_dir->GetBinContent(im+1) - hv22_dir->GetBinContent(nref_3+1)*R_assoc*Y_cent/Y_peri;

		gv22_sub_nref_3->SetPoint(im-nref_3-1, 10*im+5, v22_sub);
		gv22_sub_nref_3->SetPointError(im-nref_3-1, 0, 0);

		cout << "mult: " << im 
			<< ", R_assoc: " << R_assoc 
			<< ", Y_cent: " << Y_cent 
			<< ", Y_peri: " << Y_peri 
			<< ", R_jet: " << Y_cent/Y_peri 
			<< ", v22 cent: " << hv22_dir->GetBinContent(im+1) 
			<< ", v22 peri: " << hv22_dir->GetBinContent(nref_3+1) 
			<< ", v22_sub: " << v22_sub
			<< endl;
    g_count++;
	}

	//return;

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	SetPadStyle();
	gPad->SetBottomMargin(0.13);
	gPad->SetLeftMargin(0.16);

	htmp = (TH1D*)gPad->DrawFrame(0,-0.003,150,0.016);
	SetHistoStyle("N^{ch}","v_{22}","",0.065);
	htmp->GetXaxis()->SetTitleOffset(0.9);

  TBox *box = new TBox(0,-0.03*0.03,150,0.03*0.03);
  box->SetFillStyle(1001);
  box->SetFillColorAlpha(11,0.5);
  box->Draw("");

	TLine *line = new TLine(0,0,150,0);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->SetLineColor(1);
	line->Draw();

	hv22_dir->SetMarkerStyle(24);
	hv22_dir->SetMarkerColor(1);
	hv22_dir->SetLineColor(1);
	hv22_dir->SetLineWidth(2);
	hv22_dir->Draw("p same");

	gv22_sub_nref_1->SetLineWidth(4);
	gv22_sub_nref_1->SetLineColor(1);
	gv22_sub_nref_1->Draw("c same");

	gv22_sub_nref_2->SetLineWidth(4);
	gv22_sub_nref_2->SetLineColor(2);
	//gv22_sub_nref_2->SetLineStyle(2);
	gv22_sub_nref_2->Draw("c same");

	gv22_sub_nref_3->SetLineWidth(4);
	gv22_sub_nref_3->SetLineColor(4);
	//gv22_sub_nref_3->SetLineStyle(5);
	gv22_sub_nref_3->Draw("c same");

	{
		TLegend *leg = new TLegend(0.45,0.45,0.85,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		le = leg->AddEntry("","PYTHIA8 p+p 13 TeV","");
		le->SetTextFont(62);
		le->SetTextSize(0.05);
		le = leg->AddEntry("","Method 2","");
		le->SetTextFont(32);
		le->SetTextSize(0.05);
		le = leg->AddEntry("","0.5<p_{T}^{Trig}, p_{T}^{Assoc}<5.0 GeV/c","");
		le->SetTextSize(0.045);
		le = leg->AddEntry("","0<|#eta^{Trig}|, |#eta^{Assoc}|<2.5","");
		le->SetTextSize(0.045);
		le = leg->AddEntry(hv22_dir,"Fourier fits","P");
		le->SetTextSize(0.045);
		le = leg->AddEntry(gv22_sub_nref_1,Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_1*10,nref_1*10+10),"L");
		le->SetTextSize(0.045);
		le = leg->AddEntry(gv22_sub_nref_2,Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_2*10,nref_2*10+10),"L");
		le->SetTextSize(0.045);
		le = leg->AddEntry(gv22_sub_nref_3,Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_3*10,nref_3*10+10),"L");
		le->SetTextSize(0.045);
		leg->Draw();
	}

	{
		TLatex *tex = new TLatex(10,0.014,"(b)");
		tex->SetTextSize(0.045);
		tex->SetTextFont(62);
		tex->Draw();
	}

	//TProfile *hv22pc_mult = (TProfile*)infile->Get("hv22pc_mult");
	//hv22pc_mult->Draw("same");

	/*
	TCanvas *c4 = new TCanvas("c4","c4",1.2*500,500);
	SetPadStyle();
	gPad->SetBottomMargin(0.13);
	htmp = (TH1D*)gPad->DrawFrame(0,-0.06,150,0.015);
	SetHistoStyle("N_{TRK}^{CH}","w_{n}");
	hcoeff_dir[0]->SetMarkerStyle(20);
	hcoeff_dir[0]->SetLineColor(1);
	hcoeff_dir[0]->SetLineWidth(2);
	hcoeff_dir[0]->Draw("p same");
	hcoeff_dir[1]->SetMarkerStyle(24);
	hcoeff_dir[1]->SetLineColor(2);
	hcoeff_dir[1]->SetLineWidth(2);
	hcoeff_dir[1]->Draw("p same");

	{
		TLegend *leg = new TLegend(0.3,0.25,0.9,0.4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry("","Pythia p+p 13 TeV","");
		leg->AddEntry("","N_{TRK}^{CH}: p_{T}>0.4 GeV, |#eta|<2.5","");
		leg->AddEntry("","0.5<p_{T}<5.0 GeV, 2.0<|#Delta#eta|<5.0","");
		leg->Draw();
	}
	*/

  if ( bSAVE ){
    c1->cd();
    c1->SaveAs(Form("%s_1d.gif",dataset));

    c3->cd();
    c3->SaveAs(Form("%s_v22.gif",dataset));
  }

}

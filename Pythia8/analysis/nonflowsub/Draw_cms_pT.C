#include "Style.h"

void Draw_cms_pT(const char *dataset="pp13TeV_nondiffractive_grp0_try4", const bool bSAVE=false){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TLegendEntry *le;

  const char *system = "p+p";

	TFile *infile = new TFile(Form("outfile_hist_%s.root",dataset),"read");

	const int nmult = 5;
	const float multcut_mid[nmult+1] = {0, 10, 20, 30, 85, 500};
	const float const_pi = TMath::Pi();

  const int npt = 8;
  const float ptcut[npt+1] = {0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};

	const int nref_1 = 0;
	const int nref_2 = 1;
	const int nref_3 = 2;

	TH2D *h2d_same_dphi_deta[nmult][npt];
	TH2D *h2d_mixed_dphi_deta[nmult][npt];

	TH1D *h1d_same_dphi_long[nmult][npt][2];
	TH1D *h1d_mixed_dphi_long[nmult][npt][2];

	TH1D *h1d_same_dphi_short[nmult][npt][2];
	TH1D *h1d_mixed_dphi_short[nmult][npt][2];

	TH1D *h1d_same_dphi_sub[nmult][npt];

	TH2D *hntrig = (TH2D*)infile->Get("hntrig_mid");
	TH2D *hntrig_mixed = (TH2D*)infile->Get("hntrig_mixed_mid");

	for (int im=0; im<nmult; im++){
		for (int ipt=0; ipt<npt; ipt++){
			if ( ipt==0 ) continue;

			h2d_same_dphi_deta[im][ipt] = (TH2D*)infile->Get(Form("h2d_same_dphi_deta_m%02d_pt%02d",im,ipt));
			h2d_mixed_dphi_deta[im][ipt] = (TH2D*)infile->Get(Form("h2d_mixed_dphi_deta_m%02d_pt%02d",im,ipt));

			h2d_same_dphi_deta[im][ipt]->RebinX(3);
			h2d_mixed_dphi_deta[im][ipt]->RebinX(3);
			h2d_same_dphi_deta[im][ipt]->RebinY(2);
			h2d_mixed_dphi_deta[im][ipt]->RebinY(2);
			h2d_same_dphi_deta[im][ipt]->Sumw2();
			h2d_mixed_dphi_deta[im][ipt]->Sumw2();
			h2d_same_dphi_deta[im][ipt]->Divide(h2d_mixed_dphi_deta[im][ipt]);

			float ntrig = hntrig->GetBinContent(im+1, ipt+1);
			float ntrig_mixed = hntrig_mixed->GetBinContent(im+1, ipt+1);
			float norm = h2d_mixed_dphi_deta[im][ipt]->GetBinContent(h2d_mixed_dphi_deta[im][ipt]->FindBin(0,0));
			h2d_same_dphi_deta[im][ipt]->Scale(norm/ntrig*ntrig_mixed/ntrig);

			int etamin = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(-5.0+0.0001);
			int etamax = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(-2.0-0.0001);
			h1d_same_dphi_long[im][ipt][0] = (TH1D*)h2d_same_dphi_deta[im][ipt]->ProjectionX(Form("h1d_same_dphi_long_%d_%d_0",im,ipt),etamin,etamax);

			etamin = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(+2.0+0.0001);
			etamax = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(+5.0-0.0001);
			h1d_same_dphi_long[im][ipt][1] = (TH1D*)h2d_same_dphi_deta[im][ipt]->ProjectionX(Form("h1d_same_dphi_long_%d_%d_1",im,ipt),etamin,etamax);

			h1d_same_dphi_long[im][ipt][0]->Add(h1d_same_dphi_long[im][ipt][1]);
			h1d_same_dphi_long[im][ipt][0]->Scale(1./6.0/h1d_same_dphi_long[im][ipt][0]->GetBinWidth(1));

			etamin = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(-1.0+0.0001);
			etamax = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(+1.0-0.0001);

			h1d_same_dphi_short[im][ipt][0] = (TH1D*)h2d_same_dphi_deta[im][ipt]->ProjectionX(Form("h1d_same_dphi_short_%d_%d_0",im,ipt),etamin,etamax);
			h1d_same_dphi_short[im][ipt][0]->Scale(1./2.0/h1d_same_dphi_short[im][ipt][0]->GetBinWidth(1));

			h1d_same_dphi_sub[im][ipt] = (TH1D*)h1d_same_dphi_short[im][ipt][0]->Clone(Form("h1d_same_dphi_sub_%d_%d",im,ipt));
			h1d_same_dphi_sub[im][ipt]->Add(h1d_same_dphi_long[im][ipt][0],-1);
		}
	}//

  TCanvas *c1[nmult];
  TCanvas *c2[nmult];

	float Nassoc[nmult][npt];
	TH1D *hv22_dir[nmult];
  TGraphErrors *gv22_sub_nref_1[nmult];
  TGraphErrors *gv22_sub_nref_2[nmult];
  TGraphErrors *gv22_sub_nref_3[nmult];
  TF1 *fdir[nmult][npt];

  for (int im=0; im<nmult; im++){

    if ( !(im==nmult-1 || im==nref_1 || im==nref_2 || im==nref_3) ) continue;

    hv22_dir[im] = new TH1D(Form("hv22_dir_mult%d",im),"",npt,ptcut);
    gv22_sub_nref_1[im] = new TGraphErrors;
    gv22_sub_nref_2[im] = new TGraphErrors;
    gv22_sub_nref_3[im] = new TGraphErrors;

    c1[im] = new TCanvas(Form("c1_%d",im),Form("c1_%d",im),1.1*300*4,300*2);
    c1[im]->Divide(4,2);

    c2[im] = new TCanvas(Form("c2_%d",im),Form("c2_%d",im),1.1*200*4,200*2);
    c2[im]->Divide(4,2);

    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      c1[im]->cd(ipt);
			SetPadStyle(0.20,0.01,0.14,0.02);

      float max = TMath::Max(h1d_same_dphi_long[im][ipt][0]->GetMaximum(),h1d_same_dphi_short[im][ipt][0]->GetMaximum());
      float min = TMath::Min(h1d_same_dphi_long[im][ipt][0]->GetMinimum(),h1d_same_dphi_short[im][ipt][0]->GetMinimum());

      htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.05*(max-min), 3*const_pi/2, max+0.1*(max-min));
			SetHistoStyle("#Delta#phi (rad)","(1/N^{Trig}) dN^{Pair}/(d#Delta#phi)","",0.08);
			htmp->GetXaxis()->CenterTitle();
			htmp->GetXaxis()->SetTitleOffset(0.9);
			htmp->GetYaxis()->CenterTitle();
			htmp->GetYaxis()->SetNdivisions(5,3,0);
			htmp->GetYaxis()->SetTitleOffset(1.15);
      h1d_same_dphi_long[im][ipt][0]->SetMarkerStyle(24);
      h1d_same_dphi_long[im][ipt][0]->Draw("p same");

      h1d_same_dphi_short[im][ipt][0]->SetMarkerStyle(24);
      h1d_same_dphi_short[im][ipt][0]->SetMarkerColor(2);
      h1d_same_dphi_short[im][ipt][0]->Draw("p same");

      fdir[im][ipt] = new TF1(Form("fdir_%d_%d",im,ipt),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
      fdir[im][ipt]->SetLineWidth(3);
      fdir[im][ipt]->SetLineStyle(2);
      fdir[im][ipt]->SetLineColor(1);

      h1d_same_dphi_long[im][ipt][0]->Fit(fdir[im][ipt], "R0Q");
      fdir[im][ipt]->Draw("same");

      hv22_dir[im]->SetBinContent(ipt+1, fdir[im][ipt]->GetParameter(2));
      hv22_dir[im]->SetBinError(ipt+1, fdir[im][ipt]->GetParError(2));
      Nassoc[im][ipt] = fdir[im][ipt]->GetParameter(0)*2*const_pi;

			{
				TLatex *tex = new TLatex(4.4,max-0.2*(max-min),Form("%g<p_{T}^{Trig}<%g GeV/c",ptcut[ipt],ptcut[ipt+1]));
				tex->SetTextFont(62);
				tex->SetTextAlign(32);
				tex->SetTextSize(0.07);
				tex->Draw();
			}

			{
				if ( im==0 ){
					TLatex *tex = new TLatex(4.4,max-0.05*(max-min),Form("(%c)",97+ipt-1));
					tex->SetTextFont(62);
					tex->SetTextAlign(32);
					tex->SetTextSize(0.065);
					tex->Draw();
				}else{
					TLatex *tex = new TLatex(4.4,max-0.05*(max-min),Form("(%c)",97+ipt+7-1));
					tex->SetTextFont(62);
					tex->SetTextAlign(32);
					tex->SetTextSize(0.065);
					tex->Draw();
				}
			}

			c2[im]->cd(ipt);
      SetPadStyle();

      max = TMath::Max(h1d_same_dphi_sub[im][ipt]->GetMaximum(),h1d_same_dphi_sub[im][ipt]->GetMaximum());
      min = h1d_same_dphi_sub[im][ipt]->GetMinimum();

      htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, -1.1*fabs(min), 3*const_pi/2, 1.1*max);
      SetHistoStyle("","");

      h1d_same_dphi_sub[im][ipt]->SetMarkerStyle(24);
      h1d_same_dphi_sub[im][ipt]->Draw("p same");


    }//ipt

		c1[im]->cd(8);
		TLegend *leg = new TLegend(0.10,0.10,0.8,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		le = leg->AddEntry("",Form("h^{#pm}-h^{#pm}"),"");
		le->SetTextFont(62);
		le->SetTextSize(0.07);
		le = leg->AddEntry("",Form("PYTHIA8 %s 13 TeV",system),"");
		le->SetTextFont(62);
		le->SetTextSize(0.07);
		if ( im==nmult-1 ){
			le = leg->AddEntry("",Form("N^{ch}#geq%d",int(multcut_mid[im])),"");
		}else{
			le = leg->AddEntry("",Form("%d#leqN^{ch}<%d",int(multcut_mid[im]),int(multcut_mid[im+1])),"");
		}
		le->SetTextFont(62);
		le->SetTextSize(0.07);
		le = leg->AddEntry("",Form("0.5<p_{T}^{Assoc}<5.0 GeV/c"),"");
		le->SetTextSize(0.07);
		le = leg->AddEntry("",Form("0<|#eta|^{Trig}, |#eta|^{Assoc}<2.5"),"");
		le->SetTextSize(0.07);
		le = leg->AddEntry(h1d_same_dphi_long[im][1][0],"Long-range 2<|#Delta#eta|<5","P");
		le->SetTextSize(0.07);
		le = leg->AddEntry(h1d_same_dphi_short[im][1][0],"Short-range 0<|#Delta#eta|<1","P");
		le->SetTextSize(0.07);
		le = leg->AddEntry(fdir[im][1],"Fourier fits","L");
		le->SetTextSize(0.07);
		leg->Draw();

  }//

	//return;

	for (int im=nref_1+1; im<nmult; im++){
		if ( !(im==nmult-1) ) continue;
    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      float R_assoc = Nassoc[nref_1][ipt]/Nassoc[im][ipt];

      float Y_peri = h1d_same_dphi_sub[nref_1][ipt]->Integral(h1d_same_dphi_sub[nref_1][ipt]->FindBin(-1.2),h1d_same_dphi_sub[nref_1][ipt]->FindBin(+1.2)); 
      float Y_cent = h1d_same_dphi_sub[im][ipt]->Integral(h1d_same_dphi_sub[im][ipt]->FindBin(-1.2),h1d_same_dphi_sub[im][ipt]->FindBin(+1.2)); 

      float v22_sub = hv22_dir[im]->GetBinContent(ipt+1) - hv22_dir[nref_1]->GetBinContent(ipt+1)*R_assoc*Y_cent/Y_peri;

      gv22_sub_nref_1[im]->SetPoint(ipt-1, (ptcut[ipt]+ptcut[ipt+1])/2, v22_sub);
      gv22_sub_nref_1[im]->SetPointError(ipt-1, 0, 0);

      cout << "mult: " << im 
        << ", pt: " << ipt
        << ", R_assoc: " << R_assoc 
        << ", R_jet: " << Y_cent/Y_peri 
        << ", v22 cent: " << hv22_dir[im]->GetBinContent(ipt+1) 
        << ", v22 peri: " << hv22_dir[nref_1]->GetBinContent(ipt+1) 
        << ", v22_sub: " << v22_sub
        << endl;
    }
	}

	for (int im=nref_2+1; im<nmult; im++){
		if ( !(im==nmult-1) ) continue;
    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      float R_assoc = Nassoc[nref_2][ipt]/Nassoc[im][ipt];

      float Y_peri = h1d_same_dphi_sub[nref_2][ipt]->Integral(h1d_same_dphi_sub[nref_2][ipt]->FindBin(-1.2),h1d_same_dphi_sub[nref_2][ipt]->FindBin(+1.2)); 
      float Y_cent = h1d_same_dphi_sub[im][ipt]->Integral(h1d_same_dphi_sub[im][ipt]->FindBin(-1.2),h1d_same_dphi_sub[im][ipt]->FindBin(+1.2)); 

      float v22_sub = hv22_dir[im]->GetBinContent(ipt+1) - hv22_dir[nref_2]->GetBinContent(ipt+1)*R_assoc*Y_cent/Y_peri;

      gv22_sub_nref_2[im]->SetPoint(ipt-1, (ptcut[ipt]+ptcut[ipt+1])/2, v22_sub);
      gv22_sub_nref_2[im]->SetPointError(ipt-1, 0, 0);

      cout << "mult: " << im 
        << ", pt: " << ipt
        << ", R_assoc: " << R_assoc 
        << ", R_jet: " << Y_cent/Y_peri 
        << ", v22 cent: " << hv22_dir[im]->GetBinContent(ipt+1) 
        << ", v22 peri: " << hv22_dir[nref_2]->GetBinContent(ipt+1) 
        << ", v22_sub: " << v22_sub
        << endl;
    }
	}

	for (int im=nref_3+1; im<nmult; im++){
		if ( !(im==nmult-1) ) continue;
    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      float R_assoc = Nassoc[nref_3][ipt]/Nassoc[im][ipt];

      float Y_peri = h1d_same_dphi_sub[nref_3][ipt]->Integral(h1d_same_dphi_sub[nref_3][ipt]->FindBin(-1.2),h1d_same_dphi_sub[nref_3][ipt]->FindBin(+1.2)); 
      float Y_cent = h1d_same_dphi_sub[im][ipt]->Integral(h1d_same_dphi_sub[im][ipt]->FindBin(-1.2),h1d_same_dphi_sub[im][ipt]->FindBin(+1.2)); 

      float v22_sub = hv22_dir[im]->GetBinContent(ipt+1) - hv22_dir[nref_3]->GetBinContent(ipt+1)*R_assoc*Y_cent/Y_peri;

      gv22_sub_nref_3[im]->SetPoint(ipt-1, (ptcut[ipt]+ptcut[ipt+1])/2, v22_sub);
      gv22_sub_nref_3[im]->SetPointError(ipt-1, 0, 0);

      cout << "mult: " << im 
        << ", pt: " << ipt
        << ", R_assoc: " << R_assoc 
        << ", R_jet: " << Y_cent/Y_peri 
        << ", v22 cent: " << hv22_dir[im]->GetBinContent(ipt+1) 
        << ", v22 peri: " << hv22_dir[nref_3]->GetBinContent(ipt+1) 
        << ", v22_sub: " << v22_sub
        << endl;
    }
	}

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	SetPadStyle();
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.16);

	htmp = (TH1D*)gPad->DrawFrame(0,-0.002,5,0.015);
	SetHistoStyle("p_{T}^{Trig} (GeV/c)","v_{22}","",0.065);

	TLine *line = new TLine(0,0,5,0);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->SetLineColor(1);
	line->Draw();

	hv22_dir[nmult-1]->SetMarkerStyle(24);
	hv22_dir[nmult-1]->SetMarkerColor(1);
	hv22_dir[nmult-1]->SetLineColor(1);
	hv22_dir[nmult-1]->SetLineWidth(2);
	hv22_dir[nmult-1]->Draw("p same");

	hv22_dir[nref_1]->SetMarkerStyle(24);
	hv22_dir[nref_1]->SetMarkerColor(1);
	hv22_dir[nref_1]->SetLineColor(1);
	//hv22_dir[nref_1]->Draw("p same");

	hv22_dir[nref_2]->SetMarkerStyle(24);
	hv22_dir[nref_2]->SetMarkerColor(2);
	hv22_dir[nref_2]->SetLineColor(2);
	//hv22_dir[nref_2]->Draw("p same");

	hv22_dir[nref_3]->SetMarkerStyle(24);
	hv22_dir[nref_3]->SetMarkerColor(4);
	hv22_dir[nref_3]->SetLineColor(4);
	//hv22_dir[nref_3]->Draw("p same");

	gv22_sub_nref_1[nmult-1]->SetLineWidth(4);
	gv22_sub_nref_1[nmult-1]->SetLineColor(1);
	gv22_sub_nref_1[nmult-1]->Draw("c same");

	gv22_sub_nref_2[nmult-1]->SetLineWidth(4);
	gv22_sub_nref_2[nmult-1]->SetLineColor(2);
	gv22_sub_nref_2[nmult-1]->Draw("c same");

	gv22_sub_nref_3[nmult-1]->SetLineWidth(4);
	gv22_sub_nref_3[nmult-1]->SetLineColor(4);
	gv22_sub_nref_3[nmult-1]->Draw("c same");

	TLegend *leg = new TLegend(0.18,0.48,0.55,0.93);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	le = leg->AddEntry("","PYTHIA8 p+p 13 TeV","");
  le->SetTextFont(62);
  le->SetTextSize(0.05);
	le = leg->AddEntry("","Method 2","");
  le->SetTextFont(32);
  le->SetTextSize(0.05);
	le = leg->AddEntry("","0.5<p_{T}^{Assoc}<5.0 GeV/c","");
	le->SetTextSize(0.045);
	le = leg->AddEntry("","0<|#eta^{Trig}|, |#eta^{Assoc}|<2.5","");
	le->SetTextSize(0.045);
	le = leg->AddEntry(hv22_dir[nmult-1],"Fourier fits","P");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gv22_sub_nref_1[nmult-1],Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_1*10,nref_1*10+10),"L");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gv22_sub_nref_2[nmult-1],Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_2*10,nref_2*10+10),"L");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gv22_sub_nref_3[nmult-1],Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_3*10,nref_3*10+10),"L");
	le->SetTextSize(0.045);
	leg->Draw();

	{
		TLatex *tex = new TLatex(4,0.01,Form("N^{ch}#geq85"));
		tex->SetTextFont(62);
		tex->SetTextSize(0.05);
		tex->Draw();
	}

	{
		TLatex *tex = new TLatex(4.5,0.013,"(b)");
		tex->SetTextSize(0.045);
		tex->SetTextFont(62);
		tex->Draw();
	}

  //return;

	/*
  if ( bSAVE ){
    c1->cd();
    c1->SaveAs(Form("%s_1d.gif",dataset));

    c3->cd();
    c3->SaveAs(Form("%s_v22.gif",dataset));
  }
	*/

}

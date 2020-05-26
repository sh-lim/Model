#include "Style.h"

#define NHAR 3 
#define PI 3.1415927

TH1D *g_h_central;
TH1D *g_h_peripheral;
TH1D *g_h_template_fit;
double g_ped;

Double_t fitf(Double_t *x,Double_t *par){//par[0]=scale;par[1--NHAR]=vnn
	static TF1 *constant = new TF1("constant","[0]",-TMath::Pi()/2,1.5*PI);
	static TF1 *vnn[NHAR]={};

	if(!vnn[0]){
		std::cout<<"Fitf Initializing functions"<<std::endl;
		char name [100];
		char title[100];
		for(int i=0;i<NHAR;i++){
			sprintf(name ,"v%d%d"                ,i+2,i+2);
			sprintf(title,"[0]*(2*[1]*cos(%d*x))",i+2);
			vnn[i]= new TF1(name,title,-PI/2,1.5*PI);
		}
	}

	if(x[0]<-1.4){
		g_h_template_fit->Reset();

		g_h_template_fit->Add(g_h_peripheral);
		g_h_template_fit->Scale(par[0]);
		g_ped = (g_h_central->Integral()-g_h_template_fit->Integral())/g_h_peripheral->GetNbinsX();

		constant->SetParameter(0,g_ped);
		g_h_template_fit->Add(constant);
		for(int i=0;i<NHAR;i++){
			vnn[i]->SetParameters(g_ped,par[i+1]);
			g_h_template_fit->Add(vnn[i]);
		}
	}

	double returnval = g_h_template_fit->GetBinContent(g_h_template_fit->FindBin(x[0]));
	return returnval;
}

void MyChi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag){
	double chisq = 0;
	double delta;
	for (int i=1; i<=g_h_peripheral->GetNbinsX(); i++) {
		Double_t xvalue = g_h_central->GetBinCenter(i);
		delta  = g_h_central->GetBinContent(i)-fitf(&xvalue,par);

		chisq += delta*delta/
			(g_h_central->GetBinError(i)*g_h_central->GetBinError(i)
			 +par[0]*par[0]*g_h_peripheral->GetBinError(i)*g_h_peripheral->GetBinError(i));
	}
	f=chisq;
}

void Draw_atlas_pT(const char *dataset="pp13TeV_nondiffractive_grp0_try4", const bool bSAVE=false){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TLegendEntry *le;

  const char *system = "p+p";

	TFile *infile = new TFile(Form("outfile_hist_%s.root",dataset),"read");

	const int nmult = 5;
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
      h1d_same_dphi_long[im][ipt][0] = (TH1D*)h2d_same_dphi_deta[im][ipt]->ProjectionX(Form("h1d_same_dphi_long_%d_%d_0",ipt,im),etamin,etamax);

      etamin = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(+2.0+0.0001);
      etamax = h2d_same_dphi_deta[im][ipt]->GetYaxis()->FindBin(+5.0-0.0001);
      h1d_same_dphi_long[im][ipt][1] = (TH1D*)h2d_same_dphi_deta[im][ipt]->ProjectionX(Form("h1d_same_dphi_long_%d_%d_1",ipt,im),etamin,etamax);

      h1d_same_dphi_long[im][ipt][0]->Add(h1d_same_dphi_long[im][ipt][1]);
      h1d_same_dphi_long[im][ipt][0]->Scale(1./6.0/h1d_same_dphi_long[im][ipt][0]->GetBinWidth(1));
    }
	}//

	//return;

	TH1D *hv22_dir_1[nmult];
	TH1D *hv22_dir_2[nmult];
	TH1D *hv22_dir_3[nmult];

  TGraphErrors *gv22_sub_1[nmult];
  TGraphErrors *gv22_sub_2[nmult];
  TGraphErrors *gv22_sub_3[nmult];

	TF1 *fdir_1[nmult][npt];
	TF1 *fdir_2[nmult][npt];
	TF1 *fdir_3[nmult][npt];

  TF1 *func_1[nmult][npt];
  TF1 *func_2[nmult][npt];
  TF1 *func_3[nmult][npt];

  TF1 *f_ped_1[nmult][npt];
  TF1 *f_ped_2[nmult][npt];
  TF1 *f_ped_3[nmult][npt];

  TF1 *f_combined_1[nmult][npt];
  TF1 *f_combined_2[nmult][npt];
  TF1 *f_combined_3[nmult][npt];

  TH1D *h_rescaledperipheral_1[nmult][npt];
  TH1D *h_rescaledperipheral_2[nmult][npt];
  TH1D *h_rescaledperipheral_3[nmult][npt];

  TCanvas *c1_1[nmult];
	for (int im=nref_1; im<nmult; im++){
		if ( !(im==nmult-1) ) continue;

    hv22_dir_1[im] = new TH1D(Form("hv22_dir_1_mult%d",im),"",npt,ptcut);
    gv22_sub_1[im] = new TGraphErrors;

    c1_1[im] = new TCanvas(Form("c1_1_%d",im),Form("c1_1_%d",im),1.1*200*3,200*3);
    c1_1[im]->Divide(3,3);

    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      c1_1[im]->cd(ipt);
      SetPadStyle();

      float max = h1d_same_dphi_long[im][ipt][0]->GetMaximum();
      float min = h1d_same_dphi_long[im][ipt][0]->GetMinimum();

      htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.1*(max-min), 3*const_pi/2, max+0.1*(max-min));
      SetHistoStyle("#Delta#phi (rad)","Per trigger yield");
      h1d_same_dphi_long[im][ipt][0]->SetMarkerStyle(24);
      h1d_same_dphi_long[im][ipt][0]->Draw("p same");

      fdir_1[im][ipt] = new TF1(Form("fdir_1_%d_%d",im,ipt),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
      fdir_1[im][ipt]->SetLineWidth(3);
      fdir_1[im][ipt]->SetLineStyle(2);
      fdir_1[im][ipt]->SetLineColor(1);

      h1d_same_dphi_long[im][ipt][0]->Fit(fdir_1[im][ipt], "R0Q");
      fdir_1[im][ipt]->Draw("same");

      hv22_dir_1[im]->SetBinContent(ipt+1, fdir_1[im][ipt]->GetParameter(2));
      hv22_dir_1[im]->SetBinError(ipt+1, fdir_1[im][ipt]->GetParError(2));

      if ( im>nref_1 ){
        g_h_central = h1d_same_dphi_long[im][ipt][0];
        g_h_peripheral = h1d_same_dphi_long[nref_1][ipt][0];
        g_h_template_fit = (TH1D*)g_h_peripheral->Clone(Form("h_template_fit_1_for_mult%d_pt%d",im,ipt));
        g_h_template_fit->Reset();

        func_1[im][ipt] = new TF1("f_fit",fitf,-const_pi/2,1.5*const_pi,NHAR+1);
        func_1[im][ipt]->SetParameter(0,1);
        for (int i=1; i<NHAR+1; i++) func_1[im][ipt]->SetParameter(i,0.0);

        TH1D* l_h_central_ = (TH1D*)g_h_central->Clone("l_h_central_");
        TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);

        l_h_central_->Fit(func_1[im][ipt],"QU0");

        TF1 *f_pedestal = new TF1("f_pedestal","[0]",-0.5*const_pi,1.5*const_pi);
        f_pedestal->SetLineColor(kOrange-3);
        f_pedestal->SetLineStyle(2);
        f_pedestal->SetLineWidth(3);
        f_pedestal->SetParameter(0,g_ped);
        f_pedestal->SetParError(0,func_1[im][ipt]->GetParError(0)*g_h_peripheral->Integral()/g_h_peripheral->GetNbinsX());
        //f_pedestal->Draw("same");
        //
        //cout << "f_pedestal par: " << f_pedestal->GetParameter(0) << ", err: " << f_pedestal->GetParError(0) << endl;
        //cout << "f_scale par: " << func_1[im][ipt]->GetParameter(0) << ", err: " << func_1[im][ipt]->GetParError(0) << endl;

        h_rescaledperipheral_1[im][ipt] = (TH1D*)h1d_same_dphi_long[nref_1][ipt][0]->Clone(Form("h_rescaledperipheral_1_mult%d_pt%d",im,ipt));
        h_rescaledperipheral_1[im][ipt]->Scale(func_1[im][ipt]->GetParameter(0));
        h_rescaledperipheral_1[im][ipt]->Add(f_pedestal);
        h_rescaledperipheral_1[im][ipt]->SetMarkerStyle(31);
        h_rescaledperipheral_1[im][ipt]->SetLineWidth(3);

        f_ped_1[im][ipt] = new TF1("f_ped_1","[0]",-0.5*const_pi,1.5*const_pi);
        //f_ped_1[im][ipt]->SetParameter(0,h_rescaledperipheral_1[im][ipt]->GetBinContent(h_rescaledperipheral_1[im][ipt]->FindBin(0)));
        f_ped_1[im][ipt]->SetParameter(0,h_rescaledperipheral_1[im][ipt]->GetMinimum());
        f_ped_1[im][ipt]->SetLineColor(kOrange-3);
        f_ped_1[im][ipt]->SetLineStyle(2);
        f_ped_1[im][ipt]->SetLineWidth(3);

        f_combined_1[im][ipt] = new TF1("f_combined_1","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x)) + [4]",-0.5*const_pi,1.5*const_pi);
        f_combined_1[im][ipt]->SetParameter(0, g_ped);
        f_combined_1[im][ipt]->SetParameter(1, func_1[im][ipt]->GetParameter(1));
        f_combined_1[im][ipt]->SetParError(1, func_1[im][ipt]->GetParError(1));
        f_combined_1[im][ipt]->SetParameter(2, func_1[im][ipt]->GetParameter(2));
        f_combined_1[im][ipt]->SetParameter(3, func_1[im][ipt]->GetParameter(3));
        f_combined_1[im][ipt]->SetParameter(4, h_rescaledperipheral_1[im][ipt]->GetMinimum() - g_ped);
        f_combined_1[im][ipt]->SetLineStyle(2);
        f_combined_1[im][ipt]->SetLineColor(4);
        f_combined_1[im][ipt]->SetLineWidth(3);

        func_1[im][ipt]->Draw("same");
        f_ped_1[im][ipt]->Draw("same");
        h_rescaledperipheral_1[im][ipt]->Draw("p same");
        f_combined_1[im][ipt]->Draw("same");

        float v22 = f_combined_1[im][ipt]->GetParameter(1);
        float v22_err = f_combined_1[im][ipt]->GetParError(1);

        gv22_sub_1[im]->SetPoint(ipt-1, (ptcut[ipt]+ptcut[ipt+1])/2, v22);
        gv22_sub_1[im]->SetPointError(ipt-1, 0, 0);
      }//

			/*
      {
        TLegend *leg = new TLegend(0.15,0.73,0.7,0.93);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry("",Form("%s 200 GeV",system),"");
        leg->AddEntry("",Form("Centrality %d-%d%c",cent_fwd[im+1],cent_fwd[im],'%'),"");
        leg->AddEntry("",Form("%g<p_{T}^{Trig}<%g GeV/c, 0.2<p_{T}^{Assoc}<3.0 GeV/c",ptcut[ipt],ptcut[ipt+1]),"");
        leg->AddEntry(h1d_same_dphi_long[im][ipt][0],"1.0<|#Delta#eta|<1.8","P");
        leg->Draw();
      }
			*/

    }//ipt

		c1_1[im]->cd();
		c1_1[im]->Update();
  }//

  TCanvas *c1_2[nmult];
	for (int im=nref_2; im<nmult; im++){
		if ( !(im==nmult-1) ) continue;

    hv22_dir_2[im] = new TH1D(Form("hv22_dir_2_mult%d",im),"",npt,ptcut);
    gv22_sub_2[im] = new TGraphErrors;

    c1_2[im] = new TCanvas(Form("c1_2_%d",im),Form("c1_2_%d",im),1.1*200*3,200*3);
    c1_2[im]->Divide(3,3);

    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      c1_2[im]->cd(ipt);
      SetPadStyle();

      float max = h1d_same_dphi_long[im][ipt][0]->GetMaximum();
      float min = h1d_same_dphi_long[im][ipt][0]->GetMinimum();

      htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.1*(max-min), 3*const_pi/2, max+0.1*(max-min));
      SetHistoStyle("#Delta#phi (rad)","Per trigger yield");
      h1d_same_dphi_long[im][ipt][0]->SetMarkerStyle(24);
      h1d_same_dphi_long[im][ipt][0]->Draw("p same");

      fdir_2[im][ipt] = new TF1(Form("fdir_2_%d_%d",im,ipt),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
      fdir_2[im][ipt]->SetLineWidth(3);
      fdir_2[im][ipt]->SetLineStyle(2);
      fdir_2[im][ipt]->SetLineColor(1);

      h1d_same_dphi_long[im][ipt][0]->Fit(fdir_2[im][ipt], "R0Q");
      fdir_2[im][ipt]->Draw("same");

      hv22_dir_2[im]->SetBinContent(ipt+1, fdir_2[im][ipt]->GetParameter(2));
      hv22_dir_2[im]->SetBinError(ipt+1, fdir_2[im][ipt]->GetParError(2));

      if ( im>nref_2 ){
        g_h_central = h1d_same_dphi_long[im][ipt][0];
        g_h_peripheral = h1d_same_dphi_long[nref_2][ipt][0];
        g_h_template_fit = (TH1D*)g_h_peripheral->Clone(Form("h_template_fit_2_for_mult%d_pt%d",im,ipt));
        g_h_template_fit->Reset();

        func_2[im][ipt] = new TF1("f_fit",fitf,-const_pi/2,1.5*const_pi,NHAR+1);
        func_2[im][ipt]->SetParameter(0,1);
        for (int i=1; i<NHAR+1; i++) func_2[im][ipt]->SetParameter(i,0.0);

        TH1D* l_h_central_ = (TH1D*)g_h_central->Clone("l_h_central_");
        TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);

        l_h_central_->Fit(func_2[im][ipt],"QU0");

        TF1 *f_pedestal = new TF1("f_pedestal","[0]",-0.5*const_pi,1.5*const_pi);
        f_pedestal->SetLineColor(kOrange-3);
        f_pedestal->SetLineStyle(2);
        f_pedestal->SetLineWidth(3);
        f_pedestal->SetParameter(0,g_ped);
        f_pedestal->SetParError(0,func_2[im][ipt]->GetParError(0)*g_h_peripheral->Integral()/g_h_peripheral->GetNbinsX());
        //f_pedestal->Draw("same");
        //
        //cout << "f_pedestal par: " << f_pedestal->GetParameter(0) << ", err: " << f_pedestal->GetParError(0) << endl;
        //cout << "f_scale par: " << func_2[im][ipt]->GetParameter(0) << ", err: " << func_2[im][ipt]->GetParError(0) << endl;

        h_rescaledperipheral_2[im][ipt] = (TH1D*)h1d_same_dphi_long[nref_2][ipt][0]->Clone(Form("h_rescaledperipheral_2_mult%d_pt%d",im,ipt));
        h_rescaledperipheral_2[im][ipt]->Scale(func_2[im][ipt]->GetParameter(0));
        h_rescaledperipheral_2[im][ipt]->Add(f_pedestal);
        h_rescaledperipheral_2[im][ipt]->SetMarkerStyle(31);
        h_rescaledperipheral_2[im][ipt]->SetLineWidth(3);

        f_ped_2[im][ipt] = new TF1("f_ped_2","[0]",-0.5*const_pi,1.5*const_pi);
        //f_ped_2[im][ipt]->SetParameter(0,h_rescaledperipheral_2[im][ipt]->GetBinContent(h_rescaledperipheral_2[im][ipt]->FindBin(0)));
        f_ped_2[im][ipt]->SetParameter(0,h_rescaledperipheral_2[im][ipt]->GetMinimum());
        f_ped_2[im][ipt]->SetLineColor(kOrange-3);
        f_ped_2[im][ipt]->SetLineStyle(2);
        f_ped_2[im][ipt]->SetLineWidth(3);

        f_combined_2[im][ipt] = new TF1("f_combined_2","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x)) + [4]",-0.5*const_pi,1.5*const_pi);
        f_combined_2[im][ipt]->SetParameter(0, g_ped);
        f_combined_2[im][ipt]->SetParameter(1, func_2[im][ipt]->GetParameter(1));
        f_combined_2[im][ipt]->SetParError(1, func_2[im][ipt]->GetParError(1));
        f_combined_2[im][ipt]->SetParameter(2, func_2[im][ipt]->GetParameter(2));
        f_combined_2[im][ipt]->SetParameter(3, func_2[im][ipt]->GetParameter(3));
        f_combined_2[im][ipt]->SetParameter(4, h_rescaledperipheral_2[im][ipt]->GetMinimum() - g_ped);
        f_combined_2[im][ipt]->SetLineStyle(2);
        f_combined_2[im][ipt]->SetLineColor(4);
        f_combined_2[im][ipt]->SetLineWidth(3);


        func_2[im][ipt]->Draw("same");
        f_ped_2[im][ipt]->Draw("same");
        h_rescaledperipheral_2[im][ipt]->Draw("p same");
        f_combined_2[im][ipt]->Draw("same");

        float v22 = f_combined_2[im][ipt]->GetParameter(1);
        float v22_err = f_combined_2[im][ipt]->GetParError(1);

        gv22_sub_2[im]->SetPoint(ipt-1, (ptcut[ipt]+ptcut[ipt+1])/2, v22);
        gv22_sub_2[im]->SetPointError(ipt-1, 0, 0);
      }//

			/*
      {
        TLegend *leg = new TLegend(0.15,0.73,0.7,0.93);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry("",Form("%s 200 GeV",system),"");
        leg->AddEntry("",Form("Centrality %d-%d%c",cent_fwd[im+1],cent_fwd[im],'%'),"");
        leg->AddEntry("",Form("%g<p_{T}^{Trig}<%g GeV/c, 0.2<p_{T}^{Assoc}<3.0 GeV/c",ptcut[ipt],ptcut[ipt+1]),"");
        leg->AddEntry(h1d_same_dphi_long[im][ipt][0],"1.0<|#Delta#eta|<1.8","P");
        leg->Draw();
      }
			*/

    }//ipt

		c1_2[im]->cd();
		c1_2[im]->Update();
  }//

  TCanvas *c1_3[nmult];
	for (int im=nref_3; im<nmult; im++){
		if ( !(im==nmult-1) ) continue;

    hv22_dir_3[im] = new TH1D(Form("hv22_dir_3_mult%d",im),"",npt,ptcut);
    gv22_sub_3[im] = new TGraphErrors;

    c1_3[im] = new TCanvas(Form("c1_3_%d",im),Form("c1_3_%d",im),1.1*200*3,200*3);
    c1_3[im]->Divide(3,3);

    for (int ipt=0; ipt<npt; ipt++){

			if ( ipt==0 ) continue;

      c1_3[im]->cd(ipt);
      SetPadStyle();

      float max = h1d_same_dphi_long[im][ipt][0]->GetMaximum();
      float min = h1d_same_dphi_long[im][ipt][0]->GetMinimum();

      htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.1*(max-min), 3*const_pi/2, max+0.1*(max-min));
      SetHistoStyle("#Delta#phi (rad)","Per trigger yield");
      h1d_same_dphi_long[im][ipt][0]->SetMarkerStyle(24);
      h1d_same_dphi_long[im][ipt][0]->Draw("p same");

      fdir_3[im][ipt] = new TF1(Form("fdir_3_%d_%d",im,ipt),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
      fdir_3[im][ipt]->SetLineWidth(3);
      fdir_3[im][ipt]->SetLineStyle(2);
      fdir_3[im][ipt]->SetLineColor(1);

      h1d_same_dphi_long[im][ipt][0]->Fit(fdir_3[im][ipt], "R0Q");
      fdir_3[im][ipt]->Draw("same");

      hv22_dir_3[im]->SetBinContent(ipt+1, fdir_3[im][ipt]->GetParameter(2));
      hv22_dir_3[im]->SetBinError(ipt+1, fdir_3[im][ipt]->GetParError(2));

      if ( im>nref_3 ){
        g_h_central = h1d_same_dphi_long[im][ipt][0];
        g_h_peripheral = h1d_same_dphi_long[nref_3][ipt][0];
        g_h_template_fit = (TH1D*)g_h_peripheral->Clone(Form("h_template_fit_3_for_mult%d_pt%d",im,ipt));
        g_h_template_fit->Reset();

        func_3[im][ipt] = new TF1("f_fit",fitf,-const_pi/2,1.5*const_pi,NHAR+1);
        func_3[im][ipt]->SetParameter(0,1);
        for (int i=1; i<NHAR+1; i++) func_3[im][ipt]->SetParameter(i,0.0);

        TH1D* l_h_central_ = (TH1D*)g_h_central->Clone("l_h_central_");
        TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);

        l_h_central_->Fit(func_3[im][ipt],"QU0");

        TF1 *f_pedestal = new TF1("f_pedestal","[0]",-0.5*const_pi,1.5*const_pi);
        f_pedestal->SetLineColor(kOrange-3);
        f_pedestal->SetLineStyle(2);
        f_pedestal->SetLineWidth(3);
        f_pedestal->SetParameter(0,g_ped);
        f_pedestal->SetParError(0,func_3[im][ipt]->GetParError(0)*g_h_peripheral->Integral()/g_h_peripheral->GetNbinsX());
        //f_pedestal->Draw("same");
        //
        //cout << "f_pedestal par: " << f_pedestal->GetParameter(0) << ", err: " << f_pedestal->GetParError(0) << endl;
        //cout << "f_scale par: " << func_3[im][ipt]->GetParameter(0) << ", err: " << func_3[im][ipt]->GetParError(0) << endl;

        h_rescaledperipheral_3[im][ipt] = (TH1D*)h1d_same_dphi_long[nref_3][ipt][0]->Clone(Form("h_rescaledperipheral_3_mult%d_pt%d",im,ipt));
        h_rescaledperipheral_3[im][ipt]->Scale(func_3[im][ipt]->GetParameter(0));
        h_rescaledperipheral_3[im][ipt]->Add(f_pedestal);
        h_rescaledperipheral_3[im][ipt]->SetMarkerStyle(31);
        h_rescaledperipheral_3[im][ipt]->SetLineWidth(3);

        f_ped_3[im][ipt] = new TF1("f_ped_3","[0]",-0.5*const_pi,1.5*const_pi);
        //f_ped_3[im][ipt]->SetParameter(0,h_rescaledperipheral_3[im][ipt]->GetBinContent(h_rescaledperipheral_3[im][ipt]->FindBin(0)));
        f_ped_3[im][ipt]->SetParameter(0,h_rescaledperipheral_3[im][ipt]->GetMinimum());
        f_ped_3[im][ipt]->SetLineColor(kOrange-3);
        f_ped_3[im][ipt]->SetLineStyle(2);
        f_ped_3[im][ipt]->SetLineWidth(3);

        f_combined_3[im][ipt] = new TF1("f_combined_3","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x)) + [4]",-0.5*const_pi,1.5*const_pi);
        f_combined_3[im][ipt]->SetParameter(0, g_ped);
        f_combined_3[im][ipt]->SetParameter(1, func_3[im][ipt]->GetParameter(1));
        f_combined_3[im][ipt]->SetParError(1, func_3[im][ipt]->GetParError(1));
        f_combined_3[im][ipt]->SetParameter(2, func_3[im][ipt]->GetParameter(2));
        f_combined_3[im][ipt]->SetParameter(3, func_3[im][ipt]->GetParameter(3));
        f_combined_3[im][ipt]->SetParameter(4, h_rescaledperipheral_3[im][ipt]->GetMinimum() - g_ped);
        f_combined_3[im][ipt]->SetLineStyle(2);
        f_combined_3[im][ipt]->SetLineColor(4);
        f_combined_3[im][ipt]->SetLineWidth(3);


        func_3[im][ipt]->Draw("same");
        f_ped_3[im][ipt]->Draw("same");
        h_rescaledperipheral_3[im][ipt]->Draw("p same");
        f_combined_3[im][ipt]->Draw("same");

        float v22 = f_combined_3[im][ipt]->GetParameter(1);
        float v22_err = f_combined_3[im][ipt]->GetParError(1);

        gv22_sub_3[im]->SetPoint(ipt-1, (ptcut[ipt]+ptcut[ipt+1])/2, v22);
        gv22_sub_3[im]->SetPointError(ipt-1, 0, 0);
      }//

			/*
      {
        TLegend *leg = new TLegend(0.15,0.73,0.7,0.93);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry("",Form("%s 200 GeV",system),"");
        leg->AddEntry("",Form("Centrality %d-%d%c",cent_fwd[im+1],cent_fwd[im],'%'),"");
        leg->AddEntry("",Form("%g<p_{T}^{Trig}<%g GeV/c, 0.2<p_{T}^{Assoc}<3.0 GeV/c",ptcut[ipt],ptcut[ipt+1]),"");
        leg->AddEntry(h1d_same_dphi_long[im][ipt][0],"1.0<|#Delta#eta|<1.8","P");
        leg->Draw();
      }
			*/

    }//ipt

		c1_3[im]->cd();
		c1_3[im]->Update();
  }//

	//return;

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

	//hv22_dir_1[nref_1]->SetMarkerStyle(24);
	//hv22_dir_1[nref_1]->SetMarkerColor(1);
	//hv22_dir_1[nref_1]->SetLineColor(1);
	//hv22_dir_1[nref_1]->Draw("p same");

	hv22_dir_1[nmult-1]->SetMarkerStyle(24);
	hv22_dir_1[nmult-1]->SetMarkerColor(1);
	hv22_dir_1[nmult-1]->SetLineColor(1);
	hv22_dir_1[nmult-1]->SetLineWidth(2);
	hv22_dir_1[nmult-1]->Draw("p same");

	gv22_sub_1[nmult-1]->SetLineWidth(4);
	gv22_sub_1[nmult-1]->SetLineColor(1);
	gv22_sub_1[nmult-1]->Draw("c same");

	gv22_sub_2[nmult-1]->SetLineWidth(4);
	gv22_sub_2[nmult-1]->SetLineColor(2);
	gv22_sub_2[nmult-1]->Draw("c same");

	gv22_sub_3[nmult-1]->SetLineWidth(4);
	gv22_sub_3[nmult-1]->SetLineColor(4);
	gv22_sub_3[nmult-1]->Draw("c same");

	TLegend *leg = new TLegend(0.18,0.48,0.55,0.93);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	le = leg->AddEntry("","PYTHIA8 p+p 13 TeV","");
  le->SetTextFont(62);
  le->SetTextSize(0.05);
	le = leg->AddEntry("","Method 1","");
  le->SetTextFont(32);
  le->SetTextSize(0.05);
	le = leg->AddEntry("","0.5<p_{T}^{Assoc}<5.0 GeV/c","");
	le->SetTextSize(0.045);
	le = leg->AddEntry("","0<|#eta^{Trig}|, |#eta^{Assoc}|<2.5","");
	le->SetTextSize(0.045);
  le = leg->AddEntry(hv22_dir_1[nmult-1],"Fourier fits","P");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gv22_sub_1[nmult-1],Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_1*10,nref_1*10+10),"L");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gv22_sub_2[nmult-1],Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_2*10,nref_2*10+10),"L");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gv22_sub_3[nmult-1],Form("v_{22}^{sub}, LM %d#leqN^{ch}<%d",nref_3*10,nref_3*10+10),"L");
	le->SetTextSize(0.045);
	leg->Draw();

	{
		TLatex *tex = new TLatex(4,0.01,Form("N^{ch}#geq85"));
		tex->SetTextFont(62);
		tex->SetTextSize(0.05);
		tex->Draw();
	}

	{
		TLatex *tex = new TLatex(4.5,0.013,"(a)");
		tex->SetTextSize(0.045);
		tex->SetTextFont(62);
		tex->Draw();
	}


	/*
  if ( bSAVE ){
    //c1->cd();
    //c1->SaveAs(Form("%s_fwd_atlas_1d_ref1.gif",dataset));

    //c2->cd();
    //c2->SaveAs(Form("%s_atlas_1d_ref2.gif",dataset));

    c3->cd();
    c3->SaveAs(Form("%s_fwd_pT_atlas_v22.gif",dataset));
  }
	*/

}

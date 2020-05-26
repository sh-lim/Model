#include "Style.h"

#define NHAR 3 
#define PI 3.1415927

TH1D *g_h_central;
TH1D *g_h_peripheral;
TH1D *g_h_template_fit;
TF1 *g_f_template_fit;
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

void Draw_atlas(const char *dataset="pp13TeV_nondiffractive_grp0_try3", const bool bSAVE=false){

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

	}//

  //return;

	TH1D *hv22_dir_nref_1 = new TH1D("hv22_dir_nref_1","",20,0,200);
	TH1D *hv22_dir_nref_2 = new TH1D("hv22_dir_nref_2","",20,0,200);
	TH1D *hv22_dir_nref_3 = new TH1D("hv22_dir_nref_3","",20,0,200);

	TGraphErrors *gv22_sub_nref_1 = new TGraphErrors;
	TGraphErrors *gv22_sub_nref_2 = new TGraphErrors;
	TGraphErrors *gv22_sub_nref_3 = new TGraphErrors;

	TF1 *fdir_nref_1[nmult];
  TF1 *func_nref_1[nmult];
  TF1 *f_ped_nref_1[nmult];
  TF1 *f_combined_nref_1[nmult];
  TH1D *h_rescaledperipheral_nref_1[nmult];

	TF1 *fdir_nref_2[nmult];
  TF1 *func_nref_2[nmult];
  TF1 *f_ped_nref_2[nmult];
  TF1 *f_combined_nref_2[nmult];
  TH1D *h_rescaledperipheral_nref_2[nmult];

	TF1 *fdir_nref_3[nmult];
  TF1 *func_nref_3[nmult];
  TF1 *f_ped_nref_3[nmult];
  TF1 *f_combined_nref_3[nmult];
  TH1D *h_rescaledperipheral_nref_3[nmult];

  int g_count = 0;

	TCanvas *c1_1 = new TCanvas("c1_1","c1_1",1.1*200*5,200*3);
	c1_1->Divide(5,3);

	for (int im=0; im<nmult; im++){
		c1_1->cd(im+1);
		SetPadStyle();

		float max = h1d_same_dphi_long[im][0]->GetMaximum();
		float min = h1d_same_dphi_long[im][0]->GetMinimum();

    htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.05*(max-min), 3*const_pi/2, max+0.05*(max-min));
		SetHistoStyle("#Delta#phi (rad)","Per trigger yield");
		h1d_same_dphi_long[im][0]->SetMarkerStyle(24);
		h1d_same_dphi_long[im][0]->Draw("p same");

		fdir_nref_1[im] = new TF1(Form("fdir_%d",im),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
		fdir_nref_1[im]->SetLineWidth(3);
		fdir_nref_1[im]->SetLineStyle(2);
		fdir_nref_1[im]->SetLineColor(1);

		h1d_same_dphi_long[im][0]->Fit(fdir_nref_1[im], "R0Q");
		fdir_nref_1[im]->Draw("same");

		hv22_dir_nref_1->SetBinContent(im+1, fdir_nref_1[im]->GetParameter(2));
		hv22_dir_nref_1->SetBinError(im+1, fdir_nref_1[im]->GetParError(2));

    if ( im>nref_1+1 ){

      g_h_central = h1d_same_dphi_long[im][0];
      g_h_peripheral = h1d_same_dphi_long[nref_1][0];
			g_h_template_fit = (TH1D*)g_h_peripheral->Clone(Form("h_template_fit_for_nref_1_mult%d",im));
			g_h_template_fit->Reset();

			func_nref_1[im] = new TF1("f_fit",fitf,-const_pi/2,1.5*const_pi,NHAR+1);
			func_nref_1[im]->SetParameter(0,1);
			for (int i=1; i<NHAR+1; i++) func_nref_1[im]->SetParameter(i,0.0);

			TH1D* l_h_central_ = (TH1D*)g_h_central->Clone("l_h_central_");
			TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);

			l_h_central_->Fit(func_nref_1[im],"QU0");

			TF1 *f_pedestal = new TF1("f_pedestal","[0]",-0.5*const_pi,1.5*const_pi);
			f_pedestal->SetLineColor(kOrange-3);
			f_pedestal->SetLineStyle(2);
			f_pedestal->SetLineWidth(3);
			f_pedestal->SetParameter(0,g_ped);
			f_pedestal->SetParError(0,func_nref_1[im]->GetParError(0)*g_h_peripheral->Integral()/g_h_peripheral->GetNbinsX());
      //f_pedestal->Draw("same");
      //
			//cout << "f_pedestal par: " << f_pedestal->GetParameter(0) << ", err: " << f_pedestal->GetParError(0) << endl;
			//cout << "f_scale par: " << func_nref_1[im]->GetParameter(0) << ", err: " << func_nref_1[im]->GetParError(0) << endl;

			h_rescaledperipheral_nref_1[im] = (TH1D*)h1d_same_dphi_long[nref_1][0]->Clone(Form("h_rescaledperipheral_nref_1_mult%d",im));
			h_rescaledperipheral_nref_1[im]->Scale(func_nref_1[im]->GetParameter(0));
			h_rescaledperipheral_nref_1[im]->Add(f_pedestal);
			h_rescaledperipheral_nref_1[im]->SetMarkerStyle(31);
			h_rescaledperipheral_nref_1[im]->SetLineWidth(3);

			f_ped_nref_1[im] = new TF1("f_ped","[0]",-0.5*const_pi,1.5*const_pi);
			//f_ped_nref_1[im]->SetParameter(0,h_rescaledperipheral_nref_1[im]->GetBinContent(h_rescaledperipheral_nref_1[im]->FindBin(0)));
			f_ped_nref_1[im]->SetParameter(0,h_rescaledperipheral_nref_1[im]->GetMinimum());
			f_ped_nref_1[im]->SetLineColor(kOrange-3);
			f_ped_nref_1[im]->SetLineStyle(2);
			f_ped_nref_1[im]->SetLineWidth(3);

			f_combined_nref_1[im] = new TF1("f_combined","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x)) + [4]",-0.5*const_pi,1.5*const_pi);
			f_combined_nref_1[im]->SetParameter(0, g_ped);
			f_combined_nref_1[im]->SetParameter(1, func_nref_1[im]->GetParameter(1));
			f_combined_nref_1[im]->SetParError(1, func_nref_1[im]->GetParError(1));
			f_combined_nref_1[im]->SetParameter(2, func_nref_1[im]->GetParameter(2));
			f_combined_nref_1[im]->SetParameter(3, func_nref_1[im]->GetParameter(3));
			//f_combined_nref_1[im]->SetParameter(4, h_rescaledperipheral_nref_1[im]->GetBinContent(h_rescaledperipheral_nref_1[im]->FindBin(0)) - g_ped);
			f_combined_nref_1[im]->SetParameter(4, h_rescaledperipheral_nref_1[im]->GetMinimum() - g_ped);
			f_combined_nref_1[im]->SetLineStyle(2);
			f_combined_nref_1[im]->SetLineColor(4);
			f_combined_nref_1[im]->SetLineWidth(3);


      func_nref_1[im]->Draw("same");
			f_ped_nref_1[im]->Draw("same");
      h_rescaledperipheral_nref_1[im]->Draw("p same");
			f_combined_nref_1[im]->Draw("same");

			float v22 = f_combined_nref_1[im]->GetParameter(1);
			float v22_err = f_combined_nref_1[im]->GetParError(1);

      gv22_sub_nref_1->SetPoint(g_count, 10*im+5, v22);
      gv22_sub_nref_1->SetPointError(g_count, 0, 0);

      cout << "ref: " << nref_1 << ", mult: " << im << ", g_ped: " << g_ped << ", v22_sub: " << v22 << endl;

      g_count++;

    }

    {
      TLegend *leg = new TLegend(0.15,0.73,0.7,0.93);
      leg->SetFillStyle(0);
			leg->SetBorderSize(0);
      leg->AddEntry("",Form("%s 13 TeV",system),"");
      leg->AddEntry("",Form("%d#leqN_{TRK}^{CH}<%d",10*im,10*im+10),"");
      leg->AddEntry(h1d_same_dphi_long[im][0],"2<|#Delta#eta|<5","P");
      leg->Draw();
    }
	}//im

  c1_1->cd();
  c1_1->Update();

  //return;

  g_count = 0;
	TCanvas *c1_2 = new TCanvas("c1_2","c1_2",1.1*200*5,200*3);
	c1_2->Divide(5,3);

	for (int im=0; im<nmult; im++){
		c1_2->cd(im+1);
		SetPadStyle();

		float max = h1d_same_dphi_long[im][0]->GetMaximum();
		float min = h1d_same_dphi_long[im][0]->GetMinimum();

    htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.05*(max-min), 3*const_pi/2, max+0.05*(max-min));
		SetHistoStyle("#Delta#phi (rad)","Per trigger yield");
		h1d_same_dphi_long[im][0]->SetMarkerStyle(24);
		h1d_same_dphi_long[im][0]->Draw("p same");

		fdir_nref_2[im] = new TF1(Form("fdir_%d",im),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
		fdir_nref_2[im]->SetLineWidth(3);
		fdir_nref_2[im]->SetLineStyle(2);
		fdir_nref_2[im]->SetLineColor(1);

		h1d_same_dphi_long[im][0]->Fit(fdir_nref_2[im], "R0Q");
		fdir_nref_2[im]->Draw("same");

		hv22_dir_nref_2->SetBinContent(im+1, fdir_nref_2[im]->GetParameter(2));
		hv22_dir_nref_2->SetBinError(im+1, fdir_nref_2[im]->GetParError(2));

    if ( im>nref_2 ){

      g_h_central = h1d_same_dphi_long[im][0];
      g_h_peripheral = h1d_same_dphi_long[nref_2][0];
			g_h_template_fit = (TH1D*)g_h_peripheral->Clone(Form("h_template_fit_for_nref_2_mult%d",im));
			g_h_template_fit->Reset();

			func_nref_2[im] = new TF1("f_fit",fitf,-const_pi/2,1.5*const_pi,NHAR+1);
			func_nref_2[im]->SetParameter(0,1);
			for (int i=1; i<NHAR+1; i++) func_nref_2[im]->SetParameter(i,0.0);

			TH1D* l_h_central_ = (TH1D*)g_h_central->Clone("l_h_central_");
			TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);

			l_h_central_->Fit(func_nref_2[im],"QU0");

			TF1 *f_pedestal = new TF1("f_pedestal","[0]",-0.5*const_pi,1.5*const_pi);
			f_pedestal->SetLineColor(kOrange-3);
			f_pedestal->SetLineStyle(2);
			f_pedestal->SetLineWidth(3);
			f_pedestal->SetParameter(0,g_ped);
			f_pedestal->SetParError(0,func_nref_2[im]->GetParError(0)*g_h_peripheral->Integral()/g_h_peripheral->GetNbinsX());
      //f_pedestal->Draw("same");
      //
			//cout << "f_pedestal par: " << f_pedestal->GetParameter(0) << ", err: " << f_pedestal->GetParError(0) << endl;
			//cout << "f_scale par: " << func_nref_2[im]->GetParameter(0) << ", err: " << func_nref_2[im]->GetParError(0) << endl;

			h_rescaledperipheral_nref_2[im] = (TH1D*)h1d_same_dphi_long[nref_2][0]->Clone(Form("h_rescaledperipheral_nref_2_mult%d",im));
			h_rescaledperipheral_nref_2[im]->Scale(func_nref_2[im]->GetParameter(0));
			h_rescaledperipheral_nref_2[im]->Add(f_pedestal);
			h_rescaledperipheral_nref_2[im]->SetMarkerStyle(31);
			h_rescaledperipheral_nref_2[im]->SetLineWidth(3);

			f_ped_nref_2[im] = new TF1("f_ped","[0]",-0.5*const_pi,1.5*const_pi);
			//f_ped_nref_2[im]->SetParameter(0,h_rescaledperipheral_nref_2[im]->GetBinContent(h_rescaledperipheral_nref_2[im]->FindBin(0)));
			f_ped_nref_2[im]->SetParameter(0,h_rescaledperipheral_nref_2[im]->GetMinimum());
			f_ped_nref_2[im]->SetLineColor(kOrange-3);
			f_ped_nref_2[im]->SetLineStyle(2);
			f_ped_nref_2[im]->SetLineWidth(3);

			f_combined_nref_2[im] = new TF1("f_combined","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x)) + [4]",-0.5*const_pi,1.5*const_pi);
			f_combined_nref_2[im]->SetParameter(0, g_ped);
			f_combined_nref_2[im]->SetParameter(1, func_nref_2[im]->GetParameter(1));
			f_combined_nref_2[im]->SetParError(1, func_nref_2[im]->GetParError(1));
			f_combined_nref_2[im]->SetParameter(2, func_nref_2[im]->GetParameter(2));
			f_combined_nref_2[im]->SetParameter(3, func_nref_2[im]->GetParameter(3));
			//f_combined_nref_2[im]->SetParameter(4, h_rescaledperipheral_nref_2[im]->GetBinContent(h_rescaledperipheral_nref_2[im]->FindBin(0)) - g_ped);
			f_combined_nref_2[im]->SetParameter(4, h_rescaledperipheral_nref_2[im]->GetMinimum() - g_ped);
			f_combined_nref_2[im]->SetLineStyle(2);
			f_combined_nref_2[im]->SetLineColor(4);
			f_combined_nref_2[im]->SetLineWidth(3);


      func_nref_2[im]->Draw("same");
			f_ped_nref_2[im]->Draw("same");
      h_rescaledperipheral_nref_2[im]->Draw("p same");
			f_combined_nref_2[im]->Draw("same");

			float v22 = f_combined_nref_2[im]->GetParameter(1);
			float v22_err = f_combined_nref_2[im]->GetParError(1);

      gv22_sub_nref_2->SetPoint(g_count, 10*im+5, v22);
      gv22_sub_nref_2->SetPointError(g_count, 0, 0);

      g_count++;

    }

    {
      TLegend *leg = new TLegend(0.15,0.73,0.7,0.93);
      leg->SetFillStyle(0);
			leg->SetBorderSize(0);
      leg->AddEntry("",Form("%s 13 TeV",system),"");
      leg->AddEntry("",Form("%d#leqN_{TRK}^{CH}<%d",10*im,10*im+10),"");
      leg->AddEntry(h1d_same_dphi_long[im][0],"2<|#Delta#eta|<5","P");
      leg->Draw();
    }
	}//im

  c1_2->cd();
  c1_2->Update();

  g_count = 0;
	TCanvas *c1_3 = new TCanvas("c1_3","c1_3",1.1*200*5,200*3);
	c1_3->Divide(5,3);

	for (int im=0; im<nmult; im++){
		c1_3->cd(im+1);
		SetPadStyle();

		float max = h1d_same_dphi_long[im][0]->GetMaximum();
		float min = h1d_same_dphi_long[im][0]->GetMinimum();

    htmp = (TH1D*)gPad->DrawFrame(-const_pi/2, min-0.05*(max-min), 3*const_pi/2, max+0.05*(max-min));
		SetHistoStyle("#Delta#phi (rad)","Per trigger yield");
		h1d_same_dphi_long[im][0]->SetMarkerStyle(24);
		h1d_same_dphi_long[im][0]->Draw("p same");

		fdir_nref_3[im] = new TF1(Form("fdir_%d",im),"[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",-const_pi/2,3*const_pi/2);
		fdir_nref_3[im]->SetLineWidth(3);
		fdir_nref_3[im]->SetLineStyle(2);
		fdir_nref_3[im]->SetLineColor(1);

		h1d_same_dphi_long[im][0]->Fit(fdir_nref_3[im], "R0Q");
		fdir_nref_3[im]->Draw("same");

		hv22_dir_nref_3->SetBinContent(im+1, fdir_nref_3[im]->GetParameter(2));
		hv22_dir_nref_3->SetBinError(im+1, fdir_nref_3[im]->GetParError(2));

    if ( im>nref_3 ){

      g_h_central = h1d_same_dphi_long[im][0];
      g_h_peripheral = h1d_same_dphi_long[nref_3][0];
			g_h_template_fit = (TH1D*)g_h_peripheral->Clone(Form("h_template_fit_for_nref_3_mult%d",im));
			g_h_template_fit->Reset();

			func_nref_3[im] = new TF1("f_fit",fitf,-const_pi/2,1.5*const_pi,NHAR+1);
			func_nref_3[im]->SetParameter(0,1);
			for (int i=1; i<NHAR+1; i++) func_nref_3[im]->SetParameter(i,0.0);

			TH1D* l_h_central_ = (TH1D*)g_h_central->Clone("l_h_central_");
			TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);

			l_h_central_->Fit(func_nref_3[im],"QU0");

			TF1 *f_pedestal = new TF1("f_pedestal","[0]",-0.5*const_pi,1.5*const_pi);
			f_pedestal->SetLineColor(kOrange-3);
			f_pedestal->SetLineStyle(2);
			f_pedestal->SetLineWidth(3);
			f_pedestal->SetParameter(0,g_ped);
			f_pedestal->SetParError(0,func_nref_3[im]->GetParError(0)*g_h_peripheral->Integral()/g_h_peripheral->GetNbinsX());
      //f_pedestal->Draw("same");
      //
			//cout << "f_pedestal par: " << f_pedestal->GetParameter(0) << ", err: " << f_pedestal->GetParError(0) << endl;
			//cout << "f_scale par: " << func_nref_3[im]->GetParameter(0) << ", err: " << func_nref_3[im]->GetParError(0) << endl;

			h_rescaledperipheral_nref_3[im] = (TH1D*)h1d_same_dphi_long[nref_3][0]->Clone(Form("h_rescaledperipheral_nref_3_mult%d",im));
			h_rescaledperipheral_nref_3[im]->Scale(func_nref_3[im]->GetParameter(0));
			h_rescaledperipheral_nref_3[im]->Add(f_pedestal);
			h_rescaledperipheral_nref_3[im]->SetMarkerStyle(31);
			h_rescaledperipheral_nref_3[im]->SetLineWidth(3);

			f_ped_nref_3[im] = new TF1("f_ped","[0]",-0.5*const_pi,1.5*const_pi);
			//f_ped_nref_3[im]->SetParameter(0,h_rescaledperipheral_nref_3[im]->GetBinContent(h_rescaledperipheral_nref_3[im]->FindBin(0)));
			f_ped_nref_3[im]->SetParameter(0,h_rescaledperipheral_nref_3[im]->GetMinimum());
			f_ped_nref_3[im]->SetLineColor(kOrange-3);
			f_ped_nref_3[im]->SetLineStyle(2);
			f_ped_nref_3[im]->SetLineWidth(3);

			f_combined_nref_3[im] = new TF1("f_combined","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x)) + [4]",-0.5*const_pi,1.5*const_pi);
			f_combined_nref_3[im]->SetParameter(0, g_ped);
			f_combined_nref_3[im]->SetParameter(1, func_nref_3[im]->GetParameter(1));
			f_combined_nref_3[im]->SetParError(1, func_nref_3[im]->GetParError(1));
			f_combined_nref_3[im]->SetParameter(2, func_nref_3[im]->GetParameter(2));
			f_combined_nref_3[im]->SetParameter(3, func_nref_3[im]->GetParameter(3));
			//f_combined_nref_3[im]->SetParameter(4, h_rescaledperipheral_nref_3[im]->GetBinContent(h_rescaledperipheral_nref_3[im]->FindBin(0)) - g_ped);
			f_combined_nref_3[im]->SetParameter(4, h_rescaledperipheral_nref_3[im]->GetMinimum() - g_ped);
			f_combined_nref_3[im]->SetLineStyle(2);
			f_combined_nref_3[im]->SetLineColor(4);
			f_combined_nref_3[im]->SetLineWidth(3);


      func_nref_3[im]->Draw("same");
			f_ped_nref_3[im]->Draw("same");
      h_rescaledperipheral_nref_3[im]->Draw("p same");
			f_combined_nref_3[im]->Draw("same");

			float v22 = f_combined_nref_3[im]->GetParameter(1);
			float v22_err = f_combined_nref_3[im]->GetParError(1);

      gv22_sub_nref_3->SetPoint(g_count, 10*im+5, v22);
      gv22_sub_nref_3->SetPointError(g_count, 0, 0);
      
      g_count++;

    }

    {
      TLegend *leg = new TLegend(0.15,0.73,0.7,0.93);
      leg->SetFillStyle(0);
			leg->SetBorderSize(0);
      leg->AddEntry("",Form("%s 13 TeV",system),"");
      leg->AddEntry("",Form("%d#leqN_{TRK}^{CH}<%d",10*im,10*im+10),"");
      leg->AddEntry(h1d_same_dphi_long[im][0],"2<|#Delta#eta|<5","P");
      leg->Draw();
    }
	}//im

  c1_3->cd();
  c1_3->Update();

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

	hv22_dir_nref_1->SetMarkerStyle(24);
	hv22_dir_nref_1->SetMarkerColor(1);
	hv22_dir_nref_1->SetLineColor(1);
	hv22_dir_nref_1->SetLineWidth(2);
	hv22_dir_nref_1->Draw("p same");

	gv22_sub_nref_1->SetLineWidth(4);
	gv22_sub_nref_1->SetLineColor(1);
	gv22_sub_nref_1->Draw("c same");

	gv22_sub_nref_2->SetLineWidth(4);
	gv22_sub_nref_2->SetLineColor(2);
	gv22_sub_nref_2->Draw("c same");

	gv22_sub_nref_3->SetLineWidth(4);
	gv22_sub_nref_3->SetLineColor(4);
	gv22_sub_nref_3->Draw("c same");

	{
		TLegend *leg = new TLegend(0.45,0.45,0.85,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		le = leg->AddEntry("","PYTHIA8 p+p 13 TeV","");
		le->SetTextFont(62);
		le->SetTextSize(0.05);
		le = leg->AddEntry("","Method 1","");
		le->SetTextFont(32);
		le->SetTextSize(0.05);
		le = leg->AddEntry("","0.5<p_{T}^{Trig}, p_{T}^{Assoc}<5.0 GeV/c","");
		le->SetTextSize(0.045);
		le = leg->AddEntry("","0<|#eta^{Trig}|, |#eta^{Assoc}|<2.5","");
		le->SetTextSize(0.045);
		le = leg->AddEntry(hv22_dir_nref_1,"Fourier fits","P");
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
		TLatex *tex = new TLatex(10,0.014,"(a)");
		tex->SetTextSize(0.045);
		tex->SetTextFont(62);
		tex->Draw();
	}


  /*
  if ( bSAVE ){
    c1->cd();
    c1->SaveAs(Form("%s_atlas_1d.gif",dataset));

    c3->cd();
    c3->SaveAs(Form("%s_atlas_v22.gif",dataset));
  }
  */

}

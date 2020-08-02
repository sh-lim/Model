#include "TFile.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TString.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TColor.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TGraph.h"

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

void vnvsen(int nevents = 100, int cent = 0, double ptlow = 0.5, double pthigh = 4.0) {

  bool vnven = true; bool kvspt = true;
  
  //nevents is number of events per centrality (100 for all so far)
  //centrality integer is cent-cent+1, i.e. for 25-26% centrality enter 25
  //vnven: 1 to output graphs of vn vs en for 3 pT values, 0 to skip
  //kvspt: 1 to output graphs of k2, k3, and k3/k2 vs pT, 0 to skip
  
  char vnfilename[100];
  char enfilename[100]; 
  sprintf(enfilename,"100_Ev_for_SONIC_%d_%d_cent.txt",cent,cent+1);
  //enfilename = Form("100_ev_for_SONIC_%d_%d_cent.txt",cent,cent+1);
  int npt = 25; // number of pT steps in the vn file
  ifstream enfile;
  enfile.open(enfilename);
  if (!enfile) {
    cout << "no enfile found " << enfilename << endl;
    return;
  }
  string dummy;
  enfile >> dummy >> dummy >> dummy >> dummy;

  // storage and variables
  int event[nevents];
  double maxe2 = 0.0;
  double maxe3 = 0.0;    
  double e2[nevents];
  double e3[nevents];
  double pt[npt];
  double v2[nevents][npt];
  double v3[nevents][npt];
  double nchpart[nevents][npt];
  double ptdummy;
  double vndummy;
  double weightedv2sum[nevents];
  double weightedv3sum[nevents];
  double nchpartsum[nevents];
  double meanpt[nevents];

  TGraph *tnchpart = new TGraph();
  
  // read files, weighted sums
  for (int ievent=0; ievent<nevents; ievent++) {
    
    enfile >> event[ievent] >> e2[ievent] >> e3[ievent];
    if (e2[ievent] > maxe2) maxe2 = e2[ievent];
    if (e3[ievent] > maxe3) maxe3 = e3[ievent];    

    ifstream vnfile;
    sprintf(vnfilename,"unid_vn_fluc_%d%d/unid_vn_fluc_%d.dat",cent,cent+1,event[ievent]);    
    //    vnfilename = Form("unid_vn_fluc_%d%d/unid_vn_fluc_%d.dat",cent,cent+1,event[ievent]);
    vnfile.open(vnfilename);
    if (!vnfile) {
      cout << "no vnfile found" << vnfilename << endl;
      return;
    }
    vnfile >> dummy >> dummy >> dummy >> dummy >> dummy;
    
    for (int ipt = 0; ipt < npt; ipt++){

      if (ievent == 0) {
	vnfile >> pt[ipt];
      } else {
	vnfile >> ptdummy;
      }
      vnfile >> nchpart[ievent][ipt] >> vndummy >> vndummy >> v2[ievent][ipt] >> v3[ievent][ipt] >> vndummy >> vndummy;
      if (ipt == 0){
	weightedv2sum[ievent] = 0.0;
	weightedv3sum[ievent] = 0.0;
	nchpartsum[ievent] = 0.0;
	meanpt[ievent] = 0.0;
      }
      if (pt[ipt] >= ptlow && pt[ipt] <= pthigh){
	weightedv2sum[ievent] += (v2[ievent][ipt]*nchpart[ievent][ipt]);
	weightedv3sum[ievent] += (v3[ievent][ipt]*nchpart[ievent][ipt]);
	nchpartsum[ievent] += nchpart[ievent][ipt];
      }

      meanpt[ievent] += pt[ipt] * nchpart[ievent][ipt];
      tnchpart->SetPoint(ipt,pt[ipt],nchpart[ievent][ipt]);
      
    } // end loop over pT bins
    vnfile.close();
    cout << "meanpt = " << meanpt[ievent] / nchpartsum[ievent] << endl;
  } // end loop over events
  enfile.close();

  TCanvas *ctest = new TCanvas("ctest","ctest",10,10,600,600);
  ctest->Draw();
  tnchpart->SetMarkerStyle(20);
  tnchpart->Draw("A pl");
  tnchpart->Print();
  
  // calculate and write integrated vn 
  double intv2[nevents];
  double intv3[nevents];
  ofstream outfile;
  char* outfilename = Form("intv2v3_%d%d.txt",cent,cent+1);
  outfile.open(outfilename);
  for (int ievent=0; ievent<nevents; ievent++){
    intv2[ievent] = weightedv2sum[ievent]/nchpartsum[ievent] / e2[ievent];
    intv3[ievent] = weightedv3sum[ievent]/nchpartsum[ievent] / e3[ievent];
    outfile << intv2[ievent] << " " << intv3[ievent] << endl;
  }
  outfile.close();

  // graph setup and linear fits
  TGraph* graphs2[npt];
  TGraph* graphs3[npt];
  double parinfo2[6];
  double parinfo3[6];
  TGraphErrors* k2vpt = new TGraphErrors();
  TGraphErrors* k3vpt = new TGraphErrors();
  TGraphErrors* k3k2ratiovpt = new TGraphErrors();
  double k2;
  double k3;
  double k2err;
  double k3err;
  for (int ipt=0; ipt<npt; ipt++){
    graphs2[ipt] = new TGraph();
    graphs3[ipt] = new TGraph();
    for(int ievent=0; ievent<nevents; ievent++){
      if (!isnan(v2[ievent][ipt]) && !isnan(v3[ievent][ipt])){
        graphs2[ipt]->SetPoint(ievent,e2[ievent],v2[ievent][ipt]);
        graphs3[ipt]->SetPoint(ievent,e3[ievent],v3[ievent][ipt]);
      }
    }
    TF1* fit2 = new TF1(Form("f2_%d",ipt),"[0]*x",0.0,0.6);
    fit2->SetLineWidth(3);
    fit2->SetLineColor(kRed);
    graphs2[ipt]->Fit(fit2,"QR"); // jln - "Q" is quiet mode
    k2 = fit2->GetParameter(0);
    k2err = fit2->GetParError(0);
    k2vpt->SetPoint(ipt,pt[ipt],k2);
    k2vpt->SetPointError(ipt,0,k2err);
    TF1* fit3 = new TF1(Form("f3_%d",ipt),"[0]*x",0.0,0.6);
    fit3->SetLineWidth(3);
    fit3->SetLineColor(kRed);
    graphs3[ipt]->Fit(fit3,"QR");
    k3 = fit3->GetParameter(0);
    k3err = fit3->GetParError(0);
    k3vpt->SetPoint(ipt,pt[ipt],k3);
    k3vpt->SetPointError(ipt,0,k3err);
    if (ipt!=0) k3k2ratiovpt->SetPoint(ipt,pt[ipt],k3/k2);
    if (ipt!=0) k3k2ratiovpt->SetPointError(ipt,0,sqrt((pow(k2err,2)*pow(-1*k3/pow(k2,2),2))+(pow(k3err,2)*pow(1/k2,2))));
    if (pt[ipt]==0.56 || pt[ipt]==1.04 || pt[ipt]==1.52){
      parinfo2[(2*ipt/3)-2] = k2;
      parinfo2[(2*ipt/3)-1] = k2err;
      parinfo3[(2*ipt/3)-2] = k3;
      parinfo3[(2*ipt/3)-1] = k3err;
    }
  }

  //============================================================
  // draw vn vs en graphs
  if(vnven){

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas *c1 = new TCanvas("c1","c1",10,10,2300,1200);
    c1->Divide(3,2,1e-11,1e-11);
    double xmax = 1.0;
    double ymax = 0.35;
    if (cent==0) ymax = 0.20;
    for (int ipt=0; ipt<npt; ipt++){
      if(ipt==3||ipt==6||ipt==9){
        c1->cd(ipt/3);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.03);	
	gPad->SetBottomMargin(0.15);
	gPad->SetTicky(1);
	xmax = maxe2*1.05;	
	graphs2[ipt]->SetMarkerStyle(24);
        graphs2[ipt]->SetMarkerSize(1.3);
        graphs2[ipt]->SetMarkerColor(kBlue);
        graphs2[ipt]->SetLineColor(kBlue);
        graphs2[ipt]->GetXaxis()->SetLimits(0,xmax);
        graphs2[ipt]->GetYaxis()->SetRangeUser(0,ymax);
        graphs2[ipt]->GetXaxis()->SetTitle("#varepsilon_{2}");
        graphs2[ipt]->GetYaxis()->SetTitle("v_{2}");
	graphs2[ipt]->GetXaxis()->SetTitleSize(0.06);
	graphs2[ipt]->GetYaxis()->SetTitleSize(0.06);	
	graphs2[ipt]->GetXaxis()->SetLabelSize(0.05);
	graphs2[ipt]->GetYaxis()->SetLabelSize(0.05);
	graphs2[ipt]->GetYaxis()->SetNdivisions(5);
	if (cent==0) 	graphs2[ipt]->GetXaxis()->SetNdivisions(5);
        graphs2[ipt]->Draw("AP");

	double labelsize = 0.05;
	
        TLatex *tg2 = new TLatex();
        tg2->SetTextFont (62);
        tg2->SetTextSize(labelsize);
        tg2->SetTextColor (kBlack);
        tg2->DrawLatexNDC (0.18, 0.85, Form("Pb+Pb %d-%d%% (MAGMA I.C. + SONIC)",cent,cent+1));
	TLatex *tp2 = new TLatex();
	tp2->SetTextFont(62);
	tp2->SetTextSize(labelsize);
	tp2->SetTextColor(kBlack);
	tp2->DrawLatexNDC (0.18,0.79, Form("hadron p_{T} = %.1f GeV",pt[ipt]));
        TLatex *sl2 = new TLatex();
        sl2->SetTextFont(62);
        sl2->SetTextColor(kBlack);
        sl2->SetTextSize(labelsize);
        sl2->DrawLatexNDC(0.18,0.73,Form("#kappa_{2} = %.3f #pm %.3f",parinfo2[(2*ipt/3)-2],parinfo2[(2*ipt/3)-1]));
        TLatex *cent2 = new TLatex();
        cent2->SetTextFont(62);
        cent2->SetTextColor(kBlack);
        cent2->SetTextSize(0.04);
	//        cent2->DrawLatexNDC(0.2,0.8,Form("centrality %d - %d \%",cent,cent+1));
        TLatex *corr2 = new TLatex();
        corr2->SetTextFont(62);
        corr2->SetTextColor(kBlack);
        corr2->SetTextSize(labelsize);
        corr2->DrawLatexNDC(0.18,0.66,Form("Pearson coeff. = %.2f",graphs2[ipt]->GetCorrelationFactor()));

	xmax = maxe3*1.05;	
        c1->cd(ipt/3+3);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.03);	
	gPad->SetBottomMargin(0.15);
	gPad->SetTicky(1);
	graphs3[ipt]->SetMarkerStyle(24);
        graphs3[ipt]->SetMarkerSize(1.2);
        graphs3[ipt]->SetMarkerColor(kBlue);
        graphs3[ipt]->SetLineColor(kBlue);
        graphs3[ipt]->GetXaxis()->SetLimits(0,xmax);
        graphs3[ipt]->GetYaxis()->SetRangeUser(0,0.2);
	if (cent==0)         graphs3[ipt]->GetYaxis()->SetRangeUser(0,0.15);
        graphs3[ipt]->GetXaxis()->SetTitle("#varepsilon_{3}");
	graphs3[ipt]->GetXaxis()->SetTitleSize(0.06);
	graphs3[ipt]->GetYaxis()->SetTitleSize(0.06);	
	graphs3[ipt]->GetXaxis()->SetLabelSize(0.05);
	graphs3[ipt]->GetYaxis()->SetLabelSize(0.05);
	graphs3[ipt]->GetYaxis()->SetNdivisions(5);	
        graphs3[ipt]->GetYaxis()->SetTitle("v_{3}");
	if (cent==0) 	graphs2[ipt]->GetXaxis()->SetNdivisions(5);	
        graphs3[ipt]->Draw("AP");
        TLatex *tg3 = new TLatex();
        tg3->SetTextFont (62);
        tg3->SetTextSize(labelsize);
        tg3->SetTextColor (kBlack);
        tg3->DrawLatexNDC (0.18, 0.85, Form("Pb+Pb %d-%d%% (MAGMA I.C. + SONIC)",cent,cent+1));
	TLatex *tp3 = new TLatex();
	tp3->SetTextFont(62);
	tp3->SetTextSize(labelsize);
	tp3->SetTextColor(kBlack);
	tp3->DrawLatexNDC (0.18,0.79, Form("hadron p_{T} = %.1f GeV",pt[ipt]));
	TLatex *sl3 = new TLatex();
        sl3->SetTextFont(62);
        sl3->SetTextColor(kBlack);
        sl3->SetTextSize(labelsize);
        sl3->DrawLatexNDC(0.18,0.73,Form("#kappa_{3} = %.3f #pm %.3f",parinfo3[(2*ipt/3)-2],parinfo3[(2*ipt/3)-1]));
        TLatex *cent3 = new TLatex();
        cent3->SetTextFont(62);
        cent3->SetTextColor(kBlack);
        cent3->SetTextSize(labelsize);
	//        cent3->DrawLatexNDC(0.2,0.8,Form("centrality %d - %d #%",cent,cent+1));
        TLatex *corr3 = new TLatex();
        corr3->SetTextFont(62);
        corr3->SetTextColor(kBlack);
        corr3->SetTextSize(labelsize);
        corr3->DrawLatexNDC(0.18,0.66,Form("Pearson coeff. = %.2f",graphs3[ipt]->GetCorrelationFactor()));
      }
    } // end pt loop

    char saveout[100];
    sprintf(saveout,"figure_pbpb_%d%d_kappafits.png",cent,cent+1);
    c1->Print(saveout);
    
  } // end if (vnven)

  //============================================
  // draw kappa vs pt graphs
  if(kvspt){
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas* c3 = new TCanvas("c3","c3",10,10,700,1500);
    c3->Divide(1,3);
    c3->cd(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.02);    
    k2vpt->SetMarkerSize(1.1);
    k2vpt->SetMarkerStyle(20);    
    k2vpt->SetMarkerColor(kRed);
    k2vpt->SetLineColor(kRed);
    k2vpt->GetXaxis()->SetLimits(0,3.0);
    k2vpt->GetYaxis()->SetRangeUser(0,1.5);
    k2vpt->GetXaxis()->SetTitleSize(0.07);
    k2vpt->GetYaxis()->SetTitleSize(0.07);    
    k2vpt->GetXaxis()->SetLabelSize(0.06);
    k2vpt->GetYaxis()->SetLabelSize(0.06);    
    k2vpt->GetXaxis()->SetTitle("p_{T} (GeV)");
    k2vpt->GetYaxis()->SetTitle("#kappa_{2}");
    k2vpt->Draw("AP");
    TLatex *cen2 = new TLatex();
    cen2->SetTextFont(62);
    cen2->SetTextColor(kBlack);
    cen2->SetTextSize(0.06);
    cen2->DrawLatexNDC(0.18,0.92,Form("PbPb %d-%d%% (MAGMA I.C. + SONIC)",cent,cent+1));

    // jln - find the event average of the pT range
    double aveint = 0.0;
    for (int ievent=0;ievent<nevents;ievent++) {
      aveint += intv2[ievent];
    }
    aveint = aveint / (double) nevents;
    cout << "ave int k2 = " << aveint << endl;    
    TLine *tk2ave = new TLine(ptlow,aveint,pthigh,aveint);
    tk2ave->SetLineColor(kRed);
    tk2ave->SetLineWidth(3);
    //    tk2ave->SetLineStyle(4);
    tk2ave->Draw("l,same");
    
    c3->cd(2);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.02);        
    k3vpt->SetMarkerSize(1.1);
    k3vpt->SetMarkerStyle(20);
    k3vpt->SetMarkerColor(kBlue);
    k3vpt->SetLineColor(kBlue);
    k3vpt->GetXaxis()->SetLimits(0,3.0);
    k3vpt->GetYaxis()->SetRangeUser(0,1.5);
    k3vpt->GetXaxis()->SetTitleSize(0.07);
    k3vpt->GetYaxis()->SetTitleSize(0.07);    
    k3vpt->GetXaxis()->SetLabelSize(0.06);
    k3vpt->GetYaxis()->SetLabelSize(0.06);    
    k3vpt->GetXaxis()->SetTitle("p_{T} (GeV)");
    k3vpt->GetYaxis()->SetTitle("#kappa_{3}");
    k3vpt->Draw("AP");

    // jln - find the event average of the pT range
    aveint = 0.0;
    for (int ievent=0;ievent<nevents;ievent++) {
      aveint += intv3[ievent];
    }
    aveint = aveint / (double) nevents;
    cout << "ave int k3 = " << aveint << endl;
    TLine *tk3ave = new TLine(ptlow,aveint,pthigh,aveint);
    tk3ave->SetLineColor(kBlue);
    tk3ave->SetLineWidth(3);
    //    tk3ave->SetLineStyle(4);
    tk3ave->Draw("l,same");

    TLatex *cen3 = new TLatex();
    cen3->SetTextFont(62);
    cen3->SetTextColor(kBlack);
    cen3->SetTextSize(0.06);
    cen3->DrawLatexNDC(0.18,0.92,Form("PbPb %d-%d%% (MAGMA I.C. + SONIC)",cent,cent+1));

    c3->cd(3);
    gPad->SetTicky(1);    
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.02);        
    k3k2ratiovpt->SetMarkerColor(kGreen+2);
    k3k2ratiovpt->SetLineColor(kGreen+2);    
    k3k2ratiovpt->SetMarkerSize(1.1);    
    k3k2ratiovpt->SetMarkerStyle(20);    
    k3k2ratiovpt->GetXaxis()->SetLimits(0,3.0);
    k3k2ratiovpt->GetYaxis()->SetRangeUser(0,2.0);
    k3k2ratiovpt->GetXaxis()->SetTitleSize(0.07);
    k3k2ratiovpt->GetYaxis()->SetTitleSize(0.07);
    k3k2ratiovpt->GetXaxis()->SetLabelSize(0.06);
    k3k2ratiovpt->GetYaxis()->SetLabelSize(0.06);    
    k3k2ratiovpt->GetXaxis()->SetTitle("p_{T} (GeV)");
    k3k2ratiovpt->GetYaxis()->SetTitle("#kappa_{3} / #kappa_{2}");
    k3k2ratiovpt->Draw("AP");

    TLine *tl = new TLine(0.0,1.0,3.0,1.0);
    tl->SetLineStyle(2);
    tl->SetLineWidth(4);
    tl->Draw("l,same");
    
    TLatex *cen32 = new TLatex();
    cen32->SetTextFont(62);
    cen32->SetTextColor(kBlack);
    cen32->SetTextSize(0.06);
    cen32->DrawLatexNDC(0.18,0.92,Form("Pb+Pb %d-%d%% (MAGMA I.C. + SONIC)",cent,cent+1));

    char saveout[100];
    sprintf(saveout,"figure_pbpb_%d%d_kappaptdepend.png",cent,cent+1);
    c3->Print(saveout);

  } // end if (kvspt)
}

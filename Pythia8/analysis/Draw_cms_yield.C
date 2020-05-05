#include "Style.h"

void Draw_cms_yield(){

	gStyle->SetOptStat(0);

	TFile *infile = new TFile("../CMS-kinematics/outfile_hist_pp13TeV_set00_grp001_try001.root","read");

	const float const_pi = TMath::Pi();

	const int npt = 4;
	const float ptbin[npt+1] = {0.2, 1.0, 2.0, 3.0, 4.0};

	const int nmult = 4;
	const float multbin[nmult+1] = {0, 35, 80, 105, 150};

	const int nMarker[nmult] = {20, 21, 24, 25};
	const int nColor[nmult] = {1, 2, 4, 8};

	TH2D *h2d_same[nmult][npt];
	TH2D *h2d_mixed[nmult][npt];

	TH1D *hntrig_same[nmult];
	TH1D *hntrig_mixed[nmult];

	TH1D *h1d_deta_same[nmult][npt];
	TH1D *h1d_deta_mixed[nmult][npt];

	TH1D *h1d_dphi_same[nmult][npt];
	TH1D *h1d_dphi_mixed[nmult][npt];

	TH1D *h1d_dphi_zyam[nmult][npt];

	TF1 *f1d_dphi[nmult][npt];

	TH1D *h1d_Yassociated_pT[nmult];
	TH1D *h1d_Yassociated_mult[npt];

	for (int imult=0; imult<nmult; imult++){
		h1d_Yassociated_pT[imult] = new TH1D(Form("h1d_Yassociated_pT_mult%02d",imult),"",npt,ptbin);
	}

	for (int ipt=0; ipt<npt; ipt++){
		h1d_Yassociated_mult[ipt] = new TH1D(Form("h1d_Yassociated_mult_pt%02d",ipt),"",nmult,multbin);
	}

	for (int imult=0; imult<nmult; imult++){
		hntrig_same[imult] = (TH1D*)infile->Get(Form("hntrig_same_mult%02d",imult));
		hntrig_mixed[imult] = (TH1D*)infile->Get(Form("hntrig_mixed_mult%02d",imult));

		for (int ipt=0; ipt<npt; ipt++){

			h2d_same[imult][ipt] = (TH2D*)infile->Get(Form("h2d_same_dphi_deta_mult%02d_pt%02d",imult,ipt));
			h2d_mixed[imult][ipt] = (TH2D*)infile->Get(Form("h2d_mixed_dphi_deta_mult%02d_pt%02d",imult,ipt));

			h2d_same[imult][ipt]->RebinX(4);
			h2d_mixed[imult][ipt]->RebinX(4);

			float ntrig_same = hntrig_same[imult]->Integral(hntrig_same[imult]->FindBin(ptbin[ipt]+0.1), hntrig_same[imult]->FindBin(ptbin[ipt]+0.1));
			float ntrig_mixed = hntrig_mixed[imult]->Integral(hntrig_mixed[imult]->FindBin(ptbin[ipt]+0.1), hntrig_mixed[imult]->FindBin(ptbin[ipt]+0.1));

			float nnorm_mixed = h2d_mixed[imult][ipt]->GetBinContent(h2d_mixed[imult][ipt]->FindBin(0,0));
			nnorm_mixed /= ntrig_mixed;
			nnorm_mixed /= h2d_mixed[imult][ipt]->GetXaxis()->GetBinWidth(1);
			nnorm_mixed /= h2d_mixed[imult][ipt]->GetYaxis()->GetBinWidth(1);

			//1D projection, negative eta 
			int etabin_min = h2d_same[imult][ipt]->GetYaxis()->FindBin(-5.0+0.001);
			int etabin_max = h2d_same[imult][ipt]->GetYaxis()->FindBin(-2.0-0.001);

			h1d_dphi_same[imult][ipt] = (TH1D*)h2d_same[imult][ipt]->ProjectionX(Form("h1d_dphi_same_mult%02d_pt%02d",imult,ipt),etabin_min,etabin_max);
			h1d_dphi_mixed[imult][ipt] = (TH1D*)h2d_mixed[imult][ipt]->ProjectionX(Form("h1d_dphi_mixed_mult%02d_pt%02d",imult,ipt),etabin_min,etabin_max);

			//1D projection, positive eta 
			etabin_min = h2d_same[imult][ipt]->GetYaxis()->FindBin(+2.0+0.001);
			etabin_max = h2d_same[imult][ipt]->GetYaxis()->FindBin(+5.0-0.001);

			TH1D *htmp_same = (TH1D*)h2d_same[imult][ipt]->ProjectionX(Form("h1d_dphi_same_mult%02d_pt%02d_1",imult,ipt),etabin_min,etabin_max);
			TH1D *htmp_mixed = (TH1D*)h2d_mixed[imult][ipt]->ProjectionX(Form("h1d_dphi_mixed_mult%02d_pt%02d_1",imult,ipt),etabin_min,etabin_max);

			h1d_dphi_same[imult][ipt]->Add(htmp_same);
			h1d_dphi_mixed[imult][ipt]->Add(htmp_mixed);

			//normalization
			h2d_same[imult][ipt]->RebinY(4);
			h2d_mixed[imult][ipt]->RebinY(4);
			h2d_same[imult][ipt]->Scale(1./ntrig_same);
			h2d_mixed[imult][ipt]->Scale(1./ntrig_mixed);
			h2d_same[imult][ipt]->Divide(h2d_mixed[imult][ipt]);
			h2d_same[imult][ipt]->Scale(nnorm_mixed);

			h1d_dphi_same[imult][ipt]->Scale(1./ntrig_same);
			h1d_dphi_mixed[imult][ipt]->Scale(1./ntrig_mixed);
			h1d_dphi_same[imult][ipt]->Divide(h1d_dphi_mixed[imult][ipt]);
			h1d_dphi_same[imult][ipt]->Scale(nnorm_mixed);

			//fit w/ Fourier series
			//f1d_dphi[imult][ipt] = new TF1("f1","[0]*( 1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x))",-const_pi/2,3*const_pi/2);
			f1d_dphi[imult][ipt] = new TF1("f1","[0]*( 1 + 2*[1]*cos(x) + 2*[2]*cos(2*x))",-const_pi/2,1*const_pi/2);
			h1d_dphi_same[imult][ipt]->Fit(f1d_dphi[imult][ipt],"R0Q");

			//ZYAM subtraction
			float zyam = f1d_dphi[imult][ipt]->GetMinimum(-const_pi/2,const_pi/2);
			float zyam_x = f1d_dphi[imult][ipt]->GetMinimumX(-const_pi/2,const_pi/2);
			float Y_associated = 0.0; 
			float Y_associated_err = 0.0; 

			h1d_dphi_zyam[imult][ipt] = (TH1D*)h1d_dphi_same[imult][ipt]->Clone(Form("h1d_dphi_zyam_mult%02d_pt%02d",imult,ipt));
			for (int iphi=0; iphi<h1d_dphi_zyam[imult][ipt]->GetNbinsX(); iphi++){
				float val = h1d_dphi_zyam[imult][ipt]->GetBinContent(iphi+1);
				h1d_dphi_zyam[imult][ipt]->SetBinContent(iphi+1, val-zyam);

				//associated yield |dphi|<1.2
				float dphi = h1d_dphi_zyam[imult][ipt]->GetBinCenter(iphi+1); 
				float ddphi = h1d_dphi_zyam[imult][ipt]->GetBinWidth(iphi+1);
				//if ( fabs(dphi)<1.2 ){
				if ( fabs(dphi)<fabs(zyam_x) ){
					Y_associated += h1d_dphi_zyam[imult][ipt]->GetBinContent(iphi+1)*ddphi;
					Y_associated_err += h1d_dphi_zyam[imult][ipt]->GetBinError(iphi+1)*h1d_dphi_zyam[imult][ipt]->GetBinError(iphi+1)*ddphi;
				}
			}

			h1d_Yassociated_pT[imult]->SetBinContent(ipt+1,Y_associated); 
			h1d_Yassociated_pT[imult]->SetBinError(ipt+1,sqrt(Y_associated_err)); 

			h1d_Yassociated_mult[ipt]->SetBinContent(imult+1,Y_associated); 
			h1d_Yassociated_mult[ipt]->SetBinError(imult+1,sqrt(Y_associated_err)); 

		}//ipt
	}//imult

	TCanvas *cfig1 = new TCanvas("cfig1","cfig1",1.1*300*4,300*4);
	cfig1->Divide(4,4);

	for (int imult=0; imult<nmult; imult++){
		for (int ipt=0; ipt<npt; ipt++){
			cfig1->cd(imult*npt+ipt+1);
			SetPadStyle();
			gPad->SetLeftMargin(0.23);
			gPad->SetPhi(135);

			htmp = (TH1D*)h2d_same[imult][ipt]; 
			SetHistoStyle("","","",14,12);
			h2d_same[imult][ipt]->SetAxisRange(-4.0+0.01,4.0-0.01,"Y");
			h2d_same[imult][ipt]->GetYaxis()->SetTitle("|#Delta#eta|");
			h2d_same[imult][ipt]->GetYaxis()->CenterTitle();
			h2d_same[imult][ipt]->GetYaxis()->SetTitleOffset(5.);
			h2d_same[imult][ipt]->GetXaxis()->SetTitle("|#Delta#phi|");
			h2d_same[imult][ipt]->GetXaxis()->CenterTitle();
			h2d_same[imult][ipt]->GetXaxis()->SetTitleOffset(5.0);
			h2d_same[imult][ipt]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi}");
			h2d_same[imult][ipt]->GetZaxis()->CenterTitle();
			h2d_same[imult][ipt]->GetZaxis()->SetTitleOffset(7.0);

			float ymax = h2d_same[imult][ipt]->GetMaximum();
			float ymin = h2d_same[imult][ipt]->GetMinimum();
			h2d_same[imult][ipt]->SetMaximum(ymin + (0.8-0.2*ipt)*(ymax-ymin));
			h2d_same[imult][ipt]->Draw("surf1");

			TLegend *leg = new TLegend(0.05,0.8,0.5,0.98);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextFont(43);
			leg->SetTextSize(14);
			leg->AddEntry("","Pythia8 pp 13 TeV","h");
			leg->AddEntry("",Form("%d#leqN_{trk}<%d",int(multbin[imult]),int(multbin[imult+1])),"h");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"h");
			leg->Draw();
		}
	}

	TCanvas *cfig2_pre = new TCanvas("cfig2_pre","cfig2_pre",1.1*300*4,300*4);
	cfig2_pre->Divide(4,4);

	for (int imult=0; imult<nmult; imult++){
		for (int ipt=0; ipt<npt; ipt++){

			cfig2_pre->cd(imult*npt+ipt+1);
			SetPadStyle();
			gPad->SetLeftMargin(0.23);

			float ymax = h1d_dphi_same[imult][ipt]->GetMaximum();
			float ymin = h1d_dphi_same[imult][ipt]->GetMinimum();

			htmp = (TH1D*)gPad->DrawFrame(-const_pi/2,ymin-0.05*(ymax-ymin),const_pi*3/2,ymax+0.05*(ymax-ymin));
			SetHistoStyle("|#Delta#phi|","#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi}","",14,12);
			htmp->GetYaxis()->SetTitleOffset(7.0);
			htmp->GetXaxis()->SetTitleOffset(5.0);

			h1d_dphi_same[imult][ipt]->SetMarkerStyle(24);
			h1d_dphi_same[imult][ipt]->SetMarkerSize(1.0);
			h1d_dphi_same[imult][ipt]->Draw("same");

			f1d_dphi[imult][ipt]->SetLineWidth(3);
			f1d_dphi[imult][ipt]->SetLineStyle(2);
			f1d_dphi[imult][ipt]->Draw("same");

			float zyam_x = f1d_dphi[imult][ipt]->GetMinimumX(-const_pi/2,const_pi/2);

			TLegend *leg = new TLegend(0.25,0.73,0.65,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextFont(43);
			leg->SetTextSize(14);
			leg->AddEntry("","Pythia8 pp 13 TeV","h");
			leg->AddEntry("",Form("%d#leqN_{trk}<%d",int(multbin[imult]),int(multbin[imult+1])),"h");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"h");
			leg->AddEntry("",Form("#Delta #phi_{ZYAM}=%4.2f",fabs(zyam_x)),"h");
			leg->Draw();

		}
	}

	TCanvas *cfig2 = new TCanvas("cfig2","cfig2",1.1*300*4,300*4);
	cfig2->Divide(4,4);

	for (int imult=0; imult<nmult; imult++){
		for (int ipt=0; ipt<npt; ipt++){

			cfig2->cd(imult*npt+ipt+1);
			SetPadStyle();
			gPad->SetLeftMargin(0.23);

			float ymax = h1d_dphi_zyam[imult][ipt]->GetMaximum();
			float ymin = h1d_dphi_zyam[imult][ipt]->GetMinimum();

			htmp = (TH1D*)gPad->DrawFrame(-const_pi/2,-0.05*ymax,const_pi*3/2,1.05*ymax);
			SetHistoStyle("|#Delta#phi|","#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#etad#Delta#phi} - C_{ZYAM}","",14,12);
			htmp->GetYaxis()->SetTitleOffset(7.0);
			htmp->GetXaxis()->SetTitleOffset(5.0);

			h1d_dphi_zyam[imult][ipt]->SetMarkerStyle(24);
			h1d_dphi_zyam[imult][ipt]->SetMarkerSize(1.0);
			h1d_dphi_zyam[imult][ipt]->Draw("same");

			TLegend *leg = new TLegend(0.25,0.78,0.65,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextFont(43);
			leg->SetTextSize(14);
			leg->AddEntry("","Pythia8 pp 13 TeV","h");
			leg->AddEntry("",Form("%d#leqN_{trk}<%d",int(multbin[imult]),int(multbin[imult+1])),"h");
			leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"h");
			leg->Draw();

		}
	}

	TCanvas *cfig3 = new TCanvas("cfig3","cfig3",1.1*2*400,400);
	cfig3->Divide(2,1);

	cfig3->cd(1);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,-0.01,4,0.05);
	SetHistoStyle("p_{T} (GeV/c)","Associated yield/(GeV/c)","",20,16);

	{
		TLegend *leg = new TLegend(0.2,0.65,0.65,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");

		for (int imult=0; imult<nmult; imult++){
			h1d_Yassociated_pT[imult]->SetMarkerStyle(nMarker[imult]);
			h1d_Yassociated_pT[imult]->SetLineColor(nColor[imult]);
			h1d_Yassociated_pT[imult]->SetMarkerColor(nColor[imult]);
			h1d_Yassociated_pT[imult]->Draw("p same");
			leg->AddEntry(h1d_Yassociated_pT[imult],Form("%d#leqN_{trk}<%d",int(multbin[imult]),int(multbin[imult+1])),"PL");
		}
		leg->Draw();
	}


	cfig3->cd(2);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,-0.01,150,0.05);
	SetHistoStyle("Multiplicity","Associated yield/(GeV/c)","",20,16);

	{
		TLegend *leg = new TLegend(0.2,0.65,0.65,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry("","Pythia8 pp 13 TeV","h");

		for (int ipt=0; ipt<npt; ipt++){
			h1d_Yassociated_mult[ipt]->SetMarkerStyle(nMarker[ipt]);
			h1d_Yassociated_mult[ipt]->SetLineColor(nColor[ipt]);
			h1d_Yassociated_mult[ipt]->SetMarkerColor(nColor[ipt]);
			h1d_Yassociated_mult[ipt]->Draw("p same");
			leg->AddEntry(h1d_Yassociated_mult[ipt],Form("%g<p_{T}<%g (GeV/c)",ptbin[ipt],ptbin[ipt+1]),"PL");
		}
		leg->Draw();
	}

}

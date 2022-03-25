#include "Style.h"

void Draw_pPb5TeV00(){

	const bool bSAVE = false;

	gStyle->SetOptStat(0);

	const int nset = 2;
	const int nmult = 4;

	const int nColor[nmult] = {1, 2, 4, 6};

	const float c_alpha = 0.7;

	TFile *infile[nset];
	infile[0] = new TFile("outfile_hist_pPb5TeV_set30_grp000_try102.root","read");
	infile[1] = new TFile("outfile_hist_pPb5TeV_set31_grp000_try102.root","read");

	TH1D *hpt_pion[nset][nmult];
	TH1D *hpt_kaon[nset][nmult];
	TH1D *hpt_rho[nset][nmult];
	TH1D *hpt_kstar[nset][nmult];
	TH1D *hpt_phi[nset][nmult];

	TH1D *hmass_rho[nset];
	TH1D *hmass_phi[nset];
	TH1D *hmass_kstar[nset];

	for (int iset=0; iset<nset; iset++){
		for (int im=0; im<nmult; im++){

			hpt_pion[iset][im] = (TH1D*)infile[iset]->Get(Form("hpt_pion_m%d",im));
			hpt_kaon[iset][im] = (TH1D*)infile[iset]->Get(Form("hpt_kaon_m%d",im));
			hpt_rho[iset][im] = (TH1D*)infile[iset]->Get(Form("hpt_rho_m%d",im));
			hpt_kstar[iset][im] = (TH1D*)infile[iset]->Get(Form("hpt_kstar_m%d",im));
			hpt_phi[iset][im] = (TH1D*)infile[iset]->Get(Form("hpt_phi_m%d",im));

			hpt_pion[iset][im]->Sumw2();
			hpt_kaon[iset][im]->Sumw2();
			hpt_rho[iset][im]->Sumw2();
			hpt_kstar[iset][im]->Sumw2();
			hpt_phi[iset][im]->Sumw2();

		}//im

		hmass_rho[iset] = (TH1D*)infile[iset]->Get("hmass_rho");
		hmass_phi[iset] = (TH1D*)infile[iset]->Get("hmass_phi");
		hmass_kstar[iset] = (TH1D*)infile[iset]->Get("hmass_kstar");

		hmass_rho[iset]->SetLineColor(iset+1);
		hmass_phi[iset]->SetLineColor(iset+1);
		hmass_kstar[iset]->SetLineColor(iset+1);
	}//iset


	TH1D *hpt_kaon_pion[nset][nmult];
	TH1D *hpt_rho_pion[nset][nmult];
	TH1D *hpt_kstar_pion[nset][nmult];
	TH1D *hpt_phi_pion[nset][nmult];

	TH1D *hcp_kaon_pion[nset][nmult];
	TH1D *hcp_rho_pion[nset][nmult];
	TH1D *hcp_kstar_pion[nset][nmult];
	TH1D *hcp_phi_pion[nset][nmult];

	for (int iset=0; iset<nset; iset++){
		for (int im=0; im<nmult; im++){
			hpt_kaon_pion[iset][im] = (TH1D*)hpt_kaon[iset][im]->Clone(Form("hpt_kaon_pion_%d_%d",iset,im));
			hpt_rho_pion[iset][im] = (TH1D*)hpt_rho[iset][im]->Clone(Form("hpt_rho_pion_%d_%d",iset,im));
			hpt_kstar_pion[iset][im] = (TH1D*)hpt_kstar[iset][im]->Clone(Form("hpt_kstar_pion_%d_%d",iset,im));
			hpt_phi_pion[iset][im] = (TH1D*)hpt_phi[iset][im]->Clone(Form("hpt_phi_pion_%d_%d",iset,im));

			hpt_kaon_pion[iset][im]->Divide(hpt_pion[iset][im]);
			hpt_rho_pion[iset][im]->Divide(hpt_pion[iset][im]);
			hpt_kstar_pion[iset][im]->Divide(hpt_pion[iset][im]);
			hpt_phi_pion[iset][im]->Divide(hpt_pion[iset][im]);

			hpt_kaon_pion[iset][im]->SetMarkerColorAlpha(nColor[im],c_alpha);
			hpt_kaon_pion[iset][im]->SetLineColorAlpha(nColor[im],c_alpha);
			hpt_kaon_pion[iset][im]->SetFillColorAlpha(nColor[im],c_alpha);
			hpt_kaon_pion[iset][im]->SetFillStyle(1);

			hpt_rho_pion[iset][im]->SetMarkerColorAlpha(nColor[im],c_alpha);
			hpt_rho_pion[iset][im]->SetLineColorAlpha(nColor[im],c_alpha);
			hpt_rho_pion[iset][im]->SetFillColorAlpha(nColor[im],c_alpha);
			hpt_rho_pion[iset][im]->SetFillStyle(1);

			hpt_kstar_pion[iset][im]->SetMarkerColorAlpha(nColor[im],c_alpha);
			hpt_kstar_pion[iset][im]->SetLineColorAlpha(nColor[im],c_alpha);
			hpt_kstar_pion[iset][im]->SetFillColorAlpha(nColor[im],c_alpha);
			hpt_kstar_pion[iset][im]->SetFillStyle(1);

			hpt_phi_pion[iset][im]->SetMarkerColorAlpha(nColor[im],c_alpha);
			hpt_phi_pion[iset][im]->SetLineColorAlpha(nColor[im],c_alpha);
			hpt_phi_pion[iset][im]->SetFillColorAlpha(nColor[im],c_alpha);
			hpt_phi_pion[iset][im]->SetFillStyle(1);
		}//im
	}//iset

	for (int iset=0; iset<nset; iset++){
		for (int im=0; im<nmult; im++){
			hcp_kaon_pion[iset][im] = (TH1D*)hpt_kaon_pion[iset][im]->Clone(Form("hcp_kaon_pion_%d_%d",iset,im));
			hcp_rho_pion[iset][im] = (TH1D*)hpt_rho_pion[iset][im]->Clone(Form("hcp_rho_pion_%d_%d",iset,im));
			hcp_kstar_pion[iset][im] = (TH1D*)hpt_kstar_pion[iset][im]->Clone(Form("hcp_kstar_pion_%d_%d",iset,im));
			hcp_phi_pion[iset][im] = (TH1D*)hpt_phi_pion[iset][im]->Clone(Form("hcp_phi_pion_%d_%d",iset,im));

			hcp_kaon_pion[iset][im]->Divide(hpt_kaon_pion[iset][0]);
			hcp_rho_pion[iset][im]->Divide(hpt_rho_pion[iset][0]);
			hcp_kstar_pion[iset][im]->Divide(hpt_kstar_pion[iset][0]);
			hcp_phi_pion[iset][im]->Divide(hpt_phi_pion[iset][0]);
		}
	}

	TCanvas *c0_kaon = new TCanvas("c0_kaon","c0_kaon",1.2*2*500,500);
	c0_kaon->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c0_kaon->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,10,0.5);
		SetHistoStyle("p_{T} (GeV/c)","(K^{+}+K^{-})/(#pi^{+}+#pi^{-})","",22,20);

		for (int im=0; im<nmult; im++){
			hpt_kaon_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.2,0.9,0.2+0.05*6);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry(hpt_kaon_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_kaon_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_kaon_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_kaon_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	TCanvas *c1_kaon = new TCanvas("c1_kaon","c1_kaon",1.2*2*500,500);
	c1_kaon->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c1_kaon->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0.4,10,1.6);
		SetHistoStyle("p_{T} (GeV/c)","Ratio to 60-100%","",22,20);

		for (int im=0; im<nmult; im++){
			hcp_kaon_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.93-0.05*7,0.9,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry("","(K^{+}+K^{-})/(#pi^{+}+#pi^{-})","");
		leg->AddEntry(hpt_kaon_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_kaon_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_kaon_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_kaon_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	return;

	TCanvas *c0_rho = new TCanvas("c0_rho","c0_rho",1.2*2*500,500);
	c0_rho->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c0_rho->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,10,0.5);
		SetHistoStyle("p_{T} (GeV/c)","#rho^{0}/(#pi^{+}+#pi^{-})","",22,20);

		for (int im=0; im<nmult; im++){
			hpt_rho_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.2,0.9,0.2+0.05*6);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry(hpt_rho_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_rho_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_rho_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_rho_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	TCanvas *c1_rho = new TCanvas("c1_rho","c1_rho",1.2*2*500,500);
	c1_rho->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c1_rho->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0.4,10,1.6);
		SetHistoStyle("p_{T} (GeV/c)","Ratio to 60-100%","",22,20);

		for (int im=0; im<nmult; im++){
			hcp_rho_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.93-0.05*7,0.9,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry("","#rho^{0}/(#pi^{+}+#pi^{-})","");
		leg->AddEntry(hpt_rho_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_rho_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_rho_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_rho_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	TCanvas *c2_rho = new TCanvas("c2_rho","c2_rho",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)hmass_rho[0];
	SetHistoStyle("Mass (GeV/c^{2})","","",22,22);

	hmass_rho[0]->Draw();
	hmass_rho[1]->Draw("same");

	{
		TLegend *leg = new TLegend(0.53,0.93-0.05*5,0.95,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr","");
		leg->AddEntry("","p-Pb 5.02 TeV","");
		leg->AddEntry("","#rho^{0}#rightarrow#pi^{+}+#pi^{-}","");
		leg->AddEntry(hmass_rho[0],"Hadron rescattering off","L");
		leg->AddEntry(hmass_rho[1],"Hadron rescattering on","L");
		leg->Draw();
	}

	TCanvas *c0_kstar = new TCanvas("c0_kstar","c0_kstar",1.2*2*500,500);
	c0_kstar->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c0_kstar->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,10,0.15);
		SetHistoStyle("p_{T} (GeV/c)","(K^{*0}+#bar{K^{*0}})/(#pi^{+}+#pi^{-})","",22,20);

		for (int im=0; im<nmult; im++){
			hpt_kstar_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.2,0.9,0.2+0.05*6);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry(hpt_kstar_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_kstar_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_kstar_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_kstar_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset


	TCanvas *c1_kstar = new TCanvas("c1_kstar","c1_kstar",1.2*2*500,500);
	c1_kstar->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c1_kstar->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0.4,10,1.6);
		SetHistoStyle("p_{T} (GeV/c)","Ratio to 60-100%","",22,20);

		for (int im=0; im<nmult; im++){
			hcp_kstar_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.93-0.05*7,0.9,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry("","(K^{*0}+#bar{K^{*0}})/(#pi^{+}+#pi^{-})","");
		leg->AddEntry(hpt_kstar_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_kstar_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_kstar_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_kstar_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	TCanvas *c2_kstar = new TCanvas("c2_kstar","c2_kstar",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)hmass_kstar[0];
	SetHistoStyle("Mass (GeV/c^{2})","","",22,22);
	htmp->SetAxisRange(0.65,1.2);

	hmass_kstar[0]->Draw();
	hmass_kstar[1]->Draw("same");

	{
		TLegend *leg = new TLegend(0.53,0.93-0.05*5,0.95,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr","");
		leg->AddEntry("","p-Pb 5.02 TeV","");
		leg->AddEntry("","K^{*0}#rightarrow#pi^{-}+K^{+}","");
		leg->AddEntry(hmass_kstar[0],"Hadron rescattering off","L");
		leg->AddEntry(hmass_kstar[1],"Hadron rescattering on","L");
		leg->Draw();
	}

	//return;

	TCanvas *c0_phi = new TCanvas("c0_phi","c0_phi",1.2*2*500,500);
	c0_phi->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c0_phi->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,10,0.08);
		SetHistoStyle("p_{T} (GeV/c)","#phi/(#pi^{+}+#pi^{-})","",22,20);

		for (int im=0; im<nmult; im++){
			hpt_phi_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.2,0.9,0.2+0.05*6);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry(hpt_phi_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_phi_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_phi_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_phi_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	TCanvas *c1_phi = new TCanvas("c1_phi","c1_phi",1.2*2*500,500);
	c1_phi->Divide(2,1);
	for (int iset=0; iset<nset; iset++){
		c1_phi->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0.4,10,1.6);
		SetHistoStyle("p_{T} (GeV/c)","Ratio to 60-100%","",22,20);

		for (int im=0; im<nmult; im++){
			hcp_phi_pion[iset][im]->Draw("e3 same");
		}

		TLegend *leg = new TLegend(0.3,0.93-0.05*7,0.9,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr p-Pb 5.02 TeV","");
		if ( iset==0 ){
			leg->AddEntry("","Hadron rescattering off","");
		}else{
			leg->AddEntry("","Hadron rescattering on","");
		}
		leg->AddEntry("","#phi/(#pi^{+}+#pi^{-})","");
		leg->AddEntry(hpt_phi_pion[iset][3],"0-20%","F");
		leg->AddEntry(hpt_phi_pion[iset][2],"20-40%","F");
		leg->AddEntry(hpt_phi_pion[iset][1],"40-60%","F");
		leg->AddEntry(hpt_phi_pion[iset][0],"60-100%","F");
		leg->Draw();
	}//iset

	TCanvas *c2_phi = new TCanvas("c2_phi","c2_phi",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)hmass_phi[0];
	SetHistoStyle("Mass (GeV/c^{2})","","",22,22);
	htmp->SetAxisRange(0.65,1.2);

	hmass_phi[0]->Draw();
	hmass_phi[1]->Draw("same");

	{
		TLegend *leg = new TLegend(0.53,0.93-0.05*5,0.95,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","PYTHIA/Angantyr","");
		leg->AddEntry("","p-Pb 5.02 TeV","");
		leg->AddEntry("","#phi#rightarrowK^{+}+K^{-}","");
		leg->AddEntry(hmass_phi[0],"Hadron rescattering off","L");
		leg->AddEntry(hmass_phi[1],"Hadron rescattering on","L");
		leg->Draw();
	}

	if ( bSAVE ){
		c0_rho->cd();
		c0_rho->SaveAs("plots/PYTHIA_pPb5TeV_c0_rho.eps");

		c1_rho->cd();
		c1_rho->SaveAs("plots/PYTHIA_pPb5TeV_c1_rho.eps");

		c2_rho->cd();
		c2_rho->SaveAs("plots/PYTHIA_pPb5TeV_c2_rho.eps");

		c0_kstar->cd();
		c0_kstar->SaveAs("plots/PYTHIA_pPb5TeV_c0_kstar.eps");

		c1_kstar->cd();
		c1_kstar->SaveAs("plots/PYTHIA_pPb5TeV_c1_kstar.eps");

		c2_kstar->cd();
		c2_kstar->SaveAs("plots/PYTHIA_pPb5TeV_c2_kstar.eps");

		c0_phi->cd();
		c0_phi->SaveAs("plots/PYTHIA_pPb5TeV_c0_phi.eps");

		c1_phi->cd();
		c1_phi->SaveAs("plots/PYTHIA_pPb5TeV_c1_phi.eps");

		c2_phi->cd();
		c2_phi->SaveAs("plots/PYTHIA_pPb5TeV_c2_phi.eps");
	}

}

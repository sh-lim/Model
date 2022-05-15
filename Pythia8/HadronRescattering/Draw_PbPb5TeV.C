#include "Style.h"

void Draw_PbPb5TeV(){

	const int nset = 2;
	const int nmult = 9;
	const int cent_arr[nmult+1] = {100, 80, 70, 60, 50, 40, 30, 20, 10, 0};

	gStyle->SetOptStat(0);

	TFile *infile[nset];
	infile[0] = new TFile("outfile_hist_PbPb5TeV_set0.root","read");
	infile[1] = new TFile("outfile_hist_PbPb5TeV_set1.root","read");

	TH2D *h2d_rho[nset];
	TH2D *h2d_kstar[nset];
	TH2D *h2d_phi[nset];

	TH1D *h1d_rho[nset][nmult];
	TH1D *h1d_kstar[nset][nmult];
	TH1D *h1d_phi[nset][nmult];

	for (int iset=0; iset<nset; iset++){

		h2d_rho[iset] = (TH2D*)infile[iset]->Get("hmult_bin_mass_rho");
		h2d_kstar[iset] = (TH2D*)infile[iset]->Get("hmult_bin_mass_kstar");
		h2d_phi[iset] = (TH2D*)infile[iset]->Get("hmult_bin_mass_phi");

		for (int im=0; im<nmult; im++){

			h1d_rho[iset][im] = (TH1D*)h2d_rho[iset]->ProjectionY(Form("h1d_rho_set%d_%d",iset,im),im+1,im+1);
			h1d_kstar[iset][im] = (TH1D*)h2d_kstar[iset]->ProjectionY(Form("h1d_kstar_set%d_%d",iset,im),im+1,im+1);
			h1d_phi[iset][im] = (TH1D*)h2d_phi[iset]->ProjectionY(Form("h1d_phi_set%d_%d",iset,im),im+1,im+1);

		}//
	}//

	TCanvas *c1_rho = new TCanvas("c1_rho","c1_rho",1.2*3*300,3*300);
	c1_rho->Divide(3,3);

	for (int im=0; im<nmult; im++){

		c1_rho->cd(im+1);
		SetPadStyle();

		h1d_rho[0][im]->Scale(1./h1d_rho[0][im]->Integral());
		h1d_rho[1][im]->Scale(1./h1d_rho[1][im]->Integral());
		h1d_rho[0][im]->SetMarkerStyle(24);
		h1d_rho[0][im]->SetMarkerSize(0.6);
		h1d_rho[0][im]->SetLineColor(1);
		h1d_rho[1][im]->SetLineColor(4);

		htmp = (TH1D*)gPad->DrawFrame(0.3,0,1.5,1.1*h1d_rho[0][im]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",18,16);
		htmp->GetXaxis()->SetTitleOffset(2.5);

		h1d_rho[0][im]->Draw("same");
		h1d_rho[1][im]->Draw("same");

		{
			TLegend *leg = new TLegend(0.55,0.6,0.95,0.9);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.045);
			leg->AddEntry("","PYTHIA8 Angantyr","h");
			leg->AddEntry("","Pb-Pb 5 TeV","h");
			leg->AddEntry("","#rho#rightarrow#pi#pi","h");
			leg->AddEntry("",Form("Centrality %d%%-%d%%",cent_arr[im+1],cent_arr[im]),"h");
			leg->AddEntry(h1d_rho[0][im],"Rescattering off","L");
			leg->AddEntry(h1d_rho[1][im],"Rescattering on","L");
			leg->Draw();
		}//
	}//im

	TCanvas *c1_kstar = new TCanvas("c1_kstar","c1_kstar",1.2*3*300,3*300);
	c1_kstar->Divide(3,3);

	for (int im=0; im<nmult; im++){

		c1_kstar->cd(im+1);
		SetPadStyle();

		h1d_kstar[0][im]->Scale(1./h1d_kstar[0][im]->Integral());
		h1d_kstar[1][im]->Scale(1./h1d_kstar[1][im]->Integral());
		h1d_kstar[0][im]->SetMarkerStyle(24);
		h1d_kstar[0][im]->SetMarkerSize(0.6);
		h1d_kstar[0][im]->SetLineColor(1);
		h1d_kstar[1][im]->SetLineColor(4);

		htmp = (TH1D*)gPad->DrawFrame(0.64,0,1.2,1.1*h1d_kstar[0][im]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",18,16);
		htmp->GetXaxis()->SetTitleOffset(2.5);

		h1d_kstar[0][im]->Draw("same");
		h1d_kstar[1][im]->Draw("same");

		{
			TLegend *leg = new TLegend(0.55,0.6,0.95,0.9);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.045);
			leg->AddEntry("","PYTHIA8 Angantyr","h");
			leg->AddEntry("","Pb-Pb 5 TeV","h");
			leg->AddEntry("","K^{*0}#rightarrowK#pi","h");
			leg->AddEntry("",Form("Centrality %d%%-%d%%",cent_arr[im+1],cent_arr[im]),"h");
			leg->AddEntry(h1d_kstar[0][im],"Rescattering off","L");
			leg->AddEntry(h1d_kstar[1][im],"Rescattering on","L");
			leg->Draw();
		}//
	}//im

	TCanvas *c1_phi = new TCanvas("c1_phi","c1_phi",1.2*3*300,3*300);
	c1_phi->Divide(3,3);

	for (int im=0; im<nmult; im++){

		c1_phi->cd(im+1);
		SetPadStyle();

		h1d_phi[0][im]->Scale(1./h1d_phi[0][im]->Integral());
		h1d_phi[1][im]->Scale(1./h1d_phi[1][im]->Integral());
		h1d_phi[0][im]->SetMarkerStyle(24);
		h1d_phi[0][im]->SetMarkerSize(0.6);
		h1d_phi[0][im]->SetLineColor(1);
		h1d_phi[1][im]->SetLineColor(4);

		htmp = (TH1D*)gPad->DrawFrame(1.00,0,1.04,1.1*h1d_phi[0][im]->GetMaximum());
		SetHistoStyle("Mass (GeV/c^{2})","","",18,16);
		htmp->GetXaxis()->SetTitleOffset(2.5);

		h1d_phi[0][im]->Draw("same");
		h1d_phi[1][im]->Draw("same");

		{
			TLegend *leg = new TLegend(0.55,0.6,0.95,0.9);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.045);
			leg->AddEntry("","PYTHIA8 Angantyr","h");
			leg->AddEntry("","Pb-Pb 5 TeV","h");
			leg->AddEntry("","#phi#rightarrowKK","h");
			leg->AddEntry("",Form("Centrality %d%%-%d%%",cent_arr[im+1],cent_arr[im]),"h");
			leg->AddEntry(h1d_phi[0][im],"Rescattering off","P");
			leg->AddEntry(h1d_phi[1][im],"Rescattering on","L");
			leg->Draw();
		}//
	}//im

}

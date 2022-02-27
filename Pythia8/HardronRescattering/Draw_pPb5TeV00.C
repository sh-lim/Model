void Draw_pPb5TeV00(){

	const int nset = 2;
	const int nmult = 4;

	const int nColor[nmult] = {1, 2, 4, 6};

	TFile *infile[nset];
	infile[0] = new TFile("outfile_hist_pPb5TeV_set30_grp000_try101.root","read");
	infile[1] = new TFile("outfile_hist_pPb5TeV_set31_grp000_try101.root","read");

	TH1D *hpt_pion[nset][nmult];
	TH1D *hpt_kaon[nset][nmult];
	TH1D *hpt_rho[nset][nmult];
	TH1D *hpt_kstar[nset][nmult];
	TH1D *hpt_phi[nset][nmult];

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
	}//iset


	TH1D *hpt_pion[nset][nmult];
	TH1D *hpt_kstar_pion[nset][nmult];
	TH1D *hpt_phi_pion[nset][nmult];

	for (int iset=0; iset<nset; iset++){
		for (int im=0; im<nmult; im++){
			hpt_rho_pion[iset][im] = (TH1D*)hpt_rho[iset][im]->Clone(Form("hpt_rho_pion_%d_%d",iset,im));
			hpt_kstar_pion[iset][im] = (TH1D*)hpt_kstar[iset][im]->Clone(Form("hpt_kstar_pion_%d_%d",iset,im));
			hpt_phi_pion[iset][im] = (TH1D*)hpt_phi[iset][im]->Clone(Form("hpt_phi_pion_%d_%d",iset,im));

			hpt_rho_pion[iset][im]->Divide(hpt_pion[iset][im]);
			hpt_kstar_pion[iset][im]->Divide(hpt_pion[iset][im]);
			hpt_phi_pion[iset][im]->Divide(hpt_pion[iset][im]);
		}//im
	}//iset

	TCanvas *c0_rho = new TCanvas("c0_rho","c0_rho",1.2*2*500,500);
	for (int iset=0; iset<nset; iset++){
		c0_rho->cd(iset+1);
	}

	TCanvas *c0_kstar[nset];
	TCanvas *c0_phi[nset];



}

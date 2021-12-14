#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>
#include <TLorentzVector.h>


void SubJetPt(){

	const float const_pi = TMath::Pi();
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	float f_jet_pt[2000], f_jet_eta[2000], f_jet_phi[2000];

	float I_f_jet_pt[2000], I_f_jet_eta[2000], I_f_jet_phi[2000];

	float B_f_p_pt[2000], B_f_p_eta[2000], B_f_p_phi[2000];

	int i_np, i_njet;
	int i_p_id[2000], B_i_p_id[2000];
	int B_i_np, I_i_njet;
	float dR = 0.4;
	float area = const_pi*dR*dR;

	TH1F *njet = new TH1F("njet","",10,0,10);
	TH1F *subjetpt = new TH1F("subjetpt", "subjetpt", 100, 0, 200);
	//TH1F *jetptsub = new TH1F("jetptsub", "jetptsub", 100, 0, 200);
	TH1F *perpjetpt = new TH1F("perpjetpt", "perpjetpt", 100, 0, 200);
	TH1F *jetpt = new TH1F("jetpt", "jetpt", 100, 0, 200);
	TH1F *sigjetpt = new TH1F("sigjetpt", "sigjetpt", 100, 0, 200);
	TH1F *leadjet = new TH1F("leadjet", "leadjet", 100, 0, 200);
	//subjetpt->Sumw2(); jetpt->Sumw2(); jetptsub->Sumw2();

	TH1F *hdiffpt = new TH1F("hdiffpt","",100,-10,10);
	TH2D *hsigpt_subpt = new TH2D("hsigpt_subpt","",100,0,200,100,0,200);
	TH2D *hsigpt_incpt = new TH2D("hsigpt_incpt","",100,0,200,100,0,200);

	TFile *f = new TFile("/alice/data/shlim/PYTHIA/pp5TeV_set22_grp000/outfile_pp5TeV_set22_grp000_00000.root", "read");
	TTree *T = (TTree*)f ->Get("T");

	T->SetBranchAddress("IncjetR04_pt",I_f_jet_pt);
	T->SetBranchAddress("IncjetR04_eta",I_f_jet_eta);
	T->SetBranchAddress("IncjetR04_phi",I_f_jet_phi);
	T->SetBranchAddress("nIncjetR04",&I_i_njet);

	T->SetBranchAddress("Sig_p_pt",f_p_pt);
	T->SetBranchAddress("Sig_np",&i_np);
	T->SetBranchAddress("Sig_p_eta",f_p_eta);
	T->SetBranchAddress("Sig_p_phi",f_p_phi);
	T->SetBranchAddress("SigjetR04_pt",f_jet_pt);
	T->SetBranchAddress("SigjetR04_eta",f_jet_eta);
	T->SetBranchAddress("SigjetR04_phi",f_jet_phi);
	T->SetBranchAddress("nSigjetR04",&i_njet);
	T->SetBranchAddress("Sig_p_id",&i_p_id);

	T->SetBranchAddress("Bkg_p_pt",B_f_p_pt);
	T->SetBranchAddress("Bkg_np",&B_i_np);
	T->SetBranchAddress("Bkg_p_eta",B_f_p_eta);
	T->SetBranchAddress("Bkg_p_phi",B_f_p_phi);
	T->SetBranchAddress("Bkg_p_id",&B_i_p_id);


	int nent = T->GetEntries();
	cout << "nentries: " << nent << endl;

	for (int ien=0; ien<nent; ien++)
	{
		T->GetEntry(ien);

		if ( I_i_njet!=i_njet ) continue;

		njet->Fill(i_njet);
		//event selection
		int max = 0;
		if (I_f_jet_pt[max]<10) continue;
		leadjet->Fill(I_f_jet_pt[max]);

		//signal jets
		for (int i = 0; i < i_njet; i++)
		{
			if (fabs(f_jet_eta[i])>0.5) continue;
			//if (f_jet_pt[i]<10) continue;
			jetpt->Fill(f_jet_pt[i]);
		}

		//leading jet and rorate
		TVector3 vec_leadjet;
		vec_leadjet.SetPtEtaPhi(I_f_jet_pt[max], I_f_jet_eta[max], I_f_jet_phi[max]);
		TVector3 vOrtho(vec_leadjet);
		TVector3 perpjet;
		vOrtho.RotateZ(const_pi/2.);

		TVector3 Sigjetdummy;
		TVector3 Bkgjetdummy;
		std::vector<TVector3> Sigtracks;
		std::vector<TVector3> Bkgtracks;

		for (int itr = 0; itr < i_np; itr++)
		{
			if (fabs(f_p_eta[itr])>1) continue;
			if ( abs(i_p_id[itr])==12 || abs(i_p_id[itr])==14 || abs(i_p_id[itr])==16 ) continue;
			float deta = f_p_eta[itr] - vOrtho.Eta();
			float dphi = f_p_phi[itr] - vOrtho.Phi();

			if ( dphi < -const_pi ) dphi += 2*const_pi;
			else if ( dphi > const_pi ) dphi -= 2*const_pi;

			float newdR = sqrt(deta*deta + dphi*dphi);
			if(newdR>0.4) continue; 
			TVector3 sigtrack;
			sigtrack.SetPtEtaPhi(f_p_pt[itr], f_p_eta[itr], f_p_phi[itr]);
			Sigjetdummy += sigtrack;
			Sigtracks.push_back(sigtrack);
		}//itr  track loop
		for (int it = 0; it < B_i_np; it++)
		{
			if (fabs(B_f_p_eta[it])>1) continue;
			if ( abs(B_i_p_id[it])==12 || abs(B_i_p_id[it])==14 || abs(B_i_p_id[it])==16 ) continue;
			float deta = B_f_p_eta[it] - vOrtho.Eta();
			float dphi = B_f_p_phi[it] - vOrtho.Phi();

			if ( dphi < -const_pi ) dphi += 2*const_pi;
			else if ( dphi > const_pi ) dphi -= 2*const_pi;

			float newdR = sqrt(deta*deta + dphi*dphi);
			if(newdR>0.4) continue; 
			TVector3 bkgtrack;
			bkgtrack.SetPtEtaPhi(B_f_p_pt[it], B_f_p_eta[it], B_f_p_phi[it]);
			Bkgjetdummy += bkgtrack;
			Bkgtracks.push_back(bkgtrack);
		}//itr  track loop
		//????
		perpjet.SetPtEtaPhi(Sigjetdummy.Pt()+Bkgjetdummy.Pt(),Sigjetdummy.Eta()+Bkgjetdummy.Eta(),Sigjetdummy.Phi()+Bkgjetdummy.Phi());

		float totpt = 0;
		float totalpt = 0;
		float Sigtotpt = 0.;
		for (int j = 0; j < Sigtracks.size(); j++)
		{
			Sigtotpt += Sigtracks[j].Pt()/area;
		}
		float Bkgtotpt = 0.;
		for (int k = 0; k < Bkgtracks.size(); k++)
		{
			Bkgtotpt += Bkgtracks[k].Pt()/area;
		}

		totpt = Sigtotpt+Bkgtotpt;
		totalpt += perpjet.Pt();

		perpjetpt->Fill(totpt*area);

		for (int jet = 0; jet < i_njet; jet++)
		{
			if (fabs(I_f_jet_eta[jet])>0.5) continue;
			float subjet = 0;
			subjet = I_f_jet_pt[jet]-totpt*area;
			//if (subjet<10) continue;
			subjetpt->Fill(subjet);
			//jetptsub->Fill(I_f_jet_pt[jet]-totalpt);
		}

		if ( !(I_i_njet==1 && i_njet==1) ) continue;
		if ( fabs(I_f_jet_eta[0])>0.5 ) continue;
		if ( I_f_jet_pt[0]<10.0 ) continue;

		hsigpt_incpt->Fill(f_jet_pt[0], I_f_jet_pt[0]);
		hsigpt_subpt->Fill(f_jet_pt[0], I_f_jet_pt[0]-totpt*area);
		hdiffpt->Fill((I_f_jet_pt[0]-totpt*area) - f_jet_pt[0]);

	}//ien


	TFile *outfile = new TFile("ptcorr.root", "RECREATE");


	subjetpt->Write();
	njet->Write();
	jetpt->Write();
	perpjetpt->Write();
	leadjet->Write();
	//jetptsub->Write();
	//sigjetpt->Write();
	//hleadvsbkg->Write();
	//hmeanbkg->Write();

	hsigpt_incpt->Write();
	hsigpt_subpt->Write();
	hdiffpt->Write();

	// subjt->Write();
	// subz->Write();


	outfile->Close();

}//void

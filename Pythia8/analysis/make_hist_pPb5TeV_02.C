//Hadron rescattering effect
//rho, Kstar, phi

#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace std;

//void make_hist_pPb5TeV_02(int set = 0, int grp = 0){
void make_hist_pPb5TeV_02(int set = 30, const char *listname = "file.lst"){

	const float const_pi = TMath::Pi();

	const float mass_pion = 0.13957; 
	const float mass_kaon = 0.49367; 

	const int nmult = 4;
	float multcut[nmult] = {28, 45, 68, 500};
	if ( set==30 ){
		multcut[0] = 27;
		multcut[1] = 42;
		multcut[2] = 63;
		multcut[3] = 500;
	}

	const int npt = 15;
	const float ptbin[npt+1] = {0.00, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 7.00, 8.00, 10.0};

	char fname[500];
	//sprintf(fname, "file_pPb5TeV_set%d_grp%03d.lst", set, grp);

	ifstream flist;
	flist.open(listname);

	char ffname[300];

	int i_np;
	int i_nmult_v0a;
	int i_p_id[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	float f_p_vt[2000];

	int i_p_mom1_id[2000];
	int i_p_mom2_id[2000];
	int i_p_mom1_index[2000];
	float f_p_mom1_mass[2000];

	TH1D *hmult_v0a = new TH1D("hmult_v0a","",500,0,500);

	TH1D *hpt_pion[nmult];
	TH1D *hpt_kaon[nmult];
	TH1D *hpt_rho[nmult];
	TH1D *hpt_kstar[nmult];
	TH1D *hpt_phi[nmult];

	for (int im=0; im<nmult; im++){
		hpt_pion[im] = new TH1D(Form("hpt_pion_m%d",im),"",npt,ptbin);
		hpt_kaon[im] = new TH1D(Form("hpt_kaon_m%d",im),"",npt,ptbin);

		hpt_rho[im] = new TH1D(Form("hpt_rho_m%d",im),"",npt,ptbin);
		hpt_kstar[im] = new TH1D(Form("hpt_kstar_m%d",im),"",npt,ptbin);
		hpt_phi[im] = new TH1D(Form("hpt_phi_m%d",im),"",npt,ptbin);
	}

	TH1D *hmass_truth_rho = new TH1D("hmass_truth_rho","",100,0.30,1.50);
	TH1D *hmass_truth_kstar = new TH1D("hmass_truth_kstar","",100,0.60,1.20);
	TH1D *hmass_truth_phi = new TH1D("hmass_truth_phi","",100,1.00,1.04);

	TH1D *hmass_rho = new TH1D("hmass_rho","",100,0.30,1.50);
	TH1D *hmass_kstar = new TH1D("hmass_kstar","",100,0.60,1.20);
	TH1D *hmass_phi = new TH1D("hmass_phi","",100,1.00,1.04);

	vector<int> rho_index;
	vector<int> rho_mom_index;
	vector<int> phi_index;
	vector<int> phi_mom_index;
	vector<int> kstar_index;
	vector<int> kstar_mom_index;
	

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("nmult_v0a",&i_nmult_v0a);

		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id",i_p_id);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_pt",f_p_pt);
		T->SetBranchAddress("p_vt",f_p_vt);

		T->SetBranchAddress("p_mom1_id",i_p_mom1_id);
		T->SetBranchAddress("p_mom2_id",i_p_mom2_id);
		T->SetBranchAddress("p_mom1_index",i_p_mom1_index);
		T->SetBranchAddress("p_mom1_mass",f_p_mom1_mass);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			rho_index.clear();
			rho_mom_index.clear();
			kstar_index.clear();
			kstar_mom_index.clear();
			phi_index.clear();
			phi_mom_index.clear();

			//multiplicity
			hmult_v0a->Fill(i_nmult_v0a);

			//continue;

			int mult_index = -1;

			if ( i_nmult_v0a<multcut[0] ) mult_index = 0;
			else if ( i_nmult_v0a<multcut[1] ) mult_index = 1;
			else if ( i_nmult_v0a<multcut[2] ) mult_index = 2;
			else mult_index = 3;

			//charged pion
			for (int ip=0; ip<i_np; ip++){
				if ( abs(i_p_id[ip])!=211 ) continue;

				TLorentzVector lvec;
				lvec.SetPtEtaPhiM(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], mass_pion);

				if ( fabs(lvec.Rapidity())>0.5 ) continue;

				hpt_pion[mult_index]->Fill(f_p_pt[ip]);
			}

			//charged kaon
			for (int ip=0; ip<i_np; ip++){
				if ( abs(i_p_id[ip])!=321 ) continue;

				TLorentzVector lvec;
				lvec.SetPtEtaPhiM(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], mass_kaon);

				if ( fabs(lvec.Rapidity())>0.5 ) continue;

				hpt_kaon[mult_index]->Fill(f_p_pt[ip]);
			}

			//rho
			for (int ip=0; ip<i_np; ip++){
				if ( abs(i_p_id[ip])!=211 ) continue;
				if ( i_p_mom1_id[ip]!=113 ) continue;
				//if ( i_p_mom2_id[ip]!=90 ) continue;

				rho_index.push_back(ip);
				rho_mom_index.push_back(i_p_mom1_index[ip]);
			}

			int nd_rho = rho_index.size();

			for (int ii=0; ii<nd_rho-1; ii++){

				if ( rho_mom_index[ii]!=rho_mom_index[ii+1] ) continue;

				TLorentzVector lvec_ii, lvec_jj;
				lvec_ii.SetPtEtaPhiM(f_p_pt[rho_index[ii]], f_p_eta[rho_index[ii]], f_p_phi[rho_index[ii]], mass_pion);
				lvec_jj.SetPtEtaPhiM(f_p_pt[rho_index[ii+1]], f_p_eta[rho_index[ii+1]], f_p_phi[rho_index[ii+1]], mass_pion);

				TLorentzVector lvec_rho = lvec_ii + lvec_jj;

				//cout << lvec_rho.Pt() << " " << lvec_rho.Eta() << " " << lvec_rho.M() << endl;

				hmass_truth_rho->Fill(f_p_mom1_mass[rho_index[ii]]);
				hmass_rho->Fill(lvec_rho.M());

				if ( fabs(lvec_rho.Rapidity())>0.5 ) continue;

				hpt_rho[mult_index]->Fill(lvec_rho.Pt());
			}

			//continue;

			//k-star
			for (int ip=0; ip<i_np; ip++){
				if ( abs(i_p_id[ip])!=321 && abs(i_p_id[ip])!=211 ) continue;
				if ( i_p_mom1_id[ip]!=313 ) continue;
				//if ( i_p_mom2_id[ip]!=90 ) continue;

				kstar_index.push_back(ip);
				kstar_mom_index.push_back(i_p_mom1_index[ip]);
			}

			int nd_kstar = kstar_index.size();

			for (int ii=0; ii<nd_kstar-1; ii++){

				if ( kstar_mom_index[ii]!=kstar_mom_index[ii+1] ) continue;

				TLorentzVector lvec_ii, lvec_jj;
				if ( abs(i_p_id[kstar_index[ii]])==211 ){
					lvec_ii.SetPtEtaPhiM(f_p_pt[kstar_index[ii]], f_p_eta[kstar_index[ii]], f_p_phi[kstar_index[ii]], mass_pion);
					lvec_jj.SetPtEtaPhiM(f_p_pt[kstar_index[ii+1]], f_p_eta[kstar_index[ii+1]], f_p_phi[kstar_index[ii+1]], mass_kaon);
				}else{
					lvec_ii.SetPtEtaPhiM(f_p_pt[kstar_index[ii]], f_p_eta[kstar_index[ii]], f_p_phi[kstar_index[ii]], mass_kaon);
					lvec_jj.SetPtEtaPhiM(f_p_pt[kstar_index[ii+1]], f_p_eta[kstar_index[ii+1]], f_p_phi[kstar_index[ii+1]], mass_pion);
				}

				TLorentzVector lvec_kstar = lvec_ii + lvec_jj;

				//cout << lvec_kstar.Pt() << " " << lvec_kstar.Eta() << " " << lvec_kstar.M() << endl;

				hmass_truth_kstar->Fill(f_p_mom1_mass[kstar_index[ii]]);
				hmass_kstar->Fill(lvec_kstar.M());

				if ( fabs(lvec_kstar.Rapidity())>0.5 ) continue;

				hpt_kstar[mult_index]->Fill(lvec_kstar.Pt());
			}

			//phi
			for (int ip=0; ip<i_np; ip++){
				if ( abs(i_p_id[ip])!=321 ) continue;
				if ( i_p_mom1_id[ip]!=333 ) continue;
				//if ( i_p_mom2_id[ip]!=90 ) continue;

				phi_index.push_back(ip);
				phi_mom_index.push_back(i_p_mom1_index[ip]);
			}

			int nd_phi = phi_index.size();

			for (int ii=0; ii<nd_phi-1; ii++){

				if ( phi_mom_index[ii]!=phi_mom_index[ii+1] ) continue;

				TLorentzVector lvec_ii, lvec_jj;
				lvec_ii.SetPtEtaPhiM(f_p_pt[phi_index[ii]], f_p_eta[phi_index[ii]], f_p_phi[phi_index[ii]], mass_kaon);
				lvec_jj.SetPtEtaPhiM(f_p_pt[phi_index[ii+1]], f_p_eta[phi_index[ii+1]], f_p_phi[phi_index[ii+1]], mass_kaon);

				TLorentzVector lvec_phi = lvec_ii + lvec_jj;

				//cout << lvec_phi.Pt() << " " << lvec_phi.Eta() << " " << lvec_phi.M() << endl;

				hmass_phi->Fill(lvec_phi.M());

				if ( fabs(lvec_phi.Rapidity())>0.5 ) continue;

				hpt_phi[mult_index]->Fill(lvec_phi.Pt());
			}

		}//ien

		delete infile;

		/*
		*/

	}//

	//TFile *outfile = new TFile(Form("outfile_hist_set%d_grp%03d.root",set,grp),"recreate");
	TFile *outfile = new TFile("outfile_hist.root","recreate");

	hmult_v0a->Write();

	for (int im=0; im<nmult; im++){
		hpt_pion[im]->Write();
		hpt_kaon[im]->Write();

		hpt_rho[im]->Write();
		hpt_kstar[im]->Write();
		hpt_phi[im]->Write();
	}

	hmass_rho->Write();
	hmass_kstar->Write();
	hmass_phi->Write();

	hmass_truth_rho->Write();
	hmass_truth_kstar->Write();
	hmass_truth_phi->Write();

	/*
	hLc_pt_v0m->Write();
	hD0_pt_v0m->Write();

	htrk_eta_pt->Write();

	hjetR04_eta_pt->Write();
	hjetR04_trk_dphi_deta->Write();

	hjetR04_jt_pt->Write();
	hjetR04_zh_pt->Write();
	*/


	outfile->Close();

}

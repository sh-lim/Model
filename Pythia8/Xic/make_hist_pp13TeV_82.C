//Xic - e, Xi pair analysis

#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace std;

void make_hist_pp13TeV_82(int set = 0, int grp = 0){

	const float const_pi = TMath::Pi();

	const int nmult = 4;
	//const float cut_mult[nmult+1] = {100, 60, 20, 5, 1, 0};

	//const int npt = 5;
	//const float ptbin[npt+1] = {0.2, 1.0, 2.0, 3.0, 4.0, 5.0};

	const int npt = 7;
	const float ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TH1D *hpt = new TH1D("hpt","",npt,ptbin);

	char fname[500];
	sprintf(fname, "file_set%d_grp%03d.lst", set, grp);

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int i_np;
	int i_multV0M;
	int i_p_id[100];
	float f_p_pt[100], f_p_eta[100], f_p_phi[100], f_p_mass[100];

	int i_p1_id[100][30], i_p1_index[100][30], i_p1_status[100][30];

	TH1D *hmult_v0m = new TH1D("hmult_v0m","",500,0,500);

	TH2D *hmass_pt_sig3Xic0 = new TH2D("hmass_pt_sig3Xic0","",100,0,5,20,0,20);
	TH2D *hmass_pt_sig4Xic0 = new TH2D("hmass_pt_sig4Xic0","",100,0,5,20,0,20);
	TH2D *hmass_pt_sig4Xicp = new TH2D("hmass_pt_sig4Xicp","",100,0,5,20,0,20);

	TH2D *hmass_pt_cut_sig3Xic0 = new TH2D("hmass_pt_cut_sig3Xic0","",100,0,5,20,0,20);
	TH2D *hmass_pt_cut_sig4Xic0 = new TH2D("hmass_pt_cut_sig4Xic0","",100,0,5,20,0,20);
	TH2D *hmass_pt_cut_sig4Xicp = new TH2D("hmass_pt_cut_sig4Xicp","",100,0,5,20,0,20);

	TH2D *hmass_pt_cut2_sig3Xic0 = new TH2D("hmass_pt_cut2_sig3Xic0","",100,0,5,20,0,20);
	TH2D *hmass_pt_cut2_sig4Xic0 = new TH2D("hmass_pt_cut2_sig4Xic0","",100,0,5,20,0,20);
	TH2D *hmass_pt_cut2_sig4Xicp = new TH2D("hmass_pt_cut2_sig4Xicp","",100,0,5,20,0,20);

	TH1D *hphi_e = new TH1D("hphi_e","",90,-const_pi,const_pi);
	TH1D *hphi_Xi = new TH1D("hphi_Xi","",90,-const_pi,const_pi);

	TH2D *hmass_pt_corrJet[2];
	TH2D *hmass_pt_corrXib0[2];
	TH2D *hmass_pt_corrXibm[2];
	TH2D *hmass_pt_corrLambdab0[2];
	TH2D *hmass_pt_corrOmegabm[2];
	TH2D *hmass_pt_corrOmegac0[2];
	TH2D *hmass_pt_corrXiccp[2];
	TH2D *hmass_pt_corrB0[2];
	TH2D *hmass_pt_corrDiquark[2];
	TH2D *hmass_pt_comb[2];
	TH2D *hmass_pt_mixed1[2];
	TH2D *hmass_pt_mixed2[2];

	TH2D *hoa_pt_comb[2];
	TH2D *hoa_pt_mixed1[2];
	TH2D *hoa_pt_mixed2[2];
	TH2D *hoa_pt_sig3Xic0[2];
	TH2D *hoa_pt_sig4Xic0[2];
	TH2D *hoa_pt_sig4Xicp[2];
	TH2D *hoa_pt_corrJet[2];
	TH2D *hoa_pt_corrDiquark[2];
	TH2D *hoa_pt_corrXibm[2];

	TH2D *hoa_mass_comb = new TH2D("hoa_mass_comb","",100,0,5,90,0,180);
	TH2D *hoa_mass_corrJet = new TH2D("hoa_mass_corrJet","",100,0,5,90,0,180);
	TH2D *hoa_mass_sig3Xic0 = new TH2D("hoa_mass_sig3Xic0","",100,0,5,90,0,180);

	TH2D *hmass_pt_cut_corrJet[2];
	TH2D *hmass_pt_cut_corrXib0[2];
	TH2D *hmass_pt_cut_corrXibm[2];
	TH2D *hmass_pt_cut_corrLambdab0[2];
	TH2D *hmass_pt_cut_corrOmegabm[2];
	TH2D *hmass_pt_cut_corrOmegac0[2];
	TH2D *hmass_pt_cut_corrXiccp[2];
	TH2D *hmass_pt_cut_corrB0[2];
	TH2D *hmass_pt_cut_corrDiquark[2];
	TH2D *hmass_pt_cut_comb[2];
	TH2D *hmass_pt_cut_mixed1[2];
	TH2D *hmass_pt_cut_mixed2[2];

	TH2D *hmass_pt_cut2_corrJet[2];
	TH2D *hmass_pt_cut2_corrXib0[2];
	TH2D *hmass_pt_cut2_corrXibm[2];
	TH2D *hmass_pt_cut2_corrLambdab0[2];
	TH2D *hmass_pt_cut2_corrOmegabm[2];
	TH2D *hmass_pt_cut2_corrOmegac0[2];
	TH2D *hmass_pt_cut2_corrXiccp[2];
	TH2D *hmass_pt_cut2_corrB0[2];
	TH2D *hmass_pt_cut2_corrDiquark[2];
	TH2D *hmass_pt_cut2_comb[2];
	TH2D *hmass_pt_cut2_mixed1[2];
	TH2D *hmass_pt_cut2_mixed2[2];

	for (int ii=0; ii<2; ii++){
		hmass_pt_corrJet[ii] = new TH2D(Form("hmass_pt_corrJet_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrXibm[ii] = new TH2D(Form("hmass_pt_corrXibm_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrXib0[ii] = new TH2D(Form("hmass_pt_corrXib0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrLambdab0[ii] = new TH2D(Form("hmass_pt_corrLambda0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrOmegabm[ii] = new TH2D(Form("hmass_pt_corrOmegabm_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrOmegac0[ii] = new TH2D(Form("hmass_pt_corrOmegac0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrXiccp[ii] = new TH2D(Form("hmass_pt_corrXiccp_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrB0[ii] = new TH2D(Form("hmass_pt_corrB0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_corrDiquark[ii] = new TH2D(Form("hmass_pt_corrDiquark_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_comb[ii] = new TH2D(Form("hmass_pt_comb_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_mixed1[ii] = new TH2D(Form("hmass_pt_mixed1_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_mixed2[ii] = new TH2D(Form("hmass_pt_mixed2_chg%d",ii),"",100,0,5,20,0,20);

		hoa_pt_comb[ii] = new TH2D(Form("hoa_pt_comb_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_mixed1[ii] = new TH2D(Form("hoa_pt_mixed1_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_mixed2[ii] = new TH2D(Form("hoa_pt_mixed2_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_sig3Xic0[ii] = new TH2D(Form("hoa_pt_sig3Xic0_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_sig4Xic0[ii] = new TH2D(Form("hoa_pt_sig4Xic0_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_sig4Xicp[ii] = new TH2D(Form("hoa_pt_sig4Xicp_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_corrJet[ii] = new TH2D(Form("hoa_pt_corrJet_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_corrDiquark[ii] = new TH2D(Form("hoa_pt_corrDiquark_chg%d",ii),"",20,0,20,90,0,180);
		hoa_pt_corrXibm[ii] = new TH2D(Form("hoa_pt_corrXibm_chg%d",ii),"",20,0,20,90,0,180);

		hmass_pt_cut_corrJet[ii] = new TH2D(Form("hmass_pt_cut_corrJet_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrXibm[ii] = new TH2D(Form("hmass_pt_cut_corrXibm_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrXib0[ii] = new TH2D(Form("hmass_pt_cut_corrXib0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrLambdab0[ii] = new TH2D(Form("hmass_pt_cut_corrLambda0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrOmegabm[ii] = new TH2D(Form("hmass_pt_cut_corrOmegabm_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrOmegac0[ii] = new TH2D(Form("hmass_pt_cut_corrOmegac0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrXiccp[ii] = new TH2D(Form("hmass_pt_cut_corrXiccp_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrB0[ii] = new TH2D(Form("hmass_pt_cut_corrB0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_corrDiquark[ii] = new TH2D(Form("hmass_pt_cut_corrDiquark_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_comb[ii] = new TH2D(Form("hmass_pt_cut_comb_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_mixed1[ii] = new TH2D(Form("hmass_pt_cut_mixed1_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut_mixed2[ii] = new TH2D(Form("hmass_pt_cut_mixed2_chg%d",ii),"",100,0,5,20,0,20);

		hmass_pt_cut2_corrJet[ii] = new TH2D(Form("hmass_pt_cut2_corrJet_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrXibm[ii] = new TH2D(Form("hmass_pt_cut2_corrXibm_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrXib0[ii] = new TH2D(Form("hmass_pt_cut2_corrXib0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrLambdab0[ii] = new TH2D(Form("hmass_pt_cut2_corrLambda0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrOmegabm[ii] = new TH2D(Form("hmass_pt_cut2_corrOmegabm_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrOmegac0[ii] = new TH2D(Form("hmass_pt_cut2_corrOmegac0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrXiccp[ii] = new TH2D(Form("hmass_pt_cut2_corrXiccp_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrB0[ii] = new TH2D(Form("hmass_pt_cut2_corrB0_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_corrDiquark[ii] = new TH2D(Form("hmass_pt_cut2_corrDiquark_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_comb[ii] = new TH2D(Form("hmass_pt_cut2_comb_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_mixed1[ii] = new TH2D(Form("hmass_pt_cut2_mixed1_chg%d",ii),"",100,0,5,20,0,20);
		hmass_pt_cut2_mixed2[ii] = new TH2D(Form("hmass_pt_cut2_mixed2_chg%d",ii),"",100,0,5,20,0,20);
	}

	vector<float> vec_e_pt;
	vector<float> vec_e_eta;
	vector<float> vec_e_phi;
	vector<float> vec_e_mass;
	vector<int> vec_e_id;

	vector<float> vec_Xi_pt;
	vector<float> vec_Xi_eta;
	vector<float> vec_Xi_phi;
	vector<float> vec_Xi_mass;
	vector<int> vec_Xi_id;


	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;
		TFile *infile = new TFile(ffname,"read");

		TH1D *_hmult_v0m = (TH1D*)infile->Get("hmult_v0m");
		hmult_v0m->Add(_hmult_v0m);

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("multV0M",&i_multV0M);
		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id",i_p_id);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_pt",f_p_pt);
		T->SetBranchAddress("p_mass",f_p_mass);

		//T->SetBranchAddress("p1_id",i_p1_id);
		//T->SetBranchAddress("p1_index",i_p1_index);
		T->SetBranchAddress("p1_index",i_p1_id);
		T->SetBranchAddress("p1_id",i_p1_index);
		T->SetBranchAddress("p1_status",i_p1_status);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			//particle info
			for (int ip=0; ip<i_np; ip++){

				if ( fabs(f_p_eta[ip])>0.8 ) continue;
				if ( abs(i_p_id[ip])==11 && fabs(f_p_pt[ip])<0.5 ) continue;

				TLorentzVector ivec;
				ivec.SetPtEtaPhiM(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], f_p_mass[ip]);

				if ( abs(i_p_id[ip])==11 ){
					hphi_e->Fill(ivec.Phi());
				}else{
					hphi_Xi->Fill(ivec.Phi());
				}

				for (int jp=ip+1; jp<i_np; jp++){

					if ( fabs(f_p_eta[jp])>0.8 ) continue;
					if ( abs(i_p_id[jp])==11 && fabs(f_p_pt[jp])<0.5 ) continue;
					if ( abs(i_p_id[ip])==abs(i_p_id[jp]) ) continue;

					//11:e-, 3312:Xi-
					int ind_sign = (i_p_id[ip]*i_p_id[jp]>0) ? 1 : 0;

					TLorentzVector jvec;
					jvec.SetPtEtaPhiM(f_p_pt[jp], f_p_eta[jp], f_p_phi[jp], f_p_mass[jp]);

					TLorentzVector pvec = ivec + jvec;

					TLorentzVector evec;
					TLorentzVector Xivec;

					if ( abs(i_p_id[ip])==11 ){
						evec = ivec;
						Xivec = jvec;
					}else{
						Xivec = ivec;
						evec = jvec;
					}

					float cosoa = (evec.Px()*Xivec.Px() + evec.Py()*Xivec.Py() + evec.Pz()*Xivec.Pz())/evec.P()/Xivec.P();
					float oa = acos(cosoa)*180.0/const_pi;
					//if ( oa>90.0 ) bgood = false;
					//

					int _pt_ind = hpt->FindBin(pvec.Pt());

					int _mom_id = -1;
					int _matched_ii = -1;
					int _matched_jj = -1;

					for (int ii=0; ii<30; ii++){
						if ( i_p1_id[ip][ii]==0 ) continue;
						if ( abs(i_p1_status[ip][ii])<20 ){
							continue;
						}else if ( abs(i_p1_status[ip][ii])<30 ){
							if ( !(abs(i_p1_status[ip][ii])==23 ||  abs(i_p1_status[ip][ii])==24) ) continue;
						}else if ( abs(i_p1_status[ip][ii])<40 ){
							if ( !(abs(i_p1_status[ip][ii])==33) ) continue;
						}else if ( abs(i_p1_status[ip][ii])<50 ){
							if ( !(abs(i_p1_status[ip][ii])==43 ||  abs(i_p1_status[ip][ii])==44) ) continue;
						}

						for (int jj=0; jj<30; jj++){
							if ( i_p1_id[jp][jj]==0 ) continue;
							if ( abs(i_p1_status[jp][jj])<20 ){
								continue;
							}else if ( abs(i_p1_status[jp][jj])<30 ){
								if ( !(abs(i_p1_status[jp][jj])==23 ||  abs(i_p1_status[jp][jj])==24) ) continue;
							}else if ( abs(i_p1_status[jp][jj])<40 ){
								if ( !(abs(i_p1_status[jp][jj])==33) ) continue;
							}else if ( abs(i_p1_status[jp][jj])<50 ){
								if ( !(abs(i_p1_status[jp][jj])==43 ||  abs(i_p1_status[jp][jj])==44) ) continue;
							}

							if ( i_p1_index[ip][ii]==i_p1_index[jp][jj] ){
								_mom_id = i_p1_id[ip][ii];
								_matched_ii = ii;
								_matched_jj = jj;
								//cout << "matched: " << ii << " " << jj << " " << i_p1_id[ip][ii] << " " << i_p1_id[jp][jj] << " " << i_p1_index[ip][ii] << endl;
								break;
							}

						}//jj

						if ( _matched_ii>-1 && _matched_jj>-1 ){
							break;
						}
					}//ii

					if ( _matched_ii>-1 && _matched_jj>-1 ){

						if ( _matched_ii==0 && _matched_jj==0 && abs(_mom_id)==4132 ){
							hmass_pt_sig3Xic0->Fill(pvec.M(), pvec.Pt());
							hoa_pt_sig3Xic0[ind_sign]->Fill(pvec.Pt(), oa);
							if ( pvec.Pt()>2 && pvec.Pt()<4 ){
								hoa_mass_sig3Xic0->Fill(pvec.M(), oa);
							}
							if ( oa<90.0 ){
								hmass_pt_cut_sig3Xic0->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_sig3Xic0->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==4132 ){
							hmass_pt_sig4Xic0->Fill(pvec.M(), pvec.Pt());
							hoa_pt_sig4Xic0[ind_sign]->Fill(pvec.Pt(), oa);
							if ( oa<90.0 ){
								hmass_pt_cut_sig4Xic0->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_sig4Xic0->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==4232 ){
							hmass_pt_sig4Xicp->Fill(pvec.M(), pvec.Pt());
							hoa_pt_sig4Xicp[ind_sign]->Fill(pvec.Pt(), oa);
							if ( oa<90.0 ){
								hmass_pt_cut_sig4Xicp->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_sig4Xicp->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==4332 ){
							hmass_pt_corrOmegac0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							if ( oa<90.0 ){
								hmass_pt_cut_corrOmegac0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrOmegac0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==4412 || abs(_mom_id)==4414 ){
							hmass_pt_corrXiccp[ind_sign]->Fill(pvec.M(), pvec.Pt());
							if ( oa<90.0 ){
								hmass_pt_cut_corrXiccp[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrXiccp[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==511 ){
							hmass_pt_corrB0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							if ( oa<90.0 ){
								hmass_pt_cut_corrB0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrB0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==5122 ){
							hmass_pt_corrLambdab0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							if ( oa<90.0 ){
								hmass_pt_cut_corrLambdab0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrLambdab0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==5132 ){
							hmass_pt_corrXibm[ind_sign]->Fill(pvec.M(), pvec.Pt());
							hoa_pt_corrXibm[ind_sign]->Fill(pvec.Pt(), oa);
							if ( oa<90.0 ){
								hmass_pt_cut_corrXibm[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrXibm[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==5232 ){
							hmass_pt_corrXib0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							if ( oa<90.0 ){
								hmass_pt_cut_corrXib0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrXib0[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)==5332 ){
							hmass_pt_corrOmegabm[ind_sign]->Fill(pvec.M(), pvec.Pt());
							if ( oa<90.0 ){
								hmass_pt_cut_corrOmegabm[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrOmegabm[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else if ( abs(_mom_id)<100 ){

							if ( 0 ){
								cout << "corrJet" << endl;
								cout << i_p_id[ip] << " " 
									<< i_p1_id[ip][0] << "(" << i_p1_index[ip][0] << ") " 
									<< i_p1_id[ip][1] << "(" << i_p1_index[ip][1] << ") " 
									<< i_p1_id[ip][2] << "(" << i_p1_index[ip][2] << ") " 
									<< i_p1_id[ip][3] << "(" << i_p1_index[ip][3] << ") " 
									<< i_p1_id[ip][4] << "(" << i_p1_index[ip][4] << ") " 
									<< endl;
								cout << i_p_id[jp] << " " 
									<< i_p1_id[jp][0] << "(" << i_p1_index[jp][0] << ") " 
									<< i_p1_id[jp][1] << "(" << i_p1_index[jp][1] << ") " 
									<< i_p1_id[jp][2] << "(" << i_p1_index[jp][2] << ") " 
									<< i_p1_id[jp][3] << "(" << i_p1_index[jp][3] << ") " 
									<< i_p1_id[jp][4] << "(" << i_p1_index[jp][4] << ") " 
									<< endl;
								cout << "#######" << endl;
							}

							hmass_pt_corrJet[ind_sign]->Fill(pvec.M(), pvec.Pt());
							hoa_pt_corrJet[ind_sign]->Fill(pvec.Pt(), oa);
							if ( pvec.Pt()>2 && pvec.Pt()<4 ){
								hoa_mass_corrJet->Fill(pvec.M(), oa);
							}
							if ( oa<90.0 ){
								hmass_pt_cut_corrJet[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
							if ( oa<100.0-10*_pt_ind ){
								hmass_pt_cut2_corrJet[ind_sign]->Fill(pvec.M(), pvec.Pt());
							}
						}else{
							int _digit10 = ((_mom_id)%100)/10;
							if ( _digit10==0 ){
								hmass_pt_corrDiquark[ind_sign]->Fill(pvec.M(), pvec.Pt());
								hoa_pt_corrDiquark[ind_sign]->Fill(pvec.Pt(), oa);
								if ( oa<90.0 ){
									hmass_pt_cut_corrDiquark[ind_sign]->Fill(pvec.M(), pvec.Pt());
								}
								if ( oa<100.0-10*_pt_ind ){
									hmass_pt_cut2_corrDiquark[ind_sign]->Fill(pvec.M(), pvec.Pt());
								}
							}else{
								cout << _mom_id << " " << i_p_id[ip] << " " << i_p_id[jp] << endl;
							}
						}

					}else{

						hmass_pt_comb[ind_sign]->Fill(pvec.M(), pvec.Pt());
						hoa_pt_comb[ind_sign]->Fill(pvec.Pt(), oa);

						if ( pvec.Pt()>2 && pvec.Pt()<4 ){
							hoa_mass_comb->Fill(pvec.M(), oa);
						}
						if ( oa<90.0 ){
							hmass_pt_cut_comb[ind_sign]->Fill(pvec.M(), pvec.Pt());
						}
						if ( oa<100.0-10*_pt_ind ){
							hmass_pt_cut2_comb[ind_sign]->Fill(pvec.M(), pvec.Pt());
						}

						if ( 0 ){
							cout << "combinatorial" << endl;
							cout << i_p_id[ip] << " " 
								<< i_p1_id[ip][0] << "(" << i_p1_index[ip][0] << ") " 
								<< i_p1_id[ip][1] << "(" << i_p1_index[ip][1] << ") " 
								<< i_p1_id[ip][2] << "(" << i_p1_index[ip][2] << ") " 
								<< i_p1_id[ip][3] << "(" << i_p1_index[ip][3] << ") " 
								<< i_p1_id[ip][4] << "(" << i_p1_index[ip][4] << ") " 
								<< i_p1_id[ip][5] << "(" << i_p1_index[ip][5] << ") " 
								<< i_p1_id[ip][6] << "(" << i_p1_index[ip][6] << ") " 
								<< i_p1_id[ip][7] << "(" << i_p1_index[ip][7] << ") " 
								<< i_p1_id[ip][8] << "(" << i_p1_index[ip][8] << ") " 
								<< i_p1_id[ip][9] << "(" << i_p1_index[ip][9] << ") " 
								<< endl;
							cout << i_p_id[jp] << " " 
								<< i_p1_id[jp][0] << "(" << i_p1_index[jp][0] << ") " 
								<< i_p1_id[jp][1] << "(" << i_p1_index[jp][1] << ") " 
								<< i_p1_id[jp][2] << "(" << i_p1_index[jp][2] << ") " 
								<< i_p1_id[jp][3] << "(" << i_p1_index[jp][3] << ") " 
								<< i_p1_id[jp][4] << "(" << i_p1_index[jp][4] << ") " 
								<< i_p1_id[jp][5] << "(" << i_p1_index[jp][5] << ") " 
								<< i_p1_id[jp][6] << "(" << i_p1_index[jp][6] << ") " 
								<< i_p1_id[jp][7] << "(" << i_p1_index[jp][7] << ") " 
								<< i_p1_id[jp][8] << "(" << i_p1_index[jp][8] << ") " 
								<< i_p1_id[jp][9] << "(" << i_p1_index[jp][9] << ") " 
								<< endl;
							cout << "#######" << endl;
						}

					}//matched

				}//jp

				//mixed event
				if ( abs(i_p_id[ip])==3312 ){
					for (unsigned int kp=0; kp<vec_e_pt.size(); kp++){

						//11:e-, 3312:Xi-
						int ind_sign = (i_p_id[ip]*vec_e_id[kp]>0) ? 1 : 0;

						TLorentzVector kvec;
						kvec.SetPtEtaPhiM(vec_e_pt[kp], vec_e_eta[kp], vec_e_phi[kp], vec_e_mass[kp]);

						TLorentzVector pvec = ivec + kvec;

						hmass_pt_mixed1[ind_sign]->Fill(pvec.M(), pvec.Pt());

						float cosoa = (kvec.Px()*ivec.Px() + kvec.Py()*ivec.Py() + kvec.Pz()*ivec.Pz())/kvec.P()/ivec.P();
						float oa = acos(cosoa)*180.0/const_pi;
						//if ( oa>90.0 ) bgood = false;

						int _pt_ind = hpt->FindBin(pvec.Pt());

						hoa_pt_mixed1[ind_sign]->Fill(pvec.Pt(), oa);

						if ( oa<90.0 ){
							hmass_pt_cut_mixed1[ind_sign]->Fill(pvec.M(), pvec.Pt());
						}
						if ( oa<100.0-10*_pt_ind ){
							hmass_pt_cut2_mixed1[ind_sign]->Fill(pvec.M(), pvec.Pt());
						}

					}//kp
				}else if ( abs(i_p_id[ip])==11 ){

					for (unsigned int kp=0; kp<vec_Xi_pt.size(); kp++){

						//11:e-, 3312:Xi-
						int ind_sign = (i_p_id[ip]*vec_Xi_id[kp]>0) ? 1 : 0;

						TLorentzVector kvec;
						kvec.SetPtEtaPhiM(vec_Xi_pt[kp], vec_Xi_eta[kp], vec_Xi_phi[kp], vec_Xi_mass[kp]);

						TLorentzVector pvec = ivec + kvec;

						hmass_pt_mixed2[ind_sign]->Fill(pvec.M(), pvec.Pt());

						float cosoa = (kvec.Px()*ivec.Px() + kvec.Py()*ivec.Py() + kvec.Pz()*ivec.Pz())/kvec.P()/ivec.P();
						float oa = acos(cosoa)*180.0/const_pi;
						//if ( oa>90.0 ) bgood = false;
						
						int _pt_ind = hpt->FindBin(pvec.Pt());

						hoa_pt_mixed2[ind_sign]->Fill(pvec.Pt(), oa);

						if ( oa<90.0 ){
							hmass_pt_cut_mixed2[ind_sign]->Fill(pvec.M(), pvec.Pt());
						}
						if ( oa<100.0-10*_pt_ind ){
							hmass_pt_cut2_mixed2[ind_sign]->Fill(pvec.M(), pvec.Pt());
						}

					}//kp

				}//

			}//ip

			//mixed-event pool
			vec_e_pt.clear();
			vec_e_eta.clear();
			vec_e_phi.clear();
			vec_e_mass.clear();
			vec_e_id.clear();

			vec_Xi_pt.clear();
			vec_Xi_eta.clear();
			vec_Xi_phi.clear();
			vec_Xi_mass.clear();
			vec_Xi_id.clear();

			for (int ip=0; ip<i_np; ip++){

				if ( fabs(f_p_eta[ip])>0.8 ) continue;

				if ( abs(i_p_id[ip])==11 ){

					if ( f_p_pt[ip]<0.5 ) continue;

					vec_e_pt.push_back(f_p_pt[ip]);
					vec_e_eta.push_back(f_p_eta[ip]);
					vec_e_phi.push_back(f_p_phi[ip]);
					vec_e_mass.push_back(f_p_mass[ip]);
					vec_e_id.push_back(i_p_id[ip]);
				}else if ( abs(i_p_id[ip])==3312 ){
					vec_Xi_pt.push_back(f_p_pt[ip]);
					vec_Xi_eta.push_back(f_p_eta[ip]);
					vec_Xi_phi.push_back(f_p_phi[ip]);
					vec_Xi_mass.push_back(f_p_mass[ip]);
					vec_Xi_id.push_back(i_p_id[ip]);
				}

			}//ip

		}//ien

		delete infile;

	}//

	TFile *outfile = new TFile(Form("outfile_hist_set%d_grp%03d.root",set,grp),"recreate");

	hmass_pt_sig3Xic0->Write();
	hmass_pt_sig4Xic0->Write();
	hmass_pt_sig4Xicp->Write();

	hmass_pt_cut_sig3Xic0->Write();
	hmass_pt_cut_sig4Xic0->Write();
	hmass_pt_cut_sig4Xicp->Write();

	hmass_pt_cut2_sig3Xic0->Write();
	hmass_pt_cut2_sig4Xic0->Write();
	hmass_pt_cut2_sig4Xicp->Write();

	for (int ii=0; ii<2; ii++){
		hmass_pt_corrJet[ii]->Write();
		hmass_pt_corrXibm[ii]->Write();
		hmass_pt_corrXib0[ii]->Write();
		hmass_pt_corrLambdab0[ii]->Write();
		hmass_pt_corrOmegabm[ii]->Write();
		hmass_pt_corrOmegac0[ii]->Write();
		hmass_pt_corrXiccp[ii]->Write();
		hmass_pt_corrB0[ii]->Write();
		hmass_pt_corrDiquark[ii]->Write();
		hmass_pt_comb[ii]->Write();
		hmass_pt_mixed1[ii]->Write();
		hmass_pt_mixed2[ii]->Write();

		hoa_pt_comb[ii]->Write();
		hoa_pt_mixed1[ii]->Write();
		hoa_pt_mixed2[ii]->Write();
		hoa_pt_sig3Xic0[ii]->Write();
		hoa_pt_sig4Xic0[ii]->Write();
		hoa_pt_sig4Xicp[ii]->Write();
		hoa_pt_corrJet[ii]->Write();
		hoa_pt_corrDiquark[ii]->Write();
		hoa_pt_corrXibm[ii]->Write();

		hmass_pt_cut_corrJet[ii]->Write();
		hmass_pt_cut_corrXibm[ii]->Write();
		hmass_pt_cut_corrXib0[ii]->Write();
		hmass_pt_cut_corrLambdab0[ii]->Write();
		hmass_pt_cut_corrOmegabm[ii]->Write();
		hmass_pt_cut_corrOmegac0[ii]->Write();
		hmass_pt_cut_corrXiccp[ii]->Write();
		hmass_pt_cut_corrB0[ii]->Write();
		hmass_pt_cut_corrDiquark[ii]->Write();
		hmass_pt_cut_comb[ii]->Write();
		hmass_pt_cut_mixed1[ii]->Write();
		hmass_pt_cut_mixed2[ii]->Write();

		hmass_pt_cut2_corrJet[ii]->Write();
		hmass_pt_cut2_corrXibm[ii]->Write();
		hmass_pt_cut2_corrXib0[ii]->Write();
		hmass_pt_cut2_corrLambdab0[ii]->Write();
		hmass_pt_cut2_corrOmegabm[ii]->Write();
		hmass_pt_cut2_corrOmegac0[ii]->Write();
		hmass_pt_cut2_corrXiccp[ii]->Write();
		hmass_pt_cut2_corrB0[ii]->Write();
		hmass_pt_cut2_corrDiquark[ii]->Write();
		hmass_pt_cut2_comb[ii]->Write();
		hmass_pt_cut2_mixed1[ii]->Write();
		hmass_pt_cut2_mixed2[ii]->Write();
	}//

	hoa_mass_comb->Write();
	hoa_mass_corrJet->Write();
	hoa_mass_sig3Xic0->Write();

	hphi_e->Write();
	hphi_Xi->Write();

	/*
	htrk_eta_pt->Write();

	hjetR04_eta_pt->Write();
	hjetR04_trk_dphi_deta->Write();

	hjetR04_jt_pt->Write();
	hjetR04_zh_pt->Write();
	*/


	outfile->Close();

}

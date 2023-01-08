//Jet jT, zh

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

using namespace std;

void make_hist_pp5TeV_01(const char *fname="file.lst", bool isMB=false){

	//CMS pt range
	//const float pt_min = 0.3;
	//const float pt_max = 3.0;
	//const float eta_max = 2.5;
	//const float deta_max = 5.0;
	//ATLAS pt range
	//const float pt_min = 0.5;
	//const float pt_max = 5.0;
	//const float eta_max = 2.5;
	//const float deta_max = 5.0;
	//ALICE pt range
	const float pt_min = 0.2;
	const float pt_max = 4.0;
	const float eta_max = 0.9;
	const float deta_max = 1.8;
	
	const float deta_cut = 2.0;
	const float const_pi = TMath::Pi();

	const int nmult = 3;
	//const float cut_mult[nmult+1] = {100, 60, 20, 5, 1, 0};

	const int npt = 5;
	const float ptbin[npt+1] = {0.2, 1.0, 2.0, 3.0, 4.0, 5.0};

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int i_np, i_njetR04;

	int i_p_id[2000];
	bool b_p_chg[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];

	float f_jetR04_pt[100], f_jetR04_eta[100], f_jetR04_phi[100];

	float f_scale;


	TH1D *hevent_mult = new TH1D("hevent_mult","",200,0,200);

	TH2D *htrack_eta_pt[nmult];
	TH2D *hjetR04_eta_pt[nmult];
	TH2D *hjetR04_eta_subpt[nmult];
	TH2D *hjetR04_eta_subpt2[nmult];
	TH2D *hjetR04_eta_bkgpt[nmult];

	TProfile *hprof_pt_scale[nmult];

	for (int im=0; im<nmult; im++){
		htrack_eta_pt[im] = new TH2D(Form("htrack_eta_pt_im%d",im),"",20,-1,1,100,0,10);
		hjetR04_eta_pt[im] = new TH2D(Form("hjetR04_eta_pt_im%d",im),"",20,-1.0,1.0,100,0,100);
		hjetR04_eta_subpt[im] = new TH2D(Form("hjetR04_eta_subpt_im%d",im),"",20,-1.0,1.0,100,0,100);
		hjetR04_eta_subpt2[im] = new TH2D(Form("hjetR04_eta_subpt2_im%d",im),"",20,-1.0,1.0,100,0,100);
		hjetR04_eta_bkgpt[im] = new TH2D(Form("hjetR04_eta_bkgpt_im%d",im),"",20,-1.0,1.0,100,0,100);

		hprof_pt_scale[im] = new TProfile(Form("hprof_pt_scale_im%d",im),"",20,0,100);
	}


	TH2D *hjetR04_trk_dphi_deta = new TH2D("hjetR04_trk_dphi_deta","",96,-const_pi/2,3*const_pi/2,50,-2.5,2.5); 

	TH2D *hjetR04_jt_pt = new TH2D("hjetR04_jt_pt","",100,0,10,10,0,100);
	TH2D *hjetR04_zh_pt = new TH2D("hjetR04_zh_pt","",100,0,1,10,0,100);

	/*
	TH2D *hevent_Q_mult[nmult];
	TH2D *hevent_nMPI_mult[nmult];
	TH2D *hevent_bMPI_mult[nmult];

	TH2D *hevent_Q_nMPI_mult[nmult];
	TH2D *hevent_Q_bMPI_mult[nmult];

	for (int im=0; im<nmult; im++){
		hevent_Q_mult[im] = new TH2D(Form("hevent_Q_mult%02d",im),"",50*2,0,50,100*5,0,100);
		hevent_nMPI_mult[im] = new TH2D(Form("hevent_nMPI_mult%02d",im),"",50*2,0,50,51,-0.5,50.5);
		hevent_bMPI_mult[im] = new TH2D(Form("hevent_bMPI_mult%02d",im),"",50*2,0,50,150,0,3);

		hevent_Q_nMPI_mult[im] = new TH2D(Form("hevent_Q_nMPI_mult%02d",im),"",50*2,0,50,51,-0.5,50.5);
		hevent_Q_bMPI_mult[im] = new TH2D(Form("hevent_Q_bMPI_mult%02d",im),"",50*2,0,50,150,0,3);
	}
	*/

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");

		T->SetBranchAddress("scale",&f_scale);

		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id",i_p_id);
		T->SetBranchAddress("p_chg",b_p_chg);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_pt",f_p_pt);

		T->SetBranchAddress("njetR04",&i_njetR04);
		T->SetBranchAddress("jetR04_pt",f_jetR04_pt);
		T->SetBranchAddress("jetR04_eta",f_jetR04_eta);
		T->SetBranchAddress("jetR04_phi",f_jetR04_phi);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			int mult = 0;

			//charged particle
			for (int ip=0; ip<i_np; ip++){
				if ( !b_p_chg[ip] ) continue;
				if ( fabs(f_p_eta[ip])>1 ) continue;

				mult++;
			}

			hevent_mult->Fill(mult);

			//charged particle
			for (int ip=0; ip<i_np; ip++){
				if ( !b_p_chg[ip] ) continue;
				if ( fabs(f_p_eta[ip])>1 ) continue;

				htrack_eta_pt[0]->Fill(f_p_eta[ip], f_p_pt[ip]);

				if ( mult>=60 ){
					htrack_eta_pt[1]->Fill(f_p_eta[ip], f_p_pt[ip]);
				}

				if ( mult>=60 && mult<80 ){
					htrack_eta_pt[2]->Fill(f_p_eta[ip], f_p_pt[ip]);
				}
			}

			if ( i_njetR04==0 ) continue;

			//leading jet and rorate
			TVector3 vec_leadjet;
			vec_leadjet.SetPtEtaPhi(f_jetR04_pt[0], f_jetR04_eta[0], f_jetR04_phi[0]);
			TVector3 vOrtho(vec_leadjet);
			TVector3 perpjet;
			vOrtho.RotateZ(const_pi/2.);

			float bkg_pt = 0.0;
			float bkg_ptot = 0.0;

			//background estimation
			for (int ip=0; ip<i_np; ip++){
				if ( fabs(f_p_eta[ip])>1 ) continue;
				if ( abs(i_p_id[ip])==12 || abs(i_p_id[ip])==14 || abs(i_p_id[ip])==16 ) continue;

				float deta = f_p_eta[ip] - vOrtho.Eta();
				float dphi = f_p_phi[ip] - vOrtho.Phi();

				if ( dphi < -const_pi ) dphi += 2*const_pi;
				else if ( dphi > const_pi ) dphi -= 2*const_pi;

				float newdR = sqrt(deta*deta + dphi*dphi);
				if(newdR>0.4) continue; 

				TVector3 vec;
				vec.SetPtEtaPhi(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip]);

				bkg_pt += f_p_pt[ip];
				bkg_ptot += vec.Mag(); 
			}

			//leading jet
			float leadjet_pt = 0;
			for (int ijet=0; ijet<i_njetR04; ijet++){
				if ( fabs(f_jetR04_eta[ijet])>0.5 ) continue; 

				leadjet_pt = f_jetR04_pt[ijet]-bkg_pt;
				break;
			}

			if ( isMB ){
				hprof_pt_scale[0]->Fill(leadjet_pt, f_scale);
			}else{
				if ( mult>=60 ){
					hprof_pt_scale[1]->Fill(leadjet_pt, f_scale);
				}

				if ( mult>=60 && mult<80 ){
					hprof_pt_scale[2]->Fill(leadjet_pt, f_scale);
				}
			}


			//jetR04 info
			for (int ijet=0; ijet<i_njetR04; ijet++){
				if ( fabs(f_jetR04_eta[ijet])>0.5 ) continue; 

				TVector3 vec;
				vec.SetPtEtaPhi(f_jetR04_pt[ijet], f_jetR04_eta[ijet], f_jetR04_phi[ijet]);
				vec.SetMag(vec.Mag() - bkg_ptot);

				hjetR04_eta_pt[0]->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]);
				hjetR04_eta_subpt[0]->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]-bkg_pt);
				hjetR04_eta_subpt2[0]->Fill(f_jetR04_eta[ijet], vec.Pt());
				hjetR04_eta_bkgpt[0]->Fill(f_jetR04_eta[ijet], bkg_pt);

				if ( mult>=60 ){
					hjetR04_eta_pt[1]->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]);
					hjetR04_eta_subpt[1]->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]-bkg_pt);
					hjetR04_eta_subpt2[1]->Fill(f_jetR04_eta[ijet], vec.Pt());
					hjetR04_eta_bkgpt[1]->Fill(f_jetR04_eta[ijet], bkg_pt);
				}

				if ( mult>=60 && mult<80 ){
					hjetR04_eta_pt[2]->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]);
					hjetR04_eta_subpt[2]->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]-bkg_pt);
					hjetR04_eta_subpt2[2]->Fill(f_jetR04_eta[ijet], vec.Pt());
					hjetR04_eta_bkgpt[2]->Fill(f_jetR04_eta[ijet], bkg_pt);
				}
			}//ijet


		}//ien

		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	hevent_mult->Write();
	for (int im=0; im<nmult; im++){
		htrack_eta_pt[im]->Write();
		hjetR04_eta_pt[im]->Write();
		hjetR04_eta_subpt[im]->Write();
		hjetR04_eta_subpt2[im]->Write();
		hjetR04_eta_bkgpt[im]->Write();

		hprof_pt_scale[im]->Write();
	}

	outfile->Close();

}

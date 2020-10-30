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

void make_hist_pp13TeV(const char *fname="file.lst"){

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

	const int nmult = 4;
	//const float cut_mult[nmult+1] = {100, 60, 20, 5, 1, 0};

	const int npt = 5;
	const float ptbin[npt+1] = {0.2, 1.0, 2.0, 3.0, 4.0, 5.0};

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np, i_njetR04;
	int i_p_id[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	float f_jetR04_pt[100], f_jetR04_eta[100], f_jetR04_phi[100];


	TH1D *hevent_Q = new TH1D("hevent_Q","",1500,0,150);
	TH1D *hevent_nMPI = new TH1D("hevent_nMPI","",51,-0.5,50.5);
	TH1D *hevent_bMPI = new TH1D("hevent_bMPI","",300,0,3);
	TH2D *hevent_mult_mid_fwd = new TH2D("hevent_mult_mid_fwd","",200,0,200,200,0,200);
	TH2D *hevent_mult_fwd = new TH2D("hevent_mult_fwd","",nmult,0,nmult,200,0,200);
	TH2D *hevent_mult_mid = new TH2D("hevent_mult_mid","",nmult,0,nmult,200,0,200);

	TH2D *htrk_eta_pt = new TH2D("htrk_eta_pt","",100,-5,5,100,0,10);

	TH2D *hjetR04_eta_pt = new TH2D("hjetR04_eta_pt","",50,-2.5,2.5,100,0,200);
	TH2D *hjetR04_trk_dphi_deta = new TH2D("hjetR04_trk_dphi_deta","",96,-const_pi/2,3*const_pi/2,50,-2.5,2.5); 

	TH2D *hjetR04_jt_pt = new TH2D("hjetR04_jt_pt","",100,0,10,20,0,100);
	TH2D *hjetR04_zh_pt = new TH2D("hjetR04_zh_pt","",100,0,1,20,0,100);

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
		T->SetBranchAddress("bMPI",&f_bMPI);
		T->SetBranchAddress("nMPI",&i_nMPI);
		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id",i_p_id);
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

			hevent_Q->Fill(f_scale);
			hevent_nMPI->Fill(i_nMPI);
			hevent_bMPI->Fill(f_bMPI);

			//jetR04 info
			for (int ijet=0; ijet<i_njetR04; ijet++){
				hjetR04_eta_pt->Fill(f_jetR04_eta[ijet], f_jetR04_pt[ijet]);
			}//ijet

			//particle info
			for (int ip=0; ip<i_np; ip++){
				htrk_eta_pt->Fill(f_p_eta[ip], f_p_pt[ip]);
			}//ip

			//jetR04 loop 
			for (int ijet=0; ijet<i_njetR04; ijet++){

				if ( fabs(f_jetR04_eta[ijet])>0.4 ) continue;

				TVector3 vec_jet;
				vec_jet.SetPtEtaPhi(f_jetR04_pt[ijet], f_jetR04_eta[ijet], f_jetR04_phi[ijet]);

				//particle info
				for (int ip=0; ip<i_np; ip++){

					if ( fabs(f_p_eta[ip])>1.0 || fabs(f_p_pt[ip])<0.2 ) continue;

					float deta = f_p_eta[ip] - f_jetR04_eta[ijet];
					float dphi = f_p_phi[ip] - f_jetR04_phi[ijet];

					if ( dphi<-1./2.*const_pi ) dphi += 2*const_pi;
					else if ( dphi>5./2.*const_pi ) dphi -= 2*const_pi;

					hjetR04_trk_dphi_deta->Fill(dphi, deta);

					float dR = sqrt(deta*deta + dphi*dphi);

					if ( dR>0.4 ) continue;

					TVector3 vec_p;
					vec_p.SetPtEtaPhi(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip]);

					TVector3 vec_cross = vec_jet.Cross(vec_p); 
					float f_dot = vec_jet.Dot(vec_p);

					float f_jt = vec_cross.Mag()/vec_jet.Mag();
					hjetR04_jt_pt->Fill(f_jt, f_jetR04_pt[ijet]);

					float f_zh = f_dot/vec_jet.Mag2();
					hjetR04_zh_pt->Fill(f_zh, f_jetR04_pt[ijet]);

				}//ip

			}//ijet


		}//ien

		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	htrk_eta_pt->Write();

	hjetR04_eta_pt->Write();
	hjetR04_trk_dphi_deta->Write();

	hjetR04_jt_pt->Write();
	hjetR04_zh_pt->Write();


	outfile->Close();

}

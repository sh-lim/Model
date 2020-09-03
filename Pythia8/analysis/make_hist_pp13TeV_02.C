//CMS acceptance (vs pT)

#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>

using namespace std;

void make_hist_pp13TeV_02(const char *fname="file.lst"){

	//CMS pt range
	const float pt_min = 0.2;
	const float pt_max = 6.0;
	const float eta_max = 2.5;
	const float deta_max = 5.0;
	//ATLAS pt range
	//const float pt_min = 0.5;
	//const float pt_max = 5.0;
	//const float eta_max = 2.5;
	//const float deta_max = 5.0;
	//ALICE pt range
	//const float pt_min = 1.0;
	//const float pt_max = 2.0;
	//const float eta_max = 0.9;
	//const float deta_max = 1.8;
	
	const float deta_cut = 2.0;
	const float const_pi = TMath::Pi();

	const int nmult = 3;
	//const float cut_mult[nmult+1] = {100, 60, 20, 5, 1, 0};

	//const int npt = 4;
	//const float ptbin[npt+1] = {0.2, 1.0, 2.0, 3.0, 4.0};

	const int npt = 8;
	const float ptbin[npt+1] = {0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0};

	/*
	TFile *infile_ref = new TFile("/alice/home/shlim/work/Pythia/jobs/outfile_hist_pp13TeV_set00_grp000_try000.root","read");
	TH2D *href_mid_fwd = (TH2D*)infile_ref->Get("hevent_mult_mid_fwd")->Clone("href_mid_fwd");
	TH1D *href_fwd = (TH1D*)href_mid_fwd->ProjectionY("href_fwd");
	*/

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np;
	int i_p_id[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];

	TH1D *hntrig_same[nmult];
	TH1D *hntrig_mixed[nmult];

	TH2D *h2d_same_dphi_deta[nmult][npt];
	TH2D *h2d_mixed_dphi_deta[nmult][npt];

	for (int im=0; im<nmult; im++){
		hntrig_same[im] = new TH1D(Form("hntrig_same_mult%02d",im),"",npt,ptbin);
		hntrig_mixed[im] = new TH1D(Form("hntrig_mixed_mult%02d",im),"",npt,ptbin);

		for (int ipt=0; ipt<npt; ipt++){
			h2d_same_dphi_deta[im][ipt] = new TH2D(Form("h2d_same_dphi_deta_mult%02d_pt%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
			h2d_mixed_dphi_deta[im][ipt] = new TH2D(Form("h2d_mixed_dphi_deta_mult%02d_pt%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		}//ipt
	}//im

	TH2D *heta_pt_ana = new TH2D("heta_pt_ana","",100,-5,5,100,0,10);

	TH1D *hevent_Q = new TH1D("hevent_Q","",1500,0,150);
	TH1D *hevent_nMPI = new TH1D("hevent_nMPI","",51,-0.5,50.5);
	TH1D *hevent_bMPI = new TH1D("hevent_bMPI","",300,0,3);
	TH2D *hevent_mult_mid_fwd = new TH2D("hevent_mult_mid_fwd","",200,0,200,200,0,200);
	TH2D *hevent_mult_fwd = new TH2D("hevent_mult_fwd","",nmult,0,nmult,200,0,200);
	TH2D *hevent_mult_mid = new TH2D("hevent_mult_mid","",nmult,0,nmult,200,0,200);


	vector<float> vec_pt[nmult];
	vector<float> vec_phi[nmult];
	vector<float> vec_eta[nmult];

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

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			hevent_Q->Fill(f_scale);
			hevent_nMPI->Fill(i_nMPI);
			hevent_bMPI->Fill(f_bMPI);

			int nmult_mid_nocut = 0;
			int nmult_mid = 0, nmult_fwd = 0;

			//count multiplicity for bin
			for (int ip=0; ip<i_np; ip++){

				if ( fabs(f_p_eta[ip])<2.5 ){
					nmult_mid_nocut++;
				}

				if ( fabs(f_p_eta[ip])<2.5 && f_p_pt[ip]>0.4 ){
					nmult_mid++;
				}

				if ( fabs(f_p_eta[ip])>3.1 && fabs(f_p_eta[ip])<4.9 ){
					nmult_fwd++;
				}
			}//i_np

			if ( nmult_mid_nocut>0 ){
				hevent_mult_mid_fwd->Fill(nmult_mid, nmult_fwd);
			}

			int ind_mult = -1;

			if ( nmult_mid<10 ) ind_mult = 0;
			else if ( nmult_mid>=10 && nmult_mid<20 ) ind_mult = 1;
			else if ( nmult_mid>=120 ) ind_mult = 2;

			if ( ind_mult<0 ) continue;

			hevent_mult_mid->Fill(ind_mult+0.5, nmult_mid); 
			hevent_mult_fwd->Fill(ind_mult+0.5, nmult_fwd); 

			for (int ip=0; ip<i_np; ip++){

				//same events
				if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;
				if ( fabs(f_p_eta[ip])>eta_max ) continue;

				hntrig_same[ind_mult]->Fill(f_p_pt[ip]);

				int ind_pt_i = hntrig_same[ind_mult]->FindBin(f_p_pt[ip]) - 1;

				for (int jp=0; jp<i_np; jp++){
					if ( ip==jp ) continue;
					if ( f_p_pt[jp]<pt_min || f_p_pt[jp]>pt_max ) continue;
					if ( fabs(f_p_eta[jp])>eta_max ) continue;

					int ind_pt_j = hntrig_same[ind_mult]->FindBin(f_p_pt[jp]) - 1;

					if ( ind_pt_i!=ind_pt_j ) continue;

					float dphi = f_p_phi[jp] - f_p_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = f_p_eta[jp] - f_p_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_same_dphi_deta[ind_mult][ind_pt_i]->Fill(dphi, deta);
				}//jp

				//mixed events
				if ( vec_pt[ind_mult].size()>0 ){
					hntrig_mixed[ind_mult]->Fill(f_p_pt[ip]);
				}

				for (unsigned int jj=0; jj<vec_pt[ind_mult].size(); jj++){
					//if ( vec_pt[ind_mult][jj]<pt_min || vec_pt[ind_mult][jj]>pt_max ) continue;
					//if ( fabs(vec_eta[ind_mult][jj])>eta_max ) continue;
					
					int ind_pt_j = hntrig_same[ind_mult]->FindBin(vec_pt[ind_mult][jj]) - 1;

					if ( ind_pt_i!=ind_pt_j ) continue;

					float dphi = vec_phi[ind_mult][jj] - f_p_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_eta[ind_mult][jj] - f_p_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_mixed_dphi_deta[ind_mult][ind_pt_j]->Fill(dphi, deta);
				}//jj

			}//ip

			//fill mixed event pool
			vec_pt[ind_mult].clear();
			vec_phi[ind_mult].clear();
			vec_eta[ind_mult].clear();

			for (int ip=0; ip<i_np; ip++){
				if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;
				if ( fabs(f_p_eta[ip])>eta_max ) continue;

				vec_pt[ind_mult].push_back(f_p_pt[ip]);
				vec_phi[ind_mult].push_back(f_p_phi[ip]);
				vec_eta[ind_mult].push_back(f_p_eta[ip]);
			}

		}//ien

		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	for (int im=0; im<nmult; im++){
		hntrig_same[im]->Write();
		hntrig_mixed[im]->Write();

		for (int ipt=0; ipt<npt; ipt++){
			h2d_same_dphi_deta[im][ipt]->Write();
			h2d_mixed_dphi_deta[im][ipt]->Write();
		}
	}

	hevent_Q->Write();
	hevent_nMPI->Write();
	hevent_bMPI->Write();

	hevent_mult_mid_fwd->Write();
	hevent_mult_mid->Write();
	hevent_mult_fwd->Write();

	outfile->Close();

}

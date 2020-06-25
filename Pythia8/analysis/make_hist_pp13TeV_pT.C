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

void make_hist_pp13TeV_pT(const char *fname="file.lst"){

	//CMS pt range
	//const float pt_min = 0.3;
	//const float pt_max = 3.0;
	//ATLAS pt range
	const float pt_min = 0.5;
	const float pt_max = 5.0;
	const float eta_max = 2.5;
	const float deta_max = 5.0;
	//const float eta_max = 0.9;
	//const float deta_max = 1.8;
	
	const float deta_cut = 2.0;
	const float const_pi = TMath::Pi();

	const int nmult = 5;
	//const double multcut_mid[nmult+1] = {0, 10, 20, 30, 85, 500}; //100-70-50-35-5-0
	//const double multcut_fwd[nmult+1] = {0, 15, 25, 35, 85, 500}; //100-70-50-35-5-0
	//const double multcut_mid[nmult+1] = {0, 10, 20, 30, 70, 500}; //100-70-50-35-10-0
	//const double multcut_fwd[nmult+1] = {0, 15, 25, 35, 70, 500}; //100-70-50-35-10-0
	const double multcut_mid[nmult+1] = {0, 14, 26, 49, 69, 500}; //100-60-40-20-10-0
	const double multcut_fwd[nmult+1] = {0, 19, 31, 52, 70, 500}; //100-60-40-20-10-0

	//const int npt = 8;
	//const float ptcut[npt+1] = {0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};
	const int npt = 12;
	const double ptcut[npt+1] = 
	{
		0.5, 1.0, 1.5, 2.0, 2.5,
		3.0, 4.0, 5.0, 7.5, 10.0,
		15.0, 20.0, 30.0
	};

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int npart;
	int part_pid[1000];
	float part_eta[1000], part_phi[1000], part_pt[1000];

	//TH2D *h2d_same_dphi_deta[nmult];
	//TH2D *h2d_mixed_dphi_deta[nmult];

	TH2D *h2d_same_dphi_deta[nmult][npt];
	TH2D *h2d_mixed_dphi_deta[nmult][npt];

	TH2D *h2d_mid2_same_dphi_deta[nmult][npt];
	TH2D *h2d_mid2_mixed_dphi_deta[nmult][npt];

	TH2D *h2d_fwd_same_dphi_deta[nmult][npt];
	TH2D *h2d_fwd_mixed_dphi_deta[nmult][npt];

	TH2D *h2d_fwd2_same_dphi_deta[nmult][npt];
	TH2D *h2d_fwd2_mixed_dphi_deta[nmult][npt];

	for (int im=0; im<nmult; im++){
		for (int ipt=0; ipt<npt; ipt++){
			h2d_same_dphi_deta[im][ipt] = new TH2D(Form("h2d_same_dphi_deta_m%02d_pt%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
			h2d_mixed_dphi_deta[im][ipt] = new TH2D(Form("h2d_mixed_dphi_deta_m%02d_pt%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);

			h2d_mid2_same_dphi_deta[im][ipt] = new TH2D(Form("h2d_mid2_same_dphi_deta_m%02d_pt%02d",im,ipt),"",96,0,2*const_pi,100,-5.0,5.0);
			h2d_mid2_mixed_dphi_deta[im][ipt] = new TH2D(Form("h2d_mid2_mixed_dphi_deta_m%02d_pt%02d",im,ipt),"",96,0,2*const_pi,100,-5.0,5.0);

			h2d_fwd_same_dphi_deta[im][ipt] = new TH2D(Form("h2d_fwd_same_dphi_deta_m%02d_pt%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
			h2d_fwd_mixed_dphi_deta[im][ipt] = new TH2D(Form("h2d_fwd_mixed_dphi_deta_m%02d_pt%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);

			h2d_fwd2_same_dphi_deta[im][ipt] = new TH2D(Form("h2d_fwd2_same_dphi_deta_m%02d_pt%02d",im,ipt),"",96,0,2*const_pi,100,-5.0,5.0);
			h2d_fwd2_mixed_dphi_deta[im][ipt] = new TH2D(Form("h2d_fwd2_mixed_dphi_deta_m%02d_pt%02d",im,ipt),"",96,0,2*const_pi,100,-5.0,5.0);
		}
	}

	TH1D *hmult_all = new TH1D("hmult_all","",500,0,500);
	TH1D *hmult_mid = new TH1D("hmult_mid","",500,0,500);
	TH1D *hmult_fwd = new TH1D("hmult_fwd","",500,0,500);
	TH2D *hmult_mid_fwd = new TH2D("hmult_mid_fwd","",100,0,500,100,0,500);
	TH2D *hmult_cnt_fwd = new TH2D("hmult_cnt_fwd","",100,0,500,100,0,500);

	TH2D *hntrig_mid = new TH2D("hntrig_mid","",nmult,multcut_mid,npt,ptcut);
	TH2D *hntrig_mixed_mid = new TH2D("hntrig_mixed_mid","",nmult,multcut_mid,npt,ptcut);

	TH2D *hntrig_fwd = new TH2D("hntrig_fwd","",nmult,multcut_fwd,npt,ptcut);
	TH2D *hntrig_mixed_fwd = new TH2D("hntrig_mixed_fwd","",nmult,multcut_fwd,npt,ptcut);

	TH2D *heta_pt_all = new TH2D("heta_pt_all","",100,-5,5,300,0,30);
	TH2D *heta_pt_ana = new TH2D("heta_pt_ana","",100,-5,5,300,0,30);

	vector<float> vec_pt[nmult];
	vector<float> vec_phi[nmult];
	vector<float> vec_eta[nmult];

	vector<float> vec_fwd_pt[nmult];
	vector<float> vec_fwd_phi[nmult];
	vector<float> vec_fwd_eta[nmult];

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("npart",&npart);
		T->SetBranchAddress("part_pid",part_pid);
		T->SetBranchAddress("part_eta",part_eta);
		T->SetBranchAddress("part_phi",part_phi);
		T->SetBranchAddress("part_pt",part_pt);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
		//for (int ien=0; ien<1000; ien++){
			T->GetEntry(ien);

			hmult_all->Fill(npart);

			int nmult_mid = 0, nmult_fwd = 0, nmult_cnt = 0;

			//count multiplicity for bin
			for (int ip=0; ip<npart; ip++){
				if ( fabs(part_eta[ip])<2.5 && part_pt[ip]>0.4 ){
					nmult_mid++;
				}

				if ( fabs(part_eta[ip])<1.0 && part_pt[ip]>0.4 ){
					nmult_cnt++;
				}

				if ( fabs(part_eta[ip])>3.1 && fabs(part_eta[ip])<5.0 ){
					nmult_fwd++;
				}

				heta_pt_all->Fill(part_eta[ip], part_pt[ip]);
			}

			hmult_mid->Fill(nmult_mid);
			hmult_fwd->Fill(nmult_fwd);
			hmult_mid_fwd->Fill(nmult_mid,nmult_fwd);
			hmult_cnt_fwd->Fill(nmult_cnt,nmult_fwd);

			//continue;

			int imult = hntrig_mid->GetXaxis()->FindBin(nmult_mid) - 1; 
			int imult_fwd = hntrig_fwd->GetXaxis()->FindBin(nmult_fwd) - 1; 
			if ( imult>=nmult ) continue;

			for (int ip=0; ip<npart; ip++){
				//int ip_ptbin = int(part_pt[ip]/0.5);
				//if ( ip_ptbin>=20 ) continue;

				//if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
				if ( part_pt[ip]<ptcut[0] || part_pt[ip]>ptcut[npt] ) continue;
				if ( fabs(part_eta[ip])>eta_max ) continue;

				hntrig_mid->Fill(nmult_mid, part_pt[ip]);
				hntrig_fwd->Fill(nmult_fwd, part_pt[ip]);
				heta_pt_ana->Fill(part_eta[ip], part_pt[ip]);

				int ipt = hntrig_mid->GetYaxis()->FindBin(part_pt[ip]) - 1;

				if ( ipt>=npt ){
					cout << "pT: " << part_pt[ip] << endl;
					continue;
				}

				for (int jp=0; jp<npart; jp++){
					if ( ip==jp ) continue;
					if ( part_pt[jp]<pt_min || part_pt[jp]>pt_max ) continue;
					if ( fabs(part_eta[jp])>eta_max ) continue;

					float dphi = part_phi[jp] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = part_eta[jp] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_same_dphi_deta[imult][ipt]->Fill(dphi, deta);
					h2d_fwd_same_dphi_deta[imult_fwd][ipt]->Fill(dphi, deta);

					float dphi2 = part_phi[jp] - part_phi[ip];
					if ( dphi2<0 ) dphi2 += 2*const_pi;

					h2d_mid2_same_dphi_deta[imult][ipt]->Fill(dphi2, deta);
					h2d_fwd2_same_dphi_deta[imult_fwd][ipt]->Fill(dphi2, deta);

				}//jp

				if ( vec_pt[imult].size()>0 ){
					hntrig_mixed_mid->Fill(nmult_mid, part_pt[ip]);
				}

				//mixed events
				for (unsigned int jj=0; jj<vec_pt[imult].size(); jj++){
					if ( vec_pt[imult][jj]<pt_min || vec_pt[imult][jj]>pt_max ) continue;
					if ( fabs(vec_eta[imult][jj])>eta_max ) continue;

					float dphi = vec_phi[imult][jj] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_eta[imult][jj] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_mixed_dphi_deta[imult][ipt]->Fill(dphi, deta);

					float dphi2 = vec_phi[imult][jj] - part_phi[ip];
					if ( dphi2<0 ) dphi2 += 2*const_pi;

					h2d_mid2_mixed_dphi_deta[imult][ipt]->Fill(dphi2, deta);
				}//jj

				if ( vec_fwd_pt[imult_fwd].size()>0 ){
					hntrig_mixed_fwd->Fill(nmult_fwd, part_pt[ip]);
				}

				//mixed events
				for (unsigned int jj=0; jj<vec_fwd_pt[imult_fwd].size(); jj++){
					if ( vec_fwd_pt[imult_fwd][jj]<pt_min || vec_fwd_pt[imult_fwd][jj]>pt_max ) continue;
					if ( fabs(vec_fwd_eta[imult_fwd][jj])>eta_max ) continue;

					float dphi = vec_fwd_phi[imult_fwd][jj] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_fwd_eta[imult_fwd][jj] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_fwd_mixed_dphi_deta[imult_fwd][ipt]->Fill(dphi, deta);

					float dphi2 = vec_fwd_phi[imult_fwd][jj] - part_phi[ip];
					if ( dphi2<0 ) dphi2 += 2*const_pi;

					h2d_fwd2_mixed_dphi_deta[imult_fwd][ipt]->Fill(dphi2, deta);
				}//jj

			}//ip


			//fill mixed event pool
			vec_pt[imult].clear();
			vec_phi[imult].clear();
			vec_eta[imult].clear();

			vec_fwd_pt[imult_fwd].clear();
			vec_fwd_phi[imult_fwd].clear();
			vec_fwd_eta[imult_fwd].clear();

			for (int ip=0; ip<npart; ip++){
				if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
				if ( fabs(part_eta[ip])>eta_max ) continue;
				vec_pt[imult].push_back(part_pt[ip]);
				vec_phi[imult].push_back(part_phi[ip]);
				vec_eta[imult].push_back(part_eta[ip]);

				vec_fwd_pt[imult_fwd].push_back(part_pt[ip]);
				vec_fwd_phi[imult_fwd].push_back(part_phi[ip]);
				vec_fwd_eta[imult_fwd].push_back(part_eta[ip]);
			}

		}//ien


		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	for (int im=0; im<nmult; im++){
		for (int ipt=0; ipt<npt; ipt++){
			h2d_same_dphi_deta[im][ipt]->Write();
			h2d_mixed_dphi_deta[im][ipt]->Write();

			h2d_fwd_same_dphi_deta[im][ipt]->Write();
			h2d_fwd_mixed_dphi_deta[im][ipt]->Write();

			h2d_mid2_same_dphi_deta[im][ipt]->Write();
			h2d_mid2_mixed_dphi_deta[im][ipt]->Write();

			h2d_fwd2_same_dphi_deta[im][ipt]->Write();
			h2d_fwd2_mixed_dphi_deta[im][ipt]->Write();
		}
	}

	heta_pt_all->Write();
	heta_pt_ana->Write();
	hmult_all->Write();
	hmult_mid->Write();
	hmult_fwd->Write();
	hmult_cnt_fwd->Write();
	hmult_mid_fwd->Write();
	hntrig_mid->Write();
	hntrig_mixed_mid->Write();
	hntrig_fwd->Write();
	hntrig_mixed_fwd->Write();
	//hv22pc_mult->Write();
	outfile->Close();

}

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

void make_hist_pAu200GeV(const char *fname="file.lst"){

	//STAR pt range
	//const float pt_min = 0.2;
	//const float pt_max = 3.0;
	const float pt_min = 0.001;
	const float pt_max = 10.0;
	//STAR eta range
	//const float eta_max = 0.9;
	//const float deta_max = 1.8;
	//LHC eta range
	const float eta_max = 2.5;
	const float deta_max = 5.0;
	const float const_pi = TMath::Pi();

	const int nmult = 20;
	//const int dmult = 5;
	const int dmult = 10;
	const int nmult_fwd = 5;

	const double multcut_fwd[nmult_fwd+1] = {0, 2, 7, 14, 22, 500}; //100-86-50-20-5

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int npart;
	int part_pid[1000];
	float part_eta[1000], part_phi[1000], part_pt[1000];

	int npp;
	int pp_pid[1000];
	float pp_x[1000], pp_y[1000], pp_z[1000];

	TH1D *hpp = new TH1D("hpp","",96,0,const_pi);

	TH2D *h2d_mid_same_dphi_deta[nmult];
	TH2D *h2d_mid_mixed_dphi_deta[nmult];

	TH2D *h2d_mb_same_dphi_deta[nmult];
	TH2D *h2d_mb_mixed_dphi_deta[nmult];

	TH2D *h2d_fwd_same_dphi_deta[nmult];
	TH2D *h2d_fwd_mixed_dphi_deta[nmult];

	for (int im=0; im<nmult; im++){
		h2d_mid_same_dphi_deta[im] = new TH2D(Form("h2d_mid_same_dphi_deta_m%02d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		h2d_mid_mixed_dphi_deta[im] = new TH2D(Form("h2d_mid_mixed_dphi_deta_m%02d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);

		h2d_mb_same_dphi_deta[im] = new TH2D(Form("h2d_mb_same_dphi_deta_m%02d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		h2d_mb_mixed_dphi_deta[im] = new TH2D(Form("h2d_mb_mixed_dphi_deta_m%02d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
	}

	for (int im=0; im<nmult_fwd; im++){
		h2d_fwd_same_dphi_deta[im] = new TH2D(Form("h2d_fwd_same_dphi_deta_%d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		h2d_fwd_mixed_dphi_deta[im] = new TH2D(Form("h2d_fwd_mixed_dphi_deta_%d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
	}

	TProfile *hv2pp_mid = new TProfile("hv2pp_mid","",20,0,100);
	TProfile *hv2pp_fwd = new TProfile("hv2pp_fwd","",nmult_fwd,multcut_fwd);
	TProfile *hv2pp_sq_mid = new TProfile("hv2pp_sq_mid","",20,0,100);
	TProfile *hv2pp_sq_fwd = new TProfile("hv2pp_sq_fwd","",nmult_fwd,multcut_fwd);

	TH2D *h2d_ecc2_mult_fwd = new TH2D("h2d_ecc2_mult_fwd","",nmult_fwd,multcut_fwd,100,0,1);
	TH2D *h2d_ecc2_sq_mult_fwd = new TH2D("h2d_ecc2_sq_mult_fwd","",nmult_fwd,multcut_fwd,100,0,1);
	TProfile *hprof_ecc2_mult_fwd = new TProfile("hprof_ecc2_mult_fwd","",nmult_fwd,multcut_fwd,0,1);
	TProfile *hprof_ecc2_sq_mult_fwd = new TProfile("hprof_ecc2_sq_mult_fwd","",nmult_fwd,multcut_fwd,0,1);

	TH2D *h2d_ecc2_mult_mid = new TH2D("h2d_ecc2_mult_mid","",20,0,100,100,0,1);
	TH2D *h2d_ecc2_sq_mult_mid = new TH2D("h2d_ecc2_sq_mult_mid","",20,0,100,100,0,1);
	TProfile *hprof_ecc2_mult_mid = new TProfile("hprof_ecc2_mult_mid","",20,0,100,0,1);
	TProfile *hprof_ecc2_sq_mult_mid = new TProfile("hprof_ecc2_sq_mult_mid","",20,0,100,0,1);

	TH1D *hmult_all = new TH1D("hmult_all","",500,0,500);
	TH1D *hmult_mid = new TH1D("hmult_mid","",500,0,500);
	TH1D *hmult_mid_mb = new TH1D("hmult_mid_mb","",500,0,500);
	TH1D *hmult_fwd = new TH1D("hmult_fwd","",500,0,500);

	//TProfile *hprof_mult_mid = new TProfile("hprof_mult_mid","",nmult,multcut_mid);
	//TProfile *hprof_mult_fwd = new TProfile("hprof_mult_fwd","",nmult,multcut_fwd);

	TH1D *hntrig_mid = new TH1D("hntrig_mid","",200,0,200);
	TH1D *hntrig_mb = new TH1D("hntrig_mb","",200,0,200);
	TH1D *hntrig_fwd = new TH1D("hntrig_fwd","",nmult_fwd,multcut_fwd);

	TH1D *hntrig_mixed_mid = new TH1D("hntrig_mixed_mid","",200,0,200);
	TH1D *hntrig_mixed_mb = new TH1D("hntrig_mixed_mb","",200,0,200);
	TH1D *hntrig_mixed_fwd = new TH1D("hntrig_mixed_fwd","",nmult_fwd,multcut_fwd);

	TH2D *heta_pt_all = new TH2D("heta_pt_all","",100,-5,5,100,0,10);
	TH2D *heta_pt_ana = new TH2D("heta_pt_ana","",100,-5,5,100,0,10);
	TH2D *heta_pt_all_mult_fwd[nmult_fwd];
	for (int im=0; im<nmult_fwd; im++){
		heta_pt_all_mult_fwd[im] = new TH2D(Form("heta_pt_all_mult_fwd_%d",im),"",100,-5,5,100,0,10);
	}

	vector<float> vec_mid_pt[nmult];
	vector<float> vec_mid_phi[nmult];
	vector<float> vec_mid_eta[nmult];

	vector<float> vec_mb_pt[nmult];
	vector<float> vec_mb_phi[nmult];
	vector<float> vec_mb_eta[nmult];

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

		TTree *Tpp = (TTree*)infile->Get("Tpp");
		Tpp->SetBranchAddress("npp",&npp);
		Tpp->SetBranchAddress("pp_x",pp_x);
		Tpp->SetBranchAddress("pp_y",pp_y);
		Tpp->SetBranchAddress("pp_z",pp_z);

		int nentries = T->GetEntries();
		int nentries_pp = Tpp->GetEntries();

		if ( nentries!=nentries_pp ){
			cout << "Inconsistent number of entries" << endl;
			delete infile;
			continue;
		}

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);
			Tpp->GetEntry(ien);

			//calculate participant plane
			float cmx = 0;
			float cmy = 0;
			for (int ip=0; ip<npp; ip++){
				cmx += pp_x[ip];
				cmy += pp_y[ip];
			}
			cmx /= float(npp);
			cmy /= float(npp);

			float qx2 = 0, qy2 = 0;
			float sum_r2 = 0;
			for (int ip=0; ip<npp; ip++){
				float r2 = (pp_x[ip]-cmx)*(pp_x[ip]-cmx) + (pp_y[ip]-cmy)*(pp_y[ip]-cmy);
				float phi = atan2(pp_y[ip]-cmy, pp_x[ip]-cmx);

				qx2 += r2*cos(2*phi);
				qy2 += r2*sin(2*phi);
				sum_r2 += r2;
			}
			qx2 /= float(npp);
			qy2 /= float(npp);
			sum_r2 /= float(npp);

			float psi2 = (atan2(qy2, qx2) + const_pi)/2;
			float e2 = sqrt(qx2*qx2 + qy2*qy2) / sum_r2;
			hpp->Fill(psi2);

			//calculate event multiplicity
			hmult_all->Fill(npart);

			int nmult_mid = 0, nmult_fwd = 0, nmult_bwd = 0;

			for (int ip=0; ip<npart; ip++){
				if ( fabs(part_eta[ip])<eta_max && part_pt[ip]>pt_min ){
					nmult_mid++;
				}
				if ( fabs(part_eta[ip])>3.3 && fabs(part_eta[ip])<5.0 ){
					if ( part_eta[ip]<0 ){
						nmult_fwd++;
					}else{
						nmult_bwd++; //opposite def
					}
				}
				heta_pt_all->Fill(part_eta[ip], part_pt[ip]);
			}

			hmult_mid->Fill(nmult_mid);
			hmult_fwd->Fill(nmult_fwd);

			bool bMB = (nmult_fwd>0 && nmult_bwd>0) ? true : false;
			if ( bMB ){
				hmult_mid_mb->Fill(nmult_mid);
			}

			int imult_mid = int(nmult_mid/dmult);
			if ( imult_mid>=nmult ) continue;

			h2d_ecc2_mult_fwd->Fill(nmult_fwd,e2);
			h2d_ecc2_sq_mult_fwd->Fill(nmult_fwd,e2*e2);
			hprof_ecc2_mult_fwd->Fill(nmult_fwd,e2);
			hprof_ecc2_sq_mult_fwd->Fill(nmult_fwd,e2*e2);

			h2d_ecc2_mult_mid->Fill(nmult_mid,e2);
			h2d_ecc2_sq_mult_mid->Fill(nmult_mid,e2*e2);
			hprof_ecc2_mult_mid->Fill(nmult_mid,e2);
			hprof_ecc2_sq_mult_mid->Fill(nmult_mid,e2*e2);

			int imult_fwd = hntrig_fwd->FindBin(nmult_fwd) - 1; 

			for (int ip=0; ip<npart; ip++){

				heta_pt_all_mult_fwd[imult_fwd]->Fill(part_eta[ip], part_pt[ip]);

				if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
				if ( fabs(part_eta[ip])>eta_max ) continue;

				heta_pt_ana->Fill(part_eta[ip], part_pt[ip]);
				hntrig_mid->Fill(nmult_mid);
				hntrig_fwd->Fill(nmult_fwd);

				if ( bMB ){
					hntrig_mb->Fill(nmult_mid);
				}

				float v2pp = cos(2*(part_phi[ip] - psi2));
				int v2pp_sign = (v2pp>0) ? 1 : -1;

				hv2pp_mid->Fill(nmult_mid, v2pp);
				hv2pp_sq_mid->Fill(nmult_mid, v2pp_sign*v2pp*v2pp);
				hv2pp_fwd->Fill(nmult_fwd, v2pp);
				hv2pp_sq_fwd->Fill(nmult_fwd, v2pp_sign*v2pp*v2pp);

				for (int jp=0; jp<npart; jp++){
					if ( ip==jp ) continue;
					if ( part_pt[jp]<pt_min || part_pt[jp]>pt_max ) continue;
					if ( fabs(part_eta[jp])>eta_max ) continue;

					float dphi = part_phi[jp] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = part_eta[jp] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_mid_same_dphi_deta[imult_mid]->Fill(dphi, deta);
					h2d_fwd_same_dphi_deta[imult_fwd]->Fill(dphi, deta);

					if ( bMB ){
						h2d_mb_same_dphi_deta[imult_mid]->Fill(dphi, deta);
					}
				}//jp

				if ( vec_mid_pt[imult_mid].size()>0 ){
					hntrig_mixed_mid->Fill(nmult_mid);
				}

				if ( vec_fwd_pt[imult_fwd].size()>0 ){
					hntrig_mixed_fwd->Fill(nmult_fwd);
				}

				//mixed events mid
				for (unsigned int jj=0; jj<vec_mid_pt[imult_mid].size(); jj++){
					float dphi = vec_mid_phi[imult_mid][jj] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_mid_eta[imult_mid][jj] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_mid_mixed_dphi_deta[imult_mid]->Fill(dphi, deta);
				}//jj

				//mixed events fwd 
				for (unsigned int jj=0; jj<vec_fwd_pt[imult_fwd].size(); jj++){
					float dphi = vec_fwd_phi[imult_fwd][jj] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_fwd_eta[imult_fwd][jj] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_fwd_mixed_dphi_deta[imult_fwd]->Fill(dphi, deta);
				}//jj

				if ( bMB ){
					if ( vec_mb_pt[imult_mid].size()>0 ){
						hntrig_mixed_mb->Fill(nmult_mid);
					}

					//mixed events mid
					for (unsigned int jj=0; jj<vec_mb_pt[imult_mid].size(); jj++){
						float dphi = vec_mb_phi[imult_mid][jj] - part_phi[ip];
						if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
						else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

						float deta = vec_mb_eta[imult_mid][jj] - part_eta[ip];
						if ( fabs(deta)>deta_max ) continue;

						h2d_mb_mixed_dphi_deta[imult_mid]->Fill(dphi, deta);
					}//jj
				}//bMB

			}//ip


			//fill mixed event pool
			vec_mid_pt[imult_mid].clear();
			vec_mid_phi[imult_mid].clear();
			vec_mid_eta[imult_mid].clear();

			vec_fwd_pt[imult_fwd].clear();
			vec_fwd_phi[imult_fwd].clear();
			vec_fwd_eta[imult_fwd].clear();

			for (int ip=0; ip<npart; ip++){
				if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
				if ( fabs(part_eta[ip])>eta_max ) continue;
				vec_mid_pt[imult_mid].push_back(part_pt[ip]);
				vec_mid_phi[imult_mid].push_back(part_phi[ip]);
				vec_mid_eta[imult_mid].push_back(part_eta[ip]);

				vec_fwd_pt[imult_fwd].push_back(part_pt[ip]);
				vec_fwd_phi[imult_fwd].push_back(part_phi[ip]);
				vec_fwd_eta[imult_fwd].push_back(part_eta[ip]);
			}

			if ( bMB ){
				vec_mb_pt[imult_mid].clear();
				vec_mb_phi[imult_mid].clear();
				vec_mb_eta[imult_mid].clear();

				for (int ip=0; ip<npart; ip++){
					if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
					if ( fabs(part_eta[ip])>eta_max ) continue;
					vec_mb_pt[imult_mid].push_back(part_pt[ip]);
					vec_mb_phi[imult_mid].push_back(part_phi[ip]);
					vec_mb_eta[imult_mid].push_back(part_eta[ip]);
				}//
			}//bMB

		}//ien


		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	for (int im=0; im<nmult; im++){
		h2d_mid_same_dphi_deta[im]->Write();
		h2d_mid_mixed_dphi_deta[im]->Write();
		h2d_mb_same_dphi_deta[im]->Write();
		h2d_mb_mixed_dphi_deta[im]->Write();
	}

	for (int im=0; im<nmult_fwd; im++){
		h2d_fwd_same_dphi_deta[im]->Write();
		h2d_fwd_mixed_dphi_deta[im]->Write();
	}

	heta_pt_all->Write();
	heta_pt_ana->Write();
	for (int im=0; im<nmult_fwd; im++){
		heta_pt_all_mult_fwd[im]->Write();
	}

	h2d_ecc2_mult_fwd->Write();
	h2d_ecc2_sq_mult_fwd->Write();
	hprof_ecc2_mult_fwd->Write();
	hprof_ecc2_sq_mult_fwd->Write();

	h2d_ecc2_mult_mid->Write();
	h2d_ecc2_sq_mult_mid->Write();
	hprof_ecc2_mult_mid->Write();
	hprof_ecc2_sq_mult_mid->Write();

	hv2pp_mid->Write();
	hv2pp_sq_mid->Write();
	hv2pp_fwd->Write();
	hv2pp_sq_fwd->Write();

	hmult_all->Write();
	hmult_mid->Write();
	hmult_mid_mb->Write();
	hmult_fwd->Write();
	//hprof_mult_mid->Write();
	//hprof_mult_fwd->Write();
	hntrig_mid->Write();
	hntrig_fwd->Write();
	hntrig_mb->Write();
	hntrig_mixed_mid->Write();
	hntrig_mixed_fwd->Write();
	hntrig_mixed_mb->Write();
	outfile->Close();

}

//ALICE acceptance w/o MPI

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

void make_hist_pp13TeV_01(const char *fname="file.lst"){

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
	const float pt_min = 1.0;
	const float pt_max = 2.0;
	const float eta_max = 0.9;
	const float deta_max = 1.8;
	
	const float deta_cut = 2.0;
	const float const_pi = TMath::Pi();

	//const int nmult = 5;
	//const float cut_mult[nmult+1] = {100, 60, 20, 5, 1, 0};
	const int nmult = 1;

	const int npt = 6;

	TFile *infile_ref = new TFile("/alice/home/shlim/work/Pythia/jobs/outfile_hist_pp13TeV_set00_grp000_try000.root","read");
	TH2D *href_mid_fwd = (TH2D*)infile_ref->Get("hevent_mult_mid_fwd")->Clone("href_mid_fwd");
	TH1D *href_fwd = (TH1D*)href_mid_fwd->ProjectionY("href_fwd");

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	float f_scale;
	float f_bMPI;
	int i_nMPI;

	int i_np;
	int i_p_id[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];

	int i_njet;
	float f_jet_pt[100], f_jet_eta[100], f_jet_phi[100];

	TH1D *hntrig_same[nmult];
	TH1D *hntrig_mixed[nmult];

	TH2D *h2d_same_dphi_deta[nmult][npt];
	TH2D *h2d_mixed_dphi_deta[nmult][npt];

	for (int im=0; im<nmult; im++){
		hntrig_same[im] = new TH1D(Form("hntrig_same_mult%02d",im),"",100,0,100);
		hntrig_mixed[im] = new TH1D(Form("hntrig_mixed_mult%02d",im),"",100,0,100);

		for (int ipt=0; ipt<npt; ipt++){
			h2d_same_dphi_deta[im][ipt] = new TH2D(Form("h2d_same_dphi_deta_mult%02d_ptlead%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
			h2d_mixed_dphi_deta[im][ipt] = new TH2D(Form("h2d_mixed_dphi_deta_mult%02d_ptlead%02d",im,ipt),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		}
	}

	TH2D *h2d_jet_dphi_deta = new TH2D("h2d_jet_dphi_deta","",96,-const_pi/2,3*const_pi/2,100,-5,5);

	TH1D *hjet_ptlead = new TH1D("hjet_ptlead","",100,0,500);

	TH1D *hevent_Q = new TH1D("hevent_Q","",1500,0,150);
	TH1D *hevent_nMPI = new TH1D("hevent_nMPI","",51,-0.5,50.5);
	TH1D *hevent_bMPI = new TH1D("hevent_bMPI","",300,0,3);
	TH2D *hevent_mult_mid_fwd = new TH2D("hevent_mult_mid_fwd","",200,0,200,200,0,200);
	TH2D *hevent_mult_fwd = new TH2D("hevent_mult_fwd","",nmult,0,nmult,200,0,200);
	TH2D *hevent_mult_mid = new TH2D("hevent_mult_mid","",nmult,0,nmult,200,0,200);

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

	vector<float> vec_pt;
	vector<float> vec_phi;
	vector<float> vec_eta;

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

		T->SetBranchAddress("njet",&i_njet);
		T->SetBranchAddress("jet_eta",f_jet_eta);
		T->SetBranchAddress("jet_phi",f_jet_phi);
		T->SetBranchAddress("jet_pt",f_jet_pt);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			if ( i_njet==2 && fabs(f_jet_eta[0])<0.4 ){

				float dphi = f_jet_phi[1] - f_jet_phi[0];
				if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
				else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

				float deta = f_jet_eta[1] - f_jet_eta[0];

				h2d_jet_dphi_deta->Fill(dphi, deta);
			}

			//continue;

			int njet_good = 0;
			float pt_lead = 0; //leading jet pT
			for (int ijet=0; ijet<i_njet; ijet++){
				if ( f_jet_pt[ijet]>10.0 && fabs(f_jet_eta[ijet])<0.4 ){
				//if ( f_jet_pt[ijet]>10.0 && fabs(f_jet_eta[ijet])<1.0 ){
				//if ( f_jet_pt[ijet]>10.0 && fabs(f_jet_eta[ijet])<0.1 ){
					njet_good++;

					if ( f_jet_pt[ijet]>pt_lead ){
						pt_lead = f_jet_pt[ijet];
					}
				}
			}//

			if ( njet_good<1 ) continue;

			int nmult_mid = 0, nmult_fwd = 0;

			//count multiplicity for bin
			for (int ip=0; ip<i_np; ip++){
				if ( fabs(f_p_eta[ip])<eta_max && f_p_pt[ip]>pt_min ){
					nmult_mid++;
				}

				if ( (f_p_eta[ip]>2.8 && f_p_eta[ip]<5.1) || (f_p_eta[ip]>-3.7 && f_p_eta[ip]<-1.7) ){
					nmult_fwd++;
				}
			}

			//centrality check
			float centrality = href_fwd->Integral(1,href_fwd->FindBin(nmult_fwd))/href_fwd->Integral()*100;
			int ind_mult = -1;

			/*
			for (int im=0; im<nmult; im++){
				if ( centrality>(100-cut_mult[im]) && centrality<=(100-cut_mult[im+1]) ){
					ind_mult = im;
					break;
				}
			}//im
			*/

			if ( centrality<90.0 ) continue;
			ind_mult = 0;

			hjet_ptlead->Fill(pt_lead);

			//continue;

			hevent_Q->Fill(f_scale);
			hevent_nMPI->Fill(i_nMPI);
			hevent_bMPI->Fill(f_bMPI);


			hevent_mult_mid_fwd->Fill(nmult_mid, nmult_fwd);

			//continue;

			hevent_mult_mid->Fill(ind_mult+0.5, nmult_mid); 
			hevent_mult_fwd->Fill(ind_mult+0.5, nmult_fwd); 

			/*
				 if ( ien<1000 )
				 cout << ind_mult << " " << nmult_fwd << " " << centrality << endl;
				 */
			//continue;

			/*
			float pt_lead = 0.0;

			for (int ip=0; ip<i_np; ip++){
				if ( fabs(f_p_eta[ip])>eta_max ) continue;
				if ( f_p_pt[ip]>pt_lead ){
					pt_lead = f_p_pt[ip];
				}
			}
			*/

			/*
			hevent_Q_mult[ind_mult]->Fill(pt_lead, f_scale);
			hevent_bMPI_mult[ind_mult]->Fill(pt_lead, f_bMPI);
			hevent_nMPI_mult[ind_mult]->Fill(pt_lead, i_nMPI);
			*/

			hevent_Q_bMPI_mult[ind_mult]->Fill(f_scale, f_bMPI);
			hevent_Q_nMPI_mult[ind_mult]->Fill(f_scale, i_nMPI);

			//continue;

			for (int ip=0; ip<i_np; ip++){

				//same events
				if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;
				if ( fabs(f_p_eta[ip])>eta_max ) continue;

				hntrig_same[ind_mult]->Fill(pt_lead);

				//continue;

				for (int jp=0; jp<i_np; jp++){
					if ( ip==jp ) continue;
					if ( f_p_pt[jp]<pt_min || f_p_pt[jp]>pt_max ) continue;
					if ( fabs(f_p_eta[jp])>eta_max ) continue;

					float dphi = f_p_phi[jp] - f_p_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = f_p_eta[jp] - f_p_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					if ( pt_lead>10 ) h2d_same_dphi_deta[0][0]->Fill(dphi, deta);
					if ( pt_lead>20 ) h2d_same_dphi_deta[0][1]->Fill(dphi, deta);
					if ( pt_lead>30 ) h2d_same_dphi_deta[0][2]->Fill(dphi, deta);
					if ( pt_lead>40 ) h2d_same_dphi_deta[0][3]->Fill(dphi, deta);
					if ( pt_lead>50 ) h2d_same_dphi_deta[0][4]->Fill(dphi, deta);
					if ( pt_lead>60 ) h2d_same_dphi_deta[0][5]->Fill(dphi, deta);
				}//jp


				//mixed events
				if ( vec_pt.size()>0 ){
					hntrig_mixed[0]->Fill(pt_lead);
				}

				for (unsigned int jj=0; jj<vec_pt.size(); jj++){
					//if ( vec_pt[ind_mult][jj]<pt_min || vec_pt[ind_mult][jj]>pt_max ) continue;
					//if ( fabs(vec_eta[ind_mult][jj])>eta_max ) continue;

					float dphi = vec_phi[jj] - f_p_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_eta[jj] - f_p_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					if ( pt_lead>10 ) h2d_mixed_dphi_deta[0][0]->Fill(dphi, deta);
					if ( pt_lead>20 ) h2d_mixed_dphi_deta[0][1]->Fill(dphi, deta);
					if ( pt_lead>30 ) h2d_mixed_dphi_deta[0][2]->Fill(dphi, deta);
					if ( pt_lead>40 ) h2d_mixed_dphi_deta[0][3]->Fill(dphi, deta);
					if ( pt_lead>50 ) h2d_mixed_dphi_deta[0][4]->Fill(dphi, deta);
					if ( pt_lead>60 ) h2d_mixed_dphi_deta[0][5]->Fill(dphi, deta);
				}//jj

			}//


			//fill mixed event pool
			vec_pt.clear();
			vec_phi.clear();
			vec_eta.clear();

			for (int ip=0; ip<i_np; ip++){
				if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;
				if ( fabs(f_p_eta[ip])>eta_max ) continue;

				vec_pt.push_back(f_p_pt[ip]);
				vec_phi.push_back(f_p_phi[ip]);
				vec_eta.push_back(f_p_eta[ip]);
			}

		}//ien

		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");


	for (int im=0; im<nmult; im++){
		for (int ipt=0; ipt<npt; ipt++){
			h2d_same_dphi_deta[im][ipt]->Write();
			h2d_mixed_dphi_deta[im][ipt]->Write();
		}
	}

	for (int im=0; im<nmult; im++){
		hntrig_same[im]->Write();
		hntrig_mixed[im]->Write();
	}

	/*
	for (int im=0; im<nmult; im++){
		hevent_Q_mult[im]->Write();
		hevent_bMPI_mult[im]->Write();
		hevent_nMPI_mult[im]->Write();

		hevent_Q_bMPI_mult[im]->Write();
		hevent_Q_nMPI_mult[im]->Write();
	}
	*/

	hjet_ptlead->Write();


	hevent_Q->Write();
	hevent_nMPI->Write();
	hevent_bMPI->Write();

	hevent_mult_mid_fwd->Write();
	hevent_mult_mid->Write();
	hevent_mult_fwd->Write();

	h2d_jet_dphi_deta->Write();

	outfile->Close();

}

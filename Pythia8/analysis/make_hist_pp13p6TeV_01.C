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

void make_hist_pp13p6TeV_01(const char *dataset="pp13p6TeV", int set=30){

	//const int nmult = 20;

	ifstream flist;
	char ffname[300];
	sprintf(ffname, "file_%s_set%d.lst",dataset,set);
	flist.open(ffname);


	int i_nmult_v0a;
	int i_nmult_v0c;

	int i_np;
	int i_p_id[10000];
	int i_p_status[10000];
	float f_p_pt[10000];
	float f_p_eta[10000];
	float f_p_phi[10000];
	float f_p_m[10000];
	float f_p_vx[10000];
	float f_p_vy[10000];
	float f_p_vz[10000];
	float f_p_vt[10000];

	float ptbin[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0};
	const int nptbin = sizeof(ptbin)/sizeof(float) - 1;

	cout << "nptbin: " << nptbin << endl;

	TH1D *hpt_3322 = new TH1D("hpt_3322","",nptbin,ptbin);
	TH1D *hpt_3312 = new TH1D("hpt_3312","",nptbin,ptbin);
	TH1D *hpt_3324 = new TH1D("hpt_3324","",nptbin,ptbin);
	TH1D *hpt_3314 = new TH1D("hpt_3314","",nptbin,ptbin);

	TH1D *hrap = new TH1D("hrap","",100,-5,5);

	TH1D *hmult_fwd = new TH1D("hmult_fwd","",500,0,500);

	/*
	TH1D *hmult_all = new TH1D("hmult_all","",500,0,500);
	TH1D *hmult_mid = new TH1D("hmult_mid","",500,0,500);
	TH1D *hntrig_mid = new TH1D("hntrig_mid","",200,0,200);
	TH1D *hntrig_mixed_mid = new TH1D("hntrig_mixed_mid","",200,0,200);

	TH2D *heta_pt_all = new TH2D("heta_pt_all","",100,-5,5,100,0,10);
	TH2D *heta_pt_ana = new TH2D("heta_pt_ana","",100,-5,5,100,0,10);
	*/

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");
		//T->SetBranchAddress("nmult_v0a",&i_nmult_v0a);
		//T->SetBranchAddress("nmult_v0c",&i_nmult_v0c);

		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id", i_p_id);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_pt", f_p_pt);
		//T->SetBranchAddress("p_m",  f_p_m);
		T->SetBranchAddress("p_vt", f_p_vt);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			bool b_inelgt0 = false;

			//CHECK INEL>0
			for (int ip=0; ip<i_np; ip++){
				if ( fabs(f_p_eta[ip])<1.0 ){
					b_inelgt0 = true;
					break;
				}
			}

			if ( !b_inelgt0 ) continue;

			//hmult_fwd->Fill(i_nmult_v0a + i_nmult_v0c);

			for (int ip=0; ip<i_np; ip++){

				if ( f_p_vt[ip]>0.1 ) continue;
				if ( !(abs(i_p_id[ip])==3322 || abs(i_p_id[ip])==3312 || abs(i_p_id[ip])==3324 || abs(i_p_id[ip])==3314) ) continue;

				float mass = 1.53;
				if ( abs(i_p_id[ip])==3322 ){
					mass  = 1.31;
				}else if ( abs(i_p_id[ip])==3312 ){
					mass  = 1.32;
				}

				TLorentzVector lvec;
				lvec.SetPtEtaPhiM(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], mass);

				hrap->Fill(lvec.Rapidity());

				if ( fabs(lvec.Rapidity())>0.5 ) continue;

				if ( abs(i_p_id[ip])==3322 ){
					hpt_3322->Fill(f_p_pt[ip]);
				}else if ( abs(i_p_id[ip])==3312 ){
					hpt_3312->Fill(f_p_pt[ip]);
				}else if ( abs(i_p_id[ip])==3324 ){
					hpt_3324->Fill(f_p_pt[ip]);
				}else if ( abs(i_p_id[ip])==3314 ){
					hpt_3314->Fill(f_p_pt[ip]);
				}


			}//

		}//ien

		infile->Close();
		delete infile;

	}//

	TFile *outfile = new TFile(Form("outfile_hist_%s_set%d_00.root",dataset,set),"recreate");
	//TFile *outfile = new TFile("outfile_hist_pp13p0TeV_00.root","recreate");
	//TFile *outfile = new TFile("outfile_hist_pp2p76TeV_00.root","recreate");

	//hmult_fwd->Write();

	hrap->Write();
	hpt_3322->Write();
	hpt_3312->Write();
	hpt_3324->Write();
	hpt_3314->Write();

	outfile->Close();

}

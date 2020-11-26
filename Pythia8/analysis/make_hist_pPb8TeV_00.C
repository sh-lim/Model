//ATLAS vn-pt correlation

#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>
#include <TComplex.h>

using namespace std;

float calc2_2subevent(float Xn_a, float Yn_a, float M_a, float Xn_b, float Yn_b, float M_b)
{
	if ( M_a<1 || M_b<1 ) return -9999;
	float numerator = Xn_a*Xn_b + Yn_a*Yn_b;
	float denominator = M_a*M_b;
	return numerator/denominator;
}

float calc4_2subevent(
		float Xn_a, float Yn_a, float X2n_a, float Y2n_a, float M_a, 
		float Xn_b, float Yn_b, float X2n_b, float Y2n_b, float M_b
		)
{
	if ( M_a<2 || M_b<2 ) return -9999;

	TComplex tcn_a(Xn_a, Yn_a);
	TComplex tcn_b(Xn_b, Yn_b);
	TComplex tc2n_a(X2n_a, Y2n_a);
	TComplex tc2n_b(X2n_b, Y2n_b);

	TComplex tc_one = tcn_a*tcn_a - tc2n_a;
	TComplex tc_two = TComplex::Conjugate(tcn_b*tcn_b - tc2n_b);
	TComplex tc_numerator = tc_one*tc_two;

	float numerator = tc_numerator.Re();
	float denominator = M_a*(M_a-1)*M_b*(M_b-1);
	return numerator/denominator;
}

void make_hist_pPb8TeV_00(const char *fname="file.lst"){

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
	const float pt_min = 0.3;
	const float pt_max = 2.0;
	const float eta_max = 0.9;
	const float deta_max = 1.8;
	
	const float deta_cut = 2.0;
	const float const_pi = TMath::Pi();

	const bool bPRE = false;

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int i_np;
	int i_p_id[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000], f_p_vt[2000];

	TH2D *heta_pt_ana = new TH2D("heta_pt_ana","",100,-5,5,100,0,10);
	TH1D *hevent_mult_mid = new TH1D("hevent_mult_mid","",250,0,250);

	TProfile *hprof_meanpt = new TProfile("hprof_meanpt","",10,0,200);
	TProfile *hprof_ck = new TProfile("hprof_ck","",10,0,200);
	TProfile *hprof_corr2 = new TProfile("hprof_corr2","",10,0,200);
	TProfile *hprof_corr4 = new TProfile("hprof_corr4","",10,0,200);
	TProfile *hprof_ptcorr2 = new TProfile("hprof_ptcorr2","",10,0,200);

	TFile *infile_ref = new TFile("/alice/home/shlim/work/Pythia/jobs/outfile_hist_pPb8TeV_set00_grp000_try102.root","read");
	TProfile *hprof_ref_meanpt = (TProfile*)infile_ref->Get("hprof_meanpt");

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id",i_p_id);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_pt",f_p_pt);
		T->SetBranchAddress("p_vt",f_p_vt);

		int nentries = T->GetEntries();

		if ( bPRE ){
			for (int ien=0; ien<nentries; ien++){
				T->GetEntry(ien);

				int nmult_mid = 0;

				//count multiplicity for bin
				for (int ip=0; ip<i_np; ip++){
					if ( 
							fabs(f_p_eta[ip])<2.5 
							&& f_p_vt[ip]<0.1 
							&& f_p_pt[ip]>0.5 
							&& f_p_pt[ip]<5.0 
						 ){
						++nmult_mid;
					}
				}

				hevent_mult_mid->Fill(nmult_mid);

				float B_sumpt = 0.0;
				int B_ntrk = 0;

				for (int ip=0; ip<i_np; ip++){

					if ( f_p_vt[ip]>0.1 ) continue;
					if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;
					if ( fabs(f_p_eta[ip])>0.5 ) continue;

					B_sumpt += f_p_pt[ip];
					++B_ntrk;

				}

				if ( B_ntrk>0 ){
					hprof_meanpt->Fill(float(nmult_mid), B_sumpt/B_ntrk);
				}

			}//ien

			continue;
		}

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			int nmult_mid = 0;

			//count multiplicity for bin
			for (int ip=0; ip<i_np; ip++){
				if ( 
						fabs(f_p_eta[ip])<2.5 
						&& f_p_vt[ip]<0.1 
						&& f_p_pt[ip]>0.5 
						&& f_p_pt[ip]<5.0 
						){
					++nmult_mid;
				}
			}//ip

			float B_sumpt = 0.0;
			int B_ntrk = 0;

			for (int ip=0; ip<i_np; ip++){

				if ( f_p_vt[ip]>0.1 ) continue;
				if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;
				if ( fabs(f_p_eta[ip])>0.5 ) continue;

				B_sumpt += f_p_pt[ip];
				++B_ntrk;
			}

			float B_meanpt = hprof_ref_meanpt->GetBinContent(hprof_meanpt->FindBin(nmult_mid+0.01));
			float B_sumck = 0.0;
			int B_npair = 0;

			float Qx2_A = 0.0, Qy2_A = 0.0, Qx4_A = 0.0, Qy4_A = 0.0, Qw_A = 0.0;
			float Qx2_C = 0.0, Qy2_C = 0.0, Qx4_C = 0.0, Qy4_C = 0.0, Qw_C = 0.0;
			float Px2_A = 0.0, Py2_A = 0.0, Pw_A = 0.0;
			float Px2_C = 0.0, Py2_C = 0.0, Pw_C = 0.0;

			for (int ip=0; ip<i_np; ip++){

				if ( f_p_vt[ip]>0.1 ) continue;
				if ( f_p_pt[ip]<pt_min || f_p_pt[ip]>pt_max ) continue;

				if ( f_p_eta[ip]>-2.5 && f_p_eta[ip]<-0.75 ){
					Qx2_A += cos(2*f_p_phi[ip]);
					Qy2_A += sin(2*f_p_phi[ip]);
					Qx4_A += cos(4*f_p_phi[ip]);
					Qy4_A += sin(4*f_p_phi[ip]);
					++Qw_A;

					if ( B_ntrk>0 ){
						Px2_A += cos(2*f_p_phi[ip])*(B_sumpt/B_ntrk - B_meanpt);
						Py2_A += sin(2*f_p_phi[ip])*(B_sumpt/B_ntrk - B_meanpt);
						++Pw_A;
					}

				}else if ( f_p_eta[ip]>0.75 && f_p_eta[ip]<2.5 ){
					Qx2_C += cos(2*f_p_phi[ip]);
					Qy2_C += sin(2*f_p_phi[ip]);
					Qx4_C += cos(4*f_p_phi[ip]);
					Qy4_C += sin(4*f_p_phi[ip]);
					++Qw_C;

					if ( B_ntrk>0 ){
						Px2_C += cos(2*f_p_phi[ip])*(B_sumpt/B_ntrk - B_meanpt);
						Py2_C += sin(2*f_p_phi[ip])*(B_sumpt/B_ntrk - B_meanpt);
						++Pw_C;
					}
				}

				if ( fabs(f_p_eta[ip])>0.5 ) continue;

				for (int jp=0; jp<i_np; jp++){

					if ( ip==jp ) continue;
					if ( f_p_vt[jp]>0.1 ) continue;
					if ( f_p_pt[jp]<pt_min || f_p_pt[jp]>pt_max ) continue;
					if ( fabs(f_p_eta[jp])>0.5 ) continue;

					B_sumck += (B_meanpt - f_p_pt[ip])*(B_meanpt - f_p_pt[jp]);
					++B_npair;

				}//jp
			}//ip

			if ( B_npair>0 ){
				hprof_ck->Fill(float(nmult_mid), B_sumck/B_npair);
			}

			float corr2 = calc2_2subevent(Qx2_A,Qy2_A,Qw_A,Qx2_C,Qy2_C,Qw_C);
			float corr4 = calc4_2subevent(
					Qx2_A,Qy2_A,Qx4_A,Qy4_A,Qw_A,
					Qx2_C,Qy2_C,Qx4_C,Qy4_C,Qw_C
					);

			float ptcorr2 = calc2_2subevent(Px2_A,Py2_A,Pw_A,Px2_C,Py2_C,Pw_C);

			if ( corr2>-1 ){
				hprof_corr2->Fill(float(nmult_mid), corr2);
			}

			if ( corr4>-1 ){
				hprof_corr4->Fill(float(nmult_mid), corr4);
			}

			if ( ptcorr2>-1 ){
				hprof_ptcorr2->Fill(float(nmult_mid), ptcorr2);
			}

		}//ien

		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	hevent_mult_mid->Write();
	hprof_meanpt->Write();
	hprof_ck->Write();
	hprof_corr2->Write();
	hprof_corr4->Write();
	hprof_ptcorr2->Write();

	outfile->Close();

}

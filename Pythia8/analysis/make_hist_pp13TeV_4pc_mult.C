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

float calc2_event(float Xn, float Yn, float M)
{
	if ( M<2 ) return -9999;
	float numerator = Xn*Xn + Yn*Yn - M;
	float denominator = M*(M-1);
	return numerator/denominator;
}

float calc2_2subevent(float Xn_a, float Yn_a, float M_a, float Xn_b, float Yn_b, float M_b)
{
	if ( M_a<1 || M_b<1 ) return -9999;
	float numerator = Xn_a*Xn_b + Yn_a*Yn_b;
	float denominator = M_a*M_b;
	return numerator/denominator;
}

float calc2_3subevent(float Xn_a, float Yn_a, float M_a, float Xn_b, float Yn_b, float M_b)
{
	if ( M_a<1 || M_b<1 ) return -9999;
	float numerator = Xn_a*Xn_b + Yn_a*Yn_b;
	float denominator = M_a*M_b;
	return numerator/denominator;
}

float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

	if ( M<4 ) return -9999;

	float Qn2 = Xn*Xn+Yn*Yn;
	float Qn2d = Xn*Xn-Yn*Yn;

	float one   = Qn2*Qn2;
	float two   = X2n*X2n+Y2n*Y2n;
	float three = (2*(X2n*Qn2d + 2*Y2n*Xn*Yn));
	float four  = 2*(2*(M-2)*Qn2);
	float five  = 2*M*(M-3);

	float numerator = one + two - three - four + five;
	float denominator = M*(M-1)*(M-2)*(M-3);

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

float calc4_3subevent(
		float Xn_a, float Yn_a, float X2n_a, float Y2n_a, float M_a, 
		float Xn_b, float Yn_b, float X2n_b, float Y2n_b, float M_b,
		float Xn_c, float Yn_c, float X2n_c, float Y2n_c, float M_c
		)
{
	if ( M_a<2 || M_b<1 || M_c<1 ) return -9999;

	TComplex tcn_a(Xn_a, Yn_a);
	TComplex tcn_b(Xn_b, Yn_b);
	TComplex tcn_c(Xn_c, Yn_c);
	TComplex tc2n_a(X2n_a, Y2n_a);
	TComplex tc2n_b(X2n_b, Y2n_b);
	TComplex tc2n_c(X2n_c, Y2n_c);

	TComplex tc_one = tcn_a*tcn_a - tc2n_a;
	TComplex tc_two = TComplex::Conjugate(tcn_b)*TComplex::Conjugate(tcn_c);
	TComplex tc_numerator = tc_one*tc_two;

	float numerator = tc_numerator.Re();
	float denominator = M_a*(M_a-1)*M_b*M_c;
	return numerator/denominator;
}

void make_hist_pp13TeV_4pc_mult(const char *fname="file.lst"){

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
	const double multcut_mid[nmult+1] = {0, 10, 20, 30, 85, 500};
	const double multcut_fwd[nmult+1] = {0, 15, 25, 35, 85, 500};

	//const int npt = 8;
	//const float ptcut[npt+1] = {0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};
	const int npt = 12;
	const double ptcut[npt+1] = 
	{
		0.5, 1.0, 1.5, 2.0, 2.5,
		3.0, 4.0, 5.0, 7.5, 10.0,
		15.0, 20.0, 30.0
	};

	const int nptmax = 3;
	const double ptmax[nptmax] = {2.5, 5.0, 1000.0};

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int npart;
	int part_pid[2000];
	float part_eta[2000], part_phi[2000], part_pt[2000];

	TProfile *hprof_mid_4pc[nmult];
	TProfile *hprof_fwd_4pc[nmult];
	TProfile *hprof_mid_4pc_ref[nmult];
	TProfile *hprof_fwd_4pc_ref[nmult];
	TProfile *hprof_mid_2pc[nmult];
	TProfile *hprof_fwd_2pc[nmult];
	TProfile *hprof_mid_2pc_ref[nmult];
	TProfile *hprof_fwd_2pc_ref[nmult];

	TProfile *hprof_4pc_mult[nptmax];
	TProfile *hprof_2pc_mult[nptmax];
	TProfile *hprof_4pc_mult_sub[nptmax];
	TProfile *hprof_2pc_mult_sub[nptmax];
	TProfile *hprof_4pc_mult_3sub[nptmax];
	TProfile *hprof_2pc_mult_3sub_ab[nptmax];
	TProfile *hprof_2pc_mult_3sub_ac[nptmax];

	for (int icut=0; icut<nptmax; icut++){
		hprof_4pc_mult[icut] = new TProfile(Form("hprof_4pc_ptmax%d_mult",icut),"",20,0,200,-1,1);
		hprof_2pc_mult[icut] = new TProfile(Form("hprof_2pc_ptmax%d_mult",icut),"",20,0,200,-1,1);
		hprof_4pc_mult_sub[icut] = new TProfile(Form("hprof_4pc_ptmax%d_mult_sub",icut),"",20,0,200,-1,1);
		hprof_2pc_mult_sub[icut] = new TProfile(Form("hprof_2pc_ptmax%d_mult_sub",icut),"",20,0,200,-1,1);
		hprof_4pc_mult_3sub[icut] = new TProfile(Form("hprof_4pc_ptmax%d_mult_3sub",icut),"",20,0,200,-1,1);
		hprof_2pc_mult_3sub_ab[icut] = new TProfile(Form("hprof_2pc_ptmax%d_mult_3sub_ab",icut),"",20,0,200,-1,1);
		hprof_2pc_mult_3sub_ac[icut] = new TProfile(Form("hprof_2pc_ptmax%d_mult_3sub_ac",icut),"",20,0,200,-1,1);
	}//icut

	TProfile *hprof_4pc_mult_dir = new TProfile("hprof_4pc_mult_dir","",20,0,200,-1,1);
	TProfile *hprof_2pc_mult_dir = new TProfile("hprof_2pc_mult_dir","",20,0,200,-1,1);
	TProfile *hprof_4pc_mult_sub_dir = new TProfile("hprof_4pc_mult_sub_dir","",20,0,200,-1,1);
	TProfile *hprof_2pc_mult_sub_dir = new TProfile("hprof_2pc_mult_sub_dir","",20,0,200,-1,1);

	TProfile *hprof_4pc_mult_evt = new TProfile("hprof_4pc_mult_evt","",1,0,1,-1,1);
	TProfile *hprof_2pc_mult_evt = new TProfile("hprof_2pc_mult_evt","",1,0,1,-1,1);
	TProfile *hprof_4pc_mult_sub_evt = new TProfile("hprof_4pc_mult_sub_evt","",1,0,1,-1,1);
	TProfile *hprof_2pc_mult_sub_evt = new TProfile("hprof_2pc_mult_sub_evt","",1,0,1,-1,1);

	TH1D *hmult_all = new TH1D("hmult_all","",500,0,500);
	TH1D *hmult_mid = new TH1D("hmult_mid","",500,0,500);
	TH1D *hmult_fwd = new TH1D("hmult_fwd","",500,0,500);
	TH2D *hmult_mid_fwd = new TH2D("hmult_mid_fwd","",100,0,500,100,0,500);
	TH2D *hmult_cnt_fwd = new TH2D("hmult_cnt_fwd","",100,0,500,100,0,500);

	//TH2D *hntrig_mid = new TH2D("hntrig_mid","",nmult,multcut_mid,npt,ptcut);
	//TH2D *hntrig_fwd = new TH2D("hntrig_fwd","",nmult,multcut_fwd,npt,ptcut);
	TH2D *hntrig_mid = new TH2D("hntrig_mid","",nmult,multcut_mid,npt,ptcut);
	TH2D *hntrig_fwd = new TH2D("hntrig_fwd","",nmult,multcut_fwd,npt,ptcut);

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
		//for (int ien=0; ien<10000; ien++){
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

			//Q-vec calculation
			float qx2[nptmax] = {0.0}, qy2[nptmax] = {0.0};
			float qx3[nptmax] = {0.0}, qy3[nptmax] = {0.0};
			float qx4[nptmax] = {0.0}, qy4[nptmax] = {0.0};
			float qw[nptmax] = {0.0};

			//Q-vec for calculaton for 2-subevents
			float qx2_2a[nptmax] = {0.0}, qy2_2a[nptmax] = {0.0};
			float qx2_2b[nptmax] = {0.0}, qy2_2b[nptmax] = {0.0};
			float qx4_2a[nptmax] = {0.0}, qy4_2a[nptmax] = {0.0};
			float qx4_2b[nptmax] = {0.0}, qy4_2b[nptmax] = {0.0};
			float qw_2a[nptmax] = {0.0}, qw_2b[nptmax] = {0.0};

			//Q-vec for calculaton for 3-subevents
			float qx2_3a[nptmax] = {0.0}, qy2_3a[nptmax] = {0.0};
			float qx2_3b[nptmax] = {0.0}, qy2_3b[nptmax] = {0.0};
			float qx2_3c[nptmax] = {0.0}, qy2_3c[nptmax] = {0.0};
			float qx4_3a[nptmax] = {0.0}, qy4_3a[nptmax] = {0.0};
			float qx4_3b[nptmax] = {0.0}, qy4_3b[nptmax] = {0.0};
			float qx4_3c[nptmax] = {0.0}, qy4_3c[nptmax] = {0.0};
			float qw_3a[nptmax] = {0.0}, qw_3b[nptmax] = {0.0}, qw_3c[nptmax] = {0.0};

			for (int ip=0; ip<npart; ip++){
				if ( fabs(part_eta[ip])>2.5 ) continue;
				if ( part_pt[ip]<pt_min ) continue;

				for (int icut=0; icut<nptmax; icut++){

					if ( part_pt[ip]>ptmax[icut] ) continue;

					qx2[icut] += cos(2*part_phi[ip]);
					qy2[icut] += sin(2*part_phi[ip]);
					qx3[icut] += cos(3*part_phi[ip]);
					qy3[icut] += sin(3*part_phi[ip]);
					qx4[icut] += cos(4*part_phi[ip]);
					qy4[icut] += sin(4*part_phi[ip]);
					qw[icut]++;

					if ( part_eta[ip]<0 ){
						qx2_2a[icut] += cos(2*part_phi[ip]);
						qy2_2a[icut] += sin(2*part_phi[ip]);
						qx4_2a[icut] += cos(4*part_phi[ip]);
						qy4_2a[icut] += sin(4*part_phi[ip]);
						qw_2a[icut]++;
					}else{
						qx2_2b[icut] += cos(2*part_phi[ip]);
						qy2_2b[icut] += sin(2*part_phi[ip]);
						qx4_2b[icut] += cos(4*part_phi[ip]);
						qy4_2b[icut] += sin(4*part_phi[ip]);
						qw_2b[icut]++;
					}

					if ( part_eta[ip]<-0.8 ){
						qx2_3a[icut] += cos(2*part_phi[ip]);
						qy2_3a[icut] += sin(2*part_phi[ip]);
						qx4_3a[icut] += cos(4*part_phi[ip]);
						qy4_3a[icut] += sin(4*part_phi[ip]);
						qw_3a[icut]++;
					}else if ( part_eta[ip]<0.8 ){
						qx2_3b[icut] += cos(2*part_phi[ip]);
						qy2_3b[icut] += sin(2*part_phi[ip]);
						qx4_3b[icut] += cos(4*part_phi[ip]);
						qy4_3b[icut] += sin(4*part_phi[ip]);
						qw_3b[icut]++;
					}else{
						qx2_3c[icut] += cos(2*part_phi[ip]);
						qy2_3c[icut] += sin(2*part_phi[ip]);
						qx4_3c[icut] += cos(4*part_phi[ip]);
						qy4_3c[icut] += sin(4*part_phi[ip]);
						qw_3c[icut]++;
					}

				}//icut

			}//ipart

			for (int icut=0; icut<nptmax; icut++){
				float qq2 = calc2_event(qx2[icut], qy2[icut], qw[icut]);
				float qq4 = calc4_event(qx2[icut], qy2[icut], qx4[icut], qy4[icut], qw[icut]);

				if ( qq2>-1 ){
					hprof_2pc_mult[icut]->Fill(nmult_mid, qq2);
				}
				if ( qq4>-1 ){
					hprof_4pc_mult[icut]->Fill(nmult_mid, qq4);
				}

				float qq2_2sub = calc2_2subevent(qx2_2a[icut], qy2_2a[icut], qw_2a[icut], qx2_2b[icut], qy2_2b[icut], qw_2b[icut]);
				float qq4_2sub = calc4_2subevent(
						qx2_2a[icut], qy2_2a[icut], qx4_2a[icut], qy4_2a[icut], qw_2a[icut], 
						qx2_2b[icut], qy2_2b[icut], qx4_2b[icut], qy4_2b[icut], qw_2b[icut]
						);

				if ( qq2_2sub>-1 ){
					hprof_2pc_mult_sub[icut]->Fill(nmult_mid, qq2_2sub);
				}
				if ( qq4_2sub>-1 ){
					hprof_4pc_mult_sub[icut]->Fill(nmult_mid, qq4_2sub);
				}

				float qq2_3sub_ab = calc2_3subevent(qx2_3a[icut], qy2_3a[icut], qw_3a[icut], qx2_3b[icut], qy2_3b[icut], qw_3b[icut]);
				float qq2_3sub_ac = calc2_3subevent(qx2_3a[icut], qy2_3a[icut], qw_3a[icut], qx2_3c[icut], qy2_3c[icut], qw_3c[icut]);
				float qq4_3sub = calc4_3subevent(
						qx2_3a[icut], qy2_3a[icut], qx4_3a[icut], qy4_3a[icut], qw_3a[icut], 
						qx2_3b[icut], qy2_3b[icut], qx4_3b[icut], qy4_3b[icut], qw_3b[icut],
						qx2_3c[icut], qy2_3c[icut], qx4_3c[icut], qy4_3c[icut], qw_3c[icut]
						);

				if ( qq2_3sub_ab>-1 ){
					hprof_2pc_mult_3sub_ab[icut]->Fill(nmult_mid, qq2_3sub_ab);
				}
				if ( qq2_3sub_ac>-1 ){
					hprof_2pc_mult_3sub_ac[icut]->Fill(nmult_mid, qq2_3sub_ac);
				}
				if ( qq4_3sub>-1 ){
					hprof_4pc_mult_3sub[icut]->Fill(nmult_mid, qq4_3sub);
				}

			}//icut



			hprof_4pc_mult_evt->Reset();
			hprof_2pc_mult_evt->Reset();

			hprof_4pc_mult_sub_evt->Reset();
			hprof_2pc_mult_sub_evt->Reset();

			//cumulant
			//for (int ip=0; ip<npart; ip++){
			for (int ip=0; ip<0; ip++){ //skip for now
				if ( fabs(part_eta[ip])>2.5 ) continue;
				if ( part_pt[ip]<ptcut[0] || part_pt[ip]>ptcut[npt] ) continue;

				for (int jp=ip+1; jp<npart; jp++){
					if ( fabs(part_eta[jp])>2.5 ) continue;
					if ( part_pt[jp]<ptcut[0] || part_pt[jp]>ptcut[npt] ) continue;

					for (int kp=jp+1; kp<npart; kp++){
						if ( fabs(part_eta[kp])>2.5 ) continue;
						if ( part_pt[kp]<ptcut[0] || part_pt[kp]>ptcut[npt] ) continue;

						for (int lp=kp+1; lp<npart; lp++){
							if ( fabs(part_eta[lp])>2.5 ) continue;
							if ( part_pt[lp]<ptcut[0] || part_pt[lp]>ptcut[npt] ) continue;

							if (
									part_pt[ip]<pt_max 
									&& part_pt[jp]<pt_max 
									&& part_pt[kp]<pt_max 
									&& part_pt[lp]<pt_max 
									){

								//i-j-k-l
								float dphi = part_phi[ip] + part_phi[jp] - part_phi[kp] - part_phi[lp];
								hprof_mid_4pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
								hprof_fwd_4pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
								hprof_4pc_mult_evt->Fill(0.5,cos(2*dphi));

								if ( 
										(part_eta[ip]<0 && part_eta[jp]<0 && part_eta[kp]>0 && part_eta[lp]>0)
										|| (part_eta[ip]>0 && part_eta[jp]>0 && part_eta[kp]<0 && part_eta[lp]<0) 
									 ){
									hprof_4pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
								}


								//i-k-j-l
								dphi = part_phi[ip] + part_phi[kp] - part_phi[jp] - part_phi[lp];
								hprof_mid_4pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
								hprof_fwd_4pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
								hprof_4pc_mult_evt->Fill(0.5,cos(2*dphi));

								if ( 
										(part_eta[ip]<0 && part_eta[kp]<0 && part_eta[jp]>0 && part_eta[lp]>0)
										|| (part_eta[ip]>0 && part_eta[kp]>0 && part_eta[jp]<0 && part_eta[lp]<0) 
									 ){
									hprof_4pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
								}

								//i-l-j-k
								dphi = part_phi[ip] + part_phi[lp] - part_phi[jp] - part_phi[kp];
								hprof_mid_4pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
								hprof_fwd_4pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
								hprof_4pc_mult_evt->Fill(0.5,cos(2*dphi));

								if ( 
										(part_eta[ip]<0 && part_eta[lp]<0 && part_eta[jp]>0 && part_eta[kp]>0)
										|| (part_eta[ip]>0 && part_eta[lp]>0 && part_eta[jp]<0 && part_eta[kp]<0) 
									 ){
									hprof_4pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
								}

								//j-k-i-l
								dphi = part_phi[jp] + part_phi[kp] - part_phi[ip] - part_phi[lp];
								hprof_mid_4pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
								hprof_fwd_4pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
								hprof_4pc_mult_evt->Fill(0.5,cos(2*dphi));

								if ( 
										(part_eta[jp]<0 && part_eta[kp]<0 && part_eta[ip]>0 && part_eta[lp]>0)
										|| (part_eta[jp]>0 && part_eta[kp]>0 && part_eta[ip]<0 && part_eta[lp]<0) 
									 ){
									hprof_4pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
								}

								//j-l-i-k
								dphi = part_phi[jp] + part_phi[lp] - part_phi[ip] - part_phi[kp];
								hprof_mid_4pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
								hprof_fwd_4pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
								hprof_4pc_mult_evt->Fill(0.5,cos(2*dphi));

								if ( 
										(part_eta[jp]<0 && part_eta[lp]<0 && part_eta[ip]>0 && part_eta[kp]>0)
										|| (part_eta[jp]>0 && part_eta[lp]>0 && part_eta[ip]<0 && part_eta[kp]<0) 
									 ){
									hprof_4pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
								}

								//k-l-i-j
								dphi = part_phi[kp] + part_phi[lp] - part_phi[ip] - part_phi[jp];
								hprof_mid_4pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
								hprof_fwd_4pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
								hprof_4pc_mult_evt->Fill(0.5,cos(2*dphi));

								if ( 
										(part_eta[kp]<0 && part_eta[lp]<0 && part_eta[ip]>0 && part_eta[jp]>0)
										|| (part_eta[kp]>0 && part_eta[lp]>0 && part_eta[ip]<0 && part_eta[jp]<0) 
									 ){
									hprof_4pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
								}
							}

						}//lp
					}//kp
				}//jp
			}//ip

			//for (int ip=0; ip<npart; ip++){
			for (int ip=0; ip<0; ip++){ //skip for now
				if ( fabs(part_eta[ip])>2.5 ) continue;
				if ( part_pt[ip]<ptcut[0] || part_pt[ip]>ptcut[npt] ) continue;

				for (int jp=ip+1; jp<npart; jp++){
					if ( fabs(part_eta[jp])>2.5 ) continue;
					if ( part_pt[jp]<ptcut[0] || part_pt[jp]>ptcut[npt] ) continue;

					if ( part_pt[ip]<pt_max && part_pt[jp]<pt_max ){
						float dphi = part_phi[ip] - part_phi[jp];
						hprof_mid_2pc_ref[imult]->Fill(0.5,cos(2*dphi)); 
						hprof_fwd_2pc_ref[imult_fwd]->Fill(0.5,cos(2*dphi)); 
						hprof_2pc_mult_evt->Fill(0.5,cos(2*dphi));

						if ( part_eta[ip]*part_eta[jp]<0 ){
							hprof_2pc_mult_sub_evt->Fill(0.5,cos(2*dphi));
						}
					}//pt_max

				}//jp
			}//ip

			hprof_4pc_mult_dir->Fill(nmult_mid, hprof_4pc_mult_evt->GetBinContent(1));
			hprof_2pc_mult_dir->Fill(nmult_mid, hprof_2pc_mult_evt->GetBinContent(1));

			hprof_4pc_mult_sub_dir->Fill(nmult_mid, hprof_4pc_mult_sub_evt->GetBinContent(1));
			hprof_2pc_mult_sub_dir->Fill(nmult_mid, hprof_2pc_mult_sub_evt->GetBinContent(1));
		}//ien


		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	hprof_4pc_mult_sub_dir->Write();
	hprof_2pc_mult_sub_dir->Write();
	hprof_4pc_mult_dir->Write();
	hprof_2pc_mult_dir->Write();

	for (int icut=0; icut<nptmax; icut++){
		hprof_4pc_mult_sub[icut]->Write();
		hprof_2pc_mult_sub[icut]->Write();
		hprof_4pc_mult[icut]->Write();
		hprof_2pc_mult[icut]->Write();
		hprof_4pc_mult_3sub[icut]->Write();
		hprof_2pc_mult_3sub_ab[icut]->Write();
		hprof_2pc_mult_3sub_ac[icut]->Write();
	}//icut

	heta_pt_all->Write();
	heta_pt_ana->Write();
	hmult_all->Write();
	hmult_mid->Write();
	hmult_fwd->Write();
	hmult_cnt_fwd->Write();
	hmult_mid_fwd->Write();
	hntrig_mid->Write();
	hntrig_fwd->Write();
	outfile->Close();

}

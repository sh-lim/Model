#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>
#include <TLorentzVector.h>


void JtCalculate01(){

	double j_edge[28]= {0,0.01,0.0524807,0.0630957,0.0758577,0.0912006,0.109648,0.131825,0.15849,0.190546,0.229087,0.275423,0.331131,0.398108,0.47863,0.575439,0.691831,0.831764,0.999998,1.20226,1.44544,1.7378,2.08929,2.51189,3.01996,3.63078,4.36516,100};
	double z_edge[28] = {0};
	for (int ii=1; ii<28; ii++){
		z_edge[ii] = pow(10,-2.7+0.1*ii);
	}

	const int npt = 3;
	const float ptbin[npt+1] = {25, 40, 60, 80};

	const float const_pi = TMath::Pi();
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	int i_p_id[2000];
	bool b_p_chg[2000];

	float f_jet_pt[20], f_jet_eta[20], f_jet_phi[20];

	int i_np, i_njet;
	float dR = 0.4;
	float area = const_pi*dR*dR;

	TH1D *njet = new TH1D("njet","",10,0,10);
	TH1D *jetpt = new TH1D("jetpt", "jetpt", 100, 0, 200);

	TH1D *jeteta[npt];
	TH1D *perpeta[npt];

	TH1D *incljt[npt];
	TH1D *inclz[npt];
	TH1D *perp_incljt[npt];
	TH1D *perp_inclz[npt];

	for (int ii=0; ii<npt; ii++){
		perpeta[ii] = new TH1D(Form("perpeta_pt%d",ii), Form("perpeta_pt%d",ii), 50, -2.5, 2.5);
		jeteta[ii] = new TH1D(Form("jeteta_pt%d",ii), Form("jeteta_pt%d",ii), 50, -2.5, 2.5);

		incljt[ii] = new TH1D(Form("incljt_pt%d",ii), Form("incljt_pt%d",ii), 27, j_edge);
		inclz[ii] = new TH1D(Form("inclz_pt%d",ii), Form("inclz_pt%d",ii), 27,z_edge);
		perp_incljt[ii] = new TH1D(Form("perp_incljt_pt%d",ii), Form("perp_incljt_pt%d",ii), 27, j_edge);
		perp_inclz[ii] = new TH1D(Form("perp_inclz_pt%d",ii), Form("perp_inclz_pt%d",ii), 27,z_edge);
	}

	ifstream flist;
	flist.open("file.lst");

	string fname;

	while ( flist >> fname ){
		TFile *f = new TFile(fname.c_str(), "read");

		if ( f->IsOpen() ){
			cout << "OPEN: " << f->GetName() << endl;
		}else{
			continue;
		}
		TTree *T = (TTree*)f ->Get("T");

		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_pt",f_p_pt);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_chg",b_p_chg);
		T->SetBranchAddress("p_id",&i_p_id);

		T->SetBranchAddress("njetR04",&i_njet);
		T->SetBranchAddress("jetR04_pt",f_jet_pt);
		T->SetBranchAddress("jetR04_eta",f_jet_eta);
		T->SetBranchAddress("jetR04_phi",f_jet_phi);

		int nent = T->GetEntries();
		cout << "nentries: " << nent << endl;

		for (int ien=0; ien<nent; ien++)
		{
			T->GetEntry(ien);

			//Multiplicity calculation
			int mult = 0;

			//charged particle
			for (int ip=0; ip<i_np; ip++){
				if ( !b_p_chg[ip] ) continue;
				if ( fabs(f_p_eta[ip])>1 ) continue;

				mult++;
			}
			if ( mult<60 || mult>=80 ) continue;

			//event selection
			njet->Fill(i_njet);
			if ( i_njet==0 ) continue;

			//leading jet and rorate
			TVector3 vec_leadjet;
			vec_leadjet.SetPtEtaPhi(f_jet_pt[0], f_jet_eta[0], f_jet_phi[0]);
			TVector3 vOrtho(vec_leadjet);
			TVector3 perpjet;
			vOrtho.RotateZ(const_pi/2.);

			float bkg_pt = 0.0;
			float bkg_ptot = 0.0;

			for (int itr = 0; itr < i_np; itr++)
			{
				if (fabs(f_p_eta[itr])>1) continue;
				if ( abs(i_p_id[itr])==12 || abs(i_p_id[itr])==14 || abs(i_p_id[itr])==16 ) continue;
				float deta = f_p_eta[itr] - vOrtho.Eta();
				float dphi = f_p_phi[itr] - vOrtho.Phi();

				if ( dphi < -const_pi ) dphi += 2*const_pi;
				else if ( dphi > const_pi ) dphi -= 2*const_pi;

				float newdR = sqrt(deta*deta + dphi*dphi);
				if(newdR>0.4) continue; 

				TVector3 vec;
				vec.SetPtEtaPhi(f_p_pt[itr], f_p_eta[itr], f_p_phi[itr]);

				bkg_pt += f_p_pt[itr];
				bkg_ptot += vec.Mag(); 
			}//itr  track loop

			//signal jets
			for (int kjet = 0; kjet < i_njet; kjet++)
			{
				if (fabs(f_jet_eta[kjet])>0.5) continue;
				//bakcground subtraction
				TVector3 vec_jet;
				vec_jet.SetPtEtaPhi(f_jet_pt[kjet], f_jet_eta[kjet], f_jet_phi[kjet]);
				vec_jet.SetMag(vec_jet.Mag() - bkg_ptot);

				if (vec_jet.Pt()<25 || vec_jet.Pt()>80) continue;

				int ptbin = 0;
				if ( vec_jet.Pt()>60 ) ptbin = 2; 
				else if ( vec_jet.Pt()>40 ) ptbin = 1; 

				jetpt->Fill(vec_jet.Pt());
				jeteta[ptbin]->Fill(vec_jet.Eta());

				TVector3 vPerp(vec_jet);
				vPerp.RotateZ(const_pi/2.);

				for (int ist = 0; ist < i_np; ist++)
				{
					if (fabs(f_p_eta[ist])>1) continue;
					if (!b_p_chg[ist]) continue;
					float deta = f_p_eta[ist] - vec_jet.Eta();
					float dphi = f_p_phi[ist] - vec_jet.Phi();

					if ( dphi < -const_pi ) dphi += 2*const_pi;
					else if ( dphi > const_pi ) dphi -= 2*const_pi;

					float newdR = sqrt(deta*deta + dphi*dphi);
					if(newdR>0.4) continue;

					TVector3 constptc;
					constptc.SetPtEtaPhi(f_p_pt[ist], f_p_eta[ist], f_p_phi[ist]);
					TVector3 vec_cross = vec_jet.Cross(constptc);

					float f_jt = vec_cross.Mag()/vec_jet.Mag();
					incljt[ptbin]->Fill(f_jt);

					float f_dot = vec_jet.Dot(constptc);
					float f_zh= f_dot/vec_jet.Mag2();

					inclz[ptbin]->Fill(f_zh);

				}//ist signal track loop

				//scan all inclusive jet
				bool bBAD = false;
				for (int ijet = 0; ijet < i_njet; ijet++)
				{
					float djeteta = f_jet_eta[ijet] - vPerp.Eta();
					float djetphi = f_jet_phi[ijet] - vPerp.Phi();

					if ( djetphi < -const_pi ) djetphi += 2*const_pi;
					else if ( djetphi > const_pi ) djetphi -= 2*const_pi;
					float dR = sqrt(djeteta*djeteta+djetphi*djetphi);

					if(dR<0.8){
						bBAD = true;
						break;
					}
				}

				if ( bBAD ) continue;

				perpeta[ptbin]->Fill(vPerp.Eta());

				for (int itr = 0; itr < i_np; itr++)
				{
					if (fabs(f_p_eta[itr])>1) continue;
					if (!b_p_chg[itr]) continue;
					float deta = f_p_eta[itr] - vPerp.Eta();
					float dphi = f_p_phi[itr] - vPerp.Phi();

					if ( dphi < -const_pi ) dphi += 2*const_pi;
					else if ( dphi > const_pi ) dphi -= 2*const_pi;

					float newdR = sqrt(deta*deta + dphi*dphi);

					if(newdR>0.4) continue;

					TVector3 constptc;
					constptc.SetPtEtaPhi(f_p_pt[itr], f_p_eta[itr], f_p_phi[itr]);
					TVector3 vec_cross = vPerp.Cross(constptc);

					float f_jt = vec_cross.Mag()/vPerp.Mag();
					perp_incljt[ptbin]->Fill(f_jt);

					float f_dot = vPerp.Dot(constptc);
					float f_zh= f_dot/vPerp.Mag2();

					perp_inclz[ptbin]->Fill(f_zh);

				}//itr signal perpjet loop

			}//jet

		}//ien

		f->Close();
		delete f;

	}

	flist.close();

	TFile *outfile = new TFile("outfile_hist.root", "RECREATE");


	njet->Write();
	jetpt->Write();

	for (int ii=0; ii<npt; ii++){
		jeteta[ii]->Write();
		perpeta[ii]->Write();

		incljt[ii]->Write();
		inclz[ii]->Write();

		perp_incljt[ii]->Write();
		perp_inclz[ii]->Write();
	}

	outfile->Close();

}//void

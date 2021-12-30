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


void JtCalculate(){

	double j_edge[28]= {0,0.01,0.0524807,0.0630957,0.0758577,0.0912006,0.109648,0.131825,0.15849,0.190546,0.229087,0.275423,0.331131,0.398108,0.47863,0.575439,0.691831,0.831764,0.999998,1.20226,1.44544,1.7378,2.08929,2.51189,3.01996,3.63078,4.36516,100};
	double z_edge[28] = {0};
	for (int ii=1; ii<28; ii++){
		z_edge[ii] = pow(10,-2.7+0.1*ii);
	}

	//return;

	const float const_pi = TMath::Pi();
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	float B_f_p_pt[2000], B_f_p_eta[2000], B_f_p_phi[2000];
	int i_p_id[2000], B_i_p_id[2000];
	bool b_p_chg[2000], B_b_p_chg[2000];

	float f_jet_pt[20], f_jet_eta[20], f_jet_phi[20];
	float I_f_jet_pt[20], I_f_jet_eta[20], I_f_jet_phi[20];

	int i_np, i_njet;
	int B_i_np, I_i_njet;
	float dR = 0.4;
	float area = const_pi*dR*dR;

	TH1D *njet = new TH1D("njet","",10,0,10);

	TH1D *jetpt = new TH1D("jetpt", "jetpt", 100, 0, 200);

	TH1D *sigjt = new TH1D("sigjt", "sigjt", 27, j_edge);
	TH1D *sigz = new TH1D("sigz", "sigz", 27,z_edge);
	TH1D *incljt = new TH1D("incljt", "incljt", 27, j_edge);
	TH1D *inclz = new TH1D("inclz", "inclz", 27,z_edge);

	TH1D *perp_sigjt = new TH1D("perp_sigjt", "perp_sigjt", 27, j_edge);
	TH1D *perp_sigz = new TH1D("perp_sigz", "perp_sigz", 27,z_edge);
	TH1D *perp_incljt = new TH1D("perp_incljt", "perp_incljt", 27, j_edge);
	TH1D *perp_inclz = new TH1D("perp_inclz", "perp_inclz", 27,z_edge);

	TH1D *perpeta = new TH1D("perpeta", "perpeta", 50, -2.5, 2.5);
	TH1D *jeteta = new TH1D("jeteta", "jeteta", 50, -2.5, 2.5);


	for (int ii=0; ii<50; ii++){
		TFile *f = new TFile(Form("/alice/data/shlim/PYTHIA/pp5TeV_set22_grp000/outfile_pp5TeV_set22_grp000_%05d.root",ii), "read");

		if ( f->IsOpen() ){
			cout << "OPEN: " << f->GetName() << endl;
		}else{
			continue;
		}
		TTree *T = (TTree*)f ->Get("T");

		T->SetBranchAddress("Sig_p_pt",f_p_pt);
		T->SetBranchAddress("Sig_np",&i_np);
		T->SetBranchAddress("Sig_p_eta",f_p_eta);
		T->SetBranchAddress("Sig_p_phi",f_p_phi);
		T->SetBranchAddress("Sig_p_chg",b_p_chg);
		T->SetBranchAddress("Sig_p_id",&i_p_id);

		T->SetBranchAddress("Bkg_p_pt",B_f_p_pt);
		T->SetBranchAddress("Bkg_np",&B_i_np);
		T->SetBranchAddress("Bkg_p_eta",B_f_p_eta);
		T->SetBranchAddress("Bkg_p_phi",B_f_p_phi);
		T->SetBranchAddress("Bkg_p_chg",B_b_p_chg);
		T->SetBranchAddress("Bkg_p_id",&B_i_p_id);

		T->SetBranchAddress("nSigjetR04",&i_njet);
		T->SetBranchAddress("SigjetR04_pt",f_jet_pt);
		T->SetBranchAddress("SigjetR04_eta",f_jet_eta);
		T->SetBranchAddress("SigjetR04_phi",f_jet_phi);

		T->SetBranchAddress("nIncjetR04",&I_i_njet);
		T->SetBranchAddress("IncjetR04_pt",I_f_jet_pt);
		T->SetBranchAddress("IncjetR04_eta",I_f_jet_eta);
		T->SetBranchAddress("IncjetR04_phi",I_f_jet_phi);


		int nent = T->GetEntries();
		cout << "nentries: " << nent << endl;

		for (int ien=0; ien<nent; ien++)
		{
			T->GetEntry(ien);

			//event selection
			if ( I_i_njet!=i_njet ) continue;
			njet->Fill(i_njet);

			//if (i_njet>2) continue;

			//signal jets
			TVector3 vec_jet;
			for (int kjet = 0; kjet < i_njet; kjet++)
			{
				if (fabs(f_jet_eta[kjet])>0.5) continue;
				if (f_jet_pt[kjet]<25 || f_jet_pt[kjet]>40) continue;

				jetpt->Fill(f_jet_pt[kjet]);
				jeteta->Fill(f_jet_eta[kjet]);

				vec_jet.SetPtEtaPhi(f_jet_pt[kjet], f_jet_eta[kjet], f_jet_phi[kjet]);
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
					sigjt->Fill(f_jt);
					incljt->Fill(f_jt);

					float f_dot = vec_jet.Dot(constptc);
					float f_zh= f_dot/vec_jet.Mag2();

					sigz->Fill(f_zh);
					inclz->Fill(f_zh);

				}//ist signal track loop

				for (int ibt = 0; ibt < i_np; ibt++)
				{
					if (fabs(B_f_p_eta[ibt])>1) continue;
					if (!B_b_p_chg[ibt]) continue;
					float deta = B_f_p_eta[ibt] - vec_jet.Eta();
					float dphi = B_f_p_phi[ibt] - vec_jet.Phi();

					if ( dphi < -const_pi ) dphi += 2*const_pi;
					else if ( dphi > const_pi ) dphi -= 2*const_pi;

					float newdR = sqrt(deta*deta + dphi*dphi);
					if(newdR>0.4) continue;

					TVector3 constptc;
					constptc.SetPtEtaPhi(B_f_p_pt[ibt], B_f_p_eta[ibt], B_f_p_phi[ibt]);
					TVector3 vec_cross = vec_jet.Cross(constptc);

					float f_jt = vec_cross.Mag()/vec_jet.Mag();
					incljt->Fill(f_jt);

					float f_dot = vec_jet.Dot(constptc);
					float f_zh= f_dot/vec_jet.Mag2();

					inclz->Fill(f_zh);

				}//ist signal track loop


				//scan all inclusive jet
				bool bBAD = false;
				for (int ijet = 0; ijet < I_i_njet; ijet++)
				{
					float djeteta = I_f_jet_eta[ijet] - vPerp.Eta();
					float djetphi = I_f_jet_phi[ijet] - vPerp.Phi();

					if ( djetphi < -const_pi ) djetphi += 2*const_pi;
					else if ( djetphi > const_pi ) djetphi -= 2*const_pi;
					float dR = sqrt(djeteta*djeteta+djetphi*djetphi);

					if(dR<0.8){
						bBAD = true;
						break;
					}
				}

				if ( bBAD ) continue;

				perpeta->Fill(vPerp.Eta());

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
					perp_sigjt->Fill(f_jt);
					perp_incljt->Fill(f_jt);

					float f_dot = vPerp.Dot(constptc);
					float f_zh= f_dot/vPerp.Mag2();

					perp_sigz->Fill(f_zh);
					perp_inclz->Fill(f_zh);

				}//itr signal perpjet loop


				for (int itb = 0; itb < B_i_np; itb++)
				{
					if (fabs(B_f_p_eta[itb])>1) continue;
					if (!B_b_p_chg[itb]) continue;
					float deta = B_f_p_eta[itb] - vPerp.Eta();
					float dphi = B_f_p_phi[itb] - vPerp.Phi();

					if ( dphi < -const_pi ) dphi += 2*const_pi;
					else if ( dphi > const_pi ) dphi -= 2*const_pi;

					float newdR = sqrt(deta*deta + dphi*dphi);

					if(newdR>0.4) continue;

					TVector3 constptc;
					constptc.SetPtEtaPhi(B_f_p_pt[itb], B_f_p_eta[itb], B_f_p_phi[itb]);
					TVector3 vec_cross = vPerp.Cross(constptc);

					float f_jt = vec_cross.Mag()/vPerp.Mag();
					perp_incljt->Fill(f_jt);

					float f_dot = vPerp.Dot(constptc);
					float f_zh= f_dot/vPerp.Mag2();

					perp_inclz->Fill(f_zh);

				}//itb background perpjet loop
			}

		}//ien

		f->Close();
		delete f;

	}


	TFile *outfile = new TFile("jtresult.root", "RECREATE");


	njet->Write();

	jetpt->Write();
	jeteta->Write();
	perpeta->Write();

	incljt->Write();
	inclz->Write();
	sigjt->Write();
	sigz->Write();

	perp_incljt->Write();
	perp_inclz->Write();
	perp_sigjt->Write();
	perp_sigz->Write();

	outfile->Close();

}//void

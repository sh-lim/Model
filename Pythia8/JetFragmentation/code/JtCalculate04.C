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


void JtCalculate04(const bool isMB = true){

	//const int njt = 26;
	//double j_edge[njt+1]= {0.01,0.0524807,0.0630957,0.0758577,0.0912006,0.109648,0.131825,0.15849,0.190546,0.229087,0.275423,0.331131,0.398108,0.47863,0.575439,0.691831,0.831764,0.999998,1.20226,1.44544,1.7378,2.08929,2.51189,3.01996,3.63078,4.36516,10};

	const int njt = 64;
	double j_edge[njt+1] = {
		0.01, 0.0111397, 0.0124094, 0.0138237, 0.0153993, 0.0171544, 0.0191095, 0.0212875, 0.0237137, 0.0264165, 0.0294273, 0.0327812, 0.0365174, 0.0406794, 0.0453158, 0.0504807, 0.0562341, 0.0626434, 0.0697831, 0.0777365, 0.0865964, 0.0964662, 0.107461, 0.119709, 0.133352, 0.148551, 0.165482, 0.184342, 0.205353, 0.228757, 0.25483, 0.283874, 0.316228, 0.352269, 0.392419, 0.437144, 0.486968, 0.542469, 0.604296, 0.67317, 0.749894, 0.835363, 0.930572, 1.03663, 1.15478, 1.2864, 1.43301, 1.59634, 1.77828, 1.98096, 2.20673, 2.45824, 2.73842, 3.05053, 3.39821, 3.78552, 4.21697, 4.69759, 5.23299, 5.82942, 6.49382, 7.23394, 8.05842, 8.97687, 10
	};

	/*
	double z_edge[28] = {0};
	for (int ii=1; ii<28; ii++){
		z_edge[ii] = pow(10,-2.7+0.1*ii);
	}
	*/

	const int nmult = 1;
	const int npt = 4;
	//const float pt_edge[npt+1] = {10, 20, 40, 60, 100};
	const float pt_edge[npt+1] = {20, 40, 60, 80, 100};

	//const float eta_edge = 0.50;
	const float eta_edge = 0.25;

	const int nz = 3;

	TH1D *jetptbin = new TH1D("jetptbin","",npt,pt_edge);

	const float const_pi = TMath::Pi();
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	int i_p_id[2000];
	bool b_p_chg[2000];

	float f_jet_pt[20], f_jet_eta[20], f_jet_phi[20];

	int i_np, i_njet;
	float f_scale;
	float dR = 0.4;
	float area = const_pi*dR*dR;

	TH1D *njet[nmult];
	TH1D *jetpt[nmult];
	for (int im=0; im<nmult; im++){
		njet[im] = new TH1D(Form("njet_m%d",im),"",10,0,10);
		jetpt[im] = new TH1D(Form("jetpt_m%d",im), "", 100, 0, 200);
	}

	TH1D *scale[nmult][npt];

	TH1D *jeteta[nmult][npt];
	TH1D *perpeta[nmult][npt];

	TH1D *incljt[nmult][npt];
	TH1D *perp_incljt[nmult][npt];

	TH1D *zdepjt[nmult][npt][nz];
	TH1D *perp_zdepjt[nmult][npt][nz];

	for (int im=0; im<nmult; im++){
		for (int ii=0; ii<npt; ii++){
			scale[im][ii] = new TH1D(Form("scale_m%d_pt%d",im,ii), "", 100,0,100);

			perpeta[im][ii] = new TH1D(Form("perpeta_m%d_pt%d",im,ii), "", 50, -2.5, 2.5);
			jeteta[im][ii] = new TH1D(Form("jeteta_m%d_pt%d",im,ii), "", 50, -2.5, 2.5);

			incljt[im][ii] = new TH1D(Form("incljt_m%d_pt%d",im,ii), "", njt, j_edge);
			perp_incljt[im][ii] = new TH1D(Form("perp_incljt_m%d_pt%d",im,ii), "", njt, j_edge);

			for (int iz=0; iz<nz; iz++){
				zdepjt[im][ii][iz] = new TH1D(Form("zdepjt_m%d_pt%d_z%d",im,ii,iz), "", njt, j_edge);
				perp_zdepjt[im][ii][iz] = new TH1D(Form("perp_zdepjt_m%d_pt%d_z%d",im,ii,iz), "", njt, j_edge);
			}
		}
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

		if ( !isMB ){
			T->SetBranchAddress("scale",&f_scale);
		}

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

		//T->SetBranchAddress("nchgjetR04",&i_njet);
		//T->SetBranchAddress("chgjetR04_pt",f_jet_pt);
		//T->SetBranchAddress("chgjetR04_eta",f_jet_eta);
		//T->SetBranchAddress("chgjetR04_phi",f_jet_phi);

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

			//event selection
			if ( isMB ){
				njet[0]->Fill(i_njet);
			}
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
				if ( fabs(f_p_eta[itr])>0.9 ) continue;
				if ( !b_p_chg[itr] ) continue;
				//if ( abs(i_p_id[itr])==12 || abs(i_p_id[itr])==14 || abs(i_p_id[itr])==16 ) continue;
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
				if ( fabs(f_jet_eta[kjet])>eta_edge ) continue;
				//bakcground subtraction
				TVector3 vec_jet;
				vec_jet.SetPtEtaPhi(f_jet_pt[kjet], f_jet_eta[kjet], f_jet_phi[kjet]);
				//vec_jet.SetMag(vec_jet.Mag() - bkg_ptot);

				if (vec_jet.Pt()<pt_edge[0] || vec_jet.Pt()>pt_edge[npt]) continue;
				int ptbin = jetptbin->FindBin(vec_jet.Pt()) - 1;

				jetpt[0]->Fill(vec_jet.Pt());
				jeteta[0][ptbin]->Fill(vec_jet.Eta());

				TVector3 vPerp(vec_jet);
				vPerp.RotateZ(const_pi/2.);

				//continue;

				for (int ist = 0; ist < i_np; ist++)
				{
					if ( fabs( f_p_eta[ist])>0.9 ) continue;
					if ( f_p_pt[ist]<0.15 ) continue;
					if ( !b_p_chg[ist] ) continue;
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
					float f_dot = vec_jet.Dot(constptc);
					float f_zh = f_dot/vec_jet.Mag2();

					incljt[0][ptbin]->Fill(f_jt);

					if ( f_zh<0.2 ){
						zdepjt[0][ptbin][0]->Fill(f_jt);
					}else if ( f_zh<0.4 ){
						zdepjt[0][ptbin][1]->Fill(f_jt);
					}else{
						zdepjt[0][ptbin][2]->Fill(f_jt);
					}

				}//ist signal track loop

				//continue;

				//scan all inclusive jet
				bool bBAD = false;
				for (int ijet = 0; ijet < i_njet; ijet++)
				{
					if ( f_jet_pt[ijet]<10.0 ) continue;
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

				perpeta[0][ptbin]->Fill(vPerp.Eta());

				//continue;

				for (int itr = 0; itr < i_np; itr++)
				{
					if ( fabs(f_p_eta[itr])>0.9 ) continue;
					if ( f_p_pt[itr]<0.15 ) continue;
					if ( !b_p_chg[itr] ) continue;
					float deta = f_p_eta[itr] - vPerp.Eta();
					float dphi = f_p_phi[itr] - vPerp.Phi();

					if ( dphi < -const_pi ) dphi += 2*const_pi;
					else if ( dphi > const_pi ) dphi -= 2*const_pi;

					float newdR = sqrt(deta*deta + dphi*dphi);

					if ( newdR>0.4 ) continue;

					TVector3 constptc;
					constptc.SetPtEtaPhi(f_p_pt[itr], f_p_eta[itr], f_p_phi[itr]);
					TVector3 vec_cross = vPerp.Cross(constptc);

					float f_jt = vec_cross.Mag()/vPerp.Mag();
					float f_dot = vPerp.Dot(constptc);
					float f_zh = f_dot/vPerp.Mag2();

					perp_incljt[0][ptbin]->Fill(f_jt);

					if ( f_zh<0.2 ){
						perp_zdepjt[0][ptbin][0]->Fill(f_jt);
					}else if ( f_zh<0.4 ){
						perp_zdepjt[0][ptbin][1]->Fill(f_jt);
					}else{
						perp_zdepjt[0][ptbin][2]->Fill(f_jt);
					}

				}//itr signal perpjet loop

			}//jet

		}//ien

		f->Close();
		delete f;

	}

	flist.close();

	TFile *outfile = new TFile("outfile_JtCalculate04_hist.root", "RECREATE");


	for (int im=0; im<nmult; im++){
		njet[im]->Write();
		jetpt[im]->Write();
	}

	for (int im=0; im<nmult; im++){
		for (int ii=0; ii<npt; ii++){
			scale[im][ii]->Write();

			jeteta[im][ii]->Write();
			perpeta[im][ii]->Write();

			incljt[im][ii]->Write();
			perp_incljt[im][ii]->Write();
			for (int iz=0; iz<nz; iz++){
				zdepjt[im][ii][iz]->Write();
				perp_zdepjt[im][ii][iz]->Write();
			}
		}
	}

	outfile->Close();

}//void

void ScanVnV2(const char *dataset="pPb5TeV_grp5"){

	char fname[500];

	ifstream flist;
	//sprintf(fname, "file_%s.lst", dataset);
	sprintf(fname, "file.lst");
	flist.open(fname);

	int npart;
	int part_pid[5000];
	float part_pt[5000], part_eta[5000], part_phi[5000];

	int npp;
	int pp_pid[5000];
	float pp_x[5000], pp_y[5000], pp_z[5000];
	float pp_px[5000], pp_py[5000], pp_pz[5000];

	int nph;
	int ph_ppid[5000], ph_hpid[5000];
	float ph_px[5000], ph_py[5000], ph_pz[5000];

	TH1D *hnPsi[3][2];
	TH1D *hcosndPsi[3][2];
	TProfile *hcosndphi[3][2][4]; //[EP][order][pid]

	TH2D *hPsiCorr[3][2];

	TH1D *hnPPsi[2];
	TProfile *hcosndphiP[2][4]; //[order][pid]

	TProfile *hcosndphiInitParton[2][2]; //[order][pid]
	TProfile *hcosndphiFinalParton[2][2]; //[order][pid]
	TProfile *hcosndphiFinalPartonH[2][2]; //[order][pid]

	TH2D *h2d_eta_pt_init_parton = new TH2D("h2d_eta_pt_init_parton","",100,-5,5,100,0,10);
	TH2D *h2d_eta_pt_final_parton = new TH2D("h2d_eta_pt_final_parton","",100,-5,5,100,0,10);

	TH1D *hmult_bwd = new TH1D("hmult_bwd","",500,0,500);

	const float ptbinv2[15] = {
		0.00, 0.20, 0.40, 0.60, 0.80,
		1.00, 1.25, 1.50, 2.00, 3.00,
		4.00, 5.00, 6.00, 8.00, 10.0
	};

	const float ptbinv3[6] = {
		0.00, 0.50, 1.00, 1.50, 2.00, 3.00
	};

	//ii-event plane
	//jj-order(v2, v3)
	for (int ii=0; ii<3; ii++){
		for (int jj=0; jj<2; jj++){
			hnPsi[ii][jj] = new TH1D(Form("hnPsi_eta%d_order%d",ii,jj+2),"",36,-TMath::Pi(),TMath::Pi());
			hcosndPsi[ii][jj] = new TH1D(Form("hcosndPsi_eta%d_order%d",ii,jj+2),"",36,-1.0,1.0);

			for (int kk=0; kk<4; kk++){
				if ( jj==0 ){
					hcosndphi[ii][jj][kk] = new TProfile(Form("hcosndphi_eta%d_order%d_pid%d",ii,jj+2,kk),"",14, ptbinv2);
				}else{
					hcosndphi[ii][jj][kk] = new TProfile(Form("hcosndphi_eta%d_order%d_pid%d",ii,jj+2,kk),"",5, ptbinv3);
				}
			}

			hPsiCorr[ii][jj] = new TH2D(Form("hPsiCorr_eta%d_order%d",ii,jj+2),"",72,-TMath::Pi(),TMath::Pi(),72,-TMath::Pi(),TMath::Pi());
		}
	}

	for (int jj=0; jj<2; jj++){
		hnPPsi[jj] = new TH1D(Form("hnPPsi_order%d",jj+2),"",36,-TMath::Pi(),TMath::Pi());
		for (int kk=0; kk<4; kk++){
			if ( jj==0 ){
				hcosndphiP[jj][kk] = new TProfile(Form("hcosndphiP_order%d_pid%d",jj+2,kk),"",14, ptbinv2);
			}else{
				hcosndphiP[jj][kk] = new TProfile(Form("hcosndphiP_order%d_pid%d",jj+2,kk),"",5, ptbinv3);
			}
		}
	}

	for (int jj=0; jj<2; jj++){
		for (int kk=0; kk<2; kk++){
			if ( jj==0 ){
				hcosndphiInitParton[jj][kk] = new TProfile(Form("hcosndphiInitParton_order%d_pid%d",jj+2,kk),"",14, ptbinv2);
				hcosndphiFinalParton[jj][kk] = new TProfile(Form("hcosndphiFinalParton_order%d_pid%d",jj+2,kk),"",14, ptbinv2);
				hcosndphiFinalPartonH[jj][kk] = new TProfile(Form("hcosndphiFinalPartonH_order%d_pid%d",jj+2,kk),"",14, ptbinv2);
			}else{
				hcosndphiInitParton[jj][kk] = new TProfile(Form("hcosndphiInitParton_order%d_pid%d",jj+2,kk),"",5, ptbinv3);
				hcosndphiFinalParton[jj][kk] = new TProfile(Form("hcosndphiFinalParton_order%d_pid%d",jj+2,kk),"",5, ptbinv3);
				hcosndphiFinalPartonH[jj][kk] = new TProfile(Form("hcosndphiFinalPartonH_order%d_pid%d",jj+2,kk),"",5, ptbinv3);
			}
		}
	}

	int cut_mult_bwd = 0;
	if ( strstr(dataset,"He3Au200") ){
		cut_mult_bwd = 27;
	}else if ( strstr(dataset,"dAu200") ){
		cut_mult_bwd = 21;
	}else if ( strstr(dataset,"pAu200") ){
		cut_mult_bwd = 14;
	}else if ( strstr(dataset,"pPb5TeV") ){
		//cut_mult_bwd = 71;
		//all particle in -5<eta<-3
		cut_mult_bwd = 140;
	}else{
		cut_mult_bwd = 140;
	}

	cout << "dataset: " << dataset << endl;
	cout << "multiplicity cut for 0-5%: " << cut_mult_bwd << endl;


	while ( flist >> fname ){

		cout << "OPEN: " << fname << endl;
		TFile *infile = new TFile(fname,"read");

		TTree *T = (TTree*)infile->Get("T");
		TTree *Tpp = (TTree*)infile->Get("Tpp");
		TTree *Tph = (TTree*)infile->Get("Tph");

		if ( !T || !Tpp || !Tph ){
			infile->Close();
			delete infile;
			continue;
		}

		T->SetBranchAddress("npart",&npart);
		T->SetBranchAddress("part_pid",part_pid);
		T->SetBranchAddress("part_pt",part_pt);
		T->SetBranchAddress("part_eta",part_eta);
		T->SetBranchAddress("part_phi",part_phi);

		Tpp->SetBranchAddress("npp",&npp);
		Tpp->SetBranchAddress("pp_pid",pp_pid);
		Tpp->SetBranchAddress("pp_x",pp_x);
		Tpp->SetBranchAddress("pp_y",pp_y);
		Tpp->SetBranchAddress("pp_z",pp_z);
		Tpp->SetBranchAddress("pp_px",pp_px);
		Tpp->SetBranchAddress("pp_py",pp_py);
		Tpp->SetBranchAddress("pp_pz",pp_pz);

		Tph->SetBranchAddress("nph",&nph);
		Tph->SetBranchAddress("ph_ppid",ph_ppid);
		Tph->SetBranchAddress("ph_hpid",ph_hpid);
		Tph->SetBranchAddress("ph_px",ph_px);
		Tph->SetBranchAddress("ph_py",ph_py);
		Tph->SetBranchAddress("ph_pz",ph_pz);

		int nentries = T->GetEntries();
		int nentries_pp = Tpp->GetEntries();
		int nentries_ph = Tph->GetEntries();

		if ( nentries!=nentries_pp || nentries!=nentries_ph ){
			cout << "inconsistent entries, " << nentries << " " << nentries_pp << " " << nentries_ph << endl;
			infile->Close();
			delete infile;
			continue;
		}

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);
			Tpp->GetEntry(ien);
			Tph->GetEntry(ien);

			if ( npp>=5000 || npart>=5000 || nph>=5000 ) continue;

			int mult_bwd = 0;
			for (int ip=0; ip<npart; ip++){
				float eta = part_eta[ip];
				if ( eta>-5.0 && eta<-3.0 ){
					mult_bwd++;
				}
			}//ip

			hmult_bwd->Fill(mult_bwd);

			//multiplicity selection
			if ( mult_bwd<cut_mult_bwd ) continue;

			//EP with parton
			float Px[2] = {0.};
			float Py[2] = {0.};
			float Pw[2] = {0.};
			float PPsi[2] = {0.};

			float MeanX = 0.0, MeanY = 0.0;

			for (int ipp=0; ipp<npp; ipp++){
				MeanX += pp_x[ipp];
				MeanY += pp_y[ipp];
			}
			MeanX /= npp;
			MeanY /= npp;

			for (int ip=0; ip<npp; ip++){

				float r2 = (pp_x[ip]-MeanX)*(pp_x[ip]-MeanX) + (pp_y[ip]-MeanY)*(pp_y[ip]-MeanY);
				float phi = atan2(pp_y[ip]-MeanY, pp_x[ip]-MeanX);

				Px[0] += r2*cos(2*phi);
				Py[0] += r2*sin(2*phi);

				Px[1] += r2*cos(3*phi);
				Py[1] += r2*sin(3*phi);
			}

			for (int jj=0; jj<2; jj++){
				Px[jj] /= npp;
				Py[jj] /= npp;

				float psi = atan2(Py[jj], Px[jj]) + TMath::Pi();
				PPsi[jj] = atan2(sin(psi), cos(psi))/(jj+2);

				hnPPsi[jj]->Fill((jj+2)*PPsi[jj]);
			}

			//initial parton flow
			for (int ipp=0; ipp<npp; ipp++){

				int pid = abs(pp_pid[ipp]);

				float px = pp_px[ipp];
				float py = pp_py[ipp];
				float pz = pp_pz[ipp];

				TVector3 vec(px, py, pz);

				h2d_eta_pt_init_parton->Fill(vec.Eta(), vec.Pt());

				if ( fabs(vec.Eta())>1.0 ) continue;

				for (int jj=0; jj<2; jj++){
					float dphi = (jj+2)*(vec.Phi() - PPsi[jj]);

					if ( pid==1 || pid==2 ){
						hcosndphiInitParton[jj][0]->Fill(vec.Pt(), cos(dphi));
					}else if ( pid==3 ){
						hcosndphiInitParton[jj][1]->Fill(vec.Pt(), cos(dphi));
					}

				}//jj
			}//ipp

			//final parton flow
			for (int iph=0; iph<nph; iph++){

				int ppid = abs(ph_ppid[iph]);
				int hpid = abs(ph_hpid[iph]);

				float px = ph_px[iph];
				float py = ph_py[iph];
				float pz = ph_pz[iph];

				TVector3 vec(px, py, pz);

				h2d_eta_pt_final_parton->Fill(vec.Eta(), vec.Pt());

				if ( fabs(vec.Eta())>1.0 ) continue;

				for (int jj=0; jj<2; jj++){
					float dphi = (jj+2)*(vec.Phi() - PPsi[jj]);

					if ( ppid==1 || ppid==2 ){
						hcosndphiFinalParton[jj][0]->Fill(vec.Pt(), cos(dphi));
						if ( hpid==211 ){
							hcosndphiFinalPartonH[jj][0]->Fill(vec.Pt(), cos(dphi));
						}else if ( hpid==2212 ){
							hcosndphiFinalPartonH[jj][1]->Fill(vec.Pt(), cos(dphi));
						}
					}else if ( ppid==3 ){
						hcosndphiFinalParton[jj][1]->Fill(vec.Pt(), cos(dphi));
					}

				}//jj
			}//iph

			//EP with particle
			float Qx[3][2] = {0.};
			float Qy[3][2] = {0.};
			float Qw[3] = {0.};

			for (int ip=0; ip<npart; ip++){

				float phi = part_phi[ip];
				float eta = part_eta[ip];

				if ( eta>-1.0 && eta<1.0 ){
					Qx[0][0] += cos(2*phi);
					Qy[0][0] += sin(2*phi);

					Qx[0][1] += cos(3*phi);
					Qy[0][1] += sin(3*phi);

					Qw[0]++;
				}else if ( eta>-3.0 && eta<-1.0 ){
					Qx[1][0] += cos(2*phi);
					Qy[1][0] += sin(2*phi);

					Qx[1][1] += cos(3*phi);
					Qy[1][1] += sin(3*phi);

					Qw[1]++;
				}else if ( eta>-5.0 && eta<-3.0 ){
					Qx[2][0] += cos(2*phi);
					Qy[2][0] += sin(2*phi);

					Qx[2][1] += cos(3*phi);
					Qy[2][1] += sin(3*phi);

					Qw[2]++;
				}
			}//ip

			if ( Qw[0]<1 || Qw[1]<1 || Qw[2]<1 ){
				cout << "zero multiplicity: " << Qw[0] << ", " << Qw[1] << ", " << Qw[2] << endl;
				continue;
			}

			float Psi[3][2] = {0.};
			for (int ii=0; ii<3; ii++){
				for (int jj=0; jj<2; jj++){
					Psi[ii][jj] = atan2(Qy[ii][jj], Qx[ii][jj])/(jj+2);

					hnPsi[ii][jj]->Fill((jj+2)*Psi[ii][jj]);

					hPsiCorr[ii][jj]->Fill((jj+2)*PPsi[jj], (jj+2)*Psi[ii][jj]);
				}//jj
			}//ii


			//EP resolution
			for (int jj=0; jj<2; jj++){

				//A-B
				float ndPsi = (jj+2)*(Psi[0][jj] - Psi[1][jj]);
				hcosndPsi[0][jj]->Fill(cos(ndPsi));

				//B-C
				ndPsi = (jj+2)*(Psi[1][jj] - Psi[2][jj]);
				hcosndPsi[1][jj]->Fill(cos(ndPsi));

				//C-A
				ndPsi = (jj+2)*(Psi[2][jj] - Psi[0][jj]);
				hcosndPsi[2][jj]->Fill(cos(ndPsi));

			}//


			//particle flow
			for (int ip=0; ip<npart; ip++){

				int pid = part_pid[ip];

				float pt = part_pt[ip];
				float phi = part_phi[ip];
				float eta = part_eta[ip];

				if ( fabs(eta)>1.0 ) continue;
				if ( pt>10.0 ) continue;

				//EP with particle
				for (int ii=0; ii<3; ii++){
					for (int jj=0; jj<2; jj++){
						float dphi = (jj+2)*(phi - Psi[ii][jj]);

						hcosndphi[ii][jj][0]->Fill(pt, cos(dphi));

						if ( abs(pid)==211 ){
							hcosndphi[ii][jj][1]->Fill(pt, cos(dphi));
						}else if ( abs(pid)==321 ){
							hcosndphi[ii][jj][2]->Fill(pt, cos(dphi));
						}else if ( abs(pid)==2212 ){
							hcosndphi[ii][jj][3]->Fill(pt, cos(dphi));
						}
					}//jj
				}//ii

				for (int jj=0; jj<2; jj++){
					float dphi = (jj+2)*(phi - PPsi[jj]);

					hcosndphiP[jj][0]->Fill(pt, cos(dphi));

					if ( abs(pid)==211 ){
						hcosndphiP[jj][1]->Fill(pt, cos(dphi));
					}else if ( abs(pid)==321 ){
						hcosndphiP[jj][2]->Fill(pt, cos(dphi));
					}else if ( abs(pid)==2212 ){
						hcosndphiP[jj][3]->Fill(pt, cos(dphi));
					}
				}//jj

			}//ip


		}//ien


		infile->Close();
		delete infile;

	}//

	//TFile *outfile = new TFile(Form("outfileVnV2_%s.root",dataset),"recreate");
	TFile *outfile = new TFile("outfileVnV2.root","recreate");
	
	hmult_bwd->Write();

	for (int ii=0; ii<3; ii++){
		for (int jj=0; jj<2; jj++){

			hnPsi[ii][jj]->Write();
			hcosndPsi[ii][jj]->Write();
			hPsiCorr[ii][jj]->Write();

			for (int kk=0; kk<4; kk++){
				hcosndphi[ii][jj][kk]->Write();
			}//kk
		}
	}

	for (int jj=0; jj<2; jj++){
		hnPPsi[jj]->Write();

		for (int kk=0; kk<4; kk++){
			hcosndphiP[jj][kk]->Write();
		}
	}

	for (int jj=0; jj<2; jj++){
		for (int kk=0; kk<2; kk++){
			hcosndphiInitParton[jj][kk]->Write();
			hcosndphiFinalParton[jj][kk]->Write();
			hcosndphiFinalPartonH[jj][kk]->Write();
		}
	}

	h2d_eta_pt_init_parton->Write();
	h2d_eta_pt_final_parton->Write();

}

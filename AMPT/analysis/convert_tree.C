#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <fstream>
#include <iostream>

using namespace std;

void convert_tree(const char *fname="ana/ampt.dat", const char *fname2="ana/parton-initial-afterPropagation.dat"){

	int npart = 0;
	int part_pid[2000];
	float part_eta[2000], part_phi[2000], part_pt[2000];

	TTree *T = new TTree("T","T");
	T->Branch("npart",&npart,"npart/I");
	T->Branch("part_pid",part_pid,"part_pid[npart]/I");
	T->Branch("part_eta",part_eta,"part_eta[npart]/F");
	T->Branch("part_phi",part_phi,"part_phi[npart]/F");
	T->Branch("part_pt",part_pt,"part_pt[npart]/F");

	int    evtnumber;
	int    testnum;
	int    nlist;
	double impactpar;
	int    npartproj;
	int    nparttarg;
	int    npartprojelas;
	int    npartprojinelas;
	int    nparttargelas;
	int    nparttarginelas;
	double junk;

	ifstream fdata;
	fdata.open(fname);

	while ( fdata >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk ){

		//cout << evtnumber << " " << impactpar << endl;

		npart = 0;
		for (int ii=0; ii<2000; ii++){
			part_pid[ii] = 0;
			part_eta[ii] = part_phi[ii] = part_pt[ii] = 0.0;
		}

		for (int ipart=0; ipart<nlist; ipart++){

			int    partid;
			float  pv[3];
			float  mass;
			double space[4];

			fdata >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

			TVector3 vec(pv[0], pv[1], pv[2]);

			if ( (abs(partid)==211 || abs(partid)==321 || abs(partid)==2212) && vec.Pt()>0 ){
				part_pid[npart] = partid;
				part_eta[npart] = vec.Eta();
				part_phi[npart] = vec.Phi();
				part_pt[npart] = vec.Pt();
				npart++;
			}
		}//ipart

		T->Fill();

	}//while

	fdata.close();

	int npp = 0;
	int pp_pid[2000];
	float pp_x[2000], pp_y[2000], pp_z[2000];

	TTree *Tpp = new TTree("Tpp","Tpp");
	Tpp->Branch("npp",&npp,"npp/I");
	Tpp->Branch("pp_pid",pp_pid,"pp_pid[npp]/I");
	Tpp->Branch("pp_x",pp_x,"pp_x[npp]/F");
	Tpp->Branch("pp_y",pp_y,"pp_y[npp]/F");
	Tpp->Branch("pp_z",pp_z,"pp_z[npp]/F");

	fdata.open(fname2);

	int n_baryon_formed;
	int n_meson_formed;
	int n_inipart;
	int n_inipart_notinzpc;

	while ( fdata >> evtnumber >> testnum >> nlist >> n_baryon_formed >> n_meson_formed >> n_inipart >> n_inipart_notinzpc ){

		//cout << evtnumber << " " << nlist << endl;

		npp = 0;
		for (int ii=0; ii<2000; ii++){
			pp_pid[ii] = 0;
			pp_x[ii] = pp_y[ii] = pp_z[ii] = 0.0;
		}//

		for (int ipart=0; ipart<nlist; ipart++){

			int    partid;
			float  pv[3];
			float  mass;
			double space[4];
			int tmp_a;
			float tmp_b, tmp_c;

			//fdata >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];
			fdata >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3] >> tmp_a >> tmp_b >> tmp_c;

			/*
			cout << partid << " " << pv[0] << " " << pv[1] << " " << pv[2] << " "
				<< mass << " " << space[0] << " " << space[1] << " " << space[2] << " " << space[3] << " "
				<< tmp_a << " " << tmp_b << " " << tmp_c
				<< endl;
			*/

			pp_pid[npp] = partid;
			pp_x[npp] = space[0];
			pp_y[npp] = space[1];
			pp_z[npp] = space[2];
			npp++;
		}

		Tpp->Fill();

	}//while

	fdata.close();


	TFile *outfile = new TFile("outfile.root","recreate");
	T->Write();
	Tpp->Write();
	outfile->Close();

}

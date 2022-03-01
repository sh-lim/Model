// main03.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how different processes can be selected and studied.
// All input is specified in the main03.cmnd file.
// Also illustrated output to be plotted by Python/Matplotlib/pyplot.

#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("mainEx01a.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np;
	int i_nmult_v0a;
	int i_nmult_v0c;

	int i_p_id[10000];
	int i_p_status[10000];
	float f_p_pt[10000];
	float f_p_eta[10000];
	float f_p_phi[10000];
	float f_p_vx[10000];
	float f_p_vy[10000];
	float f_p_vz[10000];
	float f_p_vt[10000];

	//bool b_p_final[10000]; 
	int i_p_mom1_id[10000];
	int i_p_mom2_id[10000];
	int i_p_mom1_status[10000];
	int i_p_mom2_status[10000];
	int i_p_mom1_index[10000];
	int i_p_mom2_index[10000];
	float f_p_mom1_vx[10000];
	float f_p_mom1_vy[10000];
	float f_p_mom1_vz[10000];
	float f_p_mom1_vt[10000];
	float f_p_mom2_vx[10000];
	float f_p_mom2_vy[10000];
	float f_p_mom2_vz[10000];
	float f_p_mom2_vt[10000];
	float f_p_mom1_mass[10000];
	float f_p_mom2_mass[10000];

	/*
	float f_jet_pt[100];
	float f_jet_eta[100];
	float f_jet_phi[100];
	*/

	// Tree output
	auto T = new TTree("T","Pythia event");
	//T->Branch("scale",&f_scale,"scale/F");
	//T->Branch("bMPI",&f_bMPI,"bMPI/F");
	//T->Branch("nMPI",&i_nMPI,"nMPI/I");
	T->Branch("nmult_v0a",&i_nmult_v0a,"nmult_v0a/I");
	T->Branch("nmult_v0c",&i_nmult_v0c,"nmult_v0c/I");

	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_status",i_p_status,"p_status[np]/I");
	//T->Branch("p_final",b_p_final,"p_final[np]/O");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_vx",f_p_vx,"p_vx[np]/F");
	T->Branch("p_vy",f_p_vy,"p_vy[np]/F");
	T->Branch("p_vz",f_p_vz,"p_vz[np]/F");
	T->Branch("p_vt",f_p_vt,"p_vt[np]/F");

	T->Branch("p_mom1_id",i_p_mom1_id,"p_mom1_id[np]/I");
	T->Branch("p_mom2_id",i_p_mom2_id,"p_mom2_id[np]/I");
	T->Branch("p_mom1_index",i_p_mom1_index,"p_mom1_index[np]/I");
	T->Branch("p_mom2_index",i_p_mom2_index,"p_mom2_index[np]/I");
	T->Branch("p_mom1_status",i_p_mom1_status,"p_mom1_status[np]/I");
	T->Branch("p_mom2_status",i_p_mom2_status,"p_mom2_status[np]/I");
	T->Branch("p_mom1_mass",f_p_mom1_mass,"p_mom1_mass[np]/F");
	T->Branch("p_mom2_mass",f_p_mom2_mass,"p_mom2_mass[np]/F");
	T->Branch("p_mom1_vx",f_p_mom1_vx,"p_mom1_vx[np]/F");
	T->Branch("p_mom1_vy",f_p_mom1_vy,"p_mom1_vy[np]/F");
	T->Branch("p_mom1_vz",f_p_mom1_vz,"p_mom1_vz[np]/F");
	T->Branch("p_mom1_vt",f_p_mom1_vt,"p_mom1_vt[np]/F");
	T->Branch("p_mom2_vx",f_p_mom2_vx,"p_mom2_vx[np]/F");
	T->Branch("p_mom2_vy",f_p_mom2_vy,"p_mom2_vy[np]/F");
	T->Branch("p_mom2_vz",f_p_mom2_vz,"p_mom2_vz[np]/F");
	T->Branch("p_mom2_vt",f_p_mom2_vt,"p_mom2_vt[np]/F");
	//T->Branch("p_mom1_had",b_p_mom1_had,"p_mom1_had[np]/O");
	//T->Branch("p_mom2_had",b_p_mom2_had,"p_mom2_had[np]/O");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		f_scale = pythia.event.scale();
		f_bMPI = pythia.info.bMPI();
		i_nMPI = pythia.info.nMPI();

		i_nmult_v0a = 0;
		i_nmult_v0c = 0;
    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isHadron() && event[i].isFinal()) ) continue;

			int id = event[i].id();
			if ( !(abs(id)==211 || abs(id)==321 || abs(id)==2212) ) continue;

			float eta = event[i].eta();
			if ( eta>-5.1 && eta<-2.8 ){
				i_nmult_v0a++;
			}else if ( eta>1.7 && eta<3.7 ){
				i_nmult_v0c++;
			}
		}

		i_np = 0;
    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isHadron() && event[i].isFinal()) ) continue;

			float eta = event[i].eta();
			if ( fabs(eta)>1.5 ) continue;

			int status 	= event[i].status();
			int id 			= event[i].id();
			if ( !(abs(id)==211 || abs(id)==321 || abs(id)==2212) ) continue;

			float pt	= event[i].pT();
			float phi = event[i].phi();

			i_p_id[i_np] = id;
			i_p_status[i_np] = status;

			//b_p_final[i_np] = event[i].isFinal();

			f_p_pt[i_np] = pt;
			f_p_eta[i_np] = eta;
			f_p_phi[i_np] = phi;
			f_p_vx[i_np] = event[i].xProd();
			f_p_vy[i_np] = event[i].yProd();
			f_p_vz[i_np] = event[i].zProd();
			f_p_vt[i_np] = event[i].tProd();

			int index_mom1 = event[i].mother1();
			int index_mom2 = event[i].mother2();

			i_p_mom1_id[i_np] = event[index_mom1].id();
			i_p_mom2_id[i_np] = event[index_mom2].id();

			i_p_mom1_index[i_np] = index_mom1;
			i_p_mom2_index[i_np] = index_mom2;

			i_p_mom1_status[i_np] = event[index_mom1].status();
			i_p_mom2_status[i_np] = event[index_mom2].status();

			f_p_mom1_vx[i_np] = event[index_mom1].xProd();
			f_p_mom1_vy[i_np] = event[index_mom1].yProd();
			f_p_mom1_vz[i_np] = event[index_mom1].zProd();
			f_p_mom1_vt[i_np] = event[index_mom1].tProd();
			f_p_mom2_vx[i_np] = event[index_mom2].xProd();
			f_p_mom2_vy[i_np] = event[index_mom2].yProd();
			f_p_mom2_vz[i_np] = event[index_mom2].zProd();
			f_p_mom2_vt[i_np] = event[index_mom2].tProd();

			f_p_mom1_mass[i_np] = event[index_mom1].m();
			f_p_mom2_mass[i_np] = event[index_mom2].m();

			//b_p_mom1_had[i_np] = event[index_mom1].isHadron();
			//b_p_mom2_had[i_np] = event[index_mom2].isHadron();

			i_np++;

			if ( i_np>=10000 ) break;
		}

		//cout << "DONE" << endl;

		//if ( nmult_mid>=120 ){
		if ( 1 ){
			T->Fill();
		}

  }//iEvent

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}

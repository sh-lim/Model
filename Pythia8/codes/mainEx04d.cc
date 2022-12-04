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
#include "TH1.h"

//#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
	pythia.readFile("mainEx04d.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	float f_scale;
	float f_bMPI;
	int i_nMPI;

	int i_np;
	int i_mult_v0m;

	int i_p_id[100];
	float f_p_pt[100];
	float f_p_eta[100];
	float f_p_phi[100];
	float f_p_mass[100];

	int i_p1_id[100][30];
	int i_p1_index[100][30];
	int i_p1_status[100][30];

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("scale",&f_scale,"scale/F");
	T->Branch("bMPI",&f_bMPI,"bMPI/F");
	T->Branch("nMPI",&i_nMPI,"nMPI/I");
	T->Branch("multV0M",&i_mult_v0m,"multV0M/I");

	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_mass",f_p_mass,"p_mass[np]/F");

	T->Branch("p1_id",i_p1_id,"p1_id[np][30]/I");
	T->Branch("p1_index",i_p1_index,"p1_index[np][30]/I");
	T->Branch("p1_status",i_p1_status,"p1_status[np][30]/I");

	/*
	T->Branch("p2_id",i_p2_id,"p2_id[np]/I");
	T->Branch("p2_index",i_p2_index,"p2_index[np]/I");
	T->Branch("p3_id",i_p3_id,"p3_id[np]/I");
	T->Branch("p3_index",i_p3_index,"p3_index[np]/I");
	T->Branch("p4_id",i_p4_id,"p4_id[np]/I");
	T->Branch("p4_index",i_p4_index,"p4_index[np]/I");
	*/

	TH1D *hmult_v0m = new TH1D("hmult_v0m","",500,0,500);

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

		//multiplicity calculation
		i_mult_v0m = 0;

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) ) continue;

			float tmp_eta = event[i].eta();

			if ( (tmp_eta>2.8 && tmp_eta<5.1) || (tmp_eta>-3.7 && tmp_eta<-1.7) ){
				i_mult_v0m++;
			}
		}

		hmult_v0m->Fill(i_mult_v0m);

		i_np = 0;

		bool b_e = false;
		bool b_Xi = false;

		for (int i = 0; i < event.size(); i++) {
			if ( !event[i].isFinal() ) continue;

			float eta = event[i].eta();
			if ( fabs(eta)>1.2 ) continue; 

			int id = event[i].id();

			if ( abs(id)==11 ){
				b_e = true;
			}else if ( abs(id)==3312 ){
				b_Xi = true;
			}else{
				continue;
			}

			i_p_id[i_np] = event[i].id();
			f_p_pt[i_np] = event[i].pT();
			f_p_eta[i_np] = event[i].eta();
			f_p_phi[i_np] = event[i].phi();
			f_p_mass[i_np] = event[i].m();

			int p_ind_prev = event[i].mother1();

			for (int ii=0; ii<30; ii++){
				i_p1_id[i_np][ii] = p_ind_prev;
				i_p1_index[i_np][ii] = event[p_ind_prev].id();
				i_p1_status[i_np][ii] = event[p_ind_prev].status();

				p_ind_prev = event[p_ind_prev].mother1();

				if ( p_ind_prev==0 ) break;
			}

			/*
			int p1_index = event[i].mother1();
			int p1_id = event[p1_index].id();
			int p1_status = event[p1_index].status();

			int p2_index = event[p1_index].mother1();
			int p2_id = event[p2_index].id();
			int p2_status = event[p2_index].status();

			int p3_index = event[p2_index].mother1();
			int p3_id = event[p3_index].id();
			int p3_status = event[p3_index].status();

			int p4_index = event[p3_index].mother1();
			int p4_id = event[p4_index].id();
			int p4_status = event[p4_index].status();

			i_p1_id[i_np][0] = p1_id;
			i_p1_index[i_np][0] = p1_index;
			i_p1_status[i_np][0] = p1_status;

			i_p1_id[i_np][1] = p2_id;
			i_p1_index[i_np][1] = p2_index;
			i_p1_status[i_np][1] = p2_status;

			i_p1_id[i_np][2] = p3_id;
			i_p1_index[i_np][2] = p3_index;
			i_p1_status[i_np][2] = p3_status;

			i_p1_id[i_np][3] = p4_id;
			i_p1_index[i_np][3] = p4_index;
			i_p1_status[i_np][3] = p4_status;
			*/

			i_np++;

		}//i

		if ( b_e && b_Xi ){
			T->Fill();
		}

		for (int ii=0; ii<100; ii++){
			i_p_id[ii] = 0;
			f_p_pt[ii] = f_p_eta[ii] = f_p_phi[ii] = f_p_mass[ii] = -999;

			for (int jj=0; jj<30; jj++){
				i_p1_id[ii][jj] = i_p1_index[ii][jj] = i_p1_status[ii][jj] = 0;
			}
		}

  }//iEvent

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	hmult_v0m->Write();
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}

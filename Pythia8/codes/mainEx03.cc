// main03.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how different processes can be selected and studied.
// All input is specified in the main03.cmnd file.
// Also illustrated output to be plotted by Python/Matplotlib/pyplot.

// Pythia for jet fragmentation study

#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"

#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("mainEx03.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	int p_id = 0; 
	float p_pT = -999, p_eta = -999, p_phi = -999;

	int d_id[3];
	float d_pT[3], d_eta[3], d_phi[3];

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("p_id",&p_id,"p_id/I");
	T->Branch("p_pT",&p_pT,"p_pT/F");
	T->Branch("p_eta",&p_eta,"p_eta/F");
	T->Branch("p_phi",&p_phi,"p_phi/F");

	T->Branch("d_id",d_id,"d_id[3]/I");
	T->Branch("d_pT",d_pT,"d_pT[3]/F");
	T->Branch("d_eta",d_eta,"d_eta[3]/F");
	T->Branch("d_phi",d_phi,"d_phi[3]/F");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		p_id = 0; 
		p_pT = p_eta = p_phi = -999;

		for (int i = 0; i < event.size(); ++i) {
			if ( abs(event[i].id())==4132 ){
				p_id = event[i].id();
				p_pT = event[i].pT();
				p_eta = event[i].eta();
				p_phi = event[i].phi();
				cout << iEvent << " " << p_id << " " << p_pT << " " << p_eta << " " << p_phi << endl;
			}
		}

		if ( fabs(p_eta)>5.0 ) continue;

		vector<int> vec_d_id;
		vector<float> vec_d_pT;
		vector<float> vec_d_eta;
		vector<float> vec_d_phi;

		for (int i = 0; i < event.size(); ++i) {

			if ( event[event[i].mother1()].id()==p_id ){
				cout << iEvent << " " << event[i].id() << " " << event[i].pT() << " " << event[i].eta() << " " << event[i].phi() << endl;
				vec_d_id.push_back(event[i].id());
				vec_d_pT.push_back(event[i].pT());
				vec_d_eta.push_back(event[i].eta());
				vec_d_phi.push_back(event[i].phi());
			}

		}

		if ( vec_d_id.size()!=3 ) continue;

		for (int ii=0; ii<3; ii++){
			d_id[ii] = vec_d_id[ii];
			d_pT[ii] = vec_d_pT[ii];
			d_eta[ii] = vec_d_eta[ii];
			d_phi[ii] = vec_d_phi[ii];
		}

		cout << endl;

		T->Fill();

  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}

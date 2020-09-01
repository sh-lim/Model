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

#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
	pythia.readFile("mainEx04.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	int i_np;

	int i_p_id[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];
	bool f_p_charge[5000];

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_charge",f_p_charge,"p_charge[np]/O");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		i_np = 0;

		int nPSI = 0;

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) ) continue;

			int tmp_id = event[i].id();

			if ( !(abs(tmp_id)==11 || abs(tmp_id)==13) ) continue;

			int tmp_pid = event[event[i].mother1()].id();

			if ( tmp_pid==443 || tmp_pid==100443 ){
				nPSI++;
			}
		}

		if ( nPSI<2 ) continue;
		//cout << nPSI << endl;
		//if ( !bPSI ) continue;

		for (int i = 0; i < event.size(); ++i) {

			//if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			//if ( !(event[i].isFinal()) ) continue;

			int tmp_id = event[i].id();

			if ( !(event[i].isFinal() || tmp_id==443 || tmp_id==100443) ) continue;

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();
			bool tmp_charge = event[i].isCharged();

			if ( fabs(tmp_eta)>5.0 ) continue;

			i_p_id[i_np] = tmp_id;
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;
			f_p_charge[i_np] = tmp_charge;

			i_np++;

			if ( i_np>=2000 ) break;

		}

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

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

//#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main() {

	bool bDEBUG = false;

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
	pythia.readFile("mainEx05.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	int i_np;

	int i_p_id[50];
	float f_p_pt[50];
	float f_p_eta[50];
	float f_p_phi[50];
	bool f_p_bb[50];
	bool f_p_cc[50];

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_cc",f_p_cc,"p_cc[np]/O");
	T->Branch("p_bb",f_p_bb,"p_bb[np]/O");

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

		for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) ) continue;

			int pid = event[i].id();
			float pt = event[i].pT();
			float eta = event[i].eta();
			float phi = event[i].phi();

			if ( !(abs(pid)==11 || abs(pid)==13) ) continue;
			if ( fabs(eta)>1.0 ) continue;

			bool bQ = false;
			bool bC = false, bB = false;
			int pidC = 0, pidB = 0;

			int p1_ind = event[i].mother1();
			int p2_ind = event[i].mother2();
			int p1_pid = abs(event[p1_ind].id());
			int p2_pid = abs(event[p2_ind].id());

			if ( bDEBUG ){
				cout << p1_ind << " " << p1_pid << " " << p2_ind << " " << p2_pid << endl;
			}

			if ( p1_pid==4 || (p1_pid>400 && p1_pid<500) || (p1_pid>4000 && p1_pid<5000) ){
				pidC = p1_pid;
				bC = true;
			}if ( p1_pid==5 || (p1_pid>500 && p1_pid<600) || (p1_pid>5000 && p1_pid<6000) ){
				pidB = p1_pid;
				bB = true;
			}

			if ( p2_ind!=0 || (p1_ind==0 && p2_ind==0) ){
				bQ = true;
			}

			while ( !bQ ){

				int pp1_ind = event[p1_ind].mother1();
				int pp2_ind = event[p1_ind].mother2();

				int pp1_pid = abs(event[pp1_ind].id());
				int pp2_pid = abs(event[pp2_ind].id());
				
				if ( bDEBUG ){
					cout << pp1_ind << " " << pp1_pid << " " << pp2_ind << " " << pp2_pid << endl;
				}

				if ( pp1_pid==4 || pp2_pid==4 || (pp1_pid>400 && pp1_pid<500) || (pp1_pid>4000 && pp1_pid<5000) ){
					pidC = pp1_pid;
					bC = true;
				}if ( pp1_pid==5 || pp2_pid==5 || (pp1_pid>500 && pp1_pid<600) || (pp1_pid>5000 && pp1_pid<6000) ){
					pidB = pp1_pid;
					bB = true;
				}

				if ( pp2_ind!=0 || (pp1_ind==0 && pp2_ind==0) ){
					bQ = true;
				}else{
					p1_ind = pp1_ind;
					p2_ind = pp2_ind;
				}

			}//bQ

			i_p_id[i_np] = pid;
			f_p_pt[i_np] = pt;
			f_p_eta[i_np] = eta;
			f_p_phi[i_np] = phi;

			f_p_cc[i_np] = bC;
			f_p_bb[i_np] = bB;

			i_np++;

		}//iEvent

		if ( i_np>0 ){
			T->Fill();
		}

  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}

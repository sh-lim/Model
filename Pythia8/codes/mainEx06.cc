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
  pythia.readFile("mainEx06.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

  // Set up FastJet jet finder.
  //   one can use either explicitly use antikt, cambridge, etc., or
  //   just use genkt_algorithm with specification of power
  //fastjet::JetAlgorithm algorithm;
  //if (power == -1)      algorithm = fastjet::antikt_algorithm;
  //if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  //if (power ==  1)      algorithm = fastjet::kt_algorithm;
  //fastjet::JetDefinition jetDef(algorithm, R);
  // there's no need for a pointer to the jetDef (it's a fairly small object)
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, 0.4, -1);
  std::vector <fastjet::PseudoJet> fjInputs;

	int i_np;

	int i_p_id[10000];
	float f_p_pt[10000];
	float f_p_eta[10000];
	float f_p_phi[10000];

	float f_p_vx[10000];
	float f_p_vy[10000];
	float f_p_vz[10000];
	float f_p_vt[10000];

	//const float const_pt_cut = 0.2;

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");

	T->Branch("p_vx",f_p_vx,"p_vx[np]/F");
	T->Branch("p_vy",f_p_vy,"p_vy[np]/F");
	T->Branch("p_vz",f_p_vz,"p_vz[np]/F");
	T->Branch("p_vt",f_p_vt,"p_vt[np]/F");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		int nmult_mid = 0, nmult_fwd = 0;

		i_np = 0;
    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;

			i_p_id[i_np] = event[i].id();

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();

			//if ( fabs(tmp_eta)>2.5 ) continue;
			//if ( tmp_pt<0.3 ) continue;
			if ( fabs(tmp_eta)>5.0 ) continue;
			if ( fabs(tmp_eta)<2.5 && tmp_pt<0.2 ) continue;

			float xprod = event[i].xProd();
			float yprod = event[i].yProd();
			float zprod = event[i].zProd();
			float tprod = event[i].tProd();

			if ( tprod<0.1 ){
				if ( fabs(tmp_eta)<0.9 && tmp_pt>0.4 ){
					nmult_mid++;
				}

				if ( tmp_eta>-5.1 && tmp_eta<-2.8 ){
					nmult_fwd++;
				}
			}

			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			f_p_vx[i_np] = xprod;
			f_p_vy[i_np] = yprod;
			f_p_vz[i_np] = zprod;
			f_p_vt[i_np] = tprod;

			i_np++;

			if ( i_np>=10000 ) break;

		}

		//string shoving
		//if ( nmult_mid>=48 || nmult_fwd>=125 ){
		//default PYTHIA
		if ( nmult_mid>=36 || nmult_fwd>=89 ){
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

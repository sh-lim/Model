// main03.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how different processes can be selected and studied.
// All input is specified in the main03.cmnd file.
// Also illustrated output to be plotted by Python/Matplotlib/pyplot.

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
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
	pythia.readFile("mainEx00a.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	float f_b;

	int i_nAbsProj;
	int i_nDiffProj;
	int i_nElProj;

	int i_nAbsTarg;
	int i_nDiffTarg;
	int i_nElTarg;
	
	int i_nCollNDTot;
	int i_nCollND;
	int i_nCollSDP;
	int i_nCollSDT;
	int i_nCollDD;
	int i_nCollCD;
	int i_nCollEL;

	/*
	int i_p_id[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];
	float f_p_vt[5000];
	*/

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("f_b",&f_b,"f_b/F");

	T->Branch("i_nAbsProj",&i_nAbsProj,"i_nAbsProj/I");
	T->Branch("i_nDiffProj",&i_nDiffProj,"i_nDiffProj/I");
	T->Branch("i_nElProj",&i_nElProj,"i_nElProj/I");

	T->Branch("i_nAbsTarg",&i_nAbsTarg,"i_nAbsTarg/I");
	T->Branch("i_nDiffTarg",&i_nDiffTarg,"i_nDiffTarg/I");
	T->Branch("i_nElTarg",&i_nElTarg,"i_nElTarg/I");

	T->Branch("i_nCollNDTot",&i_nCollNDTot,"i_nCollNDTot/I");
	T->Branch("i_nCollND",&i_nCollND,"i_nCollND/I");
	T->Branch("i_nCollSDP",&i_nCollSDP,"i_nCollSDP/I");
	T->Branch("i_nCollSDT",&i_nCollSDT,"i_nCollSDT/I");
	T->Branch("i_nCollDD",&i_nCollDD,"i_nCollDD/I");
	T->Branch("i_nCollCD",&i_nCollCD,"i_nCollCD/I");
	T->Branch("i_nCollEL",&i_nCollEL,"i_nCollEL/I");

	/*
	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_vt",f_p_vt,"p_vt[np]/F");
	*/

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		f_b = pythia.info.hiInfo->b();

		i_nAbsProj = pythia.info.hiInfo->nAbsProj();
		i_nDiffProj = pythia.info.hiInfo->nDiffProj();
		i_nElProj = pythia.info.hiInfo->nElProj();

		i_nAbsTarg = pythia.info.hiInfo->nAbsTarg();
		i_nDiffTarg = pythia.info.hiInfo->nDiffTarg();
		i_nElTarg = pythia.info.hiInfo->nElTarg();

		i_nCollNDTot = pythia.info.hiInfo->nCollNDTot();
		i_nCollND = pythia.info.hiInfo->nCollND();
		i_nCollSDP = pythia.info.hiInfo->nCollSDP();
		i_nCollSDT = pythia.info.hiInfo->nCollSDT();
		i_nCollDD = pythia.info.hiInfo->nCollDD();
		i_nCollCD = pythia.info.hiInfo->nCollCD();
		i_nCollEL = pythia.info.hiInfo->nCollEL();

		/*
		i_np = 0;

    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) ) continue;

			int id = abs(event[i].id());

			if ( !(id==111 || id==211 || id==321 || id==2212) ) continue;

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();

			if ( fabs(tmp_eta)>5.5 ) continue;

			i_p_id[i_np] = event[i].id();
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			f_p_vt[i_np] = tprod;

			i_np++;

			if ( i_np>=2000 ) break;
		}
		*/

		T->Fill();

  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event_mainEx00a.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}

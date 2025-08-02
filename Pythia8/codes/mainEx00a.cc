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
#include "TLorentzVector.h"

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

  int i_np;
  int i_p_id[5000];
  float f_p_pt[5000];
  float f_p_eta[5000];
  float f_p_phi[5000];
  float f_p_rap[5000];
  float f_p_vt[5000];
  //int i_num_dmeson;


	TFile *outfile = new TFile("Pythia8_event.root","recreate");

  // Tree output
  auto T = new TTree("T","Pythia event");

  T->Branch("np",&i_np,"np/I");
  T->Branch("p_id",i_p_id,"p_id[np]/I");
  T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
  T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
  T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
  T->Branch("p_rap",f_p_rap,"p_rap[np]/F");
  T->Branch("p_vt",f_p_vt,"p_vt[np]/F");
  //T->Branch("num_dmeson", &i_num_dmeson, "num_dmeson/I");

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
	//i_num_dmeson=0;

    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) ) continue;
	
			//int id = abs(event[i].id());

			// Count the number of D-meson(D0, D+, D-) exising in this event
			/*
			if ( id==421 || id==411 )
			{
				i_num_dmeson++;
			}
			*/

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();
			float tmp_mass = event[i].m();
			float tmp_t = event[i].tProd();

			TLorentzVector lvec;
			lvec.SetPtEtaPhiM(tmp_pt, tmp_eta, tmp_phi, tmp_mass);

			float tmp_rap = lvec.Rapidity();


			if ( fabs(tmp_rap)>1.0 ) continue;
			if ( tmp_t>1.0 ) continue;
			//if ( fabs(tmp_eta)>6.0 ) continue;

			i_p_id[i_np] = event[i].id();
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;
			f_p_rap[i_np] = tmp_rap;
			f_p_vt[i_np] = tmp_t;

			i_np++;

			if ( i_np>5000 ) break;
		}
	
		T->Fill();

		if ((iEvent) % 10000 == 0) {
			T->AutoSave("SaveSelf");    // 중간 저장
			T->FlushBaskets();          // 버퍼 비우기
		}
  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

  //TFile *outfile = new TFile("Pythia8_event_mainEx00a.root","recreate");
  T->Write();
  outfile->Close();

  // Done.
  return 0;
}

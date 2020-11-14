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

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("mainEx01.cfg");

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
  //fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, 0.4, -1);
  //std::vector <fastjet::PseudoJet> fjInputs;

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np;
	//int i_njet;

	int i_p_id[2000];
	int i_p_status[2000];
	float f_p_pt[2000];
	float f_p_eta[2000];
	float f_p_phi[2000];
	bool b_p_final[2000]; 
	int i_p_mom1[2000];
	int i_p_mom2[2000];
	bool b_p_mom1_had[2000]; 
	bool b_p_mom2_had[2000]; 

	float f_p_vt[2000];

	/*
	float f_jet_pt[100];
	float f_jet_eta[100];
	float f_jet_phi[100];
	*/

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("scale",&f_scale,"scale/F");
	T->Branch("bMPI",&f_bMPI,"bMPI/F");
	T->Branch("nMPI",&i_nMPI,"nMPI/I");

	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_status",i_p_status,"p_status[np]/I");
	T->Branch("p_final",b_p_final,"p_final[np]/O");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_vt",f_p_vt,"p_vt[np]/F");

	T->Branch("p_mom1",i_p_mom1,"p_mom1[np]/I");
	T->Branch("p_mom2",i_p_mom2,"p_mom2[np]/I");
	T->Branch("p_mom1_had",b_p_mom1_had,"p_mom1_had[np]/O");
	T->Branch("p_mom2_had",b_p_mom2_had,"p_mom2_had[np]/O");

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

		/*
		fjInputs.resize(0);
    for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			if ( fabs(event[i].eta())>3.0 ) continue;

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp = event[i];

			fjInputs.push_back(particleTemp);
		}

		// Run Fastjet algorithm and sort jets in pT order.
		vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
		fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
		inclusiveJets = clustSeq.inclusive_jets(5.0);
		sortedJets = sorted_by_pt(inclusiveJets);
		*/

		i_np = 0;
    for (int i = 0; i < event.size(); ++i) {

			//if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			if ( !(event[i].isHadron()) || !(event[i].isCharged()) ) continue;

			int status 	= event[i].status();
			int id 			= event[i].id();

			float pt	= event[i].pT();
			float eta = event[i].eta();
			float phi = event[i].phi();
			float tprod = event[i].tProd();

			if ( fabs(eta)>5.0 ) continue;
			if ( fabs(eta)<2.5 && pt<0.2 ) continue;

			int index_mom1 = event[i].mother1();
			int index_mom2 = event[i].mother2();

			int id_mom1 = event[index_mom1].id();
			int id_mom2 = event[index_mom2].id();

			i_p_id[i_np] = id;
			i_p_status[i_np] = status;

			b_p_final[i_np] = event[i].isFinal();

			f_p_pt[i_np] = pt;
			f_p_eta[i_np] = eta;
			f_p_phi[i_np] = phi;
			f_p_vt[i_np] = tprod;

			i_p_mom1[i_np] = id_mom1;
			i_p_mom2[i_np] = id_mom2;

			b_p_mom1_had[i_np] = event[index_mom1].isHadron();
			b_p_mom2_had[i_np] = event[index_mom2].isHadron();

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

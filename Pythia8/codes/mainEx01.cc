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
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, 0.4, -1);
  std::vector <fastjet::PseudoJet> fjInputs;

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np;
	int i_njet;

	int i_p_id[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];

	float f_jet_pt[100];
	float f_jet_eta[100];
	float f_jet_phi[100];

	const float const_pt_cut = 0.2;

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("scale",&f_scale,"scale/F");
	T->Branch("bMPI",&f_bMPI,"bMPI/F");
	T->Branch("nMPI",&i_nMPI,"nMPI/I");

	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");

	T->Branch("njet",&i_njet,"njet/I");
	T->Branch("jet_pt",f_jet_pt,"jet_pt[njet]/F");
	T->Branch("jet_eta",f_jet_eta,"jet_eta[njet]/F");
	T->Branch("jet_phi",f_jet_phi,"jet_phi[njet]/F");

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
		cout 
			<< pythia.event.scale() << " " 
			<< pythia.info.bMPI() << " " 
			<< pythia.info.nMPI() << " " 
			//<< pythia.info.nISR() << " " 
			//<< pythia.info.nFSRinProc() << " " 
			//<< pythia.info.nFSRinRes() << " " 
			<< endl;
		*/

		fjInputs.resize(0);
    for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			if ( fabs(event[i].eta())>1.0 ) continue;

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp = event[i];

			fjInputs.push_back(particleTemp);
		}

		// Run Fastjet algorithm and sort jets in pT order.
		vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
		fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
		inclusiveJets = clustSeq.inclusive_jets(10.0);
		sortedJets = sorted_by_pt(inclusiveJets);

		//cout << "# of jets: " << sortedJets.size() << endl;

		//continue;
		if ( sortedJets.size()<1 ) continue;

		i_njet = 0;
		for (int i = 0; i < int(sortedJets.size()); ++i) {
			f_jet_pt[i_njet] = sortedJets[i].perp(); 
			f_jet_eta[i_njet] = sortedJets[i].rap(); 
			f_jet_phi[i_njet] = sortedJets[i].phi_std(); 
			i_njet++;
		}

		i_np = 0;
    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;

			i_p_id[i_np] = event[i].id();

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();

			if ( fabs(tmp_eta)>6.0 ) continue;
			if ( fabs(tmp_eta)<1.5 && tmp_pt<const_pt_cut ) continue;

			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			/*
			int status = event[i].status();
			int id = event[i].id();

			double xprod = event[i].xProd();
			double yprod = event[i].yProd();
			double zprod = event[i].zProd();
			double tprod = event[i].tProd();
			double rapidity_tau = (tprod-zprod)<1e-10 ? 0 : 0.5*log((tprod+zprod)/(tprod-zprod));

			double px = event[i].px();
			double py = event[i].py();
			double pz = event[i].pz();
			double ee = event[i].e();
			double rapidity = (ee-fabs(pz))<1e-10 ? 0 : 0.5*log((ee+pz)/(ee-pz));
			*/

			i_np++;

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

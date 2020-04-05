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
  pythia.readFile("mainEx02.cfg");

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
  fastjet::JetDefinition jetDefR03(fastjet::genkt_algorithm, 0.3, -1);
  fastjet::JetDefinition jetDefR04(fastjet::genkt_algorithm, 0.4, -1);
  fastjet::JetDefinition jetDefR05(fastjet::genkt_algorithm, 0.5, -1);
  std::vector <fastjet::PseudoJet> fjInputs;

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np;
	int i_njetR03;
	int i_njetR04;
	int i_njetR05;

	int i_p_id[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];

	float f_jetR03_pt[10];
	float f_jetR03_eta[10];
	float f_jetR03_phi[10];
	int f_jetR03_nconst[10];

	float f_jetR03_const_pt[10][30];
	float f_jetR03_const_eta[10][30];
	float f_jetR03_const_phi[10][30];

	float f_jetR04_pt[10];
	float f_jetR04_eta[10];
	float f_jetR04_phi[10];
	int f_jetR04_nconst[10];

	float f_jetR04_const_pt[10][30];
	float f_jetR04_const_eta[10][30];
	float f_jetR04_const_phi[10][30];

	float f_jetR05_pt[10];
	float f_jetR05_eta[10];
	float f_jetR05_phi[10];
	int f_jetR05_nconst[10];

	float f_jetR05_const_pt[10][30];
	float f_jetR05_const_eta[10][30];
	float f_jetR05_const_phi[10][30];

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

	T->Branch("njetR03",&i_njetR03,"njetR03/I");
	T->Branch("jetR03_pt",f_jetR03_pt,"jetR03_pt[njetR03]/F");
	T->Branch("jetR03_eta",f_jetR03_eta,"jetR03_eta[njetR03]/F");
	T->Branch("jetR03_phi",f_jetR03_phi,"jetR03_phi[njetR03]/F");
	T->Branch("jetR03_nconst",f_jetR03_nconst,"jetR03_nconst[njetR03]/I");
	T->Branch("jetR03_const_pt",f_jetR03_const_pt,"jetR03_const_pt[njetR03][100]/F");
	T->Branch("jetR03_const_eta",f_jetR03_const_eta,"jetR03_const_eta[njetR03][100]/F");
	T->Branch("jetR03_const_phi",f_jetR03_const_phi,"jetR03_const_phi[njetR03][100]/F");

	T->Branch("njetR04",&i_njetR04,"njetR04/I");
	T->Branch("jetR04_pt",f_jetR04_pt,"jetR04_pt[njetR04]/F");
	T->Branch("jetR04_eta",f_jetR04_eta,"jetR04_eta[njetR04]/F");
	T->Branch("jetR04_phi",f_jetR04_phi,"jetR04_phi[njetR04]/F");
	T->Branch("jetR04_nconst",f_jetR04_nconst,"jetR04_nconst[njetR04]/I");
	T->Branch("jetR04_const_pt",f_jetR04_const_pt,"jetR04_const_pt[njetR04][100]/F");
	T->Branch("jetR04_const_eta",f_jetR04_const_eta,"jetR04_const_eta[njetR04][100]/F");
	T->Branch("jetR04_const_phi",f_jetR04_const_phi,"jetR04_const_phi[njetR04][100]/F");

	T->Branch("njetR05",&i_njetR05,"njetR05/I");
	T->Branch("jetR05_pt",f_jetR05_pt,"jetR05_pt[njetR05]/F");
	T->Branch("jetR05_eta",f_jetR05_eta,"jetR05_eta[njetR05]/F");
	T->Branch("jetR05_phi",f_jetR05_phi,"jetR05_phi[njetR05]/F");
	T->Branch("jetR05_nconst",f_jetR05_nconst,"jetR05_nconst[njetR05]/I");
	T->Branch("jetR05_const_pt",f_jetR05_const_pt,"jetR05_const_pt[njetR05][100]/F");
	T->Branch("jetR05_const_eta",f_jetR05_const_eta,"jetR05_const_eta[njetR05][100]/F");
	T->Branch("jetR05_const_phi",f_jetR05_const_phi,"jetR05_const_phi[njetR05][100]/F");

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

		//multiplicity calculation 
		int i_mult_fwd = 0;
		for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			float tmp_eta = event[i].eta();

			if ( (tmp_eta>2.8 && tmp_eta<5.1) || (tmp_eta>-3.7 && tmp_eta<-1.7) ){
				i_mult_fwd++;
			}//

		}

		fjInputs.resize(0);
    for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			//ATLAS acceptance
			if ( fabs(event[i].eta())>2.0 ) continue;

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp = event[i];

			fjInputs.push_back(particleTemp);
		}

		// Run Fastjet algorithm and sort jets in pT order.
		vector <fastjet::PseudoJet> inclusiveJetsR03, sortedJetsR03;
		vector <fastjet::PseudoJet> inclusiveJetsR04, sortedJetsR04;
		vector <fastjet::PseudoJet> inclusiveJetsR05, sortedJetsR05;

		fastjet::ClusterSequence clustSeqR03(fjInputs, jetDefR03);
		inclusiveJetsR03 = clustSeqR03.inclusive_jets(20.0);
		sortedJetsR03 = sorted_by_pt(inclusiveJetsR03);

		fastjet::ClusterSequence clustSeqR04(fjInputs, jetDefR04);
		inclusiveJetsR04 = clustSeqR04.inclusive_jets(20.0);
		sortedJetsR04 = sorted_by_pt(inclusiveJetsR04);

		fastjet::ClusterSequence clustSeqR05(fjInputs, jetDefR05);
		inclusiveJetsR05 = clustSeqR05.inclusive_jets(20.0);
		sortedJetsR05 = sorted_by_pt(inclusiveJetsR05);

		//cout << "# of jets: " << sortedJets.size() << endl;

		//continue;
		if ( (sortedJetsR03.size()+sortedJetsR04.size()+sortedJetsR05.size())<1 ) continue;

		i_njetR03 = 0;
		i_njetR04 = 0;
		i_njetR05 = 0;

		for (int ii=0; ii<10; ii++){
			f_jetR03_pt[ii] = f_jetR03_eta[ii] = f_jetR03_phi[ii] = -999;
			f_jetR04_pt[ii] = f_jetR04_eta[ii] = f_jetR04_phi[ii] = -999;
			f_jetR05_pt[ii] = f_jetR05_eta[ii] = f_jetR05_phi[ii] = -999;

			f_jetR03_nconst[ii] = f_jetR04_nconst[ii] = f_jetR05_nconst[ii] = 0;

			for (int jj=0; jj<30; jj++){
				f_jetR03_const_pt[ii][jj] = f_jetR03_const_eta[ii][jj] = f_jetR03_const_phi[ii][jj] = -999;
				f_jetR04_const_pt[ii][jj] = f_jetR04_const_eta[ii][jj] = f_jetR04_const_phi[ii][jj] = -999;
				f_jetR05_const_pt[ii][jj] = f_jetR05_const_eta[ii][jj] = f_jetR05_const_phi[ii][jj] = -999;
			}
		}

		for (int i = 0; i < int(sortedJetsR03.size()); ++i) {
			f_jetR03_pt[i_njetR03] = sortedJetsR03[i].perp(); 
			f_jetR03_eta[i_njetR03] = sortedJetsR03[i].rap(); 
			f_jetR03_phi[i_njetR03] = sortedJetsR03[i].phi_std(); 

			vector<fastjet::PseudoJet> constituents = sortedJetsR03[i].constituents();

			f_jetR03_nconst[i_njetR03] = constituents.size();
			for (int j=0; j<f_jetR03_nconst[i_njetR03]; ++j) {
				if ( j>=30 ) break;
				f_jetR03_const_pt[i][j] = constituents[j].perp();
				f_jetR03_const_eta[i][j] = constituents[j].rap();
				f_jetR03_const_phi[i][j] = constituents[j].phi_std();
			}

			i_njetR03++;
		}

		for (int i = 0; i < int(sortedJetsR04.size()); ++i) {
			f_jetR04_pt[i_njetR04] = sortedJetsR04[i].perp(); 
			f_jetR04_eta[i_njetR04] = sortedJetsR04[i].rap(); 
			f_jetR04_phi[i_njetR04] = sortedJetsR04[i].phi_std(); 

			vector<fastjet::PseudoJet> constituents = sortedJetsR04[i].constituents();

			f_jetR04_nconst[i_njetR04] = constituents.size();
			for (int j=0; j<f_jetR04_nconst[i_njetR04]; ++j) {
				if ( j>=30 ) break;
				f_jetR04_const_pt[i][j] = constituents[j].perp();
				f_jetR04_const_eta[i][j] = constituents[j].rap();
				f_jetR04_const_phi[i][j] = constituents[j].phi_std();
			}

			i_njetR04++;
		}

		for (int i = 0; i < int(sortedJetsR05.size()); ++i) {
			f_jetR05_pt[i_njetR05] = sortedJetsR05[i].perp(); 
			f_jetR05_eta[i_njetR05] = sortedJetsR05[i].rap(); 
			f_jetR05_phi[i_njetR05] = sortedJetsR05[i].phi_std(); 

			vector<fastjet::PseudoJet> constituents = sortedJetsR05[i].constituents();

			f_jetR05_nconst[i_njetR05] = constituents.size();
			for (int j=0; j<f_jetR05_nconst[i_njetR05]; ++j) {
				if ( j>=30 ) break;
				f_jetR05_const_pt[i][j] = constituents[j].perp();
				f_jetR05_const_eta[i][j] = constituents[j].rap();
				f_jetR05_const_phi[i][j] = constituents[j].phi_std();
			}

			i_njetR05++;
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

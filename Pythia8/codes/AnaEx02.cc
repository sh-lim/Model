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
#include "TLorentzVector.h"

#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main(int argc, char *argv[]) {
//int main() {

	if ( argc<2 ){
		cout << "Usage: ./AnaEx02 cut_jetpt" << endl;
		return -1;
	}

	float cut_jetpt = atof(argv[1]);
	cout << "Cut Jet pT: " << cut_jetpt << endl;

	int i_np;
	int i_njetR04;
	int i_nchgjetR04;

	/*
	int i_njetR03;
	int i_njetR05;
	int i_nchgjetR03;
	int i_nchgjetR05;
	*/

	int i_p_id[5000];
	bool b_p_chg[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];
	float f_p_e[5000];
	float f_p_m[5000];

	float f_jetR04_pt[10];
	float f_jetR04_eta[10];
	float f_jetR04_phi[10];

	float f_chgjetR04_pt[10];
	float f_chgjetR04_eta[10];
	float f_chgjetR04_phi[10];

	/*
	float f_jetR03_pt[10];
	float f_jetR03_eta[10];
	float f_jetR03_phi[10];

	float f_jetR05_pt[10];
	float f_jetR05_eta[10];
	float f_jetR05_phi[10];

	float f_chgjetR03_pt[10];
	float f_chgjetR03_eta[10];
	float f_chgjetR03_phi[10];

	float f_chgjetR05_pt[10];
	float f_chgjetR05_eta[10];
	float f_chgjetR05_phi[10];
	*/

	TFile *infile = new TFile("Pythia8_event.root","read");

	TTree *T = (TTree*)infile->Get("T");

	T->SetBranchAddress("np",&i_np);
	T->SetBranchAddress("p_id",i_p_id);
	T->SetBranchAddress("p_chg",b_p_chg);
	T->SetBranchAddress("p_pt",f_p_pt);
	T->SetBranchAddress("p_eta",f_p_eta);
	T->SetBranchAddress("p_phi",f_p_phi);
	T->SetBranchAddress("p_e",f_p_e);
	T->SetBranchAddress("p_m",f_p_m);

	T->SetBranchAddress("njetR04",&i_njetR04);
	T->SetBranchAddress("jetR04_pt",f_jetR04_pt);
	T->SetBranchAddress("jetR04_eta",f_jetR04_eta);
	T->SetBranchAddress("jetR04_phi",f_jetR04_phi);

	T->SetBranchAddress("nchgjetR04",&i_nchgjetR04);
	T->SetBranchAddress("chgjetR04_pt",f_chgjetR04_pt);
	T->SetBranchAddress("chgjetR04_eta",f_chgjetR04_eta);
	T->SetBranchAddress("chgjetR04_phi",f_chgjetR04_phi);

	int nentries = T->GetEntries();

	cout << "Number of events: " << nentries << endl;


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
  fastjet::JetDefinition jetDefR06(fastjet::genkt_algorithm, 0.6, -1);
  fastjet::JetDefinition jetDefR07(fastjet::genkt_algorithm, 0.7, -1);
  std::vector <fastjet::PseudoJet> fjInputs;
  std::vector <fastjet::PseudoJet> fjInputsChg;

	for (int ien=0; ien<nentries; ien++){

		T->GetEntry(ien);

		fjInputs.resize(0);
		fjInputsChg.resize(0);

		//Jet reconstruction
		for (int ip=0; ip<i_np; ip++){

			//ALICE acceptance
			if ( fabs(f_p_eta[ip])>1.0 ) continue;

			int id = abs(i_p_id[ip]);
			if ( id==12 || id==14 || id==16 ) continue;

			TLorentzVector lvec;
			lvec.SetPtEtaPhiE(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], f_p_e[ip]);

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.E());

			fjInputs.push_back(particleTemp);

			if ( !(b_p_chg[ip]) ) continue;

			fjInputsChg.push_back(particleTemp);

		}//ip

		vector <fastjet::PseudoJet> inclusiveJetsR04, sortedJetsR04;
		vector <fastjet::PseudoJet> inclusiveChgJetsR04, sortedChgJetsR04;

		fastjet::ClusterSequence clustSeqR04(fjInputs, jetDefR04);
		inclusiveJetsR04 = clustSeqR04.inclusive_jets(cut_jetpt);
		sortedJetsR04 = sorted_by_pt(inclusiveJetsR04);

		fastjet::ClusterSequence clustSeqChgR04(fjInputsChg, jetDefR04);
		inclusiveChgJetsR04 = clustSeqChgR04.inclusive_jets(cut_jetpt);
		sortedChgJetsR04 = sorted_by_pt(inclusiveChgJetsR04);

		if ( (i_njetR04!=int(sortedJetsR04.size())) || (i_nchgjetR04!=int(sortedChgJetsR04.size())) ){
			cout << "Number of full jets: " << i_njetR04 << ", " << sortedJetsR04.size() << endl;
			cout << "Number of charged jets: " << i_nchgjetR04 << ", " << sortedChgJetsR04.size() << endl;
			continue;
		}

		for (int i = 0; i < int(sortedJetsR04.size()); ++i) {
			cout << "full jetR04: " << i << endl; 
			cout 
				<< "(" << f_jetR04_pt[i] 
				<< "," << f_jetR04_eta[i]
				<< "," << f_jetR04_phi[i]
				<< ")  "
				<< "(" << sortedJetsR04[i].perp() 
				<< "," << sortedJetsR04[i].rap()
				<< "," << sortedJetsR04[i].phi_std()
				<< ")  "
				<< endl;
		}

		for (int i = 0; i < int(sortedChgJetsR04.size()); ++i) {
			cout << "charged jetR04: " << i << endl; 
			cout 
				<< "(" << f_chgjetR04_pt[i] 
				<< "," << f_chgjetR04_eta[i]
				<< "," << f_chgjetR04_phi[i]
				<< ")  "
				<< "(" << sortedChgJetsR04[i].perp() 
				<< "," << sortedChgJetsR04[i].rap()
				<< "," << sortedChgJetsR04[i].phi_std()
				<< ")  "
				<< endl;
		}

	}//ien


	//Done
	return 0;

	/*

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


	//float f_scale;
	//float f_bMPI;
	//int i_nMPI;

	//const float const_pt_cut = 0.2;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

		// Run Fastjet algorithm and sort jets in pT order.
		vector <fastjet::PseudoJet> inclusiveJetsR03, sortedJetsR03;
		vector <fastjet::PseudoJet> inclusiveJetsR05, sortedJetsR05;

		vector <fastjet::PseudoJet> inclusiveChgJetsR03, sortedChgJetsR03;
		vector <fastjet::PseudoJet> inclusiveChgJetsR05, sortedChgJetsR05;

		fastjet::ClusterSequence clustSeqR03(fjInputs, jetDefR03);
		inclusiveJetsR03 = clustSeqR03.inclusive_jets(cut_jetpt);
		sortedJetsR03 = sorted_by_pt(inclusiveJetsR03);


		fastjet::ClusterSequence clustSeqR05(fjInputs, jetDefR05);
		inclusiveJetsR05 = clustSeqR05.inclusive_jets(cut_jetpt);
		sortedJetsR05 = sorted_by_pt(inclusiveJetsR05);

		fastjet::ClusterSequence clustSeqChgR03(fjInputsChg, jetDefR03);
		inclusiveChgJetsR03 = clustSeqChgR03.inclusive_jets(cut_jetpt);
		sortedChgJetsR03 = sorted_by_pt(inclusiveChgJetsR03);

		fastjet::ClusterSequence clustSeqChgR05(fjInputsChg, jetDefR05);
		inclusiveChgJetsR05 = clustSeqChgR05.inclusive_jets(cut_jetpt);
		sortedChgJetsR05 = sorted_by_pt(inclusiveChgJetsR05);

		//cout << "# of jets: " << sortedJets.size() << endl;

		int njet = sortedJetsR03.size() + sortedChgJetsR03.size();
		njet += sortedJetsR04.size() + sortedChgJetsR04.size();
		njet += sortedJetsR05.size() + sortedChgJetsR05.size();

		//continue;
		//if ( (sortedJetsR03.size()+sortedJetsR04.size()+sortedJetsR05.size()+sortedJetsR06.size()+sortedJetsR07.size())<1 ) continue;
		if ( njet<1 ) continue;

		i_njetR03 = 0;
		i_njetR04 = 0;
		i_njetR05 = 0;

		i_nchgjetR03 = 0;
		i_nchgjetR04 = 0;
		i_nchgjetR05 = 0;

		for (int ii=0; ii<10; ii++){
			f_jetR03_pt[ii] = f_jetR03_eta[ii] = f_jetR03_phi[ii] = -999;
			f_jetR04_pt[ii] = f_jetR04_eta[ii] = f_jetR04_phi[ii] = -999;
			f_jetR05_pt[ii] = f_jetR05_eta[ii] = f_jetR05_phi[ii] = -999;

			f_chgjetR03_pt[ii] = f_chgjetR03_eta[ii] = f_chgjetR03_phi[ii] = -999;
			f_chgjetR04_pt[ii] = f_chgjetR04_eta[ii] = f_chgjetR04_phi[ii] = -999;
			f_chgjetR05_pt[ii] = f_chgjetR05_eta[ii] = f_chgjetR05_phi[ii] = -999;
		}

		for (int i = 0; i < int(sortedJetsR03.size()); ++i) {
			f_jetR03_pt[i_njetR03] = sortedJetsR03[i].perp(); 
			f_jetR03_eta[i_njetR03] = sortedJetsR03[i].rap(); 
			f_jetR03_phi[i_njetR03] = sortedJetsR03[i].phi_std(); 

			i_njetR03++;
		}

		for (int i = 0; i < int(sortedJetsR04.size()); ++i) {
			f_jetR04_pt[i_njetR04] = sortedJetsR04[i].perp(); 
			f_jetR04_eta[i_njetR04] = sortedJetsR04[i].rap(); 
			f_jetR04_phi[i_njetR04] = sortedJetsR04[i].phi_std(); 

			i_njetR04++;
		}

		for (int i = 0; i < int(sortedJetsR05.size()); ++i) {
			f_jetR05_pt[i_njetR05] = sortedJetsR05[i].perp(); 
			f_jetR05_eta[i_njetR05] = sortedJetsR05[i].rap(); 
			f_jetR05_phi[i_njetR05] = sortedJetsR05[i].phi_std(); 

			i_njetR05++;
		}

		for (int i = 0; i < int(sortedChgJetsR03.size()); ++i) {
			f_chgjetR03_pt[i_nchgjetR03] = sortedChgJetsR03[i].perp(); 
			f_chgjetR03_eta[i_nchgjetR03] = sortedChgJetsR03[i].rap(); 
			f_chgjetR03_phi[i_nchgjetR03] = sortedChgJetsR03[i].phi_std(); 

			i_nchgjetR03++;
		}

		for (int i = 0; i < int(sortedChgJetsR04.size()); ++i) {
			f_chgjetR04_pt[i_nchgjetR04] = sortedChgJetsR04[i].perp(); 
			f_chgjetR04_eta[i_nchgjetR04] = sortedChgJetsR04[i].rap(); 
			f_chgjetR04_phi[i_nchgjetR04] = sortedChgJetsR04[i].phi_std(); 

			i_nchgjetR04++;
		}

		for (int i = 0; i < int(sortedChgJetsR05.size()); ++i) {
			f_chgjetR05_pt[i_nchgjetR05] = sortedChgJetsR05[i].perp(); 
			f_chgjetR05_eta[i_nchgjetR05] = sortedChgJetsR05[i].rap(); 
			f_chgjetR05_phi[i_nchgjetR05] = sortedChgJetsR05[i].phi_std(); 

			i_nchgjetR05++;
		}


		i_np = 0;
		for (int i = 0; i < event.size(); ++i) {

			//if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;
			if ( !(event[i].isFinal()) ) continue;

			i_p_id[i_np] = event[i].id();

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();

			if ( fabs(tmp_eta)>1.0 ) continue;
			//if ( fabs(tmp_eta)<2.5 && tmp_pt<const_pt_cut ) continue;

			i_p_id[i_np] = event[i].id();
			b_p_chg[i_np] = event[i].isCharged();
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;
			f_p_e[i_np] = event[i].e();
			f_p_m[i_np] = event[i].m();

			i_np++;

		}//i

		T->Fill();

  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;
	*/
}

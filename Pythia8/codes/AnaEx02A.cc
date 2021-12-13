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
#include "TChain.h"
#include "TLorentzVector.h"

#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main(int argc, char *argv[]) {
//int main() {

	if ( argc<3 ){
		cout << "Usage: ./AnaEx02A index cut_jetpt" << endl;
		return -1;
	}

	int index = atoi(argv[1]);
	float cut_jetpt = atof(argv[2]);
	cout << "Index: " << index << endl;
	cout << "Cut Jet pT: " << cut_jetpt << endl;

	//Input
	int i_np;
	int i_njetR04;

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

	int ii_np;

	int ii_p_id[5000];
	bool bb_p_chg[5000];
	float ff_p_pt[5000];
	float ff_p_eta[5000];
	float ff_p_phi[5000];
	float ff_p_e[5000];
	float ff_p_m[5000];

	TFile *infileSig = new TFile(Form("/alice/data/shlim/PYTHIA/pp5TeV_set20_grp000/outfile_pp5TeV_set20_grp000_%05d.root",index),"read");

	TTree *TSig = (TTree*)infileSig->Get("T");

	TSig->SetBranchAddress("np",&i_np);
	TSig->SetBranchAddress("p_id",i_p_id);
	TSig->SetBranchAddress("p_chg",b_p_chg);
	TSig->SetBranchAddress("p_pt",f_p_pt);
	TSig->SetBranchAddress("p_eta",f_p_eta);
	TSig->SetBranchAddress("p_phi",f_p_phi);
	TSig->SetBranchAddress("p_e",f_p_e);
	TSig->SetBranchAddress("p_m",f_p_m);

	TSig->SetBranchAddress("njetR04",&i_njetR04);
	TSig->SetBranchAddress("jetR04_pt",f_jetR04_pt);
	TSig->SetBranchAddress("jetR04_eta",f_jetR04_eta);
	TSig->SetBranchAddress("jetR04_phi",f_jetR04_phi);

	TChain *TBkg = new TChain("T");
	for (int ii=0; ii<5; ii++){
		TBkg->AddFile(Form("/alice/data/shlim/PYTHIA/pp5TeV_set21_grp000/outfile_pp5TeV_set21_grp000_%05d.root",5*index+ii));
	}

	TBkg->SetBranchAddress("np",&ii_np);
	TBkg->SetBranchAddress("p_id",ii_p_id);
	TBkg->SetBranchAddress("p_chg",bb_p_chg);
	TBkg->SetBranchAddress("p_pt",ff_p_pt);
	TBkg->SetBranchAddress("p_eta",ff_p_eta);
	TBkg->SetBranchAddress("p_phi",ff_p_phi);
	TBkg->SetBranchAddress("p_e",ff_p_e);
	TBkg->SetBranchAddress("p_m",ff_p_m);

	int nentriesSig = TSig->GetEntries();
	int nentriesBkg = TBkg->GetEntries();

	cout << "Number of events Sig: " << nentriesSig << " Bkg: " << nentriesBkg << endl;

	//Output
	int _i_nSigjetR04;
	int _i_nIncjetR04;

	float _f_SigjetR04_pt[10];
	float _f_SigjetR04_eta[10];
	float _f_SigjetR04_phi[10];

	float _f_IncjetR04_pt[10];
	float _f_IncjetR04_eta[10];
	float _f_IncjetR04_phi[10];

	// Tree output
	TFile *outfile = new TFile("Pythia8_JetReco.root","recreate");
	TTree *T = new TTree("T","Reco Jets");

	T->Branch("nSigjetR04",&_i_nSigjetR04,"nSigjetR04/I");
	T->Branch("SigjetR04_pt",_f_SigjetR04_pt,"SigjetR04_pt[nSigjetR04]/F");
	T->Branch("SigjetR04_eta",_f_SigjetR04_eta,"SigjetR04_eta[nSigjetR04]/F");
	T->Branch("SigjetR04_phi",_f_SigjetR04_phi,"SigjetR04_phi[nSigjetR04]/F");

	T->Branch("nIncjetR04",&_i_nIncjetR04,"nIncjetR04/I");
	T->Branch("IncjetR04_pt",_f_IncjetR04_pt,"IncjetR04_pt[nIncjetR04]/F");
	T->Branch("IncjetR04_eta",_f_IncjetR04_eta,"IncjetR04_eta[nIncjetR04]/F");
	T->Branch("IncjetR04_phi",_f_IncjetR04_phi,"IncjetR04_phi[nIncjetR04]/F");


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
	std::vector <fastjet::PseudoJet> fjInputsSig;
	std::vector <fastjet::PseudoJet> fjInputsInc;

	for (int ien=0; ien<nentriesSig; ien++){
	//for (int ien=0; ien<1000; ien++){

		TSig->GetEntry(ien);
		TBkg->GetEntry(ien);

		fjInputsSig.resize(0);
		fjInputsInc.resize(0);

		_i_nSigjetR04 = _i_nIncjetR04 = 0;
		for (int ii=0; ii<10; ii++){
			_f_SigjetR04_pt[ii] = _f_SigjetR04_eta[ii] = _f_SigjetR04_phi[ii] = -999;
			_f_IncjetR04_pt[ii] = _f_IncjetR04_eta[ii] = _f_IncjetR04_phi[ii] = -999;
		}

		//Signal
		for (int ip=0; ip<i_np; ip++){
			//ALICE acceptance
			if ( fabs(f_p_eta[ip])>1.0 ) continue;

			int id = abs(i_p_id[ip]);
			if ( id==12 || id==14 || id==16 ) continue;

			TLorentzVector lvec;
			lvec.SetPtEtaPhiE(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], f_p_e[ip]);

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.E());

			fjInputsInc.push_back(particleTemp);
			fjInputsSig.push_back(particleTemp);
		}//ip

		//Background
		for (int ip=0; ip<ii_np; ip++){
			//ALICE acceptance
			if ( fabs(ff_p_eta[ip])>1.0 ) continue;

			int id = abs(ii_p_id[ip]);
			if ( id==12 || id==14 || id==16 ) continue;

			TLorentzVector lvec;
			lvec.SetPtEtaPhiE(ff_p_pt[ip], ff_p_eta[ip], ff_p_phi[ip], ff_p_e[ip]);

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.E());

			fjInputsInc.push_back(particleTemp);
		}//ip

		//cout << "Sig: " << fjInputsSig.size() << ", Inc: " << fjInputsInc.size() << endl;

		//Jet reconstruction
		vector <fastjet::PseudoJet> IncJetsR04, SortedIncJetsR04;
		vector <fastjet::PseudoJet> SigJetsR04, SortedSigJetsR04;

		fastjet::ClusterSequence clustSeqSigR04(fjInputsSig, jetDefR04);
		SigJetsR04 = clustSeqSigR04.inclusive_jets(cut_jetpt);
		SortedSigJetsR04 = sorted_by_pt(SigJetsR04);

		fastjet::ClusterSequence clustSeqIncR04(fjInputsInc, jetDefR04);
		IncJetsR04 = clustSeqIncR04.inclusive_jets(cut_jetpt);
		SortedIncJetsR04 = sorted_by_pt(IncJetsR04);

		if ( (i_njetR04!=int(SortedSigJetsR04.size())) ){
			cout << "Number of full jets: " << i_njetR04 << ", " << SortedSigJetsR04.size() << endl;
			continue;
		}

		_i_nSigjetR04 = int(SortedSigJetsR04.size());
		for (int ii=0; ii<int(SortedSigJetsR04.size()); ++ii) {
			_f_SigjetR04_pt[ii] = SortedSigJetsR04[ii].perp(); 
			_f_SigjetR04_eta[ii] = SortedSigJetsR04[ii].rap();
			_f_SigjetR04_phi[ii] = SortedSigJetsR04[ii].phi_std();
		}

		_i_nIncjetR04 = int(SortedIncJetsR04.size());
		for (int ii=0; ii<int(SortedIncJetsR04.size()); ++ii) {
			_f_IncjetR04_pt[ii] = SortedIncJetsR04[ii].perp(); 
			_f_IncjetR04_eta[ii] = SortedIncJetsR04[ii].rap();
			_f_IncjetR04_phi[ii] = SortedIncJetsR04[ii].phi_std();
		}

		T->Fill();

	}//ien

	outfile->cd();
	T->Write();
	outfile->Close();


	//Done
	return 0;

}

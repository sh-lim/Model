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
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main(int argc, char *argv[]) {
//int main() {

	if ( argc<3 ){
		cout << "Usage: ./AnaEx02A filename cut_jetpt" << endl;
		return -1;
	}

	string filename = argv[1];
	float cut_jetpt = atof(argv[2]);
	cout << "Input file: " << filename << endl;
	cout << "Cut Jet pT: " << cut_jetpt << endl;

	TDatabasePDG *pdgdb = new TDatabasePDG();

	//fastjet
  // Set up FastJet jet finder.
  //   one can use either explicitly use antikt, cambridge, etc., or
  //   just use genkt_algorithm with specification of power
  //fastjet::JetAlgorithm algorithm;
  //if (power == -1)      algorithm = fastjet::antikt_algorithm;
  //if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  //if (power ==  1)      algorithm = fastjet::kt_algorithm;
  //fastjet::JetDefinition jetDef(algorithm, R);
  // there's no need for a pointer to the jetDef (it's a fairly small object)
  fastjet::JetDefinition jetDefR04(fastjet::genkt_algorithm, 0.4, -1);
	std::vector <fastjet::PseudoJet> fjInputsInc;

	//input
	TFile *infile = new TFile(filename.c_str(),"read");
	TTree *Tin = (TTree*)infile->Get("Events");

	int Event_numberP;
	int Particle_pid[5000];
	int Particle_status[5000];
	double Particle_ctau[5000];
	double Particle_px[5000];
	double Particle_py[5000];
	double Particle_pz[5000];
	double Particle_energy[5000];
	double Particle_mass[5000];

	Tin->SetBranchAddress("Event_numberP",&Event_numberP);
	Tin->SetBranchAddress("Particle_pid",Particle_pid);
	Tin->SetBranchAddress("Particle_status",Particle_status);
	Tin->SetBranchAddress("Particle_ctau",Particle_ctau);
	Tin->SetBranchAddress("Particle_px",Particle_px);
	Tin->SetBranchAddress("Particle_py",Particle_py);
	Tin->SetBranchAddress("Particle_pz",Particle_pz);
	Tin->SetBranchAddress("Particle_energy",Particle_energy);
	Tin->SetBranchAddress("Particle_mass",Particle_mass);
	//Tin->SetBranchAddress("",);

	int nentries = Tin->GetEntries();

	//output
	int i_njetR04;
	float f_jetR04_pt[10];
	float f_jetR04_eta[10];
	float f_jetR04_phi[10];

	int i_np;
	int i_p_id[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];
	float f_p_e[5000];
	float f_p_m[5000];


	TFile *outfile = new TFile("outfile_herwig.root","recreate");

	TTree *Tout = new TTree("T", "Reco Jets");
	Tout->Branch("njetR04",&i_njetR04,"njetR04/I");
	Tout->Branch("jetR04_pt",f_jetR04_pt,"jetR04_pt[njetR04]/F");
	Tout->Branch("jetR04_eta",f_jetR04_eta,"jetR04_eta[njetR04]/F");
	Tout->Branch("jetR04_phi",f_jetR04_phi,"jetR04_phi[njetR04]/F");

	Tout->Branch("np",&i_np,"np/I");
	Tout->Branch("p_id",i_p_id,"p_id[np]/I");
	Tout->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	Tout->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	Tout->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	Tout->Branch("p_e",f_p_e,"p_e[np]/F");
	Tout->Branch("p_m",f_p_m,"p_m[np]/F");

	for (int ien=0; ien<nentries; ien++){
	//for (int ien=0; ien<10; ien++){

		Tin->GetEntry(ien);

		for (int ip=0; ip<Event_numberP; ip++){

			if ( Particle_status[ip]!=1 ) continue;
			if ( Particle_ctau[ip]>10.0 ) continue;

			Float_t charge = pdgdb->GetParticle(Particle_pid[ip])->Charge();
			if ( fabs(charge)<0.5 ) continue;

			TLorentzVector lvec;
			lvec.SetPxPyPzE(Particle_px[ip], Particle_py[ip], Particle_pz[ip], Particle_energy[ip]);

			if ( fabs(lvec.Eta())>1.0 ) continue;

			// Create a PseudoJet from the complete Pythia particle.
			fastjet::PseudoJet particleTemp(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.E());

			fjInputsInc.push_back(particleTemp);

			i_p_id[i_np] = Particle_pid[ip];
			f_p_pt[i_np] = lvec.Pt();
			f_p_eta[i_np] = lvec.Eta();
			f_p_phi[i_np] = lvec.Phi();
			f_p_e[i_np] = Particle_energy[ip];
			f_p_m[i_np] = Particle_mass[ip];

			i_np++;
		}//ip

		//jet reconstruction
		vector <fastjet::PseudoJet> IncJetsR04, SortedIncJetsR04;

		fastjet::ClusterSequence clustSeqIncR04(fjInputsInc, jetDefR04);
		IncJetsR04 = clustSeqIncR04.inclusive_jets(cut_jetpt);
		SortedIncJetsR04 = sorted_by_pt(IncJetsR04);

		i_njetR04 = int(SortedIncJetsR04.size());
		if ( i_njetR04<1 ) continue;

		for (int ii=0; ii<int(SortedIncJetsR04.size()); ++ii) {
			f_jetR04_pt[ii] = SortedIncJetsR04[ii].perp(); 
			f_jetR04_eta[ii] = SortedIncJetsR04[ii].rap();
			f_jetR04_phi[ii] = SortedIncJetsR04[ii].phi_std();
		}

		Tout->Fill();

		//clean up
		i_np = 0, i_njetR04 = 0;
		for (int ii=0; ii<10; ii++){
			f_jetR04_pt[ii] = f_jetR04_eta[ii] = f_jetR04_phi[ii] = -999;
		}
		fjInputsInc.clear();

	}//ien

	Tout->Write();
	outfile->Close();

	//Done
	return 0;

}

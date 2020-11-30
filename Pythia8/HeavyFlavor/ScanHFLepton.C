void ScanHFLepton(){

	char fname[300];

	int np;
	int p_id[50];
	float p_pt[50], p_eta[50], p_phi[50];
	bool p_cc[50], p_bb[50];

	TH1D *hHFemu_cc = new TH1D("hHFemu_cc","",100,0,10);
	TH1D *hHFemu_bb = new TH1D("hHFemu_bb","",100,0,10);

	ifstream flist;
	flist.open("file.lst");

	while ( flist >> fname ){

		cout << "OPEN: " << fname << endl;

		TFile *infile = new TFile(fname,"read");
		TTree *T = (TTree*)infile->Get("T");

		T->SetBranchAddress("np",&np);
		T->SetBranchAddress("p_id",p_id);
		T->SetBranchAddress("p_pt",p_pt);
		T->SetBranchAddress("p_eta",p_eta);
		T->SetBranchAddress("p_phi",p_phi);
		T->SetBranchAddress("p_cc",p_cc);
		T->SetBranchAddress("p_bb",p_bb);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			for (int ip=0; ip<np; ip++){

				if ( fabs(p_eta[ip])>0.35 ) continue;

				if ( p_bb[ip] ){
					hHFemu_bb->Fill(p_pt[ip]);
				}else if ( p_cc[ip] ){
					hHFemu_cc->Fill(p_pt[ip]);
				}

			}

		}//ien


		infile->Close();
		delete infile;
	}//


	TFile *outfile = new TFile("outfile.root","recreate");
	hHFemu_bb->Write();
	hHFemu_cc->Write();
	outfile->Close();


}

void merge_hist(const char *package="pp13TeV_inelastic"){

	ifstream flist;
	flist.open(Form("file_%s.lst",package));

	char fname[500];
	int bfirst = true;

	TH1 *h[2000];
	TH1 *htmp[2000];

	TFile *infile0;

	while ( flist >> fname ){
		TString fname_str(fname);
		cout << "OPEN: " << fname_str.Data() << endl;

		if ( bfirst ){
			infile0 = new TFile(fname_str,"read");
			TList *key_list = infile0->GetListOfKeys(); 
			cout << "size: " << key_list->GetSize() << endl;
			for (int ij=0; ij<key_list->GetSize(); ij++){
				TKey *key = (TKey*)key_list->At(ij);
				h[ij] = (TH1*)key->ReadObj();
				//cout << key->GetName() << endl;
				cout << h[ij]->GetName() << endl;
			}   
			bfirst = false;
		}else{
			TFile *infile1 = new TFile(fname_str,"read");
			TList *key_list = infile1->GetListOfKeys(); 

			for (int ij=0; ij<key_list->GetSize(); ij++){
				TKey *key = (TKey*)key_list->At(ij);
				htmp[ij] = (TH1*)key->ReadObj();
				h[ij]->Add(htmp[ij]);
			}   

			infile1->Close();
			delete infile1;
		}   
	}

	TFile *outfile = new TFile(Form("outfile_hist_%s.root",package),"recreate");
	TList *key_list = infile0->GetListOfKeys(); 
	for (int ij=0; ij<key_list->GetSize(); ij++){
		h[ij]->Write();
	}

}


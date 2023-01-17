const Int_t a = 269;
Int_t eposid[a] = {1, 2, 3, 4, 5, 6, 10, 9, 12, -12, 11, -11, 14, -14, 13, -13, 16, 15, 110, 120, -120, 220, 130, -130, 230, -230, 20, -20, 330, 111, 121, -121, 221, 131, -131, 231, -231, 331, -140, 240, 1120, 1220, 2130, 1130, 1230, 2230, 1330, 2330, 1111, 1121, 1221, 2221, 1131, 2231, 1331, 2331, 3331, 2140, 17, 18, 19, 0, 99, 1112, 1113, 1114, 2222, 2223, 2224, 2224, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1132, 1133, 1134, 2232, 2233, 2234, 90, 80, 81, 85, 86, 82, 83, 84, 1200, 2300, 1300, 2400, 1400, 3400, 2500, 1500, 3500, 4500, 2200, 1200, 2300, 1300, 3300, 2400, 1400, 3400, 4400, 2500, 1500, 3500, 4500, 5500, 800000091, 800000092, 800000093, 800000094, -340, 340, -241, 241, -141, 141, -341, 341, 250, 150, 350, 450, 251, 151, 351, 451, 440, 441, 550, 551, 2240, 1240, 1140, 2241, 1241, 3240, 2340, 3140, 1340, 3340, 2341, 1341, 3341, 2440, 2441, 1440, 1441, 3440, 3441, 4441, 2250, 2150, 3250, 4250, 1250, 1150, 3150, 4150, 2350, 1350, 3350, 4350, 2450, 1450, 3450, 4450, 2550, 1550, 3550, 2251, 1151, 2351, 1351, 3351, 2451, 1451, 3451, 4451, 2551, 1551, 3551, 4551, 5551, 123, 122, 233, 232, 133, 132, 143, 132, 243, 242, 343, 342, 223, 222, 113, 112, 333, 332, 443, 442, 444, 253, 252, 153, 152, 353, 352, 453, 452, 553, 552, 124, 125, 234, 235, 134, 135, 144, 135, 244, 245, 344, 345, 224, 225, 114, 115, 334, 335, 444, 445, 254, 255, 154, 155, 354, 355, 454, 455, 554, 555, 11099, 12099, 22099, 33099, 44099, 112099, 122099, 800000110, 800000990};
Int_t pdgid[a] = {2, 1, 3, 4, 5, 6, 22, 21, 11, -11, 12, -12, 13, -13, 14, -14, 15, 16, 111, 211, -211, 221, 321, -321, 311, -311, 310, -310, 331, 113, 213, -213, 223, 323, -323, 313, -313, 333, 421, -411, 2212, 2112, 3122, 3222, 3212, 3112, 3322, 3312, 2224, 2214, 2114, 1114, 3224, 3114, 3324, 3314, 3334, 4122, 99, 99, 99, 99, 99, 32224, 12224, 12222, 31114, 11114, 11112, 21114, 12212, 2124, 32214, 2216, 12214, 22124, 11212, 13122, 3124, 23122, 13212, 23212, 53122, 13216, 13222, 23222, 13226, 13112, 23112, 13116, 23, 24, 25, 32, 33, 35, 36, 37, 2101, 3101, 3201, 4101, 4201, 4301, 5101, 5201, 5301, 5401, 1103, 1103, 3103, 3203, 3303, 4103, 4203, 4303, 4403, 5103, 5203, 5303, 5403, 5503, 91, 92, 93, 94, 431, -431, 413, -413, 423, -423, 433, -433, 511, 521, 531, 541, 513, 523, 533, 543, 441, 443, 551, 553, 4112, 4212, 4222, 4114, 4224, 4132, 4312, 4232, 4322, 4332, 4314, 4324, 4334, 4412, 4414, 4422, 4424, 4432, 4434, 4444, 5112, 5122, 5132, 5142, 5212, 5222, 5232, 5242, 5312, 5322, 5332, 5342, 5412, 5422, 5432, 5442, 5512, 5522, 5532, 5114, 5224, 5314, 5324, 5334, 5414, 5424, 5434, 5444, 5524, 5524, 5534, 5544, 5554, 10213, 10211, 10313, 10311, 10323, 10321, 10423, 10421, 10413, 10411, 10433, 10431, 10113, 10111, 10223, 10221, 10333, 10331, 10443, 10441, 10443, 10513, 10511, 10523, 10521, 10533, 10531, 10543, 10541, 10553, 10551, 20213, 215, 20313, 315, 20323, 325, 20423, 425, 20413, 415, 20433, 435, 20113, 115, 20223, 225, 20333, 335, 20443, 445, 20513, 515, 20523, 525, 20533, 535, 20543, 545, 20553, 555, 9900110, 9900210, 9900220, 9900330, 9900440, 9902210, 9902110, 110, 990}; 

Int_t EposToPdg(Int_t code)
{
	for(Int_t i=0;i<a;i++) {
		if (eposid[i]==code) return pdgid[i];
	}

	return 0;	
}


void ScanEPOS02(){

	//output file

	TFile *outfile = new TFile("outfile_epos4.root","recreate");

	int _np = 0;
	int _pid[5000];
	float _pt[5000];
	float _eta[5000];
	float _phi[5000];

	TTree *T = new TTree("T","T");
	T->Branch("np",&_np,"np/I");
	T->Branch("pid",_pid,"pid[np]/I");
	T->Branch("pt",_pt,"pt[np]/F");
	T->Branch("eta",_eta,"eta[np]/F");
	T->Branch("phi",_phi,"phi[np]/F");

	// Parameters in EPOS Tree:
	//--------------teposevent---------------
	Int_t np;       //number of particles
	Float_t bim;    //impact parameter
	std::vector<Float_t> zus; //different meaning depending on ptl type:
	// partons: presently unused
	// hadrons:  decay information :
	// -999 : hadron is decay product from decay 
	//        in cascade part (mother unknown)
	//   -1 : hadron is decay product, mother not stored
	//   >0 : hadron is decay product, mother index = zus
	//   so, zus is the number n of the particle, not the id
	//   to get the id you need id(n)
	//   -2 : no mother
	// phi 331, K+ 130
	std::vector<Float_t> px;  //px
	std::vector<Float_t> py;  //py   particle four momentum
	std::vector<Float_t> pz;  //pz
	std::vector<Float_t> e;   //energy pf particle
	std::vector<Float_t> x;   //x component of formation point
	std::vector<Float_t> y;   //y component of formation point
	std::vector<Float_t> z;   //z component of formation point
	std::vector<Float_t> t;   //formation time
	std::vector<Int_t> id;    //particle id
	std::vector<Int_t> ist;   //particle status (hadron last generation(0) or not(1))
	std::vector<Int_t> ity;   //type of particle origin (20-29 from soft strings, 30-39 from hard strings, 40-59 from remnants, 60 from fluid)
	std::vector<Int_t> ior;   //index of father (resonance decay products)
	std::vector<Int_t> jor;   //index of mother (mothers are needed for exemple for strings: the partons between ior and jor constitute the string)

	//  --------------teposhead---------------
	Int_t fIversn;  //EPOS version number
	Int_t fLaproj;  //atomic number projectile
	Int_t fMaproj;  //mass number projectile
	Int_t fLatarg;  //atomic number target
	Int_t fMatarg;  //mass number target
	Float_t fEngy;  //energy in the CMS in GeV
	Int_t fNfull;   //number of full events
	Int_t fNfreeze; //number of freeze outs per full event

	ifstream flist;
	flist.open("file.lst");

	string fname;
	flist >> fname;

	flist.close();

	TFile *infile = new TFile(fname.c_str(),"read");

	if ( infile->IsOpen() ){
		cout << "OPEN: " << infile->GetName() << endl;
	}else{
		return;
	}

	TTree *fTreeHeader = (TTree*)gDirectory->Get("teposhead");
	fTreeHeader->SetMakeClass(1);
	fTreeHeader->SetBranchAddress("iversn",&fIversn);
	fTreeHeader->SetBranchAddress("laproj",&fLaproj);
	fTreeHeader->SetBranchAddress("maproj",&fMaproj);
	fTreeHeader->SetBranchAddress("latarg",&fLatarg);
	fTreeHeader->SetBranchAddress("matarg",&fMatarg);
	fTreeHeader->SetBranchAddress("engy",&fEngy);
	fTreeHeader->SetBranchAddress("nfull",&fNfull);
	fTreeHeader->SetBranchAddress("nfreeze",&fNfull);
	fTreeHeader->GetEvent(0);

	TTree *fTreeNtuple = (TTree*)gDirectory->Get("teposevent");
	fTreeNtuple->SetMakeClass(1);

	const Int_t maxnp = fTreeNtuple->GetMaximum("np");
	cout << "Maximum np: " << maxnp << endl;

	zus.resize(maxnp);
	px.resize(maxnp);
	py.resize(maxnp);
	pz.resize(maxnp);
	e.resize(maxnp);
	x.resize(maxnp);
	y.resize(maxnp);
	z.resize(maxnp);
	t.resize(maxnp);
	id.resize(maxnp);
	ist.resize(maxnp);
	ity.resize(maxnp);
	ior.resize(maxnp);
	jor.resize(maxnp);

	fTreeNtuple->SetBranchAddress("np",&np);
	fTreeNtuple->SetBranchAddress("bim",&bim);
	fTreeNtuple->SetBranchAddress("zus",&zus[0]);
	fTreeNtuple->SetBranchAddress("px",&px[0]);
	fTreeNtuple->SetBranchAddress("py",&py[0]);
	fTreeNtuple->SetBranchAddress("pz",&pz[0]);
	fTreeNtuple->SetBranchAddress("e",&e[0]);
	fTreeNtuple->SetBranchAddress("x",&x[0]);
	fTreeNtuple->SetBranchAddress("y",&y[0]);
	fTreeNtuple->SetBranchAddress("z",&z[0]);
	fTreeNtuple->SetBranchAddress("t",&t[0]);
	fTreeNtuple->SetBranchAddress("id",&id[0]);
	fTreeNtuple->SetBranchAddress("ist",&ist[0]);
	fTreeNtuple->SetBranchAddress("ity",&ity[0]);
	fTreeNtuple->SetBranchAddress("ior",&ior[0]);
	fTreeNtuple->SetBranchAddress("jor",&jor[0]);

	int nentries = fTreeNtuple->GetEntries();

	for (int ien=0; ien<nentries; ien++){

		fTreeNtuple->GetEntry(ien);

		_np = 0;
		for (int ii=0; ii<5000; ii++){
			_pid[ii] = 0;
			_pt[ii] = _eta[ii] = _phi[ii] = 0;
		}

		for (int ip=0; ip<np; ip++){

			if ( ist[ip]!=0 ) continue;
			//if ( zus[ip]<0 ) continue;

			int _id = EposToPdg(id[ip]);

			if ( !(abs(_id)==211 || abs(_id)==321 || abs(_id)==2212) ) continue;

			TVector3 vec(px[ip], py[ip], pz[ip]);

			_pid[_np] = _id;
			_pt[_np] = vec.Pt();
			_eta[_np] = vec.Eta();
			_phi[_np] = vec.Phi();

			++_np;

		}//ip

		T->Fill();

	}//ien

	outfile->cd();
	T->Write();

	outfile->Close();

	return;

	/*
	*/

}

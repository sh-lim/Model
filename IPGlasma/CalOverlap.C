void CalOverlap(const char *set="pPb", const int grp=13){

	const int nevt = 2000;

	TH1D *h1 = new TH1D("h1","",100000,0,1000000);

	TH1D *he2 = new TH1D("he2","",100,0,1);
	TH1D *he3 = new TH1D("he3","",100,0,1);

	for (int ievt=0; ievt<nevt; ievt++){

		TFile *infile = new TFile(Form("IPGlasma_%s_grp%d/Epsilon_IPGlasma_%s_run%05d.root",set,grp,set,ievt),"read");
		cout << "OPEN: " << infile->GetName() << endl;
		TH1D *_h1t = (TH1D*)infile->Get("h1_time");
		int ntime = int(_h1t->GetEntries());
		TH2D *_h2 = (TH2D*)infile->Get(Form("h2_evt00000_t0000%d",ntime-1));

		h1->Fill(_h2->Integral());

		infile->Close();
		delete infile;
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*400,400);
	gPad->SetLogy();
	h1->Draw();

	cout << h1->Integral() << endl;

	int ncent0005 = 0;

	for (int ievt=0; ievt<nevt; ievt++){

		TFile *infile = new TFile(Form("IPGlasma_%s_grp%d/Epsilon_IPGlasma_%s_run%05d.root",set,grp,set,ievt),"read");
		cout << "OPEN: " << infile->GetName() << endl;
		TH1D *_h1t = (TH1D*)infile->Get("h1_time");
		int ntime = int(_h1t->GetEntries());
		TH2D *_h2 = (TH2D*)infile->Get(Form("h2_evt00000_t0000%d",ntime-1));

		float norm = _h2->Integral();
		int bin = h1->FindBin(norm);
		float frac = h1->Integral(1,bin-1)/nevt;

		if ( frac<0.95 ){
			infile->Close();
			delete infile;
			continue;
		}

		float meanXw = 0, meanYw = 0;
		float sumw = 0;

		float nbinsx = _h2->GetNbinsX();
		float nbinsy = _h2->GetNbinsY();

		for (int ix=0; ix<nbinsx; ix++){
			for (int iy=0; iy<nbinsy; iy++){
				float xx = _h2->GetXaxis()->GetBinCenter(ix+1);
				float yy = _h2->GetYaxis()->GetBinCenter(iy+1);

				float ww = _h2->GetBinContent(ix+1, iy+1);

				meanXw += xx*ww;
				meanYw += yy*ww;

				sumw += ww;
			}//iy
		}//ix

		meanXw /= sumw;
		meanYw /= sumw;

		double cosphi2 = 0.0, sinphi2 = 0.0, rn2 = 0.0;;
		double cosphi3 = 0.0, sinphi3 = 0.0, rn3 = 0.0;;

		for (int ix=0; ix<nbinsx; ix++){
			for (int iy=0; iy<nbinsy; iy++){
				float xx = _h2->GetXaxis()->GetBinCenter(ix+1) - meanXw;
				float yy = _h2->GetYaxis()->GetBinCenter(iy+1) - meanYw;

				float rr = sqrt(xx*xx + yy*yy);
				float phi = atan2(yy, xx);

				float ww = _h2->GetBinContent(ix+1, iy+1);

				cosphi2 += ww*pow(rr,2)*cos(2*phi);
				sinphi2 += ww*pow(rr,2)*sin(2*phi);
				rn2 += ww*pow(rr,2);

				cosphi3 += ww*pow(rr,3)*cos(3*phi);
				sinphi3 += ww*pow(rr,3)*sin(3*phi);
				rn3 += ww*pow(rr,3);
			}
		}

		float e2 = sqrt(cosphi2*cosphi2 + sinphi2*sinphi2)/rn2;
		float e3 = sqrt(cosphi3*cosphi3 + sinphi3*sinphi3)/rn3;

		he2->Fill(e2);
		he3->Fill(e3);

		ncent0005++;

		infile->Close();
		delete infile;

	}//ievt

	cout << "ncent0005: " << ncent0005 << endl;

	TCanvas *c2 = new TCanvas("c2","c2",1.1*2*400,400);
	c2->Divide(2,1);

	c2->cd(1);
	he2->Draw();

	c2->cd(2);
	he3->Draw();

	TFile *outfile = new TFile(Form("outfile_Ecc_IPGlasma_%s_grp%d.root",set,grp),"recreate");

	he2->Write();
	he3->Write();

	outfile->Close();

}

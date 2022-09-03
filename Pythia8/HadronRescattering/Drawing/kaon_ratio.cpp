void kaon_ratio(){

	gStyle -> SetOptStat(0);

	const int had = 2;
	const int nfile = 2;
	const int nmult = 6;
	const int npt = 15;

	TFile * infile[had][nfile];

	infile[0][0] = new TFile("outfile_kaon_on_0_changeptmult.root", "read");
	infile[0][1] = new TFile("outfile_kaon_on_1_changeptmult.root", "read");
	infile[1][0] = new TFile("outfile_kaon_off_0_changeptmult.root", "read");
	infile[1][1] = new TFile("outfile_kaon_off_1_changeptmult.root", "read");

	TH1D * ratio[had][nfile][nmult];

	for( int i = 0 ; i < nmult ; i++ ){
		for( int j = 0 ; j < had ; j++ ){
			for( int k = 0 ; k < nfile ; k++ ){

				ratio[j][k][i] = ( TH1D * ) infile[j][k] -> Get(Form("histo_kaon_pT_m%d", i));
			}//k
			ratio[j][0][i] -> Add(ratio[j][1][i], 1);
		}//j
		ratio[0][0][i] -> Divide(ratio[1][0][i]);
		ratio[0][0][i] -> Sumw2();
	}//i

	TCanvas * c0 = new TCanvas("c0", "c0", 1200, 1200);

	ratio[0][0][0] -> SetMarkerStyle(20);
	ratio[0][0][0] -> SetMarkerColor(kRed);
	ratio[0][0][1] -> SetMarkerStyle(21);
	ratio[0][0][1] -> SetMarkerColor(kCyan);
	ratio[0][0][2] -> SetMarkerStyle(22);
	ratio[0][0][2] -> SetMarkerColor(kGreen + 1);
	ratio[0][0][3] -> SetMarkerStyle(23);
	ratio[0][0][3] -> SetMarkerColor(kBlue);
	ratio[0][0][4] -> SetMarkerStyle(34);
	ratio[0][0][4] -> SetMarkerColor(kViolet);
	ratio[0][0][5] -> SetMarkerStyle(8);
	ratio[0][0][5] -> SetMarkerColor(kBlack);

	ratio[0][0][0] -> GetXaxis() -> SetTitle("pT (GeV/c)");

	ratio[0][0][0] -> GetYaxis() -> SetRangeUser(0.90, 1.10);

	for( int i = 0 ; i < nmult ; i++ ){
		ratio[0][0][i] -> Draw("HIST SAME P");
	}//i

	TLegend * rati = new TLegend(0.5, 0.7, 0.8, 0.85);
	rati -> SetFillStyle(0);
	rati -> SetBorderSize(0);
	rati -> SetTextFont(43);
	rati -> SetTextSize(21);
	rati -> AddEntry(ratio[0][0][5], " I ", "pl");
	rati -> AddEntry(ratio[0][0][4], " II ", "pl");
	rati -> AddEntry(ratio[0][0][3], " III ", "pl");
	rati -> AddEntry(ratio[0][0][2], " IV ", "pl");
	rati -> AddEntry(ratio[0][0][1], " V ", "pl");
	rati -> AddEntry(ratio[0][0][0], " VI ", "pl");
	rati -> Draw("same");
}

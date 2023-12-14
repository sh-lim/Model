void Draw00(){

	bool bSAVE = false;

	const int nset = 2;

	int nMarker[4] = {20, 24, 21, 25};
	int nColor[4] = {1, 2, 4, 6};

	TFile *infile[nset];
	infile[0] = new TFile("./Ana00_pp13TeV_opt1_grp140.root","read");
	infile[1] = new TFile("./Ana00_pp13TeV_opt2_grp140.root","read");
	//infile[0] = new TFile("Ana00_pp13TeV_opt3_grp110.root","read");
	//infile[3] = new TFile("Ana00_pp13TeV_opt4_grp110.root","read");

	TH1D *hmult[nset][2];
	TH1D *heta[nset][2];
	TH1D *hmult_211[nset][2];
	TH1D *hmult_321[nset][2];
	TH1D *hmult_113[nset][6];
	TH1D *hmult_313[nset][6];
	TH1D *hmult_333[nset][6];

	TH1D *hratio_333_321[nset][6];
	TH1D *hratio_313_321[nset][6];
	TH1D *hratio_113_211[nset][6];
	TH1D *hratio_321_211[nset][2];

	for (int iset=0; iset<nset; iset++){

		for (int ii=0; ii<2; ii++){

			hmult[iset][ii] = (TH1D*)infile[iset]->Get(Form("hmult_opt%d",ii));
			hmult[iset][ii]->Sumw2();
			hmult[iset][ii]->SetMarkerStyle(nMarker[ii]);
			hmult[iset][ii]->SetMarkerSize(0.6);
			hmult[iset][ii]->SetMarkerColor(nColor[iset]);
			hmult[iset][ii]->SetLineColor(nColor[iset]);
			//hmult[iset][ii]->SetLineStyle(ii+1);

			cout << iset << " " << ii << " " << hmult[iset][ii]->Integral() << " " << hmult[iset][ii]->Integral(2,-1) << endl;

			heta[iset][ii] = (TH1D*)infile[iset]->Get(Form("heta_opt%d",ii));
			heta[iset][ii]->Scale(1./(hmult[iset][ii]->Integral(2,-1)*heta[iset][ii]->GetBinWidth(1)));
			heta[iset][ii]->SetLineColor(nColor[iset]);
			heta[iset][ii]->SetMarkerColor(nColor[iset]);
			heta[iset][ii]->SetMarkerStyle(nMarker[ii]);
			heta[iset][ii]->SetMarkerSize(0.6);

			hmult_211[iset][ii] = (TH1D*)infile[iset]->Get(Form("hmult_211_opt%d",ii));
			hmult_321[iset][ii] = (TH1D*)infile[iset]->Get(Form("hmult_321_opt%d",ii));

			hmult_211[iset][ii]->Rebin(10);
			hmult_321[iset][ii]->Rebin(10);

			hmult_211[iset][ii]->Sumw2();
			hmult_321[iset][ii]->Sumw2();

			hratio_321_211[iset][ii] = (TH1D*)hmult_321[iset][ii]->Clone(Form("hratio_321_211_%d_%d",iset,ii));
			hratio_321_211[iset][ii]->Divide(hmult_211[iset][ii]);
			hratio_321_211[iset][ii]->SetMarkerStyle(nMarker[ii]);
			hratio_321_211[iset][ii]->SetMarkerColor(nColor[iset]);
			hratio_321_211[iset][ii]->SetLineColor(nColor[iset]);
			hratio_321_211[iset][ii]->SetMarkerSize(0.8);
		}//ii

		for (int ii=0; ii<6; ii++){
			hmult_113[iset][ii] = (TH1D*)infile[iset]->Get(Form("hmult_113_opt%d",ii));
			cout << "113 " << iset << " " << ii << " " << hmult_113[iset][ii]->Integral() << endl;

			hmult_113[iset][ii]->Rebin(10);
			hmult_113[iset][ii]->Sumw2();

			if ( ii==3 ){
				hmult_113[iset][ii]->Add(hmult_113[iset][ii-1]);
			}

			hratio_113_211[iset][ii] = (TH1D*)hmult_113[iset][ii]->Clone(Form("hratio_113_211_%d_%d",iset,ii));
			if ( iset==0 ){
				hratio_113_211[iset][ii]->Divide(hmult_211[iset][1]);
			}else{
				hratio_113_211[iset][ii]->Divide(hmult_211[iset][0]);
			}
			hratio_113_211[iset][ii]->SetMarkerStyle(nMarker[ii]);
			hratio_113_211[iset][ii]->SetMarkerColor(nColor[iset]);
			hratio_113_211[iset][ii]->SetLineColor(nColor[iset]);
			hratio_113_211[iset][ii]->SetMarkerSize(0.8);
			//hratio_113_211[iset][ii]->Scale(2.0);
		}

		for (int ii=0; ii<6; ii++){
			hmult_313[iset][ii] = (TH1D*)infile[iset]->Get(Form("hmult_313_opt%d",ii));
			cout << "313 " << iset << " " << ii << " " << hmult_313[iset][ii]->Integral() << endl;

			hmult_313[iset][ii]->Rebin(10);
			hmult_313[iset][ii]->Sumw2();

			if ( ii==3 ){
				hmult_313[iset][ii]->Add(hmult_313[iset][ii-1]);
			}

			hratio_313_321[iset][ii] = (TH1D*)hmult_313[iset][ii]->Clone(Form("hratio_313_321_%d_%d",iset,ii));
			if ( iset==0 ){
				hratio_313_321[iset][ii]->Divide(hmult_321[iset][1]);
			}else{
				hratio_313_321[iset][ii]->Divide(hmult_321[iset][0]);
			}
			hratio_313_321[iset][ii]->SetMarkerStyle(nMarker[ii]);
			hratio_313_321[iset][ii]->SetMarkerColor(nColor[iset]);
			hratio_313_321[iset][ii]->SetLineColor(nColor[iset]);
			hratio_313_321[iset][ii]->SetMarkerSize(0.8);
			hratio_313_321[iset][ii]->Scale(2.0/0.67);
		}//ii

		for (int ii=0; ii<6; ii++){
			hmult_333[iset][ii] = (TH1D*)infile[iset]->Get(Form("hmult_333_opt%d",ii));
			cout << "333 " << iset << " " << ii << " " << hmult_333[iset][ii]->Integral() << endl;

			hmult_333[iset][ii]->Rebin(10);
			hmult_333[iset][ii]->Sumw2();

			if ( ii==3 ){
				hmult_333[iset][ii]->Add(hmult_333[iset][ii-1]);
			}

			hratio_333_321[iset][ii] = (TH1D*)hmult_333[iset][ii]->Clone(Form("hratio_333_321_%d_%d",iset,ii));
			if ( iset==0 ){
				hratio_333_321[iset][ii]->Divide(hmult_321[iset][1]);
			}else{
				hratio_333_321[iset][ii]->Divide(hmult_321[iset][0]);
			}
			hratio_333_321[iset][ii]->SetMarkerStyle(nMarker[ii]);
			hratio_333_321[iset][ii]->SetMarkerColor(nColor[iset]);
			hratio_333_321[iset][ii]->SetLineColor(nColor[iset]);
			hratio_333_321[iset][ii]->SetMarkerSize(0.8);
			hratio_333_321[iset][ii]->Scale(2.0/0.49);
		}

		/*

		if ( iset<2 ){
			hmult_333[iset]->Scale(1./0.49);
			hmult_313[iset]->Scale(1./0.67);
		}
		hmult_211[iset]->Scale(0.5);
		hmult_321[iset]->Scale(0.5);
		hmult_313[iset]->Scale(0.5);

		*/

	}//iset

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	c1->SetMargin(0.14,0.03,0.12,0.03);
	gPad->SetTicks();

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1,120,2*hmult[0][0]->GetMaximum());
		gPad->SetLogy();

		hmult[0][0]->Draw("p same");
		hmult[0][1]->Draw("p same");
		hmult[1][0]->Draw("p same");
	}

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	c2->SetMargin(0.14,0.03,0.12,0.03);
	gPad->SetTicks();

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(-5,0,5,10);
		htmp->GetXaxis()->SetTitle("#eta");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("dN_{ch}/d#eta");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		heta[0][0]->Draw("p same");
		heta[0][1]->Draw("p same");
		heta[1][0]->Draw("p same");

		TLegend *leg = new TLegend(0.4,0.2,0.9,0.35);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(heta[0][0],Form("opt 0, before UrQMD"),"p");
		leg->AddEntry(heta[0][1],Form("opt 0, after UrQMD"),"p");
		leg->AddEntry(heta[1][0],Form("opt 1, no UrQMD"),"p");
		leg->Draw();

	}

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	c3->SetMargin(0.14,0.03,0.12,0.03);
	gPad->SetTicks();

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.2);
		htmp->GetXaxis()->SetTitle("N_{ch, |#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("K/#pi");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		hratio_321_211[0][0]->Draw("p same");
		hratio_321_211[0][1]->Draw("p same");
		hratio_321_211[1][0]->Draw("p same");

		TLegend *leg = new TLegend(0.4,0.2,0.9,0.35);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(heta[0][0],Form("opt 0, before UrQMD"),"p");
		leg->AddEntry(heta[0][1],Form("opt 0, after UrQMD"),"p");
		leg->AddEntry(heta[1][0],Form("opt 1, no UrQMD"),"p");
		leg->Draw();
	}

	TCanvas *c6 = new TCanvas("c6","c6",1.2*500,500);
	c6->SetMargin(0.14,0.03,0.12,0.03);
	gPad->SetTicks();

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.15);
		htmp->GetXaxis()->SetTitle("N_{ch, |#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("#rho^{0}/#pi^{#pm}");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		hratio_113_211[0][1]->Draw("p same");
		hratio_113_211[0][2]->Draw("p same");
		hratio_113_211[0][3]->Draw("p same");
		hratio_113_211[1][0]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.2,0.9,0.4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(hratio_113_211[0][1],Form("opt 0, after UrQMD"),"p");
		leg->AddEntry(hratio_113_211[0][3],Form("opt 0, EPOS3 tag, all resonance"),"p");
		leg->AddEntry(hratio_113_211[0][2],Form("opt 0, EPOS3 tag, survived resonance"),"p");
		leg->AddEntry(hratio_113_211[1][0],Form("opt 1, no UrQMD"),"p");
		leg->Draw();
	}

	TCanvas *c5 = new TCanvas("c5","c5",1.2*500,500);
	c5->SetMargin(0.14,0.03,0.12,0.03);
	gPad->SetTicks();

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,1.0);
		htmp->GetXaxis()->SetTitle("N_{ch, |#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("K^{*0}/K^{-}");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		hratio_313_321[0][1]->Draw("p same");
		hratio_313_321[0][2]->Draw("p same");
		hratio_313_321[0][3]->Draw("p same");
		hratio_313_321[1][0]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.2,0.9,0.4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(hratio_113_211[0][1],Form("opt 0, after UrQMD"),"p");
		leg->AddEntry(hratio_113_211[0][3],Form("opt 0, EPOS3 tag, all resonance"),"p");
		leg->AddEntry(hratio_113_211[0][2],Form("opt 0, EPOS3 tag, survived resonance"),"p");
		leg->AddEntry(hratio_113_211[1][0],Form("opt 1, no UrQMD"),"p");
		leg->Draw();
	}

	TCanvas *c4 = new TCanvas("c4","c4",1.2*500,500);
	c4->SetMargin(0.14,0.03,0.12,0.03);
	gPad->SetTicks();

	{
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,120,0.20);
		htmp->GetXaxis()->SetTitle("N_{ch, |#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitle("#phi/K^{-}");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);

		hratio_333_321[0][1]->Draw("p same");
		hratio_333_321[0][2]->Draw("p same");
		hratio_333_321[0][3]->Draw("p same");
		hratio_333_321[1][0]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.2,0.9,0.4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(hratio_113_211[0][1],Form("opt 0, after UrQMD"),"p");
		leg->AddEntry(hratio_113_211[0][3],Form("opt 0, EPOS3 tag, all resonance"),"p");
		leg->AddEntry(hratio_113_211[0][2],Form("opt 0, EPOS3 tag, survived resonance"),"p");
		leg->AddEntry(hratio_113_211[1][0],Form("opt 1, no UrQMD"),"p");
		leg->Draw();
	}

	if ( bSAVE ){
		c2->cd();
		c2->SaveAs("Figure/FigDraw00_grp140_c2.png");
		c2->SaveAs("Figure/FigDraw00_grp140_c2.pdf");
		c3->cd();
		c3->SaveAs("Figure/FigDraw00_grp140_c3.png");
		c3->SaveAs("Figure/FigDraw00_grp140_c3.pdf");
		c6->cd();
		c6->SaveAs("Figure/FigDraw00_grp140_c6.png");
		c6->SaveAs("Figure/FigDraw00_grp140_c6.pdf");
		c5->cd();
		c5->SaveAs("Figure/FigDraw00_grp140_c5.png");
		c5->SaveAs("Figure/FigDraw00_grp140_c5.pdf");
		c4->cd();
		c4->SaveAs("Figure/FigDraw00_grp140_c4.png");
		c4->SaveAs("Figure/FigDraw00_grp140_c4.pdf");
	}


	return;

	/*
	*/

	/*
	*/



}

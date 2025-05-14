void pileup_template(){
    std::string base = "Amp/minbias";
    TFile *file = new TFile(Form("%s/histograms.stabilityFix.root",base.c_str()),"read");
    TFile *foutput = new TFile(Form("%s/pileup.root",base.c_str()),"recreate");
    // TH2D *h2 = (TH2D*)file->Get("hLbAmp0")->Clone("da");
    // TH1D *h1 = h2->ProjectionY();
    TH1D *h_0 = (TH1D*)file->Get("hsumC_energy")->Clone("da");
    TH1D *h_pileup_0 = new TH1D("h_pileup_0","",200,0,25000);
    TH1D *h_1 = (TH1D*)file->Get("hsumA_energy")->Clone("da");
    TH1D *h_pileup_1 = new TH1D("h_pileup_1","",200,0,25000);

    for(int i=0; i<1000000; i++){
        float ran = h_0->GetRandom();
        float ran2 = h_0->GetRandom();
        float ran3 = h_1->GetRandom();
        float ran4 = h_1->GetRandom();
        h_pileup_0->Fill(ran + ran2);
        h_pileup_1->Fill(ran3 + ran4);
    }
    // h_pileup->Scale(0.02/h_pileup->GetEntries());
    // h1->Scale(1.0/h1->GetEntries());
    // h1->Add(h_pileup,-1);
    // h1->SetLineColor(kRed);
    // h1->Draw();
    // h_pileup->Draw("same");
    // gPad->SetLogy();
    foutput->cd();
    h_pileup_0->Write();
    h_pileup_1->Write();
    foutput->Close();
}
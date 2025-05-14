void resolution(){
    std::string base = "Amp/xor";
    TFile *file = new TFile(Form("%s/histograms.stabilityFix.root",base.c_str()),"read");
    TFile *foutput = new TFile(Form("%s/resolution.root",base.c_str()),"recreate");
    // TH2D *h2 = (TH2D*)file->Get("hLbAmp0")->Clone("da");
    // TH1D *h1 = h2->ProjectionY();
    TH1D *h_0 = (TH1D*)file->Get("hsumC_energy")->Clone("da");
    TH1D *h_res_0 = new TH1D("h_res_0","",200,0,25000);
    TH1D *h_1 = (TH1D*)file->Get("hsumA_energy")->Clone("da");
    TH1D *h_res_1 = new TH1D("h_res_1","",200,0,25000);
    // TH2D *h_corr = new TH1D("h_correlation",";Measured;Smeared",200,0,25000,200,0,25000);
    double mean = 1.0;  // Example mean value
    double sigma = 0.25 * mean; // 25% of mean

    TRandom3 randGen(0); // Random generator (seed = 0 for randomness)
    for(int i=0; i<1000000; i++){
        float ran0 = h_0->GetRandom();
        float ran1 = h_1->GetRandom();
        double value = randGen.Gaus(mean, sigma); 
        h_res_0->Fill(ran0*value);
        h_res_1->Fill(ran1*value);
    }
    // h_res_0->Scale(1.0/h_res_0->GetEntries());
    // h_0->Scale(1.0/h_0->GetEntries());
    // //h_res_1->Scale(1.0/h_res_1->GetEntries());
    // h_0->Draw();
    // h_res_0->Draw("same");
    // h_res_0->SetLineColor(kRed);
    // gPad->SetLogy();
    foutput->cd();
    h_res_0->Write();
    h_res_1->Write();
    foutput->Close();
}
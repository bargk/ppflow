void linear_func(){
    std::string data = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp");
    TFile *input = new TFile(Form("%s/histograms.root",data.c_str()),"read");
    TH2D *h0 = (TH2D*)input->Get("h0_h1h2Rat_corr_had2")->Clone("h2");
    TH1D *h1 = h0->ProfileX();
    // Define function
    TF1* line = new TF1("sideC", "[0]/x + [1]/(x*x) + [2]/(x*x*x) + [3]", 10, 4400);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);

    // Fit
    h1->Fit(line, "R");
        // Get fit parameter errors
        double par[4], err[4];
        for (int i = 0; i < 4; i++) {
            par[i] = line->GetParameter(i);
            err[i] = line->GetParError(i);
        }

        // Create upper and lower uncertainty bands (±2 sigma)
        TF1* line_upper = new TF1("sideC_upper", 
                                Form("(%f+%f*2)/x + (%f+%f*2)/(x*x) + (%f+%f*2)/(x*x*x) + (%f+%f*2)", 
                                    par[0], err[0], par[1], err[1], par[2], err[2], par[3], err[3]), 
                                10, 4400);
        line_upper->SetLineColor(kBlue);
        line_upper->SetLineStyle(2);

        TF1* line_lower = new TF1("sideC_lower", 
                                Form("(%f-%f*2)/x + (%f-%f*2)/(x*x) + (%f-%f*2)/(x*x*x) + (%f-%f*2)", 
                                    par[0], err[0], par[1], err[1], par[2], err[2], par[3], err[3]), 
                                10, 4400);
        line_lower->SetLineColor(kBlue);
        line_lower->SetLineStyle(2);

        // Draw everything
        TCanvas *c1 = new TCanvas("c1", "Fit with ±2 Sigma Uncertainty", 800, 600);
        h1->Draw();
        line->Draw("same");
        line_upper->Draw("same");
        line_lower->Draw("same");

        // Save the function
        TFile *output = new TFile("fitted_function_with_2sigma_uncertainty.root", "RECREATE");
        line->Write("sideC");
        line_upper->Write("sideC_upper");
        line_lower->Write("sideC_lower");
        output->Close();
}
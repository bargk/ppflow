//ploting figures for internal note that related to the ZDC section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide/1.5sigma";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ZDC";

TH1D * h_sum[2];
TH1D * hclone;
TFile *file_calibrated;
TFile *file_uncalibrated;

void plot_zdc(){
    file_calibrated = new TFile(Form("%s/zdc_calibrated.root",base.c_str()),"read");
    file_uncalibrated = new TFile(Form("%s/zdc_uncalib.root",base.c_str()),"read");
    TCanvas* c0 = new TCanvas("c0","",4000,2600);
    TLegend *legend0;
    SetAtlasStyle(); 

    //plot uncalibrated ZDC without fit
    //std::vector<std::string> histNames = {"h_c_opposite", "h_a_opposite"};
    std::vector<std::string> histNames = {"h_c", "h_a"};
    std::vector<int> lineColors = {kGreen + 1, kMagenta};
    std::vector<std::string> fileNames = {"side_c_uncalib_nofit.pdf", "side_a_uncalib_nofit.pdf"};
    gSystem->Exec(Form("mkdir -p %s/energy_calib",figures.c_str()));
    // Loop over histograms
    for (int i = 0; i < histNames.size(); i++) {
        // Retrieve histogram from file
        h_sum[i] = (TH1D*)file_uncalibrated->Get(histNames[i].c_str());
        // Set color and style
        h_sum[i]->SetLineColor(lineColors[i]);
        h_sum[i]->SetMarkerColor(lineColors[i]);
        Common::FormatHist(h_sum[i], Common::StandardFormat());
        h_sum[i]->SetMarkerStyle(20);
        h_sum[i]->SetMarkerSize(2);

        // Draw the histogram
        h_sum[i]->Draw();
        h_sum[i]->GetXaxis()->SetRange(3,60);
        h_sum[i]->GetYaxis()->SetTitleOffset(1.1);
        //h_sum[i]->GetFunction("fit")->SetLineWidth(10);
        h_sum[i]->GetFunction("fit")->Delete();

        // Add text annotations
        double X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", 90, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, 90, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", 90, 43); Y -= 0.05;

        // Save the canvas
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        }
    //------------------------------------------------------------------------------------------------------------------

    //plot uncalibrated ZDC with fit
    file_uncalibrated->Clear();
    //histNames = {"h_c_opposite", "h_a_opposite"};
    histNames = {"h_c", "h_a"};
    lineColors = {kGreen + 1, kMagenta};
    fileNames = {"side_c_uncalib_fit.pdf", "side_a_uncalib_fit.pdf"};
    // Loop over histograms
    for (size_t i = 0; i < histNames.size(); i++) {
        // Retrieve histogram from file
        h_sum[i] = (TH1D*)file_uncalibrated->Get(histNames[i].c_str());
        
        // Set color and style
        h_sum[i]->SetLineColor(lineColors[i]);
        h_sum[i]->SetMarkerColor(lineColors[i]);
        Common::FormatHist(h_sum[i], Common::StandardFormat());
        h_sum[i]->SetMarkerStyle(20);
        h_sum[i]->SetMarkerSize(2);

        // Draw the histogram
        h_sum[i]->Draw();
        h_sum[i]->GetXaxis()->SetRange(3,60);
        h_sum[i]->GetYaxis()->SetTitleOffset(1.1);
        //h_sum[i]->GetFunction("fit")->SetLineWidth(10);
        h_sum[i]->GetFunction("fit")->Draw("same");

        // Add text annotations
        double X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", 90, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, 90, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", 90, 43); Y -= 0.05;

        // Save the canvas
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        h_sum[i]->Clear();
        }
    //--------------------------------------------------------------------------------------------------------------
    //Draw 1n peak fit after calibration
    histNames = {"h_c", "h_a"};
    lineColors = {kGreen + 1, kMagenta};
    fileNames = {"side_c_calib_fit.pdf", "side_a_calib_fit.pdf"};
    // Loop over histograms
    for (size_t i = 0; i < histNames.size(); i++) {
        // Retrieve histogram from file
        h_sum[i] = (TH1D*)file_calibrated->Get(histNames[i].c_str());
        
        // Set color and style
        h_sum[i]->SetLineColor(lineColors[i]);
        h_sum[i]->SetMarkerColor(lineColors[i]);
        Common::FormatHist(h_sum[i], Common::StandardFormat());
        h_sum[i]->SetMarkerStyle(20);
        h_sum[i]->SetMarkerSize(2);

        // Draw the histogram
        h_sum[i]->GetXaxis()->SetRangeUser(1000,8000);
        //fit for 1n peak
        TF1 *fit_1n = new TF1("fit_1n","gaus");
        fit_1n->SetLineColor(kRed);
        fit_1n->SetLineWidth(5);
        if(i ==0 ){
            fit_1n->SetRange(2350,3150);
        }
        else{
            fit_1n->SetRange(2200,3100);    
        }
        h_sum[i]->Fit("fit_1n","Rq");
        h_sum[i]->Draw();
        

        // Add text annotations
        double X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", 90, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, 90, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", 90, 43); Y -= 0.10;
        Common::myText2(X, Y, 1, Form("#mu = : %i",(int)fit_1n->GetParameter(1)), 90, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, Form("#sigma /#mu = : %.2f",fit_1n->GetParameter(2)/fit_1n->GetParameter(1)), 90, 43); Y -= 0.05;

        // Save the canvas
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        h_sum[i]->Clear();
        }
}
//ploting figures for internal note that related to the ZDC calibration section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"


std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ZDC";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ZDC";

TH1D * h_sum[8];
THStack *hs0;
TH1D * hclone;
TFile *file_calibrated;
TFile *file_uncalibrated;
double X, Y;
int Size = 100;
void plot_ana_zdc_calibration(){
    file_calibrated = new TFile(Form("%s/zdc_calibrated.root",base.c_str()),"read");
    file_uncalibrated = new TFile(Form("%s/zdc_uncalib_itr1.root",base.c_str()),"read");
    TCanvas* c0 = new TCanvas("c0","",4000,2600);
    TLegend *legend0;
    SetAtlasStyle(); 

    //plot uncalibrated ZDC without fit
    //std::vector<std::string> histNames = {"h_c_opposite", "h_a_opposite"};
    std::vector<std::string> histNames = {"h_c", "h_a"};
    std::vector<int> lineColors = {kBlack, kMagenta};
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
        //h_sum[i]->GetXaxis()->SetRange(3,90);
        h_sum[i]->GetXaxis()->SetTitle("ZDC Sum [a.u]");
        h_sum[i]->GetYaxis()->SetTitle("Events");
        h_sum[i]->GetYaxis()->SetTitleOffset(1.1);
        //h_sum[i]->GetFunction("fit")->SetLineWidth(10);
        h_sum[i]->GetFunction("fit")->Delete();

        // Add text annotations
        X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.05;

        // Save the canvas
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        }
    //------------------------------------------------------------------------------------------------------------------
    //plot zdc sides top on each other
    h_sum[0] = (TH1D*)file_uncalibrated->Get(histNames[0].c_str());
    h_sum[1] = (TH1D*)file_uncalibrated->Get(histNames[1].c_str());
    h_sum[0]->SetLineColor(lineColors[0]);
    h_sum[1]->SetLineColor(lineColors[1]);
    Common::FormatHist(h_sum[1], Common::StandardFormat());
    hs0 = new THStack("sides_uncalib", ";ZDC Sum [a.u];Events");
    hs0->Add(h_sum[0]);
    hs0->Add(h_sum[1]);
    hs0->Draw("nostack");
    X = 0.5, Y = 0.85;
    Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
    Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
    Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.05;
    legend0 = new TLegend(0.5,0.45,0.7,0.65);
    legend0->AddEntry(h_sum[0], "side C", "l");
    legend0->AddEntry(h_sum[1], "side A", "l");
    legend0->SetBorderSize(0);
    legend0->Draw();
    c0->SaveAs(Form("%s/energy_calib/sides_uncalib.pdf", figures.c_str()));
    c0->Clear();
    legend0->Clear();
    
    //plot uncalibrated ZDC with fit
    file_uncalibrated->Clear();
    //histNames = {"h_c_opposite", "h_a_opposite"};
    histNames = {"h_c", "h_a"};
    //lineColors = {kGreen + 1, kMagenta};
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
        //h_sum[i]->GetXaxis()->SetRange(3,90);
        h_sum[i]->GetYaxis()->SetTitleOffset(1.1);
        //h_sum[i]->GetFunction("fit")->SetLineWidth(10);
        h_sum[i]->GetFunction("fit")->Draw("same");

        // Add text annotations
        X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.05;

        // Save the canvas
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        h_sum[i]->Clear();
        }
    //--------------------------------------------------------------------------------------------------------------
    //Draw 1n peak fit after calibration
    histNames = {"h_c", "h_a"};
    //lineColors = {kGreen + 1, kMagenta};
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
        h_sum[i]->GetXaxis()->SetTitle("Energy [GeV]");

        // Draw the histogram
        //h_sum[i]->GetXaxis()->SetRangeUser(1000,8000);
        //fit for 1n peak
        TF1 *fit_1n = new TF1("fit_1n","gaus");
        fit_1n->SetLineColor(kRed);
        fit_1n->SetLineWidth(5);
        if(i ==0 ){
            fit_1n->SetRange(2250,3000);
        }
        else{
            fit_1n->SetRange(2250,3000);    
        }
        h_sum[i]->Fit("fit_1n","Rq");
        h_sum[i]->Draw();
        

        // Add text annotations
        X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.10;
        Common::myText2(X, Y, 1, Form("#mu = : %i",(int)fit_1n->GetParameter(1)), Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, Form("#sigma /#mu = : %.2f",fit_1n->GetParameter(2)/fit_1n->GetParameter(1)), Size, 43); Y -= 0.05;

        // Save the canvas
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        h_sum[i]->Clear();
        }

    //--------------------------------------------------------------------------------------------------------------
    //Draw each side distribtuion after several iterations
    file_uncalibrated = new TFile(Form("%s/zdc_uncalib_itr2.root",base.c_str()),"read");
    TFile *file_uncalibrated1 = new TFile(Form("%s/zdc_uncalib_itr4.root",base.c_str()),"read");
    TFile *file_uncalibrated2 = new TFile(Form("%s/zdc_uncalib_itr7.root",base.c_str()),"read");
    TFile *file_uncalibrated3 = new TFile(Form("%s/zdc_uncalib_itr9.root",base.c_str()),"read");
    histNames = {"h_c", "h_a"};
    fileNames = {"side_c_iterations.pdf", "side_a_iterations.pdf"};
    for (size_t i = 0; i < histNames.size(); i++) {
        h_sum[0] = (TH1D*)file_uncalibrated->Get(histNames[i].c_str());
        h_sum[1] = (TH1D*)file_uncalibrated1->Get(histNames[i].c_str());
        h_sum[2] = (TH1D*)file_uncalibrated2->Get(histNames[i].c_str());
        h_sum[3] = (TH1D*)file_uncalibrated3->Get(histNames[i].c_str());
        
        // Set color and style
        for(int idx = 0; idx<4; idx++){
            h_sum[idx]->SetLineColor(kBlue+2*idx);
            h_sum[idx]->SetMarkerColor(kBlue+2*idx);
            Common::FormatHist(h_sum[idx], Common::StandardFormat());
            h_sum[idx]->SetMarkerStyle(20);
            h_sum[idx]->SetMarkerSize(2);
            h_sum[idx]->GetFunction("fit")->Delete();
        }
        h_sum[3]->Draw();
        h_sum[0]->Draw("same");
        h_sum[1]->Draw("same");
        h_sum[2]->Draw("same");

        legend0 = new TLegend(0.5,0.45,0.7,0.65);
        //legend0 = new TLegend(0.5,0.4,0.87,0.8);
        legend0->AddEntry(h_sum[0], "1 iteration", "l");
        legend0->AddEntry(h_sum[1], "3 iterations", "l");
        legend0->AddEntry(h_sum[2], "6 iterations", "l");
        legend0->AddEntry(h_sum[3], "9 iterations", "l");
        legend0->Draw();
        // Add text annotations
        X = 0.5, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.10;
        c0->SaveAs(Form("%s/energy_calib/%s", figures.c_str(), fileNames[i].c_str()));
        c0->Clear();
        h_sum[i]->Clear();
        legend0->Clear();
    }

    //zdc weights convergence plots
    TH1D *h_zdc_weights[2][4];
    std::vector<TVectorD*> weights;
    TVectorD* multiplied_weights = new TVectorD(4);
        //initialize histograms
        for(int side =0; side<2; side++){
            for(int mod =0; mod<4; mod++){
                h_zdc_weights[side][mod] = new TH1D(Form("h%i_%i",side,mod),";Iteration;Weight [#frac{GeV}{ADC}]",10,1.5,10.5);
                h_zdc_weights[side][mod]->SetMarkerColor(lineColors[side]);
                h_zdc_weights[side][mod]->SetLineColor(lineColors[side]);
                h_zdc_weights[side][mod]->SetMarkerStyle(20);
                h_zdc_weights[side][mod]->SetMarkerSize(5);
            }
        }
        for(int side =0; side<2; side++){
            for(int idx=1; idx<=10; idx++){
                file_uncalibrated = new TFile(Form("%s/zdcWeights_side%i_itr%i.root",base.c_str(),side,idx),"read");
                TVectorD* weights_itr =(TVectorD*)file_uncalibrated->Get("gains_avg");
                TVectorD* weights_itr_err =(TVectorD*)file_uncalibrated->Get("gains_std");
                for(int mod =0; mod<4; mod++){
                    h_zdc_weights[side][mod]->SetBinContent(idx, (*weights_itr)[mod]);      
                    h_zdc_weights[side][mod]->SetBinError(idx, (*weights_itr_err)[mod]);      
                }
            }
        }
        

    //drawing
    for(int mod =0; mod<4; mod++){ 
        THStack *hs = new THStack(Form("hs%i", mod), Form(";Iteration;w_{%i} [#frac{GeV}{a.u}]",mod));
        gStyle->SetErrorX(0.001);
        double maximum;
        maximum =(h_zdc_weights[0][mod]->GetMaximum() > h_zdc_weights[1][mod]->GetMaximum()) ? h_zdc_weights[0][mod]->GetMaximum() : h_zdc_weights[1][mod]->GetMaximum();
        hs->SetMaximum(1.3*maximum);
        hs->Add(h_zdc_weights[0][mod]);
        hs->Add(h_zdc_weights[1][mod]);
        hs->Draw("nostack;p;E1");
        hs->GetYaxis()->SetTitleOffset(0.95);
        hs->GetXaxis()->SetTitleOffset(1.0);
        legend0 = new TLegend(0.7,0.75,0.9,0.9);
        legend0->AddEntry(h_zdc_weights[0][mod], "side C", "lep");
        legend0->AddEntry(h_zdc_weights[1][mod], "side A", "lep");
        legend0->SetBorderSize(0);
        legend0->Draw();

        X = 0.15, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
        Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.10;
        c0->SaveAs(Form("%s/energy_calib/weights_mod%i.pdf", figures.c_str(),mod));
        c0->Clear();
    }

    //drawing calibrated sides on top of each other
    file_calibrated = new TFile(Form("%s/zdc_calibrated.root",base.c_str()),"read");
    h_sum[0] = (TH1D*)file_calibrated->Get(histNames[0].c_str());
    h_sum[1] = (TH1D*)file_calibrated->Get(histNames[1].c_str());
    h_sum[0]->SetLineColor(lineColors[0]);
    h_sum[1]->SetLineColor(lineColors[1]);
    
    Common::FormatHist(h_sum[1], Common::StandardFormat());
    THStack *hs1 = new THStack("sides_calib", ";ZDC Sum [GeV];Events");
    hs1->Add(h_sum[0]);
    hs1->Add(h_sum[1]);
    hs1->Draw("nostack");
    //hs1->GetXaxis()->SetNdivisions(505);
    c0->SetBottomMargin(0.15);  // Increase bottom margin (default is ~0.1)
    c0->SetLeftMargin(0.15);  // Increase bottom margin (default is ~0.1)
    X = 0.5, Y = 0.85;
    Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
    Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
    Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", Size, 43); Y -= 0.05;
    legend0 = new TLegend(0.5,0.45,0.7,0.65);
    legend0->AddEntry(h_sum[0], "side C", "l");
    legend0->AddEntry(h_sum[1], "side A", "l");
    legend0->SetBorderSize(0);
    legend0->Draw();
    c0->SaveAs(Form("%s/energy_calib/sides_calib.pdf", figures.c_str()));
    c0->Clear();




  
    
}


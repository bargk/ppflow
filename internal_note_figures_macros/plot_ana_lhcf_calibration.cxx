//ploting information about the zdc in lhcf run after applying calibration factors
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

//#define pileup
#define resolution
//#define sides_calib

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ZDC";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/energy_scale_studies/applying_calibration_factors";
double X, Y;
int Size = 17;

void plot_ana_lhcf_calibration(){
    SetAtlasStyle();
    std::vector<int> lineColors = {kRed+3, kBlue+3};
    std::vector<string> triggers = {"/","/xor/","/minbias/"};
    std::vector<string> trigger_name = {"AND","XOR","minbias"};
    std::vector<string> hist_name = {"hsumC_energy","hsumA_energy"};
    std::vector<string> sides_name = {"C","A"};
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TLegend *legend0;
    TCanvas *c0 = new TCanvas("c0");
    float mu = 0.02;
    //pileup drawing
    #ifdef pileup
        for(int side =0; side<2; side++){
            for(int i=0; i<trigger_name.size(); i++){
                TFile *ftrue = new TFile(Form("%s%shistograms.stabilityFix.root",base.c_str(), triggers.at(i).c_str()),"read");
                TFile *fpileup = new TFile(Form("%s%spileup.root",base.c_str(), triggers.at(i).c_str()),"read");
                TH1D *h1 = (TH1D*)ftrue->Get(Form("%s",hist_name.at(side).c_str()))->Clone(Form("clone_%i_%s",side,trigger_name.at(i).c_str()));
                TH1D *h_pileup = (TH1D*)fpileup->Get(Form("h_pileup_%i",side))->Clone(Form("clone_h_pileup_%i_%s",side,trigger_name.at(i).c_str()));
                h_pileup->Scale(mu/h_pileup->GetEntries()); h_pileup->SetMarkerColor(lineColors.at(0)); h_pileup->SetMarkerStyle(20); h_pileup->SetMarkerSize(0.8);
                h1->Scale(1.0/h1->GetEntries());            h1->SetMarkerColor(lineColors.at(1)); h1->SetMarkerStyle(4); h1->SetMarkerSize(0.8);
                THStack *hs = new THStack(Form("hs_%i_%s",side,trigger_name.at(i).c_str()),";ZDC energy [GeV];Normalized Events");
                // Define two pads: upper for histogram, lower for ratio
                TPad* pad1 = new TPad("pad1", "Main Plot", 0, 0.3, 1, 1);
                TPad* pad2 = new TPad("pad2", "Ratio Plot", 0, 0.0, 1, 0.3);
                // Set margins
                pad1->SetBottomMargin(0.01);  // Remove bottom margin for main plot
                pad2->SetTopMargin(0);     // Remove top margin for ratio plot
                pad2->SetBottomMargin(0.3); // Space for ratio axis labels

                c0->cd();
                // Draw pads
                pad1->Draw();
                pad2->Draw();
                hs->Add(h1);
                hs->Add(h_pileup);
                pad1->cd();
                hs->Draw("nostack");
                hs->GetXaxis()->SetLimits(0,20000);
                hs->GetXaxis()->SetNdivisions(505);
                // Create a small pad inside the large canvas
                TPad* pad = new TPad("pad", "Small Pad", 0.36, 0.1, 0.9, 0.7);
                pad->SetFillColor(kWhite); // Background color
                pad->Draw();
                pad->cd();  // Set focus to the small pad
                THStack *hs_small = new THStack(Form("hs_%i_%s",side,trigger_name.at(i).c_str()),";ZDC energy [GeV];Normalized Events");
                hs_small->Add(h1);
                hs_small->Add(h_pileup);
                hs_small->Draw("nostack");
                hs_small->SetMinimum(0);
                hs_small->SetMaximum(mu/20);
                hs_small->GetXaxis()->SetLimits(0,20000);
                hs_small->GetXaxis()->SetNdivisions(505);
                // Draw a vertical line at x = 6.8 TeV
                TLine* line = new TLine(6800, gPad->GetUymin(), 6800, gPad->GetUymax());
                line->SetLineColor(kBlack);
                line->SetLineStyle(2);  // Dashed line
                line->SetLineWidth(2);
                line->Draw();

                // Draw the lower pad (Ratio plot)
                TH1D *hRatio = (TH1D*)h_pileup->Clone("hRatio");
                hRatio->Divide(h1);
                pad2->cd();
                hRatio->SetMarkerColor(kBlack);
                hRatio->SetStats(0);  // No stats box
                hRatio->SetTitle(""); // No title
                hRatio->GetYaxis()->SetRangeUser(0.0, 0.5);  // Define the ratio scale
                hRatio->GetYaxis()->SetTitleSize(0.1);
                hRatio->GetYaxis()->SetTitleOffset(0.5);
                hRatio->GetYaxis()->SetLabelSize(0.08);
                hRatio->GetYaxis()->SetTitle("#frac{2 pileup}{data}");
                hRatio->GetXaxis()->SetTitle("ZDC energy [GeV]");
                hRatio->GetXaxis()->SetTitleSize(0.1);
                hRatio->GetXaxis()->SetLabelSize(0.1);
                hRatio->Draw("L");
                TLine* line1 = new TLine(6800, 0.0, 6800, 0.5);
                line1->SetLineColor(kBlack);
                line1->SetLineStyle(2);  // Dashed line
                line1->SetLineWidth(2);
                line1->Draw();
                c0->cd();
                X = 0.3, Y = 0.9;
                Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("Side %s",sides_name.at(side).c_str()), Size, 43); Y -= 0.05;
                legend0 = new TLegend(0.6,0.8,0.8,0.88);
                legend0->SetBorderSize(0);
                legend0->AddEntry(h1, "Measured", "p");
                legend0->AddEntry(h_pileup, "2 pileup", "p");
                legend0->Draw();
                c0->SaveAs(Form("%s/pileup%i_%s.pdf",figures.c_str(),side,trigger_name.at(i).c_str()));
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------
    //plot correlation Amp with units of energy
    c0->Clear();
    for(int i=0; i<trigger_name.size() -1; i++){
        TFile *input = new TFile(Form("%s%shistograms.stabilityFix.root",base.c_str(), triggers.at(i).c_str()),"read");
        TH2D *h2 = (TH2D*)input->Get(Form("hAmpCorr_energy"))->Clone(Form("clone_hAmpCorr_energy%s",trigger_name.at(i).c_str()));
        h2->GetXaxis()->SetTitle("Side C [GeV]");
        h2->GetYaxis()->SetTitle("Side A [GeV]");
        h2->GetZaxis()->SetTitle("Events");
        h2->GetXaxis()->SetRangeUser(0,5000);
        h2->GetYaxis()->SetRangeUser(0,5000);
        h2->GetXaxis()->SetNdivisions(505);
        h2->GetYaxis()->SetNdivisions(505);
        c0->cd();
        h2->Draw("colz");
        gPad->SetLogy(0);
        gPad->SetLogz();
        gPad->SetRightMargin(0.15);
        TLine* line1 = new TLine(1000, 0.0, 1000, 5000);
        TLine* line2 = new TLine(0.0, 1000, 5000, 1000);
        line1->SetLineColor(kBlack);
        line1->SetLineStyle(2);  // Dashed line
        line1->SetLineWidth(2);
        line1->Draw();
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(2);  // Dashed line
        line2->SetLineWidth(2);
        line2->Draw("same");
        X = 0.6, Y = 0.85;
        Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
        Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
        Common::myText2(X, Y, 1, Form("%s trigger", trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
        std::string figures_final = Form("%s/correlation",figures.c_str());
        gSystem->Exec(Form("mkdir -p %s",figures_final.c_str()));
        c0->SaveAs(Form("%s/sides_corr_%s.pdf",figures_final.c_str(),trigger_name.at(i).c_str()));
        h2->Reset();
    }
    //------------------------------------------------------------------------------------------------------------------
    #ifdef resolution
    
        for(int irbin=0; irbin<2; irbin++){
            for(int side =0; side<2; side++){
                for(int i=0; i<trigger_name.size(); i++){
                    TFile *ftrue = new TFile(Form("%s%shistograms.stabilityFix.root",base.c_str(), triggers.at(i).c_str()),"read");
                    TFile *fpileup = new TFile(Form("%s%sresolution.root",base.c_str(), triggers.at(i).c_str()),"read");
                    TH1D *h1 = (TH1D*)ftrue->Get(Form("%s",hist_name.at(side).c_str()))->Clone(Form("clone_%i_%s",side,trigger_name.at(i).c_str()));
                    TH1D *h_pileup = (TH1D*)fpileup->Get(Form("h_res_%i",side))->Clone(Form("clone_h_pileup_%i_%s",side,trigger_name.at(i).c_str()));
                    h_pileup->Scale(1.0/h_pileup->GetEntries()); h_pileup->SetMarkerColor(lineColors.at(0)-1); h_pileup->SetMarkerStyle(20); h_pileup->SetMarkerSize(0.8);
                    h1->Scale(1.0/h1->GetEntries());            h1->SetMarkerColor(lineColors.at(1)); h1->SetMarkerStyle(4); h1->SetMarkerSize(0.8);
                    THStack *hs = new THStack(Form("hs_%i_%s",side,trigger_name.at(i).c_str()),";ZDC energy [GeV];Normalized Events");
                    
                    if(irbin ==1){
                        h_pileup->Rebin(4);
                        h1->Rebin(4);
                    }

                    c0->cd();
                    hs->Add(h1);
                    hs->Add(h_pileup);
                    gPad->SetLogy();
                    //pad1->cd();
                    hs->Draw("nostack");
                    hs->GetXaxis()->SetLimits(0,20000);
                    hs->GetXaxis()->SetNdivisions(505);
                    double ymin = TMath::Power(10, gPad->GetUymin());  // Convert log scale to linear
                    double ymax = TMath::Power(10, gPad->GetUymax());  // Convert log scale to linear
                    TLine* line = new TLine(6800, ymin, 6800, ymax);
                    line->SetLineColor(kBlack);
                    line->SetLineStyle(2);  // Dashed line
                    line->SetLineWidth(2);
                    line->Draw();
                    X = 0.6, Y = 0.9;
                    Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                    Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, Form("Side %s",sides_name.at(side).c_str()), Size, 43); Y -= 0.05;
                    legend0 = new TLegend(0.65,0.68,0.8,0.78);
                    legend0->SetBorderSize(0);
                    legend0->SetTextSize(0.045);
                    legend0->AddEntry(h1, "Measured", "p");
                    legend0->AddEntry(h_pileup, "Smeared", "p");
                    legend0->Draw();
                    if(irbin ==0) c0->SaveAs(Form("%s/resolution%i_%s.pdf",figures.c_str(),side,trigger_name.at(i).c_str()));
                    if(irbin ==1) c0->SaveAs(Form("%s/resolution%i_%s_rebin.pdf",figures.c_str(),side,trigger_name.at(i).c_str())); 
                    c0->Clear();
                    c0->cd();
                    // Draw Ratio plot
                    TH1D *hRatio = (TH1D*)h_pileup->Clone("hRatio");
                    hRatio->Divide(h1);
                    for(int bin =0; bin < hRatio->GetNbinsX(); bin++){
                        hRatio->SetBinError(bin,0);
                    }
                    hRatio->SetLineColor(kBlack);
                    hRatio->GetXaxis()->SetRangeUser(1000, 6800); 
                    hRatio->GetYaxis()->SetTitleOffset(0.95);
                    hRatio->GetYaxis()->SetTitle("#frac{Smeared}{Measaured}");
                    hRatio->GetXaxis()->SetTitle("ZDC energy [GeV]");
                    hRatio->Draw();
                    hRatio->GetXaxis()->SetNdivisions(505);
                    hRatio->SetMaximum(2);
                    hRatio->SetMinimum(0.8);
                    gPad->SetLogy(0);
                    gPad->SetGrid();
                    if(irbin ==0) c0->SaveAs(Form("%s/resolution_ratio%i_%s.pdf",figures.c_str(),side,trigger_name.at(i).c_str()));
                    if(irbin ==1) c0->SaveAs(Form("%s/resolution_ratio%i_%s_rebin.pdf",figures.c_str(),side,trigger_name.at(i).c_str()));
                    gPad->SetGrid(0,0);
                }
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

    //plot sides after calibration
    std::vector<string> figures_name = {"sides_AND.pdf","sides_XOR.pdf","sides_minbias.pdf"};
    lineColors.clear();
    lineColors = {kBlack, kMagenta};
    #ifdef sides_calib
        for(int i =0; i<triggers.size(); i++){
            TFile *input = new TFile(Form("%s%shistograms.stabilityFix.root",base.c_str(), triggers.at(i).c_str()),"read");
            THStack *hs = new THStack("hs",";ZDC sum [GeV];Events");
            TH1D *h1_0 = (TH1D*)input->Get(Form("%s",hist_name.at(0).c_str()))->Clone("side0");
            TH1D *h1_1 = (TH1D*)input->Get(Form("%s",hist_name.at(1).c_str()))->Clone("side1");
            h1_0->SetLineColor(lineColors.at(0)); h1_0->SetLineWidth(2);  
            h1_1->SetLineColor(lineColors.at(1)); h1_1->SetLineWidth(2);
            c0->cd();
            hs->Add(h1_0);
            hs->Add(h1_1);
            hs->Draw("nostack");
            hs->GetXaxis()->SetLimits(0,20000);
            hs->GetXaxis()->SetNdivisions(505);
            gPad->SetLogy();
            X = 0.6, Y = 0.85;
            Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
            Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, Form("%s trigger",trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
            legend0 = new TLegend(0.2,0.2,0.44,0.43);
            legend0->SetBorderSize(0);
            legend0->AddEntry(h1_0, "Side C", "l");
            legend0->AddEntry(h1_1, "Side A", "l");
            legend0->Draw();
            c0->SaveAs(Form("%s/%s",figures.c_str(),figures_name.at(i).c_str()));
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

}
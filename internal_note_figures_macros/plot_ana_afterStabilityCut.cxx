//ploting figures for internal note that related stability cut
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

#define sides_uncalib_cuts 
#define sides_uncalib_cuts_ratio 
#define stability_error
#define average_adc 

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ZDC";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/energy_scale_studies/StabilityCut";
double X, Y;
int Size = 20;

void plot_ana_afterStabilityCut(){
    SetAtlasStyle();
    std::vector<int> lineColors = {kRed+3, kBlue+3};
    std::vector<string> triggers = {"/","/xor/","/minbias/"};
    std::vector<string> trigger_name = {"AND","XOR","minbias"};
    std::vector <string> side_name ={"C", "A"};
    std::vector <string> module_name ={"HAD1", "HAD2", "HAD3"};
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TLegend *legend0;
    TCanvas *c0 = new TCanvas("c0");
    //plot zdc sides before and after cut
    #ifdef sides_uncalib_cuts
        for(int side =0; side<2; side++){
            for(int i=0; i<trigger_name.size() -1; i++){
                TFile *input_before = new TFile(Form("%s%shistograms.Had2Cut.root",base.c_str(), triggers.at(i).c_str()),"read");
                TFile *input_after = new TFile(Form("%s%shistograms.stabilityFix.root",base.c_str(), triggers.at(i).c_str()),"read");
                THStack *hs = new THStack("hs",";ZDC sum [ADC];Events");
                // TH2D *h2_0 = (TH2D*)input_before->Get(Form("hLbAmp%i",side))->Clone("before");
                // TH2D *h2_1 = (TH2D*)input_after->Get(Form("hLbAmp%i",side))->Clone("after");
                TH1D *h2_0_proj =(TH1D*)input_before->Get(Form("hsum%s_uncalib",side_name.at(side).c_str())); h2_0_proj->SetLineColor(lineColors.at(0)); h2_0_proj->SetLineWidth(2);  
                TH1D *h2_1_proj =(TH1D*)input_after->Get(Form("hsum%s_uncalib",side_name.at(side).c_str())); h2_1_proj->SetLineWidth(2);
                c0->cd();
                hs->Add(h2_0_proj);
                hs->Add(h2_1_proj);
                hs->Draw("nostack");
                gPad->SetLogy();
                X = 0.6, Y = 0.85;
                Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("Side %s ",side_name.at(side).c_str()), Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("%s trigger",trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
                legend0 = new TLegend(0.2,0.2,0.44,0.43);
                legend0->SetBorderSize(0);
                legend0->AddEntry(h2_0_proj, "Before scaling", "l");
                legend0->AddEntry(h2_1_proj, "After scaling", "l");
                legend0->Draw();
                c0->SaveAs(Form("%s/side%i_%s.pdf",figures.c_str(),side,trigger_name.at(i).c_str()));
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

//plot the modules <ADC> for different triggers
lineColors = {kBlue,kRed, kOrange+3};
gPad->SetLogy(0);
#ifdef average_adc
    for(int side =0; side <2; side++){
        for(int mod =1; mod<4; mod++){
            THStack *hs = new THStack("hs",";LB;#LT ADC#GT");
            TFile *input_0 = new TFile(Form("%s/histograms.stabilityFix.root",base.c_str()),"read");
            //TFile *input_1 = new TFile(Form("%s/minbias/histograms.Had2Cut.root",base.c_str()),"read");
            TFile *input_2 = new TFile(Form("%s/xor/histograms.stabilityFix.root",base.c_str()),"read");

            TH2D *h2_00 = (TH2D*)input_0->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_AND",side,mod));
            //TH2D *h2_01 = (TH2D*)input_1->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_minbias",side,mod));
            TH2D *h2_02 = (TH2D*)input_2->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_XOR",side,mod));
            TH1D *h2_00_profx = h2_00->ProfileX(); h2_00_profx->SetMarkerColor(lineColors.at(0));   h2_00_profx->SetLineColor(kBlack);
            //TH1D *h2_01_profx = h2_01->ProfileX(); h2_01_profx->SetMarkerColor(lineColors.at(1));   h2_01_profx->SetLineColor(kBlack);
            TH1D *h2_02_profx = h2_02->ProfileX(); h2_02_profx->SetMarkerColor(lineColors.at(2));   h2_02_profx->SetLineColor(kBlack);
            h2_00_profx->SetMarkerStyle(20);
            //h2_01_profx->SetMarkerStyle(20);
            h2_02_profx->SetMarkerStyle(20);
            h2_00_profx->SetMarkerSize(0.5);
            //h2_01_profx->SetMarkerSize(0.5);
            h2_02_profx->SetMarkerSize(0.5);
            hs->Add(h2_00_profx);
            //hs->Add(h2_01_profx);
            hs->Add(h2_02_profx);
            hs->SetMinimum(0.7*h2_00_profx->GetBinContent(1500));
            hs->SetMaximum(1.3*h2_00_profx->GetBinContent(1500));
            hs->Draw("nostack");
            hs->GetXaxis()->SetLimits(1816,4755);
            X = 0.2, Y = 0.85;
            Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
            Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, Form("%s%s",module_name.at(mod-1).c_str(),side_name.at(side).c_str()), Size, 43); Y -= 0.05;
            legend0 = new TLegend(0.7,0.72,0.88,0.88);
            legend0->SetBorderSize(0);
            legend0->AddEntry(h2_00_profx, "AND trigger", "lep");
            //legend0->AddEntry(h2_01_profx, "minbias trigger", "lep");
            legend0->AddEntry(h2_02_profx, "XOR trigger", "lep");
            legend0->Draw();

            c0->SaveAs(Form("%s/HAD%i_side%i.pdf",figures.c_str(),mod,side));
        }
    }
#endif
//------------------------------------------------------------------------------------------------------------------

    //plot relative uncertainty with respect to "center" block
    TCanvas *c1 = new TCanvas("c1");
    std::vector<int> lb_center = {Bins::GetLbIndex(1816,1887), Bins::GetLbIndex(1902,1928), Bins::GetLbIndex(2812,2996), Bins::GetLbIndex(2996,3109), Bins::GetLbIndex(3218,3678)};
    TFile *flb_mean_and = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/mean_block1.root");
    TFile *flb_mean_xor = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/xor/mean_block1.root");
    #ifdef stability_error
        for(int side =0; side <2; side++){
            for(int mod =1; mod<4; mod++){
                THStack *hs = new THStack("hs",";LB;#frac{#LT ADC#GT}{#LT#LT ADC#GT#GT_{center}}");
                TFile *input_0 = new TFile(Form("%s/histograms.stabilityFix.root",base.c_str()),"read");
                //TFile *input_1 = new TFile(Form("%s/minbias/histograms.Had2Cut.root",base.c_str()),"read"); 
                TFile *input_2 = new TFile(Form("%s/xor/histograms.stabilityFix.root",base.c_str()),"read");
                
                TH2D *h2_00 = (TH2D*)input_0->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_AND",side,mod));
                TH2D *h2_02 = (TH2D*)input_2->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_XOR",side,mod));
                TH1D *h2_00_profx = h2_00->ProfileX(); h2_00_profx->SetMarkerColor(lineColors.at(0));   h2_00_profx->SetLineColor(kBlack);
                TH1D *h2_02_profx = h2_02->ProfileX(); h2_02_profx->SetMarkerColor(lineColors.at(2));   h2_02_profx->SetLineColor(kBlack);
                //get the mean of lb blocks where lhcf was at cetner
                TH1D* h_mean_and = (TH1D*)flb_mean_and->Get(Form("h_mean_side%i_mod%i",side,mod));
                TH1D* h_mean_xor = (TH1D*)flb_mean_xor->Get(Form("h_mean_side%i_mod%i",side,mod));
                double mean_center_and =0;
                double mean_center_xor =0;
                for(const auto& lumi : lb_center){
                    mean_center_and += h_mean_and->GetBinContent(lumi + 1)/lb_center.size();
                    mean_center_xor += h_mean_xor->GetBinContent(lumi + 1)/lb_center.size();
                }
                h2_00_profx->Scale(1.0/mean_center_and);
                h2_02_profx->Scale(1.0/mean_center_xor);
                hs->Add(h2_00_profx);
                hs->Add(h2_02_profx);
                hs->SetMinimum(0.95);
                hs->SetMaximum(1.1);
                c0->cd();
                hs->Draw("nostack");
                TLine *line = new TLine(1816, 1, 4755, 1); // Horizontal line at y=1
                line->SetLineColor(kBlack); // Set line color
                line->SetLineWidth(2); // Set line width
                line->SetLineStyle(2); // Dashed line
                line->Draw("same");
                hs->GetXaxis()->SetLimits(1816,4755);
                X = 0.2, Y = 0.85;
                Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("%s%s",module_name.at(mod-1).c_str(),side_name.at(side).c_str()), Size, 43); Y -= 0.05;
                legend0 = new TLegend(0.6,0.72,0.88,0.88);
                legend0->SetBorderSize(0);
                legend0->AddEntry(h2_00_profx, "AND trigger", "lep");
                legend0->AddEntry(h2_02_profx, "XOR trigger", "lep");
                legend0->Draw();
                legend0->SetTextSize(0.033);
                c0->SaveAs(Form("%s/stability_uncertainty_HAD%i_side%i.pdf",figures.c_str(),mod,side));
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

    //plot the ratio of the sides before the scale and after the scale
    #ifdef sides_uncalib_cuts_ratio
    for(int side =0; side<2; side++){
        for(int i=0; i<trigger_name.size() -1; i++){
            TFile *input_before = new TFile(Form("%s%shistograms.Had2Cut.root",base.c_str(), triggers.at(i).c_str()),"read");
            TFile *input_after = new TFile(Form("%s%shistograms.stabilityFix.root",base.c_str(), triggers.at(i).c_str()),"read");
            THStack *hs = new THStack("hs",";ZDC sum [ADC];#frac{Before scaling}{After scaling}");
            // TH2D *h2_0 = (TH2D*)input_before->Get(Form("hLbAmp%i",side))->Clone("before");
            // TH2D *h2_1 = (TH2D*)input_after->Get(Form("hLbAmp%i",side))->Clone("after");
            TH1D *h2_0_proj =(TH1D*)input_before->Get(Form("hsum%s_uncalib",side_name.at(side).c_str())); h2_0_proj->SetLineColor(lineColors.at(0)); h2_0_proj->SetLineWidth(2);  
            TH1D *h2_1_proj =(TH1D*)input_after->Get(Form("hsum%s_uncalib",side_name.at(side).c_str())); h2_1_proj->SetLineWidth(2);
            c0->cd();
            h2_0_proj->SetLineColor(kBlack);
            h2_0_proj->Divide(h2_1_proj);
            for(int ibin =0; ibin < h2_0_proj->GetNbinsX(); ibin++){
                h2_0_proj->SetBinError(ibin+1,0);
            }
            hs->Add(h2_0_proj);
            hs->Draw("nostack;L");
            gPad->SetGrid();
            gPad->SetLogy(0);
            c0->SaveAs(Form("%s/side%i_%s_ratio.pdf",figures.c_str(),side,trigger_name.at(i).c_str()));
        }
    }
#endif
//------------------------------------------------------------------------------------------------------------------
gPad->SetGrid(0,0);
}
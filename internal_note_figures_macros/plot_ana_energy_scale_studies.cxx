//ploting figures for internal note that related to the Energy scale section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"
#define average_adc 
#define stability_error 
#define sides_uncalib 
#define modules_bad
#define HAD2_frac
#define HAD2_frac_appendix

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ZDC";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/energy_scale_studies";
double X, Y;
int Size = 20;
void plot_ana_energy_scale_studies(){
    SetAtlasStyle();
    std::vector<int> lineColors = {kBlue, kRed, kOrange+3};
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TLegend *legend0;
    TCanvas *c0 = new TCanvas("c0");
    //plot the modules <ADC> for different triggers
    std::vector <string> side_name ={"C", "A"};
    std::vector <string> module_name ={"HAD1", "HAD2", "HAD3"};
    #ifdef average_adc
        for(int side =0; side <2; side++){
            for(int mod =1; mod<4; mod++){
                THStack *hs = new THStack("hs",";LB;#LT ADC#GT");
                TFile *input_0 = new TFile(Form("%s/histograms.root",base.c_str()),"read");
                TFile *input_1 = new TFile(Form("%s/minbias/histograms.root",base.c_str()),"read"); 
                TFile *input_2 = new TFile(Form("%s/xor/histograms.root",base.c_str()),"read");

                TH2D *h2_00 = (TH2D*)input_0->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_AND",side,mod));
                TH2D *h2_01 = (TH2D*)input_1->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_minbias",side,mod));
                TH2D *h2_02 = (TH2D*)input_2->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_XOR",side,mod));
                TH1D *h2_00_profx = h2_00->ProfileX(); h2_00_profx->SetMarkerColor(lineColors.at(0));   h2_00_profx->SetLineColor(kBlack);
                TH1D *h2_01_profx = h2_01->ProfileX(); h2_01_profx->SetMarkerColor(lineColors.at(1));   h2_01_profx->SetLineColor(kBlack);
                TH1D *h2_02_profx = h2_02->ProfileX(); h2_02_profx->SetMarkerColor(lineColors.at(2));   h2_02_profx->SetLineColor(kBlack);
                h2_00_profx->SetMarkerStyle(20);
                h2_01_profx->SetMarkerStyle(20);
                h2_02_profx->SetMarkerStyle(20);
                h2_00_profx->SetMarkerSize(0.5);
                h2_01_profx->SetMarkerSize(0.5);
                h2_02_profx->SetMarkerSize(0.5);
                hs->Add(h2_00_profx);
                hs->Add(h2_01_profx);
                hs->Add(h2_02_profx);
                hs->SetMinimum(0.7*h2_01_profx->GetBinContent(1500));
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
                legend0->AddEntry(h2_01_profx, "minbias trigger", "lep");
                legend0->AddEntry(h2_02_profx, "XOR trigger", "lep");
                legend0->Draw();

                c0->SaveAs(Form("%s/HAD%i_side%i.pdf",figures.c_str(),mod,side));
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

    //plot relative uncertainty with respect to "true" block  
    std::vector<int> lb_center = {Bins::GetLbIndex(1816,1887), Bins::GetLbIndex(1902,1928), Bins::GetLbIndex(2812,2996), Bins::GetLbIndex(2996,3109), Bins::GetLbIndex(3218,3678)};
    TFile *flb_mean_and = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/mean_block0.root");
    TFile *flb_mean_xor = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/xor/mean_block0.root");
    #ifdef stability_error
        for(int side =0; side <2; side++){
            for(int mod =1; mod<4; mod++){
                THStack *hs = new THStack("hs",";LB;#frac{#LT ADC#GT}{#LT#LT ADC#GT#GT_{center}}");
                TFile *input_0 = new TFile(Form("%s/histograms.root",base.c_str()),"read");
                //TFile *input_1 = new TFile(Form("%s/minbias/histograms.root",base.c_str()),"read"); 
                TFile *input_2 = new TFile(Form("%s/xor/histograms.root",base.c_str()),"read");
                
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
    
    //plot unclibrated distributions of the sides for each trigger
    std::vector<string> trigger_name = {"AND","XOR","minbias"};
    std::vector<string> triggers = {"/","/xor/","/minbias/"};
    std::vector<string> figures_name = {"sides_and.pdf","sides_xor.pdf","sides_minbias.pdf"};
    lineColors.clear();
    lineColors = {kBlack, kMagenta};
    #ifdef sides_uncalib
        for(int i =0; i<triggers.size(); i++){
            TFile *input = new TFile(Form("%s%shistograms.root",base.c_str(), triggers.at(i).c_str()),"read");
            THStack *hs = new THStack("hs",";ZDC sum [ADC];Events");
            TH2D *h2_0 = (TH2D*)input->Get(Form("hLbAmp0"))->Clone("side0");
            TH2D *h2_1 = (TH2D*)input->Get(Form("hLbAmp1"))->Clone("side1");
            TH1D *h2_0_proj = h2_0->ProjectionY(); h2_0_proj->SetLineColor(lineColors.at(0)); h2_0_proj->SetLineWidth(2);  
            TH1D *h2_1_proj = h2_1->ProjectionY(); h2_1_proj->SetLineColor(lineColors.at(1)); h2_1_proj->SetLineWidth(2);
            c0->cd();
            hs->Add(h2_0_proj);
            hs->Add(h2_1_proj);
            hs->Draw("nostack");
            gPad->SetLogy();
            X = 0.6, Y = 0.85;
            Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
            Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, Form("%s trigger",trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
            legend0 = new TLegend(0.2,0.2,0.44,0.43);
            legend0->SetBorderSize(0);
            legend0->AddEntry(h2_0_proj, "Side C", "l");
            legend0->AddEntry(h2_1_proj, "Side A", "l");
            legend0->Draw();
            c0->SaveAs(Form("%s/%s",figures.c_str(),figures_name.at(i).c_str()));
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

    //plot modules distribution under bump
    figures_name = {"modules_bad_sideC.pdf","modules_bad_sideA.pdf"};
    #ifdef modules_bad
        for(int i = 0; i < figures_name.size(); i++){
            TFile *input = new TFile(Form("%s/histograms.root",base.c_str()),"read");
            THStack *hs = new THStack("hs",";ZDC Module [ADC];Events");
            TH2D *h2_1 = (TH2D*)input->Get(Form("h_bad%i_1",i))->Clone("HAD1");
            TH2D *h2_2 = (TH2D*)input->Get(Form("h_bad%i_2",i))->Clone("HAD2");
            TH2D *h2_3 = (TH2D*)input->Get(Form("h_bad%i_3",i))->Clone("HAD3");
            TH1D *h2_1_proj = h2_1->ProjectionY(); h2_1_proj->SetLineColor(kOrange+3); h2_1_proj->SetLineWidth(2);  
            TH1D *h2_2_proj = h2_2->ProjectionY(); h2_2_proj->SetLineColor(kOrange-3); h2_2_proj->SetLineWidth(2);  
            TH1D *h2_3_proj = h2_3->ProjectionY(); h2_3_proj->SetLineColor(kOrange+9); h2_3_proj->SetLineWidth(2);  
            c0->cd();
            hs->Add(h2_1_proj);
            hs->Add(h2_2_proj);
            hs->Add(h2_3_proj);
            hs->Draw("nostack");
            gPad->SetLogy();
            X = 0.6, Y = 0.85;
            Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
            Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, Form("Side %s",side_name.at(i).c_str()), Size, 43); Y -= 0.05;
            legend0 = new TLegend(0.18,0.25,0.3,0.43);
            legend0->SetBorderSize(0);
            legend0->AddEntry(h2_1_proj, "HAD1", "l");
            legend0->AddEntry(h2_2_proj, "HAD2", "l");
            legend0->AddEntry(h2_3_proj, "HAD3", "l");
            legend0->Draw();
            c0->SaveAs(Form("%s/modules_bad_side%i.pdf",figures.c_str(),i));
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------
    
    //plot HAD2 fraction correlation for the 3 triggers
    figures_name = {"HAD2_frac_sideC","HAD2_frac_sideA"};
    #ifdef HAD2_frac
        for(int i =0; i<triggers.size(); i++){
            TFile *input = new TFile(Form("%s%shistograms.root",base.c_str(), triggers.at(i).c_str()),"read");
            for(int side =0; side<2; side++){
                TH2D *h2 = (TH2D*)input->Get(Form("h%i_h2AmpRat_Amp",side))->Clone(Form("HAD2_frac_side%i",side));
                h2->GetXaxis()->SetTitle("ZDC sum [ADC]");
                h2->GetYaxis()->SetTitle("HAD2 Fraction");
                h2->GetZaxis()->SetTitle("Events");
                c0->cd();
                h2->Draw("colz");
                gPad->SetLogy(0);
                gPad->SetLogz();
                gPad->SetRightMargin(0.15);
                X = 0.2, Y = 0.85;
                Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("%s trigger", trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("Side %s",side_name.at(side).c_str()), Size, 43); Y = 0.85;
                Common::myText2(X + 0.4, Y, 1, Form("Mean y: %.2f", h2->GetMean(2)), Size, 43); Y -= 0.05;
                Common::myText2(X + 0.4, Y, 1, Form("StdDev y: %.2f", h2->GetStdDev(2)), Size, 43); Y -= 0.05;
                c0->SaveAs(Form("%s/%s_%s.pdf",figures.c_str(),figures_name.at(side).c_str(),trigger_name.at(i).c_str()));
                h2->Reset();
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

    //plot HAD2 Fraction for the LB ranges for appendix
    std::vector<string> lb_range;
    for(int ilb=0; ilb<Bins::NLB; ilb++){
        lb_range.push_back(Bins::label_lb(ilb));
    }
    #ifdef HAD2_frac_appendix
        for(int side =0; side<2; side++){
            for(int i =0; i<trigger_name.size(); i++){
                TFile *input = new TFile(Form("%s%shistograms.root",base.c_str(), triggers.at(i).c_str()),"read");
                for(int ilb =0; ilb<Bins::NLB; ilb++){
                    TH2D *h2 = (TH2D*)input->Get(Form("h%i_h2AmpRat_Amp_ilb_%i",side,ilb));
                    h2->GetXaxis()->SetTitle("ZDC sum [ADC]");
                    h2->GetYaxis()->SetTitle("HAD2 Fraction");
                    h2->GetZaxis()->SetTitle("Events");
                    c0->cd();
                    h2->Draw("colz");
                    gPad->SetLogy(0);
                    gPad->SetLogz();
                    gPad->SetRightMargin(0.15);
                    X = 0.2; Y = 0.85;
                    Common::myText2(X, Y, 1, Form("LB range : %s",lb_range.at(ilb).c_str()), Size+5, 43);
                    X = 0.6, Y = 0.85;
                    Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                    Common::myText2(X + 0.1, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, Form("%s trigger", trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, Form("Side %s",side_name.at(side).c_str()), Size, 43); Y = 0.85;
                    Common::myText2(X + 0.4, Y, 1, Form("Mean y: %.2f", h2->GetMean(2)), Size, 43); Y -= 0.05;
                    Common::myText2(X + 0.4, Y, 1, Form("StdDev y: %.2f", h2->GetStdDev(2)), Size, 43); Y -= 0.05;
                    std::string figures_final = Form("%s/appendix/HAD2_frac/%s",figures.c_str(),trigger_name.at(i).c_str());
                    gSystem->Exec(Form("mkdir -p %s",figures_final.c_str()));
                    c0->SaveAs(Form("%s/%s_%s.pdf",figures_final.c_str(),figures_name.at(side).c_str(),lb_range.at(ilb).c_str()));
                    h2->Reset();
                }
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

}

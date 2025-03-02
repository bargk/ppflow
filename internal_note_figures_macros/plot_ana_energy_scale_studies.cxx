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
int Size = 15;
void plot_ana_energy_scale_studies(){
    SetAtlasStyle();
    std::vector<int> lineColors = {kBlue, kRed, kOrange+3};
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TLegend *legend0;
    TCanvas *c0 = new TCanvas("c0");
    //plot the modules <ADC> for different triggers
    #ifdef average_adc
        for(int side =0; side <2; side++){
            for(int mod =1; mod<4; mod++){
                THStack *hs = new THStack("hs",";LB;<a.u>");
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
                Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
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
    TCanvas *c1 = new TCanvas("c1", "", 800,600);
    c1->Divide(2,1);
    #ifdef stability_error
        for(int side =0; side <2; side++){
            for(int mod =1; mod<4; mod++){
                THStack *hs = new THStack("hs",";LB;<a.u>/Reference value");
                THStack *hs1 = new THStack("hs1",";LB;<a.u>/Reference value");
                TFile *input_0 = new TFile(Form("%s/histograms.root",base.c_str()),"read");
                TFile *input_1 = new TFile(Form("%s/minbias/histograms.root",base.c_str()),"read"); 
                TFile *input_2 = new TFile(Form("%s/xor/histograms.root",base.c_str()),"read");
                
                TH2D *h2_00 = (TH2D*)input_0->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_AND",side,mod));
                TH2D *h2_01 = (TH2D*)input_1->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_minbias",side,mod));
                TH2D *h2_02 = (TH2D*)input_2->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_XOR",side,mod));
                TH1D *h2_00_profx = h2_00->ProfileX(); h2_00_profx->SetMarkerColor(lineColors.at(0));   h2_00_profx->SetLineColor(kBlack);
                TH1D *h2_01_profx = h2_01->ProfileX(); h2_01_profx->SetMarkerColor(lineColors.at(1));   h2_01_profx->SetLineColor(kBlack);
                TH1D *h2_02_profx = h2_02->ProfileX(); h2_02_profx->SetMarkerColor(lineColors.at(2));   h2_02_profx->SetLineColor(kBlack);
                //get the mean of each true block
                std::vector<int> bin_low_0;
                std::vector<int> bin_high_0;
                std::vector<int> bin_low_1;
                std::vector<int> bin_high_1;
                std::vector<int> bin_low_2;
                std::vector<int> bin_high_2;
                bin_low_0.push_back(h2_00_profx->FindBin(1816));    bin_high_0.push_back(h2_00_profx->FindBin(1889));
                bin_low_0.push_back(h2_00_profx->FindBin(1902));    bin_high_0.push_back(h2_00_profx->FindBin(1929));
                bin_low_0.push_back(h2_00_profx->FindBin(2812));    bin_high_0.push_back(h2_00_profx->FindBin(3111));
                bin_low_0.push_back(h2_00_profx->FindBin(3218));    bin_high_0.push_back(h2_00_profx->FindBin(3680));
                bin_low_1.push_back(h2_01_profx->FindBin(1816));    bin_high_1.push_back(h2_01_profx->FindBin(1889));
                bin_low_1.push_back(h2_01_profx->FindBin(1902));    bin_high_1.push_back(h2_01_profx->FindBin(1929));
                bin_low_1.push_back(h2_01_profx->FindBin(2812));    bin_high_1.push_back(h2_01_profx->FindBin(3111));
                bin_low_1.push_back(h2_01_profx->FindBin(3218));    bin_high_1.push_back(h2_01_profx->FindBin(3680));
                bin_low_2.push_back(h2_02_profx->FindBin(1816));    bin_high_2.push_back(h2_02_profx->FindBin(1889));
                bin_low_2.push_back(h2_02_profx->FindBin(1902));    bin_high_2.push_back(h2_02_profx->FindBin(1929));
                bin_low_2.push_back(h2_02_profx->FindBin(2812));    bin_high_2.push_back(h2_02_profx->FindBin(3111));
                bin_low_2.push_back(h2_02_profx->FindBin(3218));    bin_high_2.push_back(h2_02_profx->FindBin(3680));
                double mean_0 =0;
                double mean_1 =0;
                double mean_2 =0;
                for(int i=0; i<bin_low_0.size(); i++){
                    mean_0 += h2_00_profx->Integral(bin_low_0.at(i),bin_high_0.at(i),"width");
                    mean_1 += h2_01_profx->Integral(bin_low_1.at(i),bin_high_0.at(i),"width");
                    mean_2 += h2_02_profx->Integral(bin_low_2.at(i),bin_high_2.at(i),"width");
                    //cout<< h2_00_profx->Integral(bin_low_0.at(i),bin_high_0.at(i)) / (bin_high_0.at(i)-bin_low_0.at(i)) << endl;
                }
                mean_0 = mean_0/((bin_high_0.at(0)-bin_low_0.at(0)) + (bin_high_0.at(1)-bin_low_0.at(1)) + (bin_high_0.at(2)-bin_low_0.at(2)) + (bin_high_0.at(3)-bin_low_0.at(3)));
                mean_1 = mean_1/((bin_high_0.at(0)-bin_low_0.at(0)) + (bin_high_0.at(1)-bin_low_0.at(1)) + (bin_high_0.at(2)-bin_low_0.at(2)) + (bin_high_0.at(3)-bin_low_0.at(3)));
                mean_2 = mean_2/((bin_high_0.at(0)-bin_low_0.at(0)) + (bin_high_0.at(1)-bin_low_0.at(1)) + (bin_high_0.at(2)-bin_low_0.at(2)) + (bin_high_0.at(3)-bin_low_0.at(3)));
                cout << "average: " << mean_0 << endl;
                h2_00_profx->SetMarkerStyle(20);
                h2_01_profx->SetMarkerStyle(20);
                h2_02_profx->SetMarkerStyle(20);
                h2_00_profx->SetMarkerSize(0.5);
                h2_01_profx->SetMarkerSize(0.5);
                h2_02_profx->SetMarkerSize(0.5);
                //scale histograms by their mean
                h2_00_profx->Scale(1.0/mean_0);
                h2_01_profx->Scale(1.0/mean_1);
                h2_02_profx->Scale(1.0/mean_2);
                hs->Add(h2_00_profx);
                //hs->Add(h2_01_profx);
                hs->Add(h2_02_profx);
                hs->SetMinimum(0.7);
                hs->SetMaximum(1.3);
                c1->cd(1);
                hs->Draw("nostack");
                TLine *line = new TLine(1816, 1, 4755, 1); // Horizontal line at y=1
                line->SetLineColor(kBlack); // Set line color
                line->SetLineWidth(2); // Set line width
                line->SetLineStyle(2); // Dashed line
                line->Draw("same");
                hs->GetXaxis()->SetLimits(1816,4755);
                X = 0.2, Y = 0.85;
                Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                Common::myText2(X + 0.14, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                legend0 = new TLegend(0.6,0.72,0.88,0.88);
                legend0->SetBorderSize(0);
                legend0->AddEntry(h2_00_profx, "AND trigger", "lep");
                legend0->AddEntry(h2_01_profx, "minbias trigger", "lep");
                legend0->AddEntry(h2_02_profx, "XOR trigger", "lep");
                legend0->Draw();
                legend0->SetTextSize(0.033);
                c1->cd(2);
                hs1->Add(h2_01_profx);
                hs1->Draw("nostack");
                hs1->GetXaxis()->SetLimits(1816,4755);
                line->Draw("same");
                hs1->SetMinimum(0.7);
                hs1->SetMaximum(1.3);
                c1->SaveAs(Form("%s/stability_uncertainty_HAD%i_side%i.pdf",figures.c_str(),mod,side));
                bin_low_0.clear();
                bin_high_0.clear();
                bin_low_1.clear();
                bin_high_1.clear();
                bin_low_2.clear();
                bin_high_2.clear();
            }
        }
    #endif
    //------------------------------------------------------------------------------------------------------------------

    //plot unclibrated distributions of the sides for each trigger
    std::vector<string> triggers = {"/","/xor/","/minbias/"};
    std::vector<string> figures_name = {"sides_and.pdf","sides_xor.pdf","sides_minbias.pdf"};
    lineColors.clear();
    lineColors = {kBlack, kMagenta};
    #ifdef sides_uncalib
        for(int i =0; i<triggers.size(); i++){
            TFile *input = new TFile(Form("%s%shistograms.root",base.c_str(), triggers.at(i).c_str()),"read");
            THStack *hs = new THStack("hs",";ZDC signal [a.u];Events");
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
            Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
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
            THStack *hs = new THStack("hs",";ZDC Module [a.u];Events");
            TH2D *h2_1 = (TH2D*)input->Get(Form("h_bad%i_1",i))->Clone("HAD1");
            TH2D *h2_2 = (TH2D*)input->Get(Form("h_bad%i_2",i))->Clone("HAD2");
            TH2D *h2_3 = (TH2D*)input->Get(Form("h_bad%i_3",i))->Clone("HAD3");
            TH1D *h2_1_proj = h2_1->ProjectionY(); h2_1_proj->SetLineColor(kBlue+3); h2_1_proj->SetLineWidth(2);  
            TH1D *h2_2_proj = h2_2->ProjectionY(); h2_2_proj->SetLineColor(kBlue-3); h2_2_proj->SetLineWidth(2);  
            TH1D *h2_3_proj = h2_3->ProjectionY(); h2_3_proj->SetLineColor(kBlue+4); h2_3_proj->SetLineWidth(2);  
            c0->cd();
            hs->Add(h2_1_proj);
            hs->Add(h2_2_proj);
            hs->Add(h2_3_proj);
            hs->Draw("nostack");
            gPad->SetLogy();
            X = 0.6, Y = 0.85;
            Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
            Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
            Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
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
    std::vector<string> trigger_name = {"AND","XOR","minbias"};
    #ifdef HAD2_frac
        for(int i =0; i<triggers.size(); i++){
            TFile *input = new TFile(Form("%s%shistograms.root",base.c_str(), triggers.at(i).c_str()),"read");
            for(int side =0; side<2; side++){
                TH2D *h2 = (TH2D*)input->Get(Form("h%i_h2AmpRat_Amp",side))->Clone(Form("HAD2_frac_side%i",side));
                h2->GetXaxis()->SetTitle("ZDC signal [a.u]");
                h2->GetYaxis()->SetTitle("HAD2 Fraction");
                h2->GetZaxis()->SetTitle("Events");
                c0->cd();
                h2->Draw("colz");
                gPad->SetLogy(0);
                gPad->SetLogz();
                gPad->SetRightMargin(0.15);
                X = 0.6, Y = 0.85;
                Common::myText2(X, Y, 1, "ATLAS ", Size, 73);
                Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("%s trigger", trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("Mean y: %.2f", h2->GetMean(2)), Size, 43); Y -= 0.05;
                Common::myText2(X, Y, 1, Form("StdDev y: %.2f", h2->GetStdDev(2)), Size, 43); Y -= 0.05;
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
                    TH2D *h2 = (TH2D*)input->Get(Form("h%i_h2AmpRat_Amp_LB_%s",side,lb_range.at(ilb).c_str()));
                    h2->GetXaxis()->SetTitle("ZDC signal [a.u]");
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
                    Common::myText2(X + 0.08, Y, 1, Common::Internal, Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, "#it{pp 22}, #sqrt{#it{s}} = 13.6 TeV", Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, Form("%s trigger", trigger_name.at(i).c_str()), Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, Form("Mean y: %.3f", h2->GetMean(2)), Size, 43); Y -= 0.05;
                    Common::myText2(X, Y, 1, Form("StdDev y: %.3f", h2->GetStdDev(2)), Size, 43); Y -= 0.05;
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

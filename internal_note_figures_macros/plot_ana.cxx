#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana";
TH1D *h_zdc[3];
TH1D *h1_profile[3];
TH2D *h2;
TLegend *legend0;
float X, Y;
void plot_ana(){
    int size = 20;
    SetAtlasStyle();
    TCanvas *c0 = new TCanvas("c0");
    std::vector<int> lineColors = {kBlue, kRed, kOrange+3};
    //plot events in each bin of EE for all the 3 triggers
    TFile *input = new TFile(Form("%s/histograms.root",base.c_str()));
    TFile *input1 = new TFile(Form("%s/minbias/histograms.root",base.c_str()));
    TFile *input2 = new TFile(Form("%s/xor/histograms.root",base.c_str()));
    h_zdc[0] = (TH1D*)input->Get("heff_after_cut");     h_zdc[0]->SetLineColor(lineColors.at(0));     h_zdc[0]->SetMarkerColor(lineColors.at(0));
    h_zdc[1] = (TH1D*)input1->Get("heff_after_cut");    h_zdc[1]->SetLineColor(lineColors.at(1));      h_zdc[1]->SetMarkerColor(lineColors.at(1));
    h_zdc[2] = (TH1D*)input2->Get("heff_after_cut");    h_zdc[2]->SetLineColor(lineColors.at(2));   h_zdc[2]->SetMarkerColor(lineColors.at(2));
    double maxY = 0; // Variable to store the global maximum
    THStack *hs = new THStack("hs",";Effective Energy [TeV];Events");
    
    hs->Add(h_zdc[0]);
    hs->Add(h_zdc[1]);
    hs->Add(h_zdc[2]);
    hs->Draw("nostack");
    hs->SetMaximum(10e8);
    legend0 = new TLegend(0.2,0.3,0.4,0.38);
    legend0->SetBorderSize(0);
    legend0->AddEntry(h_zdc[0],"AND trigger","l");
    legend0->AddEntry(h_zdc[1],"minbias triggers","l");
    legend0->AddEntry(h_zdc[2],"XOR triggers","l");
    legend0->Draw();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.25);
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp}22 #sqrt{#it{s}} = 13.6 TeV", size, 43);
    c0->SaveAs(Form("%s/nevents.pdf",figures.c_str()));
    //----------------------------------------------------------------------------------

    //plot 2D correlation between nrec to effective energy
    std::vector<string> triggers = {"/","/xor/","/minbias/"};
    std::vector<string> trigger_name = {"AND","XOR","minbias"};
    lineColors = {kBlue,kOrange+3,kRed};
    THStack *hs1 = new THStack("hs1",";E_{Eff} [TeV]; #LT N_{ch}^{rec}#GT");
    for(int i=0; i<triggers.size(); i++){
        input = new TFile(Form("%s%shistograms.root",base.c_str(),triggers.at(i).c_str()));
        c0->cd();
        h2 = (TH2D*)input->Get("hNtrkEff")->Clone(Form("clone_%s",trigger_name.at(i).c_str()));
        h2->GetXaxis()->SetRange(1,16);
        h2->GetXaxis()->SetTitle("E_{Eff} [TeV]");
        h2->GetYaxis()->SetTitle("N_{ch}^{rec}");
        if(i == 1){
            h2->GetXaxis()->SetRange(9,18);
        }
        if(i == 2){
            h2->GetXaxis()->SetRange(1,20);
        }
        gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
        h2->Draw("colz");
        gPad->SetLogz();
        gPad->SetLogy(0);
        h1_profile[i] = h2->ProfileX("hprofile");
        h1_profile[i]->SetMarkerStyle(20);
        h1_profile[i]->SetMarkerSize(1);
        h1_profile[i]->SetMarkerColor(lineColors.at(i));
        // h1_profile->Draw("same");
        X=0.20,Y=0.9;
        Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp}22 #sqrt{#it{s}} = 13.6 TeV", size, 43);Y -= 0.05;
        Common::myText2(X, Y, 1, Form("%s trigger", trigger_name.at(i).c_str()), size, 43); 
        c0->SaveAs(Form("%s/nch_eff_%s.pdf",figures.c_str(), trigger_name.at(i).c_str()));
        c0->Clear();
        hs1->Add(h1_profile[i]);
    }
    hs1->Draw("nostack");
    hs1->SetMinimum(14);
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp}22 #sqrt{#it{s}} = 13.6 TeV", size, 43);
    legend0 = new TLegend(0.2,0.7,0.4,0.8);
    legend0->SetBorderSize(0);
    legend0->AddEntry(h1_profile[0],"AND trigger","lep");
    legend0->AddEntry(h1_profile[1],"XOR trigger","lep");
    legend0->AddEntry(h1_profile[2],"minbias trigger","lep");
    legend0->Draw();
    c0->SaveAs(Form("%s/avg_nch.pdf",figures.c_str()));

}
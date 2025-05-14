#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
void plot_modules(){
    SetAtlasStyle();
    std::string directory = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/modules");
    std::string data = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp");
    gSystem->Exec(Form("mkdir -p %s",directory.c_str()));
    
    TFile *input_0 = new TFile(Form("%s/histograms.root",data.c_str()),"read");
    TFile *input_1 = new TFile(Form("%s/histograms.1.root",data.c_str()),"read"); //first cut based fractional HAD2 from zdc sum
    TFile *input_2 = new TFile(Form("%s/histograms.2.root",data.c_str()),"read"); //second cut based fractional HAD3 from zdc sum and the first cut

    TH2D *h2_01 = (TH2D*)input_0->Get("h0_1")->Clone("had1_AND");
    TH2D *h2_02 = (TH2D*)input_1->Get("h0_1")->Clone("had1_minbias");
    TH2D *h2_03 = (TH2D*)input_2->Get("h0_1")->Clone("had1_XOR");
    TH1D *h2_01_profx = h2_01->ProfileX(); h2_01_profx->SetMarkerColor(kGray);  h2_01_profx->SetLineColor(kGray);
    TH1D *h2_02_profx = h2_02->ProfileX(); h2_02_profx->SetMarkerColor(kGray+1);h2_02_profx->SetMarkerColor(kGray+2);
    TH1D *h2_03_profx = h2_03->ProfileX(); h2_03_profx->SetMarkerColor(kGray+2);h2_03_profx->SetMarkerColor(kGray+3);

    TCanvas *c0 = new TCanvas("c0");
    h2_01_profx->Draw();
    h2_02_profx->Draw("same");
    h2_03_profx->Draw("same");
    c0->SaveAs(Form("%s/HAD1.png",directory.c_str()));
    //TODO do it for all the modules

    //plot the modules <ADC> for different triggers
    for(int side =0; side <2; side++){
        for(int mod =1; mod<4; mod++){

            THStack *hs = new THStack("hs",";LB;<a.u>");
            TFile *input_and = new TFile(Form("%s/histograms.2.root",data.c_str()),"read");
            TFile *input_minbias = new TFile(Form("%s/minbias/histograms.2.root",data.c_str()),"read");
            TFile *input_xor = new TFile(Form("%s/xorE2/histograms.2.root",data.c_str()),"read");
            h2_01 = (TH2D*)input_and->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_AND",side,mod));
            h2_02 = (TH2D*)input_minbias->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_minbias",side,mod));
            h2_03 = (TH2D*)input_xor->Get(Form("h%i_%i",side,mod))->Clone(Form("had%i_side%i_XOR",side,mod));
            TH1D *h2_01_profx = h2_01->ProfileX();
            TH1D *h2_02_profx = h2_02->ProfileX();
            TH1D *h2_03_profx = h2_03->ProfileX();
            h2_01_profx->SetMarkerColor(kBlue);  h2_01_profx->SetLineColor(kBlue);
            h2_02_profx->SetMarkerColor(kRed);  h2_02_profx->SetLineColor(kRed);
            h2_03_profx->SetMarkerColor(kOrange+3);  h2_03_profx->SetLineColor(kOrange+3);
            hs->Add(h2_01_profx);
            hs->Add(h2_02_profx);
            hs->Add(h2_03_profx);
            hs->SetMinimum(0.7*h2_02_profx->GetBinContent(1500));
            hs->SetMaximum(1.3*h2_01_profx->GetBinContent(1500));
            hs->Draw("nostack");
            c0->SaveAs(Form("%s/HAD%i_side%i.png",directory.c_str(),mod,side));
        }
    }
    


}
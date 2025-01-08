
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base1 = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/1.5sigma";
std::string  base2 = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/2.0sigma";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/systematics";

TH1D *h[4];
void plot_ana_systematics(){
    SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TCanvas *c0 = new TCanvas("c0","",3000,2500);
    TFile *input1 =  new TFile(Form("%s/TemplateFits_vnn.root",base1.c_str()));
    TFile *input2 =  new TFile(Form("%s/TemplateFits_vnn.root",base2.c_str()));
    //TFile *input3 =  new TFile(Form("%s/sameSide/TemplateFits_vnn.root",base1.c_str()));
    TFile *input4 =  new TFile(Form("%s/sameSide/TemplateFits_vnn.root",base2.c_str()));
    for(int i=0; i<4; i++){
        h[i] = new TH1D(Form("h_v2_zdc%i",i), ";E_{Eff} [TeV];v_{2}(p_{T}^{b})", 8, Bins::CENT_LO[23], Bins::CENT_HI[30]);
        h[i]->SetLineColor(kBlue+2*i);
        h[i]->SetMarkerColor(kBlue+2*i);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(3);
    }
    std::pair<float, float> vnn_zdc1;
    std::pair<float, float> vnn_zdc2;
    int pericent =0;
    for(int icent =0; icent< 8; icent++){
                vnn_zdc1=      Bins::GetVnPtb(icent +23,13,5,5,2,1,2,pericent,0,input1,input1); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                h[0]->SetBinContent(icent+1, vnn_zdc1.first);
                h[0]->SetBinError(icent+1, vnn_zdc1.second);
                vnn_zdc2=      Bins::GetVnPtb(icent +23,13,5,5,2,1,2,pericent,0,input2,input2); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                h[1]->SetBinContent(icent+1, vnn_zdc2.first);
                h[1]->SetBinError(icent+1, vnn_zdc2.second);
                // vnn_zdc1=      Bins::GetVnPtb(icent +23,13,5,5,2,1,2,pericent,0,input3,input3); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                // h[2]->SetBinContent(icent+1, vnn_zdc1.first);
                // h[2]->SetBinError(icent+1, vnn_zdc1.second);
                vnn_zdc2=      Bins::GetVnPtb(icent +23,13,5,5,2,1,2,pericent,0,input4,input4); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                h[3]->SetBinContent(icent+1, vnn_zdc2.first);
                h[3]->SetBinError(icent+1, vnn_zdc2.second);
        }
    THStack *hs_average = new THStack("hs",";E_{Eff} [TeV]; v_{2}(p_{T}^{b})");
    hs_average->Add(h[0]);
    hs_average->Add(h[1]);
    hs_average->Add(h[3]);
    hs_average->SetMaximum(0.15);
    hs_average->SetMinimum(0);
    hs_average->Draw("nostack;E1");
    TLegend *legend0 = new TLegend(0.5,0.8,0.8,0.9);
    legend0->AddEntry(h[0],"1.5 #sigma opposite ","lep");
    legend0->AddEntry(h[1],"2.0 #sigma opposite ","lep");
    //legend0->AddEntry(h[2],"1.5 #sigma same side ","lep");
    legend0->AddEntry(h[3],"2.0 #sigma same side ","lep");
    legend0->Draw();

    c0->SaveAs(Form("%s/v2_sys.pdf",figures.c_str()));
}
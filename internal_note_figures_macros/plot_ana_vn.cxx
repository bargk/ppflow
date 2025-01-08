//ploting v2 for analaysis section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/2.0sigma/sameSide";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/vn";

void plot_ana_vn(){
    SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TH1D *h_v2_zdc = new TH1D("h_v2_zdc", ";E_{Eff} [TeV];v_{2}(p_{T}^{b})", 8, Bins::CENT_LO[23], Bins::CENT_HI[30]); 
    TH1D *h_v2_xor = new TH1D("h_v2_xor", ";E_{Eff} [TeV];", 4, Bins::CENT_LO[27], Bins::CENT_HI[30]); 
    TH1D *h_v2_minbias_1 = new TH1D("h_v2_minbias_1", ";E_{Eff} [TeV];", 1, Bins::CENT_LO[22], Bins::CENT_HI[22]); 
    TH1D *h_v2_minbias_2 = new TH1D("h_v2_minbias_2", ";E_{Eff} [TeV];", 1, Bins::CENT_LO[19], Bins::CENT_HI[19]);
    h_v2_zdc->SetMarkerColor(kBlue); h_v2_zdc->SetLineColor(kBlue);
    h_v2_minbias_1->SetMarkerColor(kRed); h_v2_minbias_1->SetLineColor(kRed);
    h_v2_minbias_2->SetMarkerColor(kRed); h_v2_minbias_2->SetLineColor(kRed);
    h_v2_xor->SetMarkerColor(kOrange); h_v2_xor->SetLineColor(kOrange);

    TFile *input =  new TFile(Form("%s/TemplateFits_vnn.root",base.c_str()));
    TFile *input1 = new TFile(Form("%s/minbias/TemplateFits_vnn.root",base.c_str()));
    TFile *input2 = new TFile(Form("%s/xorE2/TemplateFits_vnn.root",base.c_str()));

    int pericent =0;
    std::pair<float, float> vnn_zdc;
    std::pair<float, float> vnn_minbias1;
    std::pair<float, float> vnn_minbias2;
    for(int icent =0; icent< 8; icent++){
                vnn_zdc=      Bins::GetVnPtb(icent +23,13,5,5,2,1,2,pericent,0,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                h_v2_zdc->SetBinContent(icent+1, vnn_zdc.first);
                h_v2_zdc->SetBinError(icent+1, vnn_zdc.second);
        }
    for(int icent =0; icent< 4; icent++){
        vnn_zdc=      Bins::GetVnPtb(icent +27,13,5,5,2,1,2,pericent,0,input2,input2); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
        h_v2_xor->SetBinContent(icent+1, vnn_zdc.first);
        h_v2_xor->SetBinError(icent+1, vnn_zdc.second);
    }
    vnn_minbias1=      Bins::GetVnPtb(22,13,5,5,2,1,2,pericent,0,input1,input1); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
    h_v2_minbias_1->SetBinContent(1, vnn_minbias1.first);
    h_v2_minbias_1->SetBinError(1, vnn_minbias1.second);
    vnn_minbias2=      Bins::GetVnPtb(19,13,5,5,2,1,2,pericent,0,input1,input1); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
    h_v2_minbias_2->SetBinContent(1, vnn_minbias2.first);
    h_v2_minbias_2->SetBinError(1, vnn_minbias2.second);

    TCanvas *c0 = new TCanvas("c0","",3000,2500);
    THStack *hs_average = new THStack("hs",";E_{Eff} [TeV]; v_{2}(p_{T}^{b})");
    hs_average->Add(h_v2_zdc);
    hs_average->Add(h_v2_xor);
    hs_average->Add(h_v2_minbias_1);
    hs_average->Add(h_v2_minbias_2);
    hs_average->SetMaximum(0.15);
    hs_average->SetMinimum(0);
    hs_average->Draw("nostack;E1");
    float X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
    TLegend *legend0 = new TLegend(0.5,0.8,0.8,0.9);
    legend0->AddEntry(h_v2_zdc,"ZDC AND ","lep");
    legend0->AddEntry(h_v2_minbias_1,"Minimum bias","lep");
    legend0->AddEntry(h_v2_xor,"ZDC XOR","lep");
    legend0->Draw();
    c0->SaveAs(Form("%s/v2.pdf",figures.c_str())); 
}
//ploting v2 for analaysis section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/vn";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/vn";

void plot_ana_vn(){
    SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    //empty histogram just for limit of canvas
    TH1D *h_limits = new TH1D("h_limits","",1,6.1,13.61);
    h_limits->SetBinContent(1,-10);

    // const int nbins_and = 7;
    // Double_t edges_and[nbins_and + 1] = {3.6,6.1,7.1,8.1,9.1,10.1,11.1,11.6};
    const int nbins_and = 10;
    Double_t edges_and[nbins_and + 1] = {6.1,7.1,7.6,8.1,8.6,9.1,9.6,10.1,10.6,11.1,11.6};
    TH1D *hv2_and = new TH1D("hv2_and","",nbins_and,edges_and);
    const int nbins_xor = 9;
    Double_t edges_xor[nbins_xor + 1] = {7.1,8.6,9.1,9.6,10.1,10.6,11.1,11.6,12.1,12.6};
    TH1D *hv2_xor = new TH1D("hv2_xor","",nbins_xor,edges_xor);
    const int nbins_minbias = 2;
    Double_t edges_minbias[nbins_minbias + 1] = {11.6,13.1,13.61};
    TH1D *hv2_minbias = new TH1D("hv2_minbias","",nbins_minbias,edges_minbias);

    hv2_and->SetMarkerColor(kBlue); hv2_and->SetMarkerStyle(20);
    hv2_xor->SetMarkerColor(kOrange+3); hv2_xor->SetMarkerStyle(20);
    hv2_minbias->SetMarkerColor(kRed); hv2_minbias->SetMarkerStyle(20);

    TCanvas *c0 = new TCanvas("c0");
    TFile *input =  new TFile(Form("%s/TemplateFits_vnn.root",base.c_str()));
    TFile *input1 = new TFile(Form("%s/minbias/TemplateFits_vnn.root",base.c_str()));
    TFile *input2 = new TFile(Form("%s/xor/TemplateFits_vnn.root",base.c_str()));

    int pericent = 22;
    int peritrk = 5;
    int itrk =4;
    int ipt1 = 5;
    int ipt2 = 5;
    //std::pair<float, float> vnn_zdc;
    std::pair<float, float> vnn_zdc_bigbin;
    std::pair<float, float> vnn_minbias1;
    std::pair<float, float> vnn_minbias2;
    for(int icent =0; icent< nbins_and; icent++){
                int bin_idx = Bins::GetCentIndex(edges_and[icent],edges_and[icent +1]);
                std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,ipt1,ipt2,2,1,2,pericent,peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                hv2_and->SetBinContent(icent + 1, vnn_zdc.first);
                hv2_and->SetBinError(icent + 1, vnn_zdc.second);
    }

    for(int icent =0; icent< nbins_xor; icent++){
            int bin_idx = Bins::GetCentIndex(edges_xor[icent],edges_xor[icent +1]);
            std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,ipt1,ipt2,2,1,2,pericent,peritrk,input2,input2); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
            hv2_xor->SetBinContent(icent + 1, vnn_zdc.first);
            hv2_xor->SetBinError(icent + 1, vnn_zdc.second);
    }

    for(int icent =0; icent< nbins_minbias; icent++){
        int bin_idx = Bins::GetCentIndex(edges_minbias[icent],edges_minbias[icent +1]);
        std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,ipt1,ipt2,2,1,2,pericent,peritrk,input1,input1); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
        hv2_minbias->SetBinContent(icent + 1, vnn_zdc.first);
        hv2_minbias->SetBinError(icent + 1, vnn_zdc.second);
    }

    THStack *hs_average = new THStack("hs",";E_{Eff} [TeV]; v_{2}");
    hs_average->Add(h_limits);
    hs_average->Add(hv2_and);
    hs_average->Add(hv2_xor);
    hs_average->Add(hv2_minbias);
    hs_average->SetMaximum(0.09);
    hs_average->SetMinimum(0.045);
    hs_average->Draw("nostack;E1");
    //hs_average->GetXaxis()->SetLimits(3.6,13.6);
    float X=0.2,Y=0.88;
    int size =17;
    Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, Form("%s",Bins::label_ptab(ipt1,ipt2).c_str()), size, 43); Y=Y-0.05;
    Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta (Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43);
    TLegend *legend0 = new TLegend(0.45,0.7,0.88,0.88);
    legend0->SetBorderSize(0);
    legend0->AddEntry(hv2_and,"ZDC AND ","lep");
    legend0->AddEntry(hv2_minbias,"Minimum bias","lep");
    legend0->AddEntry(hv2_xor,"ZDC XOR","lep");
    legend0->Draw();
    c0->SaveAs(Form("%s/v2.pdf",figures.c_str()));
}
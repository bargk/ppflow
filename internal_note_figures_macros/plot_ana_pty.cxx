//ploting 2d pc for analaysis section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/1.5sigma";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/pty";
TH1D *h[31];
std::vector<int> bins_cent = {23,24,25,26,27,28,29,30,22,19}; 

void plot_ana_pty(){
    SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TFile *input = new TFile(Form("%s/PTY1D.root",base.c_str()));
    TFile * fileinput = new TFile(Form("%s/histograms.root",base.c_str()));
    TH2D *h2_19 = (TH2D*)fileinput->Get("hNtrkEff")->Clone();
    TH1D *h1 = h2_19->ProfileX("h1");

    //ploting manually the 2 higher bins
    h[19] = (TH1D*)input->Get(Form("PTY_cent19_trk13_pta5_ptb05_ch2_deta01"));
    TCanvas *c1 = new TCanvas("c1","",3000,2500);
    h[19]->SetMarkerSize(5);
    h[19]->SetMarkerStyle(20);
    h[19]->Draw();
    float X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X , Y, 1, "#it{pp} 22  #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    Common::myText2(X , Y , 1, Bins::label_cent(19)     , 70, 43);Y=Y-0.05;
    Common::myText2(X , Y , 1, Bins::label_ptab(5, 5), 70, 43);Y=Y-0.05;
    double nrec = h1->GetBinContent(19+1);
    Common::myText2(X , Y , 1, Form("< N_{ch}^{rec} > = %.2f",nrec), 70, 43);Y=Y-0.05;
    c1->SaveAs(Form("%s/PTY_cent19_trk13_pta5_ptb05_ch2_deta01.pdf",figures.c_str()));
    c1->Clear();

    h[22] = (TH1D*)input->Get(Form("PTY_cent22_trk13_pta5_ptb05_ch2_deta01"));
    h[22]->SetMarkerSize(5);
    h[22]->SetMarkerStyle(20);
    h[22]->Draw();
    TH2D *h2_22 = (TH2D*)fileinput->Get("hNtrkEff")->Clone();
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X , Y, 1, "#it{pp} 22  #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    Common::myText2(X , Y , 1, Bins::label_cent(22)     , 70, 43);Y=Y-0.05;
    Common::myText2(X , Y , 1, Bins::label_ptab(5, 5), 70, 43);Y=Y-0.05;
    nrec = h2_22->ProjectionY("hprjc",16,18)->GetMean();
    Common::myText2(X , Y , 1, Form("< N_{ch}^{rec} > = %.2f",nrec), 70, 43);Y=Y-0.05;
    c1->SaveAs(Form("%s/PTY_cent22_trk13_pta5_ptb05_ch2_deta01.pdf",figures.c_str()));
    c1->Clear();

    TH2D *h2 = (TH2D*)fileinput->Get("hNtrkEff")->Clone();
    h2->RebinX(2);
    h1 = h2->ProfileX("hprofile");
    for(int icent =23; icent<31; icent++){
        //h[icent -23] = (TH1D*)input->Get(Form("h_central_cent%.2i_pericent00_peritrk00_trk13_pta5_ptb05_ch2_deta01",icent)); 
        h[icent -23] = (TH1D*)input->Get(Form("PTY_cent%.2i_trk13_pta5_ptb05_ch2_deta01",icent)); 
        TCanvas *c0 = new TCanvas("c0","",3000,2500);
        h[icent -23]->SetMarkerSize(5);
        h[icent -23]->SetMarkerStyle(20);
        h[icent -23]->Draw();

        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X , Y, 1, "#it{pp} 22  #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        Common::myText2(X , Y , 1, Bins::label_cent(icent)     , 70, 43);Y=Y-0.05;
        Common::myText2(X , Y , 1, Bins::label_ptab(5, 5), 70, 43);Y=Y-0.05;

        nrec = h1->GetBinContent(icent -22);
        Common::myText2(X , Y , 1, Form("< N_{ch}^{rec} > = %.2f",nrec), 70, 43);Y=Y-0.05;
        c0->SaveAs(Form("%s/PTY_cent%.2i_trk13_pta5_ptb05_ch2_deta01.pdf",figures.c_str(),icent));
        delete c0;
    }
    //plot peripheral bin
    h[0] = (TH1D*)input->Get(Form("PTY_cent00_trk00_pta5_ptb05_ch2_deta01"));
    c1->cd(1);
    h[0]->SetMarkerSize(5);
    h[0]->SetMarkerStyle(20);
    h[0]->Draw();
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X , Y, 1, "#it{pp} 22  #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    c1->SaveAs(Form("%s/PTY_cent00_trk00_pta5_ptb05_ch2_deta01.pdf",figures.c_str()));

    
}
//ploting 2d pc for analaysis section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide/xor";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/pty";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/pty/xor";

std::vector<int> bins_cent = Bins::CentBins();
float X,Y;
void plot_ana_pty(){
    SetAtlasStyle();
    int size = 20;
    // std::vector<string> triggers = {"/","/xor/","/minbias/"};
    // std::vector<string> trigger_name = {"AND","XOR","minbias"};
    TCanvas *c0 = new TCanvas("c0");
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TFile *input = new TFile(Form("%s/PTY1D.root",base.c_str()));
    TFile * fileinput = new TFile(Form("%s/histograms.root",base.c_str()));
    for(auto &icent : bins_cent){
        TH2D *h2 = (TH2D*)fileinput->Get("hNtrkEff")->Clone();
        int bin_low = h2->GetXaxis()->FindBin(Bins::CENT_LO[icent]);
        int bin_high = h2->GetXaxis()->FindBin(Bins::CENT_HI[icent]);
        TH1D *h_proj = h2->ProjectionY(Form("%i",icent),bin_low,bin_high);
        float nrec = h_proj->GetMean();
        TH1D *h_pty =(TH1D*)input->Get(Form("PTY_cent%02i_trk04_pta5_ptb05_ch2_deta01",icent));
        h_pty->SetMarkerStyle(20);
        h_pty->GetYaxis()->SetTitle("Y(#Delta#phi)");
        h_pty->GetXaxis()->SetTitle("#Delta#phi");
        h_pty->Draw();
        float X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
        Common::myText2(X , Y, 1, "#it{pp} 22, #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
        Common::myText2(X , Y , 1, "XOR trigger", size, 43);Y=Y-0.06;
        Common::myText2(X , Y , 1, Bins::label_cent(icent)     , size, 43);Y=Y-0.05;
        Common::myText2(X , Y , 1, Bins::label_ptab(5, 5), size, 43);Y=Y-0.06;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta (Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43);Y=Y-0.06;
        Common::myText2(X , Y , 1, Form("#LT N_{ch}^{rec}#GT = %.2f",nrec), size, 43);Y=Y-0.06;
        c0->SaveAs(Form("%s/PTY_cent%02i_trk04_pta5_ptb05_ch2_deta01.pdf",figures.c_str(),icent));
    }
    //--------------------------------------------------------------------------------------------------------------------------

    //ploting peripheral bin
    bins_cent = {20,21,22,23};
    for(auto &icent : bins_cent){
        TH1D *h_pty =(TH1D*)input->Get(Form("PTY_cent%02i_trk05_pta5_ptb05_ch2_deta01",icent));
        h_pty->SetMarkerStyle(20);
        h_pty->GetYaxis()->SetTitle("Y(#Delta#phi)");
        h_pty->GetXaxis()->SetTitle("#Delta#phi");
        h_pty->Draw();
        float X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
        Common::myText2(X , Y, 1, "#it{pp} 22, #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
        Common::myText2(X , Y , 1, "XOR trigger", size, 43);Y=Y-0.06;
        Common::myText2(X , Y , 1, Bins::label_cent(icent)     , size, 43);Y=Y-0.05;
        Common::myText2(X , Y , 1, Bins::label_trk(5), size, 43);Y=Y-0.06;
        Common::myText2(X , Y , 1, Bins::label_ptab(5, 5), size, 43);Y=Y-0.06;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta (Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43);Y=Y-0.06;
        c0->SaveAs(Form("%s/peri_PTY_cent%02i_trk05_pta5_ptb05_ch2_deta01.pdf",figures.c_str(),icent));
    }
}
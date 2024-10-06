#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "bins.h"
#include "common.C"
TH1D *h_zdc[Bins::NTRK];
TH1D *h_minbias[Bins::NTRK];

/*-----------------------------------------------------------------------------
 *  Plot some basic kinematics to compare between minbias & zdc triggers
 *-----------------------------------------------------------------------------*/

 void plots_monitoring(){
    SetAtlasStyle();
    std::string outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    TFile *input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/histograms.root");
    TFile *input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/histograms.root");

    TCanvas* c0 = new TCanvas("c0","",3000,3000);
    TLegend *legend0;
    c0->Divide(4,3);
    //plot eta comparison
    for(int itrk=0; itrk<Bins::NTRK-1; itrk++){
        c0->cd(itrk+1);
        h_zdc[itrk] = (TH1D*)input_zdc->Get(Form("h_eta_itrk%.2i",itrk));
        h_minbias[itrk] = (TH1D*)input_minbias->Get(Form("h_eta_itrk%.2i",itrk));
        h_zdc[itrk]->GetYaxis()->SetTitle("Normalized counts");
        h_minbias[itrk]->GetYaxis()->SetTitle("Normalized counts");
        h_zdc[itrk]->SetLineColor(kBlue);
        h_minbias[itrk]->SetLineColor(kRed);

        //normalizing
        h_zdc[itrk]->Scale(1/h_zdc[itrk]->GetEntries());
        h_minbias[itrk]->Scale(1/h_minbias[itrk]->GetEntries());
        h_zdc[itrk]->Draw("HIST");
        h_minbias[itrk]->Draw("SAME;HIST");
        float X=0.50,Y=0.5;
        Common::myText2(X, Y, 1, Bins::label_trk (itrk) , 33, 43);
        if(itrk ==0){
            legend0 = new TLegend(0.3,0.2,0.52,0.3);
            legend0->AddEntry(h_zdc[itrk],"ZDC triggers","l");
            legend0->AddEntry(h_minbias[itrk],"minbias triggers","l");
            gStyle->SetLegendTextSize(0.02);
            legend0->Draw();
        }
    }
    c0->SaveAs(Form("%s/eta.pdf",outputpath.c_str()));
    c0->Clear();
    c0->Divide(4,3);
    //plot pt
    for(int itrk=0; itrk<Bins::NTRK-1; itrk++){
        // h_zdc[itrk]->Clear();
        // h_minbias[itrk]->Clear();
        c0->cd(itrk+1);
        h_zdc[itrk] = (TH1D*)input_zdc->Get(Form("h_pt_itrk%.2i",itrk));
        h_minbias[itrk] = (TH1D*)input_minbias->Get(Form("h_pt_itrk%.2i",itrk));
        h_zdc[itrk]->GetYaxis()->SetTitle("Normalized counts");
        h_minbias[itrk]->GetYaxis()->SetTitle("Normalized counts");
        h_zdc[itrk]->SetLineColor(kBlue);
        h_minbias[itrk]->SetLineColor(kRed);

        //normalizing
        h_zdc[itrk]->Scale(1/h_zdc[itrk]->GetEntries());
        h_minbias[itrk]->Scale(1/h_minbias[itrk]->GetEntries());
        h_zdc[itrk]->Draw("HIST");
        h_minbias[itrk]->Draw("SAME;HIST");
        float X=0.5,Y=0.5;
        Common::myText2(X, Y, 1, Bins::label_trk (itrk) , 33, 43);
        if(itrk ==0){
            legend0 = new TLegend(0.3,0.2,0.52,0.3);
            legend0->AddEntry(h_zdc[itrk],"ZDC triggers","l");
            legend0->AddEntry(h_minbias[itrk],"minbias triggers","l");
            gStyle->SetLegendTextSize(0.02);
            legend0->Draw();
        }
    }
    c0->SaveAs(Form("%s/pt.pdf",outputpath.c_str()));
    delete c0;
 }
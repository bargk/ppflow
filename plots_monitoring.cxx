#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "bins.h"
#include "common.C"
TH1D *h_zdc[Bins::NCENT];
TH1D *h_minbias[Bins::NCENT];
#define eta
#define pt
#define PTY

/*-----------------------------------------------------------------------------
 *  Plot some basic phyisics to compare between minbias & zdc triggers
 *-----------------------------------------------------------------------------*/

 void plots_monitoring(){
    //SetAtlasStyle();
    std::string outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    TFile *input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/histograms.root");
    TFile *input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/histograms.root");

    TCanvas* c0 = new TCanvas("c0","",3000,3000);
    TLegend *legend0;
    c0->Divide(4,3);
    #ifdef eta
        //plot eta zdc
        for(int icent=0; icent<12; icent++){
            c0->cd(icent+1);
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_eta_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_eta_icent%.2i",icent));
            h_zdc[icent]->GetYaxis()->SetTitle("Normalized counts");
            h_minbias[icent]->GetYaxis()->SetTitle("Normalized counts");
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);

            h_zdc[icent]->Draw("HIST");
            h_minbias[icent]->Draw("SAME;HIST");
            float X=0.50,Y=0.5;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 33, 43);
            if(icent ==0){
                legend0 = new TLegend(0.3,0.2,0.52,0.3);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.03);
                legend0->Draw();
            }
        }
        c0->SaveAs(Form("%s/eta.png",outputpath.c_str()));
        c0->Clear();
    #endif
    
    #ifdef pt
        //plot pt
        c0->Divide(4,3);
        for(int icent=0; icent<12; icent++){
            // h_zdc[icent]->Clear();
            // h_minbias[icent]->Clear();
            c0->cd(icent+1);
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_pt_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_pt_icent%.2i",icent));
            h_zdc[icent]->GetYaxis()->SetTitle("Normalized counts");
            h_minbias[icent]->GetYaxis()->SetTitle("Normalized counts");
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);

            h_zdc[icent]->Draw("HIST");
            h_minbias[icent]->Draw("SAME;HIST");
            gPad->SetLogy();
            float X=0.5,Y=0.5;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 33, 43);
            if(icent ==0){
                legend0 = new TLegend(0.3,0.7,0.52,0.8);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.03);
                legend0->Draw();
            }
        }
        c0->SaveAs(Form("%s/pt.png",outputpath.c_str()));
        #endif

    #ifdef PTY
        //plot  PTY for multiplicity
        input_zdc->Clear();
        input_minbias->Clear();
        c0->Clear();
        c0->Divide(3,2);
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/TemplateFits.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/TemplateFits.root");
        for(int itrk=3; itrk<10; itrk++){
            c0->cd(itrk-2);  
            h_zdc[itrk] = (TH1D*)input_zdc->Get(Form("h_central_cent20_pericent20_peritrk13_trk%.2i_pta5_ptb05_ch2_deta01",itrk)); 
            h_minbias[itrk] = (TH1D*)input_minbias->Get(Form("h_central_cent20_pericent20_peritrk13_trk%.2i_pta5_ptb05_ch2_deta01",itrk));

            h_zdc[itrk]->SetLineColor(kBlue);
            h_minbias[itrk]->SetLineColor(kRed);
            h_minbias[itrk]->Draw(); 
            h_zdc[itrk]->Draw("SAME");

            float X=0.2,Y=0.85;
            Common::myText2(X, Y, 1, Bins::label_trk (itrk) , 33, 43);
            if(itrk ==3){
                legend0 = new TLegend(0.2,0.7,0.52,0.8);
                legend0->AddEntry(h_zdc[itrk],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[itrk],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.03);
                legend0->Draw();
            }
        }
        c0->SaveAs(Form("%s/PTY.png",outputpath.c_str()));
        
        //plot  PTY for effective energy
        input_zdc->Clear();
        input_minbias->Clear();
        c0->Clear();
        c0->Divide(3,3);
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/TemplateFits.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/TemplateFits.root");
        for(int icent=0; icent<10; icent++){
            c0->cd(icent+1);  
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_central_cent%.2i_pericent00_peritrk13_trk14_pta5_ptb05_ch2_deta01",icent)); 
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_central_cent%.2i_pericent00_peritrk13_trk14_pta5_ptb05_ch2_deta01",icent));

            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            h_zdc[icent]->Draw();
            double maximum =h_minbias[icent]->GetMaximum();  
            h_zdc[icent]->SetMaximum(maximum+ maximum*0.02);
            h_minbias[icent]->Draw("same");

            float X=0.2,Y=0.85;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 33, 43);
            if(icent ==0){
                legend0 = new TLegend(0.2,0.7,0.52,0.8);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.03);
                legend0->Draw();
            }
        }
        c0->SaveAs(Form("%s/PTY_energy0.png",outputpath.c_str()));

        //plot  PTY for effective energy
        input_zdc->Clear();
        input_minbias->Clear();
        c0->Clear();
        c0->Divide(3,3);
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/TemplateFits.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/TemplateFits.root");
        for(int icent=10; icent<20; icent++){
            c0->cd(icent-9);  
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_central_cent%.2i_pericent00_peritrk13_trk14_pta5_ptb05_ch2_deta01",icent)); 
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_central_cent%.2i_pericent00_peritrk13_trk14_pta5_ptb05_ch2_deta01",icent));

            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            h_zdc[icent]->Draw();
            double maximum =h_minbias[icent]->GetMaximum();  
            h_zdc[icent]->SetMaximum(maximum+ maximum*0.02);
            h_minbias[icent]->Draw("same");

            float X=0.2,Y=0.85;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 33, 43);
            if(icent ==10){
                legend0 = new TLegend(0.2,0.7,0.52,0.8);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.03);
                legend0->Draw();
            }
        }
        c0->SaveAs(Form("%s/PTY_energy1.png",outputpath.c_str()));
        #endif

    //plot zdc correlation
    input_zdc->Clear();
    input_minbias->Clear();
    c0->Clear();
    c0->Divide(1,2);  
    input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/histograms.root");
    input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/histograms.root");
    TH2D *h2_minbias = (TH2D*)input_minbias->Get("hZdcCorr");
    TH2D *h2_zdc = (TH2D*)input_zdc->Get("hZdcCorr");
    c0->cd(1);
    h2_zdc->Draw("COLZ");
    Common::myText2(0.7, 0.7, 1,"zdc triggers" , 43, 43);
    gPad->SetLogz();
    c0->cd(2);
    h2_minbias->Draw("COLZ");
    Common::myText2(0.7, 0.7, 1,"minbias triggers" , 43, 43);
    gPad->SetLogz();
    gStyle->SetOptStat(0);
    c0->SaveAs(Form("%s/zdcCorr.png",outputpath.c_str()));

    //plot zdc sum
    SetAtlasStyle();
    c0->Clear();
    c0->Divide(1,1);
    h_zdc[0] = (TH1D*)input_zdc->Get("hzdc"); 
    h_minbias[0] = (TH1D*)input_minbias->Get("hzdc");
    c0->cd(1);
    h_zdc[0]->SetLineColor(kBlue);
    h_minbias[0]->SetLineColor(kRed);
    //h_minbias[0]->GetXaxis()->SetRangeUser(3000,25000);
    h_minbias[0]->Draw("HIST");
    h_zdc[0]->Draw("HIST;same");
    legend0 = new TLegend(0.5,0.83,0.82,0.93);
    legend0->AddEntry(h_zdc[0],"ZDC triggers","l");
    legend0->AddEntry(h_minbias[0],"minbias triggers","l");
    gStyle->SetLegendTextSize(0.03);
    legend0->Draw();
    gPad->SetLogy();
    c0->SaveAs(Form("%s/zdcSum.png",outputpath.c_str()));

    //plot <Nch> per effective energy bin
    c0->Clear();
    c0->Divide(1,1);
    h2_minbias = (TH2D*)input_minbias->Get("hNtrkEff");
    h2_zdc = (TH2D*)input_zdc->Get("hNtrkEff");
    c0->cd(1);
    h_minbias[0] = h2_minbias->ProfileX("h0");
    h_zdc[0] = h2_zdc->ProfileX("h1");
    h_minbias[0]->GetYaxis()->SetTitle("<N_{ch}^{rec}>");
    h_minbias[0]->SetLineColor(kRed); 
    h_zdc[0]->SetLineColor(kBlue);
    h_minbias[0]->Draw();
    h_zdc[0]->Draw("same");
    legend0 = new TLegend(0.5,0.83,0.82,0.93);
    legend0->AddEntry(h_zdc[0],"ZDC triggers","l");
    legend0->AddEntry(h_minbias[0],"minbias triggers","l");
    gStyle->SetLegendTextSize(0.03);
    legend0->Draw();
    c0->SaveAs(Form("%s/avg_multiplicity.png",outputpath.c_str()));

    //plot zdcSum before & after cut
    c0->Clear();
    c0->Divide(1,2);
    h_zdc[0] = (TH1D*)input_zdc->Get("hzdc"); 
    h_minbias[0] = (TH1D*)input_minbias->Get("hzdc");
    h_zdc[1] = (TH1D*)input_zdc->Get("hzdc_cut"); 
    h_minbias[1] = (TH1D*)input_minbias->Get("hzdc_cut");
    h_minbias[0]->SetLineColor(kBlack); 
    h_minbias[1]->SetLineColor(kGreen); 
    h_zdc[0]->SetLineColor(kBlack);
    h_zdc[1]->SetLineColor(kGreen);
    c0->cd(1);
    h_zdc[0]->Draw("HIST");
    h_zdc[1]->Draw("HIST;same");
    float entry1 = h_zdc[0]->GetEntries();
    float entry2 = h_zdc[1]->GetEntries();
    Common::myText2(0.5, 0.7, 1,Form("Total events thrown: %.2f %%",(entry1-entry2)*100/(entry1)) , 43, 43);
    gPad->SetLogy();
    c0->cd(2);
    h_minbias[0]->Draw("HIST");
    h_minbias[1]->Draw("HIST;same");
    entry1 = h_minbias[0]->GetEntries();
    entry2 = h_minbias[1]->GetEntries();
    Common::myText2(0.5, 0.7, 1,Form("Total events thrown: %.2f %%",(entry1-entry2)*100/(entry1)) , 43, 43);
    gPad->SetLogy();
    c0->SaveAs(Form("%s/zdc_sum_cut.png",outputpath.c_str()));
 }
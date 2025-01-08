#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "bins.h"
#include "common.C"

SetAtlasStyle();
TH1D *h_zdc[Bins::NCENT];
TH1D *h_minbias[Bins::NCENT];
TH1D *h_eff_no_ps_zdc;
TH1D *h_eff_no_ps_minbias;
TH1D *h_ratio;
TH1D *h_mean_zdc = new TH1D("h_mean_zdc", ";E_{Eff};", Bins::NCENT - 5, 0.0, Bins::CENT_LO[Bins::NCENT - 5]); h_mean_zdc->SetMarkerSize(1.8);
TH1D *h_mean_minbias = new TH1D("h_mean_minbias", ";E_{Eff};", Bins::NCENT - 5, 0.0, Bins::CENT_LO[Bins::NCENT - 5]); h_mean_minbias->SetMarkerSize(1.8);
THStack *hs_average; // for ploting the <..>
std::map<std::string,double> m_format=Common::StandardFormat();
float X,Y;
#define eta
#define pt
#define ntrk
#define PTY_int
#define Na
#define v2
//#define PTY

/*-----------------------------------------------------------------------------
 *  Plot some basic phyisics to compare between minbias & zdc triggers
 *-----------------------------------------------------------------------------*/

 void plots_monitoring(){
    double findMinimum(double a, double b);
    double findMaximum(double a, double b);
    double minY;
    double maxY;
    char name[100];
    std::string outputpath;
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring/ntrk";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring/eta";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring/pt";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring/pty_int";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring/Na";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring/v2";
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));
    outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/monitoring";

    TFile *input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/histograms.root");
    TFile *input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/histograms.root");
    h_eff_no_ps_zdc = (TH1D*)input_zdc->Get("h_eff_no_ps");
    h_eff_no_ps_minbias = (TH1D*)input_minbias->Get("h_eff_no_ps");
    int nevt_zdc;
    int nevt_minbias;
    

    TCanvas* c0 = new TCanvas("c0","",3000,3000);
    TCanvas* c1 = new TCanvas("c1","",3000,3000);
    TCanvas* c2 = new TCanvas("c2","",2500,3000);
    TLegend *legend0;
    c0->Divide(3,3);
    c1->Divide(3,3);
    c2->Divide(1,2);
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef eta
        //plot eta per event
        for(int icent=0; icent<10; icent++){
            c0->cd(icent+1);
            auto hs = new THStack("hs",";#eta;# per event");
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_eta_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_eta_icent%.2i",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            maxY = findMaximum(h_minbias[icent]->GetMaximum(),h_zdc[icent]->GetMaximum());
            minY = findMaximum(h_minbias[icent]->GetMinimum(),h_zdc[icent]->GetMinimum());
            hs->SetMaximum(1.2*maxY);
            // hs->SetMinimum(1.3*minY);
            hs->Draw("nostack");
             X=0.58,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                legend0 = new TLegend(0.2,0.2,0.68,0.3);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.05);
                legend0->Draw();
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

            c1->cd(icent+1);
            h_ratio = (TH1D*)h_zdc[icent]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->Divide(h_minbias[icent]);
            h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            h_ratio->GetYaxis()->SetTitle("(# per event)_{zdc}/ (# per event)_{minbias}");
            h_ratio->Draw();
            X=0.55,Y=0.9;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
                        
        }
        c0->SaveAs(Form("%s/eta/eta0.png",outputpath.c_str()));
        c1->SaveAs(Form("%s/eta/eta0_ratio.png",outputpath.c_str()));
        c0->Clear();
        c1->Clear();

        input_minbias->Clear();
        input_zdc->Clear();
        //plot eta per event
        c0->Divide(3,3);
        c1->Divide(3,3);
        for(int icent=10; icent<15; icent++){
           c0->cd(icent-9);
            auto hs = new THStack("hs",";#eta;# per event");
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_eta_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_eta_icent%.2i",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            maxY = findMaximum(h_minbias[icent]->GetMaximum(),h_zdc[icent]->GetMaximum());
            minY = findMaximum(h_minbias[icent]->GetMinimum(),h_zdc[icent]->GetMinimum());
            hs->SetMaximum(1.2*maxY);
            hs->SetMinimum(minY);
            hs->Draw("nostack");
             X=0.55,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                legend0 = new TLegend(0.2,0.2,0.68,0.3);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.05);
                legend0->Draw();
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

            c1->cd(icent-9);
            h_ratio = (TH1D*)h_zdc[icent]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->Divide(h_minbias[icent]);
            h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            h_ratio->GetYaxis()->SetTitle("(# per event)_{zdc}/ (# per event)_{minbias}");
            h_ratio->Draw();
            X=0.55,Y=0.9;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
        }
    c0->SaveAs(Form("%s/eta/eta1.png",outputpath.c_str()));
    c1->SaveAs(Form("%s/eta/eta1_ratio.png",outputpath.c_str()));
    c0->Clear();
    c1->Clear();

    //plot <eta> per EE bin
    input_minbias->Clear();
    input_zdc->Clear();
    c2->Clear();
    c2->Divide(1,2);

    for(int icent =0; icent <  Bins::NCENT - 5; icent++){
    h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_eta_icent%.2i",icent));
    h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_eta_icent%.2i",icent));
    nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
    nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
    h_zdc[icent]->Scale(1.0/nevt_zdc);
    h_minbias[icent]->Scale(1.0/nevt_minbias);
    h_mean_zdc->SetBinContent(icent +1, h_zdc[icent]->GetMean());
    h_mean_minbias->SetBinContent(icent +1, h_minbias[icent]->GetMean());
    h_mean_zdc->SetBinError(icent +1, h_zdc[icent]->GetMeanError());
    h_mean_minbias->SetBinError(icent +1, h_minbias[icent]->GetMeanError());
    }

    h_mean_zdc->SetLineColor(kBlue);
    h_mean_minbias->SetLineColor(kRed);
    hs_average = new THStack("hs",";E_{Eff} [TeV];<#eta>");
    hs_average->Add(h_mean_minbias);
    hs_average->Add(h_mean_zdc);
    minY = findMinimum(h_mean_minbias->GetMinimum(), h_mean_zdc->GetMinimum());
    maxY = findMaximum(h_mean_minbias->GetMaximum(), h_mean_zdc->GetMaximum());
    hs_average->SetMaximum(1.3*maxY);
    //hs_average->SetMinimum(1.3*minY);
    c2->cd(1);
    hs_average->Draw("nostack");
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    legend0 = new TLegend(0.7,0.2,0.9,0.3);
    legend0->AddEntry(h_mean_zdc,"ZDC triggers","l");
    legend0->AddEntry(h_mean_minbias,"minbias triggers","l");
    gStyle->SetLegendTextSize(0.035);
    legend0->Draw();
    h_ratio = (TH1D*)h_mean_zdc->Clone();
    h_ratio->SetLineColor(kBlack);
    h_ratio->Divide(h_mean_minbias);
    h_ratio->GetYaxis()->SetRangeUser(0,2.5);
    h_ratio->GetYaxis()->SetTitle("(<#eta>)_{zdc} / (<#eta>)_{minbias}");
    c2->cd(2);
    h_ratio->Draw();
    c2->SaveAs(Form("%s/eta/eta_average.png",outputpath.c_str()));
    #endif
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef pt
        c0->Clear();
        c1->Clear();
        c0->Divide(3,3);
        c1->Divide(3,3);
        input_minbias->Clear();
        input_zdc->Clear();
        //plot pt per event
        for(int icent=0; icent<10; icent++){
            c0->cd(icent+1);
            auto hs = new THStack("hs",";p_{T} [GeV];# per event");
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_pt_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_pt_icent%.2i",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            maxY = findMaximum(h_minbias[icent]->GetMaximum(),h_zdc[icent]->GetMaximum());
            minY = findMaximum(h_minbias[icent]->GetMinimum(),h_zdc[icent]->GetMinimum());
            hs->SetMaximum(1.3*maxY);
            //hs->SetMinimum(1.3*minY);
            hs->Draw("nostack");
            gPad->SetLogy();
             X=0.58,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                legend0 = new TLegend(0.2 +0.2,0.6,0.68 +0.2,0.7);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.05);
                legend0->Draw();
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

            c1->cd(icent+1);
            h_ratio = (TH1D*)h_zdc[icent]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->Divide(h_minbias[icent]);
            h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            h_ratio->GetYaxis()->SetTitle("(# per event)_{zdc}/ (# per event)_{minbias}");
            h_ratio->Draw();
            X=0.55,Y=0.9;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
                        
        }
        c0->SaveAs(Form("%s/pt/pt0.png",outputpath.c_str()));
        c1->SaveAs(Form("%s/pt/pt0_ratio.png",outputpath.c_str()));
        c0->Clear();
        c1->Clear();

        input_minbias->Clear();
        input_zdc->Clear();
        //plot pt per event
        c0->Divide(3,3);
        c1->Divide(3,3);
        for(int icent=10; icent<15; icent++){
           c0->cd(icent-9);
            auto hs = new THStack("hs",";p_{T} [GeV];# per event");
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_pt_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_pt_icent%.2i",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            maxY = findMaximum(h_minbias[icent]->GetMaximum(),h_zdc[icent]->GetMaximum());
            minY = findMaximum(h_minbias[icent]->GetMinimum(),h_zdc[icent]->GetMinimum());
            hs->SetMaximum(1.3*maxY);
            hs->SetMinimum(1.4*minY);
            hs->Draw("nostack");
            gPad->SetLogy();
             X=0.55,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                legend0 = new TLegend(0.2 +0.2,0.6,0.68 +0.2,0.7);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.05);
                legend0->Draw();
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

            c1->cd(icent-9);
            h_ratio = (TH1D*)h_zdc[icent]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->Divide(h_minbias[icent]);
            h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            h_ratio->GetYaxis()->SetTitle("(# per event)_{zdc}/ (# per event)_{minbias}");
            h_ratio->Draw();
            X=0.55,Y=0.9;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
        }
        c0->SaveAs(Form("%s/pt/pt1.png",outputpath.c_str()));
        c1->SaveAs(Form("%s/pt/pt1_ratio.png",outputpath.c_str()));
        c0->Clear();
        c1->Clear();

        //plot <pt> per EE bin
        input_minbias->Clear();
        input_zdc->Clear();
        c2->Clear();
        c2->Divide(1,2);

        for(int icent =0; icent <  Bins::NCENT - 5; icent++){
        h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_pt_icent%.2i",icent));
        h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_pt_icent%.2i",icent));
        nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
        nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
        h_zdc[icent]->Scale(1.0/nevt_zdc);
        h_minbias[icent]->Scale(1.0/nevt_minbias);
        h_mean_zdc->SetBinContent(icent +1, h_zdc[icent]->GetMean());
        h_mean_minbias->SetBinContent(icent +1, h_minbias[icent]->GetMean());
        h_mean_zdc->SetBinError(icent +1, h_zdc[icent]->GetMeanError());
        h_mean_minbias->SetBinError(icent +1, h_minbias[icent]->GetMeanError());
        }

        h_mean_zdc->SetLineColor(kBlue);
        h_mean_minbias->SetLineColor(kRed);
        hs_average = new THStack("hs",";E_{Eff} [TeV];<p_{T}> [GeV]");
        hs_average->Add(h_mean_minbias);
        hs_average->Add(h_mean_zdc);
        minY = findMinimum(h_mean_minbias->GetMinimum(), h_mean_zdc->GetMinimum());
        maxY = findMaximum(h_mean_minbias->GetMaximum(), h_mean_zdc->GetMaximum());
        hs_average->SetMaximum(1.05);
        hs_average->SetMinimum(0.9);
        c2->cd(1);
        hs_average->Draw("nostack");
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        legend0 = new TLegend(0.7,0.2,0.9,0.3);
        legend0->AddEntry(h_mean_zdc,"ZDC triggers","l");
        legend0->AddEntry(h_mean_minbias,"minbias triggers","l");
        gStyle->SetLegendTextSize(0.035);
        legend0->Draw();
        h_ratio = (TH1D*)h_mean_zdc->Clone();
        h_ratio->GetXaxis()->SetTitle("E_{Eff} [TeV]");
        h_ratio->SetLineColor(kBlack);
        h_ratio->Divide(h_mean_minbias);
        //h_ratio->GetYaxis()->SetRangeUser(0,2);
        h_ratio->GetYaxis()->SetTitle("(<p_{T}>)_{zdc} / (<p_{T}>)_{minbias}");
        c2->cd(2);
        //m_format["YTitleOffset"]=1.6;
        //Common::FormatHist(h_ratio,m_format);
        h_ratio->GetYaxis()->SetTitleOffset(1.6);
        h_ratio->Draw();
        c2->SaveAs(Form("%s/pt/pt_average.png",outputpath.c_str()));
    #endif
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef PTY
        //plot  PTY for effective energy
        input_zdc->Clear();
        input_minbias->Clear();
        c0->Clear();
        c0->Divide(3,3);
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/TemplateFits.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/TemplateFits.root");
        for(int icent=0; icent<10; icent++){
            c0->cd(icent+1);
            THStack *hs = new THStack(Form("hs%i",icent), ";#Delta#phi;Y(#Delta#phi)");  
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_central_cent%.2i_pericent00_peritrk00_trk13_pta5_ptb05_ch2_deta01",icent)); 
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_central_cent%.2i_pericent00_peritrk00_trk13_pta5_ptb05_ch2_deta01",icent));

            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            double minimum_minbias =h_minbias[icent]->GetMinimum();  
            double minimum_zdc =h_zdc[icent]->GetMinimum();
            double minimum = findMinimum(minimum_minbias, minimum_zdc);  
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            hs->SetMinimum(0.9*minimum);
            hs->Draw("NOSTACK");

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
            THStack *hs = new THStack(Form("hs%i",icent), ";#Delta#phi;Y(#Delta#phi)");  
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_central_cent%.2i_pericent00_peritrk00_trk14_pta5_ptb05_ch2_deta01",icent)); 
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_central_cent%.2i_pericent00_peritrk00_trk14_pta5_ptb05_ch2_deta01",icent));

            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            double minimum_minbias =h_minbias[icent]->GetMinimum();  
            double minimum_zdc =h_zdc[icent]->GetMinimum();
            double minimum = findMinimum(minimum_minbias, minimum_zdc);  
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            hs->SetMinimum(0.9*minimum);
            hs->Draw("NOSTACK");

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


    
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef ntrk
        c0->Clear();
        c1->Clear();
        c0->Divide(3,3);
        c1->Divide(3,3);
        input_minbias->Clear();
        input_zdc->Clear();
        //plot nch per event
        for(int icent=0; icent<10; icent++){
            c0->cd(icent+1);
            auto hs = new THStack("hs",";N_{ch}^{rec};# per event");
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_ntrk_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_ntrk_icent%.2i",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            maxY = findMaximum(h_minbias[icent]->GetMaximum(),h_zdc[icent]->GetMaximum());
            minY = findMaximum(h_minbias[icent]->GetMinimum(),h_zdc[icent]->GetMinimum());
            // hs->SetMaximum(1.3*maxY);
            // hs->SetMinimum(1.3*minY);
            hs->Draw("nostack");
            gPad->SetLogy();
             X=0.58,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                legend0 = new TLegend(0.2,0.2,0.6,0.3);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.035);
                legend0->Draw();
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

            c1->cd(icent+1);
            h_ratio = (TH1D*)h_zdc[icent]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->Divide(h_minbias[icent]);
            //h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            h_ratio->GetXaxis()->SetRangeUser(0.0,90.0);
            h_ratio->GetYaxis()->SetTitle("(# per event)_{zdc}/ (# per event)_{minbias}");
            h_ratio->Draw();
            X=0.55,Y=0.9;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
                        
        }
        c0->SaveAs(Form("%s/ntrk/ntrk0.png",outputpath.c_str()));
        c1->SaveAs(Form("%s/ntrk/ntrk0_ratio.png",outputpath.c_str()));
        c0->Clear();
        c1->Clear();

        input_minbias->Clear();
        input_zdc->Clear();
        //plot nch per event
        c0->Divide(3,3);
        c1->Divide(3,3);
        for(int icent=10; icent<15; icent++){
           c0->cd(icent-9);
            auto hs = new THStack("hs",";N_{ch}^{rec};# per event");
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_ntrk_icent%.2i",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_ntrk_icent%.2i",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            h_zdc[icent]->SetLineColor(kBlue);
            h_minbias[icent]->SetLineColor(kRed);
            hs->Add(h_zdc[icent]);
            hs->Add(h_minbias[icent]);
            maxY = findMaximum(h_minbias[icent]->GetMaximum(),h_zdc[icent]->GetMaximum());
            minY = findMaximum(h_minbias[icent]->GetMinimum(),h_zdc[icent]->GetMinimum());
            //hs->SetMaximum(1.3*maxY);
            //hs->SetMinimum(1.3*minY);
            hs->Draw("nostack");
            gPad->SetLogy();
             X=0.55,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                legend0 = new TLegend(0.2,0.2,0.6,0.3);
                legend0->AddEntry(h_zdc[icent],"ZDC triggers","l");
                legend0->AddEntry(h_minbias[icent],"minbias triggers","l");
                gStyle->SetLegendTextSize(0.035);
                legend0->Draw();
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

            c1->cd(icent-9);
            h_ratio = (TH1D*)h_zdc[icent]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->Divide(h_minbias[icent]);
            //h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            h_ratio->GetXaxis()->SetRangeUser(0.0,80.0);
            h_ratio->GetYaxis()->SetTitle("(# per event)_{zdc}/ (# per event)_{minbias}");
            h_ratio->Draw();
            X=0.55,Y=0.9;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
        }
        c0->SaveAs(Form("%s/ntrk/ntrk1.png",outputpath.c_str()));
        c1->SaveAs(Form("%s/ntrk/ntrk1_ratio.png",outputpath.c_str()));
        c0->Clear();
        c1->Clear();

        //plot <nch> per EE bin
        input_minbias->Clear();
        input_zdc->Clear();
        c2->Clear();
        c2->Divide(1,2);

        for(int icent =0; icent <  Bins::NCENT - 5; icent++){
        h_zdc[icent] = (TH1D*)input_zdc->Get(Form("h_ntrk_icent%.2i",icent));
        h_minbias[icent] = (TH1D*)input_minbias->Get(Form("h_ntrk_icent%.2i",icent));
        nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
        nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
        h_zdc[icent]->Scale(1.0/nevt_zdc);
        h_minbias[icent]->Scale(1.0/nevt_minbias);
        h_mean_zdc->SetBinContent(icent +1, h_zdc[icent]->GetMean());
        h_mean_minbias->SetBinContent(icent +1, h_minbias[icent]->GetMean());
        h_mean_zdc->SetBinError(icent +1, h_zdc[icent]->GetMeanError());
        h_mean_minbias->SetBinError(icent +1, h_minbias[icent]->GetMeanError());
        }

        h_mean_zdc->SetLineColor(kBlue);
        h_mean_minbias->SetLineColor(kRed);
        hs_average = new THStack("hs",";E_{Eff} [TeV];<N_{ch}^{rec}>");
        hs_average->Add(h_mean_minbias);
        hs_average->Add(h_mean_zdc);
        minY = findMinimum(h_mean_minbias->GetMinimum(), h_mean_zdc->GetMinimum());
        maxY = findMaximum(h_mean_minbias->GetMaximum(), h_mean_zdc->GetMaximum());
        // hs_average->SetMaximum(1.3*maxY);
        hs_average->SetMinimum(5.0);
        c2->cd(1);
        hs_average->Draw("nostack");
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        legend0 = new TLegend(0.7,0.2,0.9,0.3);
        legend0->AddEntry(h_mean_zdc,"ZDC triggers","l");
        legend0->AddEntry(h_mean_minbias,"minbias triggers","l");
        gStyle->SetLegendTextSize(0.035);
        legend0->Draw();
        h_ratio = (TH1D*)h_mean_zdc->Clone();
        h_ratio->SetLineColor(kBlack);
        h_ratio->Divide(h_mean_minbias);
        //h_ratio->GetYaxis()->SetRangeUser(0,2);
        h_ratio->GetXaxis()->SetTitle("E_{Eff} [TeV]");
        h_ratio->GetYaxis()->SetTitle("(<N_{ch}^{rec}>)_{zdc} / (<N_{ch}^{rec}>)_{minbias}");
        c2->cd(2);
        h_ratio->GetYaxis()->SetTitleOffset(1.6);
        h_ratio->Draw();
        c2->SaveAs(Form("%s/ntrk/ntrk_average.png",outputpath.c_str()));
    #endif
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef PTY_int
        //plot PTY integral per effective energy bin
        input_minbias->Clear();
        input_zdc->Clear();
        c2->Clear();
        c2->Divide(1,2);
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/PTY1D.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/PTY1D.root");
        const int n = Bins::NCENT -4;
        TH1D *h_pty_zdc = new TH1D("h_pty_zdc", ";E_{Eff} [TeV];", n,0.0,Bins::CENT_LO[n]);   h_pty_zdc->Sumw2();
        TH1D *h_pty_mb = new TH1D("h_pty_mb", ";E_{Eff} [TeV];", n,0.0,Bins::CENT_LO[n]);     h_pty_mb->Sumw2();
        double integral[2][n] = {0.0};
        double integralError[2][n] = {0.0};


        for(int icent = 0; icent<n; icent ++){
            h_zdc[0] = (TH1D*)input_zdc->Get(Form("PTY_cent%.2i_trk13_pta5_ptb05_ch2_deta01",icent)); 
            h_minbias[0] = (TH1D*)input_minbias->Get(Form("PTY_cent%.2i_trk13_pta5_ptb05_ch2_deta01",icent));
            integral[0][icent] = h_zdc[0]->IntegralAndError(1, h_zdc[0]->GetNbinsX(), integralError[0][icent],"width");
            integral[1][icent] = h_minbias[0]->IntegralAndError(1, h_minbias[0]->GetNbinsX(), integralError[1][icent],"width");
            float norm = (float)h_zdc[0]->GetNbinsX()*h_zdc[0]->GetBinWidth(1);
            h_pty_zdc->SetBinContent(icent + 1, integral[0][icent]/norm);
            h_pty_zdc->SetBinError(icent + 1, integralError[0][icent]/norm);
            h_pty_mb->SetBinContent(icent + 1, integral[1][icent]/norm);
            h_pty_mb->SetBinError(icent + 1, integralError[1][icent]/norm);
            
        }
        hs_average = new THStack("hs",";E_{Eff} [TeV]; # of pairs per N_{a}");
        h_pty_zdc->SetLineColor(kBlue);
        h_pty_mb->SetLineColor(kRed);
        hs_average->Add(h_pty_zdc);
        hs_average->Add(h_pty_mb);
        c2->cd(1);
        hs_average->Draw("nostack");
        hs_average->SetMinimum(0.3);
        hs_average->SetMaximum(0.7);
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS"         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        legend0 = new TLegend(0.7,0.2,0.9,0.3);
        legend0->AddEntry(h_pty_zdc,"ZDC triggers","l");
        legend0->AddEntry(h_pty_mb,"minbias triggers","l");
        gStyle->SetLegendTextSize(0.035);
        legend0->Draw();
        c2->cd(2);
        h_ratio = (TH1D*)h_pty_zdc->Clone();
        h_ratio->SetLineColor(kBlack);
        h_ratio->Divide(h_pty_mb);
        h_ratio->GetYaxis()->SetTitle("#frac{(# of pairs per N_{a})_{zdc}}{(# of pairs per N_{a})_{minbias}}");
        h_ratio->Draw(); 
        h_ratio->GetYaxis()->SetRangeUser(0.87,1.05);
        c2->SaveAs(Form("%s/pty_int/pty_integral.png",outputpath.c_str()));
    #endif
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef Na
        //Na vs Effective energy
        input_minbias->Clear();
        input_zdc->Clear();
        c2->Clear();
        c2->Divide(1,2);
        c0->Clear();
        c0->Divide(1,1);
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/RebinTrk.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/RebinTrk.root");
        h_mean_zdc->Reset();
        h_mean_minbias->Reset();
        for(int icent =0; icent< Bins::NCENT -4; icent++){
            h_zdc[icent] = (TH1D*)input_zdc->Get(Form("N_trigger_cent%.2i_trk13",icent));
            h_minbias[icent] = (TH1D*)input_minbias->Get(Form("N_trigger_cent%.2i_trk13",icent));
            nevt_zdc = h_eff_no_ps_zdc->GetBinContent(icent +1);
            nevt_minbias = h_eff_no_ps_minbias->GetBinContent(icent +1);
            h_zdc[icent]->Scale(1.0/nevt_zdc);
            h_minbias[icent]->Scale(1.0/nevt_minbias);
            int bin_lo  =h_zdc[icent]->FindBin(Bins::PT1_LO[5]+.0001);
            int bin_high=h_zdc[icent]->FindBin(Bins::PT1_HI[5]-.0001);
            double Na_zdc_err;
            double Na_minbias_err;
            double Na_zdc= h_zdc[icent]->IntegralAndError(bin_lo,bin_high,Na_zdc_err);
            double Na_minbias= h_minbias[icent]->IntegralAndError(bin_lo,bin_high,Na_minbias_err);
            h_mean_zdc->SetBinContent(icent +1 , Na_zdc);
            h_mean_minbias->SetBinContent(icent +1 , Na_minbias);
            h_mean_zdc->SetBinError(icent +1 , Na_zdc_err);
            h_mean_minbias->SetBinError(icent +1 , Na_minbias_err);
        }
    //h_mean_zdc->Divide(h_mean_minbias);
        hs_average = new THStack("hs",";E_{Eff} [TeV]; N_{a} per event");
        h_mean_zdc->SetLineColor(kBlue);
        h_mean_minbias->SetLineColor(kRed);
        hs_average->Add(h_mean_zdc);
        hs_average->Add(h_mean_minbias);
        c2->cd(1);
        hs_average->Draw("nostack");
        hs_average->SetMinimum(4.0);
        // hs_average->SetMaximum(0.7);
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS"         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        legend0 = new TLegend(0.7,0.2,0.9,0.3);
        legend0->AddEntry(h_mean_zdc,"ZDC triggers","l");
        legend0->AddEntry(h_mean_minbias,"minbias triggers","l");
        gStyle->SetLegendTextSize(0.035);
        legend0->Draw();
        c2->cd(2);
        h_ratio = (TH1D*)h_mean_zdc->Clone();
        h_ratio->SetLineColor(kBlack);
        h_ratio->Divide(h_mean_minbias);
        h_ratio->GetYaxis()->SetTitle("#frac{(N_{a} per event)_{zdc}}{(N_{a} per event)_{minbias}}");
        h_ratio->Draw(); 
        c2->SaveAs(Form("%s/Na/Na.png",outputpath.c_str()));
    #endif
//--------------------------------------------------------------------------------------------------------------------------------------- 
    #ifdef v2
        //ratio for v2
        input_minbias->Clear();
        input_zdc->Clear();
        c2->Clear();
        c2->Divide(1,2);
        h_mean_zdc->Reset();
        h_mean_minbias->Reset();
        input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/TemplateFits_vnn.root");
        input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/TemplateFits_vnn.root");
        for(int icent =0; icent< Bins::NCENT -4; icent++){
            std::pair<float, float> vnn_zdc=    Bins::GetVnPtb(icent,13,5,5,2,1,2,0,0,input_zdc,input_zdc); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
            std::pair<float, float> vnn_minbias=Bins::GetVnPtb(icent,13,5,5,2,1,2,0,0,input_minbias,input_minbias); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
            //if(vnn_minbias.first <0.0) continue;
            h_mean_zdc->SetBinContent(icent+1, vnn_zdc.first);
            h_mean_zdc->SetBinError(icent+1, vnn_zdc.second);
            h_mean_minbias->SetBinContent(icent+1, vnn_minbias.first);
            h_mean_minbias->SetBinError(icent+1, vnn_minbias.second);
        }
        hs_average = new THStack("hs",";E_{Eff} [TeV]; v_{2}");
        h_mean_zdc->SetLineColor(kBlue);
        h_mean_minbias->SetLineColor(kRed);
        hs_average->Add(h_mean_zdc);
        hs_average->Add(h_mean_minbias);
        c2->cd(1);
        hs_average->Draw("nostack");
        hs_average->SetMaximum(0.2);
        hs_average->SetMinimum(0);
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS"         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        legend0 = new TLegend(0.7,0.8,0.9,0.9);
        legend0->AddEntry(h_mean_zdc,"ZDC triggers","l");
        legend0->AddEntry(h_mean_minbias,"minbias triggers","l");
        gStyle->SetLegendTextSize(0.035);
        legend0->Draw();
        c2->cd(2);
        h_ratio = (TH1D*)h_mean_zdc->Clone();
        h_ratio->SetLineColor(kBlack);
        h_ratio->Divide(h_mean_minbias);
        h_ratio->GetYaxis()->SetRangeUser(0,3);
        h_ratio->GetXaxis()->SetTitle("E_{Eff} [TeV]");
        h_ratio->GetYaxis()->SetTitle("#frac{(v_{2})_{zdc}}{(v_{2})_{minbias}}");
        h_ratio->Draw(); 
        c2->SaveAs(Form("%s/v2/v2.png",outputpath.c_str()));
    #endif
//---------------------------------------------------------------------------------------------------------------------------------------   
// information on zdc

    // //plot zdc sum
    // input_minbias->Clear();
    // input_zdc->Clear();
    // c0->Clear();
    // c0->Divide(1,1);
    // h_zdc[0] = (TH1D*)input_zdc->Get("hzdc"); 
    // h_minbias[0] = (TH1D*)input_minbias->Get("hzdc");
    // c0->cd(1);
    // h_zdc[0]->SetLineColor(kBlue);
    // h_minbias[0]->SetLineColor(kRed);
    // //h_minbias[0]->GetXaxis()->SetRangeUser(3000,25000);
    // h_minbias[0]->Draw("HIST");
    // h_zdc[0]->Draw("HIST;same");
    // legend0 = new TLegend(0.5,0.83,0.82,0.93);
    // legend0->AddEntry(h_zdc[0],"ZDC triggers","l");
    // legend0->AddEntry(h_minbias[0],"minbias triggers","l");
    // gStyle->SetLegendTextSize(0.03);
    // legend0->Draw();
    // gPad->SetLogy();
    // c0->SaveAs(Form("%s/zdcSum.png",outputpath.c_str()));

    // //plot zdcSum before & after cut
    // input_minbias->Clear();
    // input_zdc->Clear();
    // c0->Clear();
    // c0->Divide(1,2);
    // h_zdc[0] = (TH1D*)input_zdc->Get("hzdc"); 
    // h_minbias[0] = (TH1D*)input_minbias->Get("hzdc");
    // h_zdc[1] = (TH1D*)input_zdc->Get("hzdc_cut"); 
    // h_minbias[1] = (TH1D*)input_minbias->Get("hzdc_cut");
    // h_minbias[0]->SetLineColor(kBlack); 
    // h_minbias[1]->SetLineColor(kGreen); 
    // h_zdc[0]->SetLineColor(kBlack);
    // h_zdc[1]->SetLineColor(kGreen);
    // c0->cd(1);
    // h_zdc[0]->Draw("HIST");
    // h_zdc[1]->Draw("HIST;same");
    // float entry1 = h_zdc[0]->GetEntries();
    // float entry2 = h_zdc[1]->GetEntries();
    // Common::myText2(0.5, 0.7, 1,Form("Total events thrown: %.2f %%",(entry1-entry2)*100/(entry1)) , 43, 43);
    // gPad->SetLogy();
    // c0->cd(2);
    // h_minbias[0]->Draw("HIST");
    // h_minbias[1]->Draw("HIST;same");
    // entry1 = h_minbias[0]->GetEntries();
    // entry2 = h_minbias[1]->GetEntries();
    // Common::myText2(0.5, 0.7, 1,Form("Total events thrown: %.2f %%",(entry1-entry2)*100/(entry1)) , 43, 43);
    // gPad->SetLogy();
    // c0->SaveAs(Form("%s/zdc_sum_cut.png",outputpath.c_str()));

    //plot zdc correlation
    // input_zdc->Clear();
    // input_minbias->Clear();
    // c0->Clear();
    // c0->Divide(1,2);  
    // input_zdc = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/histograms.root");
    // input_minbias = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias/histograms.root");
    // TH2D *h2_minbias = (TH2D*)input_minbias->Get("hZdcCorr");
    // TH2D *h2_zdc = (TH2D*)input_zdc->Get("hZdcCorr");
    // c0->cd(1);
    // h2_zdc->Draw("COLZ");
    // Common::myText2(0.7, 0.7, 1,"zdc triggers" , 43, 43);
    // gPad->SetLogz();
    // c0->cd(2);
    // h2_minbias->Draw("COLZ");
    // Common::myText2(0.7, 0.7, 1,"minbias triggers" , 43, 43);
    // gPad->SetLogz();
    // gStyle->SetOptStat(0);
    // c0->SaveAs(Form("%s/zdcCorr.png",outputpath.c_str()));

    // c0->Clear();
    // c0->Divide(3,3);
 }


 double findMinimum(double a, double b) {
    return (a < b) ? a : b;
}

 double findMaximum(double a, double b) {
    return (a > b) ? a : b;
}
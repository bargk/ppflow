
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "bins.h"
#include "common.C"
TH1D *h_zdc[Bins::NCENT];
TH1D *h1[Bins::NCENT];
THStack *hs_average;
TH1D *h1_profile;
TH2D *h2;
TH1D *h_mean_zdc = new TH1D("h_mean_zdc", ";E_{Eff} [TeV];", Bins::NCENT - 4, 0.0, Bins::CENT_LO[Bins::NCENT - 4]); 
TH1D *h_v2_zdc = new TH1D("h_v2_zdc", ";E_{Eff} [TeV];", 8, 0.0, Bins::CENT_HI[30]); 
TH1D *h_v2_xor = new TH1D("h_v2_xor", ";E_{Eff} [TeV];", 4, Bins::CENT_LO[27], Bins::CENT_HI[30]); 
TH1D *h_v2_minbias_1 = new TH1D("h_v2_minbias_1", ";E_{Eff} [TeV];", 1, Bins::CENT_LO[22], Bins::CENT_HI[22]); 
TH1D *h_v2_minbias_2 = new TH1D("h_v2_minbias_2", ";E_{Eff} [TeV];", 1, Bins::CENT_LO[19], Bins::CENT_HI[19]); 
h_v2_minbias_1->SetMarkerStyle(20);
h_v2_minbias_2->SetMarkerStyle(20);
h_v2_zdc->SetMarkerStyle(20);
h_v2_xor->SetMarkerStyle(20);

h_v2_minbias_1->SetMarkerSize(2.1);
h_v2_minbias_2->SetMarkerSize(2.1);
h_v2_zdc->SetMarkerSize(2.1);
h_v2_xor->SetMarkerSize(2.1);
TFile *input;
TFile *input1;
TFile *input2;
TFile *input3;
h_mean_zdc->SetMarkerSize(2);
h_mean_zdc->SetMarkerStyle(20);
float X,Y;

//#define pt


/*-----------------------------------------------------------------------------
 *  Plot some basic phyisics for analysis (not including Template)
 *-----------------------------------------------------------------------------*/
 std::string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/1.5sigma/sameSide";
void plots_analysis(int Trig1 =0){
    int Trig = Trig1; 
    SetAtlasStyle();
    std::string inputpath;
    std::string outputpath;
    if( Trig == 0){
        std::cout << "Working on AND trigger " << std::endl;
        outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/plots_analysis";
        inputpath = Form("%s",base.c_str());
    }
    else if( Trig == 1){
        std::cout << "Working on minbias trigger " << std::endl;
        outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/plots_analysis/minbias";
        inputpath = Form("%s/minbias",base.c_str());
    }
    else if( Trig == 2){
        std::cout << "Working on XOR trigger " << std::endl;
        outputpath = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/plots_analysis/xorE2";
        inputpath = Form("%s/xorE2",base.c_str());
    }
    else{
        std::cerr << "Error: no valid trigger index provided!" <<std::endl;
        exit(-1);
    }
    gSystem->Exec(Form("mkdir -p %s",outputpath.c_str()));

   
    input = new TFile(Form("%s/histograms.root",inputpath.c_str()));

    TLegend *legend0;
    TCanvas* c0 = new TCanvas("c0","",3000,3000);
    TCanvas* c1 = new TCanvas("c1","",2500,2500);
    TCanvas* c2 = new TCanvas("c2","",3000,2500);
    c0->Divide(3,3);
    //c1->Divide(3,3);
    //c2->Divide(1,2);

    #ifdef pt
    //plot pt
        c0->Clear();
        c1->Clear();
        c0->Divide(3,3);
        c1->Divide(3,3);
        input->Clear();
        //plot pt per event
        for(int icent=0; icent<10; icent++){
            c0->cd(icent+1);
            auto hs = new THStack("hs",";p_{T} [GeV];Counts");
            h_zdc[icent] = (TH1D*)input->Get(Form("h_pt_icent%.2i",icent));
            h_zdc[icent]->SetLineColor(kBlue);
            hs->Add(h_zdc[icent]);
            hs->Draw("nostack");
            gPad->SetLogy();
             X=0.58,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 0){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }

        }
        c0->SaveAs(Form("%s/pt0.png",outputpath.c_str()));
        c0->Clear();

        c0->Divide(3,3);
        for(int icent =10; icent<15; icent++){
            c0->cd(icent-9);
            auto hs = new THStack("hs",";p_{T} [GeV];Counts");
            h_zdc[icent] = (TH1D*)input->Get(Form("h_pt_icent%.2i",icent));
            h_zdc[icent]->SetLineColor(kBlue);
            hs->Add(h_zdc[icent]);
            hs->Draw("nostack");
            gPad->SetLogy();
             X=0.58,Y=0.88;
            Common::myText2(X, Y, 1, Bins::label_cent (icent) , 40, 43);
            if(icent == 10){
                X=0.20,Y=0.88;
                Common::myText2(X     ,Y,1,"ATLAS "         ,40,73);
                Common::myText2(X+0.16,Y,1,Common::Internal,40,43);Y=Y-0.05;
                Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 40, 43);Y=Y-0.05;
            }
        }
        c0->SaveAs(Form("%s/pt1.png",outputpath.c_str()));

        // plot <pt>
        c1->Clear();

        for(int icent =0; icent <  Bins::NCENT - 5; icent++){
        h_zdc[icent] = (TH1D*)input->Get(Form("h_pt_icent%.2i",icent));
        h_mean_zdc->SetBinContent(icent +1, h_zdc[icent]->GetMean());
        h_mean_zdc->SetBinError(icent +1, h_zdc[icent]->GetMeanError());
        }

        h_mean_zdc->SetLineColor(kBlue);
        h_mean_zdc->GetYaxis()->SetTitle("<p_{T}> [GeV]");
        c1->cd(1);
        h_mean_zdc->Draw();
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        c1->SaveAs(Form("%s/pt_avg.png",outputpath.c_str()));
    #endif

    //plot Nch vs Eeff and <nch> vs Eeff
    c1->Clear();
    c1->Divide(1,1);
    h2 = (TH2D*)input->Get("hNtrkEff");
    h2->GetXaxis()->SetRange(1,15);
    h2->GetYaxis()->SetRange(1,13);
    h2->GetXaxis()->SetTitle("E_{Eff} [TeV]");
    if(Trig == 2){
        h2->GetXaxis()->SetRange(10,17);
    }
    gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
    h2->Draw("colz");
    gPad->SetLogz();
    h1_profile = h2->ProfileX("hprofile");
    h1_profile->GetYaxis()->SetTitle("<N_{ch}>");
    h1_profile->SetMarkerStyle(20);
    h1_profile->SetMarkerSize(2);
    h1_profile->Draw("same");
    X=0.20,Y=0.9;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    c1->SaveAs(Form("%s/nch_eff.png",outputpath.c_str()));
    c1->Clear();
    h1_profile->Draw();
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    c1->SaveAs(Form("%s/avg_nch.png",outputpath.c_str()));

    // //plot PTY peripheral bin
    // c1->Clear();
    // input = new TFile(Form("%s/PTY1D.root",inputpath.c_str()));
    // h_zdc[0] = (TH1D*)input->Get("PTY_cent00_trk00_pta5_ptb05_ch2_deta01");
    // h_zdc[0]->GetXaxis()->SetTitle("#Delta#phi");
    // h_zdc[0]->GetYaxis()->SetTitle("Y^{peri}(#Delta#phi)");
    // h_zdc[0]->SetMarkerStyle(20);
    // h_zdc[0]->SetMarkerSize(2);
    // float X=0.2,Y=0.9;
    // h_zdc[0]->Draw();
    // X=0.20,Y=0.88;
    // Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    // Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    // Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    // Common::myText2(X, Y, 1, "0.5 < P^{a,b}_{T} < 5 GeV" , 80, 43);Y=Y-0.05;
    // std::string peri_label=Bins::label_cent_peri(0) + " , " + Bins::label_trk_peri(0);
    // Common::myText2(X, Y, 1, "E_{Eff}<0.68 TeV" , 80, 43);Y=Y-0.05;
    // Common::myText2(X, Y, 1, "N_{ch}^{rec} < 10" , 80, 43);Y=Y-0.05;

    // c1->SaveAs(Form("%s/peri_bin.png",outputpath.c_str()));

    //plot events in each bin of EE for all the 3 triggers
    c1->Clear();
    input = new TFile(Form("%s/histograms.root",base.c_str()));
    input1 = new TFile(Form("%s/minbias/histograms.root",base.c_str()));
    input2 = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
    h_zdc[0] = (TH1D*)input->Get("h_eff_no_ps");     h_zdc[0]->SetLineColor(kBlue);     h_zdc[0]->SetMarkerColor(kBlue);
    h_zdc[1] = (TH1D*)input1->Get("h_eff_no_ps");    h_zdc[1]->SetLineColor(kRed);      h_zdc[1]->SetMarkerColor(kRed);
    h_zdc[2] = (TH1D*)input2->Get("h_eff_no_ps");    h_zdc[2]->SetLineColor(kBlack);   h_zdc[2]->SetMarkerColor(kBlack);
    double maxY = 0; // Variable to store the global maximum

    // First loop to determine the maximum Y value
    for (int i = 0; i < 3; i++) {
        double localMax = h_zdc[i]->GetMaximum();
        if (localMax > maxY) {
            maxY = localMax;
        }
    }
    c2->Divide(1,1);
    c2->cd();
    for(int i =0; i<3; i++){
        h_zdc[i]->GetXaxis()->SetLabelSize(0.04); // Adjust size as needed
        h_zdc[i]->GetXaxis()->LabelsOption("v"); // Rotate labels vertically for readability
        //h_zdc[i]->GetXaxis()->SetRange(1,15); 
        h_zdc[i]->GetXaxis()->SetTitle("E_{Eff} [TeV]"); 
        h_zdc[i]->GetYaxis()->SetTitle("Events"); 
        h_zdc[i]->SetMarkerStyle(20);
        h_zdc[i]->SetMarkerSize(2.1);
        h_zdc[i]->SetMaximum(maxY * 10.0); // Set Y-axis maximum with some margin
        if (i == 0) h_zdc[i]->Draw("E1");
        else{h_zdc[i]->Draw("E1;SAME");}
    }
    h_zdc[0]->GetXaxis()->SetTitleOffset(2.25); // Default is ~0.01; increase for more space
    legend0 = new TLegend(0.4,0.25,0.7,0.4);
    legend0->AddEntry(h_zdc[0],"Two sided ZDC trigger","l");
    legend0->AddEntry(h_zdc[1],"Minimum bias","l");
    legend0->AddEntry(h_zdc[2],"One side only ZDC triggers","l");
    legend0->Draw();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.25);
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
    c2->SaveAs(Form("%s/nevents.png",outputpath.c_str()));
//--------------------------------------------------------------------------------------------------
    c2->Clear();
    c2->Divide(1,1);
    input = new TFile(Form("%s/histograms.root",inputpath.c_str()));
    h_zdc[0] = (TH1D*)input->Get("heff");     h_zdc[0]->SetLineColor(kBlue);     h_zdc[0]->SetMarkerColor(kBlue);
    h_zdc[0]->Draw("E1");
    h_zdc[0]->GetXaxis()->SetTitle("E_{Eff} [TeV]");
    h_zdc[0]->GetYaxis()->SetTitle("Events");
    gPad->SetLogy(0);
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
    c2->SaveAs(Form("%s/effetive_energy.png",outputpath.c_str()));

//     //plot v2 with zdc trigger and minbias trigger
//     c2->Clear();
//     c2->Divide(1,1);
//    input = new TFile(Form("%s/TemplateFits_vnn.root",base.c_str()));
//     input1 = new TFile(Form("%s/minbias/TemplateFits_vnn.root",base.c_str()));
//     input2 = new TFile(Form("%s/xorE2/TemplateFits_vnn.root",base.c_str()));
//     if(Trig ==2){
//         input = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/xorE2/TemplateFits_vnn.root");
//     }
//     std::pair<float, float> vnn_zdc;
//     std::pair<float, float> vnn_minbias1;
//     std::pair<float, float> vnn_minbias2;
//     int pericent= 0; //defualt peripheral bin for zdc trigger
    
//     h_mean_zdc->Reset();
//     for(int icent =0; icent< 8; icent++){
//                 vnn_zdc=      Bins::GetVnPtb(icent +23,13,5,5,2,1,2,pericent,0,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
//                 h_v2_zdc->SetBinContent(icent+1, vnn_zdc.first);
//                 h_v2_zdc->SetBinError(icent+1, vnn_zdc.second);
//         }
//     for(int icent =0; icent< 4; icent++){
//                 vnn_zdc=      Bins::GetVnPtb(icent +27,13,5,5,2,1,2,pericent,0,input2,input2); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
//                 h_v2_xor->SetBinContent(icent+1, vnn_zdc.first);
//                 h_v2_xor->SetBinError(icent+1, vnn_zdc.second);
//         }
//         vnn_minbias1=      Bins::GetVnPtb(22,13,5,5,2,1,2,pericent,0,input1,input1); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
//         h_v2_minbias_1->SetBinContent(1, vnn_minbias1.first);
//         h_v2_minbias_1->SetBinError(1, vnn_minbias1.second);
//         vnn_minbias2=      Bins::GetVnPtb(19,13,5,5,2,1,2,pericent,0,input1,input1); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
//         h_v2_minbias_2->SetBinContent(1, vnn_minbias2.first);
//         h_v2_minbias_2->SetBinError(1, vnn_minbias2.second);
//         hs_average = new THStack("hs",";E_{Eff} [TeV]; v_{2}(p_{T}^{b})");
//         h_v2_zdc->SetLineColor(kBlue);
//         h_v2_zdc->SetMarkerColor(kBlue);
//         h_v2_xor->SetLineColor(kOrange+4);
//         h_v2_xor->SetMarkerColor(kOrange+4);

//         h_v2_minbias_1->SetLineColor(kRed);
//         h_v2_minbias_2->SetLineColor(kRed);
//         h_v2_minbias_1->SetMarkerColor(kRed);
//         h_v2_minbias_2->SetMarkerColor(kRed);

//         h_v2_xor->GetXaxis()->SetRange(10,17);
//         hs_average->Add(h_v2_zdc);
//         hs_average->Add(h_v2_xor);
//         if(Trig ==0){ //take minbias only when in AND trigger
//             hs_average->Add(h_v2_minbias_1);
//             hs_average->Add(h_v2_minbias_2);
//         }
//         hs_average->SetMaximum(0.1);
//         hs_average->SetMinimum(0);
//         c2->cd(1);
//         hs_average->Draw("nostack;E1");
//          if(Trig ==2){
//             hs_average->GetXaxis()->SetLimits(5.44,11.56);
//         }
//         X=0.20,Y=0.88;
//         Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
//         Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
//         Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
//         legend0->Clear();
//         legend0 = new TLegend(0.5,0.8,0.8,0.9);
//         if(Trig ==0){
//             legend0->AddEntry(h_v2_zdc,"ZDC AND ","lep");
//             legend0->AddEntry(h_v2_minbias_1,"Minimum bias","lep");
//             legend0->AddEntry(h_v2_xor,"ZDC XOR","lep");
//             legend0->Draw();
//         }
//         c2->SaveAs(Form("%s/v2_all_range.png",outputpath.c_str()));

        //plot each side of zdc energy
        if(Trig ==2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input =new TFile(Form("%s/histograms.root",base.c_str()));
        }
        c2->Clear();
        c2->Divide(1,1);
        h_zdc[0] = (TH1D*)input->Get("hzdc_C_without_pileup");     h_zdc[0]->SetLineColor(kGreen+1);     h_zdc[0]->SetMarkerColor(kGreen+1);
        h_zdc[1] = (TH1D*)input->Get("hzdc_A_without_pileup");     h_zdc[1]->SetLineColor(kMagenta);     h_zdc[1]->SetMarkerColor(kMagenta);
        h_zdc[0]->SetMarkerSize(2);
        h_zdc[1]->SetMarkerSize(2);
        h_zdc[0]->SetMarkerStyle(20);
        h_zdc[1]->SetMarkerStyle(20);
        h_zdc[0]->GetXaxis()->SetTitle("Energy [GeV]");
        h_zdc[0]->Draw();
        h_zdc[1]->Draw("same");
        gPad->SetLogy();
        X=0.60,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
        legend0->Clear();
        legend0 = new TLegend(0.5,0.3,0.8,0.4);
        legend0->AddEntry(h_zdc[0],"side C","lep");
        legend0->AddEntry(h_zdc[1],"side A","lep");
        legend0->Draw();
        c2->SaveAs(Form("%s/zdc_sides.png",outputpath.c_str()));

        h_zdc[0] = (TH1D*)input->Get("hzdc_A_without_pileup");     h_zdc[0]->SetLineColor(kRed);     h_zdc[0]->SetMarkerColor(kRed);
        h_zdc[1] = (TH1D*)input->Get("hzdc_A_with_pileup");     h_zdc[0]->SetLineColor(kGreen);     h_zdc[0]->SetMarkerColor(kGreen);
        h_zdc[0]->SetMarkerSize(2);
        h_zdc[1]->SetMarkerSize(2);
        h_zdc[0]->SetMarkerStyle(20);
        h_zdc[1]->SetMarkerStyle(20);
        h_zdc[0]->GetXaxis()->SetTitle("Energy [GeV]");
        h_zdc[0]->GetXaxis()->SetRangeUser(1600,10000);
        h_zdc[1]->GetXaxis()->SetRangeUser(1600,10000);
        h_zdc[0]->Draw();
        h_zdc[1]->Draw("same");
        gPad->SetLogy();
        X=0.60,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
        legend0->Clear();
        legend0 = new TLegend(0.5,0.3,0.8,0.4);
        legend0->AddEntry(h_zdc[1],"before pileup cut","lep");
        legend0->AddEntry(h_zdc[0],"after pileup cut","lep");
        legend0->Draw();
        c2->SaveAs(Form("%s/zdc_pileup.png",outputpath.c_str()));

        //plot zdc correlation
        c1->Clear();
        c1->Divide(1,1);
        if(Trig == 2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input = new TFile(Form("%s/histograms.root",base.c_str()));
        }

        h2 = (TH2D*)input->Get("hZdcCorr");
        Common::FormatHist(h2, Common::StandardFormat());
        h2->Draw("colz");
        h2->GetYaxis()->SetTitleOffset(1.5);
        h2->GetXaxis()->SetTitle("side C [GeV]");
        h2->GetYaxis()->SetTitle("side A [GeV]");
        // h2->GetYaxis()->SetRangeUser(0,13000);
        // h2->GetXaxis()->SetRangeUser(0,13000);
        gPad->SetLogz();

        X=0.60,Y=0.92;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
        c1->SaveAs(Form("%s/zdc_corr.png",outputpath.c_str()));

        //plot zdc distribtuion
        if(Trig ==2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input = new TFile(Form("%s/histograms.root",base.c_str()));
        }
        c2->Clear();
        c2->Divide(1,1);
        h_zdc[0] = (TH1D*)input->Get("hzdc");     h_zdc[0]->SetLineColor(kBlue);     h_zdc[0]->SetMarkerColor(kBlue);
        h_zdc[0]->SetMarkerSize(2);
        h_zdc[0]->SetMarkerStyle(20);
        Common::FormatHist(h_zdc[0], Common::StandardFormat());
        h_zdc[0]->Draw();
        gPad->SetLogy();
        h_zdc[0]->GetYaxis()->SetTitleOffset(1.5);
        h_zdc[0]->GetXaxis()->SetRange(1,h_zdc[0]->FindBin(40000));
        X=0.60,Y=0.92;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        double int1 = h_zdc[0]->Integral(h_zdc[0]->FindBin(13600) +1, h_zdc[0]->FindBin(50000));
        double int2 = h_zdc[0]->Integral();
        Common::myText2(X       , Y, 1, Form("Events thrown ~ %.2f %%", int1*100/int2), 70, 43);
        c2->SaveAs(Form("%s/zdc.png",outputpath.c_str()));

        //plot lucrod correlation 
        c2->Clear();
        c2->Divide(1,1);
        if(Trig == 2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input = new TFile(Form("%s/histograms.root",base.c_str()));
        }

        h2 = (TH2D*)input->Get("hLucrodCorr");
        Common::FormatHist(h2, Common::StandardFormat());
        h2->Draw("colz");
        h2->GetYaxis()->SetTitleOffset(1.5);
        h2->GetYaxis()->SetRangeUser(0,500);
        h2->GetXaxis()->SetRangeUser(0,500);
        if(Trig ==0){
            h2->GetYaxis()->SetRangeUser(0,4300);
            h2->GetXaxis()->SetRangeUser(0,4300);
        }
        gPad->SetLogz();
        gPad->SetLogy(0);
        gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
        X=0.60,Y=0.92;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        c2->SaveAs(Form("%s/lucrod_corr.png",outputpath.c_str()));

        //plot lucrod vs zdc side C
        c2->Clear();
        c2->Divide(1,1);
        if(Trig == 2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input = new TFile(Form("%s/histograms.root",base.c_str()));
        }

        h2 = (TH2D*)input->Get("hLucrodZdcCorr0");
        Common::FormatHist(h2, Common::StandardFormat());
        h2->Draw("colz");
        h2->GetYaxis()->SetTitleOffset(1.5);
        h2->GetYaxis()->SetRangeUser(0,6000);
        h2->GetXaxis()->SetRangeUser(0,500);
        if(Trig ==0){
            h2->GetYaxis()->SetRangeUser(0,8000);
            h2->GetXaxis()->SetRangeUser(0,2000);
        }
        gPad->SetLogz();
        //gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
        X=0.60,Y=0.92;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        c2->SaveAs(Form("%s/lucrodZdc_corr0.png",outputpath.c_str()));

        //plot lucrod vs zdc side C
        c2->Clear();
        c2->Divide(1,1);
       if(Trig == 2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input = new TFile(Form("%s/histograms.root",base.c_str()));
        }


        h2 = (TH2D*)input->Get("hLucrodZdcCorr1");
        Common::FormatHist(h2, Common::StandardFormat());
        h2->Draw("colz");
        h2->GetYaxis()->SetTitleOffset(1.5);
        h2->GetYaxis()->SetRangeUser(0,6000);
        h2->GetXaxis()->SetRangeUser(0,500);
        if(Trig ==0){
            h2->GetYaxis()->SetRangeUser(0,8000);
            h2->GetXaxis()->SetRangeUser(0,2000);
        }
        gPad->SetLogz();
        //gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
        X=0.60,Y=0.92;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        c2->SaveAs(Form("%s/lucrodZdc_corr1.png",outputpath.c_str()));

        //plot zdc side correlation - uncalibrated
        c2->Clear();
        c2->Divide(1,1);
        if(Trig == 2){
            input = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
        }
        else{
            input = new TFile(Form("%s/histograms.root",base.c_str()));
        }


        h2 = (TH2D*)input->Get("hAmpCorr");
        Common::FormatHist(h2, Common::StandardFormat());
        h2->Draw("colz");
        h2->GetYaxis()->SetTitleOffset(1.5);
        if(Trig ==2){
            h2->GetYaxis()->SetRangeUser(0,6000);
            h2->GetXaxis()->SetRangeUser(0,6000);
        }
        else{
            h2->GetYaxis()->SetRangeUser(0,8000);
            h2->GetXaxis()->SetRangeUser(0,8000);
        }
        gPad->SetLogz();
        //gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
        X=0.60,Y=0.92;
        Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
        c2->SaveAs(Form("%s/zdc_amp_corr.png",outputpath.c_str()));
}


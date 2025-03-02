#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana";
TH1D *h_zdc[3];
TH1D *h1_profile;
TH2D *h2;
TLegend *legend0;
void plot_ana(){
    SetAtlasStyle();
    TCanvas *c0 = new TCanvas("c0","",3000,2500);
    //plot events in each bin of EE for all the 3 triggers
    TFile *input = new TFile(Form("%s/histograms.root",base.c_str()));
    TFile *input1 = new TFile(Form("%s/minbias/histograms.root",base.c_str()));
    TFile *input2 = new TFile(Form("%s/xorE2/histograms.root",base.c_str()));
    h_zdc[0] = (TH1D*)input->Get("h_eff_no_ps");     h_zdc[0]->SetLineColor(kBlue);     h_zdc[0]->SetMarkerColor(kBlue);
    h_zdc[1] = (TH1D*)input1->Get("h_eff_no_ps");    h_zdc[1]->SetLineColor(kRed);      h_zdc[1]->SetMarkerColor(kRed);
    h_zdc[2] = (TH1D*)input2->Get("h_eff_no_ps");    h_zdc[2]->SetLineColor(kOrange+3);   h_zdc[2]->SetMarkerColor(kOrange+3);
    double maxY = 0; // Variable to store the global maximum

    // First loop to determine the maximum Y value
    for (int i = 0; i < 3; i++) {
        double localMax = h_zdc[i]->GetMaximum();
        if (localMax > maxY) {
            maxY = localMax;
        }
    }
    for(int i =0; i<3; i++){
        h_zdc[i]->GetXaxis()->SetLabelSize(0.04); // Adjust size as needed
        h_zdc[i]->GetXaxis()->LabelsOption("v"); // Rotate labels vertically for readability
        //h_zdc[i]->GetXaxis()->SetRange(1,15); 
        h_zdc[i]->GetXaxis()->SetTitle("E_{Eff} [TeV]"); 
        h_zdc[i]->GetYaxis()->SetTitle("Events"); 
        h_zdc[i]->SetMarkerStyle(20);
        h_zdc[i]->SetMarkerSize(3);
        h_zdc[i]->SetMaximum(maxY * 100.0); // Set Y-axis maximum with some margin
        if (i == 0) h_zdc[i]->Draw("E1");
        else{h_zdc[i]->Draw("E1;SAME");}
    }
    h_zdc[0]->GetXaxis()->SetTitleOffset(2.25); // Default is ~0.01; increase for more space
    legend0 = new TLegend(0.2,0.65,0.5,0.8);
    legend0->AddEntry(h_zdc[0],"Two sided ZDC trigger","lep");
    legend0->AddEntry(h_zdc[1],"Minimum bias triggers","lep");
    legend0->AddEntry(h_zdc[2],"One side only ZDC triggers","lep");
    legend0->Draw();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.25);
    float X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);
    c0->SaveAs(Form("%s/nevents.pdf",figures.c_str()));

    base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide";
    figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana";
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    //plot 2D correlation between nrec to effective energy
    input = new TFile(Form("%s/histograms.root",base.c_str()));
    c0->Clear();
    c0->Divide(1,1);
    h2 = (TH2D*)input->Get("hNtrkEff");
    h2->GetXaxis()->SetRange(1,15);
    h2->GetYaxis()->SetRange(1,13);
    h2->GetXaxis()->SetTitle("E_{Eff} [TeV]");
    // if(Trig == 2){
    //     h2->GetXaxis()->SetRange(10,17);
    // }
    gPad->SetRightMargin(0.15);  // Default is 0.1, increase for more space
    h2->Draw("colz");
    gPad->SetLogz();
    gPad->SetLogy(0);
    h1_profile = h2->ProfileX("hprofile");
    h1_profile->GetYaxis()->SetTitle("<N^{rec}_{ch}>");
    h1_profile->SetMarkerStyle(20);
    h1_profile->SetMarkerSize(4);
    // h1_profile->Draw("same");
    X=0.20,Y=0.9;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    c0->SaveAs(Form("%s/nch_eff.pdf",figures.c_str()));
    c0->Clear();
    h1_profile->Draw();
    X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 70, 43);Y=Y-0.05;
    c0->SaveAs(Form("%s/avg_nch.pdf",figures.c_str()));

}
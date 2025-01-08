#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/1.5sigma";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana";
TH1D *h_zdc[3];
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
    h_zdc[2] = (TH1D*)input2->Get("h_eff_no_ps");    h_zdc[2]->SetLineColor(kOrange);   h_zdc[2]->SetMarkerColor(kOrange);
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
        h_zdc[i]->SetMaximum(maxY * 10.0); // Set Y-axis maximum with some margin
        if (i == 0) h_zdc[i]->Draw("E1");
        else{h_zdc[i]->Draw("E1;SAME");}
    }
    h_zdc[0]->GetXaxis()->SetTitleOffset(2.25); // Default is ~0.01; increase for more space
    legend0 = new TLegend(0.45,0.25,0.75,0.4);
    legend0->AddEntry(h_zdc[0],"Two sided ZDC trigger","lep");
    legend0->AddEntry(h_zdc[1],"Minimum bias triggers","lep");
    legend0->AddEntry(h_zdc[2],"One side only ZDC triggers","lep");
    legend0->Draw();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.25);
    float X=0.20,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp} 22 #sqrt{#it{s}} = 13.6 TeV", 70, 43);
    c0->SaveAs(Form("%s/nevents.pdf",figures.c_str()));
}
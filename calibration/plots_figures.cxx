
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/1.5sigma";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/figures";

TH1D * h_sum[2];
TFile *file_calibrated;
TFile *file_uncalibrated;

float X,Y;
void plots_figures(){
    file_calibrated = new TFile(Form("%s/zdc_calibrated.root",base.c_str()),"read");
    file_uncalibrated = new TFile(Form("%s/zdc_uncalib.root",base.c_str()),"read");
    TCanvas* c0 = new TCanvas("c0","",3000,2500);
    TLegend *legend0;
    SetAtlasStyle();    


    //plot uncalibrated sum
std::vector<std::string> histNames = {"h_c_opposite", "h_a_opposite"};
std::vector<int> lineColors = {kGreen + 1, kMagenta};
std::vector<std::string> fileNames = {"side_c_uncalib.png", "side_a_uncalib.png"};

// Loop over histograms
for (size_t i = 0; i < histNames.size(); i++) {
    h_sum[i] = (TH1D*)file_uncalibrated->Get(histNames[i].c_str());

    h_sum[i]->SetLineColor(lineColors[i]);
    h_sum[i]->SetMarkerColor(lineColors[i]);
    Common::FormatHist(h_sum[i], Common::StandardFormat());
    h_sum[i]->SetMarkerStyle(20);
    h_sum[i]->SetMarkerSize(2);

    h_sum[i]->Draw();
    h_sum[i]->GetFunction("fit")->SetLineWidth(10);
    h_sum[i]->GetFunction("fit")->Draw("same");

    double X = 0.5, Y = 0.85;
    Common::myText2(X, Y, 1, "ATLAS ", 70, 73);
    Common::myText2(X + 0.08, Y, 1, Common::Internal, 70, 43); Y -= 0.05;
    Common::myText2(X, Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", 70, 43); Y -= 0.05;
    c0->SaveAs(Form("%s/%s", figures.c_str(), fileNames[i].c_str()));
    c0->Clear();
    }

    //plot after calibration
    c0->Clear();
    c0->Divide(1,1);
    c0->cd();
    h_sum[0] = (TH1D*)file_calibrated->Get("h_c");
    h_sum[1] = (TH1D*)file_calibrated->Get("h_a");

    Common::FormatHist(h_sum[0], Common::StandardFormat());
    Common::FormatHist(h_sum[1], Common::StandardFormat());
    h_sum[0]->SetLineColor(kGreen + 1);
    h_sum[1]->SetLineColor(kMagenta);
    h_sum[0]->SetMarkerColor(kGreen + 1);
    h_sum[1]->SetMarkerColor(kMagenta);
    h_sum[1]->SetMarkerSize(2);
    h_sum[0]->SetMarkerSize(2);
    h_sum[1]->SetMarkerStyle(20);
    h_sum[0]->SetMarkerStyle(20);
    h_sum[1]->GetXaxis()->SetRange(1,100);
    h_sum[1]->GetXaxis()->SetTitle("Energy [GeV]");
    h_sum[1]->Draw("");
    h_sum[0]->Draw("same");
    legend0 = new TLegend(0.78,0.7,0.88,0.8);
    legend0->AddEntry(h_sum[0],"side C","lep");
    legend0->AddEntry(h_sum[1],"side A","lep");
    legend0->Draw();
    h_sum[1]->GetYaxis()->SetTitleOffset(1.0);
    X=0.45,Y=0.85;
    Common::myText2(X     ,Y,1,"ATLAS "         ,70,73);
    Common::myText2(X+0.08,Y,1,Common::Internal,70,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{PbPb 23}, #sqrt{#it{s_{NN}}} = 5.36 TeV", 70, 43);Y=Y-0.05;

    c0->SaveAs(Form("%s/zdc_calib.png",figures.c_str()));
    c0->Clear();
    //c0->SaveAs(Form("%s/side_a_calib.png",figures.c_str()));
    
}
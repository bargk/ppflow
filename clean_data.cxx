//generating list of LB when the module is overflow

TH1D *h1[2][4];
TH2D *h2[2][4];
void clean_data(){
    TFile *fout = new TFile("badLB.root", "UPDATE");
    TFile *finput1 = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/Amp/histograms.root","read");
    TH1D *h_lb[2][4];
    for (int side=0; side<2; side++){
        for(int mod =0; mod<4; mod++){
            h_lb[side][mod] = new TH1D(Form("h_lb_%i%i",side,mod),";LB;",3800,1000,4800);
        }
    }
    int sideAna =0;
    int modAna = 2;
    h2[sideAna][modAna] = (TH2D*)finput1->Get("h0_2")->Clone("h01");
    h1[sideAna][modAna] = (TH1D*)h2[sideAna][modAna]->ProfileX("h01profile");
    float threshold = 435.0;
    int lowlim = h1[sideAna][modAna]->GetBinLowEdge(h1[sideAna][modAna]->FindBin(1000));
    int highlim = h1[sideAna][modAna]->GetBinLowEdge(h1[sideAna][modAna]->FindBin(1900));
    TF1 *constFit = new TF1("constFit", "[0]", lowlim, highlim);
    constFit->SetParameter(0, 410); // Initial guess for the constant value
//     for(int bin =lowbin; bin< highbin +1; bin++){
//         float amp = h1[sideAna][modAna]->GetBinContent(bin);
//         if(amp > threshold){
//             cout << bin << endl;
//             h_lb[sideAna][modAna]->SetBinContent(bin,1);
//             h_lb[sideAna][modAna]->SetBinError(bin,0);
//         } 
//     }
// fout->cd();
// h_lb[sideAna][modAna]->Write();
// fout->Close();
}
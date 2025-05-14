//generate the scale factor between the weighted average to the average of each LB block in <ADC> 
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"
std::vector<double> MyIntegralAndBins(TH1D* h1, int bin_lo, int bin_hi);
void parameters(){
    std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp";
    //std::vector<string> triggers = {"/","/xor/","/minbias/"};
    //std::vector<string> trigger_name = {"AND","XOR","minbias"};
    std::vector<string> triggers = {"/","/xor/"};
    std::vector<string> trigger_name = {"AND", "XOR"};
    
    for(int i=0; i<trigger_name.size(); i++){
        TFile *input = new TFile(Form("%s%shistograms.Had2Cut.root",base.c_str(), triggers.at(i).c_str()),"read");
        //TFile *input = new TFile(Form("%s%shistograms.root",base.c_str(), triggers.at(i).c_str()),"read");
        TFile *output = new TFile(Form("%s%smean_block1.root",base.c_str(), triggers.at(i).c_str()),"recreate");
        for(int side =0; side<2; side++){
            for(int mod =0; mod <4; mod++){
                TH1D *h_mean = new TH1D(Form("h_mean_side%i_mod%i",side,mod),";LB bin; #LT ADC#GT",Bins::NLB,0,Bins::NLB);
                for(int ilb =0; ilb < Bins::NLB; ilb ++){
                    TH2D *h2 = (TH2D*)input->Get(Form("h%i_%i",side,mod))->Clone(Form("clone_h%i_%i_%s",side,mod,triggers.at(i).c_str()));
                    int bin_low = h2->GetXaxis()->FindBin(Bins::LB_LO[ilb]);
                    int bin_high = h2->GetXaxis()->FindBin(Bins::LB_HI[ilb]);
                    TH1D *h2_proj = h2->ProjectionY(Form("%i%i%i",side,mod,ilb),bin_low,bin_high);
                    double mean = h2_proj->GetMean();
                    h_mean->SetBinContent(ilb+1, mean);
                }
                output->cd();
                h_mean->Write();
            }
        }
        output->Close();
    }
}

//getting list of LB when the module is overflow
// 'BRL' stands for "bad running list" 
//Generating txt file with BRL
TH1D *h1[2][4];
TH1D *h_GRL[2][4];
TH2D *h2[2][4];
void clean_data(){
    TFile *fout = new TFile("BRL.root", "RECREATE");
    TFile *finput1 = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/Amp/histograms.root","read");
    TH1D *h_lb[2][4];
    for (int side=0; side<2; side++){
        for(int mod =0; mod<4; mod++){
            h_lb[side][mod] = new TH1D(Form("h_lb_%i%i",side,mod),";LB;",3800,1000,4800);
            h_GRL[side][mod] = new TH1D(Form("h_GRL%i_%i",side,mod),";LB;",3800,1000,4800);

        }
    }
    
    std::vector<pair<int,int>> slices;
    slices.push_back(make_pair(1000,1800));
    slices.push_back(make_pair(1801,1900));
    slices.push_back(make_pair(1901,2200));
    slices.push_back(make_pair(2201,2800));
    slices.push_back(make_pair(2801,3600));
    slices.push_back(make_pair(3601,4800));

    for(int sideAna=0; sideAna<2; sideAna++){
        for(int modAna =0; modAna<4; modAna++){
        h2[sideAna][modAna] = (TH2D*)finput1->Get(Form("h%i_%i",sideAna, modAna))->Clone("h01");
        h1[sideAna][modAna] = (TH1D*)h2[sideAna][modAna]->ProfileX("h01profile");
        // Create a cleaned histogram
        
       
        for(auto s: slices){
            TH1D *h_clean0 = (TH1D *)h1[sideAna][modAna]->Clone(Form("h_clean%i%i_slice",sideAna,modAna));
            int rangeLow = s.first;
            int rangeHigh = s.second;
        
            // Compute the mean and standard deviation in the fit range
            double sum = 0, sumSq = 0;
            int nBins = 0;

            for (int bin = h1[sideAna][modAna]->FindBin(rangeLow); bin <= h1[sideAna][modAna]->FindBin(rangeHigh); ++bin) {
                double content = h1[sideAna][modAna]->GetBinContent(bin);
                if(content ==0) continue;// dont take into account zeros
                sum += content;
                sumSq += content * content;
                nBins++;
            }

            double mean = sum / nBins;
            double variance = (sumSq / nBins) - (mean * mean);
            double sigma = sqrt(variance);
            cout << mean << " " << sigma <<endl;

            // Loop through bins and save events beyond 3sigma
            for (int bin = h1[sideAna][modAna]->FindBin(rangeLow); bin <= h1[sideAna][modAna]->FindBin(rangeHigh); bin++) {
                double content = h1[sideAna][modAna]->GetBinContent(bin);
                if(content == 0) continue;
                if (fabs(content - mean) > 3 * sigma) {
                    h_GRL[sideAna][modAna]->SetBinContent(bin, 1); 
                    h_GRL[sideAna][modAna]->SetBinError(bin, 0);   
                }
            }
        }
        fout->cd();
        h_GRL[sideAna][modAna]->Write(); 
    }
}

fout->Close();

std::cout << "Generating txt file of bad lumiblocks" << std::endl;
std::ofstream outFile("LBlist.txt");
std::set<int> uniqueBins;

for (int sideAna = 0; sideAna < 2; sideAna++) {
    for (int modAna = 0; modAna < 4; modAna++) {
        // Loop over bins and add bin numbers with value 1 to the set
        for (int bin = 1; bin <= h_GRL[sideAna][modAna]->GetNbinsX(); bin++) {
            if (h_GRL[sideAna][modAna]->GetBinContent(bin) == 1) {
                int lb = h_GRL[sideAna][modAna]->GetBinLowEdge(bin);
                uniqueBins.insert(lb);
            }
        }
    }
}

// Write unique bin numbers to the file
for (const auto& lb : uniqueBins) {
    outFile << lb << std::endl;
}
 
}
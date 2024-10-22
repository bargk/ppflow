std::vector<float> no_booster = {0.54, 1.00, 0.94, 0.79,1.47,1.02,0.87,0.54}; //first 4 lower bits is side C

void initHistos();
int triggerIndex(std::vector<bool> trigger);
bool passTrigger(std::vector<bool> trigger);
float sumZdc(int side,unsigned int module_mask, float module_amplitude[]);
bool isBitSet(int x, int s);
void SaveHistos();
TVectorD* ConcatenateTVectorD(const TVectorD* vec1, const TVectorD* vec2);

//histograms
TH1D* h0[2]; 
TH1D* h1[2]; 
TH1D* h_cut[2]; 
TH1D* h_module[8]; 
TH2D* h_module_corr[4];
TH2D* h0_corr; 
TH2D* h1_corr;
TFile *output; 

void LG_graphs(float number, bool fine_tune = false){
    TVectorD *gains81;
    TVectorD *gains12;
    TVectorD *weights;
    std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/%.1fsigma",number);
    std::string path = "/gpfs0/citron/users/bargl/ZDC/user.bglik.data23_hi.00463315.calibration_ZDCCalib.merge.AOD.c1535_m2248.ANALYSIS_EXT0/";

    //load zdc weights
    if(fine_tune){
        TFile *file12 = new TFile(Form("%s/zdcWeights_side0_fineTune.root",base.c_str()),"READ");
        TFile *file81 = new TFile(Form("%s/zdcWeights_side1_fineTune.root",base.c_str()),"READ");
        gains81 = (TVectorD*)file81->Get("gains_avg");
        gains12 = (TVectorD*)file12->Get("gains_avg");
        weights = ConcatenateTVectorD(gains12,gains81);
    }
    else{
        TFile *file12 = new TFile(Form("%s/zdcWeights_side0.root",base.c_str()),"READ");
        TFile *file81 = new TFile(Form("%s/zdcWeights_side1.root",base.c_str()),"READ");
        gains81 = (TVectorD*)file81->Get("gains_avg");
        gains12 = (TVectorD*)file12->Get("gains_avg");
        weights = ConcatenateTVectorD(gains12,gains81);
    }



    char name[100];
    char name1[100];
    sprintf(name,"%s/zdc_calibrated.root", base.c_str());
    if(fine_tune){
    sprintf(name,"%s/zdc_calibrated_fineTune.root", base.c_str());
    }
    output = new TFile(name,"recreate");
    initHistos();

     //load files
    TChain *tc = new TChain("zdcTree");
    tc->Add(Form("%s*",path.c_str()),0);
    //tc->Add(Form("%suser.bglik.41707724.EXT0._000005.ZDCNT.root",path.c_str()),0);
    if (!tc || tc->IsZombie() || !tc->GetNbranches()) {
    std::cerr << "Error opening file" << std::endl;
    exit(-1);
    }

    //load branches
    TTreeReader myreader(tc);
    TTreeReaderValue<unsigned int> lumiblock(myreader,"lumiBlock");
    TTreeReaderValue<bool> HLT_noalg_ZDCPEB_L1ZDC_C(myreader,"HLT_noalg_ZDCPEB_L1ZDC_C");
    TTreeReaderValue<bool> HLT_noalg_ZDCPEB_L1ZDC_A(myreader,"HLT_noalg_ZDCPEB_L1ZDC_A");
    TTreeReaderValue<bool> HLT_noalg_ZDCPEB_L1ZDC_A_C(myreader,"HLT_noalg_ZDCPEB_L1ZDC_A_C");
    TTreeReaderValue<bool> HLT_noalg_ZDCPEB_L1ZDC_OR(myreader,"HLT_noalg_ZDCPEB_L1ZDC_OR");
    TTreeReaderValue<float> ps_HLT_noalg_ZDCPEB_L1ZDC_A_C(myreader,"ps_HLT_noalg_ZDCPEB_L1ZDC_A_C");
    TTreeReaderValue<float> ps_HLT_noalg_ZDCPEB_L1ZDC_OR(myreader,"ps_HLT_noalg_ZDCPEB_L1ZDC_OR");
    TTreeReaderValue<float> ps_HLT_noalg_ZDCPEB_L1ZDC_C(myreader,"ps_HLT_noalg_ZDCPEB_L1ZDC_C");
    TTreeReaderValue<float> ps_HLT_noalg_ZDCPEB_L1ZDC_A(myreader,"ps_HLT_noalg_ZDCPEB_L1ZDC_A");
    TTreeReaderArray<float> Amp(myreader, "zdc_ZdcAmp");
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp");

    int nentries= myreader.GetEntries();
    int istart=0,iend= nentries; //by default loop over whole dataset
    cout << "Total events : "<< nentries << endl;

    while (myreader.Next()){
        int ievt = myreader.GetCurrentEntry();
        if(ievt%100000==0) cout<<"proccesed "<<ievt<<" / "<<iend<<" events "<<" "<<tc->GetFile()->GetName()<<endl;
        if(ievt > iend) break;
        std::vector<bool> m_trig;
        std::vector<float> m_trig_ps;
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_C);            m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_C);  
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A);            m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A);
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A_C);          m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A_C);
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_OR);           m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_OR);

        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        int trig_index = triggerIndex(m_trig); // retrive the index of the relevant trigger
        float prescale = m_trig_ps.at(trig_index);
        float module_amp[8]; 
        for(int i=0; i<8; i++){
            module_amp[i] = (*weights)[i]*ModAmp[i]/no_booster.at(i); //now with energy
            h_module[i]->Fill(module_amp[i], prescale);
            if(i<4) h_module_corr[i]->Fill(module_amp[i],module_amp[i+4], prescale); //x: side C, y: side A
        }
        h1[0]->Fill(sumZdc(0, (*BitMask), module_amp), prescale);
        h1[1]->Fill(sumZdc(1, (*BitMask), module_amp), prescale);

        for(int side =0; side<2; side++){
            float total_sum =0;

            // opposite selection
            if(side ==0 && (*HLT_noalg_ZDCPEB_L1ZDC_A)){
                total_sum = sumZdc(side, (*BitMask), module_amp);
                prescale = *ps_HLT_noalg_ZDCPEB_L1ZDC_A;
                h0[side]->Fill(total_sum,prescale);
            }
            else if(side ==1 && (*HLT_noalg_ZDCPEB_L1ZDC_C)){
                total_sum = sumZdc(side, (*BitMask), module_amp);
                prescale = *ps_HLT_noalg_ZDCPEB_L1ZDC_C;
                h0[side]->Fill(total_sum,prescale);
            }
            else {continue;}
        }
    }
    if(fine_tune){
        output->cd();
        output->Write();
        output->Close();
        cout << "Saving Finished" << endl;
    }
        
    if(!fine_tune){ //create sigma cut for fine tuning in case fine tune files not exist already
        //doing sigma cut
        float chi2_arr[2];
        float ndf_arr[2];
        float constant_arr[2];
        float mean_arr[2]; 
        float sigma_arr[2];
        float constErr_arr[2];
        float meanErr_arr[2];
        float sigmaErr_arr[2];

        std::cout << "Performing sigma cut!" <<std::endl;
        cout <<h0[0]->GetMaximumBin()<<endl;
        for (int side =0; side<2; side++){
            int maximum = h0[side]->GetMaximumBin();
            double low_lim = 0.5*h0[side]->GetXaxis()->GetBinCenter(maximum);
            double high_lim = 1.5*h0[side]->GetXaxis()->GetBinCenter(maximum);

            // Perform iterative fit
            TF1 *fit = new TF1("fit","gaus");
            fit->SetLineColor(kRed);
            fit->SetLineWidth(5);
            fit->SetRange(low_lim,high_lim);
            h0[side]->Fit("fit","Rqn");
            double mean = fit->GetParameter(1);
            double sigma = fit->GetParameter(2);
            double mu_first = mean;
            double mu_second = 0; 
            int stop =0;
            while (true){
                fit->SetRange(mu_first -number*sigma, mu_first + number*sigma);
                h0[side]->Fit("fit", "Rqn");
                mu_second = fit->GetParameter(1);
                cout << fabs(mu_first -mu_second) << std::endl;
                if(fabs(mu_first - mu_second)<0.1){
                    break;
                }
                mu_first = mu_second;
                sigma = fit->GetParameter(2);
                stop +=1;
                if( stop == 20){
                    break;
                }
            }
            //update parameters
            fit->SetRange(mu_first -number*sigma, mu_first +number*sigma);
            h0[side]->Fit("fit","R");
            chi2_arr[side] = fit->GetChisquare();
            ndf_arr[side] = fit->GetNDF();
            constant_arr[side] = fit->GetParameter(0);
            mean_arr[side]   = fit->GetParameter(1);
            sigma_arr[side]  = fit->GetParameter(2);
            constErr_arr[side] = fit->GetParError(0);
            meanErr_arr[side] = fit->GetParError(1);
            sigmaErr_arr[side]= fit->GetParError(2);
        }

    output->cd();
    output->Write();
    output->Close();
    cout << "Saving Finished" << endl;
    

    //Now create histograms with sigma cut
    std::cout << "Generating matrix elements" << std::endl;
    float M0; float M1; float M2; float M3;

    TFile *file = new TFile(Form("%s/matrix_elements_fineTune.root",base.c_str()), "RECREATE");
    TTree *tree_a = new TTree("mytree_a","tree_title");
    TTree *tree_c = new TTree("mytree_c","tree_title");
    tree_a->Branch("M0", &M0, "M0/F");
    tree_a->Branch("M1", &M1, "M1/F");
    tree_a->Branch("M2", &M2, "M2/F");
    tree_a->Branch("M3", &M3, "M3/F");
    tree_c->Branch("M0", &M0, "M0/F");
    tree_c->Branch("M1", &M1, "M1/F");
    tree_c->Branch("M2", &M2, "M2/F");
    tree_c->Branch("M3", &M3, "M3/F");

    myreader.SetEntry(0);
    while (myreader.Next()){
        float module_amp[8];
        float total_sum=0;
        for(int i=0; i<8; i++){
                module_amp[i] = ModAmp[i]/no_booster.at(i);
            }
            for(int side =0; side<2; side++){
                // opposite selection
                if(side ==0 && (*HLT_noalg_ZDCPEB_L1ZDC_A)){
                    total_sum = sumZdc(side, (*BitMask), module_amp);
                    if (total_sum < mean_arr[side] - number * sigma_arr[side] || total_sum > mean_arr[side] + number * sigma_arr[side]) continue; // do sigma selection
                    M0 = module_amp[0];
                    M1 = module_amp[1];
                    M2 = module_amp[2];
                    M3 = module_amp[3];
                    tree_c->Fill();
                }
                else if(side ==1 && (*HLT_noalg_ZDCPEB_L1ZDC_C)){
                    total_sum = sumZdc(side, (*BitMask), module_amp);
                    if (total_sum < mean_arr[side] - number * sigma_arr[side] || total_sum > mean_arr[side] + number * sigma_arr[side]) continue; // do sigma selection
                    M0 = module_amp[0+4];
                    M1 = module_amp[1+4];
                    M2 = module_amp[2+4];
                    M3 = module_amp[3+4];
                    tree_a->Fill();
                }
            }
            
    }
    file->Write();
    file->Close();
    }

    char letter[100];
    std::cout << "finished" << std::endl;
    sprintf(letter,"kill -9 %d",gSystem->GetPid());
    std::cout<<letter<<std::endl;
    gSystem->Exec(letter);

    
}

void initHistos(){
    std::cout << "initHistos()" << endl;
    h0[0] = new TH1D("h_c_opposite",";#sum [GeV];Counts",500,500,15000); 
    h0[1] = new TH1D("h_a_opposite",";#sum [GeV];Counts",500,500,15000);
    h_cut[0] = new TH1D("h_c_opposite_cut",";#sum [GeV];Counts",500,500,10000); 
    h_cut[1] = new TH1D("h_a_opposite_cut",";#sum [GeV];Counts",500,500,10000);
    h1[0] = new TH1D("h_c",";#sum [GeV];Counts",500,500,15000);
    h1[1] = new TH1D("h_a",";#sum [GeV];Counts",500,500,15000);
    h0_corr = new TH2D("h_corr_opposite",";#sum [GeV] (side C);#sum [GeV] (side A)",500,100,10000,500,100,10000);
    h1_corr = new TH2D("h_corr",";#sum [GeV] (side C);#sum [GeV] (side A)",500,100,10000,500,100,10000);

    //histogram for each module
    for(int i=0; i<8; i++){
        h_module[i] = new TH1D(Form("module_%i",i),";#sum [GeV];Counts",400,0,4000);
    }
    //histograms for correlation
    for(int i=0; i<4; i++){
        h_module_corr[i] = new TH2D(Form("module_%i_corr",i),";#sum [GeV] (side C);#sum [GeV] (side A)",400,0,4000,400,0,4000);
    }
 
}

 bool passTrigger(std::vector<bool> trigger){
        bool pass = false;
        for(const auto &itrg : trigger){
            if(itrg){
            pass = true;
            break;
            }
        }
        return pass;
 }

int triggerIndex(std::vector<bool> trigger){
    int index =-1;
    for(int i=0; i<trigger.size(); i++){
        if(trigger.at(i)){
            index = i;
            break;
        }
    }
    if(index == -1){
        std::cout << "Error : trigger index is not valid" << std::endl;
        throw std::exception();
    }
    return index;
}

float sumZdc(int side,unsigned int module_mask, float module_amplitude[]){
    float sum_zdc = 0;
    if (side ==0){
        for(int i=0; i<4; i++){
            if(isBitSet(module_mask,i)) sum_zdc+= module_amplitude[i];
        }
    }
        else if(side == 1){
            for(int i=4; i<8; i++){
                if(isBitSet(module_mask,i)) sum_zdc+= module_amplitude[i];
        }
        }
        else{
            std::cerr << "Error: Unknown side value" << std::endl;
            exit(-1);
        }
        return sum_zdc;
    }

    bool isBitSet(int x, int s){
        int mask = x >> s;
        return mask % 2;
    }

    TVectorD* ConcatenateTVectorD(const TVectorD* vec1, const TVectorD* vec2) {
    // Get the sizes of the two vectors
    int size1 = vec1->GetNrows();
    int size2 = vec2->GetNrows();
    
    // Create a new TVectorD to hold the concatenated result
    TVectorD* concatenated = new TVectorD(size1 + size2);

    // Copy elements from vec1 into the new vector
    for (int i = 0; i < size1; ++i) {
        (*concatenated)[i] = (*vec1)[i];
    }

    // Copy elements from vec2 into the new vector, starting after vec1's elements
    for (int i = 0; i < size2; ++i) {
        (*concatenated)[i + size1] = (*vec2)[i];
    }

    return concatenated;  // Return the pointer to the concatenated vector
}
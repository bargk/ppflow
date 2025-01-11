std::vector<float> no_booster = {0.54, 1.00, 0.94, 0.79,1.47,1.02,0.87,0.54}; //first 4 lower bits is side C

void initHistos();
int triggerIndex(std::vector<bool> trigger);
bool passTrigger(std::vector<bool> trigger);
float sumZdc(int side,unsigned int module_mask, float module_amplitude[]);
bool isBitSet(int x, int s);
void SaveHistos();
TVectorD* ConcatenateTVectorD(const TVectorD* vec1, const TVectorD* vec2);
TVectorD* MultiplyMultipleTVectorD(const std::vector<TVectorD*>& vectors);

//histograms
TH1D* h0[2]; 
TH1D* h1[2]; 
TH1D* h_cut[2]; 
TH1D* h_module[8]; 
TH2D* h_module_corr[4];
TH2D* h0_corr; 
TH2D* h1_corr;
TFile *output; 

void LG_graphs(int itr){
    std::vector<TVectorD*> weights;
    TVectorD* multiplied_weights = new TVectorD(8);
    std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide");
    std::string path = "/gpfs0/citron/users/bargl/ZDC/user.bglik.data23_hi.00463315.calibration_ZDCCalib.merge.AOD.c1535_m2248.ANALYSIS_EXT0/";

    // Loop over all iterations
    for (int idx = 1; idx <= itr; idx++) {
        // Open files for the current iteration
        TFile* file12 = new TFile(Form("%s/zdcWeights_side0_itr%i.root", base.c_str(), idx), "READ");
        TFile* file81 = new TFile(Form("%s/zdcWeights_side1_itr%i.root", base.c_str(), idx), "READ");

        // Retrieve gain vectors for both sides
        TVectorD* gains12 = (TVectorD*)file12->Get("gains_avg");
        TVectorD* gains81 = (TVectorD*)file81->Get("gains_avg");

        // Concatenate vectors for the current iteration
        TVectorD* weights_itr = ConcatenateTVectorD(gains12, gains81);
        weights.push_back(weights_itr);

        // Clean up
        file12->Close();
        file81->Close();
    }

    // Perform element-wise multiplication across all iterations
    multiplied_weights = MultiplyMultipleTVectorD(weights);
    if (multiplied_weights) {
        std::cout << "Result of element-wise multiplication across " << itr << " iterations:" << std::endl;
        multiplied_weights->Print();
    }

    //save the weights for each side
    TVectorD* gains12 = new TVectorD(4);
    TVectorD* gains81 = new TVectorD(4);
    for(int i =0; i<4; i++){
        (*gains12)[i] = (*multiplied_weights)[i];
        (*gains81)[i] = (*multiplied_weights)[i+4];
    }
    TFile* fweights12 = new TFile(Form("%s/zdcWeights_side0.root", base.c_str()), "recreate");
    TFile* fweights81 = new TFile(Form("%s/zdcWeights_side1.root", base.c_str()), "recreate");
    fweights12->cd();
    gains12->Write("gains");
    fweights12->Close();
    fweights81->cd();
    gains81->Write("gains");
    fweights81->Close();

    char name[100];
    char name1[100];
    sprintf(name,"%s/zdc_calibrated.root", base.c_str());
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
    int istart=0,iend= 0.01*nentries; //by default loop over whole dataset
    cout << "Total events : "<< nentries << endl;

    while (myreader.Next()){
        int ievt = myreader.GetCurrentEntry();
        if(ievt%100000==0) cout<<"proccesed "<<ievt<<" / "<<iend<<" events "<<" "<<tc->GetFile()->GetName()<<endl;
        if(ievt > iend) break;
        std::vector<bool> m_trig;
        std::vector<float> m_trig_ps;
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_C);            m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_C);  
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A);            m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A);
        // m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A_C);          m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A_C);
        // m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_OR);           m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_OR);

        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        int trig_index = triggerIndex(m_trig); // retrive the index of the relevant trigger
        float prescale = m_trig_ps.at(trig_index);
        float module_amp[8]; 
        for(int i=0; i<8; i++){
            module_amp[i] = (*multiplied_weights)[i]*ModAmp[i]/no_booster.at(i); //now with energy
            h_module[i]->Fill(module_amp[i]);
        }
        if((*HLT_noalg_ZDCPEB_L1ZDC_C)){
            h1[0]->Fill(sumZdc(0, (*BitMask), module_amp)); // oposite side
            h0[1]->Fill(sumZdc(1, (*BitMask), module_amp)); // same side
        }
        if((*HLT_noalg_ZDCPEB_L1ZDC_A)){
           h1[1]->Fill(sumZdc(1, (*BitMask), module_amp));
           h0[0]->Fill(sumZdc(0, (*BitMask), module_amp));
        }
    }
    output->Write();
    output->Close();
    char letter[100];
    std::cout << "finished" << std::endl;
    sprintf(letter,"kill -9 %d",gSystem->GetPid());
    std::cout<<letter<<std::endl;
    gSystem->Exec(letter);

    
}

void initHistos(){
    std::cout << "initHistos()" << endl;
    double LimLow =0;
    double LimHigh =10000;
    h0[0] = new TH1D("h_c_opposite",";#sum [GeV];Counts",200,LimLow,LimHigh); 
    h0[1] = new TH1D("h_a_opposite",";#sum [GeV];Counts",200,LimLow,LimHigh);
    h_cut[0] = new TH1D("h_c_opposite_cut",";#sum [GeV];Counts",200,LimLow,LimHigh); 
    h_cut[1] = new TH1D("h_a_opposite_cut",";#sum [GeV];Counts",200,LimLow,LimHigh);
    h1[0] = new TH1D("h_c",";#sum [GeV];Counts",200,LimLow,LimHigh);
    h1[1] = new TH1D("h_a",";#sum [GeV];Counts",200,LimLow,LimHigh);
    // h0_corr = new TH2D("h_corr_opposite",";#sum GeV (side C);#sum GeV (side A)",500,LimLow,LimHigh,500,LimLow,LimHigh);
    // h1_corr = new TH2D("h_corr",";#sum GeV (side C);#sum GeV (side A)",500,LimLow,LimHigh,500,LimLow,LimHigh);

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

// Function to multiply multiple TVectorD objects element-wise
TVectorD* MultiplyMultipleTVectorD(const std::vector<TVectorD*>& vectors) {
    if (vectors.empty()) {
        std::cerr << "Error: No vectors provided for multiplication!" << std::endl;
        return nullptr;
    }

    // Check that all vectors have the same size
    size_t size = vectors[0]->GetNoElements();
    for (const auto& vec : vectors) {
        if (vec->GetNoElements() != size) {
            std::cerr << "Error: Vectors have different sizes!" << std::endl;
            return nullptr;
        }
    }

    // Perform element-wise multiplication
    TVectorD* result = new TVectorD(size);
    for (size_t i = 0; i < size; ++i) {
        double product = 1.0;
        for (const auto& vec : vectors) {
            product *= (*vec)[i];
        }
        (*result)[i] = product;
    }

    return result;
}
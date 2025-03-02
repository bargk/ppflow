    #include "calibration.h"

    void initHistos();
    int triggerIndex(std::vector<bool> trigger);
    bool passTrigger(std::vector<bool> trigger);
    float sumZdc(int side,unsigned int module_mask, float module_amplitude[]);
    bool isBitSet(int x, int s);
    void SaveHistos();
    TVectorD* ConcatenateTVectorD(const TVectorD* vec1, const TVectorD* vec2);
    TVectorD* MultiplyMultipleTVectorD(const std::vector<TVectorD*>& vectors);

    void data_gen(float number, int itr){

    gSystem->Exec(Form("mkdir -p %s",base.c_str()));
    char name[100];
    char name1[100];
    sprintf(name,"%s/zdc_uncalib_itr%d.root", base.c_str(),itr);
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

    //load weights from previous iteration
    TVectorD *weights = new TVectorD(8);
    TVectorD *initial_weight_12 = new TVectorD(4);
    TVectorD *initial_weight_81 = new TVectorD(4);
    int itr_prev = itr -1;
   if(itr_prev > 0){
    TFile* file12 = new TFile(Form("%s/zdcWeights_side0_itr%i.root", base.c_str(), itr_prev), "READ");
    TFile* file81 = new TFile(Form("%s/zdcWeights_side1_itr%i.root", base.c_str(), itr_prev), "READ");

    // Retrieve gain vectors for both sides
    TVectorD* gains12 = (TVectorD*)file12->Get("gains_avg");
    TVectorD* gains81 = (TVectorD*)file81->Get("gains_avg");

    // Concatenate vectors for the current iteration
    weights = ConcatenateTVectorD(gains12, gains81);

    // Clean up
    file12->Close();
    file81->Close();
   }
    else{ //initial guess is with hv gains
        for (int i = 0; i < 8; ++i) {
            //(*multiplied_weights)[i] = hv_gain.at(i)/no_booster.at(i);
            (*weights)[i] = hv_gain.at(i)/no_booster.at(i);
        }
    }
        std::cout << "Result of element-wise multiplication across" << itr << " iterations:" << std::endl;
        weights->Print();




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
        //m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A_C);          m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A_C);
        //m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_OR);           m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_OR);

        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        int trig_index = triggerIndex(m_trig); // retrive the index of the relevant trigger
        float prescale = m_trig_ps.at(trig_index);
        int lumiBlock = *lumiblock;
        if(!((lumiBlock > 268 && lumiBlock < 675))) continue; //stable beams
        float module_amp[8];
        for(int i=0; i<8; i++){
            module_amp[i] = (*weights)[i]*ModAmp[i];
        }

        if((*HLT_noalg_ZDCPEB_L1ZDC_C)){
            h1[0]->Fill(sumZdc(0, (*BitMask), module_amp)); // oposite side
            h0[1]->Fill(sumZdc(1, (*BitMask), module_amp)); // same side
            for(int i=0; i<4; i++){
                if(isBitSet((*BitMask),i)) h_module[i]->Fill(lumiBlock,module_amp[i]);
            }
        }
        if((*HLT_noalg_ZDCPEB_L1ZDC_A)){
           h1[1]->Fill(sumZdc(1, (*BitMask), module_amp));
           h0[0]->Fill(sumZdc(0, (*BitMask), module_amp));
           for(int i=4; i<8; i++){
                if(isBitSet((*BitMask),i)) h_module[i]->Fill(lumiBlock,module_amp[i]);
            }
        }
    }

    //doing sigma cut
    std::cout << "Performing sigma cut!" <<std::endl;

    for (int side =0; side<2; side++){
        TH1D *htmp = (TH1D*)h1[side]->Clone(Form("%i",side));
        int bin = htmp->FindBin(300);
        htmp->GetXaxis()->SetRange(bin,htmp->GetNbinsX());
        int maximum = htmp->GetMaximumBin();
        double low_lim = 0.8*h1[side]->GetXaxis()->GetBinCenter(maximum);
        double high_lim = 1.2*h1[side]->GetXaxis()->GetBinCenter(maximum);
        cout << "Maximum bin : " << maximum<<endl;
        // Perform iterative fit
        TF1 *fit = new TF1("fit","gaus");
        fit->SetLineColor(kRed);
        fit->SetLineWidth(5);
        fit->SetRange(low_lim,high_lim);
        h1[side]->Fit("fit","Rqn");

        double mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        double mu_first = mean;
        double mu_second = 0; 
        int stop =0;
        while (true){
            fit->SetRange(mu_first -number*sigma, mu_first + number*sigma);
            h1[side]->Fit("fit", "Rqn");
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
        h1[side]->Fit("fit","R");
        chi2_arr[side] = fit->GetChisquare();
        ndf_arr[side] = fit->GetNDF();
        constant_arr[side] = fit->GetParameter(0);
        mean_arr[side]   = fit->GetParameter(1);
        sigma_arr[side]  = fit->GetParameter(2);
        constErr_arr[side] = fit->GetParError(0);
        meanErr_arr[side] = fit->GetParError(1);
        sigmaErr_arr[side]= fit->GetParError(2);
    }
    //Now create histograms with sigma cut
    std::cout << "Generating matrix elements" << std::endl;
    float M0; float M1; float M2; float M3;

    TFile *file = new TFile(Form("%s/matrix_elements_itr%i.root",base.c_str(),itr), "RECREATE");
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
        int ievt = myreader.GetCurrentEntry();
        if(ievt%100000==0) cout<<"proccesed "<<ievt<<" / "<<iend<<" events "<<" "<<tc->GetFile()->GetName()<<endl;
        if(ievt > iend) break;
        std::vector<bool> m_trig;
        std::vector<float> m_trig_ps;
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_C);            m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_C);  
        m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A);            m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A);
        //m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_A_C);          m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_A_C);
        //m_trig.push_back(*HLT_noalg_ZDCPEB_L1ZDC_OR);           m_trig_ps.push_back(*ps_HLT_noalg_ZDCPEB_L1ZDC_OR);

        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        int lumiBlock = *lumiblock;
        if(!((lumiBlock > 268 && lumiBlock < 675))) continue; //stable beams
        float module_amp[8];
        float total_sum=0;
        for(int i=0; i<8; i++){
                module_amp[i] = (*weights)[i]*ModAmp[i];
            }
        if((*HLT_noalg_ZDCPEB_L1ZDC_C)){
            total_sum = sumZdc(0, (*BitMask), module_amp);
            if ( (total_sum < mean_arr[0] - number * sigma_arr[0] ) || (total_sum > (mean_arr[0] + number * sigma_arr[0]))) continue; // do sigma selection
            M0 = module_amp[0];
            M1 = module_amp[1];
            M2 = module_amp[2];
            M3 = module_amp[3];
            tree_c->Fill();
            h_cut[0]->Fill(total_sum);
        }

        if((*HLT_noalg_ZDCPEB_L1ZDC_A)){
            total_sum = sumZdc(1, (*BitMask), module_amp);
            if (( total_sum < (mean_arr[1] - number * sigma_arr[1]) ) || (total_sum > (mean_arr[1] + number * sigma_arr[1]))) continue; // do sigma selection
            M0 = module_amp[0+4];
            M1 = module_amp[1+4];
            M2 = module_amp[2+4];
            M3 = module_amp[3+4];
            tree_a->Fill();
            h_cut[1]->Fill(total_sum);
        }        
    }
    file->cd();
    file->Write();
    file->Close();
    output->cd();
    output->Write();
    output->Close();
    cout << "Saving Finished" << endl;

    char letter[100];
    std::cout << "finished" << std::endl;
    sprintf(letter,"kill -9 %d",gSystem->GetPid());
    std::cout<<letter<<std::endl;
    gSystem->Exec(letter);
    }

    void SaveHistos(){
        output->cd();
        output->Write();
        cout << "Saving Finished" << endl;
    }

void initHistos(){
    std::cout << "initHistos()" << endl;
    double LimLow =0;
    double LimHigh =10000;
    h0[0] = new TH1D("h_c_opposite",";#sum ADC;Counts",200,LimLow,LimHigh); 
    h0[1] = new TH1D("h_a_opposite",";#sum ADC;Counts",200,LimLow,LimHigh);
    h_cut[0] = new TH1D("h_c_cut",";#sum ADC;Counts",200,LimLow,LimHigh); 
    h_cut[1] = new TH1D("h_a_cut",";#sum ADC;Counts",200,LimLow,LimHigh);
    h1[0] = new TH1D("h_c",";#sum ADC;Counts",200,LimLow,LimHigh);
    h1[1] = new TH1D("h_a",";#sum ADC;Counts",200,LimLow,LimHigh);
    // h0_corr = new TH2D("h_corr_opposite",";#sum ADC (side C);#sum ADC (side A)",500,LimLow,LimHigh,500,LimLow,LimHigh);
    // h1_corr = new TH2D("h_corr",";#sum ADC (side C);#sum ADC (side A)",500,LimLow,LimHigh,500,LimLow,LimHigh);

    //histogram for each module
    for(int i=0; i<8; i++){
        h_module[i] = new TH2D(Form("module_%i",i),";LB;ADC",500,200,700,500,0,15000);
    }
    //histograms for correlation
    for(int i=0; i<4; i++){
        h_module_corr[i] = new TH2D(Form("module_%i_corr",i),";#sum ADC (side C);#sum ADC (side A)",400,0,4000,400,0,4000);
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
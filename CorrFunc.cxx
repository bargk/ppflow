#include "Event.h"
#include "Defs.h"
#include "bins.h"
//#define check




void CorrFunc(const int a ,const char* fileList, bool minbias1 = 0){
    minbias = minbias1;
    sprintf(output_name,"histograms_%d.root",a);
    if(minbias){
        cout << "Working on minbias triggers!" << endl;
    }
    else{
        cout << "working on ZDC triggers!" << endl;
    }
    // Split the file paths into a vector
    std::vector<std::string> files;
    std::istringstream ss(fileList);
    std::string file;
    while (std::getline(ss, file, ',')) {
        files.push_back(file);
    }
    
    if (files.empty()) {
        std::cerr << "No files provided!" << std::endl;
        return;
    }
    // loading data
    fChain = new TChain("zdcTree");
    char base[200] ="/gpfs0/citron/users/bargl/ZDC/lhcf22/user.bglik.data22_13p6TeV.00435229.ANALAYSIS_1_noLhcf_EXT0.566478889.566478889";
    // Open the files
    for (const auto file : files) {
        fChain->Add(Form("%s/%s",base,file.c_str()),0);
        if (!fChain || fChain->IsZombie()) {
            std::cerr << "Error opening file: " << file <<" Check files.txt content"<< std::endl;
            //continue;
        }
        std::cout << "Opened file: " << file << std::endl;
    }

    //fChain->Add(Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/user.steinber.data22_13p6TeV.00435229.physics_MinBias.merge.AOD.r14470_p5587.4zdc_EXT0/%s",s.c_str()),0);
    //fChain->Add(Form("%s/*.root",path.c_str()),0);
    //fChain->Add(Form("%s/user.steinber.39675801.EXT0._000981.ZDCNT.root",path.c_str())); // zdc
    if (!fChain || fChain->IsZombie() || !fChain->GetNbranches()) {
      std::cerr << "Error opening file" << std::endl;
      exit(-1);
    }

    
    
    load_weights();
    // getting branches
    TTreeReader myreader(fChain);
    TTreeReaderValue<std::vector<float>> trk_pt(myreader, "trk_pt"); // units of MeV
    TTreeReaderValue<std::vector<float>> trk_eta(myreader, "trk_eta");
    TTreeReaderValue<std::vector<float>> trk_phi(myreader, "trk_phi");
    TTreeReaderValue<std::vector<signed char>> trk_charge(myreader, "trk_charge");
    TTreeReaderValue<std::vector<short>> trk_quality(myreader, "trk_quality");
    TTreeReaderValue<std::vector<float>> vtx_z(myreader, "vtx_z");
    TTreeReaderValue<unsigned int> ntrk_ptr(myreader, "ntrk");
    TTreeReaderValue<int> nvtx(myreader, "t_nvtx");
    TTreeReaderValue<unsigned int> lumiblock(myreader, "lumiBlock");
    //zdc triggers
    TTreeReaderValue<bool> HLT_L1_ZDC_XOR_E1_E3(myreader, "HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_L1_ZDC_XOR_E2(myreader, "HLT_noalg_L1ZDC_XOR_E2");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_A_AND_C(myreader, "HLT_noalg_L1ZDC_A_AND_C");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_XOR_E2(myreader,"HLT_mb_sptrk_L1ZDC_XOR_E2");
	TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_XOR_E1_E3(myreader,"HLT_mb_sptrk_L1ZDC_XOR_E1_E3");
	TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_A_AND_C(myreader, "HLT_mb_sptrk_L1ZDC_A_AND_C");
	TTreeReaderValue<bool> HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2(myreader,"HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2");
	TTreeReaderValue<bool> HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3(myreader, "HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3");
	TTreeReaderValue<bool> HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C(myreader, "HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C");
    TTreeReaderValue<float> ps_HLT_L1_ZDC_XOR_E1_E3(myreader, "ps_HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<float> ps_HLT_L1_ZDC_XOR_E2(myreader, "ps_HLT_noalg_L1ZDC_XOR_E2");
    TTreeReaderValue<float> ps_HLT_noalg_L1ZDC_A_AND_C(myreader, "ps_HLT_noalg_L1ZDC_A_AND_C");
    TTreeReaderValue<float> ps_HLT_mb_sptrk_L1ZDC_XOR_E2(myreader,"ps_HLT_mb_sptrk_L1ZDC_XOR_E2");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_L1ZDC_XOR_E1_E3(myreader,"ps_HLT_mb_sptrk_L1ZDC_XOR_E1_E3");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_L1ZDC_A_AND_C(myreader, "ps_HLT_mb_sptrk_L1ZDC_A_AND_C");
	TTreeReaderValue<float> ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2(myreader,"ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2");
	TTreeReaderValue<float> ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3(myreader, "ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3");
	TTreeReaderValue<float> ps_HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C(myreader, "ps_HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C");


    //minbias triggers
    TTreeReaderValue<bool> HLT_noalg_L1MBTS_1(myreader, "HLT_noalg_L1MBTS_1");
    TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_1(myreader, "HLT_noalg_mb_L1MBTS_1");
	TTreeReaderValue<bool> HLT_noalg_L1MBTS_1_1(myreader, "HLT_noalg_L1MBTS_1_1");
	TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_1_1(myreader, "HLT_noalg_mb_L1MBTS_1_1");
	TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_2(myreader, "HLT_noalg_mb_L1MBTS_2");
	TTreeReaderValue<bool> HLT_noalg_L1MBTS_A(myreader, "HLT_noalg_L1MBTS_A");
	TTreeReaderValue<bool> HLT_noalg_L1MBTS_C(myreader, "HLT_noalg_L1MBTS_C");
    TTreeReaderValue<bool> HLT_mb_sptrk_pt2_L1MBTS_2(myreader, "HLT_mb_sptrk_pt2_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_sptrk_pt4_L1MBTS_2(myreader, "HLT_mb_sptrk_pt4_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_sptrk_pt6_L1MBTS_2(myreader, "HLT_mb_sptrk_pt6_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_sptrk_pt8_L1MBTS_2(myreader, "HLT_mb_sptrk_pt8_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_1(myreader, "HLT_mb_mbts_L1MBTS_1");
	TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_1_1(myreader, "HLT_mb_mbts_L1MBTS_1_1");
	TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_2(myreader, "HLT_mb_mbts_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_mbts_L1RD0_FILLED(myreader, "HLT_mb_mbts_L1RD0_FILLED");
    TTreeReaderValue<bool> HLT_mb_sptrk_vetombts2in_L1RD0_FILLED(myreader, "HLT_mb_sptrk_vetombts2in_L1RD0_FILLED");
	TTreeReaderValue<bool> HLT_mb_mbts_all_L1MBTS_A(myreader, "HLT_mb_mbts_all_L1MBTS_A");
	TTreeReaderValue<bool> HLT_mb_mbts_all_L1MBTS_C(myreader, "HLT_mb_mbts_all_L1MBTS_C");
    TTreeReaderValue<float> ps_HLT_noalg_L1MBTS_1(myreader, "ps_HLT_noalg_L1MBTS_1");
    TTreeReaderValue<float> ps_HLT_noalg_mb_L1MBTS_1(myreader, "ps_HLT_noalg_mb_L1MBTS_1");
	TTreeReaderValue<float> ps_HLT_noalg_L1MBTS_1_1(myreader, "ps_HLT_noalg_L1MBTS_1_1");
	TTreeReaderValue<float> ps_HLT_noalg_mb_L1MBTS_1_1(myreader, "ps_HLT_noalg_mb_L1MBTS_1_1");
	TTreeReaderValue<float> ps_HLT_noalg_mb_L1MBTS_2(myreader, "ps_HLT_noalg_mb_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_noalg_L1MBTS_A(myreader, "ps_HLT_noalg_L1MBTS_A");
	TTreeReaderValue<float> ps_HLT_noalg_L1MBTS_C(myreader, "ps_HLT_noalg_L1MBTS_C");
    TTreeReaderValue<float> ps_HLT_mb_sptrk_pt2_L1MBTS_2(myreader, "ps_HLT_mb_sptrk_pt2_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_pt4_L1MBTS_2(myreader, "ps_HLT_mb_sptrk_pt4_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_pt6_L1MBTS_2(myreader, "ps_HLT_mb_sptrk_pt6_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_pt8_L1MBTS_2(myreader, "ps_HLT_mb_sptrk_pt8_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_mb_mbts_L1MBTS_1(myreader, "ps_HLT_mb_mbts_L1MBTS_1");
	TTreeReaderValue<float> ps_HLT_mb_mbts_L1MBTS_1_1(myreader, "ps_HLT_mb_mbts_L1MBTS_1_1");
	TTreeReaderValue<float> ps_HLT_mb_mbts_L1MBTS_2(myreader, "ps_HLT_mb_mbts_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_mb_mbts_L1RD0_FILLED(myreader, "ps_HLT_mb_mbts_L1RD0_FILLED");
    TTreeReaderValue<float> ps_HLT_mb_sptrk_vetombts2in_L1RD0_FILLED(myreader, "ps_HLT_mb_sptrk_vetombts2in_L1RD0_FILLED");
	TTreeReaderValue<float> ps_HLT_mb_mbts_all_L1MBTS_A(myreader, "ps_HLT_mb_mbts_all_L1MBTS_A");
	TTreeReaderValue<float> ps_HLT_mb_mbts_all_L1MBTS_C(myreader, "ps_HLT_mb_mbts_all_L1MBTS_C");
    
    //zdc
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp"); // First 4 lower bits corressponds to side C
    TTreeReaderArray<float> Amp(myreader, "zdc_ZdcAmp"); 
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");

    InitHistos();

    int nentries= myreader.GetEntries();
    int istart=0,iend= nentries; //by default loop over whole dataset
    cout << "Total events : "<< nentries << endl;
    
      while (myreader.Next()){
        //retrive trigger information
        std::vector<bool> m_trig;
        std::vector<float> m_trig_ps;
        if(minbias){
            m_trig.push_back(*HLT_noalg_L1MBTS_1);                                  m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_1);
            m_trig.push_back(*HLT_noalg_mb_L1MBTS_1);                               m_trig_ps.push_back(*ps_HLT_noalg_mb_L1MBTS_1);
            m_trig.push_back(*HLT_noalg_L1MBTS_1_1);                                m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_1_1);
            m_trig.push_back(*HLT_noalg_mb_L1MBTS_1_1);                             m_trig_ps.push_back(*ps_HLT_noalg_mb_L1MBTS_1_1);
            m_trig.push_back(*HLT_noalg_mb_L1MBTS_2);                               m_trig_ps.push_back(*ps_HLT_noalg_mb_L1MBTS_2);
            m_trig.push_back(*HLT_noalg_L1MBTS_A);                                  m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_A);
            m_trig.push_back(*HLT_noalg_L1MBTS_C);                                  m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_C);
            //m_trig.push_back(*HLT_mb_sptrk_pt2_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt2_L1MBTS_2);
            //m_trig.push_back(*HLT_mb_sptrk_pt4_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt4_L1MBTS_2);
            //m_trig.push_back(*HLT_mb_sptrk_pt6_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt6_L1MBTS_2);
            //m_trig.push_back(*HLT_mb_sptrk_pt8_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt8_L1MBTS_2);
            m_trig.push_back(*HLT_mb_mbts_L1MBTS_1);                                m_trig_ps.push_back(*ps_HLT_mb_mbts_L1MBTS_1);
            m_trig.push_back(*HLT_mb_mbts_L1MBTS_1_1);                              m_trig_ps.push_back(*ps_HLT_mb_mbts_L1MBTS_1_1);
            m_trig.push_back(*HLT_mb_mbts_L1MBTS_2);                                m_trig_ps.push_back(*ps_HLT_mb_mbts_L1MBTS_2);
            m_trig.push_back(*HLT_mb_mbts_L1RD0_FILLED);                            m_trig_ps.push_back(*ps_HLT_mb_mbts_L1RD0_FILLED);
            m_trig.push_back(*HLT_mb_sptrk_vetombts2in_L1RD0_FILLED);               m_trig_ps.push_back(*ps_HLT_mb_sptrk_vetombts2in_L1RD0_FILLED);
            m_trig.push_back(*HLT_mb_mbts_all_L1MBTS_A);                            m_trig_ps.push_back(*ps_HLT_mb_mbts_all_L1MBTS_A);
            m_trig.push_back(*HLT_mb_mbts_all_L1MBTS_C);                            m_trig_ps.push_back(*ps_HLT_mb_mbts_all_L1MBTS_C);
        }
        else{
            m_trig.push_back(*HLT_L1_ZDC_XOR_E1_E3);                                m_trig_ps.push_back(*ps_HLT_L1_ZDC_XOR_E1_E3);
            m_trig.push_back(*HLT_L1_ZDC_XOR_E2);                                   m_trig_ps.push_back(*ps_HLT_L1_ZDC_XOR_E2);
            m_trig.push_back(*HLT_noalg_L1ZDC_A_AND_C);                             m_trig_ps.push_back(*ps_HLT_noalg_L1ZDC_A_AND_C);
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_XOR_E2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_L1ZDC_XOR_E2);
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_XOR_E1_E3);                        m_trig_ps.push_back(*ps_HLT_mb_sptrk_L1ZDC_XOR_E1_E3);
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_A_AND_C);                          m_trig_ps.push_back(*ps_HLT_mb_sptrk_L1ZDC_A_AND_C);
            m_trig.push_back(*HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2);                 m_trig_ps.push_back(*ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2);
            m_trig.push_back(*HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3);              m_trig_ps.push_back(*ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3);
            m_trig.push_back(*HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C);                m_trig_ps.push_back(*ps_HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C);
        }

        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        trig_index = triggerIndex(m_trig); // retrive the index of the relevant trigger
        prescale = m_trig_ps.at(trig_index);
        lumiBlock = *lumiblock;
        ntrk = trk_pt->size();

        int i = myreader.GetCurrentEntry();  
        if(i>iend) {cout<<"Event loop finished! Analyzed "<< iend<<" events"<<endl;break;}
        if(i%100000==0) cout<<"proccesed "<<i<<" / "<<iend<<" events "<<" "<<fChain->GetFile()->GetName()<<endl;
        
        if(!((lumiBlock>537 && lumiBlock <859) || (lumiBlock>1014 && lumiBlock <4755))) continue; //stable beams
        if((*nvtx) != 2) continue; // right now i make sure only one primary vertex

        #ifdef check
        cout<<"pt="<<trk_pt->size()<<"  eta="<<trk_eta->size()<<"  phi="<<trk_phi->size()<<"  qop="<<trk_charge->size()
            <<"  quality="<<trk_quality->size()<<endl;
        #endif
        float sum =0;
        float sumA =0;
        float sumC =0;
        //sum the energy in the zdc
        for(int i = 0; i<4; i++){
            if(isBitSet(*BitMask,i)) sumC+= zdcWei.at(i)*ModAmp[i];
        }
        for(int i = 4; i<8; i++){
            if(isBitSet(*BitMask,i)) sumA+= zdcWei.at(i)*ModAmp[i];
        }
        //if(sumC >0.01) continue;
        sum = sumA + sumC;
        hzdc->Fill(sum,prescale);
        hZdcCorr->Fill(sumA,sumC,prescale);

        //determine effective energy bin
        float m_eff_energy = sqrt_s - sum;
        m_eff_energy = m_eff_energy/1000; //convert to TeV
        if(m_eff_energy < 0.0 || m_eff_energy > 13.6001) continue;
        m_cent_i          =Bins::GetCentBin(m_eff_energy);
        if(m_cent_i <0 || m_cent_i >=Bins::NCENT) continue;

        hzdc_after_cut->Fill(sum,prescale);
        //get multiplicity bin
        nbin = Bins::GetTrkBin(ntrk);
        if(nbin <0 || nbin >=Bins::NTRK) continue;
        //Z-vtx cuts
        m_zvtx=vtx_z->at(0);
        int zbin   = get_zPool(m_zvtx);
        if(zbin<0) continue;
        /*-----------------------------------------------------------------------------
         *  Fill some Monitor histograms
         *-----------------------------------------------------------------------------*/
        h_Zvtx->Fill(m_zvtx);
        hNtrk->Fill(ntrk);
        hNtrk_no_cut->Fill((*ntrk_ptr));
        //heff->Fill(m_cent_i);
        hNtrkEff->Fill(m_eff_energy,ntrk);


        /*-----------------------------------------------------------------------------
         *  Information about tracks in the event
         *-----------------------------------------------------------------------------*/
        Event *ev=new Event(i, m_cent_i, m_zvtx);
        for(int j=0; j<ntrk; j++){
            float pt=0,eta=10,phi=10;
            int   charge=-1;
            //if( (trk_qual->at(j)&m_track_quality)!=m_track_quality) continue; //TODO understand this
            pt =trk_pt->at(j)/1000.0;//Convert to GeV
            phi=trk_phi->at(j);
            eta=trk_eta->at(j);
            if(trk_charge->at(j)>0.0) charge=1;

            if(pt< 0.5) continue;
            if(pt>20.0) continue;//upper cutoff is to remove very high pT Tracks
            if (fabs(eta) > 2.4999) continue;

            float trk_eff=1;//eventually this stores 1/eff. TODO get this value

            int ptbin1 = Bins::GetPtBin1(pt);
            int ptbin2 = Bins::GetPtBin2(pt);
            if(ptbin1==-1 && ptbin2==-1) continue;
            N_trigger[m_cent_i][nbin]->Fill(pt,trk_eff);
            h_pt[m_cent_i]->Fill(pt);
            h_eta[m_cent_i]->Fill(eta);
            //N_ntrk[m_cent_i]->Fill(ntrk,trk_eff);
            if(ptbin1>=0) h_EtaPhi[ptbin1]->Fill(eta,phi);

            //push track into event
            ev->AddTrack(pt,eta,phi,charge,trk_eff,ptbin1,ptbin2);
        }
        Fill(ev,ev,0);
        FillMixed(ev);
     }
     //Normalize the histograms per event
    for(int icent=0; icent<Bins::NCENT; icent++){
        h_pt[icent]->Scale(1.0/iend);
        h_eta[icent]->Scale(1.0/iend);
    }
    SaveHistos();
    
}

void InitHistos(){
    //create output file
    
    std::string directory;
    directory = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    if(minbias) directory = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias";
    gSystem->Exec(Form("mkdir -p %s",directory.c_str()));
    tmpf = new TFile(Form("%s/%s",directory.c_str(),output_name),"recreate");

    
    /*-----------------------------------------------------------------------------
     *  Global monitor histograms
     *-----------------------------------------------------------------------------*/
    h_Zvtx = new TH1D("hzvtx", "hzvtx" , 300 , -300 , 300); h_Zvtx ->Sumw2();
    hzdc = new TH1D("hzdc", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_after_cut = new TH1D("hzdc_cut", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hNtrk = new TH1D("hNtrk", ";N_{ch}^{rec}; Events" , 200 , 0 , 200);hNtrk->Sumw2();
    hNtrk_no_cut = new TH1D("hNtrkNoCut", ";N_{ch}; Events" , 200 , 0 , 200);hNtrk_no_cut->Sumw2();
    hNtrkEff  = new TH2D("hNtrkEff", ";Effective energy [TeV];N_{ch}^{rec}" , Bins::NCENT, 0 , 13.6, Bins::NTRK,0,140);
    hZdcCorr  = new TH2D("hZdcCorr", ";side A;side C" , 200, 0 , 25000, 200,0,25000);
    
    
    //eta map of tracks
    for(int icent=0; icent<Bins::NCENT; icent++){
        sprintf(histname,"h_eta_icent%.2d",icent);
        h_eta[icent]=new TH1D(histname,";#eta;counts per event",50,-3.0,3.0);
        h_eta[icent]->Sumw2();
    }

    //pt map of tracks
    for(int icent=0; icent<Bins::NCENT; icent++){
        sprintf(histname,"h_pt_icent%.2d",icent);
        h_pt[icent]=new TH1D(histname,";p_{T} [GeV];counts per event",50,0,5);
        h_pt[icent]->Sumw2();
    }
    

    //Eta-phi map of tracks
    for(int ipt=0;ipt<Bins::NPT1;ipt++){
        sprintf(histname ,"h_EtaPhi_pt%d",ipt);
        sprintf(histtitle,"h_EtaPhi;#eta;#phi;N_{Tracks};Events");
        h_EtaPhi[ipt]=new TH2D(histname,histtitle,50,-2.5,2.5,64,-acos(-1.0),acos(-1.0));
    }

    //main fg and bg distributions
    int nBinsPhi = 36;
    int nBinsEta = 50;
    for(int icent=0; icent<Bins::NCENT; icent++){
        for(int itrk =0; itrk<Bins::NTRK; itrk++){
            for(int ipt1=0; ipt1<Bins::NPT1; ipt1++){
                for(int ipt2=0; ipt2<Bins::NPT2; ipt2++){
                    for(int it=0; it<Bins::NCH; it++){
                    sprintf(histname,"fg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d",icent,itrk,ipt1,ipt2,it);
                    fg[icent][itrk][ipt1][ipt2][it] = new TH2D(histname,histname,36,-PI/2,1.5*PI,50,0,5.0);
                    fg[icent][itrk][ipt1][ipt2][it]->Sumw2();

                    sprintf(histname,"bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d",icent,itrk,ipt1,ipt2,it);
                    bg[icent][itrk][ipt1][ipt2][it] = new TH2D(histname,histname,36,-PI/2,1.5*PI,50,0,5.0);
                    bg[icent][itrk][ipt1][ipt2][it]->Sumw2();
                    }
                }
            }
        }
    }

    //Trigger particle counter (for PTY)
    for(int icent=0; icent<Bins::NCENT; icent++){
        for(int itrk=0; itrk<Bins::NTRK; itrk++){
            sprintf(histname,"N_trigger_cent%.2d_trk%.2d",icent,itrk);
            N_trigger[icent][itrk] = new TH1D(histname,"N_trigger;pT_trigger;",200,0,20);
            N_trigger[icent][itrk]->Sumw2();
        }
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



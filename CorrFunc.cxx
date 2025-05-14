#include "Event.h"
#include "CorrFunc.h"
#include "bins.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/step_func.h"
//#define check


std::string weights_path;
int Trig;
bool same_side = true; // use same or opposite side calibration method
weights_path = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles";
if(same_side) weights_path = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide"; 
void CorrFunc(const int a ,const char* fileList, int Trig1 = 0){
    Trig = Trig1;
    sprintf(output_name,"histograms_%d.root",a);

    if(Trig ==0){
        cout << "working on ZDC trigger!" << endl;
    }
    else if(Trig ==1){
        cout << "Working on minbias trigger!" << endl;
    }
    else if(Trig == 2){
        cout << "Working on xor trigger!" << endl;
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
    char base[200] ="/gpfs0/citron/users/bargl/ZDC/lhcf22/user.bglik.data22_13p6TeV.00435229.ANALYSIS.2.noLhcf_EXT0.585561170.585561170";
    //char base[200] ="/gpfs0/citron/users/bargl/ZDC/lhcf22/user.bglik.data22_13p6TeV.00435229.ANALYSIS3.noLhcf_EXT0.596048724.596048724";
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
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_XOR_E1_E3                    (myreader, "HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_XOR_E1_E3                 (myreader, "HLT_mb_sptrk_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_XOR_E2                       (myreader, "HLT_noalg_L1ZDC_XOR_E2");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_XOR_E2                    (myreader, "HLT_mb_sptrk_L1ZDC_XOR_E2");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_A_AND_C                      (myreader, "HLT_noalg_L1ZDC_A_AND_C");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_A_AND_C                   (myreader, "HLT_mb_sptrk_L1ZDC_A_AND_C");
	// TTreeReaderValue<bool> HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2          (myreader,"HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2");
	// TTreeReaderValue<bool> HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3       (myreader, "HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3");
	//TTreeReaderValue<bool> HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C         (myreader, "HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C");
    TTreeReaderValue<float> ps_HLT_noalg_L1ZDC_XOR_E1_E3                     (myreader, "ps_HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<float> ps_HLT_noalg_L1ZDC_XOR_E2    (myreader, "ps_HLT_noalg_L1ZDC_XOR_E2");
    TTreeReaderValue<float> ps_HLT_noalg_L1ZDC_A_AND_C                  (myreader, "ps_HLT_noalg_L1ZDC_A_AND_C");
    TTreeReaderValue<float> ps_HLT_mb_sptrk_L1ZDC_XOR_E2                (myreader,"ps_HLT_mb_sptrk_L1ZDC_XOR_E2");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_L1ZDC_XOR_E1_E3             (myreader,"ps_HLT_mb_sptrk_L1ZDC_XOR_E1_E3");
	TTreeReaderValue<float> ps_HLT_mb_sptrk_L1ZDC_A_AND_C               (myreader, "ps_HLT_mb_sptrk_L1ZDC_A_AND_C");
	// TTreeReaderValue<float> ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2      (myreader,"ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E2");
	// TTreeReaderValue<float> ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3   (myreader, "ps_HLT_mb_sp100_trk30_hmt_L1ZDC_XOR_E1_E3");
	//TTreeReaderValue<float> ps_HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C     (myreader, "ps_HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C");


    //minbias triggers
    TTreeReaderValue<bool> HLT_noalg_L1MBTS_1(myreader, "HLT_noalg_L1MBTS_1");
    //TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_1(myreader, "HLT_noalg_mb_L1MBTS_1");
	// TTreeReaderValue<bool> HLT_noalg_L1MBTS_1_1(myreader, "HLT_noalg_L1MBTS_1_1");
	// TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_1_1(myreader, "HLT_noalg_mb_L1MBTS_1_1");
	TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_2(myreader, "HLT_noalg_mb_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_2(myreader, "HLT_mb_mbts_L1MBTS_2");
	// TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_1(myreader, "HLT_mb_mbts_L1MBTS_1");
	// TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_1_1(myreader, "HLT_mb_mbts_L1MBTS_1_1");
	
    TTreeReaderValue<float> ps_HLT_noalg_L1MBTS_1(myreader, "ps_HLT_noalg_L1MBTS_1");
    //TTreeReaderValue<float> ps_HLT_noalg_mb_L1MBTS_1(myreader, "ps_HLT_noalg_mb_L1MBTS_1");
	//TTreeReaderValue<float> ps_HLT_noalg_L1MBTS_1_1(myreader, "ps_HLT_noalg_L1MBTS_1_1");
	//TTreeReaderValue<float> ps_HLT_noalg_mb_L1MBTS_1_1(myreader, "ps_HLT_noalg_mb_L1MBTS_1_1");
	TTreeReaderValue<float> ps_HLT_noalg_mb_L1MBTS_2(myreader, "ps_HLT_noalg_mb_L1MBTS_2");
	TTreeReaderValue<float> ps_HLT_mb_mbts_L1MBTS_2(myreader, "ps_HLT_mb_mbts_L1MBTS_2");
    
    //zdc
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp"); // First 4 lower bits corressponds to side C
    TTreeReaderArray<float> Amp(myreader, "zdc_ZdcAmp"); 
    TTreeReaderArray<short> Status(myreader, "zdc_ZdcStatus");
    TTreeReaderArray<short> LucrodTriggerSideAmp(myreader, "zdc_ZdcLucrodTriggerSideAmp"); 
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");
    TTreeReaderArray<unsigned int> ModStatus(myreader, "zdc_ZdcModuleStatus");
    
    InitHistos();

    //----------------------------------coorections for LHCf -------------------------------------------
    //reading moments for late shower cut. DONT FORGET TO UNCOMMENT THE RELEVANT PART IN THE EVENT LOOP!
    std::vector<std::vector<float>> moments_C;
    std::vector<std::vector<float>> moments_A;
    for(int ilb =0; ilb<Bins::NLB; ilb++){
        moments_C.push_back(GetMoments(ilb,0, Trig));
        moments_A.push_back(GetMoments(ilb,1, Trig));
    }

    //getting scale factors for stability correction
    float scale_stability[2][4][Bins::NLB];
    for(int side =0; side<2; side++){
        std::string path_par;
        path_par = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp";
        std::vector<int> lb_center = {Bins::GetLbIndex(1816,1887), Bins::GetLbIndex(1902,1928), Bins::GetLbIndex(2812,2996), Bins::GetLbIndex(2996,3109),Bins::GetLbIndex(3218,3678)};
        TFile *flb_mean_and = new TFile(Form("%s/mean_block1.root",path_par.c_str()),"read");
        for(int mod=0; mod<4; mod++){
            TH1D* h_mean_and = (TH1D*)flb_mean_and->Get(Form("h_mean_side%i_mod%i",side,mod));
            double mean_center_and =0;
            for(const auto& lumi : lb_center){
                mean_center_and += h_mean_and->GetBinContent(lumi + 1)/lb_center.size();
            }
            for(int ilb=0; ilb<Bins::NLB; ilb++){
                scale_stability[side][mod][ilb] = mean_center_and/h_mean_and->GetBinContent(ilb+1);
                //scale_stability[side][mod][ilb] = 1.0;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------

    int nentries= myreader.GetEntries();
    int istart=0,iend= nentries; //by default loop over whole dataset
    cout << "Analyzing : "<< iend << " events" << endl;
    
      while (myreader.Next()){
        //retrive trigger information
        std::vector<bool> m_trig;
        std::vector<float> m_trig_ps;
        if(Trig == 0){
            m_trig.push_back(*HLT_noalg_L1ZDC_A_AND_C);                             m_trig_ps.push_back(*ps_HLT_noalg_L1ZDC_A_AND_C);
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_A_AND_C);                             m_trig_ps.push_back(*ps_HLT_mb_sptrk_L1ZDC_A_AND_C);
            //m_trig.push_back(*HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C);                             m_trig_ps.push_back(*ps_HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C);
        }
        else if(Trig == 1){
            m_trig.push_back(*HLT_noalg_L1MBTS_1);                                  m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_1);
            //m_trig.push_back(*HLT_noalg_mb_L1MBTS_1);                               m_trig_ps.push_back(*ps_HLT_noalg_mb_L1MBTS_1);
            //m_trig.push_back(*HLT_noalg_L1MBTS_1_1);                                m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_1_1);
            //m_trig.push_back(*HLT_noalg_mb_L1MBTS_1_1);                             m_trig_ps.push_back(*ps_HLT_noalg_mb_L1MBTS_1_1);
            m_trig.push_back(*HLT_noalg_mb_L1MBTS_2);                               m_trig_ps.push_back(*ps_HLT_noalg_mb_L1MBTS_2);
            //m_trig.push_back(*HLT_noalg_L1MBTS_A);                                  m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_A);
            //m_trig.push_back(*HLT_noalg_L1MBTS_C);                                  m_trig_ps.push_back(*ps_HLT_noalg_L1MBTS_C);
            // m_trig.push_back(*HLT_mb_sptrk_pt2_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt2_L1MBTS_2);
            // m_trig.push_back(*HLT_mb_sptrk_pt4_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt4_L1MBTS_2);
            // m_trig.push_back(*HLT_mb_sptrk_pt6_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt6_L1MBTS_2);
            // m_trig.push_back(*HLT_mb_sptrk_pt8_L1MBTS_2);                           m_trig_ps.push_back(*ps_HLT_mb_sptrk_pt8_L1MBTS_2);
            //m_trig.push_back(*HLT_mb_mbts_L1MBTS_1);                                m_trig_ps.push_back(*ps_HLT_mb_mbts_L1MBTS_1);
            //m_trig.push_back(*HLT_mb_mbts_L1MBTS_1_1);                              m_trig_ps.push_back(*ps_HLT_mb_mbts_L1MBTS_1_1);
            m_trig.push_back(*HLT_mb_mbts_L1MBTS_2);                                m_trig_ps.push_back(*ps_HLT_mb_mbts_L1MBTS_2);
            //m_trig.push_back(*HLT_mb_mbts_L1RD0_FILLED);                            m_trig_ps.push_back(*ps_HLT_mb_mbts_L1RD0_FILLED);
            //m_trig.push_back(*HLT_mb_sptrk_vetombts2in_L1RD0_FILLED);               m_trig_ps.push_back(*ps_HLT_mb_sptrk_vetombts2in_L1RD0_FILLED);
            //m_trig.push_back(*HLT_mb_mbts_all_L1MBTS_A);                            m_trig_ps.push_back(*ps_HLT_mb_mbts_all_L1MBTS_A);
            //m_trig.push_back(*HLT_mb_mbts_all_L1MBTS_C);                            m_trig_ps.push_back(*ps_HLT_mb_mbts_all_L1MBTS_C);
        }
        else if(Trig ==2){
            m_trig.push_back(*HLT_noalg_L1ZDC_XOR_E2);                                m_trig_ps.push_back(*ps_HLT_noalg_L1ZDC_XOR_E2);                            
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_XOR_E2);                             m_trig_ps.push_back(*ps_HLT_mb_sptrk_L1ZDC_XOR_E2);                                                           
            m_trig.push_back(*HLT_noalg_L1ZDC_XOR_E1_E3);                             m_trig_ps.push_back(*ps_HLT_noalg_L1ZDC_XOR_E1_E3);                                
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_XOR_E1_E3);                          m_trig_ps.push_back(*ps_HLT_mb_sptrk_L1ZDC_XOR_E1_E3);  
        }
        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        trig_index = triggerIndex(m_trig); // retrive the index of the relevant trigger
        //prescale = m_trig_ps.at(trig_index);
        prescale = 1;
        //N_evt_sampled += prescale;
        lumiBlock = *lumiblock;
        ntrk = trk_pt->size();
        if(*nvtx ==1) continue;
        if(ntrk <= 1) continue;
        int i = myreader.GetCurrentEntry();
         
        if(i>iend) {cout<<"Event loop finished! Analyzed "<< iend<<" events"<<endl;break;}
        if(i%100000==0) cout<<"proccesed "<<i<<" / "<<iend<<" events "<<" "<<fChain->GetFile()->GetName()<<endl;
        
        if(!((lumiBlock>1014 && lumiBlock <4755))) continue; //stable beams
        if(lumiBlock < 1816) continue; //Timing adjust to Had1C for trigger
        if(lumiBlock >=1889 && lumiBlock <= 1901) continue;
        if(lumiBlock >=1930 && lumiBlock <= 1935) continue;
        if(lumiBlock >=2147 && lumiBlock <= 2253) continue;
        if(lumiBlock >=2807 && lumiBlock <= 2811) continue;
        if(lumiBlock >=3111 && lumiBlock <= 3217) continue;
        if(lumiBlock >=3680 && lumiBlock <= 3686) continue;
        if(lumiBlock >=3766 && lumiBlock <= 3866) continue;
        if(lumiBlock >=4563 && lumiBlock <= 4654) continue;
          
        #ifdef check
        cout<<"pt="<<trk_pt->size()<<"  eta="<<trk_eta->size()<<"  phi="<<trk_phi->size()<<"  qop="<<trk_charge->size()
            <<"  quality="<<trk_quality->size()<<endl;
        #endif
        if(isBitSet(Status[0],2) || isBitSet(Status[1],2)) continue; // reject failbits
        
        //--------------------------getting relevant moments for late shower cut-------------
        int ilb_bin = Bins::GetLbBin(lumiBlock);
        if(ilb_bin == -1) continue;
        float mean_C = moments_C[ilb_bin].at(0);
        float stdev_C = moments_C[ilb_bin].at(1);
        float mean_A = moments_A[ilb_bin].at(0);
        float stdev_A = moments_A[ilb_bin].at(1);
        if(isBitSet((*BitMask),2) && (ModAmp[2]/Amp[0] > (mean_C + 2*stdev_C))) continue;
        if(isBitSet((*BitMask),2+4) && (ModAmp[2+4]/Amp[1] > (mean_A + 2*stdev_A))) continue;
        //----------------------------------------------------------------------------------

        /*-----------------------------------------------------------------------------
        *  Fill ZDC techinicals
        *-----------------------------------------------------------------------------*/
        hLucrodCorr->Fill(LucrodTriggerSideAmp[0],LucrodTriggerSideAmp[1]);
        hLucrodZdcCorr[0]->Fill(LucrodTriggerSideAmp[0], Amp[0]);
        hLucrodZdcCorr[1]->Fill(LucrodTriggerSideAmp[1], Amp[1]);
        hAmpCorr->Fill(Amp[0], Amp[1]);
        //------------------------------------------------------------------------------
        float sum =0;
        float sumA =0;
        float sumC =0;
        //sum the energy in the zdc
        bool LGoverflow = false;
        for(int i = 1; i<4; i++){
            if(isBitSet(*BitMask,i)) sumC+= zdcWei.at(i)*ModAmp[i]*scale_stability[0][i][ilb_bin];
            if(isBitSet(ModStatus[i],6)) LGoverflow = true;
        }
        for(int i = 5; i<8; i++){
            if(isBitSet(*BitMask,i)) sumA+= zdcWei.at(i)*ModAmp[i]*scale_stability[1][i-4][ilb_bin];
            if(isBitSet(ModStatus[i],6)) LGoverflow = true;
        }
        if(LGoverflow) continue; //reject overflowed events
        // //take into account only hadronic process for zdc triggers
        // if(Trig == 0 || Trig == 2){
        //     if(!isElectroMagnetic(0,(*BitMask))){
        //         sum+=sumC;
        //     }
        //     if(!isElectroMagnetic(1,(*BitMask))){
        //         sum+=sumA;
        //     }
        // }
        // else{
        //     sum = sumA + sumC;
        // }
        sum = sumA + sumC;
        //reject events when both sides were purely EM for zdc triggers
        //if(!(Trig==1) && sum == 0) continue;
        //if on "AND" trigger reject events when both sides didnt catch the leading baryon
        //if(Trig == 0 && (isElectroMagnetic(0,(*BitMask)) || isElectroMagnetic(1,(*BitMask)))) continue;
        hzdc->Fill(sum,prescale);
        hZdcCorr->Fill(sumC,sumA);
        if((*nvtx) == 2){
            hzdc_A_without_pileup_noOfflineCut->Fill(sumA,prescale);
            hzdc_C_without_pileup_noOfflineCut->Fill(sumC,prescale);
        } 
        /*-----------------------------------------------------------------------------
        *  Offline cuts
        *-----------------------------------------------------------------------------*/
        if(sumA > 5000 || sumC > 5000) continue;
        bool pass_offline = true;
        if(Trig ==0){
            if(!(sumA >1000 && sumC >1000)) pass_offline = false;
        }
        if(Trig ==2){
            if(!((sumA > 1000 && sumC <1000) || (sumA <1000 && sumC >1000))) pass_offline = false;
        }
        if(!pass_offline) continue;
        //-----------------------------------------------------------------------------------------
        hzdc_A_with_pileup->Fill(sumA,prescale);
        hzdc_C_with_pileup->Fill(sumC,prescale);
        if((*nvtx) != 2) continue; //pile up rejection
        hzdc_A_without_pileup->Fill(sumA,prescale);
        hzdc_C_without_pileup->Fill(sumC,prescale);
        
        //determine effective energy bin
        float m_eff_energy = sqrt_s - sum;
        m_eff_energy = m_eff_energy/1000; //convert to TeV
        heff->Fill(m_eff_energy);
        // if(m_eff_energy < 0.0 || m_eff_energy > 13.6001) continue;
        m_cent_i          =Bins::GetCentBin(m_eff_energy);
        if(m_cent_i <0 || m_cent_i >=Bins::NCENT) continue;
        //get multiplicity bin
        nbin = Bins::GetTrkBin(ntrk);
        if(nbin <0 || nbin >=Bins::NTRK) continue;
        //Z-vtx cuts
        m_zvtx=vtx_z->at(0);
        int zbin   = get_zPool(m_zvtx);
        if(zbin<0) continue;
        good_events[m_cent_i] +=1;
        /*-----------------------------------------------------------------------------
        *  Fill some Monitor histograms
        *-----------------------------------------------------------------------------*/
        hzdc_after_cut->Fill(sum);
        heff_after_cut->Fill(m_eff_energy);
        h_Zvtx->Fill(m_zvtx);
        hNtrk[m_cent_i]->Fill(ntrk);
        hNtrk_no_cut->Fill((*ntrk_ptr));
        //h_eff_no_ps->Fill(m_cent_i);
        //h_Trig->Fill(trig_index);
        hNtrkEff->Fill(m_eff_energy,ntrk);
        hZdcCorr_after_cut->Fill(sumC,sumA);
        


        /*-----------------------------------------------------------------------------
         *  Information about tracks in the event
         *-----------------------------------------------------------------------------*/
         int trk_with_pt = 0; //for testing
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
            h_pt[m_cent_i]->Fill(pt);
            h_eta[m_cent_i]->Fill(eta);
            N_trigger[m_cent_i][nbin]->Fill(pt,trk_eff);
            //N_ntrk[m_cent_i]->Fill(ntrk,trk_eff);
            if(ptbin1>=0) h_EtaPhi[ptbin1]->Fill(eta,phi);
            if (ptbin1 == 2 && m_cent_i == 9) trk_with_pt +=1; // for testing

            //push track into event
            ev->AddTrack(pt,eta,phi,charge,trk_eff,ptbin1,ptbin2);
        }
        h_trk_test->Fill(trk_with_pt); // for testing
        Fill(ev,ev,0);
        FillMixed(ev);
     }

    SaveHistos();
    
}

void InitHistos(){
    //create output file
    
    std::string directory;
    directory = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles");
    if (same_side) directory = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide");
    if(Trig == 1) directory = Form("%s/minbias",directory.c_str());
    if(Trig == 2) directory = Form("%s/xor",directory.c_str());
    gSystem->Exec(Form("mkdir -p %s",directory.c_str()));
    tmpf = new TFile(Form("%s/%s",directory.c_str(),output_name),"recreate");
    //tree = new TTree("tree", "A simple tree");
    //tree->Branch("N_evt_sampled", &N_evt_sampled, "N_evt_sampled/I"); // Create a branch with the variable
    
    //ZDC technical histograms (ADC)
    hLucrodCorr = new TH2D("hLucrodCorr",";zdc_ZdcLucrodTriggerSideAmp(0) [a.u];zdc_ZdcLucrodTriggerSideAmp(1) [a.u]", 100,0,5000, 100,0,5000);
    hAmpCorr = new TH2D("hAmpCorr",";zdc_ZdcAmp(0) [ADC];zdc_ZdcAmp(1) [ADC]", 80,0,8000, 80,0,8000);
    hLucrodZdcCorr[0] = new TH2D("hLucrodZdcCorr0",";zdc_ZdcLucrodTriggerSideAmp(0) [a.u];zdc_ZdcAmp(0) [ADC]", 100,0,5000, 80,0,8000);
    hLucrodZdcCorr[1] = new TH2D("hLucrodZdcCorr1",";zdc_ZdcLucrodTriggerSideAmp(1) [a.u];zdc_ZdcAmp(1) [ADC]", 100,0,5000, 80,0,8000);
    

    /*-----------------------------------------------------------------------------
     *  Global monitor histograms
     *-----------------------------------------------------------------------------*/
    h_Zvtx = new TH1D("hzvtx", "hzvtx" , 300 , -300 , 300); h_Zvtx ->Sumw2();
    hzdc = new TH1D("hzdc", ";ZDC energy [GeV]; Counts" , 200 , 0 , 50000);
    hzdc_A_without_pileup_noOfflineCut = new TH1D("hzdc_A_without_pileup_noOfflineCut", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_C_without_pileup_noOfflineCut = new TH1D("hzdc_C_without_pileup_noOfflineCut", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_A_with_pileup = new TH1D("hzdc_A_with_pileup", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_C_with_pileup = new TH1D("hzdc_C_with_pileup", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_A_without_pileup = new TH1D("hzdc_A_without_pileup", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_C_without_pileup = new TH1D("hzdc_C_without_pileup", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hzdc_after_cut = new TH1D("hzdc_cut", ";ZDC energy [GeV]; Counts" , 200 , 0 , 25000);
    hNtrk_no_cut = new TH1D("hNtrkNoCut", ";N_{ch}; Events" , 200 , 0 , 200);hNtrk_no_cut->Sumw2();
    hNtrkEff  = new TH2D("hNtrkEff", ";Effective energy [TeV];N_{ch}^{rec}" , Bins::NCENT, Bins::CENT_LO[0] , Bins::CENT_HI[Bins::NCENT-1], 15,0,150); hNtrkEff->Sumw2();
    hZdcCorr  = new TH2D("hZdcCorr", ";side C;side A" , 200, 0 , 25000, 200,0,25000);
    hZdcCorr_after_cut  = new TH2D("hZdcCorr_after_cut", ";side C;side A" , 200, 0 , 25000, 200,0,25000);
    heff_after_cut  = new TH1D("heff_after_cut", ";Eff energy [TeV];" , Bins::NCENT, Bins::CENT_LO[0] , Bins::CENT_HI[Bins::NCENT-1]);
    heff  = new TH1D("heff", ";Eff energy [TeV];" ,Bins::NCENT, 0 , 13.61);
    
    
    //eta map of tracks
    for(int icent=0; icent<Bins::NCENT; icent++){
        sprintf(histname,"h_eta_icent%.2d",icent);
        h_eta[icent]=new TH1D(histname,";#eta;counts per event",20,-2.5,2.5);
        h_eta[icent]->Sumw2();
    }

    //pt map of tracks
    for(int icent=0; icent<Bins::NCENT; icent++){
        sprintf(histname,"h_pt_icent%.2d",icent);
        h_pt[icent]=new TH1D(histname,";p_{T} [GeV];counts per event",39,0.5,20);
        h_pt[icent]->Sumw2();
    }
    
    //ntrks
    for(int icent=0; icent<Bins::NCENT; icent++){
        sprintf(histname,"h_ntrk_icent%.2d",icent);
        hNtrk[icent] = new TH1D(histname, ";N_{ch}^{rec}; Events" ,Bins::NTRK,0,Bins::TRK_HI[Bins::NTRK -1]);
        hNtrk[icent]->Sumw2();
    }

    //Eta-phi map of tracks
    for(int ipt=0;ipt<Bins::NPT1;ipt++){
        sprintf(histname ,"h_EtaPhi_pt%d",ipt);
        sprintf(histtitle,"h_EtaPhi;#eta;#phi;N_{Tracks};Events");
        h_EtaPhi[ipt]=new TH2D(histname,histtitle,50,-2.5,2.5,36,-acos(-1.0),acos(-1.0));
        h_EtaPhi[ipt]->Sumw2();
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
                    fg[icent][itrk][ipt1][ipt2][it] = new TH2D(histname,histname,nBinsPhi,-PI/2,1.5*PI,nBinsEta,0,5.0);
                    fg[icent][itrk][ipt1][ipt2][it]->Sumw2();

                    sprintf(histname,"bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d",icent,itrk,ipt1,ipt2,it);
                    bg[icent][itrk][ipt1][ipt2][it] = new TH2D(histname,histname,nBinsPhi,-PI/2,1.5*PI,nBinsEta,0,5.0);
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

    //for testing
    h_trk_test = new TH1D("h_trk_test", ";N_{ch}; Events" , 200 , 0 , 200);

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

void load_weights(){
    std::cout << "Initialize ZDC Weights!" << std::endl;
    std::string zdc_C; 
    std::string zdc_A; 
    zdc_C = Form("%s/zdcWeights_side0_itr10.root",weights_path.c_str());
    zdc_A = Form("%s/zdcWeights_side1_itr10.root",weights_path.c_str());
    std::cout << zdc_C << std::endl;
    std::cout << zdc_A << std::endl;
    //converting TvectorD to std::vector object
    TFile *fileC = new TFile(Form("%s",zdc_C.c_str()));
    TFile *fileA = new TFile(Form("%s",zdc_A.c_str()));
    TVectorD *tvec_C = nullptr;
    TVectorD *tvec_A = nullptr;
    fileC->GetObject("gains_avg", tvec_C);
    fileA->GetObject("gains_avg", tvec_A);
    if (!tvec_C || !tvec_A) {
        std::cerr << "Error: TVectorD object not found in file!" << std::endl;
        exit(-1);
    }
    for (int i = 0; i < 8; ++i) {
        if(i<4) zdcWei.push_back((*tvec_C)[i]);
        //if(i<4) zdcWei.push_back(1.0);
        else{
                zdcWei.push_back((*tvec_A)[i-4]);
                //zdcWei.push_back(1.0);
        }
        if(i ==0 || i==4) zdcWei.at(i) = 0.0;
        cout << zdcWei.at(i) << endl;   
    }
}

bool isBitSet(int x, int s){
  int mask = x >> s;
  return mask % 2;
}

void SaveHistos(){
  tmpf->cd();
  tmpf->Write();
  cout << "Saving Finished" << endl;
}

int get_zPool(float z){
  int bin=-1;
  if(fabs(z)>=ZMAX) return -1;
  bin =(int)  (((z+ZMAX)/(2.0*ZMAX))*nz);  // pools
  if(bin<0||bin>=nz) return -1;
  return bin;
}

void Fill(Event* event1, Event* event2, int mixtype) {

    int icent = (m_cent_i);
    if (icent < 0 || icent >= Bins::NCENT) return;

    int itrk = (nbin);
    //int itrk = Bins::GetTrkBin(event1->get_npart());
      if (itrk < 0 || itrk >= Bins::NTRK) return;
    
    int ntrk1 = event1->get_npart();
    int ntrk2 = event2->get_npart();

    float phi01, eta1, qop1, eff1; //,pt1;
    float phi02, eta2, qop2, eff2; //,pt2;
    int   ptbin1, ptbin2;

    int sig1 = (int)(ran->Rndm() * 2);
    float wei = 1.0 / nmix;

    for (int i = 0; i < ntrk1; i++) {
        Track* trk1 = event1->GetTrack(i);
        eta1   = (trk1)->get_eta();
        phi01  = (trk1)->get_phi0();
        qop1   = (trk1)->get_charge();
        ptbin1 = (int)((trk1)->get_ptbin1());//trigger bin
        eff1   = (trk1)->get_eff();
        if (ptbin1 < 0) continue;

        int j = 0;
        while (j < ntrk2)  {
            if (mixtype == TWOPCTYPE::FG_HADRON_HADRON && j == i) {j++; continue;}
            Track* trk2 = event2->GetTrack(j);
            j++;
            eta2   = (trk2)->get_eta();
            phi02  = (trk2)->get_phi0();
            qop2   = (trk2)->get_charge();
            ptbin2 = (int) ((trk2)->get_ptbin2());
            eff2   = (trk2)->get_eff();
            if (ptbin2 < 0) continue;

            float dphi0 = phi01 - phi02;
            float deta  = fabs(eta1 - eta2);
            //if (sig1) dphi0 = -dphi0;
            if      (dphi0 > 1.5 * Common::PI)  dphi0 -= 2 * Common::PI;
            else if (dphi0 < -0.5 * Common::PI) dphi0 += 2 * Common::PI;

            if (mixtype == TWOPCTYPE::FG_HADRON_HADRON) {
                if ((qop1 * qop2) > 0) fg[icent][itrk][ptbin1][ptbin2][0]->Fill(dphi0, deta, eff1 * eff2);
                else                   fg[icent][itrk][ptbin1][ptbin2][1]->Fill(dphi0, deta, eff1 * eff2);
            }
            else if (mixtype == TWOPCTYPE::BG_HADRON_HADRON) {
                if ((qop1 * qop2) > 0) bg[icent][itrk][ptbin1][ptbin2][0]->Fill(dphi0, deta, wei * eff1 * eff2);
                else                   bg[icent][itrk][ptbin1][ptbin2][1]->Fill(dphi0, deta, wei * eff1 * eff2);
            }
        }
        //cout << j <<" " <<ev->Get << endl;
    }
}

bool FillMixed(EVENT_PTR event) {
    int zbin    = get_zPool   (m_zvtx);
    int indx    = zbin + nbin*nz + m_cent_i * nz*n_multiplicity; //a * n_B * n_C + b * n_C + c
    //int dep     = depth[m_cent_i];

    if (zbin == -1) {
        std::cout << "AAAAAAAA " << zbin << "  " << m_zvtx << std::endl;
        delete event;
        return false;
    }

    nmix = pool[indx].size();
    if (nmix > 0 && nmix <= dep) {
        for (int k = 0; k < nmix; k++) {
            Fill(pool[indx][k], event    , TWOPCTYPE::BG_HADRON_HADRON); //;mixing
        }
    }


    if (nmix >= dep) {
        delete pool[indx].at(0);
        pool[indx].erase( (pool[indx]).begin() );
    }
    pool[indx].push_back( event );

    return true;
}

bool isElectroMagnetic(int side, int mask){
   if (side != 0 && side != 1) {
        std::cerr << "Error: isElectroMagnetic() - Invalid side value provided. Must be 0 or 1." << std::endl;
        throw std::invalid_argument("Invalid side value provided");
    }

    // Adjust the bit positions based on the side
    int offset = (side == 1) ? 4 : 0;

    // Check the hadronic bits
    bool had1 = isBitSet(mask, 1 + offset);
    bool had2 = isBitSet(mask, 2 + offset);
    bool had3 = isBitSet(mask, 3 + offset);

    // Return true if had1 is set, but had2 and had3 are not
    return had1 && !had2 && !had3;
}

std::vector<float> GetMoments(int ilb_bin, int side, int Trig){
    TFile *input_histo; 
    std::vector<string> lb_range;
    for(int ilb=0; ilb<Bins::NLB; ilb++){
        lb_range.push_back(Bins::label_lb(ilb));
    }
    if(Trig == 0 ) input_histo = new TFile(Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/histograms.root"),"read");
    if(Trig == 1 ) input_histo = new TFile(Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/minbias/histograms.root"),"read");
    if(Trig == 2 ) input_histo = new TFile(Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp/xor/histograms.root"),"read");
    TH2D *h2 = (TH2D*)input_histo->Get(Form("h%i_h2AmpRat_Amp_ilb_%i",side,ilb_bin));
    //h2->Print();
    std::vector<float> mom;
    mom.push_back(h2->GetMean(2));
    mom.push_back(h2->GetStdDev(2));
    return mom;
}
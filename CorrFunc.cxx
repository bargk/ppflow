#include "Event.h"
#include "Defs.h"
#include "bins.h"
//#define check




void CorrFunc(const int a ,const char* fileList){
    sprintf(output_name,"histograms_%d.root",a);

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
    char base[200] ="/gpfs0/citron/users/bargl/ZDC/lhcf22/user.steinber.data22_13p6TeV.00435229.physics_MinBias.merge.AOD.r14470_p5587.4zdc_EXT0";
    // Open the files
    for (const auto file : files) {
        fChain->Add(Form("%s/%s",base,file.c_str()),0);
        if (!fChain || fChain->IsZombie()) {
            std::cerr << "Error opening file: " << file << std::endl;
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
    TTreeReaderValue<bool> L1_ZDC_XOR_E1_E3(myreader, "HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_A_AND_C(myreader, "HLT_noalg_L1ZDC_A_AND_C");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_OR(myreader, "HLT_noalg_L1ZDC_OR");
    TTreeReaderValue<bool> trk30(myreader, "HLT_mb_sp100_trk30_hmt_L1ZDC_A_AND_C");
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp"); // First 4 lower bits corressponds to side C
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");

    InitHistos();

    int nentries= myreader.GetEntries();
    int istart=0,iend= nentries; //by default loop over wholedataset
    cout << "Total events : "<< nentries << endl;
    
      while (myreader.Next()){
        
        lumiBlock = *lumiblock;
        ntrk = *ntrk_ptr;
        int i = myreader.GetCurrentEntry();
        
        if(i>iend) {cout<<"Event loop finished! Analyzed "<< iend<<" events"<<endl;break;}
        if(i%100000==0) cout<<"proccesed "<<i<<" / "<<iend<<" events "<<" "<<fChain->GetFile()->GetName()<<endl;
        

        //if(!(*HLT_noalg_L1ZDC_A_AND_C)) continue;
        //if(!(*HLT_noalg_L1ZDC_OR)) continue;
        if(!((lumiBlock>537 && lumiBlock <859) || (lumiBlock>1014 && lumiBlock <4755))) continue; //stable beams
        if((*nvtx) != 2) continue; // right now i make sure only one primary vertex

        #ifdef check
        cout<<"pt="<<trk_pt->size()<<"  eta="<<trk_eta->size()<<"  phi="<<trk_phi->size()<<"  qop="<<trk_charge->size()
            <<"  quality="<<trk_quality->size()<<endl;
        #endif
        int sum =0;
        //sum the energy in the zdc
        for(int i = 0; i<8; i++){
                if(isBitSet(*BitMask,i)){ sum+= zdcWei.at(i)*ModAmp[i];}
        }

        //determine effective energy bin
        float m_eff_energy_percent = GetEffectiveEnergy(sum);
        if(m_eff_energy_percent < 0 || m_eff_energy_percent > 100) continue;
        m_cent_i          =Bins::GetCentBin(m_eff_energy_percent);
        //m_cent_i          =Bins::GetCentBin(float(ntrk));
        if(m_cent_i <0 || m_cent_i >=Bins::NCENT) continue;


        

        //Z-vtx cuts
        m_zvtx=vtx_z->at(0);
        int zbin   = get_zPool(m_zvtx);
        if(zbin<0) continue;

        /*-----------------------------------------------------------------------------
         *  Fill some Monitor histograms
         *-----------------------------------------------------------------------------*/
        h_Zvtx->Fill(m_zvtx);
        //hNtrk->Fill(ntrk);
        //heff->Fill(m_cent_i);
        hNtrkEff->Fill(m_eff_energy_percent,ntrk/0.8);

        //TODO write number of good tracks based on bitword in trk_quality

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

            float trk_eff=1/0.8;//eventually this stores 1/eff. TODO get this value

            
            int ptbin1 = Bins::GetPtBin1(pt);
            int ptbin2 = Bins::GetPtBin2(pt);
            if(ptbin1==-1 && ptbin2==-1) continue;
            N_trigger[m_cent_i]->Fill(pt,trk_eff);
            h_eta[m_cent_i]->Fill(eta,trk_eff);
            //N_ntrk[m_cent_i]->Fill(ntrk,trk_eff);
            if(ptbin1>=0) h_EtaPhi[ptbin1]->Fill(eta,phi);

            //push track into event
            ev->AddTrack(pt,eta,phi,charge,trk_eff,ptbin1,ptbin2);
        }

        Fill(ev,ev,0);
        FillMixed(ev);
     }
    SaveHistos();
    
}




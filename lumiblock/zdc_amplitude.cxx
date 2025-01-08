//Ploting the amplitudes of the modules.
//The code has same structrue as CorrFunc.cxx but without the 2pc anaysis.

#include "bins.h"
std::vector<bool> m_trig;

bool isBitSet(int x, int s){
  int mask = x >> s;
  return mask % 2;
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


void zdc_amplitude(const int a ,const char* fileList, int Trig1 = 0){
    //load the bad lumiblocks
    TFile *finput = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/badLB.root",read);
    int Trig = Trig1;
    char output_name[100];
    sprintf(output_name,"histograms_%d.root",a);

    if(Trig ==0){
        cout << "working on ZDC trigger!" << endl;
    }
    else if(Trig ==1){
        cout << "Working on minbias trigger!" << endl;
    }
    else if(Trig == 2){
        cout << "Working on xorE2 trigger!" << endl;
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
    TChain *fChain = new TChain("zdcTree");
    char base[200] ="/eos/project-p/ppflow2022/user.bglik.data22_13p6TeV.00435229.ANALAYSIS_1_noLhcf_EXT0";
    // Open the files
    for (const auto file : files) {
        fChain->Add(Form("%s/%s",base,file.c_str()),0);
        if (!fChain || fChain->IsZombie()) {
            std::cerr << "Error opening file: " << file <<" Check files.txt content"<< std::endl;
            //continue;
        }
        std::cout << "Opened file: " << file << std::endl;
    }


    std::string directory;
    directory = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/Amp");
    if(Trig == 1) directory = Form("%s/minbias",directory.c_str());
    if(Trig == 2) directory = Form("%s/xorE2",directory.c_str());
    gSystem->Exec(Form("mkdir -p %s",directory.c_str()));
    TFile *tmpf = new TFile(Form("%s/%s",directory.c_str(),output_name),"recreate");


    TH2D *h_module[2][4];
    for(int side =0; side<2; side++){
        for(int mod =0; mod<4; mod++){
            h_module[side][mod]= new TH2D(Form("h%i_%i",side,mod),Form(";LB;HAD[%i][%i]",side,mod),3800,1000,4800,100,0,20000);
        }
    }
    
    //zdc
    TTreeReader myreader(fChain);
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp"); // First 4 lower bits corressponds to side C
    TTreeReaderArray<float> Amp(myreader, "zdc_ZdcAmp"); 
    TTreeReaderArray<short> LucrodTriggerSideAmp(myreader, "zdc_ZdcLucrodTriggerSideAmp"); 
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");
    TTreeReaderValue<unsigned int> lumiblock(myreader, "lumiBlock");
    //zdc triggers
    TTreeReaderValue<bool> HLT_L1_ZDC_XOR_E1_E3                         (myreader, "HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_L1_ZDC_XOR_E2                            (myreader, "HLT_noalg_L1ZDC_XOR_E2");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_A_AND_C                      (myreader, "HLT_noalg_L1ZDC_A_AND_C");

    while(myreader.Next()){
        int lumiBlock = *lumiblock;
        if(!((lumiBlock>1014 && lumiBlock <4755))) continue; //stable beams
        if(Trig == 0){
            m_trig.push_back(*HLT_noalg_L1ZDC_A_AND_C);                            
        }

        else if(Trig ==2){
            m_trig.push_back(*HLT_L1_ZDC_XOR_E2);                                 
        }
        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        for(int mod =0; mod<4; mod++){
            if(isBitSet((*BitMask),mod)) h_module[0][mod]->Fill(lumiBlock,ModAmp[mod]);
            if(isBitSet((*BitMask),mod+4)) h_module[1][mod]->Fill(lumiBlock,ModAmp[mod +4]);
        }     
    }
    tmpf->Write();

}
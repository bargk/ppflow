//Ploting the amplitudes of the modules.
//The code has same structrue as CorrFunc.cxx but without the 2pc anaysis.
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/lumiblock/step_func.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::vector<float> zdcWei;
void load_weights();
double lineFunction(double x) {
    return -1.16 * x + 3500;
}

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

std::vector<float> GetMoments(int ilb_bin, int side, int Trig);

std::string weights_path = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide";
void zdc_amplitude(const int a ,const char* fileList, int Trig1 = 0){
    

    int Trig = Trig1;
    load_weights();
    char output_name[100];
    sprintf(output_name,"histograms_%d.root",a);

    if(Trig ==0){
        cout << "working on ZDC AND trigger!" << endl;
    }
    else if(Trig ==1){
        cout << "Working on minbias trigger!" << endl;
    }
    else if(Trig == 2){
        cout << "Working on ZDC xor trigger!" << endl;
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
    char base[200] ="/gpfs0/citron/users/bargl/ZDC/lhcf22/user.bglik.data22_13p6TeV.00435229.ANALYSIS.2.noLhcf_EXT0.585561170.585561170";
    //char base[200] ="/gpfs0/citron/users/bargl/ZDC/lhcf22/user.bglik.data22_13p6TeV.00435229.reprocZDC.False.ANALYSIS.noLhcf_EXT0.588509852.588509852";
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
    if(Trig == 2) directory = Form("%s/xor",directory.c_str());
    gSystem->Exec(Form("mkdir -p %s",directory.c_str()));

    //----------------------------------coorections for LHCf -------------------------------------------
    //reading moments for late shower cut. DONT FORGET TO UNCOMMENT THE RELEVANT PART IN THE EVENT LOOP!
    std::vector<std::vector<float>> moments_C;
    std::vector<std::vector<float>> moments_A;
    for(int ilb =0; ilb<Bins::NLB; ilb++){
        moments_C.push_back(GetMoments(ilb,0, Trig));
        moments_A.push_back(GetMoments(ilb,1, Trig));
        // cout<<moments_C[ilb].at(0) <<" " <<moments_C[ilb].at(1) <<endl;
        // cout<<moments_A[ilb].at(0) <<" " <<moments_A[ilb].at(1) <<endl;
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
    TFile *tmpf = new TFile(Form("%s/%s",directory.c_str(),output_name),"recreate");
    

    TH2D *hAmpLB[2]; 
    hAmpLB[0]= new TH2D("hLbAmp0",";LB;zdc_ZdcAmp(0) [ADC]", 3800,1000,4800,80,0,8000);
    hAmpLB[1]= new TH2D("hLbAmp1",";LB;zdc_ZdcAmp(1) [ADC]", 3800,1000,4800,80,0,8000);
    TH2D *hAmpCorr= new TH2D("hAmpCorr",";zdc_ZdcAmp(0) [ADC];zdc_ZdcAmp(1) [ADC]", 80,0,8000,80,0,8000);
    TH2D *h_module[2][4];
    TH2D *h_module_bad[2][4];
    TH2D *h_module_maxADC_bad[2][4];
    TH2D *h_had1_corr[2][4];
    TH2D *h_had2_corr_had3[2];
    TH2D *h_h2AmpRat_corr_Amp[2];
    TH2D *h_h2AmpRat_corr_Amp_lb[2][Bins::NLB];
    TH2D *h_h3AmpRat_corr_Amp[2];
    TH2D *h_h1h2Rat_corr_had1[2];
    TH2D *h_h1h2Rat_corr_had2[2];
    TH2D *h_h2h3Rat_corr_had3[2];
    TH2D *h_h1h2Rat_corr_had1_bad[2];
    TH2D *h_h1h2Rat_corr_had2_bad[2];
    TH2D *h_h2h3Rat_corr_had3_bad[2];
    for(int side =0; side<2; side++){
        h_had2_corr_had3[side]= new TH2D(Form("h%i_had2_had_3",side),Form(";zdc_ZdcModuleFitAmp[%i][2];zdc_ZdcModuleFitAmp[%i][3]",side,side),200,0,4600,200,0,4600);
        h_h2AmpRat_corr_Amp[side]= new TH2D(Form("h%i_h2AmpRat_Amp",side),Form(";zdc_ZdcAmp[%i];HAD2/zdc_ZdcAmp[%i]",side,side),80,0,8000,150,0,1.5);
        h_h3AmpRat_corr_Amp[side]= new TH2D(Form("h%i_h3AmpRat_Amp",side),Form(";zdc_ZdcAmp[%i];HAD3/zdc_ZdcAmp[%i]",side,side),80,0,8000,150,0,1.5);
        h_h1h2Rat_corr_had1[side]= new TH2D(Form("h%i_h1h2Rat_corr_had1",side),Form(";zdc_ZdcModuleFitAmp[%i][1];HAD1/HAD2",side),200,0,4600,1600,0,400);
        h_h1h2Rat_corr_had2[side]= new TH2D(Form("h%i_h1h2Rat_corr_had2",side),Form(";zdc_ZdcModuleFitAmp[%i][2];HAD1/HAD2",side),200,0,4600,1600,0,400);
        h_h2h3Rat_corr_had3[side]= new TH2D(Form("h%i_h2h3Rat_corr_had3",side),Form(";zdc_ZdcModuleFitAmp[%i][3];HAD2/HAD3",side),200,0,4600,1600,0,400);
        h_h1h2Rat_corr_had1_bad[side]= new TH2D(Form("h%i_h1h2Rat_corr_had1_bad",side),Form(";zdc_ZdcModuleFitAmp[%i][1];HAD1/HAD2",side),200,0,4600,1600,0,400);
        h_h1h2Rat_corr_had2_bad[side]= new TH2D(Form("h%i_h1h2Rat_corr_had2_bad",side),Form(";zdc_ZdcModuleFitAmp[%i][2];HAD1/HAD2",side),200,0,4600,1600,0,400);
        h_h2h3Rat_corr_had3_bad[side]= new TH2D(Form("h%i_h2h3Rat_corr_had3_bad",side),Form(";zdc_ZdcModuleFitAmp[%i][3];HAD2/HAD3",side),200,0,4600,1600,0,400);
        for(int mod =0; mod<4; mod++){
            h_module[side][mod]= new TH2D(Form("h%i_%i",side,mod),Form(";LB;zdc_ZdcModuleFitAmp[%i][%i]",side,mod),3800,1000,4800,200,0,4600);
            h_module_bad[side][mod]= new TH2D(Form("h_bad%i_%i",side,mod),Form(";LB;zdc_ZdcModuleFitAmp[%i][%i]",side,mod),3800,1000,4800,200,0,4600);
            h_module_maxADC_bad[side][mod]= new TH2D(Form("h_maxADC_bad%i_%i",side,mod),Form(";LB;zdc_ZdcModuleMaxADC[%i][%i]",side,mod),3800,1000,4800,200,0,4600);
            h_had1_corr[side][mod]= new TH2D(Form("h%i_had1_had_%i",side,mod),Form(";zdc_ZdcModuleFitAmp[%i][1];zdc_ZdcModuleFitAmp[%i][%i]",side,side,mod),200,0,4600,200,0,4600);
        }
    }

    for(int side =0; side<2; side++){
        for(int ilb =0; ilb<Bins::NLB; ilb++){
            h_h2AmpRat_corr_Amp_lb[side][ilb]= new TH2D(Form("h%i_h2AmpRat_Amp_ilb_%i",side,ilb),Form(";zdc_ZdcAmp[%i];HAD2/zdc_ZdcAmp[%i]",side,side),80,0,8000,150,0,1.5);
        }
    }

    TH1D *hRatio[2];
    hRatio[0] = new TH1D("hratio_0", ";HAD3/HAD2;Events",1000,0,1000);
    hRatio[1] = new TH1D("hratio_1", ";HAD3/HAD2;Events",1000,0,1000);
    TH1D *hRatio_bad[2];
    hRatio_bad[0] = new TH1D("hratio_bad0", ";HAD3/HAD2;Events",1000,0,1000);
    hRatio_bad[1] = new TH1D("hratio_bad1", ";HAD3/HAD2;Events",1000,0,1000);
    //calibrated histograms
    TH1D *h_energy[8];
    for(int i=0; i<8; i++){
        h_energy[i] = new TH1D(Form("h_energy_%i",i), ";[GeV];Events",200,0,15000); 
    }
    TH1D *hsumA_energy = new TH1D("hsumA_energy", ";[GeV];Events",200,0,25000);
    TH1D *hsumC_energy = new TH1D("hsumC_energy", ";[GeV];Events",200,0,25000);
    TH2D *hAmpCorr_energy = new TH2D("hAmpCorr_energy",";[GeV];Events",200,0,25000,200,0,25000);
    TH1D *hsumA_uncalib = new TH1D("hsumA_uncalib", ";[ADC];Events",80,0,8000);
    TH1D *hsumC_uncalib = new TH1D("hsumC_uncalib", ";[ADC];Events",80,0,8000);
    
    
    //zdc
    TTreeReader myreader(fChain);
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp"); // First 4 lower bits corressponds to side C
    TTreeReaderArray<float> Amp(myreader, "zdc_ZdcAmp"); 
    TTreeReaderArray<float> ModADC(myreader, "zdc_ZdcModuleMaxADC"); 
    TTreeReaderArray<short> LucrodTriggerSideAmp(myreader, "zdc_ZdcLucrodTriggerSideAmp"); 
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");
    TTreeReaderArray<unsigned int> ModStatus(myreader, "zdc_ZdcModuleStatus");
    TTreeReaderValue<unsigned int> lumiblock(myreader, "lumiBlock");
    TTreeReaderValue<int> nvtx(myreader, "t_nvtx");
    //zdc triggers
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_XOR_E1_E3                         (myreader, "HLT_noalg_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_XOR_E1_E3                         (myreader, "HLT_mb_sptrk_L1ZDC_XOR_E1_E3");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_XOR_E2                            (myreader, "HLT_noalg_L1ZDC_XOR_E2");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_XOR_E2                            (myreader, "HLT_mb_sptrk_L1ZDC_XOR_E2");
    TTreeReaderValue<bool> HLT_noalg_L1ZDC_A_AND_C                      (myreader, "HLT_noalg_L1ZDC_A_AND_C");
    TTreeReaderValue<bool> HLT_mb_sptrk_L1ZDC_A_AND_C                   (myreader, "HLT_mb_sptrk_L1ZDC_A_AND_C");

    //minbias triggers
    TTreeReaderValue<bool> HLT_noalg_L1MBTS_1(myreader, "HLT_noalg_L1MBTS_1");
	TTreeReaderValue<bool> HLT_noalg_mb_L1MBTS_2(myreader, "HLT_noalg_mb_L1MBTS_2");
	TTreeReaderValue<bool> HLT_mb_mbts_L1MBTS_2(myreader, "HLT_mb_mbts_L1MBTS_2");

    int LGOverflowBit = 6;
    int HGOverflowBit = 3;
    int nentries= myreader.GetEntries();
    int istart=0,iend= nentries; //by default loop over whole dataset
    while(myreader.Next()){
        std::vector<bool> m_trig;
        int i = myreader.GetCurrentEntry(); 
        if(i>iend) {cout<<"Event loop finished! Analyzed "<< iend<<" events"<<endl;break;}
        int lumiBlock = *lumiblock;
        if(*nvtx ==1) continue;
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


        // if(!(ModAmp[1]/ModAmp[2] > 0.9*func_value && ModAmp[1]/ModAmp[2] < 1.1*func_value)) continue;
        if(Trig == 0){
            m_trig.push_back(*HLT_noalg_L1ZDC_A_AND_C);                            
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_A_AND_C);                            
        } 
        else if(Trig ==1){
            m_trig.push_back(*HLT_noalg_L1MBTS_1);                                                            
            m_trig.push_back(*HLT_noalg_mb_L1MBTS_2);                                                       
             m_trig.push_back(*HLT_mb_mbts_L1MBTS_2);                                                                     
        }
        else if(Trig ==2){
            m_trig.push_back(*HLT_noalg_L1ZDC_XOR_E2);                                 
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_XOR_E2);                                 
            m_trig.push_back(*HLT_noalg_L1ZDC_XOR_E1_E3);                                 
            m_trig.push_back(*HLT_mb_sptrk_L1ZDC_XOR_E1_E3);                                 
        }
        if(!passTrigger(m_trig)) continue; // check if event passed selected triggers
        
        bool overflow = false;
        for(int i=1; i<4; i++){
            if(isBitSet(ModStatus[i],6)){
                overflow = true;
                break;
            }
        }
        if(overflow) continue;
        for(int i=5; i<8; i++){
            if(isBitSet(ModStatus[i],6)){
                overflow = true;
                break;
            }
        }
        if(overflow) continue;

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

        hAmpLB[0]->Fill(lumiBlock, Amp[0]);
        hAmpLB[1]->Fill(lumiBlock, Amp[1]);
        if(isBitSet((*BitMask),2)){
            h_h2AmpRat_corr_Amp[0]->Fill(Amp[0],ModAmp[2]/Amp[0]);
            h_h2AmpRat_corr_Amp_lb[0][ilb_bin]->Fill(Amp[0],ModAmp[2]/Amp[0]);
        }
        if(isBitSet((*BitMask),3)) h_h3AmpRat_corr_Amp[0]->Fill(Amp[0],ModAmp[3]/Amp[0]);
        if(isBitSet((*BitMask),2+4)) {
            h_h2AmpRat_corr_Amp[1]->Fill(Amp[1],ModAmp[2+4]/Amp[1]);
            h_h2AmpRat_corr_Amp_lb[1][ilb_bin]->Fill(Amp[1],ModAmp[2+4]/Amp[1]);
        }
        if(isBitSet((*BitMask),3+4)) h_h3AmpRat_corr_Amp[1]->Fill(Amp[1],ModAmp[3+4]/Amp[1]);

        if(isBitSet((*BitMask),3) && isBitSet((*BitMask),2)) hRatio[0]->Fill(ModAmp[3]/ModAmp[2]);
        if(isBitSet((*BitMask),7) && isBitSet((*BitMask),6)) hRatio[1]->Fill(ModAmp[7]/ModAmp[6]);
        hAmpCorr->Fill(Amp[0], Amp[1]);
        for(int mod =0; mod<4; mod++){
            if(isBitSet((*BitMask),mod)) h_module[0][mod]->Fill(lumiBlock,ModAmp[mod]*scale_stability[0][mod][ilb_bin]);
            if(isBitSet((*BitMask),mod+4)) h_module[1][mod]->Fill(lumiBlock,ModAmp[mod +4]*scale_stability[1][mod][ilb_bin]);
        }
        if(Amp[0] > 4000 && Amp[0] < 5000){
            if(isBitSet((*BitMask),1) && isBitSet((*BitMask),2)){
                h_h1h2Rat_corr_had1_bad[0]->Fill(ModAmp[1],ModAmp[1]/ModAmp[2]);
                h_h1h2Rat_corr_had2_bad[0]->Fill(ModAmp[2],ModAmp[1]/ModAmp[2]);
            }  
            if(isBitSet((*BitMask),1+4) && isBitSet((*BitMask),2+4)){
                h_h1h2Rat_corr_had1_bad[1]->Fill(ModAmp[1+4],ModAmp[1+4]/ModAmp[2+4]);
                h_h1h2Rat_corr_had2_bad[1]->Fill(ModAmp[2+4],ModAmp[1+4]/ModAmp[2+4]);
            }  
            if(isBitSet((*BitMask),3) && isBitSet((*BitMask),2)){
                hRatio_bad[0]->Fill(ModAmp[3]/ModAmp[2]);
                h_h2h3Rat_corr_had3_bad[0]->Fill(ModAmp[3],ModAmp[2]/ModAmp[3]);
            }
            if(isBitSet((*BitMask),7) && isBitSet((*BitMask),6)) {
                hRatio_bad[1]->Fill(ModAmp[7]/ModAmp[6]);
                h_h2h3Rat_corr_had3_bad[1]->Fill(ModAmp[7],ModAmp[6]/ModAmp[7]);
            }
           for(int mod =0; mod<4; mod++){
            if(isBitSet((*BitMask),mod) ) h_module_bad[0][mod]->Fill(lumiBlock,ModAmp[mod]);
            if(isBitSet((*BitMask),mod+4) ) h_module_bad[1][mod]->Fill(lumiBlock,ModAmp[mod +4]);
            if(isBitSet((*BitMask),mod) ) h_module_maxADC_bad[0][mod]->Fill(lumiBlock,ModADC[mod]);
            if(isBitSet((*BitMask),mod+4) ) h_module_maxADC_bad[1][mod]->Fill(lumiBlock,ModADC[mod +4]);
            } 
        }
        if(isBitSet((*BitMask),1) && isBitSet((*BitMask),2)){
            h_had1_corr[0][2]->Fill(ModAmp[1], ModAmp[2]);
            h_h1h2Rat_corr_had1[0]->Fill(ModAmp[1],ModAmp[1]/ModAmp[2]);
            h_h1h2Rat_corr_had2[0]->Fill(ModAmp[2],ModAmp[1]/ModAmp[2]);
        }      
        if(isBitSet((*BitMask),1) && isBitSet((*BitMask),3))h_had1_corr[0][3]->Fill(ModAmp[1], ModAmp[3]);      
        if(isBitSet((*BitMask),1+4) && isBitSet((*BitMask),2+4)){
            h_had1_corr[1][2]->Fill(ModAmp[1+4], ModAmp[2+4]);
            h_h1h2Rat_corr_had1[1]->Fill(ModAmp[1+4],ModAmp[1+4]/ModAmp[2+4]);
            h_h1h2Rat_corr_had2[1]->Fill(ModAmp[2+4],ModAmp[1+4]/ModAmp[2+4]);
        }      
        if(isBitSet((*BitMask),1+4) && isBitSet((*BitMask),3+4))h_had1_corr[1][3]->Fill(ModAmp[1+4], ModAmp[3+4]);
        
        if(isBitSet((*BitMask),2) && isBitSet((*BitMask),3)){
            h_had2_corr_had3[0]->Fill(ModAmp[2], ModAmp[3]);
            h_h2h3Rat_corr_had3[0]->Fill(ModAmp[3],ModAmp[2]/ModAmp[3]);
        }      
        if(isBitSet((*BitMask),6) && isBitSet((*BitMask),7)){
            h_had2_corr_had3[1]->Fill(ModAmp[6], ModAmp[7]);
            h_h2h3Rat_corr_had3[1]->Fill(ModAmp[7],ModAmp[6]/ModAmp[7]); 
        }     
        float sumA =0;
        float sumC =0;
        float sumA_uncalib =0;
        float sumC_uncalib =0;
        bool fillC = false;
        bool fillA = false;
        for(int i=0; i<4; i++){
            if(isBitSet(*BitMask,i)){
                sumC += ModAmp[i]*zdcWei.at(i)*scale_stability[1][i][ilb_bin];
                sumC_uncalib += ModAmp[i]*scale_stability[1][i][ilb_bin];
                h_energy[i]->Fill(ModAmp[i]*zdcWei.at(i)*scale_stability[0][i][ilb_bin]);
                fillC = true;
            }  
            if(isBitSet(*BitMask,i+4)){
                sumA += ModAmp[i+4]*zdcWei.at(i+4)*scale_stability[1][i][ilb_bin];
                sumA_uncalib += ModAmp[i+4]*scale_stability[1][i][ilb_bin];
                h_energy[i+4]->Fill(ModAmp[i+4]*zdcWei.at(i+4)*scale_stability[1][i][ilb_bin]);
                fillA = true;
            }  
        }
        if(fillA) {
            hsumA_energy->Fill(sumA);
            hsumA_uncalib->Fill(sumA_uncalib);
        }
        if(fillC) {
            hsumC_energy->Fill(sumC);
            hsumC_uncalib->Fill(sumC_uncalib);
        }
        hAmpCorr_energy->Fill(sumC,sumA);
        
    }
    tmpf->Write();

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





#include "calibration.h"

void centroid_LB(){
    output = new TFile("centroid.root","recreate");
    //load files
    TChain *tc = new TChain("zdcTree");
    tc->Add(Form("%s*",path.c_str()),0);
    //tc->Add(Form("%suser.bglik.41707724.EXT0._000005.ZDCNT.root",path.c_str()),0);
    if (!tc || tc->IsZombie() || !tc->GetNbranches()) {
    std::cerr << "Error opening file" << std::endl;
    exit(-1);
    }

    TH2D *h0_centroidY_lb = new TH2D("h0_centroid_y", ";LB;centroid_y", 500,200,700, 100,-50,50);
    TH2D *h1_centroidY_lb = new TH2D("h1_centroid_y", ";LB;centroid_y", 500,200,700, 100,-50,50);
    TH2D *h0_Amp_lb = new TH2D("h0_Amp", ";LB;Zdc_zdcAmp", 500,200,700, 200,0,10000);
    TH2D *h1_Amp_lb = new TH2D("h1_Amp", ";LB;Zdc_zdcAmp", 500,200,700, 200,0,10000);
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
    TTreeReaderArray<float> CentroidY(myreader, "zdc_yCentroid");
    TTreeReaderArray<float> Amp(myreader, "zdc_ZdcAmp");
    TTreeReaderValue<unsigned int> BitMask(myreader, "zdc_ZdcModuleMask");
    TTreeReaderArray<float> ModAmp(myreader, "zdc_ZdcModuleFitAmp");


    int nentries= myreader.GetEntries();
    int istart=0,iend= nentries; //by default loop over whole dataset
    cout << "Total events : "<< nentries << endl;
    
    while (myreader.Next()){
        int lumiBlock = *lumiblock;
        if(!((lumiBlock > 350 && lumiBlock < 675))) continue; //stable beams
        int ievt = myreader.GetCurrentEntry();
        if(ievt%100000==0) cout<<"proccesed "<<ievt<<" / "<<iend<<" events "<<" "<<tc->GetFile()->GetName()<<endl;
        if(*HLT_noalg_ZDCPEB_L1ZDC_C){
            h0_centroidY_lb->Fill(lumiBlock,CentroidY[0]);
            h0_Amp_lb->Fill(lumiBlock,Amp[0]);
        }

        if(*HLT_noalg_ZDCPEB_L1ZDC_A) {
            h1_centroidY_lb->Fill(lumiBlock,CentroidY[1]);
            h1_Amp_lb->Fill(lumiBlock,Amp[1]);
        }
    }
    output->Write();
    output->Close();
}
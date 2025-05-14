#include "bins.h"

TH1D *fg[Bins::NCENT+Bins::NCENT_ADD][Bins::NTRK + Bins::NTRK_ADD][Bins::NPT1+Bins::NPT1_ADD][Bins::NPT2+Bins::NPT2_ADD][Bins::NCH+Bins::NCH_ADD][Bins::NDETA];
TH1D *bg[Bins::NCENT+Bins::NCENT_ADD][Bins::NTRK + Bins::NTRK_ADD][Bins::NPT1+Bins::NPT1_ADD][Bins::NPT2+Bins::NPT2_ADD][Bins::NCH+Bins::NCH_ADD][Bins::NDETA];
double NTrigs[Bins::NCENT+Bins::NCENT_ADD][Bins::NTRK + Bins::NTRK_ADD][Bins::NPT1+Bins::NPT1_ADD];

/*-----------------------------------------------------------------------------
 *  Makes PTY distributions from the 1D-pair distributions
 *-----------------------------------------------------------------------------*/
void S06_PTY_1D(int Trig =0){
    std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles");
    if(Trig== 0){
        std::cout << "Working on AND trigger!" << std::endl;
        if(Bins::same_side1) base = Form("%s/sameSide",base.c_str());
    }
    else if(Trig == 1){
        std::cout << "Working on Minbias trigger!" << std::endl;
        if(Bins::same_side1) base = Form("%s/sameSide/minbias",base.c_str());
        else{
            base = Form("%s/minbias",base.c_str());
        }
    }
    else if(Trig ==2){
        std::cout << "Working on XOR trigger!" << std::endl;
        if(Bins::same_side1) base = Form("%s/sameSide/xor",base.c_str());
        else{
            base = Form("%s/xor",base.c_str());
        }
    }
    else{
         std::cerr << "Error: no valid trigger index provided" << std::endl;
         exit(-1);
    }
  std::cout << base << endl;
   char name [600];
   char name1[600];
   sprintf(name,"%s/RebinTrk.root",base.c_str());
   TFile *Cent=new TFile(name,"read");

   for (int icent=0; icent<Bins::NCENT+Bins::NCENT_ADD; icent++){
    for (int itrk=0; itrk<Bins::NTRK+Bins::NTRK_ADD; itrk++){
      sprintf(name,"N_trigger_cent%.2d_trk%.2d",icent, itrk);
      TH1D* h_NTrigs=(TH1D*)Cent->Get(name);
      for (int ipt1=0; ipt1<Bins::NPT1+Bins::NPT1_ADD; ipt1++){
        int bin_lo  =h_NTrigs->FindBin(Bins::PT1_LO[ipt1]+.0001);
        int bin_high=h_NTrigs->FindBin(Bins::PT1_HI[ipt1]-.0001);
        std::cout<<bin_lo<<"  "<<bin_high    <<"  "
        <<h_NTrigs->GetBinLowEdge(bin_lo)    <<"  "
        <<h_NTrigs->GetBinLowEdge(bin_high+1)<<"  "
        <<Bins::PT1_LO[ipt1]                 <<"  "
        <<Bins::PT1_HI[ipt1]<<std::endl;
        if(fabs(h_NTrigs->GetBinLowEdge(bin_lo) - Bins::PT1_LO[ipt1])>0.001){
          std::cout<<"1  "<<h_NTrigs->GetBinLowEdge(bin_lo    )<<"  "
                    <<Bins::PT1_LO[ipt1]<<std::endl;
          throw std::exception();
        }
        if(fabs(h_NTrigs->GetBinLowEdge(bin_high+1) - Bins::PT1_HI[ipt1])>0.001){
          std::cout<<"2  "<<h_NTrigs->GetBinLowEdge(bin_high+1)<<"  "
                    <<Bins::PT1_HI[ipt1]<<std::endl;
          throw std::exception();
        }
        NTrigs[icent][itrk][ipt1]=h_NTrigs->Integral(bin_lo,bin_high);
        //cout<< "Ntriggers: " << NTrigs[icent][itrk][ipt1] <<endl;
      }
    }
   }


   sprintf(name ,"%s/ProjectionX.root",base.c_str());
   sprintf(name1,"%s/PTY1D.root"      ,base.c_str());
   TFile *input =new TFile(name );
   TFile *output=new TFile(name1,"recreate");

  char PjXfgname[100],PjXbgname[100],PjXconame[100];

  const std::vector<int> cent_bins=Bins::CentBins();
  const std::vector<int> trk_bins=Bins::TrkBins();
  const std::vector<int> pt1_bins =Bins::PtaBins ();
  const std::vector<int> pt2_bins =Bins::PtbBins ();
  const std::vector<int> ch_bins  =Bins::ChBins  ();
  const std::vector<int> deta_bins=Bins::DetaBins();

  //Read In
  input->ReadAll();
  TIter next(input->GetList());
  TObject *obj;
  for(int icent:cent_bins){
    for(int itrk:trk_bins){
      for(int ipt1:pt1_bins){
        for(int ipt2:pt2_bins){
          //std::cout<<icent<<"  "<<ipt1<<"  "<<ipt2<<"  ::"<<std::endl;
          for(int ich:ch_bins){
            for(int ieta:deta_bins){
              sprintf(PjXfgname,"PjX_fg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,itrk,ipt1,ipt2,ich,ieta);
              sprintf(PjXbgname,"PjX_bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,itrk,ipt1,ipt2,ich,ieta);
              
              obj=next();
              fg[icent][itrk][ipt1][ipt2][ich][ieta] = (TH1D*)input->Get(PjXfgname);
              if(!Common::CheckObject(obj,PjXfgname)) throw std::exception();
              //if(icent ==11 &&ipt1 == 5 && ipt2 == 5&& itrk==14) cout<< fg[icent][itrk][ipt1][ipt2][ich][ieta]->GetBinError(1)  <<" " <<PjXfgname<<endl;
              obj=next();
              bg[icent][itrk][ipt1][ipt2][ich][ieta] = (TH1D*)obj;
              if(!Common::CheckObject(obj,PjXbgname)) throw std::exception();


              double integral1=fg[icent][itrk][ipt1][ipt2][ich][ieta]->Integral("width");
              fg[icent][itrk][ipt1][ipt2][ich][ieta]->Divide(bg[icent][itrk][ipt1][ipt2][ich][ieta]);
              double integral2=fg[icent][itrk][ipt1][ipt2][ich][ieta]->Integral("width");
              double width    =fg[icent][itrk][ipt1][ipt2][ich][ieta]->GetBinWidth(1);
    

              double scale_correct=integral1/integral2/width/NTrigs[icent][itrk][ipt1];
              double scale_approx =bg[icent][itrk][ipt1][ipt2][ich][ieta]->Integral()/
                                  bg[icent][itrk][ipt1][ipt2][ich][ieta]->GetNbinsX()/width/NTrigs[icent][itrk][ipt1];
              //fg[icent][ipt1][ipt2][ich][ieta]->Scale(integral1/integral2/width/NTrigs[icent][ipt1]);
              fg[icent][itrk][ipt1][ipt2][ich][ieta]->Scale(scale_correct);
              

              //check for nan
              if(fg[icent][itrk][ipt1][ipt2][ich][ieta]->Integral() != fg[icent][itrk][ipt1][ipt2][ich][ieta]->Integral()){
                fg[icent][itrk][ipt1][ipt2][ich][ieta]->Reset();
              }//clear nan histogram (redundant?)

              sprintf(name,"PTY_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,itrk,ipt1,ipt2,ich,ieta);
              fg[icent][itrk][ipt1][ipt2][ich][ieta]->SetName(name);
              fg[icent][itrk][ipt1][ipt2][ich][ieta]->Write();
            }
          }
        }
      }
    }
  }
  output->Close();
  std::cout << "finished" << std::endl;
    sprintf(name,"kill -9 %d",gSystem->GetPid());
    std::cout<<name<<std::endl;
    gSystem->Exec(name);
}

#include "bins.h"

TH1D *fg[Bins::NCENT+Bins::NCENT_ADD][Bins::NPT1+Bins::NPT1_ADD][Bins::NPT2+Bins::NPT2_ADD][Bins::NCH+Bins::NCH_ADD][Bins::NDETA];
TH1D *bg[Bins::NCENT+Bins::NCENT_ADD][Bins::NPT1+Bins::NPT1_ADD][Bins::NPT2+Bins::NPT2_ADD][Bins::NCH+Bins::NCH_ADD][Bins::NDETA];
double NTrigs[Bins::NCENT+Bins::NCENT_ADD][Bins::NPT1+Bins::NPT1_ADD];

/*-----------------------------------------------------------------------------
 *  Makes PTY distributions from the 1D-pair distributions
 *-----------------------------------------------------------------------------*/
void S05_PTY_1D(int m_use_multiplicity = 0){
   string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
   if(m_use_multiplicity == 1) base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/multiplicity";
   char name [600];
   char name1[600];
   sprintf(name,"%s/RebinEff.root",base.c_str());
   TFile *Cent=new TFile(name,"read");

   for (int icent=0; icent<Bins::NCENT+Bins::NCENT_ADD; icent++){
     sprintf(name,"N_trigger_cent%.2d",icent);
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
       NTrigs[icent][ipt1]=h_NTrigs->Integral(bin_lo,bin_high);
     }
   }


   sprintf(name ,"%s/ProjectionX.root",base.c_str());
   sprintf(name1,"%s/PTY1D.root"      ,base.c_str());
   TFile *input =new TFile(name );
   TFile *output=new TFile(name1,"recreate");

   char PjXfgname[100],PjXbgname[100],PjXconame[100];

  const std::vector<int> cent_bins=Bins::CentBins();
  const std::vector<int> pt1_bins =Bins::PtaBins ();
  const std::vector<int> pt2_bins =Bins::PtbBins ();
  const std::vector<int> ch_bins  =Bins::ChBins  ();
  const std::vector<int> deta_bins=Bins::DetaBins();

  //Read In
  input->ReadAll();
  TIter next(input->GetList());
  TObject *obj;
  for(int icent:cent_bins){
    for(int ipt1:pt1_bins){
      for(int ipt2:pt2_bins){
        std::cout<<icent<<"  "<<ipt1<<"  "<<ipt2<<"  ::"<<std::endl;
        for(int ich:ch_bins){
          for(int ieta:deta_bins){
            sprintf(PjXfgname,"PjX_fg_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,ipt1,ipt2,ich,ieta);
            sprintf(PjXbgname,"PjX_bg_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,ipt1,ipt2,ich,ieta);
          //sprintf(PjXconame,"PjX_co_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,ipt1,ipt2,ich,ieta);

            obj=next();
            fg[icent][ipt1][ipt2][ich][ieta] = (TH1D*)obj;
            if(!Common::CheckObject(obj,PjXfgname)) throw std::exception();

            obj=next();
            bg[icent][ipt1][ipt2][ich][ieta] = (TH1D*)obj;
            if(!Common::CheckObject(obj,PjXbgname)) throw std::exception();

          //obj=next();
          //TH1D* htemp = (TH1D*)obj;
          //if(!Common::CheckObject(obj,PjXconame)) throw std::exception();

            double integral1=fg[icent][ipt1][ipt2][ich][ieta]->Integral();
                             fg[icent][ipt1][ipt2][ich][ieta]->Divide(bg[icent][ipt1][ipt2][ich][ieta]);
            double integral2=fg[icent][ipt1][ipt2][ich][ieta]->Integral();
            double width    =fg[icent][ipt1][ipt2][ich][ieta]->GetBinWidth(1);

            double scale_correct=integral1/integral2/width/NTrigs[icent][ipt1];
            double scale_approx =bg[icent][ipt1][ipt2][ich][ieta]->Integral()/
                                 bg[icent][ipt1][ipt2][ich][ieta]->GetNbinsX()/width/NTrigs[icent][ipt1];
            //fg[icent][ipt1][ipt2][ich][ieta]->Scale(integral1/integral2/width/NTrigs[icent][ipt1]);
            fg[icent][ipt1][ipt2][ich][ieta]->Scale(scale_approx);
            if(icent>=2 && (icent<=13 || icent >=18) && ieta==1 && ich==2){
              std::cout<<"A "<<icent<<"  "<<ipt1<<"  "<<ipt2<<" "<<scale_correct/scale_approx<<std::endl;
            }

            //check for nan
            if(fg[icent][ipt1][ipt2][ich][ieta]->Integral() != fg[icent][ipt1][ipt2][ich][ieta]->Integral()){
              fg[icent][ipt1][ipt2][ich][ieta]->Reset();
            }//clear nan histogram (redundant?)

            sprintf(name,"PTY_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,ipt1,ipt2,ich,ieta);
            fg[icent][ipt1][ipt2][ich][ieta]->SetName(name);
            fg[icent][ipt1][ipt2][ich][ieta]->Write();
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

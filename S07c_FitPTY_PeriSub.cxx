#include "bins.h"

/*-----------------------------------------------------------------------------
 *  Evaluates the peripheral-subtracted vnn
 *-----------------------------------------------------------------------------*/
void S07c_FitPTY_PeriSub(int m_use_multiplicity =0){
    string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    if(m_use_multiplicity ==1) base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/multiplicity";

    const std::vector<int> cent_bins=Bins::CentBins();
    const std::vector<int> pt1_bins =Bins::PtaBins ();
    const std::vector<int> pt2_bins =Bins::PtbBins ();
    const std::vector<int> ch_bins  =Bins::ChBins  ();
    const std::vector<int> deta_bins=Bins::DetaBins();
    const std::vector<int> centbins_peripheral=Bins::CentBinsPeriph();

    char name [600];
    char name1[600];
    char name2[600];



    //Do the Peripheral subtraction to get the subtracted PTY
    //------------------------------------------------------------------------------
    {
    sprintf(name ,"%s/PTY1D.root"           ,base.c_str());
    sprintf(name1,"%s/ZYAM1D.root"          ,base.c_str());
    sprintf(name2,"%s/PeripheralScales.root",base.c_str());
    TFile *InFileCentral    = new TFile(name ,"read");
    TFile *InFilePeripheral = new TFile(name1,"read");
    TFile *InFileScales     = new TFile(name2,"read");

    sprintf(name ,"%s/PTYPeripheralSub.root"     ,base.c_str());
    TFile *OutFile1 = new TFile(name ,"recreate");
    OutFile1->cd();
    std::cout<<name<<"  "<<name1<<std::endl;


    for(int ipt1:pt1_bins){
      for(int ipt2:pt2_bins){
        for(int ich:ch_bins){
          for(int ideta:deta_bins){
           std::cout<<ipt1<<"  "<<ipt2<<"  "<<ich<<"  "<<ideta<<std::endl;
            for(auto icent2:centbins_peripheral){

              sprintf(name,"PTY_ZYAM_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent2,ipt1,ipt2,ich,ideta);
              const TH1D* h_peripheral   = (TH1D*)InFilePeripheral->Get(name);
              if(!h_peripheral) {std::cout<<name<<" Not Found"<<std::endl;throw std::exception();}

              sprintf(name,"h_scale_pericent%.2d_pta%d_ptb%.2d_ch%d",icent2,ipt1,ipt2,ich);
              const TH1D* h_scale   = (TH1D*)InFileScales->Get(name);
              if(!h_scale) {std::cout<<name<<" Not Found"<<std::endl;throw std::exception();}

              for(int icent1:cent_bins){
                float scale=h_scale->GetBinContent(icent1+1);

                sprintf(name,"PTY_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d"     ,icent1,ipt1,ipt2,ich,ideta);
                const TH1D* h_central= (TH1D*)InFileCentral->Get(name);
                if(!h_central){std::cout<<name<<" Not Found"<<std::endl;throw std::exception();}


                sprintf(name,"PTY_periSub_cent%.2d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent1,icent2,ipt1,ipt2,ich,ideta);
                TH1D* h_PeriSub=(TH1D*)h_central->Clone(name);
                h_PeriSub->Add(h_peripheral,-scale);


                h_PeriSub->Write();
              }
            }
          }
        }
      }
    }
    InFileCentral   ->Close();
    InFilePeripheral->Close();
    InFileScales    ->Close();
    OutFile1->Close();
    }
    //------------------------------------------------------------------------------






    //Obtain the vnn from the subtracted PTY
    //------------------------------------------------------------------------------
    {
    sprintf(name ,"%s/PTYPeripheralSub.root"  ,base.c_str());
    sprintf(name1,"%s/PeripheralSub_vnn.root" ,base.c_str());
    TFile *InFile1  = new TFile(name ,"read");
    TFile *OutFile2 = new TFile(name1,"recreate");

    TF1* FullFit=new TF1("FullFit","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))",-Common::PI/2,3*Common::PI/2.0);

    TH1D* h_v22=new TH1D("h_v22",";Centbin;v22",Bins::NCENT+Bins::NCENT_ADD,0,Bins::NCENT+Bins::NCENT_ADD);
    TH1D* h_v33=new TH1D("h_v33",";Centbin;v33",Bins::NCENT+Bins::NCENT_ADD,0,Bins::NCENT+Bins::NCENT_ADD);
    TH1D* h_v44=new TH1D("h_v44",";Centbin;v44",Bins::NCENT+Bins::NCENT_ADD,0,Bins::NCENT+Bins::NCENT_ADD);

#ifdef OLD
    InFile1->ReadAll();
    std::cout<<"Reading  Done"<<std::endl;
    TIter next(InFile1->GetList());
    TObject *obj;
#endif
    for(int ipt1:pt1_bins){
      for(int ipt2:pt2_bins){
        std::cout<<"Fourier  "<<"  "<<ipt1<<"  "<<ipt2<<"  ::"<<std::endl;
        for(int ich:ch_bins){
          for(int ideta:deta_bins){
            for(auto icent2:centbins_peripheral){
              h_v22->Reset();
              h_v33->Reset();
              h_v44->Reset();
              for(int icent:cent_bins){
                sprintf(name,"PTY_periSub_cent%.2d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,icent2,ipt1,ipt2,ich,ideta);
#ifdef OLD
                obj=next();
                if(!Common::CheckObject(obj,name)) throw std::exception();
                TH1D* h_PTYPeriSub = (TH1D*)obj;
#else
                TH1D* h_PTYPeriSub=(TH1D*)InFile1->Get(name);
                if(!Common::CheckObject(h_PTYPeriSub,name)) throw std::exception();
#endif

                if(h_PTYPeriSub->Integral()<0.0001) continue;
                if(h_PTYPeriSub->Integral()!=h_PTYPeriSub->Integral()) continue;

                FullFit->SetParameter(0,h_PTYPeriSub->GetBinContent(1));
                h_PTYPeriSub->Fit(FullFit,"Q");

                double v2    =FullFit->GetParameter(2);
                double v2_err=FullFit->GetParError (2);

                double v3    =FullFit->GetParameter(3);
                double v3_err=FullFit->GetParError (3);

                double v4    =FullFit->GetParameter(4);
                double v4_err=FullFit->GetParError (4);

                h_v22->SetBinContent(icent+1,v2    );
                h_v22->SetBinError  (icent+1,v2_err);
                h_v33->SetBinContent(icent+1,v3    );
                h_v33->SetBinError  (icent+1,v3_err);
                h_v44->SetBinContent(icent+1,v4    );
                h_v44->SetBinError  (icent+1,v4_err);
              }
              sprintf(name,"h_v22_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent2,ipt1,ipt2,ich,ideta);
              h_v22->SetName(name);
              h_v22->Write();
              sprintf(name,"h_v33_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent2,ipt1,ipt2,ich,ideta);
              h_v33->SetName(name);
              h_v33->Write();
              sprintf(name,"h_v44_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent2,ipt1,ipt2,ich,ideta);
              h_v44->SetName(name);
              h_v44->Write();
            }
          }
        }
      }
    }
    OutFile2->Close();
    }
    //------------------------------------------------------------------------------



/*
    //obtain the vn from the vnn
    //The vn is vn(cent,pta) obtained by dividing vnn(cent,pta,ptb) by sqrt(vn(cent,ptb,ptb))
    //------------------------------------------------------------------------------
    {
    sprintf(name1,"01RootFiles/%s_PeripheralSub_vnn.root",base.c_str());
    TFile *outfile = new TFile(name1,"update");
    outfile->cd();

    for(int ihar:{2,3,4}){
      for(int ich:ch_bins){
        for(int ideta:deta_bins){
          for(int ipt2:pt2_bins){
            for(auto icent2:centbins_peripheral){

              int ipt1_=Bins::GetPtaIndexForPtbIndex(ipt2);//code exits if pt1 bin corresponding to pt2 bin is not found
              sprintf(name,"h_v%d%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",ihar,ihar,icent2,ipt1_,ipt2,ich,ideta);
              sprintf(name1,"h_v%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d" ,ihar     ,icent2,ipt1_,ipt2,ich,ideta);
              TH1* h_vnn_   =(TH1*)(outfile->Get(name));
              Common::CheckObject(h_vnn_,name);
              TH1* h_vn_diag=(TH1*)(h_vnn_)->Clone(name1);
              Common::Take_Sqrt(h_vn_diag);
              h_vn_diag->Write();

              for(int ipt1:pt1_bins){
                if(ipt1==ipt1_) continue;
                sprintf(name,"h_v%d%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",ihar,ihar,icent2,ipt1,ipt2,ich,ideta);
                sprintf(name1,"h_v%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d" ,ihar     ,icent2,ipt1,ipt2,ich,ideta);
                TH1* h_vn=(TH1*)(outfile->Get(name))->Clone(name1);
                h_vn->Divide(h_vn_diag);
                h_vn->Write();
              }
            }
          }
        }
      }
    }
    outfile->Close();
    }
    //------------------------------------------------------------------------------
*/
    std::cout << "finished" << std::endl;
    sprintf(name,"kill -9 %d",gSystem->GetPid());
    std::cout<<name<<std::endl;
    gSystem->Exec(name);
}

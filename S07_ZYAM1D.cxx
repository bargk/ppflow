
#include "bins.h"
#include "common.C"

#define CORRECT_METHOD //Proper ZYAM error propagation
#define QUADRATURE_ADD_ZYAM
#define REDUCE_ZYAM_ERROR //Reduce ZYAM error by factor of 2 as for unfolded PTY it is counted twice

/*-----------------------------------------------------------------------------
 *  Makes ZYAM subtraction on the PTYs
 *-----------------------------------------------------------------------------*/
void S07_ZYAM1D(){
   string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
   char name [600];

   sprintf(name,"%s/PTY1D.root"         ,base.c_str());
   TFile *input=new TFile(name);
   sprintf(name,"%s/PTY1D_FITS.root"    ,base.c_str());
   TFile *output1=new TFile(name,"recreate");
   sprintf(name,"%s/ZYAM1D.root"        ,base.c_str());
   TFile *output2=new TFile(name,"recreate");
   sprintf(name,"%s/ZYAM1D_Summary.root",base.c_str());
   TFile *output3=new TFile(name,"recreate");


  //Read In
  input->ReadAll();
  TIter next(input->GetList());
  TObject *obj;

  char name2[600];
  sprintf(name2,"[0]*(1");
  const int NHAR_ZYAM_FIT=5;
  for(int i=0;i<NHAR_ZYAM_FIT;i++){
    sprintf(name ,"%s +2*[%d]*cos(%d*x) + 2*[%d]*sin(%d*x)",name2,2*i+1,i+1,2*i+2,i+1);
    sprintf(name2,"%s",name);
  }
  // name of the full fit :[0]*(1 + 2*[1]*cos(1*x) + 2*[2]*sin(1*x) + 2*[3]*cos(2*x) + 2*[4]*sin(2*x) + ... + 2*[9]*cos(5*x) + 2*[10]*sin(5*x))
  sprintf(name2,"%s)",name);
  std::cout<<name2<<std::endl;

  TF1 *FullFit=new TF1("FullFit",name2,-Common::PI/2.0,3*Common::PI/2.0);
  TF1 *MinVal =new TF1("MinVal" ,"[0]",-Common::PI/2.0,3*Common::PI/2.0);
  FullFit->SetLineColor(2        );
  FullFit->SetLineWidth(2        );
  MinVal ->SetLineColor(kOrange-3);
  MinVal ->SetLineWidth(2        );
  MinVal ->SetLineStyle(2        );

  TH1D *h_pedestal=0;// histogram used for pedestal substraction
  TH1D *h_ZYAM_Summary=new TH1D("h_ZYAM_Summary","",Bins::BINS_SUMMARY,0,Bins::BINS_SUMMARY);//Integrated near-side yields
  for(int i=1;i<=Bins::BINS_SUMMARY;i++) {h_ZYAM_Summary->GetXaxis()->SetBinLabel(i,Bins::BINS_SUMMARY_LABELS[i].c_str());}

  const std::vector<int> cent_bins=Bins::CentBins();
  const std::vector<int> trk_bins=Bins::TrkBins();
  const std::vector<int> pt1_bins =Bins::PtaBins ();
  const std::vector<int> pt2_bins =Bins::PtbBins ();
  const std::vector<int> ch_bins  =Bins::ChBins  ();
  const std::vector<int> deta_bins=Bins::DetaBins();

  for (int icent:cent_bins){
    for (int itrk:trk_bins){
      for (int ipt1:pt1_bins){
        for (int ipt2:pt2_bins){
          std::cout<<icent<<"  "<<itrk<<"  "<<ipt1<<"  "<<ipt2<<"  ::"<<std::endl;
          for (int ich:ch_bins){
              for (int ieta:deta_bins){
                sprintf(name,"PTY_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,itrk,ipt1,ipt2,ich,ieta);
                bool b_write_dummy=false;//Flag to write empty or dummy histograms if PTY histogram is empty

                //retrieve PTY object
                obj=next();
                TH1D* PTY = (TH1D*)obj;
                if(!Common::CheckObject(PTY,name)) throw std::exception();

                //Check if PTY histogram is empty;
                if(PTY->Integral()<0.0001) {
                  std::cout<<name<<"is empty, writing dummy values"<<std::endl;
                  b_write_dummy=true;
                }

                if(!h_pedestal){
                  h_pedestal=(TH1D*)obj->Clone("h_pedestal");
                  h_pedestal->Reset();
                  #ifdef CORRECT_METHOD
                  int NBins=h_pedestal->GetNbinsX();
                  for(int ibin=1;ibin<=NBins;ibin++) {h_pedestal->SetBinContent(ibin,1.0);h_pedestal->SetBinError(ibin,0.0);}
                  #endif
                }



                //-------------------------------------------------------------------------
                //TH1D and Fit before ZYAM
                FullFit->SetParameter(0,PTY->GetBinContent(1));
                for(int i=0;i<NHAR_ZYAM_FIT;i++) {
                  FullFit->SetParameter(2*i+1,0);
                  FullFit->SetParameter(2*i+2,0);
                }
                const int LOW=0,HIGH=1;//lables for Lower,Higher values of X at which min is reached
                double min[2]={0,0},min_err[2]={0,0},min_X[2]={0,0};

                if(b_write_dummy==true){
                  PTY->GetListOfFunctions()->Add((TF1*)FullFit->Clone());
                }
                else{
                  TFitResultPtr r=PTY->Fit(FullFit,"QS");
                  min  [HIGH]  =FullFit->GetMinimum (-.02,Common::PI);
                  min_X[HIGH]  =FullFit->GetMinimumX(-.02,Common::PI);
                  min  [LOW ]  =FullFit->GetMinimum (-Common::PI,0 );
                  min_X[LOW ]  =FullFit->GetMinimumX(-Common::PI,0 );
                  r->GetConfidenceIntervals(2, 1, 1, min_X, min_err, 0.683, false);
                }

                double min_        =(min    [0]+min    [1])/2.0;
                double min_err_    =sqrt(min_err[0]*min_err[0] + min_err[1]*min_err[1]);
                #ifdef REDUCE_ZYAM_ERROR
                min_err_/=2;
                #endif
                MinVal->SetParameter(0,min_    );
                MinVal->SetParError (0,min_err_);
                PTY->GetListOfFunctions()->Add((TF1*)MinVal->Clone());
                output1->cd();
                PTY->Write();
                //-------------------------------------------------------------------------



                //-------------------------------------------------------------------------
                //TH1D after ZYAM
                //if histogram is not empty then fit again after ZYAM
                //if it was empty the old fits are identical before and after ZYAM
                if(b_write_dummy==false){
                  #ifdef CORRECT_METHOD
                  PTY->Add(h_pedestal,-min_);
                  #else
                  h_pedestal->Reset();
                  int NBins=h_pedestal->GetNbinsX();
                  for(int ibin=1;ibin<=NBins;ibin++) {h_pedestal->SetBinContent(ibin,min_);h_pedestal->SetBinError(ibin,min_err_);}
                  PTY->Add(h_pedestal,-1.0);
                  #endif

                  PTY->GetListOfFunctions()->Clear();
                  FullFit->SetParameter(0,0);
                  PTY->Fit(FullFit,"QS");
                }

                sprintf(name,"PTY_ZYAM_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,itrk,ipt1,ipt2,ich,ieta);
                PTY->SetName(name);
                output2->cd();
                PTY->Write();
                //-------------------------------------------------------------------------


                /*
                //-------------------------------------------------------------------------
                //TF1 after ZYAM
                sprintf(name,"FitResult_ZYAM_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,ipt1,ipt2,ich,ieta);
                FullFit->SetName(name);
                FullFit->Write();
                //-------------------------------------------------------------------------
                */



                //-------------------------------------------------------------------------
                //Integrated near-side yields and other summary values
                h_ZYAM_Summary->Reset();
                if(b_write_dummy==false){
                  int bin1=PTY->FindBin(min_X[LOW ]);
                  int bin2=PTY->FindBin(min_X[HIGH]);
                  double hist_err     =0.0;//set in line below
                  double hist_integral=PTY->IntegralAndError(bin1,bin2,hist_err,"width");
                  double func_err     =FullFit->IntegralError(min_X[LOW],min_X[HIGH]);
                  double func_integral=FullFit->Integral     (min_X[LOW],min_X[HIGH]);
                  #ifdef CORRECT_METHOD
                  double par=0.0;
                  double cov=min_err_*min_err_;
                  double zyam_err_func=fabs(MinVal->IntegralError(min_X[LOW],min_X[HIGH],&par,&cov));
                  double bin1_lo=PTY->GetBinLowEdge(bin1  );
                  double bin2_hi=PTY->GetBinLowEdge(bin2+1);
                  double zyam_err_hist=fabs(MinVal->IntegralError(bin1_lo,bin2_hi,&par,&cov));
                  #ifdef QUADRATURE_ADD_ZYAM
                    func_err=sqrt(func_err*func_err +(zyam_err_func*zyam_err_func));
                    hist_err=sqrt(hist_err*hist_err +(zyam_err_hist*zyam_err_hist));
                  #else
                    func_err=fabs(func_err)+zyam_err_func;
                    hist_err=fabs(hist_err)+zyam_err_hist;
                  #endif
                  #endif
                  h_ZYAM_Summary->SetBinContent(Bins::BIN_ZYAM_MIN_ERR  ,min_err_     );
                  h_ZYAM_Summary->SetBinError  (Bins::BIN_ZYAM_MIN_ERR  ,0.0          );

                  h_ZYAM_Summary->SetBinContent(Bins::BIN_ZYAM_MIN_VAL  ,min_         );
                  h_ZYAM_Summary->SetBinError  (Bins::BIN_ZYAM_MIN_VAL  ,0.0          );

                  h_ZYAM_Summary->SetBinContent(Bins::BIN_ZYAM_ERR_HIST ,zyam_err_hist);
                  h_ZYAM_Summary->SetBinError  (Bins::BIN_ZYAM_ERR_HIST ,0.0          );
                  h_ZYAM_Summary->SetBinContent(Bins::BIN_ZYAM_ERR_FUNC ,zyam_err_func);
                  h_ZYAM_Summary->SetBinError  (Bins::BIN_ZYAM_ERR_FUNC ,0.0          );

                  h_ZYAM_Summary->SetBinContent(Bins::BIN_FUNC_INTEGRAL ,func_integral);
                  h_ZYAM_Summary->SetBinError  (Bins::BIN_FUNC_INTEGRAL ,func_err     );
                  h_ZYAM_Summary->SetBinContent(Bins::BIN_HIST_INTEGRAL ,hist_integral);
                  h_ZYAM_Summary->SetBinError  (Bins::BIN_HIST_INTEGRAL ,hist_err     );
                  h_ZYAM_Summary->SetBinContent(Bins::BIN_FUNC_CHI2_NDOF,FullFit->GetChisquare()/FullFit->GetNDF());
                }


                sprintf(name,"h_ZYAM_Summary_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,itrk,ipt1,ipt2,ich,ieta);
                h_ZYAM_Summary->SetName(name);
                output3->cd();
                h_ZYAM_Summary->Write();
                //-------------------------------------------------------------------------
            }
          }
        }
      }
    }
   }
   output1->Close();
   output2->Close();
   output3->Close();
   std::cout << "finished" << std::endl;
   sprintf(name,"kill -9 %d",gSystem->GetPid());
   std::cout<<name<<std::endl;
   gSystem->Exec(name);
}

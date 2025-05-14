#include "bins.h"

namespace TemplateFitting{

const TH1D*  g_h_central;
const TH1D*  g_h_peripheral;
TH1D*  g_h_template_fit;
double g_ped;

#define FITS_ALLHARS
//#define FITS_HAR2_ONLY
#ifdef FITS_HAR2_ONLY
const int NHAR =1     ;//v2 only
const int NPARS=NHAR+2;
  const int F_INDEX  =1;//index for fit parameter "F"
  const int V22_INDEX=2;//index for fit parameter "v22"
  const int G_INDEX  =NPARS;//index for parameter "G"
#endif

#ifdef FITS_HAR2_AND3_ONLY
const int NHAR=2;//v2-v3
const int NPARS=NHAR+2;
  const int F_INDEX  =1;//index for fit parameter "F"
  const int V22_INDEX=2;//index for fit parameter "v22"
  const int V33_INDEX=3;//index for fit parameter "v33"
  const int G_INDEX  =NPARS;//index for parameter "G"
#endif

#ifdef FITS_ALLHARS
const int NHAR=4;//v2-v5
const int NPARS=NHAR+2;
  const int F_INDEX  =1;//index for fit parameter "F"
  const int V22_INDEX=2;//index for fit parameter "v22"
  const int V33_INDEX=3;//index for fit parameter "v33"
  const int V44_INDEX=4;//index for fit parameter "v44"
  const int G_INDEX  =NPARS;//index for parameter "G"
#endif

Double_t fitf(Double_t *x,Double_t *par){//par[0]=scale;par[1--NHAR]=vnn
  static TF1 *constant = new TF1("constant","[0]",-Common::PI/2,1.5*Common::PI);
  static TF1 *vnn[NHAR]={};

  if(!vnn[0]){
    std::cout<<"Fitf Initializing functions"<<std::endl;
    char name [100];
    char title[100];
    for(int i=0;i<NHAR;i++){
      sprintf(name ,"v%d%d"                ,i+2,i+2);
      sprintf(title,"[0]*(2*[1]*cos(%d*x))",i+2);
      vnn[i]= new TF1(name,title,-Common::PI/2,1.5*Common::PI);
    }
  }

  if(x[0]<-1.4){
    g_h_template_fit->Reset();

    g_h_template_fit->Add(g_h_peripheral,par[0]);
    g_ped = (g_h_central->Integral()-g_h_template_fit->Integral())/g_h_peripheral->GetNbinsX();

    constant->SetParameter(0,g_ped);
    g_h_template_fit->Add(constant);
    for(int i=0;i<NHAR;i++){
      vnn[i]->SetParameters(g_ped,par[i+1]);
      g_h_template_fit->Add(vnn[i]);
    }
  }

  double returnval = g_h_template_fit->GetBinContent(g_h_template_fit->FindBin(x[0]));
  return returnval;
}


void MyChi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag){
  double chisq = 0;
  double delta;
  for (int i=1; i<=g_h_peripheral->GetNbinsX(); i++) {
    Double_t xvalue = g_h_central->GetBinCenter(i);
    delta  = g_h_central->GetBinContent(i)-fitf(&xvalue,par);

    chisq += delta*delta/
            (g_h_central->GetBinError(i)*g_h_central->GetBinError(i)
             +par[0]*par[0]*g_h_peripheral->GetBinError(i)*g_h_peripheral->GetBinError(i));
  }
  int dof = g_h_central->GetNbinsX() - NHAR; 
  f=chisq;
}











struct Fitting{
  TH1* h_central           =NULL;
  TH1* h_rescaledperipheral=NULL;
  TH1* h_fit_func          =NULL;
  TF1* f_pedestal          =NULL;
  TF1* f_vnn_combined      =NULL;
  TF1* f_vnn[NHAR]         ={};
  TCanvas* c1              =NULL;
  TH1* h_pars              =NULL;

  double    parms    [NHAR+2]={};//scale_for_peripheral,vnn[NHAR],pedestal
  Double_t  parms_err[NHAR+2]={};
  Double_t  chi2;


  void cd(){c1->cd();}

  void GetVnnAndError(double& vnn, double& vnn_err,int ihar){
    vnn=-10;
    vnn_err=-10;
    if(ihar<0) return;
    if(ihar>=NHAR) return;
    vnn    =parms    [ihar+1];
    vnn_err=parms_err[ihar+1];
  }

  void GetScaleAndError(double& scale, double& scale_err){
    scale     = parms    [0];
    scale_err = parms_err[0];
  }

  void Clear(){
    if(h_central )           delete h_central;
    if(h_rescaledperipheral) delete h_rescaledperipheral;
    if(h_fit_func)           delete h_fit_func;
    if(f_pedestal)           delete f_pedestal;
    if(f_vnn_combined)       delete f_vnn_combined;
    for(int i=0;i<NHAR;i++){
      if(f_vnn[i]) delete f_vnn[i];
    }
    if(c1)                   delete c1;
    if(h_pars)               delete h_pars;
  }



  void SetName(std::string NewName){
    char name[600];
    sprintf(name,"h_central_%s"           ,NewName.c_str());h_central           ->SetName(name); h_central->SetTitle(" ");
    sprintf(name,"h_rescaledperipheral_%s",NewName.c_str());h_rescaledperipheral->SetName(name);
    sprintf(name,"h_fit_func_%s"          ,NewName.c_str());h_fit_func          ->SetName(name);
    sprintf(name,"f_pedestal_%s"          ,NewName.c_str());f_pedestal          ->SetName(name);
    sprintf(name,"f_vnn_combined_%s"      ,NewName.c_str());f_vnn_combined      ->SetName(name);
    for(int i=0;i<NHAR;i++){
      sprintf(name,"f_v%d%d_%s",i+2,i+2,NewName.c_str());
      f_vnn[i] ->SetName(name);
    }
    sprintf(name,"can_fits_%s"            ,NewName.c_str());c1                  ->SetName(name);
    c1->SetTitle("");
    sprintf(name,"h_pars_%s"              ,NewName.c_str());h_pars              ->SetName(name);
  }



  void Write(){
    h_central           ->Write();
    h_rescaledperipheral->Write();
    h_fit_func          ->Write();
    f_pedestal          ->Write();
    f_vnn_combined      ->Write();
    for(int i=0;i<NHAR;i++) f_vnn[i]->Write();
    c1                  ->Write();
    h_pars              ->Write();
  }




  //Constructor:Create the objects to be drawn
  Fitting(const TH1* l_h_central){
    char name [100];
    char title[100];

      h_central            = (TH1*)l_h_central->Clone("h_central")          ;h_central->Reset();
      h_rescaledperipheral = (TH1*)h_central  ->Clone("h_peripheral_scaled");
      h_fit_func           = (TH1*)h_central  ->Clone("h_fit")              ;

      Common::format(h_central           ,1,20);
      Common::format(h_fit_func          ,2, 0);
      Common::format(h_rescaledperipheral,1,24);
      h_fit_func->SetLineWidth(2);

      f_pedestal           = new TF1("f_pedestal","[0]",-Common::PI/2,1.5*Common::PI);
      f_pedestal->SetLineColor(kOrange-3);
      f_pedestal->SetLineStyle(2);
      f_pedestal->SetLineWidth(2);

      sprintf(title,"[0]*(1");
      for(int i=0;i<NHAR;i++){
        sprintf(name ,"%s+2*[%d]*cos(%d*x)",title,i+1,i+2);
        sprintf(title,"%s",name );
        sprintf(name ,"%s",title);
      }
      sprintf(title,"%s)",name);
      f_vnn_combined=new TF1("f_vnn_combined",title,-Common::PI/2,1.5*Common::PI);
      f_vnn_combined->SetLineColor(4);
      f_vnn_combined->SetLineStyle(7);
      f_vnn_combined->SetLineWidth(2);


      for(int i=0;i<NHAR;i++){
        sprintf(name ,"f_v%d%d",i+2,i+2);
        sprintf(title,"[0]*(1+2*[1]*cos(%d*x))",i+2);
        f_vnn[i]= new TF1(name,title,-Common::PI/2,1.5*Common::PI);
        f_vnn[i]->SetLineColor(3);
        f_vnn[i]->SetLineStyle(i+2);
      }

      sprintf(name,"can_fits");
      c1= new TCanvas(name,name,800,600);
      c1->SetLeftMargin (0.15);
      c1->SetTopMargin  (0.02);
      c1->SetBottomMargin(0.12);
      c1->SetRightMargin(0.02);
      c1->cd();

      h_pars=new TH1D("h_pars",";Parameters;",NHAR+3,0,NHAR+3); //added extra bin for chi2
  }


  void Update(TF1* func,const TH1* l_h_central,const TH1* l_h_peripheral){
    //update the parameters
    for(int i=0;i<NHAR+1;i++){
      parms    [i]=func->GetParameter(i);
      parms_err[i]=func->GetParError (i);
    }
    parms    [NHAR+1]= g_ped;
    parms_err[NHAR+1]= parms_err[0]*(l_h_peripheral->Integral()/l_h_peripheral->GetNbinsX());
    int npar=NHAR+1;
    MyChi2(npar,0,chi2,func->GetParameters(),0);

    for(int i=0;i<=NHAR+1;i++){
      h_pars->SetBinContent(i+1,parms    [i]);
      h_pars->SetBinError  (i+1,parms_err[i]);
    }
    //write chi2
    h_pars->SetBinContent(h_pars->GetNbinsX(),chi2);
    h_pars->SetBinError(h_pars->GetNbinsX(),0);
    //Set parameters for objects to be drawn
    h_central->Reset();
    h_central->Add(l_h_central);

    f_pedestal->SetParameter (0,parms    [NHAR+1]);//must be done before h_rescaledperipheral
    f_pedestal->SetParError  (0,parms_err[NHAR+1]);//must be done before h_rescaledperipheral

    h_rescaledperipheral->Reset();
    h_rescaledperipheral->Add(l_h_peripheral);
    h_rescaledperipheral->Scale(parms[0]);
    h_rescaledperipheral->Add(f_pedestal);

    h_fit_func ->Reset();
    h_fit_func ->Add(func);

    for(int i=0;i<NHAR;i++){
      f_vnn[i]      ->SetParameters(parms[NHAR+1],parms    [i+1]);
      f_vnn[i]      ->SetParError  (0            ,parms_err[NHAR+1]);
      f_vnn[i]      ->SetParError  (1            ,parms_err[i+1]);

      f_vnn_combined->SetParameter (i+1          ,parms    [i+1]);
      f_vnn_combined->SetParError  (i+1          ,parms_err[i+1]);
    }
    f_vnn_combined->SetParameter(0,parms    [NHAR+1]);
    f_vnn_combined->SetParError (0,parms_err[NHAR+1]);


    //Draw all the objects
    c1->cd();
    //h_central ->SetMaximum(1 + 3 * (1 - parms[NHAR + 1] * (1.0 - 2.5 * parms[1])));
    h_central ->SetMinimum(0.95*parms[NHAR+1]*(1.0-2.5*parms[1]));
    h_central           ->Draw();
    h_rescaledperipheral->Draw("same");
    h_fit_func          ->Draw("same");
    f_pedestal          ->Draw("same");
    f_vnn_combined      ->Draw("same");
    //for(int i=0;i<NHAR;i++) f_vnn[i]->Draw("same");
  }

};





Fitting* TemplateFit(const TH1D* h_central, const TH1D* h_peripheral){
  g_h_central      = h_central;
  g_h_peripheral   = h_peripheral;
  g_h_template_fit = (TH1D*)g_h_peripheral->Clone("h_template_fit");
  g_h_template_fit->Reset();


  static Fitting* fitresult = new Fitting(g_h_central);

  static TF1 *func = new TF1("f_fit",fitf,-Common::PI/2,1.5*Common::PI,NHAR+1);
  func->SetParameter(0,1);
  for(int i=1;i<NHAR+1;i++) func->SetParameter(i,0.0);

  fitresult->cd();
  TH1D* l_h_central_=(TH1D*)g_h_central->Clone("l_h_central_");
  TVirtualFitter::Fitter(l_h_central_)->SetFCN(MyChi2);
  //l_h_central_->Fit(func,"BQU");
  if(l_h_central_->Integral() != l_h_central_->Integral() || //check for nan
    l_h_central_->Integral()<=0.0){                          //check for empty histogram
    l_h_central_->Reset();
    for(int i=0;i<NHAR+1;i++) {func->SetParameter(i,0.0);func->SetParError(i,0.0);}
    g_ped=0.0;
  }
  else{
    l_h_central_->Fit(func,"QU");
  }

  fitresult->Update(func,g_h_central,g_h_peripheral);

  delete g_h_template_fit;
  delete l_h_central_;

  return fitresult;
}
}

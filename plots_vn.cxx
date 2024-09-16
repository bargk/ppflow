#include "bins.h"
#include "common.C"
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"


std::vector<TCanvas*>        m_can_vec;
std::map<std::string,double> m_format=Common::StandardFormat();

enum{
 FEW_BINS=0,
 MORE_BINS=1,
};

TStyle* AtlasStyle();
void SetAtlasStyle();

void Peripheraldep_Centdep(int itrk,int ipt1,int ipt2,int ich, int ideta, std::vector<int> centbins_peripheral,std::vector<int> trkbins_peripheral);

TFile *InFile1;
TFile *InFile2;
TFile *OutFile;

int m_data_type;
int m_har;
//int m_correlation_type;
int m_vn_type;
bool m_do_shifts;




void plots_vn(int  ihar        = 2,//ihar=2,3,4 = v2,v3,v4
              int  l_vn_type   = Bins::VN_TEMPLATE,
              bool l_do_shifts = false){
     SetAtlasStyle();             
    std::string base1 = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    //m_correlation_type=l_correlation_type;

    std::string base2 =base1;
    // if(m_correlation_type!=DataSetEnums::HADRON_HADRON_CORRELATIONS){
    //   int l_trig_type       =DataSetEnums::NTRK_CUT_FOR_HMT          ;//overrides outer scope variable
    //   int l_correlation_type=DataSetEnums::HADRON_HADRON_CORRELATIONS;//overrides outer scope variable
    //   base2=BASENAME;
    // }

    m_har             =ihar;
    m_vn_type         =l_vn_type;
    m_do_shifts       =l_do_shifts;
 

    m_format["YTitleOffset"]=1.3;
    m_format["XTitleOffset"]=1.3;
    m_format["XNdivisions"] =505;
    m_format["YNdivisions"] =505;

    gStyle->SetErrorX(0.001);

    char name1[600];
    char name2[600];
    if(m_vn_type==Bins::VN_TEMPLATE){
      sprintf(name1,"%s/TemplateFits_vnn.root" ,base1.c_str());
      sprintf(name2,"%s/TemplateFits_vnn.root" ,base2.c_str());
    }
    else if(m_vn_type==Bins::VN_TEMPLATE_PEDESTAL){
      sprintf(name1,"01RootFiles/%s_TemplateFits_pedestal_vnn.root" ,base1.c_str());
      sprintf(name2,"01RootFiles/%s_TemplateFits_pedestal_vnn.root" ,base2.c_str());
    }
    else if(m_vn_type==Bins::VN_DIRECT){
      sprintf(name1,"01RootFiles/%s_Direct_vnn.root"        ,base1.c_str());
      sprintf(name2,"01RootFiles/%s_Direct_vnn.root"        ,base2.c_str());
    }
    else{
      std::cout<<__PRETTY_FUNCTION__<<" Unsupported vn_type="<<m_vn_type<<std::endl;
      throw std::exception();
    }


                     InFile1=new TFile(name1,"read");Common::CheckFile(InFile1,name1);
    if(base2!=base1){InFile2=new TFile(name2,"read");Common::CheckFile(InFile2,name2);}
    else             InFile2=InFile1;
    OutFile = new TFile("Rootfiles/Plots.root","update");
    std::vector<int>cent_periph=Bins::CentBinsPeriph();
    std::vector<int>trk_periph=Bins::TrkBinsPeriph();
    Peripheraldep_Centdep(Bins::GetTrkIndex(0,1000),Bins::GetPtaIndex(0.5,5.0),Bins::GetPtbIndex(0.5,5.0),2,Bins::GetDetaIndex(2.0,5.0),cent_periph,trk_periph);
    Peripheraldep_Centdep(Bins::GetTrkIndex(30,40),Bins::GetPtaIndex(0.5,5.0),Bins::GetPtbIndex(0.5,5.0),2,Bins::GetDetaIndex(2.0,5.0),cent_periph,trk_periph);
    Peripheraldep_Centdep(Bins::GetTrkIndex(100,110),Bins::GetPtaIndex(0.5,5.0),Bins::GetPtbIndex(0.5,5.0),2,Bins::GetDetaIndex(2.0,5.0),cent_periph,trk_periph);
    base1 = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/vn_plots";
    Common::SaveCanvas(m_can_vec,base1);
}

//Peripheraldep vs Ntrk
void Peripheraldep_Centdep(int itrk,int ipt1,int ipt2,int ich, int ideta, std::vector<int>centbins_peripheral, std::vector<int>trkbins_peripheral){
  if(m_vn_type==Bins::VN_DIRECT) {
    std::cout<<" Function doesnot make sense for this m_vn_type="<<m_vn_type<<std::endl;
    return ;
  }
  std::vector<int> centbins_central;
  for(int i=1; i<10; i++){
  centbins_central.push_back(Bins::GetCentIndex(1.36*i,1.36*(i+1)));
  }


  //std::vector<int> centbins_central=Bins::GetCentIndex({1.36,1.36*2,1.36*3,0.68*10});
  //std::vector<int> centbins_central=Bins::GetCentIndex(points);

  char name [600];
  char name1[600];
  char name2[600];
  sprintf(name1,"can_v%d%d_peripheralbindep_centdep",m_har,m_har);
  sprintf(name2,"trk%.2d_pta%d_ptb%d_ch%d_deta%.2d",itrk,ipt1,ipt2,ich,ideta);
  if     (m_vn_type==Bins::VN_TEMPLATE         ) sprintf(name,"%s_template_%s"        ,name1,name2);
  else if(m_vn_type==Bins::VN_TEMPLATE_PEDESTAL) sprintf(name,"%s_templatepedestal_%s",name1,name2);
  else if(m_vn_type==Bins::VN_PERISUB          ) sprintf(name,"%s_perisub_%s"         ,name1,name2);

  TCanvas *Can1=new TCanvas(name,name,1200,900);
  //TCanvas *Can1=new TCanvas(name,name,1300,600);
  Can1->Divide(2,2);
  m_can_vec.push_back(Can1);


  for(int vnn_vn:{0,1}){
     int idraw=0;
     TH1D* h_denominator=nullptr;
     TH1D* h_vnn_vn_main=nullptr;
     TH1D* h_ratio_main =nullptr;
     float max=-1e3,min=0,max_ratios=-1e3,min_ratios=1e3;

     for(auto icent2:centbins_peripheral){
        static int icent0=icent2;//the first centrality which is used for all ratios
      for(auto itrk2:trkbins_peripheral){
        TH1D* h_vnn_vn=Bins::CentdepHist(centbins_central,Common::UniqueName().c_str()); 
        if(vnn_vn==0) sprintf(name,"v_{%d,%d}(p_{T}^{a},p_{T}^{b})",m_har,m_har);
        else          sprintf(name,"v_{%d}(p_{T}^{b})"             ,m_har      );
        h_vnn_vn->GetYaxis()->SetTitle(name);
        for(unsigned int icent_bin=0;icent_bin<centbins_central.size();icent_bin++){
          //skip centralities smaller than cent_periph
          if((h_vnn_vn->GetBinLowEdge(icent_bin+2)-0.001) < Bins::GetCentVals(icent2).second) continue;
          if((h_vnn_vn->GetBinLowEdge(icent_bin+2)-0.001) < Bins::GetCentVals(icent0).second) continue;

          std::pair<float, float> vnn=Bins::GetVnn  (centbins_central[icent_bin],itrk,ipt1,ipt2,ich,ideta,m_har,icent2,itrk2, InFile1);
          if(vnn.first<0) continue; //remove negative vnn values
          if(vnn_vn==1)           vnn=Bins::GetVnPtb(centbins_central[icent_bin],itrk,ipt1,ipt2,ich,ideta,m_har,icent2,itrk2,InFile1,InFile2);
          h_vnn_vn->SetBinContent(icent_bin+1,vnn.first );
          h_vnn_vn->SetBinError  (icent_bin+1,vnn.second);
        }
        if(h_vnn_vn->GetMaximum()>max) max=1.3*h_vnn_vn->GetMaximum();
        if(h_vnn_vn->GetMinimum()<min) min=h_vnn_vn->GetMinimum();



        int sty[]={20,24,21,25,29     ,30};
        int col[]={ 1, 2, 4, 6,kOrange,40};
        Common::FormatHist(h_vnn_vn,m_format);
        Common::format    (h_vnn_vn,col[idraw],sty[idraw]);

        float shift[]={0,1,-1};
        TH1D* htemp=(TH1D*)h_vnn_vn->Clone(Common::UniqueName().c_str());
        if(m_do_shifts) Common::ShiftXaxis(htemp,shift[idraw]);

        Can1->cd(vnn_vn+1);
        if(idraw==0){
          htemp->Draw();
          h_vnn_vn_main=(TH1D*)htemp;
          h_denominator=(TH1D*)h_vnn_vn;// For ratios
        }
        else{
          htemp->Draw("same");

          TH1D* h_ratio=(TH1D*)h_vnn_vn->Clone(Common::UniqueName().c_str());
          h_ratio->Divide(h_denominator);
          std::string s1=h_ratio->GetYaxis()->GetTitle();s1=s1+" Ratio";
          h_ratio->GetYaxis()->SetTitle(s1.c_str());


          //For ratios, Find the maximum and minimum only over the nchrec>20 and beyond the given peripheral bin region
          for(int ibin=1;ibin<=h_ratio->GetNbinsX();ibin++){
            double min_cent =Bins::GetCentVals(centbins_peripheral[0]).second;
            double min_cent2=Bins::GetCentVals(icent2                ).second;
            if((h_ratio->GetBinLowEdge(ibin)+0.001)>min_cent && (h_ratio->GetBinLowEdge(ibin)+0.001)>min_cent2){
              double val=h_ratio->GetBinContent(ibin);
              if(val>max_ratios) max_ratios=val;
              if(val<min_ratios) min_ratios=val;
            }
            if(h_ratio->GetBinContent(ibin)!=0){
              h_ratio->SetBinError(ibin,0);//remove errorbar from ratio hists
            }
          }

          Can1->cd(vnn_vn+3);
          Can1->GetPad(vnn_vn+3)->SetGrid();
          if(idraw==1) {h_ratio->Draw("");h_ratio_main=h_ratio;}
          else         h_ratio->Draw("same");
        }

          //Draw Legend
         if(vnn_vn==0){
           Can1->cd(1);
           std::string legend=Bins::label_cent_peri(icent2);
           if(idraw<3) Common::myMarkerText(0.51,0.90-idraw*0.06    ,col[idraw],sty[idraw],legend,1,0.040);
           else        Common::myMarkerText(0.73,0.49-(idraw-3)*0.06,col[idraw],sty[idraw],legend,1,0.040);
         }
       idraw++;
      }
     }

     //Set Ranges for Y-axis
     float diff1=(max-min)*0.05;
     float diff2=(max_ratios-min_ratios)*0.05;
     h_vnn_vn_main->GetYaxis()->SetRangeUser(min       -diff1,max       +10*diff1);
     //h_ratio_main ->GetYaxis()->SetRangeUser(min_ratios-diff2,max_ratios+diff2);
     std::cout<<min<<"   "<<max<<"   "<<min_ratios<<"   "<<max_ratios<<endl;

  }


  Can1->cd(1);
  float X=0.20,Y=0.88;
  Common::myText2(X     ,Y,1,"ATLAS"         ,18,73);
  Common::myText2(X+0.12,Y,1,Common::Internal,18,43);Y=Y-0.05;
  Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", 18, 43);Y=Y-0.05;
  if(Bins::GetPtbIndexForPtaIndex(ipt1)==ipt2){
    Common::myText2(X     ,Y,1,Bins::label_ptab(ipt1,ipt2) ,18,43);Y=Y-0.05;
  }
  else{
    Common::myText2(X     ,Y,1,Bins::label_pta (ipt1) ,18,43);Y=Y-0.05;
    Common::myText2(X     ,Y,1,Bins::label_ptb (ipt2) ,18,43);Y=Y-0.05;
  }
  Common::myText2(X     ,Y,1,Bins::label_eta (ideta)       ,18,43);Y=Y-0.05;
  //Common::myText2(X     ,Y,1,Bins::label_vn_type(m_vn_type),18,43);
}

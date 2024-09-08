#include"TemplateFitting.h"


//#define TEST // ONLY do few bins, for debugging

std::vector<TCanvas*>        m_can_vec;


/*-----------------------------------------------------------------------------
 *  Does the template fits and stores the fits as well as the vnn
 *-----------------------------------------------------------------------------*/
void S07a_FitPTY_Template(int l_use_peripheral_pp = 0) {
    bool no_ZYAM = false;
    string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    char name [600];
    char name1[600];

    sprintf(name , "%s/PTY1D.root" , base.c_str());
    sprintf(name1, "%s/ZYAM1D.root", base.c_str());

    if (l_use_peripheral_pp == 1) {
        
        if (no_ZYAM) sprintf(name1, "01RootFiles/Data_5TeV_eff1_trig2_PTY1D.root");
        else         sprintf(name1, "%s/multiplicity/ZYAM1D.root", base.c_str());
    }

    TFile *InFileCentral    = new TFile(name , "read");
    TFile *InFilePeripheral = new TFile(name1, "read");

    sprintf(name , "%s/TemplateFits.root"     , base.c_str());
    sprintf(name1, "%s/TemplateFits_vnn.root" , base.c_str());
    if (l_use_peripheral_pp == 1) {
        if (no_ZYAM) {
            sprintf(name , "01RootFiles/%s_TemplateFits_ppPeriph_pedestal.root"     , base.c_str());
            sprintf(name1, "01RootFiles/%s_TemplateFits_ppPeriph_pedestal_vnn.root" , base.c_str());
        }
        else {
            sprintf(name , "%s/TemplateFits_ntrkPeriph.root"     , base.c_str());
            sprintf(name1, "%s/TemplateFits_ntrkPeriph_vnn.root" , base.c_str());
        }
    }
    TFile *OutFile1 = new TFile(name , "recreate"); OutFile1->cd();
    TFile *OutFile2 = new TFile(name1, "recreate");
    std::cout << name << "  " << name1 << std::endl;




#ifdef TEST
    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> pt1_bins  = {Bins::GetPtaIndex(0.5, 5.0)};
    const std::vector<int> pt2_bins  = {Bins::GetPtbIndex(0.5, 5.0)};
    const std::vector<int> ch_bins   = {2};
    const std::vector<int> deta_bins = {Bins::GetDetaIndex(2.0, 5.0)};
    const std::vector<int> centbins_peripheral = Bins::CentBinsPeriph();
#else
    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> pt1_bins  = Bins::PtaBins ();
    // const std::vector<int> pt1_bins  = {Bins::GetPtaIndex(0.5, 5.0)};
    const std::vector<int> pt2_bins  = Bins::PtbBins ();
    const std::vector<int> ch_bins   = Bins::ChBins  ();
    const std::vector<int> deta_bins = Bins::DetaBins();
    //const std::vector<int> centbins_peripheral = {52, 55, 62, 63};
    const std::vector<int> centbins_peripheral = Bins::CentBinsPeriph();
#endif

    TH1D* h_v22   = new TH1D("h_v22", ";Centbin;v22"    , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_v33   = new TH1D("h_v33", ";Centbin;v33"    , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_v44   = new TH1D("h_v44", ";Centbin;v44"    , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_v55   = new TH1D("h_v55"  , ";Centbin;v55"  , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_scale = new TH1D("h_scale", ";Centbin;Scale", Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);

    //------------------------------------------------------------------------------
    for (int ipt1 : pt1_bins) {
        for (int ipt2 : pt2_bins) {
            for (int ich : ch_bins) {
                for (int ideta : deta_bins) {
                    std::cout << "Vnn " << ipt1 << "  " << ipt2 << "  " << ich << "  " << ideta << std::endl;
                    for (auto icent2 : centbins_peripheral) {
                        for (int icent1 : cent_bins) {

                            char hCentname[100], hPeriname[100];
                            sprintf(hCentname, "PTY_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d"     , icent1, ipt1, ipt2, ich, ideta);
                            sprintf(hPeriname, "PTY_ZYAM_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);

                            if (l_use_peripheral_pp == 1) {
                                if (no_ZYAM) {
                                    sprintf(hPeriname, "PTY_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d"      , icent2, ipt1, ipt2, ich, ideta);
                                }
                                else sprintf(hPeriname, "PTY_ZYAM_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);
                            }
                            TH1D* h_central      = (TH1D*)InFileCentral   ->Get(hCentname)->Clone(Common::UniqueName().c_str());
                            TH1D* h_peripheral   = (TH1D*)InFilePeripheral->Get(hPeriname)->Clone(Common::UniqueName().c_str());
                            if (!h_central   ) {std::cout << hCentname << " Not Found" << std::endl; throw std::exception();}
                            if (!h_peripheral) {std::cout << hPeriname << " Not Found" << std::endl; throw std::exception();}

                            h_central   ->GetXaxis()->SetTitleOffset(1.2);
                            h_central   ->GetXaxis()->SetTitle("#Delta#phi");
                            h_central   ->GetYaxis()->SetTitle("Y(#Delta#phi)");
                            h_peripheral->GetXaxis()->SetTitle("#Delta#phi");
                            h_peripheral->GetYaxis()->SetTitle("Y(#Delta#phi)");

                            /*-----------------------------------------------------------------------------
                             *  Normalize as correlation function (new)
                             *-----------------------------------------------------------------------------*/
                            // {
                            //     double scale1 = h_central->Integral() / h_central->GetNbinsX();
                            //     if (scale1 != 0) {
                            //         h_central   ->Scale(1.0 / scale1);
                            //         h_peripheral->Scale(1.0 / scale1);
                            //     }
                            //     h_central   ->GetYaxis()->SetTitle("C(#Delta#phi)");
                            //     h_peripheral->GetYaxis()->SetTitle("C(#Delta#phi)");
                            // }

                            TemplateFitting::Fitting* fitresult = TemplateFitting::TemplateFit(h_central, h_peripheral);
                            sprintf(name, "cent%.2d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent1, icent2, ipt1, ipt2, ich, ideta);
                            std::string new_name = name;
                            fitresult->SetName(new_name);

                            double v22, v22_err;
                            fitresult->GetVnnAndError(v22, v22_err, 0);
                            h_v22->SetBinContent(icent1 + 1, v22);
                            h_v22->SetBinError  (icent1 + 1, v22_err);

                            double v33, v33_err;
                            fitresult->GetVnnAndError(v33, v33_err, 1); //returns -10 if v33 is not included in fit
                            h_v33->SetBinContent(icent1 + 1, v33);
                            h_v33->SetBinError  (icent1 + 1, v33_err);

                            double v44, v44_err;
                            fitresult->GetVnnAndError(v44, v44_err, 2); //returns -10 if v44 is not included in fit
                            h_v44->SetBinContent(icent1 + 1, v44);
                            h_v44->SetBinError  (icent1 + 1, v44_err);

                            double v55, v55_err;
                            fitresult->GetVnnAndError(v55, v55_err, 3); //returns -10 if v55 is not included in fit
                            h_v55->SetBinContent(icent1 + 1, v55);
                            h_v55->SetBinError  (icent1 + 1, v55_err);

                            double scale, scale_err;
                            fitresult->GetScaleAndError(scale, scale_err); //the scale factor for the peripheral bin
                            if (h_central->Integral() > 0) {
                                scale    *= h_peripheral->Integral() / h_central->Integral();
                                scale_err *= h_peripheral->Integral() / h_central->Integral();
                            }
                            h_scale->SetBinContent(icent1 + 1, scale);
                            h_scale->SetBinError  (icent1 + 1, scale_err);

                            OutFile1->cd();
                            fitresult->Write();

                            delete h_central;
                            delete h_peripheral;
                        }
                        OutFile2->cd();
                        sprintf(name, "h_v22_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);
                        //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v22_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                        h_v22->SetName(name);
                        h_v22->Write();
                        sprintf(name, "h_v33_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);
                        //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v33_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                        h_v33->SetName(name);
                        h_v33->Write();
                        sprintf(name, "h_v44_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);
                        //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v44_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                        h_v44->SetName(name);
                        h_v44->Write();
                        sprintf(name, "h_v55_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);
                        //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v55_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                        h_v55->SetName(name);
                        h_v55->Write();
                        sprintf(name, "h_v00_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);
                        //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v00_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                        h_scale->SetName(name);
                        h_scale->Write();
                    }
                }
            }
        }
    }
    InFileCentral   ->Close();
    InFilePeripheral->Close();
    OutFile1->Close();
    OutFile2->Close();
    


    /* 
        //obtain the vn from the vnn
        //The vn is vn(cent,pta) obtained by dividing vnn(cent,pta,ptb) by sqrt(vn(cent,ptb,ptb))
        sprintf(name1,"%s/TemplateFits_vnn.root" ,base.c_str());
        
        OutFile2 = new TFile(name1,"update");
        OutFile2->cd();

        for(int ihar:{2,3,4}){
        for(int ich:ch_bins){
          for(int ideta:deta_bins){
            for(auto icent2:centbins_peripheral){
              for(int ipt2:pt2_bins){

                int ipt1_=Bins::GetPtaIndexForPtbIndex(ipt2);//code exits if pt1 bin corresponding to pt2 bin is not found
                sprintf(name,"h_v%d%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",ihar,ihar,icent2,ipt1_,ipt2,ich,ideta);
                sprintf(name1,"h_v%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d" ,ihar     ,icent2,ipt1_,ipt2,ich,ideta);
                TH1* h_vnn_   =(TH1*)(OutFile2->Get(name));
                if(!h_vnn_) {std::cout<<name<<std::endl;throw std::exception();}
                TH1* h_vn_diag=(TH1*)(h_vnn_)->Clone(name1);
                Common::Take_Sqrt(h_vn_diag);
                h_vn_diag->Write();

                for(int ipt1:pt1_bins){
                  std::cout<<"Vn "<<ipt1<<"  "<<ipt2<<"  "<<ich<<"  "<<ideta<<std::endl;
                  if(ipt1==ipt1_) continue;
                  sprintf(name,"h_v%d%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",ihar,ihar,icent2,ipt1,ipt2,ich,ideta);
                  sprintf(name1,"h_v%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d" ,ihar     ,icent2,ipt1,ipt2,ich,ideta);
                  TH1* h_vn=(TH1*)(OutFile2->Get(name))->Clone(name1);
                  h_vn->Divide(h_vn_diag);
                  h_vn->Write();
                }
              }
            }
          }
        }
        }
       OutFile2->Close();
    */
    std::cout << "finished" << std::endl;
    sprintf(name,"kill -9 %d",gSystem->GetPid());
    std::cout<<name<<std::endl;
    gSystem->Exec(name);
}

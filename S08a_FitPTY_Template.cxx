#include"TemplateFitting.h"


#define TEST // ONLY do few bins, for debugging

std::vector<TCanvas*>        m_can_vec;
TemplateFitting::Fitting* fitresult;

/*-----------------------------------------------------------------------------
 *  Does the template fits and stores the fits as well as the vnn
 *-----------------------------------------------------------------------------*/
void S08a_FitPTY_Template(int Trig =0) {
    bool no_ZYAM = false;
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
        std::cout << "Working on XOR_E2 trigger!" << std::endl;
        if(Bins::same_side1) base = Form("%s/sameSide/xorE2",base.c_str());
        else{
            base = Form("%s/xorE2",base.c_str());
        }
    }
    else{
         std::cerr << "Error: no valid trigger index provided" << std::endl;
         exit(-1);
    }
    std::cout << base << endl;

    char name [600];
    char name1[600];
    sprintf(name, "%s/PTY1D.root", base.c_str());
    sprintf(name1 , "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/%.1fsigma/PTY1D.root",Bins::sigma1);
    if(Bins::same_side1) sprintf(name1 , "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/%.1fsigma/sameSide/PTY1D.root",Bins::sigma1);
    //if(Trig ==2) sprintf(name1 , "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/xorE2/PTY1D.root"); // if work on xor change peri bin
    //sprintf(name1, "%s/PTY1D.root", base.c_str()); // don't do ZYAM substruction
    //sprintf(name1, "%s/ZYAM1D.root", base.c_str());

    

    TFile *InFileCentral    = new TFile(name , "read");
    TFile *InFilePeripheral = new TFile(name1, "read");

    sprintf(name , "%s/TemplateFits.root"     , base.c_str());
    sprintf(name1, "%s/TemplateFits_vnn.root" , base.c_str());
    
    TFile *OutFile1 = new TFile(name , "recreate"); 
    TFile *OutFile2 = new TFile(name1, "recreate");
    OutFile1->cd();
    std::cout << name << "  " << name1 << std::endl;




#ifdef TEST
    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> trk_bins = {13};
    const std::vector<int> pt1_bins  = {Bins::GetPtaIndex(0.5, 5.0)};
    const std::vector<int> pt2_bins  = {Bins::GetPtbIndex(0.5, 5.0)};
    const std::vector<int> ch_bins   = {2};
    const std::vector<int> deta_bins = {Bins::GetDetaIndex(2.0, 5.0)};
    const std::vector<int> centbins_peripheral = {0,21};
    const std::vector<int> trkbins_peripheral =  {0};
#else
    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> trk_bins = Bins::TrkBins();
    const std::vector<int> pt1_bins  = Bins::PtaBins ();
    // const std::vector<int> pt1_bins  = {Bins::GetPtaIndex(0.5, 5.0)};
    const std::vector<int> pt2_bins  = Bins::PtbBins ();
    const std::vector<int> ch_bins   = Bins::ChBins  ();
    const std::vector<int> deta_bins = Bins::DetaBins();
    //const std::vector<int> centbins_peripheral = {52, 55, 62, 63};
    const std::vector<int> centbins_peripheral = Bins::CentBinsPeriph();
    const std::vector<int> trkbins_peripheral = Bins::TrkBinsPeriph();
#endif

    TH1D* h_v22   = new TH1D("h_v22", ";Centbin;v22"    , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_v33   = new TH1D("h_v33", ";Centbin;v33"    , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_v44   = new TH1D("h_v44", ";Centbin;v44"    , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_v55   = new TH1D("h_v55"  , ";Centbin;v55"  , Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);
    TH1D* h_scale = new TH1D("h_scale", ";Centbin;Scale", Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);

    //storing the original central and peripheral histograms
    TH1D* h_central;
    TH1D* h_peripheral;
    //histograms for fluctations
    TH1D *h_v22_fluc; //for storing v22 values from fit each after each fluctuation
    TH1D *h_fluctuate_peri;
    TH1D *h_fluctuate_central;

    //------------------------------------------------------------------------------
    for (int ipt1 : pt1_bins) {
        for (int ipt2 : pt2_bins) {
            for (int ich : ch_bins) {
                for (int ideta : deta_bins) {
                    std::cout << "Vnn " << ipt1 << "  " << ipt2 << "  " << ich << "  " << ideta << std::endl;
                    for (auto itrk2 : trkbins_peripheral) {
                        for (auto icent2 : centbins_peripheral) {
                            for (int itrk1 : trk_bins) {
                                h_v22->Reset();
                                h_v33->Reset();
                                h_v44->Reset();
                                h_v55->Reset();
                                for (int icent1 : cent_bins) {
                                    char hCentname[100], hPeriname[100];
                                    sprintf(hCentname, "PTY_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent1, itrk1, ipt1, ipt2, ich, ideta);
                                    sprintf(hPeriname, "PTY_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, itrk2, ipt1, ipt2, ich, ideta);
                                    h_central      = (TH1D*)InFileCentral   ->Get(hCentname)->Clone(Common::UniqueName().c_str());
                                    h_peripheral   = (TH1D*)InFilePeripheral->Get(hPeriname)->Clone(Common::UniqueName().c_str());
                                    if (!h_central   ) {std::cout << hCentname << " Not Found" << std::endl; throw std::exception();}
                                    if (!h_peripheral) {std::cout << hPeriname << " Not Found" << std::endl; throw std::exception();}
                                    h_v22_fluc = new TH1D(Form("h_v22_fluc_%s_pericent%.2i_peritrk%.2i",hCentname,icent2,itrk2), ";v22;counts",1000,-1,1);
                                    for(int ifluc =0; ifluc <1; ifluc++){ //no need to do toy mc 
                                        TRandom3 randGen(0); // Seed 0 for random initialization
                                        h_fluctuate_peri = (TH1D*)h_peripheral->Clone("h_peri_fluctuated");
                                        h_fluctuate_central = (TH1D*)h_central->Clone("h_peri_fluctuated");
                                        h_fluctuate_peri->Reset();
                                        h_fluctuate_central->Reset();
                                        for (int i = 1; i <= h_peripheral->GetNbinsX(); ++i) {
                                            double binContent = h_peripheral->GetBinContent(i);
                                            double binError = h_peripheral->GetBinError(i);
                                            double binContent_cent = h_central->GetBinContent(i);
                                            double binError_cent = h_central->GetBinError(i);
                                            // Generate fluctuated value
                                            double fluctuatedContent = randGen.Gaus(binContent, binError);
                                            double fluctuatedContent_cent = randGen.Gaus(binContent_cent, binError_cent);
                                            // Set the new bin content
                                            h_fluctuate_peri->SetBinContent(i, fluctuatedContent);
                                            h_fluctuate_central->SetBinContent(i, fluctuatedContent_cent);
                                            h_fluctuate_peri->SetBinError(i, binError); // Keep the same error
                                            h_fluctuate_central->SetBinError(i, binError_cent);

                                        }
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

                                        if(ifluc == 0) {
                                            fitresult = TemplateFitting::TemplateFit(h_central, h_peripheral); //first iteration on original histograms
                                            fitresult->Write(); //Save canvas of template fit only for the original histograms
                                        } 
                                        else{
                                            fitresult = TemplateFitting::TemplateFit(h_fluctuate_central, h_fluctuate_peri);
                                            }
                                        double v22, v22_err;
                                        fitresult->GetVnnAndError(v22, v22_err, 0);
                                        h_v22_fluc->Fill(v22);
                                    }
                                    sprintf(name, "cent%.2d_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent1,icent2, itrk2,itrk1, ipt1, ipt2, ich, ideta);
                                    std::string new_name = name;
                                    fitresult->SetName(new_name);
                                    double v22, v22_err;
                                    // v22 = h_v22_fluc->GetMean();
                                    // v22_err = h_v22_fluc->GetStdDev();
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

                                    if(icent1 == cent_bins.size() -1){
                                        fitresult = TemplateFitting::TemplateFit(h_central, h_peripheral); //from some reason reulsts is not saved for the last iteration
                                        fitresult->Write();                                                 // of cent bins so save it manually . TODO fix it
                                    }
                                    delete h_central;
                                    delete h_peripheral;
                                    // OutFile2->cd();
                                    // h_v22_fluc->Write();
                                }
                                OutFile2->cd();
                                sprintf(name, "h_v22_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, itrk2,itrk1,ipt1, ipt2, ich, ideta);
                                //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v22_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                                h_v22->SetName(name);
                                h_v22->Write();
                                sprintf(name, "h_v33_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, itrk2, itrk1, ipt1, ipt2, ich, ideta);
                                //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v33_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                                h_v33->SetName(name);
                                h_v33->Write();
                                sprintf(name, "h_v44_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, itrk2,itrk1,ipt1, ipt2, ich, ideta);
                                //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v44_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                                h_v44->SetName(name);
                                h_v44->Write();
                                sprintf(name, "h_v55_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, itrk2,itrk1,ipt1, ipt2, ich, ideta);
                                //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v55_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                                h_v55->SetName(name);
                                h_v55->Write();
                                sprintf(name, "h_v00_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, itrk2,itrk1,ipt1, ipt2, ich, ideta);
                                //if (l_use_peripheral_pp == 1) {sprintf(name, "h_v00_PPperiph%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent2, ipt1, ipt2, ich, ideta);}
                                h_scale->SetName(name);
                                h_scale->Write();

                                h_v22_fluc->Reset();
                            }
                        }
                    }
                }
            }
        }
    }
    InFileCentral   ->Close();
    InFilePeripheral->Close();
    OutFile1->Close();
    OutFile2->Close();
    


    // /* 
    //     //obtain the vn from the vnn
    //     //The vn is vn(cent,pta) obtained by dividing vnn(cent,pta,ptb) by sqrt(vn(cent,ptb,ptb))
    //     sprintf(name1,"%s/TemplateFits_vnn.root" ,base.c_str());
        
    //     OutFile2 = new TFile(name1,"update");
    //     OutFile2->cd();

    //     for(int ihar:{2,3,4}){
    //     for(int ich:ch_bins){
    //       for(int ideta:deta_bins){
    //         for(auto icent2:centbins_peripheral){
    //           for(int ipt2:pt2_bins){

    //             int ipt1_=Bins::GetPtaIndexForPtbIndex(ipt2);//code exits if pt1 bin corresponding to pt2 bin is not found
    //             sprintf(name,"h_v%d%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",ihar,ihar,icent2,ipt1_,ipt2,ich,ideta);
    //             sprintf(name1,"h_v%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d" ,ihar     ,icent2,ipt1_,ipt2,ich,ideta);
    //             TH1* h_vnn_   =(TH1*)(OutFile2->Get(name));
    //             if(!h_vnn_) {std::cout<<name<<std::endl;throw std::exception();}
    //             TH1* h_vn_diag=(TH1*)(h_vnn_)->Clone(name1);
    //             Common::Take_Sqrt(h_vn_diag);
    //             h_vn_diag->Write();

    //             for(int ipt1:pt1_bins){
    //               std::cout<<"Vn "<<ipt1<<"  "<<ipt2<<"  "<<ich<<"  "<<ideta<<std::endl;
    //               if(ipt1==ipt1_) continue;
    //               sprintf(name,"h_v%d%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",ihar,ihar,icent2,ipt1,ipt2,ich,ideta);
    //               sprintf(name1,"h_v%d_pericent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d" ,ihar     ,icent2,ipt1,ipt2,ich,ideta);
    //               TH1* h_vn=(TH1*)(OutFile2->Get(name))->Clone(name1);
    //               h_vn->Divide(h_vn_diag);
    //               h_vn->Write();
    //             }
    //           }
    //         }
    //       }
    //     }
    //     }
    //    OutFile2->Close();
    // */
    std::cout << "finished" << std::endl;
    sprintf(name,"kill -9 %d",gSystem->GetPid());
    std::cout<<name<<std::endl;
    gSystem->Exec(name);
}

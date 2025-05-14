//plot template fits for analysis section

#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"
using namespace Bins;
//using namespace DataSetEnums;


std::vector<TCanvas*> m_can_vec;
std::vector<TCanvas*> m_can_vec_cuts; // for same multiplicity and pT cuts
std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide/xor";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/template";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/template/xor";

void plot_ana_templateFit() {
    //SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    char name [600];
    char name1[600];

    TFile *InFile  = new TFile(Form("%s/TemplateFits.root",base.c_str()) , "read");

    // TFile * fileinput = new TFile("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/histograms.root");
    // TH2D *h2 = (TH2D*)fileinput->Get("hNtrkEff")->Clone();
    // TH1D *h1 = h2->ProfileX("h1");

    //------------------------------------------------------------------------------
    TCanvas *c1 = NULL;
    TH1D *h_central, *h_rescaledperipheral, *h_fit_func, *h_pars; 
    TF1  *f_pedestal, *f_vnn_combined;

    // std::vector<int> centbins_peripheral = Bins::CentBinsPeriph();
    // std::vector<int> trkbins_peripheral = Bins::TrkBinsPeriph();
    std::vector<int> centbins_peripheral = {22};
    std::vector<int> trkbins_peripheral ={5};
    std::vector<int> pt2_bins = Bins::PtbBins ();
    //std::vector<int> trk_bins = Bins::TrkBins ();
    std::vector<int> trk_bins = {4};
    std::vector<int> cent_bins = Bins::CentBins();
    std::vector<int> deta_bins = Bins::DetaBins();

    for (int itrk1 : trk_bins){
        for (int ipt1 : {Bins::GetPtaIndex(0.5, 5.0)}) {
            for (int ipt2 : {Bins::GetPtaIndex(0.5, 5.0)}) {
                for (int ich = 2; ich < 3; ich++) {
                    for (int ideta : deta_bins) {
                        std::cout << ipt1 << "  " << ipt2 << "  " << ich << "  " << ideta << std::endl;
                        for (auto icent2 : centbins_peripheral) {
                            for (auto itrk2 : trkbins_peripheral) {
                                  for (auto icent1 : cent_bins) {
                                //for (int icent1 = NCENT; icent1 < NCENT+ NCENT_ADD; icent1++) {
                                    std::cout << ipt1 << "  " << ipt2 << "  " << ich << "  " << ideta << "  " << icent2 << "  " << icent1 << std::endl;

                                    sprintf(name , "cent%.2d_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent1, icent2,itrk2,itrk1, ipt1, ipt2, ich, ideta);
                                    sprintf(name1, "h_central_%s"           , name); h_central           = (TH1D*)InFile->Get(name1); Common::CheckObject2(h_central           , name1, InFile);
                                    sprintf(name1, "h_rescaledperipheral_%s", name); h_rescaledperipheral= (TH1D*)InFile->Get(name1); Common::CheckObject2(h_rescaledperipheral, name1, InFile);
                                    sprintf(name1, "h_fit_func_%s"          , name); h_fit_func          = (TH1D*)InFile->Get(name1); Common::CheckObject2(h_fit_func          , name1, InFile);
                                    sprintf(name1, "h_pars_%s"              , name); h_pars              = (TH1D*)InFile->Get(name1); Common::CheckObject2(h_pars              , name1, InFile);
                                    sprintf(name1, "f_pedestal_%s"          , name); f_pedestal          = (TF1* )InFile->Get(name1); Common::CheckObject2(f_pedestal          , name1, InFile);
                                    sprintf(name1, "f_vnn_combined_%s"      , name); f_vnn_combined      = (TF1* )InFile->Get(name1); Common::CheckObject2(f_vnn_combined      , name1, InFile);

                                    if (h_central->Integral() == 0) {
                                        std::cout << "Skipping " << name << "as integral is 0" << std::endl;
                                        continue;
                                    }
                                    sprintf(name1, "can_fits_%s"            , name); c1                  = new TCanvas(name1, "", 800, 600);

                                    c1->cd();
                                    c1->SetLeftMargin (0.15);
                                    c1->SetTopMargin  (0.02);
                                    c1->SetBottomMargin(0.12);
                                    c1->SetRightMargin(0.02);

                                    float max = h_central->GetMaximum();
                                    float min = h_central->GetMinimum();
                                    h_central->SetMaximum(max + (max - min) / 6);
                                    h_central->SetMinimum(min - (max - min) / 6);


                                    h_central           ->Draw();
                                    h_rescaledperipheral->Draw("same");
                                    h_fit_func          ->Draw("same");
                                    f_pedestal          ->Draw("same");
                                    f_vnn_combined      ->Draw("same");

                                    TLegend *leg = new TLegend(0.18, 0.36, 0.50, 0.72);
                                    leg->SetTextSize(0.04);
                                    leg->SetBorderSize(0);
                                    leg->SetFillStyle(0);
                                    std::string label_central   = "Y(#Delta#phi)  "; //+label_cent(icent1);
                                  
                                        //std::string label_peripheral = "FY^{periph}(#Delta#phi) + G"; //+label_cent(icent2);
                                        std::string label_peripheral = "FY^{periph}(#Delta#phi)"; //+label_cent(icent2);
                                        leg->AddEntry(h_central           , label_central.c_str()    , "p");
                                        leg->AddEntry(h_rescaledperipheral, label_peripheral.c_str() , "p");
                                        leg->AddEntry(h_fit_func          , "Y^{templ}(#Delta#phi)"  , "l");
                                        //leg->AddEntry(f_vnn_combined      , "Y^{ridge}(#Delta#phi) +FY^{periph}(0)", "l");
                                        leg->AddEntry(f_vnn_combined      , "Y^{ridge}(#Delta#phi)", "l");
                                        //leg->AddEntry(f_pedestal          , "G + FY^{periph}(0)"       , "l");
                                        leg->AddEntry(f_pedestal          , "G "       , "l");
                                   
                                    leg->Draw();


                                    float X = 0.20, Y = 0.88, SIZE = 27;

                                    c1->cd(); X = 0.185; Y = 0.91; SIZE = 27;
                                    Common::myText2(X       , Y, 1, "ATLAS"         , SIZE, 73);
                                    Common::myText2(X + 0.12, Y, 1, Common::Internal, SIZE, 43); Y = Y - 0.06;
                                    Common::myText2(X       , Y, 1, "#it{pp}#sqrt{#it{s}} = 13.6 TeV", SIZE, 43); Y = Y - 0.06;

                                    if (Bins::GetPtbIndexForPtaIndex(ipt1) != ipt2) {
                                        Common::myText2(X, Y, 1, label_pta(ipt1) + ",  " + label_ptb(ipt2), SIZE, 43); Y = Y - 0.06;
                                    }
                                    else {
                                        Common::myText2(X, Y, 1, label_ptab(ipt1, ipt2), SIZE, 43); Y = Y - 0.06;
                                    }
                                    //double nrec = h1->GetBinContent(icent1 +1);
                                    X = 0.55 + 0.02;
                                    Y = 0.6;
                                    Common::myText2(X       , Y     , 1, label_cent     (icent1), SIZE, 43); Y = Y - 0.08;
                                    double chi2 = h_pars->GetBinContent(h_pars->GetNbinsX());
                                    int dof = h_central->GetNbinsX() - (h_pars->GetNbinsX()-1);
                                    chi2 = chi2/dof;
                                    Common::myText2(X       , Y     , 1, Form("#frac{#chi^{2}}{DOF} = %.2f",chi2), SIZE, 43); Y = Y - 0.06;
                                    //Common::myText2(X       , Y     , 1, Form("<N_{ch}^{rec}> = %.2f",nrec),SIZE, 43); Y = Y - 0.06;
                                    //Common::myText2(X       , Y     , 1, label_trk     (itrk1), SIZE, 43); Y = Y - 0.06;
                                    //Common::myText2(X       ,Y     ,1,label_cent_peri(icent2),SIZE,43);

                                    c1->cd();
                                    double pedestal = f_vnn_combined->GetParameter(0);
                                    X = 0.68; Y = 0.60; SIZE = 22;
                                    int nhars = f_vnn_combined->GetNpar() - 1;
                                    for (int ihar = 0; ihar < nhars; ihar++) {
                                        double vnn     = f_vnn_combined->GetParameter(ihar + 1);
                                        double vnn_err = f_vnn_combined->GetParError (ihar + 1);
                                        //sprintf(name,"Cv_{%d,%d}=%0.2e"   ,ihar+2,ihar+2,vnn    *pedestal); Common::myText2(X ,Y     ,1,name,SIZE,43);
                                        //sprintf(name,"       #pm %0.2e"                 ,vnn_err*pedestal); Common::myText2(X ,Y-0.05,1,name,SIZE,43);
                                        //Y = Y - 0.1;
                                    }
                                    gStyle->SetOptStat(0);
                                    m_can_vec.push_back(c1);
                                    c1->SaveAs(Form("%s/%s.pdf",figures.c_str(),name));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // std::cout << "Saving" << std::endl;
    // Common::SaveCanvas(m_can_vec, figures);


    InFile->Close();
    //OutFile->Close();
}

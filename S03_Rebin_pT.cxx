#include "bins.h"


TH2D *fg[Bins::NCENT + Bins::NCENT_ADD][Bins::NTRK + Bins::NTRK_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH];
TH2D *bg[Bins::NCENT + Bins::NCENT_ADD][Bins::NTRK + Bins::NTRK_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH];


/*-----------------------------------------------------------------------------
 *  Adds more pTa and pTb bins
 *-----------------------------------------------------------------------------*/
void S03_Rebin_pT(bool minbias =0) {
    std::string base;
    if(minbias){
        std::cout << "Working on Minbias triggers!" << std::endl;
         base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/minbias";
    }
    else{
        std::cout << "Working on ZDC triggers!" << std::endl;
        base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    }
    char name [600];
    char name1[600];
    sprintf(name , "%s/RebinTrk.root", base.c_str());
    sprintf(name1, "%s/Rebin_pT.root"       , base.c_str());
    TFile *input = new TFile(name );
    TFile *output = new TFile(name1, "recreate");
    output->cd();
    char fgname[100], bgname[100];

    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> trk_bins = Bins::TrkBins();
    const std::vector<int> pt1_bins  = Bins::PtaBins ();
    const std::vector<int> pt2_bins  = Bins::PtbBins ();


    //Read in All Histograms
    for (int icent : cent_bins) {
        for (int itrk : trk_bins) {
            for (int ipt1 = 0; ipt1 < Bins::NPT1; ipt1++) {
                std::cout << "Reading " << icent << "  " << ipt1 << std::endl;
                for (int ipt2 = 0; ipt2 < Bins::NPT2; ipt2++) {
                    for (int ich = 0; ich < Bins::NCH; ich++) {
                        sprintf(fgname, "fg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent, itrk, ipt1, ipt2, ich);
                        sprintf(bgname, "bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent, itrk, ipt1, ipt2, ich);
                        fg[icent][itrk][ipt1][ipt2][ich] = (TH2D*)input->Get(fgname);
                        bg[icent][itrk][ipt1][ipt2][ich] = (TH2D*)input->Get(bgname);

                        Common::CheckObject2(fg[icent][itrk][ipt1][ipt2][ich], fgname);
                        Common::CheckObject2(bg[icent][itrk][ipt1][ipt2][ich], bgname);
                    }
                }
            }
        }
    }

    //Add more pTa bins
    for (int icent : cent_bins) {
        for (int itrk : trk_bins) {
            for (int ipt1 = Bins::NPT1; ipt1 < Bins::NPT1 + Bins::NPT1_ADD; ipt1++) {
                std::cout << "Adding pta bins " << icent << "  " << itrk <<"  " <<ipt1 << std::endl;
                for (int ipt2 = 0; ipt2 < Bins::NPT2; ipt2++) {
                    for (int ich = 0; ich < Bins::NCH; ich++) {
                        for (int ipt1_add = Bins::pt1_add_lo[ipt1 - Bins::NPT1]; ipt1_add < Bins::pt1_add_up[ipt1 - Bins::NPT1]; ipt1_add++) {
                            TH2D* fg_add = fg[icent][itrk][ipt1_add][ipt2][ich];
                            TH2D* bg_add = bg[icent][itrk][ipt1_add][ipt2][ich];

                            if (ipt1_add == Bins::pt1_add_lo[ipt1 - Bins::NPT1]) {
                                sprintf(fgname, "fg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent, itrk, ipt1, ipt2, ich);
                                sprintf(bgname, "bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent, itrk, ipt1, ipt2, ich);
                                fg[icent][itrk][ipt1][ipt2][ich] = (TH2D*)fg_add->Clone(fgname);
                                bg[icent][itrk][ipt1][ipt2][ich] = (TH2D*)bg_add->Clone(bgname);
                                fg[icent][itrk][ipt1][ipt2][ich]->SetTitle(fgname);
                                bg[icent][itrk][ipt1][ipt2][ich]->SetTitle(bgname);
                            }
                            else {
                                fg[icent][itrk][ipt1][ipt2][ich]->Add(fg_add);
                                bg[icent][itrk][ipt1][ipt2][ich]->Add(bg_add);
                            }
                        }
                    }
                }
            }
        }
    }

    //Add more pTb bins
    for (int icent : cent_bins) {
        for (int itrk : trk_bins) {
            for (int ipt1 = 0; ipt1 < Bins::NPT1 + Bins::NPT1_ADD; ipt1++) {
                std::cout << "Adding ptb bins " << icent << "  " << itrk <<"  " <<ipt1 << std::endl;
                for (int ipt2 = Bins::NPT2; ipt2 < Bins::NPT2 + Bins::NPT2_ADD; ipt2++) {
                    for (int ich = 0; ich < Bins::NCH; ich++) {
                        for (int ipt2_add = Bins::pt2_add_lo[ipt2 - Bins::NPT2]; ipt2_add < Bins::pt2_add_up[ipt2 - Bins::NPT2]; ipt2_add++) {
                            TH2D* fg_add = fg[icent][itrk][ipt1][ipt2_add][ich];
                            TH2D* bg_add = bg[icent][itrk][ipt1][ipt2_add][ich];

                            if (ipt2_add == Bins::pt2_add_lo[ipt2 - Bins::NPT2]) {
                                sprintf(fgname, "fg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent,itrk, ipt1, ipt2, ich);
                                sprintf(bgname, "bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent,itrk, ipt1, ipt2, ich);
                                fg[icent][itrk][ipt1][ipt2][ich] = (TH2D*)fg_add->Clone(fgname);
                                bg[icent][itrk][ipt1][ipt2][ich] = (TH2D*)bg_add->Clone(bgname);
                                fg[icent][itrk][ipt1][ipt2][ich]->SetTitle(fgname);
                                bg[icent][itrk][ipt1][ipt2][ich]->SetTitle(bgname);
                            }
                            else {
                                fg[icent][itrk][ipt1][ipt2][ich]->Add(fg_add);
                                bg[icent][itrk][ipt1][ipt2][ich]->Add(bg_add);
                            }
                        }
                    }
                }
            }
        }
    }


    output->cd();
    //Write All Histograms
    for (int icent : cent_bins) {
        for (int itrk : trk_bins) {
            for (int ipt1 : pt1_bins) {
                std::cout << "Writing " << icent << "  " << itrk <<"  " <<ipt1 << std::endl;
                for (int ipt2 : pt2_bins) {
                    for (int ich = 0; ich < Bins::NCH; ich++) {
                        fg[icent][itrk][ipt1][ipt2][ich]->Write();
                        bg[icent][itrk][ipt1][ipt2][ich]->Write();
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

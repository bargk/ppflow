#include "bins.h"


TH2D *fg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH + Bins::NCH_ADD];
TH2D *bg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH + Bins::NCH_ADD];

/*-----------------------------------------------------------------------------
 *  Adds together the same-charge and opposite-charge histograms
 *  to get the combined charge histograms
 *-----------------------------------------------------------------------------*/
void S04_RebinCharge(int m_use_multiplicity = 0) {
    string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    if(m_use_multiplicity == 1) base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/multiplicity";
    char name [600];
    char name1[600];
    sprintf(name , "%s/Rebin_pT.root"   , base.c_str());
    sprintf(name1, "%s/RebinCharge.root", base.c_str());
    TFile *input = new TFile(name);
    TFile *output = new TFile(name1, "recreate");


    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> pt1_bins = Bins::PtaBins ();
    const std::vector<int> pt2_bins = Bins::PtbBins ();
    const std::vector<int> ch_bins  = Bins::ChBins  ();



    char fgname[100], bgname[100];
    for (int icent : cent_bins) {
        for (int pt1 : pt1_bins) {
            std::cout << "Reading " << icent << "  " << pt1 << std::endl;
            for (int pt2 : pt2_bins) {
                for (int ich = 0; ich < Bins::NCH; ich++) {
                    sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, pt1, pt2, ich);
                    sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, pt1, pt2, ich);
                    fg[icent][pt1][pt2][ich] = (TH2D*)input->Get(fgname);
                    bg[icent][pt1][pt2][ich] = (TH2D*)input->Get(bgname);
                }
            }
        }
    }

 for (int icent : cent_bins) {
    for (int pt1 : pt1_bins) {
        std::cout << "Rebinning " << icent << "  " << pt1 << std::endl;

        for (int pt2 : pt2_bins) {
            for (int ich = Bins::NCH; ich < Bins::NCH + Bins::NCH_ADD; ich++) {
                // Generate names for foreground and background histograms
                char fgname[256], bgname[256];
                sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, pt1, pt2, ich);
                sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, pt1, pt2, ich);

                // Rebinning process
                for (int ich_add = Bins::ch_add_lo[ich - Bins::NCH]; ich_add < Bins::ch_add_up[ich - Bins::NCH]; ich_add++) {
                    if (ich_add == Bins::ch_add_lo[ich - Bins::NCH]) {
                        // Clone histograms for the first additional bin
                        fg[icent][pt1][pt2][ich] = (TH2D*)fg[icent][pt1][pt2][ich_add]->Clone(fgname);
                        bg[icent][pt1][pt2][ich] = (TH2D*)bg[icent][pt1][pt2][ich_add]->Clone(bgname);
                        fg[icent][pt1][pt2][ich]->SetTitle(fgname);
                        bg[icent][pt1][pt2][ich]->SetTitle(bgname);
                    } else {
                        // Add histograms for subsequent additional bins
                        fg[icent][pt1][pt2][ich]->Add(fg[icent][pt1][pt2][ich_add]);
                        bg[icent][pt1][pt2][ich]->Add(bg[icent][pt1][pt2][ich_add]);
                    }
                }
            }
        }
    }
}


    output->cd();
    for (int icent : cent_bins) {
        for (int pt1 : pt1_bins) {
            std::cout << "Writing " << icent << "  " << pt1 << std::endl;
            for (int pt2 : pt2_bins) {
                for (int ich : ch_bins) {
                    fg[icent][pt1][pt2][ich]->Write();
                    bg[icent][pt1][pt2][ich]->Write();
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


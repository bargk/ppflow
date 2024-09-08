#include "bins.h"
#include "Defs.h"

TH2D *fg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH + Bins::NCH_ADD];
TH2D *bg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH + Bins::NCH_ADD];


/*-----------------------------------------------------------------------------
 *  Makes 1D Pair distributions for several Deta ranges from the 2D pair distributions
 *-----------------------------------------------------------------------------*/
void S04_ProjectionX() {
    string base = directory;
    char name [600];
    char name1[600];

    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> pt1_bins = Bins::PtaBins ();
    const std::vector<int> pt2_bins = Bins::PtbBins ();
    const std::vector<int> ch_bins  = Bins::ChBins  ();
    const std::vector<int> deta_bins = Bins::DetaBins();


    sprintf(name , "%s/RebinCharge.root", base.c_str());
    sprintf(name1, "%s/ProjectionX.root", base.c_str());
    TFile *input = new TFile(name );
    TFile *output = new TFile(name1, "recreate");


    //Read In
    std::cout << "Reading Starting" << std::endl;
    char fgname[100], bgname[100], PjXfgname[100], PjXbgname[100], PjXconame[100];
    input->ReadAll();
    TIter next(input->GetList());
    TObject *obj;
    for (int icent : cent_bins) {
        for (int ipt1 : pt1_bins) {
            for (int ipt2 : pt2_bins) {
                std::cout << "Reading " << icent << "  " << ipt1 << "  " << ipt2 << std::endl;
                for (int ich : ch_bins) {
                    sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);
                    sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);

                    obj = next();
                    fg[icent][ipt1][ipt2][ich] = (TH2D*)obj;
                    if (!Common::CheckObject(obj, fgname)) throw std::exception();

                    obj = next();
                    bg[icent][ipt1][ipt2][ich] = (TH2D*)obj;
                    if (!Common::CheckObject(obj, bgname)) throw std::exception();
                }
            }
        }
    }
    std::cout << "Reading Done" << std::endl;


    int bins_lo[Bins::NDETA] = {0};
    int bins_up[Bins::NDETA] = {0};
    TH2D* hist = (TH2D*)obj;
    for (int ieta : deta_bins) {
        bins_lo[ieta] = hist->GetYaxis()->FindBin(Bins::DETA_LO[ieta] + 0.001);
        bins_up[ieta] = hist->GetYaxis()->FindBin(Bins::DETA_HI[ieta] - 0.001);
    }

    //Rebin
    std::cout << "Projection Starting" << std::endl;
    int count = 1;
    for (int icent : cent_bins) {
        for (int ipt1 : pt1_bins) {
            for (int ipt2 : pt2_bins) {
                std::cout << icent << "  " << ipt1 << "  " << ipt2 << "  ::" << count << std::endl;
                for (int ich : ch_bins) {
                    for (int ieta : deta_bins) {
                        sprintf(PjXfgname, "PjX_fg_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent, ipt1, ipt2, ich, ieta);
                        sprintf(PjXbgname, "PjX_bg_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", icent, ipt1, ipt2, ich, ieta);
                        //sprintf(PjXconame,"PjX_co_cent%.2d_pta%d_ptb%.2d_ch%d_deta%.2d",icent,ipt1,ipt2,ich,ieta);
                        TH1D* PjX_fg  = (TH1D*)fg[icent][ipt1][ipt2][ich]->ProjectionX(PjXfgname, bins_lo[ieta], bins_up[ieta]);
                        TH1D* PjX_bg  = (TH1D*)bg[icent][ipt1][ipt2][ich]->ProjectionX(PjXbgname, bins_lo[ieta], bins_up[ieta]);
                        Common::Symmetrize_1D(PjX_fg);
                        Common::Symmetrize_1D(PjX_bg);

                        //PjX_fg ->SetTitle(PjXfgname);
                        //PjX_bg ->SetTitle(PjXbgname);
                        //TH1D* PjX_co  = (TH1D*)PjX_fg->Clone(PjXconame);
                        //PjX_co ->Divide(PjX_bg);
                        //PjX_co ->SetTitle(PjXconame);

                        count++;
                        //std::cout<<icent<<"  "<<ipt1<<"  "<<ipt2<<"  "<<ich<<"  "<<ieta<<"  :: "<<count<<std::endl;

                        PjX_fg->Write();
                        PjX_bg->Write();
                        //PjX_co->Write();
                        delete PjX_fg;
                        delete PjX_bg;
                        //delete PjX_co;
                    }
                }
            }
        }
    }
    std::cout << "Projection Done" << std::endl;
}

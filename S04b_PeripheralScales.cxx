#include "bins.h"
#include "common.C"
#include "Defs.h"

TH2D *fg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH + Bins::NCH_ADD];
TH2D *bg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD][Bins::NPT2 + Bins::NPT2_ADD][Bins::NCH + Bins::NCH_ADD];
double NTrigs[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1 + Bins::NPT1_ADD] = {{0}};


/*-----------------------------------------------------------------------------
 *  Determines the Peripheral scales to be used while performing the
 *  old-style peripheral subtraction.
 *  The peripheral scales are determined by matching the near-side peak
 *  in the 2D PTYs for the peripheral reference and signal centrality
 *-----------------------------------------------------------------------------*/
void S04b_PeripheralScales() {
    string base = directory;
    char name [600];
    char name1[600];
    char name2[600];

    const std::vector<int> cent_bins = Bins::CentBins();
    const std::vector<int> pt1_bins = Bins::PtaBins ();
    const std::vector<int> pt2_bins = Bins::PtbBins ();
    const std::vector<int> ch_bins  = Bins::ChBins  ();
    const std::vector<int> cent_periph = Bins::CentBinsPeriph();

    //Read In Ntrigs
    sprintf(name, "%s/RebinEff.root", base.c_str());
    TFile *Cent = new TFile(name, "read");
    for (int icent = 0; icent < Bins::NCENT + Bins::NCENT_ADD; icent++) {
        sprintf(name, "N_trigger_cent%.2d", icent);
        TH1D* h_NTrigs = (TH1D*)Cent->Get(name);
        if (!Common::CheckObject(h_NTrigs, name)) throw std::exception();
        for (int ipt1 = 0; ipt1 < Bins::NPT1 + Bins::NPT1_ADD; ipt1++) {
            int bin_lo  = h_NTrigs->FindBin(Bins::PT1_LO[ipt1] + .0001);
            int bin_high = h_NTrigs->FindBin(Bins::PT1_HI[ipt1] - .0001);
            std::cout << bin_lo << "  " << bin_high    << "  "
                      << h_NTrigs->GetBinLowEdge(bin_lo)    << "  "
                      << h_NTrigs->GetBinLowEdge(bin_high + 1) << "  "
                      << Bins::PT1_LO[ipt1]                 << "  "
                      << Bins::PT1_HI[ipt1] << std::endl;
            if (fabs(h_NTrigs->GetBinLowEdge(bin_lo) - Bins::PT1_LO[ipt1]) > 0.001) {
                cout << "1  " << h_NTrigs->GetBinLowEdge(bin_lo) << "  " << Bins::PT1_LO[ipt1] << endl;
                throw std::exception();
            }
            if (fabs(h_NTrigs->GetBinLowEdge(bin_high + 1) - Bins::PT1_HI[ipt1]) > 0.001) {
                cout << "2  " << h_NTrigs->GetBinLowEdge(bin_high + 1) << "  " << Bins::PT1_HI[ipt1] << endl;
                throw std::exception();
            }
            NTrigs[icent][ipt1] = h_NTrigs->Integral(bin_lo, bin_high);
        }
    }
    Cent->Close();



    //Read In 2D Corrs
    sprintf(name , "%s/RebinCharge.root"     , base.c_str());
    sprintf(name1, "%s/PeripheralScales.root", base.c_str());
    sprintf(name2, "%s/PTY2D.root"           , base.c_str());
    TFile *input   = new TFile(name );
    TFile *output  = new TFile(name1, "recreate");
    TFile *output2 = new TFile(name2, "recreate");
    char fgname[100], bgname[100];
    input->ReadAll();
    TIter next(input->GetList());
    TObject *obj;
    for (int icent : cent_bins) {
        for (int ipt1 : pt1_bins) {
            for (int ipt2 : pt2_bins) {
                for (int ich : ch_bins) {
                    cout << "Reading " << icent << "  " << ipt1 << "  " << ipt2 << endl;
                    sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);
                    sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);

                    obj = next();
                    fg[icent][ipt1][ipt2][ich] = (TH2D*)obj;
                    if (!Common::CheckObject(obj, fgname)) throw std::exception();

                    obj = next();
                    bg[icent][ipt1][ipt2][ich] = (TH2D*)obj;
                    if (!Common::CheckObject(obj, bgname)) throw std::exception();

                    double integral = bg[icent][ipt1][ipt2][ich]->Integral();
                    fg[icent][ipt1][ipt2][ich]->Divide(bg[icent][ipt1][ipt2][ich]);
                    fg[icent][ipt1][ipt2][ich]->Scale(integral / NTrigs[icent][ipt1]);
                }
            }
        }
    }




    /*-----------------------------------------------------------------------------
     * Evaluate Scales
     *-----------------------------------------------------------------------------*/
    int count = 1;
    for (int ipt1 : pt1_bins) {
        for (int ipt2 : pt2_bins) {
            cout << "Evaluating " << ipt1 << "  " << ipt2 << "  ::" << count << endl;
            for (int ich : ch_bins) {
                for (int icent_periph : cent_periph) {
                    sprintf(fgname, "h_scale_pericent%.2d_pta%d_ptb%.2d_ch%d", icent_periph, ipt1, ipt2, ich);
                    TH1D*hist = new TH1D(fgname, ";centbin;", Bins::NCENT + Bins::NCENT_ADD, 0, Bins::NCENT + Bins::NCENT_ADD);

                    fg[icent_periph][ipt1][ipt2][ich]->GetYaxis()->SetRangeUser(0 , 1);
                    fg[icent_periph][ipt1][ipt2][ich]->GetXaxis()->SetRangeUser(-1, 1);
                    double Unsub_jet_periph = fg[icent_periph][ipt1][ipt2][ich]->Integral();


                    fg[icent_periph][ipt1][ipt2][ich]->GetYaxis()->SetRangeUser(2 , 5);
                    fg[icent_periph][ipt1][ipt2][ich]->GetXaxis()->SetRangeUser(-1, 1);
                    double long_range_periph = fg[icent_periph][ipt1][ipt2][ich]->Integral() / 3.0;

                    double sub_jet_periph = (Unsub_jet_periph - long_range_periph);

                    for (int icent : cent_bins) {
                        fg[icent][ipt1][ipt2][ich]->GetYaxis()->SetRangeUser(0 , 1);
                        fg[icent][ipt1][ipt2][ich]->GetXaxis()->SetRangeUser(-1, 1);
                        double Unsub_jet = fg[icent][ipt1][ipt2][ich]->Integral();

                        fg[icent][ipt1][ipt2][ich]->GetYaxis()->SetRangeUser(2 , 5);
                        fg[icent][ipt1][ipt2][ich]->GetXaxis()->SetRangeUser(-1, 1);
                        double long_range = fg[icent][ipt1][ipt2][ich]->Integral() / 3.0;

                        double sub_jet = (Unsub_jet - long_range);

                        hist->SetBinContent(icent + 1, sub_jet / sub_jet_periph);
                        count++;

                        fg[icent][ipt1][ipt2][ich]->GetYaxis()->SetRangeUser(0 , 5);
                        fg[icent][ipt1][ipt2][ich]->GetXaxis()->SetRangeUser(-Common::PI / 2.0, 3 * Common::PI / 2.0);

                        output2->cd();
                        fg[icent][ipt1][ipt2][ich]->Write();
                    }
                    output->cd();
                    hist->Write();
                }
            }
        }
    }
}

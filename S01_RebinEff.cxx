#include "bins.h"
#include "Defs.h"

TH2D *fg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1][Bins::NPT2][Bins::NCH];
TH2D *bg[Bins::NCENT + Bins::NCENT_ADD][Bins::NPT1][Bins::NPT2][Bins::NCH];
TH1D *N_trigger[Bins::NCENT + Bins::NCENT_ADD];


/*-----------------------------------------------------------------------------
 *  Reads in original 2D histograms and adds more effective energy bins
 *------------------------------------------------------- ---------------------*/
void S01_RebinEff() {
    string base = directory;
    char name [600];
    char name1[600];
    TFile *input;
    TFile *output;
    sprintf(name , "%s/histograms.root"   , base.c_str());
    sprintf(name1, "%s/RebinEff.root", base.c_str());
    input  = new TFile(name);
    output = new TFile(name1, "recreate");
    output->cd();
    if (input->IsZombie()) {std::cout << name << " Not Found" << std::endl; throw std::exception();}


    char fgname[100], bgname[100], n_tname[100];

//Read in All Histograms
    cout << "Reading histograms" << endl;
    for (int icent = 0; icent < Bins::NCENT; icent++) {
        for (int ipt1 = 0; ipt1 < Bins::NPT1; ipt1++) {
            for (int ipt2 = 0; ipt2 < Bins::NPT2; ipt2++) {
                for (int ich = 0; ich < Bins::NCH; ich++) {
                    sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);
                    sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);
                    fg[icent][ipt1][ipt2][ich] = (TH2D*)input->Get(fgname);
                    bg[icent][ipt1][ipt2][ich] = (TH2D*)input->Get(bgname);
                    if (!fg[icent][ipt1][ipt2][ich]) {std::cout << fgname << " Not found" << std::endl; throw std::exception();}
                    if (!bg[icent][ipt1][ipt2][ich]) {std::cout << bgname << " Not found" << std::endl; throw std::exception();}
                }
            }
        }
    }

//Add more effective energy bins
cout << "Adding more effective energy bins" << endl;
    for (int icent = Bins::NCENT; icent < Bins::NCENT + Bins::NCENT_ADD; icent++) {
        for (int ipt1 = 0; ipt1 < Bins::NPT1; ipt1++) {
            for (int ipt2 = 0; ipt2 < Bins::NPT2; ipt2++) {
                for (int ich = 0; ich < Bins::NCH; ich++) {
                    for (int icent_add = Bins::cent_add_lo[icent - Bins::NCENT];
                            icent_add < Bins::cent_add_up[icent - Bins::NCENT];
                            icent_add++) {
                        TH2D* fg_add = fg[icent_add][ipt1][ipt2][ich];
                        TH2D* bg_add = bg[icent_add][ipt1][ipt2][ich];

                        if (icent_add == Bins::cent_add_lo[icent - Bins::NCENT]) {
                            sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);
                            sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, ipt1, ipt2, ich);
                            fg[icent][ipt1][ipt2][ich] = (TH2D*)fg_add->Clone(fgname);
                            bg[icent][ipt1][ipt2][ich] = (TH2D*)bg_add->Clone(bgname);
                            fg[icent][ipt1][ipt2][ich]->SetTitle(fgname);
                            bg[icent][ipt1][ipt2][ich]->SetTitle(bgname);
                        }
                        else {
                            fg[icent][ipt1][ipt2][ich]->Add(fg_add);
                            bg[icent][ipt1][ipt2][ich]->Add(bg_add);
                        }
                    }
                }
            }
        }
    }

    output->cd();
//Write All Histograms
cout << "Writing new histograms" << endl;
    for (int icent = 0; icent < Bins::NCENT + Bins::NCENT_ADD; icent++) {
        for (int ipt1 = 0; ipt1 < Bins::NPT1; ipt1++) {
            for (int ipt2 = 0; ipt2 < Bins::NPT2; ipt2++) {
                for (int ich = 0; ich < Bins::NCH; ich++) {
                    fg[icent][ipt1][ipt2][ich]->Write();
                    bg[icent][ipt1][ipt2][ich]->Write();
                }
            }
        }
    }



//Trigger particle spectra
    for (int icent = 0; icent < Bins::NCENT; icent++) {
        sprintf(n_tname, "N_trigger_cent%.2d", icent);
        N_trigger[icent] = (TH1D*)input->Get(n_tname);
        N_trigger[icent]->Write();
    }


    for (int icent = Bins::NCENT; icent < Bins::NCENT + Bins::NCENT_ADD; icent++) {
        for (int icent_add = Bins::cent_add_lo[icent - Bins::NCENT];
                icent_add < Bins::cent_add_up[icent - Bins::NCENT];
                icent_add++) {
            TH1D *n_t = N_trigger[icent_add];

            if (icent_add == Bins::cent_add_lo[icent - Bins::NCENT]) {
                sprintf(n_tname, "N_trigger_cent%.2d", icent);
                N_trigger[icent] = (TH1D*)n_t->Clone(n_tname);
                N_trigger[icent]->SetTitle(n_tname);
            }
            else N_trigger[icent]->Add(n_t);
        }
        N_trigger[icent]->Write();
    }

    output->Close();
    std::cout << "finished" << std::endl;
    sprintf(name,"kill -9 %d",gSystem->GetPid());
    std::cout<<name<<std::endl;
    gSystem->Exec(name);
}

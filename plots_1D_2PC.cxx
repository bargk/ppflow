#include "common.C"

char name  [600];
char name1 [600];
char fgname[600];
char bgname[600];

// This plots C(dphi) and not Y(dphi)
void plots_1D_2PC(){
    char directory[600] = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/measure_1D_2PC";
    char base[600] = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";
    gSystem->Exec(Form("mkdir -p %s",directory));
    //open projection
    TFile *input1D = TFile::Open(Form("%s/ProjectionX.root",base),"READ");

    int ieff_energy=0, inch=0, it=0; //TODO make it as function argument
    for(int ieff_energy =0; ieff_energy<1; ieff_energy++){
        for(int inch =0; inch<10; inch++){
            for(int it =0; it<2; it++){
                sprintf(fgname, "PjX_fg_%d_%d_%d",ieff_energy,inch,it);
                sprintf(bgname, "PjX_bg_%d_%d_%d",ieff_energy,inch,it);


                TH1D* fg_tmp = (TH1D*)(input1D->Get(fgname));
                TH1D* bg_tmp = (TH1D*)(input1D->Get(bgname));
                Common::CheckObject2(fg_tmp, fgname);
                Common::CheckObject2(bg_tmp, bgname);
                TH1D* fg = (TH1D*)(fg_tmp->Clone(Common::UniqueName().c_str()));
                TH1D* bg = (TH1D*)(bg_tmp->Clone(Common::UniqueName().c_str()));

                #ifndef NORM1 //normalize bg to same pairs as fg
                double scale = fg->Integral() / bg->Integral();
                bg->Scale(scale);
            #endif
                // if (m_rebin_X > 1) {
                //     fg->Rebin(m_rebin_X);
                //     bg->Rebin(m_rebin_X);
                // }
                fg->Divide(bg);
            #ifdef NORM1 //normalize to average of 1.0
                fg->Scale(fg->GetNbinsX() / fg->Integral());
            #endif

                fg->GetXaxis()->SetTitle("#Delta#phi");
                fg->GetYaxis()->SetTitle("C(#Delta#phi)");

                // if (fabs(m_range_up + 1000) > 0.001) fg->SetMaximum(m_range_up);
                // if (fabs(m_range_lo + 1000) > 0.001) fg->SetMinimum(m_range_lo);

                sprintf(name, "Can_1D_%d_%d_%d",ieff_energy,inch,it);
                TCanvas *Can = new TCanvas(name, name, 600, 450);
                //m_can_vec.push_back(Can);
                Can->SetLeftMargin (0.12);
                Can->SetTopMargin  (0.05);
                Can->SetRightMargin(0.05);
                fg->Draw();
                Can->SaveAs(Form("%s/%s.pdf",directory,name));

            }
        }
    }
    
}
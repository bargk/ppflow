
#include "common.C"
#define SYMMETRIZE
//#define RENORMALIZE

float m_range_up = -1000;
float m_range_lo = -1000;
void Draw2D(int icent,int itrk, int pt1, int pt2, int ich);

void plots_2D_2PC(){


    // for (int ieff_energy =0; ieff_energy<1; ieff_energy++){
    //     for (int inch=0; inch<10; inch++){
    //         for(int it=0; it<2; it++){
    //             Draw2D(ieff_energy, inch, it);
    //         }
    //     }
    // }

    //m_range_up = 1.4; m_range_lo = 0.9; Draw2D(1,5,39,2);
    // m_range_up = 1.35; m_range_lo = 0.9; Draw2D(0,0,1);
    //m_range_up = 1.2; m_range_lo = 0.9; Draw2D(0,1,0);
    // // m_range_up = 1.05; m_range_lo = 0.9; Draw2D(0,1,1);
    // // m_range_up = 1.05; m_range_lo = 0.9; Draw2D(0,2,0);
    // m_range_up = 1.1; m_range_lo = 0.94; Draw2D(0,4,0);
    // m_range_up = 1.1; m_range_lo = 0.93; Draw2D(0,9,0);
    Draw2D(05,13,5,05,2);
}
//TODO modify it to make rebin charge plots
void Draw2D(int icent, int pt1, int pt2, int ich){
    char directory[600] = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/measure_2D_2PC";
    char directory_root[600] = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles";

    gSystem->Exec(Form("mkdir -p %s",directory));
    int   m_rebin_X=2 , m_rebin_Y=2;

    char fgname[600];
    char bgname[600];
    TFile *input2D = TFile::Open(Form("%s/RebinCharge.root",directory_root),"READ");

    sprintf(fgname, "fg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, pt1, pt2, ich);
    sprintf(bgname, "bg_cent%.2d_pta%d_ptb%.2d_ch%d", icent, pt1, pt2, ich);

    TH2D* fg_temp = (TH2D*)(input2D->Get(fgname));
    TH2D* bg_temp = (TH2D*)(input2D->Get(bgname));
    Common::CheckObject2(fg_temp, fgname, input2D);
    Common::CheckObject2(bg_temp, bgname, input2D);
    TH2D* fg = (TH2D*)(fg_temp->Clone(Common::UniqueName().c_str()));
    TH2D* bg = (TH2D*)(bg_temp->Clone(Common::UniqueName().c_str()));
    
    bool ret = false;
    if (fg->Integral() == 0) {std::cout << fgname << " has 0 integral" << std::endl; ret = true;}
    if (bg->Integral() == 0) {std::cout << bgname << " has 0 integral" << std::endl; ret = true;}
    if (ret) return;

//the scale ensures that at large deta the corr has mean value=1
    TH1* fg1d = (TH1D*) fg->ProjectionX(Common::UniqueName().c_str(), 21, 50);
    TH1* bg1d = (TH1D*) bg->ProjectionX(Common::UniqueName().c_str(), 21, 50);
    double scale = bg1d->Integral() / fg1d->Integral();

#ifdef SYMMETRIZE
    fg = Common::Symmetrize_2D(fg);
    bg = Common::Symmetrize_2D(bg);
#endif
    // fg->Rebin2D(m_rebin_X, m_rebin_Y);
    // bg->Rebin2D(m_rebin_X, m_rebin_Y);


    #ifdef RENORMALIZE
    {
        int BinsDeta = fg->GetNbinsY();
        int BinsDphi = fg->GetNbinsX();
        for (int ibin_deta = 1; ibin_deta <= BinsDeta; ibin_deta++) {
            double count_fg = 0, count_bg = 0;
            for (int ibin_phi = 1; ibin_phi <= BinsDphi; ibin_phi++) {
                count_fg += fg->GetBinContent(ibin_phi, ibin_deta);
                count_bg += bg->GetBinContent(ibin_phi, ibin_deta);
            }
            double scale_ = count_fg / count_bg;
            for (int ibin_phi = 1; ibin_phi <= BinsDphi; ibin_phi++) {
                bg->SetBinContent(ibin_phi, ibin_deta, bg->GetBinContent(ibin_phi, ibin_deta)*scale_);
            }
        }
    }
    fg->Divide(bg);
#else
    fg->Divide(bg);
    fg->Scale(scale);
#endif


    // m_format["XTitleSize"  ] = 0.05;
    // m_format["YTitleSize"  ] = 0.05;
    // m_format["XTitleOffset"] = 1.2;
    // m_format["YTitleOffset"] = 1.2;
    // Common::FormatHist(fg, m_format);
    fg->GetXaxis()->SetTitle("#Delta#phi");
    fg->GetYaxis()->SetTitle("|#Delta#eta|");
    fg->GetZaxis()->SetTitle("C(|#Delta#eta|,#Delta#phi)");
    fg->GetZaxis()->SetNdivisions (505);
    fg->GetZaxis()->CenterTitle   ();
    fg->GetZaxis()->SetTitleSize  (0.05);
    fg->GetZaxis()->SetLabelSize  (0.05);
    fg->GetZaxis()->SetTitleOffset(1.5 );
    fg->SetLineColor(1);
    fg->SetLineWidth(2);
    Common::format_hist(fg);

    if (m_range_up != -1000) fg->SetMaximum(m_range_up);
    if (m_range_lo != -1000) fg->SetMinimum(m_range_lo);
#ifdef SYMMETRIZE
    fg->GetYaxis()->SetTitle("#Delta#eta");
    fg->GetZaxis()->SetTitle("C(#Delta#eta,#Delta#phi)");
    fg->GetYaxis()->SetRangeUser(-4.59, 4.59);
    // if (m_data_type == DataSetEnums::DATA_5TEV) {
    //     fg->GetYaxis()->SetRangeUser(-4, 4);
    // }
#else
    fg->GetYaxis()->SetRangeUser(0.0, 4.4);
#endif

    char name[600];
    sprintf(name, "Can_2Dcent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent,itrk, pt1, pt2, ich);

    TCanvas *Can = new TCanvas(name, name, 600, 600);
    Can->SetTheta(60);
    Can->SetPhi  (130);
    gStyle->SetPalette(kRainBow);
    fg->Draw("SURF1FB");
    //m_can_vec.push_back(Can);
    Can->SetLeftMargin (0.17);
    Can->SetTopMargin  (0.05);
    Can->SetRightMargin(0.03);

    Can->SaveAs(Form("%s/%s.pdf",base,name));

}
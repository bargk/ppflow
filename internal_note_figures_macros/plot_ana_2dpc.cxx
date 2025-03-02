//ploting 2d pc for analaysis section
#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

#define SYMMETRIZE  
std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/2dpc";
float m_range_lo = -1000;
float m_range_up = -1000;
void Draw2D(int icent,int itrk, int pt1, int pt2, int ich);
TFile *input0;
std::map<std::string, double> m_format = Common::StandardFormat();

void plot_ana_2dpc(){
    SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    int m_rebin_X =2; int m_rebin_Y =2;

    m_range_up = 1.1;  m_range_lo = 0.95; Draw2D(22,13, 5, 05, 2);
    m_range_up = 1.1;  m_range_lo = 0.95; Draw2D(19,13, 5, 05, 2); //those 2 is for the minbias range

    m_range_up = 1.2;  m_range_lo = 0.93; Draw2D(23,13, 5, 05, 2);
    m_range_up = 1.2;  m_range_lo = 0.93; Draw2D(24,13, 5, 05, 2);
    m_range_up = 1.15; m_range_lo = 0.93; Draw2D(25,13, 5, 05, 2);
    m_range_up = 1.15; m_range_lo = 0.93; Draw2D(26,13, 5, 05, 2);
    m_range_up = 1.15; m_range_lo = 0.93; Draw2D(27,13, 5, 05, 2);
    m_range_up = 1.15; m_range_lo = 0.93; Draw2D(28,13, 5, 05, 2);
    m_range_up = 1.15; m_range_lo = 0.93; Draw2D(29,13, 5, 05, 2);
    m_range_up = 1.13; m_range_lo = 0.93; Draw2D(30,13, 5, 05, 2);
}

void Draw2D(int icent,int itrk, int pt1, int pt2, int ich){
   
    int   m_rebin_X=2 , m_rebin_Y=2;

    char fgname[600];
    char bgname[600];
    TFile *input2D = TFile::Open(Form("%s/RebinCharge.root",base.c_str()),"READ");

    sprintf(fgname, "fg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent,itrk, pt1, pt2, ich);
    sprintf(bgname, "bg_cent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent,itrk, pt1, pt2, ich);

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
    fg->Rebin2D(m_rebin_X, m_rebin_Y);
    bg->Rebin2D(m_rebin_X, m_rebin_Y);


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


    m_format["XTitleSize"  ] = 0.05;
    m_format["YTitleSize"  ] = 0.05;
    m_format["XTitleOffset"] = 1.2;
    m_format["YTitleOffset"] = 1.2;
    Common::FormatHist(fg, m_format);
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
    sprintf(name, "Can_2Dcent%.2d_trk%.2d_pta%d_ptb%.2d_ch%d", icent, itrk, pt1, pt2, ich);

    TCanvas *Can = new TCanvas(name, name, 600, 600);
    Can->SetTheta(60);
    Can->SetPhi  (130);
    gStyle->SetPalette(kRainBow);
    fg->Draw("SURF1FB");
    //m_can_vec.push_back(Can);
    Can->SetLeftMargin (0.17);
    Can->SetTopMargin  (0.05);
    Can->SetRightMargin(0.03);

    float X = 0.05, Y = 0.95;
    int SIZE = 20;
    Common::myText2(X, Y, 1, "ATLAS ", SIZE, 73);
    Common::myText2(X + 0.15, Y, 1, Common::Internal, SIZE, 43); Y -= 0.05;
    Common::myText2(X, Y, 1, "#it{pp}, #sqrt{#it{s}} = 13.6 TeV", SIZE, 43); Y -= 0.05;

    X = 0.63; Y = 0.95;
    Common::myText2(X , Y , 1, Bins::label_cent(icent)     , SIZE, 43);
    Common::myText2(X + 0.06       , Y - 0.08      , 1, Bins::label_ptab(pt1, pt2), SIZE, 43);


    Can->SaveAs(Form("%s/%s.pdf",figures.c_str(),name));

}

#ifndef __COMMON_C__
#define __COMMON_C__


namespace Common{
    double PI= acos(-1.0);
    //std::string Internal = "";
    string Internal="Internal";


    void SaveCanvas(std::vector<TCanvas*> &can_vec, std::string base) {
    char name[600];
    gSystem->Exec(Form("mkdir -p %s",base.c_str()));
    for (auto &can : can_vec) {
        sprintf(name,"%s/%s.png",base.c_str(),can->GetName());
        can->SaveAs(name);
        // sprintf(name, "figs/%s_%s.pdf", can->GetName(), base.c_str());
        // can->SaveAs(name);
        //sprintf(name,"figs/%s_%s.png",can->GetName(),base.c_str());
        //can->SaveAs(name);
        //sprintf(name,"figs/%s_%s.C",can->GetName(),base.c_str());
        // can->SaveAs(name);
        can->Write();
    }
    can_vec.clear();
    }

    void myText(float x, float y, Color_t color, std::string text, float tsize = 0.06) {
    TLatex l; //l.SetTextAlign(12);
    l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x, y, text.c_str());
    }

    void myText2(float x, float y, Color_t color, std::string text, int size, int font) {
        TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
        l.SetNDC();
        l.SetTextColor(color);
        l.SetTextSize(size);
        l.SetTextFont(font);
        l.DrawLatex(x, y, text.c_str());
    }


    void ShiftXaxis(TH1* h, float shift, int excluded_lower_bins = 0) {
      double X[100];
      int NBins = h->GetNbinsX();
      if (NBins > 100) {std::cout << "void ShiftXaxis(TH1* h) Too many bins" << std::endl; throw std::exception();}

      for (int ibin = 0; ibin < excluded_lower_bins; ibin++) {
          X[ibin] = h->GetBinLowEdge(ibin + 1);
      }
      for (int ibin = excluded_lower_bins; ibin <= NBins; ibin++) {
          X[ibin] = h->GetBinLowEdge(ibin + 1) + shift;
      }

      h->GetXaxis()->Set(NBins, X);
    }

    void myMarkerText(Double_t x, Double_t y, Int_t color, Int_t mstyle, std::string text, Float_t msize, Double_t tsize = 0.06)
    {
        TMarker *marker = new TMarker(x - (0.4 * tsize), y, 8);
        marker->SetMarkerColor(color);  marker->SetNDC();
        marker->SetMarkerStyle(mstyle);
        marker->SetMarkerSize(msize);
        marker->Draw();

        TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
        l.SetNDC();
        l.DrawLatex(x, y, text.c_str());
    }


    void CheckFile(const TFile *file, std::string name) {
      if (file->IsZombie()) {
          std::cout << name << " File is Zombie" << std::endl;
          throw std::exception();
      }
    }

    bool CheckObject(const TObject *obj, std::string name, TFile *file = nullptr) {
        std::string file_description = "";
        if (file) {
            file_description = " in TFile ";
            file_description += file->GetName();
        }
        if (!obj) {
            std::cout << name << " Object is invalid" << file_description << std::endl;
            return false;
        }
        if (strcmp(obj->GetName(), name.c_str())) {
            std::cout << "Expected=" << obj->GetName() << "  Found=" << name << file_description << std::endl;
            return false;
        }
        return true;
    }


    bool CheckObject2(const TObject *obj, std::string name, TFile *file = nullptr) {
        if (!CheckObject(obj, name, file)) throw std::exception();
        else return true;
    }

    std::string UniqueName() {
    static int UniqueNameCounter = 0;
    std::string str = "UniqueName";
    str += std::to_string(UniqueNameCounter);
    UniqueNameCounter++;
    return str;
}

//Symmetrize 1D correlation in dphi
//Only works for a particular binning
void Symmetrize_1D(TH1 *hist) {
  static TH1* hist2 = 0;

  int NBins = hist->GetNbinsX();
  if (NBins != 36) {std::cout << "Error in Symmetrize()" << std::endl; throw std::exception();}

  if (!hist2) hist2 = (TH1*)hist->Clone("hist_temp_symmetrize");

  hist2->Reset();
  for (int ibin = 1; ibin <= 9; ibin++) {
    double val = hist->GetBinContent(ibin);
    double err = hist->GetBinError  (ibin);
    hist2->SetBinContent(18 + 1 - ibin, val);
    hist2->SetBinError  (18 + 1 - ibin, err);

    val = hist->GetBinContent(18 + 1 - ibin);
    err = hist->GetBinError  (18 + 1 - ibin);
    hist2->SetBinContent(ibin, val);
    hist2->SetBinError  (ibin, err);

    val = hist->GetBinContent(18 + ibin);
    err = hist->GetBinError  (18 + ibin);
    hist2->SetBinContent(36 + 1 - ibin, val);
    hist2->SetBinError  (36 + 1 - ibin, err);

    val = hist->GetBinContent(36 + 1 - ibin);
    err = hist->GetBinError  (36 + 1 - ibin);
    hist2->SetBinContent(18 + ibin, val);
    hist2->SetBinError  (18 + ibin, err);
  }
  hist->Add(hist2);
  double sqrt2 = sqrt(2.0);
  for (int ibin = 1; ibin <= 36; ibin++) {
    double val = hist->GetBinContent(ibin);
    double err = hist->GetBinError  (ibin);
    hist->SetBinContent(ibin, val / 2.0);
    hist->SetBinError  (ibin, err / sqrt2);
  }
}


//Symmetrize 2D correlation (both in deta and dphi)
//Only works for a particular binning
TH2D* Symmetrize_2D(TH2D *hist) {

    int NBins_phi = hist->GetNbinsX();
    int NBins_deta = hist->GetNbinsY();
    if (NBins_phi != 36) {std::cout << "Error in Symmetrize()" << std::endl; throw std::exception();}
    if (NBins_deta != 50) {std::cout << "Error in Symmetrize()" << std::endl; throw std::exception();}

    char histname[600];
    sprintf(histname, "%s_sym", hist->GetName());
    TH2D* hist2 = new TH2D(histname, histname, 36, -PI / 2, 1.5 * PI, 100, -5.0, 5.0);

    for (int ibin_deta = 1; ibin_deta <= 50; ibin_deta++) {
        for (int ibin_phi = 1; ibin_phi <= 9; ibin_phi++) {
            double val1 = hist->GetBinContent(ibin_phi     , ibin_deta);
            double err1 = hist->GetBinError  (ibin_phi     , ibin_deta);
            double val2 = hist->GetBinContent(18 + 1 - ibin_phi, ibin_deta);
            double err2 = hist->GetBinError  (18 + 1 - ibin_phi, ibin_deta);
            double val = val1 + val2;
            double err = sqrt(err1 * err1 + err2 * err2);

            hist2->SetBinContent(ibin_phi     , 50 + ibin_deta, val);
            hist2->SetBinError  (ibin_phi     , 50 + ibin_deta, err);
            hist2->SetBinContent(18 + 1 - ibin_phi, 50 + ibin_deta, val);
            hist2->SetBinError  (18 + 1 - ibin_phi, 50 + ibin_deta, err);
            hist2->SetBinContent(ibin_phi     , 51 - ibin_deta, val);
            hist2->SetBinError  (ibin_phi     , 51 - ibin_deta, err);
            hist2->SetBinContent(18 + 1 - ibin_phi, 51 - ibin_deta, val);
            hist2->SetBinError  (18 + 1 - ibin_phi, 51 - ibin_deta, err);


            val1 = hist->GetBinContent(18 +  ibin_phi, ibin_deta);
            err1 = hist->GetBinError  (18 +  ibin_phi, ibin_deta);
            val2 = hist->GetBinContent(36 + 1 - ibin_phi, ibin_deta);
            err2 = hist->GetBinError  (36 + 1 - ibin_phi, ibin_deta);
            val = val1 + val2;
            err = sqrt(err1 * err1 + err2 * err2);

            hist2->SetBinContent(18 +  ibin_phi, 50 + ibin_deta, val);
            hist2->SetBinError  (18 +  ibin_phi, 50 + ibin_deta, err);
            hist2->SetBinContent(36 + 1 - ibin_phi, 50 + ibin_deta, val);
            hist2->SetBinError  (36 + 1 - ibin_phi, 50 + ibin_deta, err);
            hist2->SetBinContent(18 +  ibin_phi, 51 - ibin_deta, val);
            hist2->SetBinError  (18 +  ibin_phi, 51 - ibin_deta, err);
            hist2->SetBinContent(36 + 1 - ibin_phi, 51 - ibin_deta, val);
            hist2->SetBinError  (36 + 1 - ibin_phi, 51 - ibin_deta, err);
        }
    }
    return hist2;
}

void format(TH1* hist, int col = 1, int sty = 20) {
    //hist->GetYaxis()->CenterTitle();
    //hist->GetXaxis()->CenterTitle();

    hist->SetMarkerColor(col);
    hist->SetLineColor  (col);
    hist->SetMarkerStyle(sty);
    // hist->SetMarkerSize(0.5);
}

void format(TGraph* hist, int col = 1, int sty = 20) {
    //hist->GetYaxis()->CenterTitle();
    //hist->GetXaxis()->CenterTitle();

    hist->SetMarkerColor(col);
    hist->SetLineColor  (col);
    hist->SetMarkerStyle(sty);
}

void format_hist(TH1* hist) {
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
}

void FormatHist(TH1* hist1, std::map<std::string, double> f) {
    if (f.find("XLabelSize"  ) != f.end()) hist1->GetXaxis()->SetLabelSize  (f["XLabelSize"]);
    if (f.find("YLabelSize"  ) != f.end()) hist1->GetYaxis()->SetLabelSize  (f["YLabelSize"]);
    if (f.find("XLabelOffset") != f.end()) hist1->GetXaxis()->SetLabelOffset(f["XLabelOffset"]);
    if (f.find("YLabelOffset") != f.end()) hist1->GetYaxis()->SetLabelOffset(f["YLabelOffset"]);
    if (f.find("XTitleOffset") != f.end()) hist1->GetXaxis()->SetTitleOffset(f["XTitleOffset"]);
    if (f.find("YTitleOffset") != f.end()) hist1->GetYaxis()->SetTitleOffset(f["YTitleOffset"]);
    if (f.find("XTitleSize"  ) != f.end()) hist1->GetXaxis()->SetTitleSize  (f["XTitleSize"]);
    if (f.find("YTitleSize"  ) != f.end()) hist1->GetYaxis()->SetTitleSize  (f["YTitleSize"]);
    if (f.find("XNdivisions")  != f.end()) hist1->GetXaxis()->SetNdivisions (int(f["XNdivisions"]));
    if (f.find("YNdivisions")  != f.end()) hist1->GetYaxis()->SetNdivisions (int(f["YNdivisions"]));

    //hist1->GetXaxis()->CenterTitle();
    //hist1->GetYaxis()->CenterTitle();
    hist1->SetStats(0);
    hist1->SetTitle("");
}

TH1* Take_Sqrt( TH1* MyHist, int flag = 0) { //Symmetric errors!
  TH1* NewHist;
  if (flag) { //create new one
    NewHist = (TH1*)MyHist->Clone();
    NewHist->Reset();
  } else {
    NewHist = MyHist;
  }
  int N = NewHist->GetNbinsX();
  for (int I = 1; I <= N; I++) {
    float val1     = MyHist->GetBinContent(I);
    float val = (val1 > 0) ? sqrt(val1) : -sqrt(-val1);
    {
      float val_err = fabs(MyHist->GetBinError  (I));
      float val_up = val1 + val_err;
      float val_dn = val1 - val_err;
      val_up  = (val_up > 0) ? sqrt(val_up) : -sqrt(-val_up);
      val_dn  = (val_dn > 0) ? sqrt(val_dn) : -sqrt(-val_dn);

      float val_err0 = fabs(val_up - val_dn) / 2.0;
      float val_err1 = fabs(val_up - val);
      float val_err2 = fabs(val_dn - val);
      val_err = val_err0;
      if (val_err1 > val_err) val_err = val_err1;
      if (val_err2 > val_err) val_err = val_err2;

      NewHist->SetBinContent(I, val    );
      NewHist->SetBinError  (I, val_err);
    }
  }
  return NewHist;
}

std::map<std::string, double> StandardFormat() {
    std::map<std::string, double> _format = {
        {"XTitleSize"  , 0.06},
        {"YTitleSize"  , 0.06},
        {"XLabelSize"  , 0.05},
        {"YLabelSize"  , 0.05},
        {"XNdivisions" , 505},
        {"YNdivisions" , 505},
        {"YTitleOffset", 1.0},
        {"XTitleOffset", 1.0},
    };
    return _format;
}

}
#endif
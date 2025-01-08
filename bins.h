#ifndef __BINS_H__
#define __BINS_H__
#include <string>
#include <map>
#include <utility>
#include <Riostream.h>

#include "common.C"

namespace Bins{

    std::string HH_LABEL  = "#it{h}-#it{h} Correlations";
    float sigma1 = 1.5;
    bool same_side1 = true;
    int  Initialize();
    void Initialize_CentAdd();
    void Initialize_Pt1Add();
    void Initialize_Pt2Add();

    std::string label_cent      (int icent);
//     std::string label_cent_short(int icent);
//     std::string label_cent_peri (int icent);
//     int GetCentIndex              (int mult_low, int mult_high);
//     std::vector<int> GetCentIndex (std::vector<std::pair<int, int>> input);
//     std::vector<int> GetCentIndex (std::vector<double> input);
//     std::pair<double, double> GetCentVals(int icent);


//     std::string label_pta(int ipt);
//     int GetPtaIndex(double pta_low, double pta_high);
//     std::vector<int> GetPtaIndex (std::vector<std::pair<double, double>> input);
//     std::pair<double, double> GetPtaVals(int ipt);

//     std::string label_ptb(int ipt);
    int GetPtbIndex(double ptb_low, double ptb_high);
    std::vector<int> GetPtbIndex (std::vector<std::pair<double, double>> input);
    //std::vector<int> GetPtbIndex (std::vector<float> input);
//     std::pair<double, double> GetPtbVals(int ipt);

//     std::string label_ptab(int ipta, int iptb);
//     int GetPtbIndexForPtaIndex(int ipta);

//     std::string label_eta(int ideta);
//     int GetDetaIndex(double deta_low, double deta_high);
//     std::vector<int> GetDetaIndex (std::vector<std::pair<double, double>> input);


    int GetPtBin1(float pt);
    int GetPtBin2(float pt);
//     //TODO modify to multiplicity
//     int GetCentrality(float FCal_Et, int data_type);
    int GetCentBin(float cent_onepercent);
    int GetTrktBin(float cent_onepercent);

//     TH1D* PtadepHist(std::vector<int>ipt1_vec, std::string name);
//     TH1D* PtbdepHist(std::vector<int>ipt2_vec, std::string name);
//     TH1D* DetadepHist(std::vector<int>deta_vec, std::string name);

    const int NO_PERIPHERAL_BIN = -1;
//     std::pair<float, float> GetVnn  (int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph, TFile *FourierFile);
//     std::pair<float, float> GetVnPtb(int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph, TFile *FourierFile1, TFile *FourierFile2);

//     std::pair<float, float> GetVnn  (TFile *TemplateFile, int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph);
//     std::pair<float, float> GetVnPtb(TFile *TemplateFile, int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph);


//     int GetPtBin1(float pt);
//     int GetPtBin2(float pt);


//     TH1D* CentdepHist(std::vector<int>centbins, int l_centrality_type, std::string name);
//     TH1D* PtadepHist(std::vector<int>ipt1_vec, std::string name);
//     TH1D* PtbdepHist(std::vector<int>ipt2_vec, std::string name);




//     std::vector<int> CentBins();
//     std::vector<int> CentBinsPeriph();
//     std::vector<int> PtaBins ();
//     std::vector<int> PtbBins ();
//     std::vector<int> DetaBins();
//     std::vector<int> ChBins  ();

    enum {
        BINS_SUMMARY       = 6,
        BIN_FUNC_INTEGRAL  = 1,
        BIN_HIST_INTEGRAL  = 2,
        BIN_FUNC_CHI2_NDOF = 3,
        BIN_ZYAM_MIN_ERR   = 4,
        BIN_ZYAM_ERR_HIST  = 5,
        BIN_ZYAM_ERR_FUNC  = 6,
        BIN_ZYAM_MIN_VAL   = 7,
    };
    std::map<int, std::string> BINS_SUMMARY_LABELS = {
    {BIN_ZYAM_MIN_ERR  , "Error of ZYAM Min"        },
    {BIN_ZYAM_MIN_VAL  , "Pedestal of ZYAM"         },
    {BIN_FUNC_INTEGRAL , "PTY Near Side (from Fit)" },
    {BIN_HIST_INTEGRAL , "PTY Near Side (from Hist)"},
    {BIN_FUNC_CHI2_NDOF, "CHI2/NDOF of Fit"         },
    {BIN_ZYAM_ERR_HIST , "ZYAM Error (hist)"        },
    {BIN_ZYAM_ERR_FUNC , "ZYAM Error (function)"    },
    };

    enum {
        VN_TYPE                         = 7,
        VN_TEMPLATE                     = 0, //OLD with ZYAM
        VN_TEMPLATE_PEDESTAL            = 1, //NEW without ZYAM
        VN_DIRECT                       = 2,
        VN_PERISUB                      = 3,
        VN_TEMPLATE_PEDESTAL_HJHHperiph = 4,
        VN_TEMPLATE_HJHHperiph          = 5,
        VN_TEMPLATE_PEDESTALMINBIAS_HJHHperiph = 6,
    };

    std::string label_vn_type(int vn_type) {
    if      (vn_type == VN_TEMPLATE                    ) return "Template fit (ZYAM)"  ;
    else if (vn_type == VN_TEMPLATE_PEDESTAL           ) return "Template fit"         ;
    else if (vn_type == VN_DIRECT                      ) return "Fourier Transform"     ;
    else if (vn_type == VN_PERISUB                     ) return "Peripheral Subtraction";
    else if (vn_type == VN_TEMPLATE_PEDESTAL_HJHHperiph) return "Template fit";
    else if (vn_type == VN_TEMPLATE_HJHHperiph)          return "Template fit (ZYAM)";
    else if (vn_type == VN_TEMPLATE_PEDESTALMINBIAS_HJHHperiph)          return "Template Fits (jet ZYAM)";
    else {
        std::cout << " " << ":: Unknown vn type" << std::endl;
        throw std::exception();
    }
    }

//     //-------------------------------------------------------------------------------------



    enum {
    NCENT = 20,
    NTRK = 13,
    NPT1 = 5,
    NPT2 = 5,
    NCH  = 2,
    SAME_CHARGE = 0,
    OPPOSITE_CHARGE = 1,

    NCENT_ADD = 11,
    NTRK_ADD = 1,

    NPT1_ADD = 1,
    NPT2_ADD = 2,
    NCH_ADD = 1,
    COMBINED_CHARGE = 2,
    NDETA   = 5,
    };


    char label[600];
    // CENTRALITY BINS
//---------------------------------------------------------------------------------------------------------------------

//LABELS                              0,    1,    2,    3,    4,    5,    6,    7,    8,     9,   10,   11,   12,   13,   14    ,15,   16,     17,    18,    19    
float CENT_LO[NCENT + NCENT_ADD] = { 0.0,  0.68, 1.36, 2.04, 2.72, 3.40, 4.08, 4.76, 5.44, 6.12, 6.80, 7.48, 8.16, 8.84, 9.52, 10.20, 10.88,  11.56, 12.24, 12.92 };  // bins for effective energy 
float CENT_HI[NCENT + NCENT_ADD] = { 0.68, 1.36, 2.04, 2.72, 3.40, 4.08, 4.76, 5.44, 6.12, 6.80, 7.48, 8.16, 8.84, 9.52, 10.20, 10.88, 11.56, 12.24, 12.92, 13.6001 }; //in TeV

float cent_add_lo[NCENT_ADD]     = {0.0};
float cent_add_up[NCENT_ADD]     = {0.0};
void Initialize_CentAdd() {
  std::vector<std::pair<double, double>> new_cent_bins = {
    //20,           21,          22,            23,           24,          25,          26,             27,            28,             29,             30,        35,       36,
    { 0.0, 13.6}, {0.0, 1.36}, {10.88, 12.92}, {0.0, 1.36}, {1.36, 2.72}, {2.72, 4.08}, {4.08, 5.44}, {5.44, 6.8}, {6.8, 8.16}, {8.16, 9.52}, {9.52, 10.88}//, {10.88, 12.24}, {12.24, 13.6}
    //, {0, 60}, { 0, 10},
    //37,       38,       39,       40,       41,       42,       43,       44,       
    //{10, 20}, {30, 40}, {50, 60}, {60, 70}, {70, 80}, {80, 90}, {90, 100}, {85, 95}
  };
    if (new_cent_bins.size() != NCENT_ADD) {
    std::cout << "Initialize_CentAdd()::new_cent_bins.size()!=NCENT_ADD " << new_cent_bins.size() << "  " << NCENT_ADD << std::endl;
    throw std::exception();
  }

  int ibin = 0;
  for (auto new_bin : new_cent_bins) {
    int low = -1, high = -1;
    for (int icent = 0; icent < NCENT; icent++) {
      if (fabs(CENT_LO[icent] - new_bin.first ) < 0.001 ) low = icent;
      if (fabs(CENT_HI[icent] - new_bin.second) < 0.001 ) high = icent + 1;
    }
    //cout << low << " " << high <<endl;
    if (low == -1 || high == -1 || low >= high) {
      std::cout << "Initialize_CentAdd():: Problem adding new bin" << std::endl;
      throw std::exception();
    }

    cent_add_lo[ibin] = low;
    cent_add_up[ibin] = high;
    CENT_LO[NCENT + ibin] = new_bin.first;
    CENT_HI[NCENT + ibin] = new_bin.second;
    ibin++;
  }
}

std::string label_cent(int icent) {
  sprintf(label, "%.2f < E_{Eff} (TeV) < %.2f", CENT_LO[icent], CENT_HI[icent]);
  std::string ret = label;
  return ret;
}
int GetCentIndex(double cent_low, double cent_high) {
  for (int index = 0; index < NCENT + NCENT_ADD; index++) {
    if ( fabs(CENT_LO[index] - cent_low) < 0.001 && fabs(CENT_HI[index] - cent_high) < 0.001 ) return index;
  }
  std::cout << "This Centrality does not exist " << cent_low << "," << cent_high << std::endl;
  throw std::exception();
}
std::vector<int> GetCentIndex(std::vector<std::pair<float, float>> input) {
  std::vector<int> ret;
  for (auto& mult_bin : input) {
    ret.push_back(GetCentIndex(mult_bin.first, mult_bin.second));
  }
  return ret;
}
std::vector<int> GetCentIndex(std::vector<float> input) {
  std::vector<int> ret;
  for (int i = 0; i < (int)(input.size()) - 1; i++) {
    ret.push_back(GetCentIndex(input[i], input[i + 1]));
  }
  return ret;
}
std::pair<double, double> GetCentVals(int icent) {
  if (icent >= (NCENT + NCENT_ADD) || icent < 0) {
    std::cout << "This Centrality index doesnot exist :" << icent << std::endl;
    throw std::exception();
  }
  return std::make_pair(CENT_LO[icent], CENT_HI[icent]);
}
//---------------------------------------------------------------------------------------------------------------------

// Multiplicity BINS
//---------------------------------------------------------------------------------------------------------------------
//LABELS                        0,  1,  2,  3,  4,  5,  6,  7,  8,  9,   10,  11,  12
int TRK_LO[NTRK + NTRK_ADD] = { 0,  10, 20, 30, 40, 50, 60, 70, 80, 90 , 100, 110, 120  }; 
int TRK_HI[NTRK + NTRK_ADD] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130  };
int trk_add_lo[NTRK_ADD]     = {0};
int trk_add_up[NTRK_ADD]     = {0};


void Initialize_TrkAdd() {
  std::vector<std::pair<double, double>> new_trk_bins = {
    //13,      14,       14,        15,        28,       29,      30,       31,       32,       33,       34,        35,       36,
    { 0, 130}//, {0, 130}, //{ 0, 100},{0, 20}, {0, 20}, {20, 40}, {40, 60}, {60, 80}, {80, 100}, {20, 100}, {0, 60}, { 0, 10},
    //37,       38,       39,       40,       41,       42,       43,       44,       
    //{10, 20}, {30, 40}, {50, 60}, {60, 70}, {70, 80}, {80, 90}, {90, 100}, {85, 95}
  };
    if (new_trk_bins.size() != NTRK_ADD) {
    std::cout << "Initialize_TrkAdd()::new_trk_bins.size()!=NTRK_ADD " << new_trk_bins.size() << "  " << NTRK_ADD << std::endl;
    throw std::exception();
  }

  int ibin = 0;
  for (auto new_bin : new_trk_bins) {
    int low = -1, high = -1;
    for (int itrk = 0; itrk < NTRK; itrk++) {
      if (fabs(TRK_LO[itrk] - new_bin.first ) < 0.001 ) low = itrk;
      if (fabs(TRK_HI[itrk] - new_bin.second) < 0.001 ) high = itrk + 1;
    }
    
    if (low == -1 || high == -1 || low >= high) {
      std::cout << "Initialize_TrkAdd():: Problem adding new bin" << std::endl;
      throw std::exception();
    }

    trk_add_lo[ibin] = low;
    trk_add_up[ibin] = high;
    TRK_LO[NTRK + ibin] = new_bin.first;
    TRK_HI[NTRK + ibin] = new_bin.second;
    ibin++;
  }
}

std::string label_trk(int itrk) {
  sprintf(label, "%d #leq N_{ch}^{rec} <%d", TRK_LO[itrk], TRK_HI[itrk]);
  std::string ret = label;
  return ret;
}
int GetTrkIndex(double trk_low, double trk_high) {
  for (int index = 0; index < NTRK + NTRK_ADD; index++) {
    if ( fabs(TRK_LO[index] - trk_low) < 0.001 && fabs(TRK_HI[index] - trk_high) < 0.001 ) return index;
  }
  std::cout << "This multiplicity does not exist " << trk_low << "," << trk_high << std::endl;
  throw std::exception();
}
std::vector<int> GetTrkIndex(std::vector<std::pair<int, int>> input) {
  std::vector<int> ret;
  for (auto& mult_bin : input) {
    ret.push_back(GetTrkIndex(mult_bin.first, mult_bin.second));
  }
  return ret;
}
std::vector<int> GetTrkIndex(std::vector<int> input) {
  std::vector<int> ret;
  for (int i = 0; i < (int)(input.size()) - 1; i++) {
    ret.push_back(GetTrkIndex(input[i], input[i + 1]));
  }
  return ret;
}
std::pair<double, double> GetTrkVals(int itrk) {
  if (itrk >= (NTRK + NTRK_ADD) || itrk < 0) {
    std::cout << "This multiplicity index doesnot exist :" << itrk << std::endl;
    throw std::exception();
  }
  return std::make_pair(TRK_LO[itrk], TRK_HI[itrk]);
}
//---------------------------------------------------------------------------------------------------------------------

// PTa BINS
//LABELS                           0,   1,   2,   3,   4,   
float PT1_LO[NPT1 + NPT1_ADD] = {0.5, 1.0, 2.0, 3.0, 4.0,};
float PT1_HI[NPT1 + NPT1_ADD] = {1.0, 2.0, 3.0, 4.0, 5.0,};
//BINS
int pt1_add_lo[NPT1_ADD] = {0};
int pt1_add_up[NPT1_ADD] = {0};
void Initialize_Pt1Add() {
  std::vector<std::pair<double, double>> new_pt1_bins = {
    //  5          6          7           8
    {0.5, 5.0}//, {1.0, 5.0}, {2.0, 5.0}, {0.5, 3.0},
  };

  if (new_pt1_bins.size() != NPT1_ADD) {
    std::cout << "Initialize_Pt1Add()::new_pt1_bins.size()!=NPT1_ADD " << new_pt1_bins.size() << "  " << NPT1_ADD << std::endl;
    throw std::exception();
  }

  int ibin = 0;
  for (auto new_bin : new_pt1_bins) {
    int low = -1, high = -1;
    for (int ipt1 = 0; ipt1 < NPT1; ipt1++) {
      if (fabs(PT1_LO[ipt1] - new_bin.first ) < 0.001 ) low = ipt1;
      if (fabs(PT1_HI[ipt1] - new_bin.second) < 0.001 ) high = ipt1 + 1;
    }
    if (low == -1 || high == -1 || low >= high) {
      std::cout << "Initialize_Pt1Add():: Problem adding new bin" << std::endl;
      throw std::exception();
    }

    pt1_add_lo[ibin] = low;
    pt1_add_up[ibin] = high;
    PT1_LO[NPT1 + ibin] = new_bin.first;
    PT1_HI[NPT1 + ibin] = new_bin.second;
    ibin++;
  }
}

std::string label_pta(int ipt) {
  sprintf(label, "%.1f<#it{p}_{T}^{#it{a}}(GeV)<%.1f", PT1_LO[ipt], PT1_HI[ipt]);
  if (PT1_LO[ipt] == int(PT1_LO[ipt]))   sprintf(label, "%d<#it{p}_{T}^{#it{a}}(GeV)<%.1f", int(PT1_LO[ipt]),   (PT1_HI[ipt]));
  if (PT1_HI[ipt] == int(PT1_HI[ipt]))   sprintf(label, "%.1f<#it{p}_{T}^{#it{a}}(GeV)<%d",   (PT1_LO[ipt]), int(PT1_HI[ipt]));
  if (PT1_LO[ipt] == int(PT1_LO[ipt])
      && PT1_HI[ipt] == int(PT1_HI[ipt]))   sprintf(label, "%d<#it{p}_{T}^{#it{a}}(GeV)<%d"  , int(PT1_LO[ipt]), int(PT1_HI[ipt]));

  std::string ret = label;
  return ret;
}
int GetPtaIndex(double pta_low, double pta_high) {
  for (int index = 0; index < NPT1 + NPT1_ADD; index++) {
    if ( fabs(PT1_LO[index] - pta_low) < 0.001 && fabs(PT1_HI[index] - pta_high) < 0.001 ) return index;
  }
  std::cout << "This pta doesnot exist " << pta_low << "," << pta_high << std::endl;
  throw std::exception();
}
std::vector<int> GetPtaIndex (std::vector<std::pair<double, double>> input) {
  std::vector<int> ret;
  for (auto& pta_bin : input) {
    ret.push_back(GetPtaIndex(pta_bin.first, pta_bin.second));
  }
  return ret;
}
std::pair<double, double> GetPtaVals(int ipt) {
  return {PT1_LO[ipt], PT1_HI[ipt]};
}
//---------------------------------------------------------------------------------------------------------------------

//PTb Bins
//---------------------------------------------------------------------------------------------------------------------
// 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,   19,  20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,   33,   34,   35,   36,   37,   38
float PT2_LO[NPT2 + NPT2_ADD] =
{0.5, 1.0, 2.0, 3.0, 4.0//, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,  10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5
};
float PT2_HI[NPT2 + NPT2_ADD] =
{1.0, 2.0, 3.0, 4.0, 5.0//, 3.5, 4.0, 4.5, 5.0, //5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0
};
int pt2_add_lo[NPT2_ADD] = {0};
int pt2_add_up[NPT2_ADD] = {0};

void Initialize_Pt2Add() {
  std::vector<std::pair<double, double>> new_pt2_bins = {
    // 5         6           7           8          49          50          51          52          53
    {0.5, 5.0}, {3.0, 5.0}//, {2.0, 3.0}, {3.0, 5.0}, //{0.5, 5.0}, {1.0, 5.0}, {2.0, 5.0}, {0.4, 3.0}, {1.0, 3.0},
    // 54         55          56          57          58          59          60          61            62
    //{1.4, 1.6}, {2.4, 2.6}, {1.0, 1.5}, {1.5, 2.0}, {3.0, 4.0}, {4.0, 6.0}, {6.0, 8.0}, {8.0, 10.0}, {10.0, 15},
    // 63         64          65          66          67          68          69          70          71        72
    //{3.0, 3.5}, {3.5, 4.0}, {4.0, 4.5}, {4.5, 5.0}, {5.0, 5.5}, {5.5, 6.0}, {6.0, 7.0}, {7.0, 8.0},  {12, 15}, {15, 20},
    // 73         74          75          76          77         78          79          80          81        82
    //{3.4, 3.6}, {4.4, 4.6}, {5.4, 5.9}, {5.9, 6.5}, {6.5, 7.2}, {10, 14}  , {14, 18}  , {18, 25}  ,  {5, 6}, {1.2, 1.5},
    // 83         84           85          86
    //{4.0, 5.0}, {5.0, 8.0}, {2.0, 2.5}, {2.5, 3.0}, {8.0, 20.0}, {10.0, 20.0}, {6.0, 10.0}, {4.0, 8.0}
  };

  if (new_pt2_bins.size() != NPT2_ADD) {
    std::cout << "Initialize_Pt2Add()::new_pt2_bins.size()!=NPT2_ADD " << new_pt2_bins.size() << "  " << NPT2_ADD << std::endl;
    throw std::exception();
  }

  int ibin = 0;
  for (auto new_bin : new_pt2_bins) {
    int low = -1, high = -1;
    for (int ipt2 = 0; ipt2 < NPT2; ipt2++) {
      if (fabs(PT2_LO[ipt2] - new_bin.first ) < 0.001 ) low = ipt2;
      if (fabs(PT2_HI[ipt2] - new_bin.second) < 0.001 ) high = ipt2 + 1;
    }
    if (low == -1 || high == -1 || low >= high) {
      std::cout << "Initialize_Pt2Add():: Problem adding new bin" << std::endl;
      throw std::exception();
    }

    pt2_add_lo[ibin] = low;
    pt2_add_up[ibin] = high;
    PT2_LO[NPT2 + ibin] = new_bin.first;
    PT2_HI[NPT2 + ibin] = new_bin.second;
    ibin++;
  }

  //checks
  //make sure edge of PT2_LO matches edge of previous PT2_HI for original NPT2 bins
  for (int ipt2 = 1; ipt2 < NPT2; ipt2++) {
    if (PT2_LO[ipt2] != PT2_HI[ipt2 - 1]) {
      std::cout << "Error in (PT2_LO,PT2_HI) original at index " << ipt2 << std::endl;
      throw std::exception();
    }
  }
  //make sure that no duplicate pt2 bins are produced
  for (int ipt2 = NPT2; ipt2 < NPT2 + NPT2_ADD; ipt2++) {
    if (GetPtbIndex(PT2_LO[ipt2], PT2_HI[ipt2]) != ipt2) {
      std::cout << "Duplicate PT2 bin found at index " << ipt2 << std::endl;
      throw std::exception();
    }
  }
}

std::string label_ptb(int ipt) {
  sprintf(label, "%.1f<#it{p}_{T}^{#it{b}}(GeV)<%.1f", PT2_LO[ipt], PT2_HI[ipt]);
  if (PT2_LO[ipt] == int(PT2_LO[ipt]))   sprintf(label, "%d<#it{p}_{T}^{#it{b}}(GeV)<%.1f", int(PT2_LO[ipt]),   (PT2_HI[ipt]));
  if (PT2_HI[ipt] == int(PT2_HI[ipt]))   sprintf(label, "%.1f<#it{p}_{T}^{#it{b}}(GeV)<%d",   (PT2_LO[ipt]), int(PT2_HI[ipt]));
  if (PT2_LO[ipt] == int(PT2_LO[ipt])
      && PT2_HI[ipt] == int(PT2_HI[ipt]))   sprintf(label, "%d<#it{p}_{T}^{#it{b}}(GeV)<%d"  , int(PT2_LO[ipt]), int(PT2_HI[ipt]));
  std::string ret = label;
  return ret;
}
int GetPtbIndex(double ptb_low, double ptb_high) {
  for (int index = 0; index < NPT2 + NPT2_ADD; index++) {
    if ( fabs(PT2_LO[index] - ptb_low) < 0.001 && fabs(PT2_HI[index] - ptb_high) < 0.001 ) return index;
  }
  std::cout << "This ptb doesnot exist " << ptb_low << "," << ptb_high << std::endl;
  throw std::exception();
}
std::vector<int> GetPtbIndex (std::vector<std::pair<double, double>> input) {
  std::vector<int> ret;
  for (auto& ptb_bin : input) {
    ret.push_back(GetPtbIndex(ptb_bin.first, ptb_bin.second));
  }
  return ret;
}
std::vector<int> GetPtbIndex(std::vector<float> input) {
  std::vector<int> ret;
  for (int i = 0; i < (int)(input.size()) - 1; i++) {
    ret.push_back(GetPtbIndex(input[i], input[i + 1]));
  }
  return ret;
}
std::pair<double, double> GetPtbVals(int ipt) {
  return {PT2_LO[ipt], PT2_HI[ipt]};
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
std::string label_ptab(int ipta, int iptb) {
  if (PT1_LO[ipta] == PT2_LO[iptb] && PT1_HI[ipta] == PT2_HI[iptb]) {
    sprintf(label, "%.1f<#it{p}_{T}^{#it{a},#it{b}}(GeV)<%.1f", PT2_LO[iptb], PT2_HI[iptb]);
    if (PT1_LO[ipta] == int(PT1_LO[ipta]))   sprintf(label, "%d<#it{p}_{T}^{#it{a},#it{b}}(GeV)<%.1f", int(PT1_LO[ipta]),   (PT1_HI[ipta]));
    if (PT1_HI[ipta] == int(PT1_HI[ipta]))   sprintf(label, "%.1f<#it{p}_{T}^{#it{a},#it{b}}(GeV)<%d",   (PT1_LO[ipta]), int(PT1_HI[ipta]));
    if (PT1_LO[ipta] == int(PT1_LO[ipta])
        && PT1_HI[ipta] == int(PT1_HI[ipta]))   sprintf(label, "%d<#it{p}_{T}^{#it{a},#it{b}}(GeV)<%d"  , int(PT1_LO[ipta]), int(PT1_HI[ipta]));
  }
  else {
    std::cout << "label_ptab() called with unpmatched pta and ptb bins" << ipta << "  " << iptb << std::endl; throw std::exception();
  }
  std::string ret = label;
  return ret;
}
int GetPtbIndexForPtaIndex(int ipta) {
  return GetPtbIndex(PT1_LO[ipta], PT1_HI[ipta]);
}
//---------------------------------------------------------------------------------------------------------------------

//Charge bins
//---------------------------------------------------------------------------------------------------------------------
//BINS
int ch_add_lo[NCH_ADD] = {0};
int ch_add_up[NCH_ADD] = {2};
//---------------------------------------------------------------------------------------------------------------------

//Delta eta
//LABELS
//---------------------------------------------------------------------------------------------------------------------
float DETA_LO[NDETA] = {0.0, 2.0, 3.0, 0.0, 4.0};
float DETA_HI[NDETA] = {2.0, 5.0, 5.0, 1.0, 5.0};
//BINS
std::string label_eta(int ideta) {
  sprintf(label, "%.1f<|#Delta#eta|<%.1f", DETA_LO[ideta], DETA_HI[ideta]);
  if (DETA_LO[ideta] == int(DETA_LO[ideta]) && DETA_HI[ideta] == int(DETA_HI[ideta])) {
    sprintf(label, "%d<|#Delta#eta|<%d", int(DETA_LO[ideta]), int(DETA_HI[ideta]));
  }
  std::string ret = label;
  return ret;
}
int GetDetaIndex(double deta_low, double deta_high) {
  for (int index = 0; index < NDETA; index++) {
    if ( fabs(DETA_LO[index] - deta_low) < 0.001 && fabs(DETA_HI[index] - deta_high) < 0.001 ) return index;
  }
  std::cout << "This deta doesnot exist " << deta_low << "," << deta_high << std::endl;
  throw std::exception();
}
std::vector<int> GetDetaIndex (std::vector<std::pair<double, double>> input) {
  std::vector<int> ret;
  for (auto& deta_bin : input) {
    ret.push_back(GetDetaIndex(deta_bin.first, deta_bin.second));
  }
  return ret;
}
//---------------------------------------------------------------------------------------------------------------------

// //---------------------------------------------------------------------------------------------------------------------
// enum {
//   NDET = 2, //(FullCal, FCal)
//   TOTAL_CALORIMETER  = 0,
//   FORWARD_CALORIMETER = 1,
//   NSIDE = 3, //(+ve,-ve,combined)
//   POSITIVE = 0,
//   NEGATIVE = 1,
//   COMBINED = 2,
//   NHAR = 2, //V2-V3

//   NHAR_DFT = 6, //v1-v6
//   v1 = 1,
//   v2 = 2,
//   v3 = 3,
//   v4 = 4,
//   v5 = 5,
//   v6 = 6,
// };


//used in Anacorr code only
int GetPtBin1(float pt) {
  for (int i = 0; i < NPT1; i++) {
    if (pt > PT1_LO[i] && pt < PT1_HI[i]) return i;
  }
  return -1;
}

//used in Anacorr code only
int GetPtBin2(float pt) {
  for (int i = 0; i < NPT2; i++) {
    if (pt > PT2_LO[i] && pt < PT2_HI[i]) return i;
  }
  return -1;
}


int GetCentBin(float cent_onepercent) {
  for (int i = 0; i < NCENT; i++) {
    if (cent_onepercent >= CENT_LO[i] && cent_onepercent < CENT_HI[i]) return i;
  }
  return -1;
}

int GetTrkBin(int ntrk) {
  for (int i = 0; i < NTRK; i++) {
    if (ntrk >= TRK_LO[i] && ntrk < TRK_HI[i]) return i;
  }
  return -1;
}


//---------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Create a centrality dependence histogram from vector of bin-numbers
TH1D* CentdepHist(std::vector<int>centbins, std::string name) {
  double Xbins[1000] = {0.0};
  int NXbins = centbins.size();

  for (int ibin = 0; ibin < NXbins; ibin++) {
    Xbins[ibin] = CENT_LO[centbins[ibin]];

    if (ibin > 0) {
      if (Xbins[ibin] != CENT_HI[centbins[ibin - 1]]) {
        std::cout << " " << " Error higher edge of lower bin doesnot match lower edge of next bin "
                  << ibin << "  " << centbins[ibin] << "  " << centbins[ibin - 1] << std::endl;
        throw std::exception();
      }
    }
  }
  Xbins[NXbins] = CENT_HI[centbins[NXbins - 1]];

  TH1D* hist = new TH1D(name.c_str(), "", NXbins, Xbins);
  hist->GetXaxis()->SetTitle("E_{Eff} [TeV]");
  //hist->GetXaxis()->SetTitle("N_{ch}");

  return hist;
}

//Create a centrality dependence histogram from vector of bin-numbers
TH1D* TrkdepHist(std::vector<int>trkbins, std::string name) {
  double Xbins[1000] = {0.0};
  int NXbins = trkbins.size();

  for (int ibin = 0; ibin < NXbins; ibin++) {
    Xbins[ibin] = TRK_LO[trkbins[ibin]];

    if (ibin > 0) {
      if (Xbins[ibin] != TRK_HI[trkbins[ibin - 1]]) {
        std::cout << " " << " Error higher edge of lower bin doesnot match lower edge of next bin "
                  << ibin << "  " << trkbins[ibin] << "  " << trkbins[ibin - 1] << std::endl;
        throw std::exception();
      }
    }
  }
  Xbins[NXbins] = TRK_HI[trkbins[NXbins - 1]];

  TH1D* hist = new TH1D(name.c_str(), "", NXbins, Xbins);
  hist->GetXaxis()->SetTitle("N_{ch}^{rec}");

  return hist;
}


// //Create a pta dependence histogram from vector of bin-numbers
// TH1D* PtadepHist(std::vector<int>ipt1_vec, std::string name) {
//   double Xbins[1000] = {0.0};
//   int NXbins = ipt1_vec.size();

//   for (int ibin = 0; ibin < NXbins; ibin++) {
//     Xbins[ibin] = PT1_LO[ipt1_vec[ibin]];
//     if (ibin > 0) {
//       if (Xbins[ibin] != PT1_HI[ipt1_vec[ibin - 1]]) {
//         std::cout << " " << " Error higher edge of lower bin doesnot match lower edge of next bin "
//                   << ibin << "  " << ipt1_vec[ibin] << "  " << ipt1_vec[ibin - 1] << std::endl;
//         throw std::exception();
//       }
//     }
//   }
//   Xbins[NXbins] = PT1_HI[ipt1_vec[NXbins - 1]];
//   TH1D* hist = new TH1D(name.c_str(), "", NXbins, Xbins);
//   hist->GetXaxis()->SetTitle("#it{p}_{T}^{#it{a}} [GeV]");
//   return hist;
// }


// //Create a ptb dependence histogram from vector of bin-numbers
// TH1D* PtbdepHist(std::vector<int>ipt2_vec, std::string name) {
//   double Xbins[1000] = {0.0};
//   int NXbins = ipt2_vec.size();

//   for (int ibin = 0; ibin < NXbins; ibin++) {
//     Xbins[ibin] = PT2_LO[ipt2_vec[ibin]];
//     if (ibin > 0) {
//       if (Xbins[ibin] != PT2_HI[ipt2_vec[ibin - 1]]) {
//         std::cout << " " << " Error higher edge of lower bin doesnot match lower edge of next bin "
//                   << ibin << "  " << ipt2_vec[ibin] << "  " << ipt2_vec[ibin - 1] << std::endl;
//         throw std::exception();
//       }
//     }
//   }
//   Xbins[NXbins] = PT2_HI[ipt2_vec[NXbins - 1]];
//   TH1D* hist = new TH1D(name.c_str(), "", NXbins, Xbins);
//   hist->GetXaxis()->SetTitle("#it{p}_{T}^{#it{b}} [GeV]");
//   return hist;
// }


// //Create a Deta dependence histogram from vector of bin-numbers
// TH1D* DetadepHist(std::vector<int>deta_vec, std::string name) {
//   double Xbins[1000] = {0.0};
//   int NXbins = deta_vec.size();

//   for (int ibin = 0; ibin < NXbins; ibin++) {
//     Xbins[ibin] = DETA_LO[deta_vec[ibin]];
//     if (ibin > 0) {
//       if (Xbins[ibin] != DETA_HI[deta_vec[ibin - 1]]) {
//         std::cout << "Error DetadepHist(): higher edge of lower bin doesnot match lower edge of next bin "
//                   << ibin << "  " << deta_vec[ibin] << "  " << deta_vec[ibin - 1] << std::endl;
//         throw std::exception();
//       }
//     }
//   }
//   Xbins[NXbins] = DETA_HI[deta_vec[NXbins - 1]];
//   TH1D* hist = new TH1D(name.c_str(), "", NXbins, Xbins);
//   hist->GetXaxis()->SetTitle("#Delta#eta");
//   return hist;
// }
// //------------------------------------------------------------------------------


//Returns the vnn and its stat error for a given index
std::pair<float, float> GetVnn(int icent,int itrk, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph,int itrk_periph, TFile *FourierFile) {
  std::pair<float, float> ret;
  char name[600];

  if (icent_periph == Bins::NO_PERIPHERAL_BIN) sprintf(name, "h_v%d%d_pta%d_ptb%.2d_ch%d_deta%.2d"             , ihar, ihar,               ipt1, ipt2, ich, ideta);
  else                                         sprintf(name, "h_v%d%d_pericent%.2d_peritrk%.2d_trk%.2d_pta%d_ptb%.2d_ch%d_deta%.2d", ihar, ihar, icent_periph,itrk_periph,itrk, ipt1, ipt2, ich, ideta);
  TH1* hVnn = (TH1*)FourierFile->Get(name);
  if (!Common::CheckObject(hVnn, name, FourierFile)) throw std::exception();

  ret.first    = hVnn->GetBinContent(icent + 1);
  ret.second   = hVnn->GetBinError  (icent + 1);
#ifdef DOUBLE_ERROR
  if (PT1_LO[ipt1] == PT2_LO[ipt2] && PT1_HI[ipt1] == PT2_HI[ipt2]) ret.second *= sqrt(2.0);
#endif
  return ret;
}


//Returns the vn{ptb} and its stat error for a given index
std::pair<float, float> GetVnPtb(int icent,int itrk, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph,int itrk_periph, TFile *FourierFile1, TFile *FourierFile2) {
  static TH1D *h_vnn = nullptr, *h_vnn_diag = nullptr;
  if (!h_vnn) {
    h_vnn     = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
    h_vnn_diag = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
  }


  int ipt2_for_ipt1 = GetPtbIndex(PT1_LO[ipt1], PT1_HI[ipt1]);
  std::pair<float, float> vnn_diag = GetVnn(icent,itrk, ipt1, ipt2_for_ipt1, ich, ideta, ihar, icent_periph,itrk_periph, FourierFile2);
  h_vnn_diag->SetBinContent(1, vnn_diag.first);
  h_vnn_diag->SetBinError  (1, vnn_diag.second);
  Common::Take_Sqrt(h_vnn_diag);


  //for diagonal input return sqrt
  if (ipt2 == ipt2_for_ipt1 && FourierFile2 == FourierFile1) {
    std::pair<float, float> ret_diag;
    ret_diag.first    = h_vnn_diag->GetBinContent(1);
    ret_diag.second   = h_vnn_diag->GetBinError  (1);
    return ret_diag;
  }
  //for off-diagonal input
  else {
    std::pair<float, float> vnn = GetVnn(icent,itrk, ipt1, ipt2, ich, ideta, ihar, icent_periph,itrk_periph, FourierFile1);
    h_vnn->SetBinContent(1, vnn.first);
    h_vnn->SetBinError  (1, vnn.second);


    std::pair<float, float> ret;
    h_vnn->Divide(h_vnn_diag);
    ret.first    = h_vnn->GetBinContent(1);
    ret.second   = h_vnn->GetBinError  (1);
    return ret;
  }
}

std::pair<float, float> GetVnn(TFile *TemplateFile, int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph) {
  std::pair<float, float> ret;
  char name[600];
  sprintf(name, "h_v%d%d_PPperiph%d_pta%d_ptb%.2d_ch%d_deta%.2d", ihar, ihar, icent_periph, ipt1, ipt2, ich, ideta);
  // std::cout << "------> " << name << std::endl;
  TH1* hVnn = (TH1*)TemplateFile->Get(name);
  if (!Common::CheckObject(hVnn, name)) {
    std::cout << " " << std::endl;
    throw std::exception();
  }
  ret.first    = hVnn->GetBinContent(icent + 1);
  ret.second   = hVnn->GetBinError  (icent + 1);
#ifdef DOUBLE_ERROR
  if (PT1_LO[ipt1] == PT2_LO[ipt2] && PT1_HI[ipt1] == PT2_HI[ipt2]) ret.second *= sqrt(2.0);
#endif
  return ret;
}

std::pair<float, float> GetVnPtb(TFile *TemplateFile, int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph=0) {
  static TH1D *h_vnn = nullptr, *h_vnn_diag = nullptr;
  if (!h_vnn) {
    h_vnn      = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
    h_vnn_diag = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
  }


  int ipt2_for_ipt1 = GetPtbIndex(PT1_LO[ipt1], PT1_HI[ipt1]);
  std::pair<float, float> vnn_diag = GetVnn(TemplateFile, icent, ipt1, ipt2_for_ipt1, ich, ideta, ihar, icent_periph);
  h_vnn_diag->SetBinContent(1, vnn_diag.first);
  h_vnn_diag->SetBinError  (1, vnn_diag.second);
  Common::Take_Sqrt(h_vnn_diag);


  //for diagonal input return sqrt
  if (ipt2 == ipt2_for_ipt1) {
    std::pair<float, float> ret_diag;
    ret_diag.first    = h_vnn_diag->GetBinContent(1);
    ret_diag.second   = h_vnn_diag->GetBinError  (1);
    return ret_diag;
  }
  //for off-diagonal input
  else {
    std::pair<float, float> vnn = GetVnn(TemplateFile, icent, ipt1, ipt2, ich, ideta, ihar, icent_periph);
    h_vnn->SetBinContent(1, vnn.first);
    h_vnn->SetBinError  (1, vnn.second);


    std::pair<float, float> ret;
    h_vnn->Divide(h_vnn_diag);
    ret.first    = h_vnn->GetBinContent(1);
    ret.second   = h_vnn->GetBinError  (1);
    return ret;
  }
}

std::string label_cent_peri(int icent) {
  //if (NTRACK_PPB_HI[icent] < 500) sprintf(label, "%d#leq#it{N}_{ ch}^{ pp,periph}<%d", NTRACK_PPB_LO[icent], NTRACK_PPB_HI[icent]);
  //else                            
  sprintf(label, "%.2f#leq E_{Eff} < %.2f TeV"   , CENT_LO[icent],CENT_HI[icent]);
  std::string ret = label;
  return ret;
}
std::string label_trk_peri(int icent) {
  //if (NTRACK_PPB_HI[icent] < 500) sprintf(label, "%d#leq#it{N}_{ ch}^{ pp,periph}<%d", NTRACK_PPB_LO[icent], NTRACK_PPB_HI[icent]);
  //else                            
  sprintf(label, "%d < N_{ch}^{rec} < %d "   , TRK_LO[icent],TRK_HI[icent]);
  if(icent == NTRK) sprintf(label, "%d < N_{ch}^{rec}"   , TRK_LO[icent]); //highest peripheral nch bin
  std::string ret = label;
  return ret;
}

// std::pair<float, float> GetVnn_selfNtrk(TFile *TemplateFile, int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph) {
//   std::pair<float, float> ret;
//   char name[600];
//   sprintf(name, "h_v%d%d_selfNtrkperiph%d_pta%d_ptb%.2d_ch%d_deta%.2d", ihar, ihar, icent_periph, ipt1, ipt2, ich, ideta);
//   // std::cout << "------> " << name << std::endl;
//   TH1* hVnn = (TH1*)TemplateFile->Get(name);
//   if (!Common::CheckObject(hVnn, name)) {
//     std::cout << " " << std::endl;
//     throw std::exception();
//   }
//   ret.first    = hVnn->GetBinContent(icent + 1);
//   ret.second   = hVnn->GetBinError  (icent + 1);
// #ifdef DOUBLE_ERROR
//   if (PT1_LO[ipt1] == PT2_LO[ipt2] && PT1_HI[ipt1] == PT2_HI[ipt2]) ret.second *= sqrt(2.0);
// #endif
//   return ret;
// }

// std::pair<float, float> GetVnPtb_selfNtrk(TFile *TemplateFile, int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, int icent_periph) {
//   static TH1D *h_vnn = nullptr, *h_vnn_diag = nullptr;
//   if (!h_vnn) {
//     h_vnn      = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
//     h_vnn_diag = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
//   }


//   int ipt2_for_ipt1 = GetPtbIndex(PT1_LO[ipt1], PT1_HI[ipt1]);
//   std::pair<float, float> vnn_diag = GetVnn_selfNtrk(TemplateFile, icent, ipt1, ipt2_for_ipt1, ich, ideta, ihar, icent_periph);
//   h_vnn_diag->SetBinContent(1, vnn_diag.first);
//   h_vnn_diag->SetBinError  (1, vnn_diag.second);
//   Common::Take_Sqrt(h_vnn_diag);


//   //for diagonal input return sqrt
//   if (ipt2 == ipt2_for_ipt1) {
//     std::pair<float, float> ret_diag;
//     ret_diag.first    = h_vnn_diag->GetBinContent(1);
//     ret_diag.second   = h_vnn_diag->GetBinError  (1);
//     return ret_diag;
//   }
//   //for off-diagonal input
//   else {
//     std::pair<float, float> vnn = GetVnn(TemplateFile, icent, ipt1, ipt2, ich, ideta, ihar, icent_periph);
//     h_vnn->SetBinContent(1, vnn.first);
//     h_vnn->SetBinError  (1, vnn.second);


//     std::pair<float, float> ret;
//     h_vnn->Divide(h_vnn_diag);
//     ret.first    = h_vnn->GetBinContent(1);
//     ret.second   = h_vnn->GetBinError  (1);
//     return ret;
//   }
// }


// //Returns the vn{pta} and its stat error for a given index
// std::pair<float, float> GetVnPta(int icent, int ipt1, int ipt2, int ich, int ideta, int ihar, TFile *FourierFile, int icent_periph) {
//   static TH1D *h_vnn = nullptr, *h_vnn_diag = nullptr;
//   if (!h_vnn) {
//     h_vnn     = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
//     h_vnn_diag = new TH1D(Common::UniqueName().c_str(), "", 1, 0, 1);
//   }


//   int ipt1_for_ipt2 = GetPtaIndex(PT2_LO[ipt2], PT2_HI[ipt2]);
//   std::pair<float, float> vnn_diag = GetVnn(icent, ipt1_for_ipt2, ipt2, ich, ideta, ihar, FourierFile, icent_periph);
//   h_vnn_diag->SetBinContent(1, vnn_diag.first);
//   h_vnn_diag->SetBinError  (1, vnn_diag.second);
//   Common::Take_Sqrt(h_vnn_diag);


//   //for diagonal input return sqrt
//   if (ipt1 == ipt1_for_ipt2) {
//     std::pair<float, float> ret_diag;
//     ret_diag.first    = h_vnn_diag->GetBinContent(1);
//     ret_diag.second   = h_vnn_diag->GetBinError  (1);
//     return ret_diag;
//   }
//   //for off-diagonal input
//   else {
//     std::pair<float, float> vnn = GetVnn(icent, ipt1, ipt2, ich, ideta, ihar, FourierFile, icent_periph);
//     h_vnn->SetBinContent(1, vnn.first);
//     h_vnn->SetBinError  (1, vnn.second);


//     std::pair<float, float> ret;
//     h_vnn->Divide(h_vnn_diag);
//     ret.first    = h_vnn->GetBinContent(1);
//     ret.second   = h_vnn->GetBinError  (1);
//     return ret;
//   }
// }
// //------------------------------------------------------------------------------
    int GetPtaIndexForPtbIndex(int iptb){
      return GetPtaIndex(PT2_LO[iptb],PT2_HI[iptb]);
    }

// //Returns a histogram storing the Fourier vn(pTb)
// TH1D* PtbdepHistVnFilled(std::vector<int> ptb_bins, std::string hist_name, int icent, int ipt1, int ich, int ideta, int ihar, TFile* InFile) {
//   TH1D* h_vn = PtbdepHist(ptb_bins, hist_name);
//   char name[600];
//   sprintf(name, "#it{v}_{%d}", ihar);
//   h_vn->GetYaxis()->SetTitle(name);
//   for (int iptb_bin = 0; iptb_bin < (int)(ptb_bins.size()); iptb_bin++) {
//     int ipt2 = ptb_bins[iptb_bin];
//     std::pair<float, float> vn = Bins::GetVnPtb(icent, ipt1, ipt2, ich, ideta, ihar, InFile);
//     h_vn->SetBinContent(iptb_bin + 1, vn.first);
//     h_vn->SetBinError  (iptb_bin + 1, vn.second);
//   }
//   return h_vn;
// }

// //Returns a histogram storing the Template vn(pTb)
// TH1D* PtbdepHistVnFilled(TFile* TemplateFile, std::vector<int> ptb_bins, std::string hist_name, int icent, int ipt1, int ich, int ideta, int ihar, int iperiph) {
//   TH1D* h_vn = PtbdepHist(ptb_bins, hist_name);
//   char name[600];
//   sprintf(name, "#it{v}_{%d}", ihar);
//   h_vn->GetYaxis()->SetTitle(name);
//   for (int iptb_bin = 0; iptb_bin < (int)(ptb_bins.size()); iptb_bin++) {
//     int ipt2 = ptb_bins[iptb_bin];
//     std::pair<float, float> vn = Bins::GetVnPtb(TemplateFile, icent, ipt1, ipt2, ich, ideta, ihar, iperiph);
//     h_vn->SetBinContent(iptb_bin + 1, vn.first);
//     h_vn->SetBinError  (iptb_bin + 1, vn.second);
//   }
//   return h_vn;
// }

// //Returns a histogram storing the Fourier vn(cent)
// TH1D* CentdepHistVnFilled(std::vector<int> centbins, std::string hist_name, int ipt1, int ipt2, int ich, int ideta, int ihar, TFile* InFile) {
//   TH1D* h_vn = Bins::CentdepHist(centbins, hist_name);
//   char name[600];
//   sprintf(name, "#it{v}_{%d}", ihar);
//   h_vn->GetYaxis()->SetTitle(name);
//   for (int icent_bin = 0; icent_bin < (int)(centbins.size()); icent_bin++) {
//     int icent = centbins[icent_bin];
//     std::pair<float, float> vn = Bins::GetVnPtb(icent, ipt1, ipt2, ich, ideta, ihar, InFile);
//     h_vn->SetBinContent(icent_bin + 1, vn.first);
//     h_vn->SetBinError  (icent_bin + 1, vn.second);
//   }
//   return h_vn;
// }

// //------------------------------------------------------------------------------
// //These functions are used to restruict the values of the different deta,pt,centrality...
// //bins that we loop over in the main code
// //For the defaults cases we should typically loop over all bins
// //but for most cross-checks only looping over a few bins is necessary

std::vector<int> CentBins() {
  std::vector<int> full_set_cent;
  std::vector<int> min_set_cent;
  for (int icent = 0; icent < NCENT + NCENT_ADD; icent++) full_set_cent.push_back(icent);
  return full_set_cent;
}
std::vector<int> TrkBins() {
  std::vector<int> full_set_cent;
  std::vector<int> min_set_cent;
  for (int itrk = 0; itrk < NTRK + NTRK_ADD; itrk++) full_set_cent.push_back(itrk);
  return full_set_cent;
}

std::vector<int> CentBinsPeriph() {
  std::vector<int> full_set_cent;
  full_set_cent.push_back(GetCentIndex(0.0, 0.68 ));
  //full_set_cent.push_back(GetCentIndex(0.0, 1.36 ));
  full_set_cent.push_back(GetCentIndex(0.0, 13.6 ));
  return full_set_cent;
}
std::vector<int> TrkBinsPeriph() {
  std::vector<int> full_set_trk;
  full_set_trk.push_back(GetTrkIndex(0, 10 ));
  full_set_trk.push_back(GetTrkIndex(0, 20 ));
  //full_set_trk.push_back(GetTrkIndex(0, 1000 ));
  return full_set_trk;
}


std::vector<int> PtaBins () {
  std::vector<int> full_set_pta;
  std::vector<int> min_set_pta;
  for (int ipt1 = 0; ipt1 < NPT1 + NPT1_ADD; ipt1++) full_set_pta.push_back(ipt1);
  //return specific bins
  // min_set_pta.push_back(GetPtaIndex(0.5, 5.0));
  // min_set_pta.push_back(GetPtaIndex(2.0, 5.0));


  // //Default case; return all pta bins
  // if ((l_data_type       == DataSetEnums::DATA_XeXe || l_data_type == DataSetEnums::DATA) &&
  //     l_apply_efficiency == DataSetEnums::DEFAULT_EFFICIENCY                              &&
  //     l_trig_type        == DataSetEnums::ALL_TRIGS                                       &&
  //     l_track_quality    == DataSetEnums::HI_TIGHT                                        &&
  //     l_nc               == 20                                                            &&
  //     l_nz               == 100                                                             ) return full_set_pta;

  //return min_set_pta;
  return full_set_pta;
}

std::vector<int> PtbBins () {
  std::vector<int> full_set_ptb;
  std::vector<int> min_set_ptb;
  for (int ipt2 = 0; ipt2 < NPT2 + NPT2_ADD; ipt2++) full_set_ptb.push_back(ipt2);
  //return specific bins
  // min_set_ptb.push_back(GetPtbIndex(0.5, 5.0));
  // min_set_ptb.push_back(GetPtbIndex(0.5, 1.0));
  // min_set_ptb.push_back(GetPtbIndex(1.0, 2.0));
  // min_set_ptb.push_back(GetPtbIndex(2.0, 3.0));
  // min_set_ptb.push_back(GetPtbIndex(3.0, 4.0));
  // min_set_ptb.push_back(GetPtbIndex(4.0, 6.0));
  // min_set_ptb.push_back(GetPtbIndex(6.0, 8.0));
  // min_set_ptb.push_back(GetPtbIndex(8.0, 10.0));
  // min_set_ptb.push_back(GetPtbIndex(10.0, 15.0));
  return full_set_ptb;
}

std::vector<int> DetaBins() {
  std::vector<int> min_set_deta = {GetDetaIndex(2.0, 5.0)};
  std::vector<int> full_set_deta;
  for (int ideta = 0; ideta < NDETA; ideta++) full_set_deta.push_back(ideta);

  return min_set_deta;
  //return full_set_deta;
}

std::vector<int> ChBins() {
  std::vector<int> min_set_ch = {COMBINED_CHARGE};
  std::vector<int> full_set_ch = {SAME_CHARGE, OPPOSITE_CHARGE, COMBINED_CHARGE};

  return min_set_ch;
  //return full_set_ch;
}
// bool PassAnalysisBin(int icent, int ipt1, int ipt2, int ich, int ideta) {
//   if (ich == COMBINED_CHARGE && ideta == GetDetaIndex(2.0, 5.0)) return true;
//   return false;
// }
// //------------------------------------------------------------------------------

int Initialize() {
  std::cout << "Called Bins::Initialize()" << std::endl;
  Initialize_Pt1Add();
  Initialize_Pt2Add();
  Initialize_CentAdd();
  Initialize_TrkAdd();
  return 1;
}
int m_is_initialized = Initialize();




}


#endif

#ifndef _Defs_HH_
#define _Defs_HH_


#include <vector>
#include "Riostream.h"
#include "TChain.h"
#include <string>
#include <fstream>

#include "common.C"
#include "bins.h"

using namespace std;

std::string path = "/gpfs0/citron/users/bargl/ZDC/lhcf22/user.steinber.data22_13p6TeV.00435229.physics_MinBias.merge.AOD.r14470_p5587.4zdc_EXT0"; //nutples


TChain *fChain        =nullptr;
TTree *tree;
int nmix;
TFile *tmpf;
char histtitle[100];
char histname[100];
char output_name[100];

const double PI = acos(-1.0);

TRandom *ran=new TRandom;
enum TWOPCTYPE {
    FG_HADRON_HADRON = 0,
    BG_HADRON_HADRON = 1,
};

// pooling information for bg distribution
enum {
  //pool bins:
  nc  = Bins::NCENT,
  n_multiplicity = Bins::NTRK,
  nz  =40,          
  ZMAX=200,
  plnbr = nc*nz*n_multiplicity,
};
vector< EVENT_PTR > pool[plnbr];
int dep =20; //each pool have size of 20


//histograms
TH1* h_trk_test; // for test
TH1* h_Trig;
TH1* h_Zvtx;
TH1* hNtrk[Bins::NCENT];
TH1* hNtrk_no_cut;
TH1* h_eff_no_ps;
TH1* heff; //for Eff distribution
TH1* heff_after_cut; //for Eff distribution
TH1* hzdc;
TH1* hzdc_A_with_pileup;
TH1* hzdc_C_with_pileup;
TH1* hzdc_A_without_pileup;
TH1* hzdc_C_without_pileup;
TH1* hzdc_A_without_pileup_noOfflineCut;
TH1* hzdc_C_without_pileup_noOfflineCut;
TH1* hzdc_after_cut;
TH1* N_trigger[Bins::NCENT][Bins::NTRK];
TH1* h_pt[Bins::NCENT];
TH1* h_eta[Bins::NCENT];
TH2* hLucrodCorr;
TH2* hAmpCorr;
TH2* hLucrodZdcCorr[2];
TH2* hNtrkEff;
TH2* hZdcCorr;
TH2* hZdcCorr_after_cut;
TH2* h_EtaPhi[Bins::NPT1];
TH2* fg [Bins::NCENT][Bins::NTRK][Bins::NPT1][Bins::NPT2][Bins::NCH];
TH2* bg [Bins::NCENT][Bins::NTRK][Bins::NPT1][Bins::NPT2][Bins::NCH];

//variables used in events analysis
float m_zvtx;
unsigned int    lumiBlock;
int             nvtx;
unsigned int     ntrk;
int             m_cent_i;
int             nbin;
double sqrt_s = 13600.0; //13.6 TeV
bool minbias;
int prescale;
int trig_index;
//int N_evt_sampled;
int good_events[Bins::NCENT]; //good events for normalizing pt and eta histograms
std::vector<std::string> m_TrigNames;
std::vector<std::string> m_EffEnergyNames;

//TODO define it as part of class

//ZDC weights
//vector<float> zdcWei ={0,2.78128,2.37211,3.31743,0,3.80388,3.26592,5.26164}; //EM on both sides set to 0 //OLD- BEFORE REPORECESSED!
//vector<float> zdcWei_a ={0,3.03828,3.47985,2.85721,0,3.8468,3.01163,3.14491}; //EM on both sides set to 0 
//std::vector<float> no_booster = {0.54, 1.00, 0.94, 0.79,1.47,1.02,0.87,0.54};
//std::vector<float> no_booster = {1.0/2.67, 1.0/1.64, 1.0/0.91,1.0/1.33,1.0/1.11,1.0/1.85,1.0/1.29,1.0/2.27}; // this is the HV gains from pb23
std::vector<float> no_booster = {1.0, 1.00, 1.0, 1.0,1.0, 1.00, 1.0, 1.0};
vector<float> zdcWei;

void load_weights();
void InitHistos();
bool passTrigger(std::vector<bool> trigger);
int triggerIndex(std::vector<bool> trigger);
bool isBitSet(int x, int s);
void SaveHistos();
int get_zPool(float z);
void Fill(Event* event1, Event* event2, int mixtype);
bool FillMixed(EVENT_PTR event);
bool isElectroMagnetic(int side, int mask);
std::vector<float> GetMoments(int ilb_bin, int side, int Trig);
#endif

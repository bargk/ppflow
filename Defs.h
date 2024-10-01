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

std::string path = "/gpfs0/citron/users/bargl/ZDC/lhcf22/user.steinber.data22_13p6TeV.00435229.physics_MinBias.merge.AOD.r14470_p5587.4zdc_EXT0"; //zdc
bool m_run_on_grid    =false;
TChain *fChain        =nullptr;
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


//global information
TH1* h_Zvtx;
TH1* hNtrk;
TH1* heff;
TH1* hzdc;
TH1* N_trigger[Bins::NCENT][Bins::NTRK];
TH1* h_pt[Bins::NTRK];
//TH1* N_ntrk[Bins::NCENT];
TH1* h_eta[Bins::NTRK];
TH2* hNtrkEff;
TH2* h_EtaPhi[Bins::NPT1];
TH2* fg [Bins::NCENT][Bins::NTRK][Bins::NPT1][Bins::NPT2][Bins::NCH];
TH2* bg [Bins::NCENT][Bins::NTRK][Bins::NPT1][Bins::NPT2][Bins::NCH];

float m_zvtx;
int m_nz              =20; // each zvtx bin has width 10 mm


unsigned int    lumiBlock;
int             nvtx;
unsigned int     ntrk;
int             m_cent_i;
int             nbin;
double sqrt_s = 13600.0; //13.6 TeV
bool minbias;

// std::vector<TTreeReaderValue<bool>> m_trig_minbias;
// std::vector<TTreeReaderValue<bool>> m_trig_zdc;

//TODO define it as part of class

//ZDC weights
vector<float> zdcWei ={0,2.78128,2.37211,3.31743,0,3.80388,3.26592,5.26164}; //EM on both sides set to 0

void InitHistos();
bool passTrigger(std::vector<bool> trigger);
bool isBitSet(int x, int s){
  int mask = x >> s;
  return mask % 2;
}



void SaveHistos(){
  tmpf->cd();
  tmpf->Write();
  cout << "Saving Finished" << endl;
}

int get_zPool(float z){
  int bin=-1;
  if(fabs(z)>=ZMAX) return -1;
  bin =(int)  (((z+ZMAX)/(2.0*ZMAX))*nz);  // pools
  if(bin<0||bin>=nz) return -1;
  return bin;
}

void Fill(Event* event1, Event* event2, int mixtype) {

    int icent = (m_cent_i);
    if (icent < 0 || icent >= Bins::NCENT) return;

    int itrk = (nbin);
    if (itrk < 0 || itrk >= Bins::NTRK) return;

    int ntrk1 = event1->get_npart();
    int ntrk2 = event2->get_npart();

    float phi01, eta1, qop1, eff1; //,pt1;
    float phi02, eta2, qop2, eff2; //,pt2;
    int   ptbin1, ptbin2;

    int sig1 = (int)(ran->Rndm() * 2);
    float wei = 1.0 / nmix;

    for (int i = 0; i < ntrk1; i++) {
        Track* trk1 = event1->GetTrack(i);
        eta1   = (trk1)->get_eta();
        phi01  = (trk1)->get_phi0();
        qop1   = (trk1)->get_charge();
        ptbin1 = (int)((trk1)->get_ptbin1());//trigger bin
        eff1   = (trk1)->get_eff();
        if (ptbin1 < 0) continue;

        int j = 0;
        while (j < ntrk2)  {
            if (mixtype == TWOPCTYPE::FG_HADRON_HADRON && j == i) {j++; continue;}
            Track* trk2 = event2->GetTrack(j);
            j++;
            eta2   = (trk2)->get_eta();
            phi02  = (trk2)->get_phi0();
            qop2   = (trk2)->get_charge();
            ptbin2 = (int) ((trk2)->get_ptbin2());
            eff2   = (trk2)->get_eff();
            if (ptbin2 < 0) continue;

            float dphi0 = phi01 - phi02;
            float deta  = fabs(eta1 - eta2);
            if (sig1) dphi0 = -dphi0;
            if      (dphi0 > 1.5 * Common::PI)  dphi0 -= 2 * Common::PI;
            else if (dphi0 < -0.5 * Common::PI) dphi0 += 2 * Common::PI;

            if (mixtype == TWOPCTYPE::FG_HADRON_HADRON) {
                if ((qop1 * qop2) > 0) fg[icent][itrk][ptbin1][ptbin2][0]->Fill(dphi0, deta, eff1 * eff2);
                else                   fg[icent][itrk][ptbin1][ptbin2][1]->Fill(dphi0, deta, eff1 * eff2);
            }
            else if (mixtype == TWOPCTYPE::BG_HADRON_HADRON) {
                if ((qop1 * qop2) > 0) bg[icent][itrk][ptbin1][ptbin2][0]->Fill(dphi0, deta, wei * eff1 * eff2);
                else                   bg[icent][itrk][ptbin1][ptbin2][1]->Fill(dphi0, deta, wei * eff1 * eff2);
            }
        }
    }
}



bool FillMixed(EVENT_PTR event) {
    int zbin    = get_zPool   (m_zvtx);
    int indx    = zbin + nbin*nz + m_cent_i * nz*n_multiplicity; //a * n_B * n_C + b * n_C + c
    //int dep     = depth[m_cent_i];

    if (zbin == -1) {
        std::cout << "AAAAAAAA " << zbin << "  " << m_zvtx << std::endl;
        delete event;
        return false;
    }

    nmix = pool[indx].size();
    if (nmix > 0 && nmix <= dep) {
        for (int k = 0; k < nmix; k++) {
            Fill(pool[indx][k], event    , TWOPCTYPE::BG_HADRON_HADRON); //;mixing
        }
    }


    if (nmix >= dep) {
        delete pool[indx].at(0);
        pool[indx].erase( (pool[indx]).begin() );
    }
    pool[indx].push_back( event );

    return true;
}




#endif

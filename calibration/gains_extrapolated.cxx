#include "calibration.h"
void gains_extrapolated(){



    //load weights for side C
    std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide");
    TVectorD *gains12;
    TVectorD *gains12_err;
    TVectorD *gains12_corrected;
    TVectorD *gains12_err_corrected;
    TFile *input = new TFile(Form("%s/zdcWeights_side0.root",base.c_str()));
    gains12 = (TVectorD*)input->Get("gains_avg");
    //gains12_err = (TVectorD*)input->Get("gains_std");
    gains12_corrected = (TVectorD*)input->Get("gains_avg");
    //gains12_err_corrected = (TVectorD*)input->Get("gains_std");

    //extapolate HAD2 side C
    (*gains12_corrected)[2]= ((had2_12->Eval(1500))/had2_12->Eval(1450))*((*gains12)[2]); //12HAD2
    // TODO add error propagation

    TFile *output = new TFile(Form("%s/zdcWeights_side0_corrected.root",base.c_str()),"RECREATE");
    output->cd();
    gains12_corrected->Write("gains");
    //gains12_err_corrected->Write("gains_std");
    output->Close();

    //print the gains HV for run 463315: f(nominal)/f(463315)

    cout << "HV gains for run 463315" << endl;
    cout << "EM C " << em12_hv << endl;
    cout << "HAD1 C " << had1_12_hv << endl;
    cout << "HAD2 C " << had2_12_hv << endl;
    cout << "HAD3 C " << had3_12_hv << endl;
    cout << "EM A " << em81_hv << endl;
    cout << "HAD1 A " << had1_81_hv << endl;
    cout << "HAD2 A " << had2_81_hv << endl;
    cout << "HAD3 A " << had3_81_hv << endl; 
}
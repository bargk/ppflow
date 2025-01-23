
void gains_extrapolated(){

    //pmts curves for HV response

    //Side 81
    //...............................81EM.......................................
    TF1 *em_81 = new TF1("f81EM", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    em_81->SetParameter(0,-9.94);
    em_81->SetParameter(1,1.38e-2);
    em_81->SetParameter(2,-2.32e-6);

    //...............................81HAD1.......................................
    TF1 *had1_81 = new TF1("f81HAD1", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    had1_81->SetParameter(0,-9.73);
    had1_81->SetParameter(1,1.33e-2);
    had1_81->SetParameter(2,-2.25e-6);

    TF1 *had1_81_err = new TF1("f81HAD1_err", "TMath::Sqrt((((2*[0]*x + [1])*f81HAD1)**2))",1200,1800);
    had1_81_err->SetParameter(0,-2.25e-6);
    had1_81_err->SetParameter(1,1.33e-2);
    // cout<< had1_81->Eval(1500) << endl;
    //had1_81_err->Draw();
    //...............................81HAD2.......................................
    TF1 *had2_81 = new TF1("f81HAD2", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    had2_81->SetParameter(0,-10.84);
    had2_81->SetParameter(1,1.44e-2);
    had2_81->SetParameter(2,-2.66e-6);

    TF1 *had2_81_err = new TF1("f81HAD2_err", "TMath::Sqrt((((2*[0]*x + [1])*f81HAD2)**2))",1200,1800);
    had2_81_err->SetParameter(0,-2.66e-6);
    had2_81_err->SetParameter(1,1.44e-2);
    //...............................81HAD3.......................................
    TF1 *had3_81 = new TF1("f81HAD3", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    had3_81->SetParameter(0,-10.44);
    had3_81->SetParameter(1,1.33e-2);
    had3_81->SetParameter(2,-2.3e-6);

    TF1 *had3_81_err = new TF1("f81HAD3_err", "TMath::Sqrt((((2*[0]*x + [1])*f81HAD3)**2))",1200,1800);
    had3_81_err->SetParameter(0,-2.3e-6);
    had3_81_err->SetParameter(1,1.33e-2);
    //............................................................................

    //Side 12
    //...............................12EM.......................................
    TF1 *em_12 = new TF1("f12EM", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    em_12->SetParameter(0,-9.56);
    em_12->SetParameter(1,1.25e-2);
    em_12->SetParameter(2,-2.05e-6);
    //...............................12HAD1.......................................
    TF1 *had1_12 = new TF1("f12HAD1", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    had1_12->SetParameter(0,-10.35);
    had1_12->SetParameter(1,1.39e-2);
    had1_12->SetParameter(2,-2.43e-6);

    TF1 *had1_12_err = new TF1("f12HAD1_err", "TMath::Sqrt((((2*[0]*x + [1])*f12HAD1)**2))",1200,1800);
    had1_12_err->SetParameter(0,-2.43e-6);
    had1_12_err->SetParameter(1,1.39e-2);

    //...............................12HAD2.......................................
    TF1 *had2_12 = new TF1("f12HAD2", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    had2_12->SetParameter(0,-9.69);
    had2_12->SetParameter(1,1.32e-2);
    had2_12->SetParameter(2,-2.26e-6);

    TF1 *had2_12_err = new TF1("f12HAD2_err", "TMath::Sqrt((((2*[0]*x + [1])*f12HAD2)**2))",1200,1800);
    had2_12_err->SetParameter(0,-2.26e-6);
    had2_12_err->SetParameter(1,1.32e-2);

    //...............................12HAD3.......................................
    TF1 *had3_12 = new TF1("f12HAD3", "TMath::Exp([0]+[1]*x+[2]*x**2)",1200,1800);
    had3_12->SetParameter(0,-11.02);
    had3_12->SetParameter(1,1.49e-2);
    had3_12->SetParameter(2,-2.74e-6);

    TF1 *had3_12_err = new TF1("f12HAD3_err", "TMath::Sqrt((((2*[0]*x + [1])*f12HAD3)**2))",1200,1800);
    had3_12_err->SetParameter(0,-2.74e-6);
    had3_12_err->SetParameter(1,1.49e-2);

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
    float em12_hv = em_12->Eval(1525.5)/em_12->Eval(1375.5);
    float had1_12_hv = had1_12->Eval(1446.0)/had1_12->Eval(1375.0);
    float had2_12_hv = had2_12->Eval(1486.0)/had2_12->Eval(1500.0);
    float had3_12_hv = had3_12->Eval(1544.0)/had3_12->Eval(1500.0);

    float em81_hv = em_81->Eval(1390.0)/em_81->Eval(1375.0);
    float had1_81_hv = had1_81->Eval(1463.0)/had1_81->Eval(1374.0);
    float had2_81_hv = had2_81->Eval(1541.0)/had2_81->Eval(1500.0);
    float had3_81_hv = had3_81->Eval(1635.0)/had3_81->Eval(1500.0);

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
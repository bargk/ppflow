void gains_extrapolated(){
    std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide");
    //pmts curves for HV response

    // Side 81
    //...............................81EM.......................................
    TF1 *em_81 = new TF1("f81EM", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    em_81->SetParameter(0, -9.94);
    em_81->SetParameter(1, 1.38e-2);
    em_81->SetParameter(2, -2.32e-6);

    TF1 *em_81_derivative = new TF1("f81EM_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    em_81_derivative->SetParameter(0, -9.94);
    em_81_derivative->SetParameter(1, 1.38e-2);
    em_81_derivative->SetParameter(2, -2.32e-6);

    //...............................81HAD1.......................................
    TF1 *had1_81 = new TF1("f81HAD1", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had1_81->SetParameter(0, -9.73);
    had1_81->SetParameter(1, 1.33e-2);
    had1_81->SetParameter(2, -2.25e-6);

    TF1 *had1_81_derivative = new TF1("f81HAD1_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had1_81_derivative->SetParameter(0, -9.73);
    had1_81_derivative->SetParameter(1, 1.33e-2);
    had1_81_derivative->SetParameter(2, -2.25e-6);

    //...............................81HAD2.......................................
    TF1 *had2_81 = new TF1("f81HAD2", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had2_81->SetParameter(0, -10.84);
    had2_81->SetParameter(1, 1.44e-2);
    had2_81->SetParameter(2, -2.66e-6);

    TF1 *had2_81_derivative = new TF1("f81HAD2_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had2_81_derivative->SetParameter(0, -10.84);
    had2_81_derivative->SetParameter(1, 1.44e-2);
    had2_81_derivative->SetParameter(2, -2.66e-6);

    //...............................81HAD3.......................................
    TF1 *had3_81 = new TF1("f81HAD3", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had3_81->SetParameter(0, -10.44);
    had3_81->SetParameter(1, 1.33e-2);
    had3_81->SetParameter(2, -2.3e-6);

    TF1 *had3_81_derivative = new TF1("f81HAD3_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had3_81_derivative->SetParameter(0, -10.44);
    had3_81_derivative->SetParameter(1, 1.33e-2);
    had3_81_derivative->SetParameter(2, -2.3e-6);

    //............................................................................

    // Side 12
    //...............................12EM.......................................
    TF1 *em_12 = new TF1("f12EM", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    em_12->SetParameter(0, -9.56);
    em_12->SetParameter(1, 1.25e-2);
    em_12->SetParameter(2, -2.05e-6);

    TF1 *em_12_derivative = new TF1("f12EM_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    em_12_derivative->SetParameter(0, -9.56);
    em_12_derivative->SetParameter(1, 1.25e-2);
    em_12_derivative->SetParameter(2, -2.05e-6);

    //...............................12HAD1.......................................
    TF1 *had1_12 = new TF1("f12HAD1", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had1_12->SetParameter(0, -10.35);
    had1_12->SetParameter(1, 1.39e-2);
    had1_12->SetParameter(2, -2.43e-6);

    TF1 *had1_12_derivative = new TF1("f12HAD1_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had1_12_derivative->SetParameter(0, -10.35);
    had1_12_derivative->SetParameter(1, 1.39e-2);
    had1_12_derivative->SetParameter(2, -2.43e-6);

    //...............................12HAD2.......................................
    TF1 *had2_12 = new TF1("f12HAD2", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had2_12->SetParameter(0, -9.69);
    had2_12->SetParameter(1, 1.32e-2);
    had2_12->SetParameter(2, -2.26e-6);

    TF1 *had2_12_derivative = new TF1("f12HAD2_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had2_12_derivative->SetParameter(0, -9.69);
    had2_12_derivative->SetParameter(1, 1.32e-2);
    had2_12_derivative->SetParameter(2, -2.26e-6);

    //...............................12HAD3.......................................
    TF1 *had3_12 = new TF1("f12HAD3", "TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had3_12->SetParameter(0, -11.02);
    had3_12->SetParameter(1, 1.49e-2);
    had3_12->SetParameter(2, -2.74e-6);

    TF1 *had3_12_derivative = new TF1("f12HAD3_derivative", "(2*[2]*x + [1]) * TMath::Exp([0]+[1]*x+[2]*x**2)", 1200, 1800);
    had3_12_derivative->SetParameter(0, -11.02);
    had3_12_derivative->SetParameter(1, 1.49e-2);
    had3_12_derivative->SetParameter(2, -2.74e-6);




    //load weights for side C
    TVectorD *gains12;
    TVectorD *gains12_err;
    TVectorD *gains12_corrected;
    TVectorD *gains12_err_corrected;
    TFile *input = new TFile(Form("%s/zdcWeights_side0.root",base.c_str()));
    gains12 = (TVectorD*)input->Get("gains_avg");
    gains12_err = (TVectorD*)input->Get("gains_std");
    gains12_corrected = (TVectorD*)input->Get("gains_avg");
    gains12_err_corrected = (TVectorD*)input->Get("gains_std");

    //extapolate HAD1C
    (*gains12_corrected)[1]= ((had1_12->Eval(1375))/had1_12->Eval(1350))*((*gains12)[1]); 

    double rel_wpb1_c = (*gains12_err)[1]/(*gains12)[1];
    double rel_vpb_had1_c = had1_12_derivative->Eval(1375)/had1_12->Eval(1375);  //assume 1volt uncertainty 
    double rel_vlhcf_had1_c = had1_12_derivative->Eval(1350)/had1_12->Eval(1350); //assume 1volt uncertainty 
    (*gains12_err_corrected)[1] = (*gains12_corrected)[1]*(TMath::Sqrt(rel_wpb1_c*rel_wpb1_c + rel_vpb_had1_c*rel_vpb_had1_c + rel_vlhcf_had1_c*rel_vlhcf_had1_c));
    TFile *output = new TFile(Form("%s/zdcWeights_side0_corrected.root",base.c_str()),"RECREATE");
    output->cd();
    gains12_corrected->Write("gains");
    gains12_err_corrected->Write("gains_std");
    output->Close();

    // //load weights for side A
    TVectorD *gains81;
    TVectorD *gains81_err;
    TVectorD *gains81_corrected;
    TVectorD *gains81_err_corrected;
    input = new TFile(Form("%s/zdcWeights_side1.root",base.c_str()));
    gains81 = (TVectorD*)input->Get("gains_avg");
    gains81_err = (TVectorD*)input->Get("gains_std");
    gains81_corrected = (TVectorD*)input->Get("gains_avg");
    gains81_err_corrected = (TVectorD*)input->Get("gains_std");

    //extapolate HAD1A
    (*gains81_corrected)[1]= ((had1_81->Eval(1375))/had1_81->Eval(1350))*((*gains81)[1]); 
    double rel_wpb1_a = (*gains81_err)[1]/(*gains81)[1];
    double rel_vpb_had1_a = had1_81_derivative->Eval(1375)/had1_81->Eval(1375);  //assume 1volt uncertainty 
    double rel_vlhcf_had1_a = had1_81_derivative->Eval(1350)/had1_81->Eval(1350); //assume 1volt uncertainty 
    (*gains81_err_corrected)[1] = (*gains81_corrected)[1]*(TMath::Sqrt(rel_wpb1_a*rel_wpb1_a + rel_vpb_had1_a*rel_vpb_had1_a + rel_vlhcf_had1_a*rel_vlhcf_had1_a));

    //extapolate HAD2A
    (*gains81_corrected)[2]= ((had2_81->Eval(1500))/had2_81->Eval(1450))*((*gains81)[2]); 
    double rel_wpb2_a = (*gains81_err)[2]/(*gains81)[2];
    double rel_vpb_had2_a = had2_81_derivative->Eval(1500)/had2_81->Eval(1500);  //assume 1volt uncertainty 
    double rel_vlhcf_had2_a = had2_81_derivative->Eval(1450)/had2_81->Eval(1450); //assume 1volt uncertainty 
    (*gains81_err_corrected)[2] = (*gains81_corrected)[2]*(TMath::Sqrt(rel_wpb2_a*rel_wpb2_a + rel_vpb_had2_a*rel_vpb_had2_a + rel_vlhcf_had2_a*rel_vlhcf_had2_a));
    TFile *output1 = new TFile(Form("%s/zdcWeights_side1_corrected.root",base.c_str()),"RECREATE");
    output1->cd();
    gains81_corrected->Write("gains");
    gains81_err_corrected->Write("gains_std");
    output1->Close();

    // //print the gains HV for run 463315: f(nominal)/f(463315)
    // float em12_hv = em_12->Eval(1525.5)/em_12->Eval(1375.5);
    // float had1_12_hv = had1_12->Eval(1446.0)/had1_12->Eval(1375.0);
    // float had2_12_hv = had2_12->Eval(1486.0)/had2_12->Eval(1500.0);
    // float had3_12_hv = had3_12->Eval(1544.0)/had3_12->Eval(1500.0);

    // float em81_hv = em_81->Eval(1390.0)/em_81->Eval(1375.0);
    // float had1_81_hv = had1_81->Eval(1463.0)/had1_81->Eval(1374.0);
    // float had2_81_hv = had2_81->Eval(1541.0)/had2_81->Eval(1500.0);
    // float had3_81_hv = had3_81->Eval(1635.0)/had3_81->Eval(1500.0);

    // cout << "HV gains for run 463315" << endl;
    // cout << "EM C " << em12_hv << endl;
    // cout << "HAD1 C " << had1_12_hv << endl;
    // cout << "HAD2 C " << had2_12_hv << endl;
    // cout << "HAD3 C " << had3_12_hv << endl;
    // cout << "EM A " << em81_hv << endl;
    // cout << "HAD1 A " << had1_81_hv << endl;
    // cout << "HAD2 A " << had2_81_hv << endl;
    // cout << "HAD3 A " << had3_81_hv << endl; 
}




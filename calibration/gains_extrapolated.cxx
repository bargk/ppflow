void gains_extrapolated(){
    //Side 81
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
    std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles");
    TVectorD *gains12;
    TVectorD *gains12_err;
    TVectorD *gains12_corrected;
    TVectorD *gains12_err_corrected;
    TFile *input = new TFile(Form("%s/2.0sigma/zdcWeights_side0.root",base.c_str()));
    gains12 = (TVectorD*)input->Get("gains_avg");
    gains12_err = (TVectorD*)input->Get("gains_std");
    gains12_corrected = (TVectorD*)input->Get("gains_avg");
    gains12_err_corrected = (TVectorD*)input->Get("gains_std");

    //extapolate HAD2 side C
    (*gains12_corrected)[2]= ((had2_12->Eval(1500))/had2_12->Eval(1450))*((*gains12)[2]); //12HAD2
    // TODO add error propagation

    TFile *output = new TFile(Form("%s/2.0sigma/zdcWeights_side0_corrected.root",base.c_str()),"RECREATE");
    output->cd();
    gains12_corrected->Write("gains_avg");
    gains12_err_corrected->Write("gains_std");
    output->Close();

    // //load gain factors from TB2023 and create new root file with the new weights
    // std::vector<int> mods ={12,81};
    // std::vector<string> names ={"hamamatsu", "no_booster", "booster", "no_pmt"};
    // TVectorT<double> *zdcWeights = new TVectorT<double>(3);
    // TVectorT<double> *zdcWeights_err = new TVectorT<double>(3);
    // TFile *file = new TFile("/afs/cern.ch/user/b/bglik/ZDC/lhcf22/calibration/zdcWeights_V2.root", "RECREATE");
    // for(const auto& mod:mods){
    //     for(const auto& name:names){
    //         TFile *fopen = new TFile(Form("/afs/cern.ch/user/b/bglik/ZDC/TB2023/no_lucrod/Version2/LG_method/1.5sigma/%i/%s/gains_avg.root",mod,name.c_str()),"READ");
    //         TVectorT<double>* gains = dynamic_cast<TVectorT<double>*>(fopen->Get("gains_avg"));
    //         TVectorT<double>* gains_err = dynamic_cast<TVectorT<double>*>(fopen->Get("gains_std"));
    //         //cout<< (*gains)[0] << endl;
    //         if(mod ==12){
    //             (*zdcWeights)[0]= ((had1_12->Eval(1500))/had1_12->Eval(1375))*((*gains)[1]); //12HAD1
    //             double rel_err_TB_had1 = (had1_12_err->Eval(1500))/(had1_12->Eval(1500)); //relative error 
    //             double rel_err_pp_had1 = (had1_12_err->Eval(1375))/(had1_12->Eval(1375));
    //             (*zdcWeights_err)[0] = TMath::Sqrt(std::pow(rel_err_pp_had1,2) + std::pow(rel_err_TB_had1,2) + std::pow((*gains_err)[1],2));
            
    //             (*zdcWeights)[1]= ((had2_12->Eval(1500))/had2_12->Eval(1450))*((*gains)[2]); //12HAD2
    //             double rel_err_TB_had2 = (had2_12_err->Eval(1500))/(had2_12->Eval(1500)); //relative error 
    //             double rel_err_pp_had2 = (had2_12_err->Eval(1450))/(had2_12->Eval(1450));
    //             (*zdcWeights_err)[1] = TMath::Sqrt(std::pow(rel_err_pp_had2,2)+ std::pow(rel_err_TB_had2,2) + std::pow((*gains_err)[2],2));

    //             (*zdcWeights)[2] = (*gains)[3];
    //             (*zdcWeights_err)[2] = (*gains_err)[3];        
    //         }
    //         else{
    //             (*zdcWeights)[0]= ((had1_81->Eval(1500))/had1_81->Eval(1375))*((*gains)[1]); //81HAD1
    //             double rel_err_TB_had1 = (had1_81_err->Eval(1500))/(had1_81->Eval(1500)); //relative error 
    //             double rel_err_pp_had1 = (had1_81_err->Eval(1375))/(had1_81->Eval(1375));
    //             (*zdcWeights_err)[0] = TMath::Sqrt(std::pow(rel_err_pp_had1,2) + std::pow(rel_err_TB_had1,2) + std::pow((*gains_err)[1],2));

    //             for(int i=1; i<3; i++){
    //                 (*zdcWeights)[i] = (*gains)[i+1];
    //                 (*zdcWeights_err)[i] = (*gains_err)[i+1];
    //             }
    //         }
    //     file->cd(); //writing to the zdcweitghts file
    //     zdcWeights->Write(Form("zdcWeights%i_%s", mod, name.c_str()));
    //     zdcWeights_err->Write(Form("zdcWeights%i_err_%s", mod, name.c_str()));
    //     }
    // }
    // file->Close();
}
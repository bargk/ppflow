
std::string path = "/gpfs0/citron/users/bargl/ZDC/user.bglik.data23_hi.00463315.calibration_ZDCCalib.merge.AOD.c1535_m2248.ANALYSIS_EXT0/";
std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide");

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


//histograms
TH1D* h0[2]; 
TH1D* h1[2]; 
TH1D* h_cut[2]; 
TH1D* h_module[8]; 
TH2D* h_module_corr[4];
TH2D* h0_corr; 
TH2D* h1_corr;
TFile *output; 

//fit parameters
float chi2_arr[2];
float ndf_arr[2];
float constant_arr[2];
float mean_arr[2]; 
float sigma_arr[2];
float constErr_arr[2];
float meanErr_arr[2];
float sigmaErr_arr[2];


// std::vector<float> no_booster = {0.54, 1.00, 0.94, 0.79,1.47,1.02,0.87,0.54}; //first 4 lower bits is side C 
std::vector<float> no_booster = {1.00, 1.00, 1.00, 1.00,1.00,1.00,1.00,1.00}; // during pbpb23 the nominal HV took into account the PMT relavtive sensitivty
float em12_hv = em_12->Eval(1525.5)/em_12->Eval(1375.5);
float had1_12_hv = had1_12->Eval(1446.0)/had1_12->Eval(1375.0);
float had2_12_hv = had2_12->Eval(1486.0)/had2_12->Eval(1500.0);
float had3_12_hv = had3_12->Eval(1544.0)/had3_12->Eval(1500.0);

float em81_hv = em_81->Eval(1390.0)/em_81->Eval(1375.0);
float had1_81_hv = had1_81->Eval(1463.0)/had1_81->Eval(1374.0);
float had2_81_hv = had2_81->Eval(1541.0)/had2_81->Eval(1500.0);
float had3_81_hv = had3_81->Eval(1635.0)/had3_81->Eval(1500.0);

std::vector<float> hv_gain = {2.67,1.64,0.91,1.33,1.11,1.85,1.29,2.27}; //HV gains to correct the PMT's voltage response for run 463315


// std::vector<float> no_booster = {0.54, 1.00, 0.94, 0.79,1.47,1.02,0.87,0.54}; //first 4 lower bits is side C 
std::vector<float> no_booster = {1.00, 1.00, 1.00, 1.00,1.00,1.00,1.00,1.00}; // during pbpb23 the nominal HV took into account the PMT relavtive sensitivty
std::vector<float> hv_gain = {2.67,1.64,0.91,1.33,1.11,1.85,1.29,2.27}; //HV gains to correct the PMT's voltage response for run 463315

std::string path = "/gpfs0/citron/users/bargl/ZDC/user.bglik.data23_hi.00463315.calibration_ZDCCalib.merge.AOD.c1535_m2248.ANALYSIS_EXT0/";
std::string base = Form("/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles/sameSide");



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



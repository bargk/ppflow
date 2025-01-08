
TH1D *sum[2];
float chi2_arr[2] = {0.0, 0.0};
float ndf_arr[2] = {0.0, 0.0};
float constant1_arr[2] = {0.0, 0.0};
float constant2_arr[2] = {0.0, 0.0};
float constant3_arr[2] = {0.0, 0.0};
float mean1_arr[2] = {0.0, 0.0};
float sigma1_arr[2] = {0.0, 0.0};
float sigma2_arr[2] = {0.0, 0.0};
float sigma3_arr[2] = {0.0, 0.0};

void plots_sumGaussians(){
    std::string base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/sumGaussiansPlots";
    std::string data = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/calibration/RootFiles";

    gSystem->Exec(Form("mkdir -p %s",base.c_str()));
    TFile *input1 = new TFile(Form("%s/sameSide/1.5sigma/zdc_calibrated.root",data.c_str()),"READ");
    sum[0] = (TH1D*)input1->Get("h_c");
    sum[1] = (TH1D*)input1->Get("h_a");
    sum[0]->GetXaxis()->SetRangeUser(500,8000);
    sum[1]->GetXaxis()->SetRangeUser(500,8000);
    int side = 0;
    TF1* gFit = new TF1("gFit","gaus");
        gFit->SetRange(2400,2900);
        sum[side]->Fit("gFit","Rqn");



        TF1 *fit = new TF1("fit", "[0]*TMath::Gaus(x, [1], [2]) + [3]*[0]*TMath::Gaus(x, 2*[1], [4]) + [5]*[0]*TMath::Gaus(x, 3*[1], [6])");
        
        //fit->SetParameters(6500,hist_mean,400,0.3,900,0.1,1500);
        fit->SetParameter(0,gFit->GetParameter(0));
        fit->SetParameter(1,gFit->GetParameter(1));
        fit->SetParameter(2,gFit->GetParameter(2));
        fit->SetParameter(3,0.2*gFit->GetParameter(0));
        fit->SetParameter(4,1.4*gFit->GetParameter(2));
        fit->SetParameter(5,0.04*gFit->GetParameter(0));
        fit->SetParameter(6,2*gFit->GetParameter(2));
        //fit->SetParameter(7,0.035*gFit->GetParameter(0));
        //fit->SetParameter(8,2*gFit->GetParameter(2));
        //fit->SetParameter(7,1e-5);
        //fit->SetParameter(8,1e-5);

        double low_lim = 0.85*gFit->GetParameter(1);
        double high_lim =2.7*gFit->GetParameter(1);
        fit->SetRange(low_lim,high_lim);
        //fit->SetParLimits(0,5e5,7e5);
        //fit->SetParLimits(1,hist_mean-5,hist_mean+5);
        fit->SetParLimits(2,0,1000);
        fit->SetParLimits(3,0,0.7);
        fit->SetParLimits(4,0,1000);
        fit->SetParLimits(5,0,0.05);
        fit->SetParLimits(6,0,3000);
        //fit->SetParLimits(7,1e-8,1e2);
        //fit->SetParLimits(8,1e-6,1e2);

        fit->SetLineColor(kRed);
        sum[side]->Fit("fit","R");
        chi2_arr[side] = fit->GetChisquare();
        ndf_arr[side] = fit->GetNDF();
        constant1_arr[side] = fit->GetParameter(0);
        constant2_arr[side] = fit->GetParameter(3);
        constant3_arr[side] = fit->GetParameter(5);
        mean1_arr[side] = fit->GetParameter(1);
        sigma1_arr[side] = fit->GetParameter(2);
        sigma2_arr[side] = fit->GetParameter(4);
        sigma3_arr[side] = fit->GetParameter(6);
        sum[side]->Draw();


        // //fit for 1n peak
        // TF1 *fit_1n = new TF1("fit_1n","gaus");
        // fit_1n->SetLineColor(kRed);
        // fit_1n->SetLineWidth(5);
        // fit_1n->SetRange(2400,3100);
        // sum[side]->Fit("fit_1n","R");
        // sum[side]->Draw();


}
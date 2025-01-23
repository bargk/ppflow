#include "calibration.h"

TFile *fopen2;
TVectorD* gains_old;




void Lagrange_method(int itr){
    int m_samples = 3; // how many sub "subsets"
    std::cout << "Dividing data into " << m_samples << " different samples!" << std::endl;
    TTree *tree = new TTree("gains","gains");
    std::vector<float> gains;
    int be_energy = 2680; // 1n peak energy [GeV]

    char name[100];
    char name1[100];
    char name2[100];
    sprintf(name,"%s/matrix_elements_itr%d.root",base.c_str(),itr);
    TChain *mychain_c = new TChain("mytree_c");
    TChain *mychain_a = new TChain("mytree_a");
    mychain_c->Add(name,0);
    mychain_a->Add(name,0);
    if(!mychain_a || !mychain_c) {std::cerr << "Error with loading data" << endl; exit(-1);}
    for (int side =0; side<2; side++){
        sprintf(name1,"%s/zdcWeights_side%i_itr%d.root",base.c_str(),side, itr);
        TTreeReader myReader;
        if(side ==0) {myReader.SetTree(mychain_c);}
        else { myReader.SetTree(mychain_a);}

        TTreeReaderValue<float> M0_branch(myReader,"M0");
        TTreeReaderValue<float> M1_branch(myReader,"M1");
        TTreeReaderValue<float> M2_branch(myReader,"M2");
        TTreeReaderValue<float> M3_branch(myReader,"M3");
       
        Long64_t counter = myReader.GetEntries();
        
        
        while(counter%m_samples) counter -= 1;
        Long64_t nevt = counter/m_samples;  

        // vectors for taking the averages
        std::vector<std::vector<double>> amp; 
        std::vector<std::vector<double>> amp2;
        std::vector<std::vector<double>> m0mj;
        std::vector<std::vector<double>> m1mj;
        std::vector<std::vector<double>> m2mj;
        std::vector<std::vector<double>> g;

        for(int i=0; i<m_samples; i++){

            //Initialize matrix elements (remember to add 0 to the element 5_5)
            float M0 = 0;
            float M0M0 = 0;
            float M0M1 = 0;
            float M0M2 = 0;
            float M0M3 = 0;
            float M1 = 0;
            float M1M1 = 0;
            float M1M2 = 0;
            float M1M3 = 0;
            float M2 = 0;
            float M2M2 = 0;
            float M2M3 = 0;
            float M3 = 0;
            float M3M3 = 0;
            double MaXADC;
            int trial =0;
            //Event loop
            myReader.SetEntriesRange(i*nevt, (i+1)*nevt); // Set range of event loop
            while(myReader.Next()){            
                M0 += *M0_branch;
                M0M0 += (*M0_branch) * (*M0_branch);
                M0M1 += (*M0_branch) * (*M1_branch);
                M0M2 += (*M0_branch) * (*M2_branch);
                M0M3 += (*M0_branch) * (*M3_branch);

                M1 += *M1_branch;
                M1M1 += (*M1_branch)* (*M1_branch);
                M1M2 += (*M1_branch)* (*M2_branch);
                M1M3 += (*M1_branch)* (*M3_branch);

                M2 += *M2_branch;
                M2M2 += (*M2_branch)* (*M2_branch);
                M2M3 += (*M2_branch)* (*M3_branch);

                M3 += *M3_branch;
                M3M3 += (*M3_branch)* (*M3_branch);
            }

            std::vector<double> temp_amp = {M0/nevt,M1/nevt, M2/nevt,M3/nevt,0};
            std::vector<double> temp_amp2 = {M0M0/nevt,M1M1/nevt, M2M2/nevt,M3M3/nevt,0};
            std::vector<double> temp_m0mj = {M0M1/nevt,M0M2/nevt, M0M3/nevt,0,0};
            std::vector<double> temp_m1mj = {M1M2/nevt,M1M3/nevt, 0,0,0};
            std::vector<double> temp_m2mj = {M2M3/nevt,0,0,0,0};

            std::vector<std::vector<double>> a = {{2*M0M0,2*M0M1,2*M0M2,2*M0M3,M0},   //Matrix elements (its symmetric)
                                                {2*M0M1,2*M1M1,2*M1M2,2*M1M3,M1},
                                                {2*M0M2,2*M1M2,2*M2M2,2*M2M3,M2},
                                                {2*M0M3,2*M1M3,2*M2M3,2*M3M3,M3},
                                                {M0,M1,M2,M3,0}};  
            for (auto& row : a) {         // Getting average of each module : <MiMj>
                for (auto& elem : row) { 
                    elem /= nevt;                     
                    }
                }
            TMatrixD A(5, 5);               //A(rows,colums)
            for (Int_t i = 0; i < 5; i++) {
                for (Int_t j = 0; j < 5; j++) {
                    A(i, j) = a[i][j];
                    }
                }        
            Double_t b[5] = {0,0,0,0,double(be_energy)}; //Initalize the constants vector
            TVectorD B(5,b);
            // Solve the system of equations using the TMatrixD::Solve method
            TDecompLU lu(A);
            Bool_t ok;
            TVectorD X = lu.Solve(B,ok);
            cout << "ok : " << ok << endl;
            cout << "gain factors g0,g1,g2,g3, lambda for side " << side <<" : ";
            std::vector<double> temp_g;
            for (Int_t i = 0; i < X.GetNrows() ; i++) {
                cout << X(i) << " ";
                if(X(i)<=0 && i<4){
                    std::cerr << "Error : calibration factor cannot be negative! " <<endl;
                    exit(-1);
                }
                temp_g.push_back(X(i));
                
            }
            cout << endl;
            cout << "Total events : " << nevt << endl;

            amp.push_back(temp_amp);
            amp2.push_back(temp_amp2);
            m0mj.push_back(temp_m0mj);
            m1mj.push_back(temp_m1mj);
            m2mj.push_back(temp_m2mj);
            g.push_back(temp_g);



            // Getting gain factors into root
            TFile *file_0 = new TFile(Form("%s/gains_side%i_%i.root",base.c_str(),side,i),"RECREATE");
            X.Write("gains"); //g0,g1,g2,g3, lambda
            file_0->Write();
            file_0->Close();    
        }
        // calculate avg constants

        std::vector<double> amp_avg; //<M>
        std::vector<double> amp2_avg; // <M^2>
        std::vector<double> m0mj_avg; //<M0Mj>
        std::vector<double> m1mj_avg; //<M1Mj>
        std::vector<double> m2mj_avg; //<M2Mj>
        std::vector<double> g_avg; // average g from m_samples
        for(int element=0; element<5; element++){
            double sum0 =0; double sum1 =0; double sum2 =0; double sum3 =0; double sum4 =0; double sum5 =0;
            for(int itr_idx =0; itr_idx<m_samples; itr_idx++){
                sum0 += amp[itr_idx].at(element);
                sum1 += amp2[itr_idx].at(element);
                sum2 += m0mj[itr_idx].at(element);
                sum3 += m1mj[itr_idx].at(element);
                sum4 += m2mj[itr_idx].at(element);
                sum5+= g[itr_idx].at(element);
            }

            amp_avg.push_back(sum0/m_samples);
            amp2_avg.push_back(sum1/m_samples);
            m0mj_avg.push_back(sum2/m_samples);
            m1mj_avg.push_back(sum3/m_samples);
            m2mj_avg.push_back(sum4/m_samples);
            g_avg.push_back(sum5/m_samples);
        }

        std::vector<TH1F*> histograms;
        for(int i=0; i<4; i++){
            TH1F *h0 = new TH1F(Form("h%i%i",side,i),"",10,0,10);
            histograms.push_back(h0);
        }

        for(int i=0; i<m_samples; i++){
            TFile* fout = new TFile(Form("%s/gains_side%i_%i.root",base.c_str(),side,i), "READ");
            TVectorT<double>* gain = dynamic_cast<TVectorT<double>*>(fout->Get("gains"));
            for(int j=0; j<4; j++){
                 histograms.at(j)->Fill((*gain)[j]);
             }
            fout->Close(); 
            gSystem->Exec(Form("rm %s/gains_side%i_%i.root",base.c_str(),side,i)); //remove the ith file of gains. 
        }

        TFile *fopen = new TFile(name1, "RECREATE");
        TVectorD gains =  TVectorD(4);
        TVectorD gains_err = TVectorD(4);
        for(int i=0; i<4; i++){
            gains[i] = histograms.at(i)->GetMean();
            gains_err[i] = histograms.at(i)->GetStdDev();
        }

        //before wrting the weights from current iteration, multiply it with weights from previous iteration
        int prev_itr = itr -1;
        //if don't have previous lagrange weights
        if(prev_itr == 0 && side == 0){
            for(int i=0; i<4; i++){
                gains[i] = gains[i] * hv_gain.at(i)/ no_booster.at(i);
            }
        }
        else if(prev_itr == 0 && side == 1){
            for(int i=0; i<4; i++){
                gains[i] = gains[i] * hv_gain.at(i+4)/ no_booster.at(i+4);
            }
        }

        else{
            TFile* file_weight = new TFile(Form("%s/zdcWeights_side%i_itr%i.root", base.c_str(), side,prev_itr), "READ");
            TVectorD* weights_prev_itr = (TVectorD*)file_weight->Get("gains_avg");
            for(int i=0; i<4; i++){
               gains[i] = gains[i] * (*weights_prev_itr)[i];
            }
        }

        std::cout << "Final weights after multiplication : " << std::endl;
        gains.Print();

        fopen->cd();
        gains.Write("gains_avg");
        gains_err.Write("gains_std");
        fopen->Write();
        fopen->Close();
        for(int i=0; i<4; i++){
            histograms.at(i)->Reset();
        }       
    }
    char letter[100];
    std::cout << "finished" << std::endl;
    sprintf(letter,"kill -9 %d",gSystem->GetPid());
    std::cout<<letter<<std::endl;
    gSystem->Exec(letter);

}


#include "/gpfs0/citron/users/bargl/ZDC/HI23/AtlasStyle.h"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/common.C"
#include "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/bins.h"

std::string  base = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/Rootfiles/sameSide";
//std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/internal-note/fig_pool/ana/systematics/vn";
std::string  figures = "/gpfs0/citron/users/bargl/ZDC/lhcf22/thesis/fig_pool/ana/systematics/vn";

int pericent = 22;
int peritrk = 5;
int itrk =4;
int ipt1 = 5;
int ipt2 = 5;

void plot_ana_systematics(){
    float X , Y;
    int size = 17;
    SetAtlasStyle();
    gSystem->Exec(Form("mkdir -p %s",figures.c_str()));
    TCanvas *c0 = new TCanvas("c0");
    TLegend *legend0;
    TFile *input; 
    THStack *hs = new THStack("hs",";E_{Eff} [TeV];v_{2}");
    THStack *hs_ratio = new THStack("hs",";E_{Eff} [TeV];v_{2} ratio");
    std::vector<int> peri_cent = {22,20,21};
    std::vector<int> peri_trk = {5,0,1};
    std::vector<int> style = {20,24,21};
    TH1D *h_zdc[4];
    TH1D *h_denominator;
    TH1D *h_ratio_1;
    TH1D *h_ratio_2;
    TLine *line;
    TH1D *h_systematics_low[3];
    TH1D *h_systematics_high[3];
    TH1D *hv2[3];
    TBox *box[4];

    legend0 = new TLegend(0.45,0.7,0.88,0.88);
    legend0->SetBorderSize(0);

    const int nbins_and = 10;
    Double_t edges_and[nbins_and + 1] = {6.1,7.1,7.6,8.1,8.6,9.1,9.6,10.1,10.6,11.1,11.6};

    const int nbins_xor = 9;
    Double_t edges_xor[nbins_xor + 1] = {7.1,8.6,9.1,9.6,10.1,10.6,11.1,11.6,12.1,12.6};

    const int nbins_minbias = 2;
    Double_t edges_minbias[nbins_minbias + 1] = {11.6,13.1,13.61};

    //bins for pt histogram
    const int nbins_pt = 5;
    Double_t edges_pt[nbins_pt + 1] = {0.5, 1.0, 2.0, 3.0, 4.0,5.0};

    std::vector<string> triggers = {"/","/minbias/","/xor/"};
    std::vector<string> trigger_name = {"AND","minbias","XOR"};

    TH1D *h_limits = new TH1D("h_limits","",1,3.6,13.61);
    h_limits->SetBinContent(1,-10);
    hs->Add(h_limits);
    hs_ratio->Add(h_limits);

    //plot vn for different effective energy peripheral bin
    for(int itrig =0; itrig<triggers.size(); itrig ++){
        legend0->Clear();
        input =  new TFile(Form("%s%sTemplateFits_vnn.root",base.c_str(),triggers.at(itrig).c_str()));
        for(int i=0; i<peri_cent.size(); i++){
            if(itrig ==0) h_zdc[i] = new TH1D(Form("hv2_and%i",i),"",nbins_and,edges_and);
            if(itrig ==1) h_zdc[i] = new TH1D(Form("hv2_minbias%i",i),"",nbins_minbias,edges_minbias);
            if(itrig ==2) h_zdc[i] = new TH1D(Form("hv2_xor%i",i),"",nbins_xor,edges_xor);
            
            h_zdc[i]->SetMarkerStyle(style.at(i));
            h_zdc[i]->SetMarkerColor(i +1);
            h_zdc[i]->SetLineColor(i +1);
    
            if(itrig ==0){
                for(int icent =0; icent< nbins_and; icent++){
                    int bin_idx = Bins::GetCentIndex(edges_and[icent],edges_and[icent +1]);
                    std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,ipt1,ipt2,2,1,2,peri_cent.at(i),peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                    h_zdc[i]->SetBinContent(icent + 1, vnn_zdc.first);
                    h_zdc[i]->SetBinError(icent + 1, vnn_zdc.second);
                }
            }
            if(itrig==1){
                for(int icent =0; icent< nbins_minbias; icent++){
                    int bin_idx = Bins::GetCentIndex(edges_minbias[icent],edges_minbias[icent +1]);
                    std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,ipt1,ipt2,2,1,2,peri_cent.at(i),peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                    h_zdc[i]->SetBinContent(icent + 1, vnn_zdc.first);
                    h_zdc[i]->SetBinError(icent + 1, vnn_zdc.second);
                }
            }
            if(itrig==2){
                for(int icent =0; icent< nbins_xor; icent++){
                    int bin_idx = Bins::GetCentIndex(edges_xor[icent],edges_xor[icent +1]);
                    std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,ipt1,ipt2,2,1,2,peri_cent.at(i),peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                    h_zdc[i]->SetBinContent(icent + 1, vnn_zdc.first);
                    h_zdc[i]->SetBinError(icent + 1, vnn_zdc.second);
                }
            }
    
            hs->Add(h_zdc[i]);
            legend0->AddEntry(h_zdc[i],Form("%s , %s ",Bins::label_cent_peri(peri_cent.at(i)).c_str(),Bins::label_trk_peri(peritrk).c_str()),"lep");
        }

        c0->cd();
        hs->Draw("nostack");
        if(itrig ==0) hs->GetXaxis()->SetLimits(edges_and[0],edges_and[nbins_and]);
        if(itrig ==1) hs->GetXaxis()->SetLimits(edges_minbias[0],edges_minbias[nbins_minbias]);
        if(itrig ==2) hs->GetXaxis()->SetLimits(edges_xor[0],edges_xor[nbins_xor]);
        
        hs->SetMaximum(0.09);
        hs->SetMinimum(0.045);
        legend0->Draw();
        X=0.20,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_ptab(ipt1,ipt2).c_str()), size, 43); Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta(Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s trigger",trigger_name.at(itrig).c_str()), size, 43);
    
        c0->SaveAs(Form("%s/v2_peri_%s.pdf",figures.c_str(),trigger_name.at(itrig).c_str()));
        c0->Clear();
        h_denominator = (TH1D*)h_zdc[0];// For ratios
        h_ratio_1 = (TH1D*)h_zdc[1]->Clone("h_ratio_1");
        h_ratio_2 = (TH1D*)h_zdc[2]->Clone("h_ratio_2");
        h_ratio_1->Divide(h_denominator);
        h_ratio_2->Divide(h_denominator);
    
        for(int ibin =1; ibin<h_ratio_1->GetNbinsX()+1; ibin++){
            h_ratio_1->SetBinError(ibin,0);
            h_ratio_2->SetBinError(ibin,0);
        }
        //for systematics
        h_systematics_low[itrig] = (TH1D*)h_ratio_1->Clone(Form("clone_low_trig%i",itrig));
        h_systematics_high[itrig] = (TH1D*)h_ratio_2->Clone(Form("clone_high_trig%i",itrig));
        hv2[itrig] = (TH1D*)h_zdc[0]->Clone(Form("clone_v2_trig%i",itrig));
        hs_ratio->Add(h_ratio_1);
        hs_ratio->Add(h_ratio_2);
        
        hs_ratio->SetMaximum(1.03);
        hs_ratio->SetMinimum(0.94);
        hs_ratio->Draw("nostack;L");

        if(itrig ==0) hs_ratio->GetXaxis()->SetLimits(edges_and[0],edges_and[nbins_and]);
        if(itrig ==1) hs_ratio->GetXaxis()->SetLimits(edges_minbias[0],edges_minbias[nbins_minbias]);
        if(itrig ==2) hs_ratio->GetXaxis()->SetLimits(edges_xor[0],edges_xor[nbins_xor]);
        gPad->SetGrid();
        c0->SaveAs(Form("%s/v2_peri_ratio_%s.pdf",figures.c_str(),trigger_name.at(itrig).c_str()));
        gPad->SetGrid(0,0);
        hs->GetHists()->Clear();
        hs_ratio->GetHists()->Clear();
        c0->Clear();

    }
    //-----------------------------------------------------------------------------------
    //now draw the original v2 but with systematics uncertainties
        // Add systematic uncertainties as shaded regions
        hv2[0]->SetMarkerColor(kBlue);
        hv2[1]->SetMarkerColor(kRed);
        hv2[2]->SetMarkerColor(kOrange+3);
        hs->Add(hv2[0]);
        hs->Add(hv2[1]);
        hs->Add(hv2[2]);
        hs->Draw("nostack");
        hs->GetXaxis()->SetLimits(edges_and[0],13.6);
        for(int itrig =0; itrig<3; itrig ++){
            if(itrig == 0 ){
                for (int j = 0; j < nbins_and; j++) {
                    double loww = hv2[0]->GetBinContent(j+1) * h_systematics_low[itrig]->GetBinContent(j+1);
                    double highh = hv2[0]->GetBinContent(j+1) * h_systematics_high[itrig]->GetBinContent(j+1);
                    box[itrig] = new TBox(edges_and[j], loww, edges_and[j+1], highh);
                    box[itrig]->SetFillColorAlpha(kBlue, 0.2);  // Grayish for syst error
                    box[itrig]->SetLineColor(0);
                    box[itrig]->Draw("same");
                }
            }
            if(itrig == 1){
                for (int j = 0; j < nbins_minbias; j++) {
                    double loww = hv2[itrig]->GetBinContent(j+1) * h_systematics_low[itrig]->GetBinContent(j+1);
                    double highh = hv2[itrig]->GetBinContent(j+1) * h_systematics_high[itrig]->GetBinContent(j+1);
                    box[itrig] = new TBox(edges_minbias[j], loww, edges_minbias[j+1], highh);
                    box[itrig]->SetFillColorAlpha(kRed, 0.2);  // Grayish for syst error
                    box[itrig]->SetLineColor(0);
                    box[itrig]->Draw("same");
                }
            }
            if(itrig == 2){
                for (int j = 0; j < nbins_xor; j++) {
                    double loww = hv2[itrig]->GetBinContent(j+1) * h_systematics_low[itrig]->GetBinContent(j+1);
                    double highh = hv2[itrig]->GetBinContent(j+1) * h_systematics_high[itrig]->GetBinContent(j+1);
                    box[itrig] = new TBox(edges_xor[j], loww, edges_xor[j+1], highh);
                    box[itrig]->SetFillColorAlpha(kOrange +3, 0.2);  // Grayish for syst error
                    box[itrig]->SetLineColor(0);
                    box[itrig]->Draw("same");
                }
            }
        }
        X=0.2,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp 22} #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_ptab(ipt1,ipt2).c_str()), size, 43); Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta (Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43);
        legend0 = new TLegend(0.45,0.7,0.88,0.88);
        legend0->SetBorderSize(0);
        legend0->AddEntry(hv2[0],"ZDC AND ","lep");
        legend0->AddEntry(hv2[1],"Minimum bias","lep");
        legend0->AddEntry(hv2[2],"ZDC XOR","lep");
        legend0->Draw();
        c0->SaveAs(Form("%s/v2_systematics.pdf",figures.c_str()));
        hs->GetHists()->Clear();
        legend0->Clear();
    //------------------------------------------------------------------------------------------------

    //plot chi2 test for AND & XOR
    c0->Clear();
    for(int index =0; index<3; index++){
        int nbinsx = hv2[index]->GetNbinsX();
        hv2[index]->Draw();
        TF1 *f_const = new TF1("f_const", "[0]", hv2[index]->GetXaxis()->GetXmin(), hv2[index]->GetXaxis()->GetXmax());
        //TF1 *f_const = new TF1("f_const", "[0]", hv2[index]->GetXaxis()->GetBinCenter(1), hv2[index]->GetXaxis()->GetBinCenter(nbinsx));
        f_const->SetParameter(0, 0.06); // Initial guess for the constant value
        hv2[index]->Fit(f_const, "QR"); // "Q" for quiet mode (avoid printouts)
        hv2[index]->GetYaxis()->SetTitle("v_{2}");
        hv2[index]->GetXaxis()->SetTitle("E_{Eff} [TeV]");
        hv2[index]->GetYaxis()->SetRangeUser(0.045,0.09);
        double chi2 = f_const->GetChisquare();
        int ndof = f_const->GetNDF();
        std::cout << "Chi2: " << chi2 << ", NDF: " << ndof << ", Chi2/NDF: " << chi2 / ndof << std::endl;
        X=0.2,Y=0.88;
        Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
        Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, "#it{pp 22} #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_ptab(ipt1,ipt2).c_str()), size, 43); Y=Y-0.05;
        Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta (Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43);
        X=0.6,Y=0.88;
        if (index == 0){ Common::myText2(X     ,Y,1,Form("Fit results: %.3f #pm %.3f",f_const->GetParameter(0),f_const->GetParError(0)),size,43);Y=Y-0.05;}
        if (index == 2){ Common::myText2(X     ,Y,1,Form("Fit results: %.4f #pm %.4f",f_const->GetParameter(0),f_const->GetParError(0)),size,43);Y=Y-0.05;}
        Common::myText2(X     ,Y,1,Form("#frac{#chi^{2}}{NDOF} = %.2f ",chi2 / ndof),size,43);Y=Y-0.05;
        c0->SaveAs(Form("%s/v2_and_chi2_%i.pdf",figures.c_str(),index));
    }
    //-----------------------------------------------------------------------------------
    //print out chi2 test for trigger consistent to each other.
    TH1D *hv2_and_rebin = (TH1D*)hv2[2]->Clone("clone_hv2_xor_trigger");
    hv2_and_rebin->SetBinContent(1,(hv2[0]->GetBinContent(2) + hv2[0]->GetBinContent(3) + hv2[0]->GetBinContent(4))/3);
    hv2_and_rebin->SetBinError(1,(hv2[0]->GetBinError(3)));
    for(int i=5; i<11; i++){
        hv2_and_rebin->SetBinContent(i-3, hv2[0]->GetBinContent(i));
        hv2_and_rebin->SetBinError(i-3,(hv2[0]->GetBinError(i))); 
    } 
    //set the last bin to 0 since AND trigger have no data there
    int last_bin = hv2_and_rebin->GetNbinsX();
    hv2_and_rebin->SetBinContent(last_bin,0);
    hv2_and_rebin->SetBinError(last_bin,0);
    hv2[2]->SetBinContent(last_bin,0);
    hv2[2]->SetBinError(last_bin,0);
    double chi2_trig = hv2_and_rebin->Chi2Test(hv2[2], "CHI2/NDF");  // Returns Chi2/NDF
    cout << "chi2 test for consistent between AND & XOR trigger: " << chi2_trig << endl;
    //plot pt variation
    c0->Clear();
    legend0->Clear();
    // legend0 = new TLegend(0.6,0.8,0.88,0.9);
    legend0->SetBorderSize(0);
    std::vector<int> pt1_val = {5,1,0};
    std::vector<int> pt2_val = {5,0};
    // style = {20,24,21};
    for(int itrig =0; itrig<triggers.size(); itrig++){
        input =  new TFile(Form("%s%sTemplateFits_vnn.root",base.c_str(),triggers.at(itrig).c_str()));
        for(int idxpt2 =0; idxpt2 < pt2_val.size(); idxpt2++){
            THStack *hs_pta = new THStack("hs_pta",";E_{Eff} [TeV];v_{2}");
            for(int i=0; i<pt1_val.size(); i++){
                if(itrig ==0) h_zdc[i] = new TH1D(Form("hv2_and%i_%i_%i",i,itrig,idxpt2),"",nbins_and,edges_and);
            else{
                h_zdc[i] = new TH1D(Form("hv2_xor%i_%i_%i",i,itrig,idxpt2),"",nbins_xor,edges_xor);
            }
                h_zdc[i]->SetMarkerStyle(style.at(i));
                h_zdc[i]->SetMarkerColor(kBlue + 2*i);
                h_zdc[i]->SetLineColor(kBlue + 2*i);
                if(itrig ==0){
                    for(int icent =0; icent< nbins_and ; icent++){
                        int bin_idx = Bins::GetCentIndex(edges_and[icent],edges_and[icent +1]);
                            std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,pt1_val.at(i),pt2_val.at(idxpt2),2,1,2,pericent,peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                            h_zdc[i]->SetBinContent(icent + 1, vnn_zdc.first);
                            h_zdc[i]->SetBinError(icent + 1, vnn_zdc.second);
                    }
                }
                else{
                    for(int icent =0; icent< nbins_xor ; icent++){
                        int bin_idx = Bins::GetCentIndex(edges_xor[icent],edges_xor[icent +1]);
                            std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(bin_idx,itrk,pt1_val.at(i),pt2_val.at(idxpt2),2,1,2,pericent,peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
                            h_zdc[i]->SetBinContent(icent + 1, vnn_zdc.first);
                            h_zdc[i]->SetBinError(icent + 1, vnn_zdc.second);
                    }
                }
                hs->Add(h_zdc[i]);
                legend0->AddEntry(h_zdc[i],Form("%s",Bins::label_pta(pt1_val.at(i)).c_str()),"lep");
            }
            hs->Draw("nostack");
            hs->SetMaximum(0.1);
            hs->SetMinimum(0.03);
            if(itrig ==0) hs->GetXaxis()->SetLimits(edges_and[0],edges_and[nbins_and]);
            else{
                hs->GetXaxis()->SetLimits(edges_xor[0],edges_xor[nbins_xor]);
            }
            legend0->Draw();
            X=0.20,Y=0.88;
            Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
            Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
            Common::myText2(X       , Y, 1, "#it{pp} #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
            Common::myText2(X       , Y, 1, Form("%s",Bins::label_ptb(pt2_val.at(idxpt2)).c_str()), size, 43); Y=Y-0.05;
            Common::myText2(X       , Y, 1, Form("%s",Bins::label_eta (Bins::GetDetaIndex(2.0,5.0)).c_str()), size, 43); Y=Y-0.05;
            Common::myText2(X       , Y, 1, Form("%s trigger",trigger_name.at(itrig).c_str()), size, 43);
            c0->SaveAs(Form("%s/v2_pta_ipt2_%i_%s.pdf",figures.c_str(),idxpt2,trigger_name.at(itrig).c_str()));
            c0->Clear();
            legend0->Clear();
            hs->GetHists()->Clear();
    }
    }
    //-----------------------------------------------------------------------------------

    //plot v2 vs PT for several effective energy bins
    hs->GetHists()->Clear();
    c0->Clear();
    hv2[0]= new TH1D(Form("hv2_pt"),"",nbins_pt,edges_pt);
    std::string trigger_s = "/";
    input =  new TFile(Form("%s%sTemplateFits_vnn.root",base.c_str(),trigger_s.c_str()));
    int icent= 30;
    for(int iptt1 =0; iptt1< nbins_pt ; iptt1++){
        //int bin_idx = Bins::GetCentIndex(icent,icent +1);
        std::pair<float, float> vnn_zdc=      Bins::GetVnPtb(icent,itrk,5,iptt1,2,1,2,pericent,peritrk,input,input); //icent_bin,itrk,ipt1,ipt2,ich,ideta,m_har,pericent,peritrk
        hv2[0]->SetBinContent(iptt1 + 1, vnn_zdc.first);
        hv2[0]->SetBinError(iptt1 + 1, vnn_zdc.second);
    }
    hv2[0]->Draw();
    // hs->Add(hv2[0]);
    // hs->Draw("nostack");
    // hs->GetXaxis()->SetTitle("p_T^{a} [GeV]");
    c0->SaveAs(Form("%s/v2_pt_and.pdf",figures.c_str()));
    //------------------------------------------------------------------------------------------------------------------------------------
    //pt distribution
    c0->Clear();
    trigger_s = "/";
    TH1D *h_pt_integraged[3];
    std::vector<int> colors = {kBlue,kRed,kOrange+3};
    THStack *hs_pt = new THStack("hs",";p_{T} [GeV];#frac{1}{N} #frac{dN}{dp_{T}N}");
    for(int itrig=0; itrig<3; itrig++){
        if(itrig == 1) continue;
        input =  new TFile(Form("%s%shistograms.root",base.c_str(),triggers.at(itrig).c_str()));
        h_pt_integraged[itrig] = (TH1D*)input->Get("h_pt_icent11")->Clone("clone_h_pt_icent00");
        h_pt_integraged[itrig]->SetMarkerColor(colors.at(itrig));
        h_pt_integraged[itrig]->SetLineColor(colors.at(itrig));
        //h_pt_integraged->GetYaxis()->SetTitle("#frac{#frac{dN}{dp_{T}N}}{N_{events}}");
        for(int icent =12; icent<15; icent++){
            TH1D *h_pt_temp = (TH1D*)input->Get(Form("h_pt_icent%02i",icent))->Clone(Form("clone_h_pt_icent%01i",icent));
            h_pt_integraged[itrig]->Add(h_pt_temp);
        }
        int bin_lo  =h_pt_integraged[itrig]->FindBin(Bins::PT1_LO[0]+.0001);
        int bin_high=h_pt_integraged[itrig]->FindBin(Bins::PT1_HI[4]-.0001);
        h_pt_integraged[itrig]->Scale(1.0/(h_pt_integraged[itrig]->Integral(bin_lo,bin_high)));
        hs_pt->Add(h_pt_integraged[itrig]);
    }
    hs_pt->Draw("nostack");
    hs_pt->GetXaxis()->SetLimits(0.5,5);
    hs_pt->SetMaximum(1.1);
    hs_pt->SetMinimum(1e-3);
    gPad->SetLogy();
    X=0.7,Y=0.88;
    Common::myText2(X     ,Y,1,"ATLAS "         ,size,73);
    Common::myText2(X+0.1,Y,1,Common::Internal,size,43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "#it{pp 22} #sqrt{#it{s}} = 13.6 TeV", size, 43);Y=Y-0.05;
    Common::myText2(X       , Y, 1, "E_{Eff} #in [9.1,11.6) TeV", size, 43);Y=Y-0.05;
    legend0 = new TLegend(0.2,0.3,0.5,0.5);
    legend0->SetBorderSize(0);
    legend0->AddEntry(h_pt_integraged[0],"ZDC AND ","lep");
    //legend0->AddEntry(h_pt_integraged[1],"Minimum bias","lep");
    legend0->AddEntry(h_pt_integraged[2],"ZDC XOR","lep");
    legend0->Draw();
    c0->SaveAs(Form("%s/pt_distribution.pdf",figures.c_str()));
}
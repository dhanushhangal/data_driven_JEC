#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include "TStopwatch.h"
#include "TEnv.h"

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>

#include <algorithm>
#include <vector>

const int nCbins = 1;
const int nbins_alpha = 6;

void draw_response(bool is_pp=0){

    double xbins_pt[] = {40.0,50.0,70.0,95.0,120.0,150.0,200.0,300.0,500.0}; 
    const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

    double xbins_alpha[] = {0.,0.05,0.1,0.2,0.3,0.4,0.5};
    const int nbins_alpha = sizeof(xbins_alpha)/sizeof(double)-1;

    double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 4.013, 5.191};
    const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;

    TString pt_tag[] = {"40<p_{T}<50", "50<p_{T}<70", "70<p_{T}<95", "95<p_{T}<120", "120<p_{T}<150", "150<p_{T}<200", "200<p_{T}<300", "300<p_{T}<500"};
    TString alpha_tag[] = {"#alpha<0.05", "#alpha<0.1", "#alpha<0.2", "#alpha<0.3", "#alpha<0.4", "#alpha<0.5",};
    TString eta_tag[] = {"|#eta| < 0.261", "0.261 < |#eta| < 0.522", "0.522 < |#eta| < 0.783", "0.783 < |#eta| < 1.044", "1.044 < |#eta| < 1.305"};

    TH1D *h_jetpt[nCbins];
    TH1D *h_jeteta[nCbins];
    TH1D *h_jetphi[nCbins];
  
    TH1D *h_phopt[nCbins];
    TH1D *h_phoeta[nCbins];
    TH1D *h_phophi[nCbins];

    TH1D *h_metpt[nCbins];
    TH1D *h_meteta[nCbins];
    TH1D *h_metphi[nCbins];

    TH2D *h_R_jtpt[nCbins][nbins_alpha];
    TH2D *h_R_mpf[nCbins][nbins_alpha];

    TH1D *h_R_jtpt_proj[nCbins][nbins_alpha][nbins_pt];
    TH1D *h_R_jtpt_data_proj[nCbins][nbins_alpha][nbins_pt];

    TH1D *h_R_mpf_proj[nCbins][nbins_alpha][nbins_pt];
    TH1D *h_R_mpf_data_proj[nCbins][nbins_alpha][nbins_pt];

    TH2D *h_R_jtpt_eta[nCbins][nbins_alpha];
    TH2D *h_R_mpf_eta[nCbins][nbins_alpha];

    TH2D *h_R_jtpt_data_eta[nCbins][nbins_alpha];
    TH2D *h_R_mpf_data_eta[nCbins][nbins_alpha];

    TH1D *h_R_jtpt_eta_proj[nCbins][nbins_alpha][nbins_eta];
    TH1D *h_R_jtpt_data_eta_proj[nCbins][nbins_alpha][nbins_eta];

    TH1D *h_R_mpf_eta_proj[nCbins][nbins_alpha][nbins_eta];
    TH1D *h_R_mpf_data_eta_proj[nCbins][nbins_alpha][nbins_eta];

    TH1D *h_R_jtpt_px[nCbins][nbins_alpha];
    TH1D *h_R_mpf_px[nCbins][nbins_alpha];

    TH1D *h_R_jtpt_alpha[nCbins][nbins_pt];
    TH1D *h_R_mpf_alpha[nCbins][nbins_pt];

    TH2D *h_R_jtpt_data[nCbins][nbins_alpha];
    TH2D *h_R_mpf_data[nCbins][nbins_alpha];

    TH1D *h_R_jtpt_data_px[nCbins][nbins_alpha];
    TH1D *h_R_mpf_data_px[nCbins][nbins_alpha];

    TH1D *h_R_jtpt_data_alpha[nCbins][nbins_pt];
    TH1D *h_R_mpf_data_alpha[nCbins][nbins_pt];

    TH1D *h_R_mpf_eta_alpha[nCbins][nbins_pt];
    TH1D *h_R_mpf_eta_data_alpha[nCbins][nbins_pt];

    TH1D *h_R_jtpt_eta_alpha[nCbins][nbins_pt];
    TH1D *h_R_jtpt_eta_data_alpha[nCbins][nbins_pt];

    TF1 *f_R_jtpt_proj[nCbins][nbins_alpha][nbins_pt];
    TF1 *f_R_jtpt_data_proj[nCbins][nbins_alpha][nbins_pt];

    TF1 *f_R_mpf_proj[nCbins][nbins_alpha][nbins_pt];
    TF1 *f_R_mpf_data_proj[nCbins][nbins_alpha][nbins_pt];

    TF1 *f_R_jtpt_eta_proj[nCbins][nbins_alpha][nbins_eta];
    TF1 *f_R_jtpt_data_eta_proj[nCbins][nbins_alpha][nbins_eta];

    TF1 *f_R_mpf_eta_proj[nCbins][nbins_alpha][nbins_eta];
    TF1 *f_R_mpf_data_eta_proj[nCbins][nbins_alpha][nbins_eta];

    TF1 *f_R_mpf_alpha[nCbins][nbins_pt];
    TF1 *f_R_jtpt_alpha[nCbins][nbins_pt];

    TF1 *f_R_jtpt_eta_alpha[nCbins][nbins_eta];
    TF1 *f_R_mpf_eta_alpha[nCbins][nbins_eta];

    TF1 *f_response_MPF = new TF1("f_response_MPF","[0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395)))",50.,500.);
    f_response_MPF->SetLineColor(kAzure+1);
    TF1 *f_response_MPF_two = new TF1("f_response_MPF_two","[0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))",50.,500.);
    f_response_MPF_two->SetLineColor(kViolet+1);
    TF1 *f_response_MPF_lin = new TF1("f_response_MPF_lin","[0]+[1]*log(x)",50.,500.);
    f_response_MPF_lin->SetLineColor(kGreen-2);
    TF1 *f_response_MPF_quad = new TF1("f_response_MPF_quad","[0]+[1]*log(x)+[2]*log(x)*log(x)",50.,500.);
    f_response_MPF_quad->SetLineColor(kOrange+1);
    TF1 *f_response_jtpt = new TF1("f_response_MPF","1./([0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))",50.,500.);

    TFile *ppMC_file;
    TFile *ppdata_file;

    if(is_pp){
        ppMC_file = TFile::Open("ppMC2017_L3Res_Sep2_histos.root");  
        ppdata_file = TFile::Open("ppdata2017_L3Res_Sep2_histos.root");  
    }
    else{
        ppMC_file = TFile::Open("PbPbMC2018_L3Res_Aug28_histos.root");  
        ppdata_file = TFile::Open("PbPbdata2018_L3Res_Aug28_histos.root");  
    }

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=0; ibin2<nbins_alpha;ibin2++){

            h_R_mpf_eta[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_mpf_eta_cent%d_alpha%d",ibin,ibin2));
            h_R_jtpt_eta[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_jtpt_eta_cent%d_alpha%d",ibin,ibin2));

            if(is_pp){
                h_R_mpf[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_mpf_cent%d_alpha%d",ibin,ibin2));
                h_R_jtpt[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_jtpt_cent%d_alpha%d",ibin,ibin2));

                h_R_mpf_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_mpf_cent%d_alpha%d",ibin,ibin2))->Clone(Form("h_R_mpf_data_cent%d_alpha%d",ibin,ibin2));
                h_R_jtpt_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_jtpt_cent%d_alpha%d",ibin,ibin2))->Clone(Form("h_R_jtpt_data_cent%d_alpha%d",ibin,ibin2));

                h_R_mpf_data_eta[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_mpf_eta_cent%d_alpha%d",ibin,ibin2))->Clone(Form("h_R_mpf_eta_data_cent%d_alpha%d",ibin,ibin2));
                h_R_jtpt_data_eta[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_jtpt_eta_cent%d_alpha%d",ibin,ibin2))->Clone(Form("h_R_jtpt_eta_data_cent%d_alpha%d",ibin,ibin2));

                if(ibin2>0) h_R_mpf_eta[ibin][ibin2]->Add(h_R_mpf_eta[ibin][ibin2-1]);
                if(ibin2>0) h_R_jtpt_eta[ibin][ibin2]->Add(h_R_jtpt_eta[ibin][ibin2-1]);

                if(ibin2>0) h_R_mpf_data_eta[ibin][ibin2]->Add(h_R_mpf_data_eta[ibin][ibin2-1]);
                if(ibin2>0) h_R_jtpt_data_eta[ibin][ibin2]->Add(h_R_jtpt_data_eta[ibin][ibin2-1]);
            }
            else{
                h_R_mpf[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_mpf_cent%d_alpha%d",3,ibin2));
                h_R_jtpt[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_jtpt_cent%d_alpha%d",3,ibin2));

                h_R_mpf_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_mpf_cent%d_alpha%d",3,ibin2))->Clone(Form("h_R_mpf_data_cent%d_alpha%d",ibin,ibin2));
                h_R_jtpt_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_jtpt_cent%d_alpha%d",3,ibin2))->Clone(Form("h_R_jtpt_data_cent%d_alpha%d",ibin,ibin2));                
            }

 //h_R_jtpt[ibin][ibin2]->RebinX(8);
 //h_R_jtpt_data[ibin][ibin2]->RebinX(8);

            if(ibin2>0) h_R_mpf[ibin][ibin2]->Add(h_R_mpf[ibin][ibin2-1]);
            if(ibin2>0) h_R_jtpt[ibin][ibin2]->Add(h_R_jtpt[ibin][ibin2-1]);

            if(ibin2>0) h_R_mpf_data[ibin][ibin2]->Add(h_R_mpf_data[ibin][ibin2-1]);
            if(ibin2>0) h_R_jtpt_data[ibin][ibin2]->Add(h_R_jtpt_data[ibin][ibin2-1]);
/*

            h_R_mpf_px[ibin][ibin2] = h_R_mpf[ibin][ibin2]->ProjectionX(Form("h_R_mpf_px_%d_%d",ibin,ibin2),1,-1,"");
            h_R_jtpt_px[ibin][ibin2] = h_R_jtpt[ibin][ibin2]->ProjectionX(Form("h_R_jtpt_px_%d_%d",ibin,ibin2),1,-1,"");

            h_R_mpf_data_px[ibin][ibin2] = h_R_mpf_data[ibin][ibin2]->ProjectionX(Form("h_R_mpf_data_px_%d_%d",ibin,ibin2),1,-1,"");
            h_R_jtpt_data_px[ibin][ibin2] = h_R_jtpt_data[ibin][ibin2]->ProjectionX(Form("h_R_jtpt_data_px_%d_%d",ibin,ibin2),1,-1,"");
*/

            h_R_mpf_px[ibin][ibin2] = new TH1D(Form("h_R_mpf_px_cent%d_alpha%d",ibin,ibin2),"",nbins_pt, xbins_pt);
            h_R_mpf_px[ibin][ibin2]->Sumw2();
            h_R_mpf_data_px[ibin][ibin2] = new TH1D(Form("h_R_mpf_data_px_cent%d_alpha%d",ibin,ibin2),"",nbins_pt, xbins_pt);
            h_R_mpf_data_px[ibin][ibin2]->Sumw2();

            h_R_jtpt_px[ibin][ibin2] = new TH1D(Form("h_R_jtpt_px_cent%d_alpha%d",ibin,ibin2),"",nbins_pt, xbins_pt);
            h_R_jtpt_px[ibin][ibin2]->Sumw2();
            h_R_jtpt_data_px[ibin][ibin2] = new TH1D(Form("h_R_jtpt_data_px_cent%d_alpha%d",ibin,ibin2),"",nbins_pt, xbins_pt);
            h_R_jtpt_data_px[ibin][ibin2]->Sumw2();

            for(int ibin3=0; ibin3<h_R_mpf[ibin][ibin2]->GetNbinsX();ibin3++){
                h_R_mpf_proj[ibin][ibin2][ibin3] = h_R_mpf[ibin][ibin2]->ProjectionY(Form("h_R_mpf_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               
                h_R_mpf_data_proj[ibin][ibin2][ibin3] = h_R_mpf_data[ibin][ibin2]->ProjectionY(Form("h_R_mpf_data_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               

                h_R_jtpt_proj[ibin][ibin2][ibin3] = h_R_jtpt[ibin][ibin2]->ProjectionY(Form("h_R_jtpt_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               
                h_R_jtpt_data_proj[ibin][ibin2][ibin3] = h_R_jtpt_data[ibin][ibin2]->ProjectionY(Form("h_R_jtpt_data_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               

                f_R_mpf_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_mpf_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                f_R_mpf_proj[ibin][ibin2][ibin3]->SetLineColor(kRed);
                f_R_mpf_data_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_mpf_data_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                f_R_mpf_data_proj[ibin][ibin2][ibin3]->SetLineColor(kBlue);

                f_R_jtpt_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_jtpt_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                f_R_jtpt_proj[ibin][ibin2][ibin3]->SetLineColor(kRed);
                f_R_jtpt_data_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_jtpt_data_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                f_R_jtpt_data_proj[ibin][ibin2][ibin3]->SetLineColor(kBlue);
            }
            
            if(is_pp){
                for(int ibin3=0; ibin3<h_R_mpf_eta[ibin][ibin2]->GetNbinsX();ibin3++){
                    h_R_mpf_eta_proj[ibin][ibin2][ibin3] = h_R_mpf_eta[ibin][ibin2]->ProjectionY(Form("h_R_mpf_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               
                    h_R_mpf_data_eta_proj[ibin][ibin2][ibin3] = h_R_mpf_data_eta[ibin][ibin2]->ProjectionY(Form("h_R_mpf_data_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               

                    h_R_jtpt_eta_proj[ibin][ibin2][ibin3] = h_R_jtpt_eta[ibin][ibin2]->ProjectionY(Form("h_R_jtpt_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");               
                    h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3] = h_R_jtpt_data_eta[ibin][ibin2]->ProjectionY(Form("h_R_jtpt_eta_data_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");                               

                    f_R_mpf_eta_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_mpf_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                    f_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kRed);
                    f_R_mpf_data_eta_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_mpf_data_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                    f_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kBlue);

                    f_R_jtpt_eta_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_jtpt_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                    f_R_jtpt_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kRed);
                    f_R_jtpt_data_eta_proj[ibin][ibin2][ibin3] = new TF1(Form("f_R_jtpt_data_eta_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                    f_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kBlue);
                }
            }
	    }

        //vs. alpha
        for(int ibin2=0; ibin2<h_R_mpf[0][0]->GetNbinsX();ibin2++){
            h_R_mpf_alpha[ibin][ibin2] = new TH1D(Form("h_R_mpf_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_mpf_alpha[ibin][ibin2]->Sumw2();
            h_R_mpf_data_alpha[ibin][ibin2] = new TH1D(Form("h_R_mpf_data_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_mpf_data_alpha[ibin][ibin2]->Sumw2();

            h_R_jtpt_alpha[ibin][ibin2] = new TH1D(Form("h_R_jtpt_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_jtpt_alpha[ibin][ibin2]->Sumw2();
            h_R_jtpt_data_alpha[ibin][ibin2] = new TH1D(Form("h_R_jtpt_data_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_jtpt_data_alpha[ibin][ibin2]->Sumw2();
        }
        //vs. alpha
        for(int ibin2=0; ibin2<5;ibin2++){
            h_R_jtpt_eta_alpha[ibin][ibin2] = new TH1D(Form("h_R_jtpt_eta_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_jtpt_eta_alpha[ibin][ibin2]->Sumw2();
            h_R_jtpt_eta_data_alpha[ibin][ibin2] = new TH1D(Form("h_R_jtpt_eta_data_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_jtpt_eta_data_alpha[ibin][ibin2]->Sumw2();

            h_R_mpf_eta_alpha[ibin][ibin2] = new TH1D(Form("h_R_mpf_eta_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_mpf_eta_alpha[ibin][ibin2]->Sumw2();
            h_R_mpf_eta_data_alpha[ibin][ibin2] = new TH1D(Form("h_R_mpf_eta_data_alpha_cent%d_alpha%d",ibin,ibin2),"",nbins_alpha, xbins_alpha);
            h_R_mpf_eta_data_alpha[ibin][ibin2]->Sumw2();
        }
    }
/*
    TCanvas *c_sample = new TCanvas("c_sample","c_sample",500,500);
    c_sample->cd();
    h_R_mpf_eta[0][0]->Draw("COLZ");
*/
    TLine *l1 = new TLine(50.,1.,500.,1.);
    TLine *l2 = new TLine(0.,1.,0.5,1.);
    TLine *l3 = new TLine(40.,0.9,200.,0.9);
    TLine *l4 = new TLine(0.,1.,1.5,1.);
    l1->SetLineStyle(2); l2->SetLineStyle(2); l3->SetLineStyle(2); l4->SetLineStyle(2);

    TLatex *tx1;
    tx1 = new TLatex(); tx1->SetTextSize(.09);
    TLatex *tx2;
    tx2 = new TLatex(); tx2->SetTextSize(.08);

    TLegend *legendc = new TLegend(0.25,0.6,0.65,0.85);
    legendc ->SetLineColor(kWhite);
    legendc ->AddEntry(h_R_mpf_proj[0][3][1], "Pythia8", "lepf");
    if(is_pp) legendc ->AddEntry(h_R_mpf_data_proj[0][3][1], "pp data", "lepf");
    else legendc ->AddEntry(h_R_mpf_data_proj[0][3][1], "PbPb data (70-90%)", "lepf");

    TCanvas *c_mpf_fits[nbins_alpha];
    
    for(int ibin2=0; ibin2<nbins_alpha;ibin2++){    
        c_mpf_fits[ibin2] = new TCanvas(Form("c_mpf_fits_%d",ibin2),Form("c_mpf_fits_%d",ibin2),1000,500);
        c_mpf_fits[ibin2]->Divide(4,2,0);
        for(int ibin=0; ibin<nCbins;ibin++){
            for(int ibin3=0; ibin3<h_R_mpf[ibin][ibin2]->GetNbinsX();ibin3++){
                c_mpf_fits[ibin2]->cd(ibin3+1);
                h_R_mpf_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_mpf_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_mpf_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_mpf_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_mpf_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.55);
                h_R_mpf_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);
                h_R_mpf_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1));

                h_R_mpf_data_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_mpf_data_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_data_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_proj[ibin][ibin2][ibin3],"Q M R","sames",0.7*f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1));

                if(ibin3==0){
                    legendc->Draw("same");
                }
                tx2->DrawLatexNDC(0.37,0.92, pt_tag[ibin3]);                
            }
        }
    }    
/*
    TCanvas *c_mpf_fits_eta[nbins_alpha];

    for(int ibin2=0; ibin2<nbins_alpha;ibin2++){    
        c_mpf_fits_eta[ibin2] = new TCanvas(Form("c_mpf_fits_eta_%d",ibin2),Form("c_mpf_fits_eta_%d",ibin2),1000,500);
        c_mpf_fits_eta[ibin2]->Divide(3,2,0);

        for(int ibin=0; ibin<nCbins;ibin++){
            for(int ibin3=0; ibin3<5;ibin3++){

                c_mpf_fits_eta[ibin2]->cd(ibin3+1);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.55);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                f_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetParameter(0,0.15);
                f_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetParameter(1,1.);
                f_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetParameter(2,0.25);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,2.);
                //if(ibin3!=3 || ibin2>1) h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetParameter(1));

                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,2.);
                //h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->GetParameter(1));

                if(ibin3==0){
                    legendc->Draw("same");
                }
                tx2->DrawLatexNDC(0.37,0.92, eta_tag[ibin3]);                

            }
        }
    }    
*/
    TCanvas *c_jtpt_fits[nbins_alpha];

    for(int ibin2=0; ibin2<nbins_alpha;ibin2++){    
        c_jtpt_fits[ibin2] = new TCanvas(Form("c_jtpt_fits_%d",ibin2),Form("c_jtpt_fits_%d",ibin2),1000,500);
        c_jtpt_fits[ibin2]->Divide(4,2,0);
        for(int ibin=0; ibin<nCbins;ibin++){
            for(int ibin3=0; ibin3<h_R_jtpt[ibin][ibin2]->GetNbinsX();ibin3++){
                c_jtpt_fits[ibin2]->cd(ibin3+1);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_jtpt_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_jtpt_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Scale(1./h_R_jtpt_proj[ibin][ibin2][ibin3]->Integral());
                h_R_jtpt_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.55);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_proj[ibin][ibin2][ibin3],"Q N M R","sames",0.,3.);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1));

                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_jtpt_data_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Scale(1./h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Integral());
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_data_proj[ibin][ibin2][ibin3],"Q N M R","sames",0.,3.);
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_data_proj[ibin][ibin2][ibin3],"Q M R","sames",0.7*f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1));

                if(ibin3==0){
                    legendc->Draw("same");
                }
                tx2->DrawLatexNDC(0.37,0.92, pt_tag[ibin3]);
            }
        }
    }    

    TCanvas *c_jtpt_fits_eta[nbins_alpha];

    for(int ibin2=0; ibin2<nbins_alpha;ibin2++){    
        c_jtpt_fits_eta[ibin2] = new TCanvas(Form("c_jtpt_fits_eta_%d",ibin2),Form("c_jtpt_fits_eta_%d",ibin2),1000,500);
        c_jtpt_fits_eta[ibin2]->Divide(3,2,0);

        for(int ibin=0; ibin<nCbins;ibin++){
            for(int ibin3=0; ibin3<5;ibin3++){

                c_jtpt_fits_eta[ibin2]->cd(ibin3+1);
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->Scale(1./h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->Integral());
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.55);
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);
                h_R_jtpt_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1));

                h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->Scale(1./h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->Integral());
                h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_data_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);
                h_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_data_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.7*f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1));

                if(ibin3==0){
                    legendc->Draw("same");
                }
                tx2->DrawLatexNDC(0.37,0.92, eta_tag[ibin3]);                

            }
        }
    }

    TCanvas *c_mpf_fits_eta[nbins_alpha];

    for(int ibin2=0; ibin2<nbins_alpha;ibin2++){    
        c_mpf_fits_eta[ibin2] = new TCanvas(Form("c_mpf_fits_eta_%d",ibin2),Form("c_mpf_fits_eta_%d",ibin2),1000,500);
        c_mpf_fits_eta[ibin2]->Divide(3,2,0);

        for(int ibin=0; ibin<nCbins;ibin++){
            for(int ibin3=0; ibin3<5;ibin3++){

                c_mpf_fits_eta[ibin2]->cd(ibin3+1);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_mpf_eta_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.55);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);
                h_R_mpf_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1));

                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);
                h_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_eta_proj[ibin][ibin2][ibin3],"Q M R","sames",0.7*f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1));

                if(ibin3==0){
                    legendc->Draw("same");
                }
                tx2->DrawLatexNDC(0.37,0.92, eta_tag[ibin3]);                

            }
        }
    }

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=0; ibin2<nbins_alpha;ibin2++){
            for(int ibin3=0; ibin3<h_R_jtpt[ibin][ibin2]->GetNbinsX();ibin3++){
                h_R_mpf_px[ibin][ibin2]->SetBinContent(ibin3+1,f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_px[ibin][ibin2]->SetBinError(ibin3+1,f_R_mpf_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_px[ibin][ibin2]->SetBinContent(ibin3+1,f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_px[ibin][ibin2]->SetBinError(ibin3+1,f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_mpf_data_px[ibin][ibin2]->SetBinContent(ibin3+1,f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_data_px[ibin][ibin2]->SetBinError(ibin3+1,f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_data_px[ibin][ibin2]->SetBinContent(ibin3+1,f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_data_px[ibin][ibin2]->SetBinError(ibin3+1,f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParError(1));

                //vs. alpha
                h_R_mpf_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_mpf_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_mpf_data_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_data_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_data_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_data_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParError(1));
            }

            for(int ibin3=0; ibin3<5;ibin3++){
                //vs. alpha
                h_R_jtpt_eta_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_jtpt_eta_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_eta_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_jtpt_eta_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_eta_data_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_eta_data_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_jtpt_data_eta_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_mpf_eta_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_eta_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_mpf_eta_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_mpf_eta_data_alpha[ibin][ibin3]->SetBinContent(ibin2+1,f_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_eta_data_alpha[ibin][ibin3]->SetBinError(ibin2+1,f_R_mpf_data_eta_proj[ibin][ibin2][ibin3]->GetParError(1));
            }

            h_R_mpf_px[ibin][ibin2]->SetMarkerColor(kRed+2); h_R_mpf_px[0][ibin2]->SetLineColor(kRed+2);
            h_R_mpf_data_px[ibin][ibin2]->SetMarkerColor(kRed+2); h_R_mpf_data_px[0][ibin2]->SetLineColor(kRed+2);

            h_R_jtpt_px[ibin][ibin2]->SetMarkerColor(kGreen-2); h_R_jtpt_px[0][ibin2]->SetLineColor(kGreen-2);
            h_R_jtpt_data_px[ibin][ibin2]->SetMarkerColor(kGreen-2); h_R_jtpt_data_px[0][ibin2]->SetLineColor(kGreen-2);
        }
    }    
cout<<"here"<<endl;
    TCanvas *c_R_mpf = new TCanvas("c_R_mpf","c_R_mpf",900,300);
    c_R_mpf->Divide(3,1);
    gStyle->SetOptStat(0);

    TLegend *legenda = new TLegend(0.25,0.6,0.65,0.85);
    legenda ->SetLineColor(kWhite);
    legenda ->AddEntry(h_R_mpf_px[0][3], "MPF ; Pythia8", "lepf");
    if(is_pp) legenda ->AddEntry(h_R_mpf_data_px[0][3], "MPF ; pp data", "lepf");
    else legenda ->AddEntry(h_R_mpf_data_px[0][3], "MPF ; PbPb data (70-90%)", "lepf");

    TLegend *legendb = new TLegend(0.25,0.6,0.65,0.85);
    legendb ->SetLineColor(kWhite);
    legendb ->AddEntry(h_R_jtpt_px[0][3], "p_{T bal} ; Pythia8", "lepf");
    if(is_pp) legendb ->AddEntry(h_R_jtpt_data_px[0][3], "p_{T bal} ; pp data", "lepf");
    else  legendb ->AddEntry(h_R_jtpt_data_px[0][3], "p_{T bal} ; PbPb data (70-90%)", "lepf");

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=2; ibin2<nbins_alpha-1;ibin2++){
            c_R_mpf->cd(ibin2+1-2);
            gPad->SetLogx();
    	    h_R_mpf_px[ibin][ibin2]->GetYaxis()->SetTitle("Absolute Response");
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetTitleSize(0.06);
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.6);
            h_R_mpf_px[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_R_mpf_px[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
    	    h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetRangeUser(50.,500.);
    	    h_R_mpf_px[ibin][ibin2]->SetMarkerStyle(20);        
            h_R_mpf_data_px[ibin][ibin2]->SetMarkerStyle(4);        
            h_R_jtpt_px[ibin][ibin2]->SetMarkerStyle(20);        
            h_R_jtpt_data_px[ibin][ibin2]->SetMarkerStyle(4);        
            //h_R_mpf_px[ibin][ibin2]->SetMarkerStyle(20);        
    	    h_R_mpf_px[ibin][ibin2]->Draw("e0 same");
            h_R_mpf_data_px[ibin][ibin2]->Draw("e0 same");
            //h_R_jtpt_px[ibin][ibin2]->Draw("e0 same");
            //h_R_jtpt_data_px[ibin][ibin2]->Draw("e0 same");
            l1->Draw("same"); //l2->Draw("same"); l3->Draw("same");

            if(ibin==0 and ibin2==3) legenda->Draw("same");
            tx1->DrawLatexNDC(0.37,0.92, alpha_tag[ibin2]);
        }
	}

    TCanvas *c_R_jtpt = new TCanvas("c_R_jtpt","c_R_jtpt",900,300);
    c_R_jtpt->Divide(3,1);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=2; ibin2<nbins_alpha-1;ibin2++){
            c_R_jtpt->cd(ibin2+1-2);
            gPad->SetLogx();
            h_R_jtpt_px[ibin][ibin2]->GetYaxis()->SetTitle("Absolute Response");
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetTitleSize(0.06);
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.6);
            h_R_jtpt_px[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_R_jtpt_px[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetRangeUser(50.,500.);
            //h_R_jtpt_px[ibin][ibin2]->SetMarkerStyle(4);        
            h_R_jtpt_px[ibin][ibin2]->SetMarkerStyle(20);        
            h_R_jtpt_data_px[ibin][ibin2]->SetMarkerStyle(4);        
            h_R_jtpt_px[ibin][ibin2]->Draw("e0 same");
            h_R_jtpt_data_px[ibin][ibin2]->Draw("e0 same");
            l1->Draw("same"); //l2->Draw("same"); l3->Draw("same");

            if(ibin==0 and ibin2==3) legendb->Draw("same");
            tx1->DrawLatexNDC(0.37,0.92, alpha_tag[ibin2]);
        }
    }

    TH1D *h_R_jtpt_px_data_MC[nCbins][nbins_alpha];
    TH1D *h_R_mpf_px_data_MC[nCbins][nbins_alpha];

    TLegend *legende = new TLegend(0.25,0.6,0.65,0.85);
    legende ->SetLineColor(kWhite);
    legende ->AddEntry(h_R_jtpt_data_px[0][3], "p_{T bal}", "lepf");
    legende ->AddEntry(h_R_mpf_data_px[0][3], "MPF", "lepf");

    TCanvas *c_R_data_mc = new TCanvas("c_R_data_mc","c_R_data_mc",900,300);
    c_R_data_mc->Divide(3,1);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=0; ibin2<nbins_alpha;ibin2++){

            h_R_jtpt_px_data_MC[ibin][ibin2] = (TH1D*)h_R_jtpt_data_px[ibin][ibin2]->Clone(Form("h_R_jtpt_px_data_MC_%d_%d",ibin,ibin2));
            h_R_jtpt_px_data_MC[ibin][ibin2]->Divide(h_R_jtpt_px[ibin][ibin2]);
            h_R_mpf_px_data_MC[ibin][ibin2] = (TH1D*)h_R_mpf_data_px[ibin][ibin2]->Clone(Form("h_R_mpf_px_data_MC_%d_%d",ibin,ibin2));
            h_R_mpf_px_data_MC[ibin][ibin2]->Divide(h_R_mpf_px[ibin][ibin2]);
            
            if(ibin2 < 2 || ibin2 ==5) continue;
            c_R_data_mc->cd(ibin2+1-2);
            gPad->SetLogx();

            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->SetTitle("Relative Response (Data/MC)");
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetXaxis()->SetTitleSize(0.06);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.6);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.95);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_R_jtpt_px_data_MC[ibin][ibin2]->GetXaxis()->SetRangeUser(50.,500.);
            h_R_jtpt_px_data_MC[ibin][ibin2]->Draw("same");

            h_R_mpf_px_data_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.9,1.1);
            h_R_mpf_px_data_MC[ibin][ibin2]->Draw("same");
            l1->Draw("same");

            if(ibin==0 and ibin2==3) legende->Draw("same");
            tx1->DrawLatexNDC(0.37,0.92, alpha_tag[ibin2]);
        }
    }

    //vs. alpha
    TH1D *h_R_jtpt_alpha_data_MC[nCbins][nbins_alpha];
    TH1D *h_R_mpf_alpha_data_MC[nCbins][nbins_alpha];

    TLegend *legendf = new TLegend(0.25,0.6,0.65,0.85);
    legendf ->SetLineColor(kWhite);
    legendf ->AddEntry(h_R_jtpt_data_alpha[0][3], "p_{T bal}", "lepf");
    legendf ->AddEntry(h_R_mpf_data_alpha[0][3], "MPF", "lepf");

    TCanvas *c_R_data_mc_alpha = new TCanvas("c_R_data_mc_alpha","c_R_data_mc_alpha",1200,600);
    c_R_data_mc_alpha->Divide(4,2);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=0; ibin2<nbins_pt;ibin2++){

            f_R_mpf_alpha[ibin][ibin2] = new TF1(Form("f_R_mpf_alpha_%d_%d",ibin,ibin2),"pol1",0.,0.5);
            f_R_jtpt_alpha[ibin][ibin2] = new TF1(Form("f_R_jtpt_alpha_%d_%d",ibin,ibin2),"pol1",0.,0.5);
            f_R_mpf_alpha[ibin][ibin2]->SetLineColor(kRed+2); 
            f_R_jtpt_alpha[ibin][ibin2]->SetLineColor(kGreen-2); 

            c_R_data_mc_alpha->cd(ibin2+1);

            h_R_mpf_data_alpha[ibin][ibin2]->SetMarkerStyle(4); h_R_mpf_data_alpha[ibin][ibin2]->SetMarkerColor(kRed+2); h_R_mpf_data_alpha[ibin][ibin2]->SetLineColor(kRed+2);
            h_R_jtpt_data_alpha[ibin][ibin2]->SetMarkerStyle(4); h_R_jtpt_data_alpha[ibin][ibin2]->SetLineColor(kGreen-2);

            h_R_jtpt_alpha_data_MC[ibin][ibin2] = (TH1D*)h_R_jtpt_data_alpha[ibin][ibin2]->Clone(Form("h_R_jtpt_alpha_data_MC_%d_%d",ibin,ibin2));
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->Divide(h_R_jtpt_alpha[ibin][ibin2]);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetTitle("Relative Response (Data/MC)");
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->SetTitleSize(0.06);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.6);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.95);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->SetNdivisions(505);            
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->SetTitle("#alpha");
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->CenterTitle();
            //h_R_jtpt_alpha_data_MC[ibin][ibin2]->GetXaxis()->SetRangeUser(50.,500.);
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->Draw("same");
            h_R_jtpt_alpha_data_MC[ibin][ibin2]->Fit(f_R_jtpt_alpha[ibin][ibin2],"Q M R","sames",0.,0.5);

            h_R_mpf_alpha_data_MC[ibin][ibin2] = (TH1D*)h_R_mpf_data_alpha[ibin][ibin2]->Clone(Form("h_R_mpf_alpha_data_MC_%d_%d",ibin,ibin2));
            h_R_mpf_alpha_data_MC[ibin][ibin2]->Divide(h_R_mpf_alpha[ibin][ibin2]);
            h_R_mpf_alpha_data_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.9,1.1);
            h_R_mpf_alpha_data_MC[ibin][ibin2]->Draw("same");
            h_R_mpf_alpha_data_MC[ibin][ibin2]->Fit(f_R_mpf_alpha[ibin][ibin2],"Q M R","sames",0.,0.5);
            l2->Draw("same");

            if(ibin==0 and ibin2==3) legendf->Draw("same");
            tx2->DrawLatexNDC(0.37,0.92, pt_tag[ibin2]);
        }
    }

    //vs. alpha
    TH1D *h_R_jtpt_alpha_eta_data_MC[nCbins][nbins_alpha];
    TH1D *h_R_mpf_alpha_eta_data_MC[nCbins][nbins_alpha];

    TCanvas *c_R_eta_data_mc_alpha = new TCanvas("c_R_eta_data_mc_alpha","c_R_eta_data_mc_alpha",1200,600);
    c_R_eta_data_mc_alpha->Divide(3,2);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=0; ibin2<5;ibin2++){

            f_R_jtpt_eta_alpha[ibin][ibin2] = new TF1(Form("f_R_jtpt_eta_alpha_%d_%d",ibin,ibin2),"pol1",0.,0.5);
            f_R_jtpt_eta_alpha[ibin][ibin2]->SetLineColor(kGreen-2); 

            f_R_mpf_eta_alpha[ibin][ibin2] = new TF1(Form("f_R_mpf_eta_alpha_%d_%d",ibin,ibin2),"pol1",0.,0.5);
            f_R_mpf_eta_alpha[ibin][ibin2]->SetLineColor(kGreen-2); 

            c_R_eta_data_mc_alpha->cd(ibin2+1);

            h_R_mpf_eta_data_alpha[ibin][ibin2]->SetMarkerStyle(4); h_R_mpf_eta_data_alpha[ibin][ibin2]->SetLineColor(kRed-2);
            h_R_jtpt_eta_data_alpha[ibin][ibin2]->SetMarkerStyle(4); h_R_jtpt_eta_data_alpha[ibin][ibin2]->SetLineColor(kGreen-2);

            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2] = (TH1D*)h_R_jtpt_eta_data_alpha[ibin][ibin2]->Clone(Form("h_R_jtpt_alpha_eta_data_MC_%d_%d",ibin,ibin2));
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->Divide(h_R_jtpt_eta_alpha[ibin][ibin2]);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetTitle("Relative Response (Data/MC)");
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->SetTitleSize(0.06);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.6);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.95);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->SetNdivisions(505);            
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->SetTitle("#alpha");
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->CenterTitle();
            //h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->GetXaxis()->SetRangeUser(50.,500.);
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->Draw("same");
            h_R_jtpt_alpha_eta_data_MC[ibin][ibin2]->Fit(f_R_jtpt_eta_alpha[ibin][ibin2],"Q M R","sames",0.,0.5);

            h_R_mpf_alpha_eta_data_MC[ibin][ibin2] = (TH1D*)h_R_mpf_eta_data_alpha[ibin][ibin2]->Clone(Form("h_R_mpf_alpha_eta_data_MC_%d_%d",ibin,ibin2));
            h_R_mpf_alpha_eta_data_MC[ibin][ibin2]->Divide(h_R_mpf_eta_alpha[ibin][ibin2]);
            h_R_mpf_alpha_eta_data_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_R_mpf_alpha_eta_data_MC[ibin][ibin2]->Draw("same");
            h_R_mpf_alpha_eta_data_MC[ibin][ibin2]->Fit(f_R_mpf_eta_alpha[ibin][ibin2],"Q M R","sames",0.,0.5);

            l2->Draw("same");

            if(ibin==0 and ibin2==3) legendf->Draw("same");
            tx2->DrawLatexNDC(0.2,0.92, eta_tag[ibin2]);
        }
    }

    TH1D *h_R_jtpt_px_data_MC_alpha[nCbins];
    TH1D *h_R_mpf_px_data_MC_alpha[nCbins];

    TH1D *h_R_jtpt_px_eta_data_MC_alpha[nCbins];
    TH1D *h_R_mpf_px_eta_data_MC_alpha[nCbins];

    for(int ibin=0; ibin<nCbins;ibin++){
        h_R_mpf_px_data_MC_alpha[ibin] = new TH1D(Form("h_R_mpf_px_data_MC_alpha_cent%d",ibin),"",nbins_pt, xbins_pt);
        h_R_mpf_px_data_MC_alpha[ibin]->Sumw2();
        h_R_jtpt_px_data_MC_alpha[ibin] = new TH1D(Form("h_R_jtpt_px_data_MC_alpha_cent%d",ibin),"",nbins_pt, xbins_pt);
        h_R_jtpt_px_data_MC_alpha[ibin]->Sumw2();        
        h_R_jtpt_px_eta_data_MC_alpha[ibin] = new TH1D(Form("h_R_jtpt_px_eta_data_MC_alpha_cent%d",ibin),"",nbins_eta, xbins_eta);
        h_R_jtpt_px_eta_data_MC_alpha[ibin]->Sumw2();        
        h_R_mpf_px_eta_data_MC_alpha[ibin] = new TH1D(Form("h_R_mpf_px_eta_data_MC_alpha_cent%d",ibin),"",nbins_eta, xbins_eta);
        h_R_mpf_px_eta_data_MC_alpha[ibin]->Sumw2();        

        for(int ibin2=0; ibin2<nbins_pt;ibin2++){
            h_R_mpf_px_data_MC_alpha[ibin]->SetBinContent(ibin2+1,f_R_mpf_alpha[ibin][ibin2]->Eval(0.));
            h_R_mpf_px_data_MC_alpha[ibin]->SetBinError(ibin2+1,f_R_mpf_alpha[ibin][ibin2]->GetParError(0));
            h_R_jtpt_px_data_MC_alpha[ibin]->SetBinContent(ibin2+1,f_R_jtpt_alpha[ibin][ibin2]->Eval(0.));
            h_R_jtpt_px_data_MC_alpha[ibin]->SetBinError(ibin2+1,f_R_jtpt_alpha[ibin][ibin2]->GetParError(0));
        }
        for(int ibin2=0; ibin2<5;ibin2++){
            h_R_jtpt_px_eta_data_MC_alpha[ibin]->SetBinContent(ibin2+1,f_R_jtpt_eta_alpha[ibin][ibin2]->Eval(0.));
            h_R_jtpt_px_eta_data_MC_alpha[ibin]->SetBinError(ibin2+1,f_R_jtpt_eta_alpha[ibin][ibin2]->GetParError(0));            

            h_R_mpf_px_eta_data_MC_alpha[ibin]->SetBinContent(ibin2+1,f_R_mpf_eta_alpha[ibin][ibin2]->Eval(0.));
            h_R_mpf_px_eta_data_MC_alpha[ibin]->SetBinError(ibin2+1,f_R_mpf_eta_alpha[ibin][ibin2]->GetParError(0));            
        }
    }

    TCanvas *c_R_data_mc_alpha0 = new TCanvas("c_R_data_mc_alpha0","c_R_data_mc_alpha0",500,500);
    gStyle->SetOptStat(0);
    gPad->SetLogx();

    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->SetTitle("Relative Response (Data/MC)");
    h_R_jtpt_px_data_MC_alpha[0]->GetXaxis()->SetTitleSize(0.06);
    h_R_jtpt_px_data_MC_alpha[0]->GetXaxis()->SetTitleOffset(0.6);
    h_R_jtpt_px_data_MC_alpha[0]->GetXaxis()->SetLabelSize(0.05);
    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->SetTitleSize(0.05);
    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->SetTitleOffset(0.95);
    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->SetLabelSize(0.05);
    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->CenterTitle();
    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->SetNdivisions(505);
    h_R_jtpt_px_data_MC_alpha[0]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
    h_R_jtpt_px_data_MC_alpha[0]->GetXaxis()->CenterTitle();
    h_R_jtpt_px_data_MC_alpha[0]->GetXaxis()->SetRangeUser(50.,500.);
    h_R_jtpt_px_data_MC_alpha[0]->GetYaxis()->SetRangeUser(0.8,1.2);
    h_R_mpf_px_data_MC_alpha[0]->SetMarkerStyle(4);h_R_mpf_px_data_MC_alpha[0]->SetMarkerColor(kRed+2);h_R_mpf_px_data_MC_alpha[0]->SetLineColor(kRed+2); 
    h_R_jtpt_px_data_MC_alpha[0]->SetMarkerStyle(4);h_R_jtpt_px_data_MC_alpha[0]->SetMarkerColor(kGreen-2);h_R_jtpt_px_data_MC_alpha[0]->SetLineColor(kGreen-2); 
    h_R_jtpt_px_data_MC_alpha[0]->Draw("e0 same");
    h_R_mpf_px_data_MC_alpha[0]->Draw("e0 same");
/*    
    h_R_mpf_px_data_MC_alpha[0]->Fit(f_response_MPF,"M R","sames",50.,500.);
    cout<<f_response_MPF->GetChisquare()/f_response_MPF->GetNDF()<<endl;
    h_R_mpf_px_data_MC_alpha[0]->Fit(f_response_MPF_two,"M R","sames",50.,500.);
    cout<<f_response_MPF_two->GetChisquare()/f_response_MPF_two->GetNDF()<<endl;
    h_R_mpf_px_data_MC_alpha[0]->Fit(f_response_MPF_lin,"M R","sames",50.,500.);
    cout<<f_response_MPF_lin->GetChisquare()/f_response_MPF_lin->GetNDF()<<endl;
    h_R_mpf_px_data_MC_alpha[0]->Fit(f_response_MPF_quad,"M R","sames",50.,500.);
    cout<<f_response_MPF_quad->GetChisquare()/f_response_MPF_quad->GetNDF()<<endl;    
    f_response_MPF->Draw("same");
    f_response_MPF_two->Draw("same");
    f_response_MPF_lin->Draw("same");
    f_response_MPF_quad->Draw("same");
*/
    //h_R_jtpt_px_data_MC_alpha[0]->Fit(f_response_jtpt,"M R","sames",50.,500.);
    legende->Draw("same");
    l1->Draw("same");


    TCanvas *c_R_eta_data_mc_alpha0 = new TCanvas("c_R_eta_data_mc_alpha0","c_R_eta_data_mc_alpha0",500,500);
    gStyle->SetOptStat(0);

    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->SetTitle("Relative Response (Data/MC)");
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetXaxis()->SetTitleSize(0.06);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetXaxis()->SetTitleOffset(0.6);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetXaxis()->SetLabelSize(0.05);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->SetTitleSize(0.05);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->SetTitleOffset(0.95);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->SetLabelSize(0.05);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->CenterTitle();
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->SetNdivisions(505);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetXaxis()->SetTitle("#eta");
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetXaxis()->CenterTitle();
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetXaxis()->SetRangeUser(0.,1.3);
    h_R_jtpt_px_eta_data_MC_alpha[0]->GetYaxis()->SetRangeUser(0.8,1.2);
    h_R_mpf_px_eta_data_MC_alpha[0]->SetMarkerStyle(4);h_R_mpf_px_eta_data_MC_alpha[0]->SetMarkerColor(kRed+2);h_R_mpf_px_eta_data_MC_alpha[0]->SetLineColor(kRed+2); 
    h_R_jtpt_px_eta_data_MC_alpha[0]->SetMarkerStyle(4);h_R_jtpt_px_eta_data_MC_alpha[0]->SetMarkerColor(kGreen-2);h_R_jtpt_px_eta_data_MC_alpha[0]->SetLineColor(kGreen-2); 
    h_R_jtpt_px_eta_data_MC_alpha[0]->Draw("e0 same");
    h_R_mpf_px_eta_data_MC_alpha[0]->Draw("e0 same");
    legende->Draw("same");
    l4->Draw("same");

/*
    TFile *output_file = new TFile("outer.root", "RECREATE");
    output_file->cd();

    for(int ibin2=0; ibin2<nbins_alpha;ibin2++){
        h_R_mpf_data_px[0][ibin2]->Write();
        h_R_mpf_px[0][ibin2]->Write();
        h_R_jtpt_data_px[0][ibin2]->Write();
        h_R_jtpt_px[0][ibin2]->Write();
        h_R_mpf_px_data_MC[0][ibin2]->Write();
        h_R_jtpt_px_data_MC[0][ibin2]->Write();
    }

    h_R_mpf_px_data_MC_alpha[0]->Write();
    h_R_jtpt_px_data_MC_alpha[0]->Write();    

    output_file->Close();
*/








}
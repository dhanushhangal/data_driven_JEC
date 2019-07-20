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

void draw_response(bool is_pp=1){

    double xbins_pt[] = {40.0,50.0,70.0,95.0,120.0,150.0,200.0,300.0,500.0}; 
    const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

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

    TH1D *h_R_jtpt_px[nCbins][nbins_alpha];
    TH1D *h_R_mpf_px[nCbins][nbins_alpha];

    TH2D *h_R_jtpt_data[nCbins][nbins_alpha];
    TH2D *h_R_mpf_data[nCbins][nbins_alpha];

    TH1D *h_R_jtpt_data_px[nCbins][nbins_alpha];
    TH1D *h_R_mpf_data_px[nCbins][nbins_alpha];

    TF1 *f_R_jtpt_proj[nCbins][nbins_alpha][nbins_pt];
    TF1 *f_R_jtpt_data_proj[nCbins][nbins_alpha][nbins_pt];

    TF1 *f_R_mpf_proj[nCbins][nbins_alpha][nbins_pt];
    TF1 *f_R_mpf_data_proj[nCbins][nbins_alpha][nbins_pt];

    TFile *ppMC_file;
    TFile *ppdata_file;

    if(is_pp){
        ppMC_file = TFile::Open("ppMC2017_L3Res_Jul20_histos.root");  
        ppdata_file = TFile::Open("ppdata2017_L3Res_Jul20_histos.root");  
    }
    else{
        ppMC_file = TFile::Open("PbPbMC2018_L3Res_Jul17_histos_old.root");  
        ppdata_file = TFile::Open("PbPbdata2018_L3Res_Jul17_histos_old.root");  
    }

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=0; ibin2<nbins_alpha;ibin2++){

            h_R_mpf[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_mpf_cent%d_alpha%d",ibin,ibin2));
	    	h_R_jtpt[ibin][ibin2] = (TH2D*)ppMC_file->Get(Form("h_R_jtpt_cent%d_alpha%d",ibin,ibin2));

            if(is_pp){
                h_R_mpf_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_mpf_cent%d_alpha%d",ibin,ibin2))->Clone(Form("h_R_mpf_data_cent%d_alpha%d",ibin,ibin2));
                h_R_jtpt_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_jtpt_cent%d_alpha%d",ibin,ibin2))->Clone(Form("h_R_jtpt_data_cent%d_alpha%d",ibin,ibin2));
            }
            else{
                h_R_mpf_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_mpf_cent%d_alpha%d",3,ibin2))->Clone(Form("h_R_mpf_data_cent%d_alpha%d",ibin,ibin2));
                h_R_jtpt_data[ibin][ibin2] = (TH2D*)ppdata_file->Get(Form("h_R_jtpt_cent%d_alpha%d",3,ibin2))->Clone(Form("h_R_jtpt_data_cent%d_alpha%d",ibin,ibin2));                
            }

            if(ibin2>0 && ibin2!=5) h_R_mpf[ibin][ibin2]->Add(h_R_mpf[ibin][ibin2-1]);
            if(ibin2>0 && ibin2!=5) h_R_jtpt[ibin][ibin2]->Add(h_R_jtpt[ibin][ibin2-1]);

            if(ibin2>0 && ibin2!=5) h_R_mpf_data[ibin][ibin2]->Add(h_R_mpf_data[ibin][ibin2-1]);
            if(ibin2>0 && ibin2!=5) h_R_jtpt_data[ibin][ibin2]->Add(h_R_jtpt_data[ibin][ibin2-1]);

            h_R_mpf_px[ibin][ibin2] = h_R_mpf[ibin][ibin2]->ProjectionX(Form("h_R_mpf_px_%d_%d",ibin,ibin2),1,-1,"");
            h_R_jtpt_px[ibin][ibin2] = h_R_jtpt[ibin][ibin2]->ProjectionX(Form("h_R_jtpt_px_%d_%d",ibin,ibin2),1,-1,"");

            h_R_mpf_data_px[ibin][ibin2] = h_R_mpf_data[ibin][ibin2]->ProjectionX(Form("h_R_mpf_data_px_%d_%d",ibin,ibin2),1,-1,"");
            h_R_jtpt_data_px[ibin][ibin2] = h_R_jtpt_data[ibin][ibin2]->ProjectionX(Form("h_R_jtpt_data_px_%d_%d",ibin,ibin2),1,-1,"");

            for(int ibin3=1; ibin3<h_R_mpf[ibin][ibin2]->GetNbinsX()+1;ibin3++){
                h_R_mpf_proj[ibin][ibin2][ibin3] = h_R_mpf[ibin][ibin2]->ProjectionY(Form("h_R_mpf_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3,ibin3,"");               
                h_R_mpf_data_proj[ibin][ibin2][ibin3] = h_R_mpf_data[ibin][ibin2]->ProjectionY(Form("h_R_mpf_data_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3,ibin3,"");               

                h_R_jtpt_proj[ibin][ibin2][ibin3] = h_R_jtpt[ibin][ibin2]->ProjectionY(Form("h_R_jtpt_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3,ibin3,"");               
                h_R_jtpt_data_proj[ibin][ibin2][ibin3] = h_R_jtpt_data[ibin][ibin2]->ProjectionY(Form("h_R_jtpt_data_proj_%d_%d_%d",ibin,ibin2,ibin3),ibin3,ibin3,"");               

                f_R_mpf_proj[ibin][ibin2][ibin3] = new TF1(Form("h_R_mpf_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                f_R_mpf_data_proj[ibin][ibin2][ibin3] = new TF1(Form("h_R_mpf_data_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 

                f_R_jtpt_proj[ibin][ibin2][ibin3] = new TF1(Form("h_R_jtpt_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
                f_R_jtpt_data_proj[ibin][ibin2][ibin3] = new TF1(Form("h_R_jtpt_data_proj_%d_%d_%d",ibin,ibin2,ibin3),"gaus",0.,3.); 
            }
	    }
    }

    TLine *l1 = new TLine(40.,1.,200.,1.);
    TLine *l2 = new TLine(40.,1.1,200.,1.1);
    TLine *l3 = new TLine(40.,0.9,200.,0.9);
    l1->SetLineStyle(2); l2->SetLineStyle(2); l3->SetLineStyle(2);

    TLegend *legendc = new TLegend(0.25,0.6,0.65,0.98);
    legendc ->SetLineColor(kWhite);
    legendc ->AddEntry(h_R_mpf_proj[0][3][1], "Pythia8", "lepf");
    if(is_pp) legendc ->AddEntry(h_R_mpf_data_proj[0][3][1], "pp data", "lepf");
    else legendc ->AddEntry(h_R_mpf_data_proj[0][3][1], "PbPb data (50-100%)", "lepf");

    TCanvas *c_mpf_fits = new TCanvas("c_mpf_fits","c_mpf_fits",1000,500);
    c_mpf_fits->Divide(4,2,0);
    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=3; ibin2<4;ibin2++){    
            for(int ibin3=1; ibin3<h_R_mpf[ibin][ibin2]->GetNbinsX()+1;ibin3++){
                c_mpf_fits->cd(ibin3);
                h_R_mpf_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_mpf_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_mpf_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_mpf_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_mpf_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.5);
                h_R_mpf_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);

                h_R_mpf_data_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_mpf_data_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Scale(1./h_R_mpf_data_proj[ibin][ibin2][ibin3]->Integral());
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_mpf_data_proj[ibin][ibin2][ibin3]->Fit(f_R_mpf_data_proj[ibin][ibin2][ibin3],"Q M R","sames",0.,3.);

                if(ibin3==1){
                    legendc->Draw("same");
                }
            }
        }
    }    

    TCanvas *c_jtpt_fits = new TCanvas("c_jtpt_fits","c_jtpt_fits",1000,500);
    c_jtpt_fits->Divide(4,2,0);
    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=3; ibin2<4;ibin2++){    
            for(int ibin3=1; ibin3<h_R_jtpt[ibin][ibin2]->GetNbinsX()+1;ibin3++){
                c_jtpt_fits->cd(ibin3);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->SetLineColor(kRed); h_R_jtpt_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed); h_R_jtpt_proj[ibin][ibin2][ibin3]->SetMarkerStyle(4);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Scale(1./h_R_jtpt_proj[ibin][ibin2][ibin3]->Integral());
                h_R_jtpt_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.,0.5);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_proj[ibin][ibin2][ibin3],"Q N M R","sames",0.,3.);
                h_R_jtpt_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_proj[ibin][ibin2][ibin3],"Q M R","sames",0.6*f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1));

                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlue); h_R_jtpt_data_proj[ibin][ibin2][ibin3]->SetMarkerStyle(5);
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Scale(1./h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Integral());
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Draw("e0 same");
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_data_proj[ibin][ibin2][ibin3],"Q N M R","sames",0.,3.);
                h_R_jtpt_data_proj[ibin][ibin2][ibin3]->Fit(f_R_jtpt_data_proj[ibin][ibin2][ibin3],"Q M R","sames",0.7*f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1),1.8*f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1));

                if(ibin3==1){
                    legendc->Draw("same");
                }
            }
        }
    }    

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=3; ibin2<4;ibin2++){
            for(int ibin3=1; ibin3<h_R_jtpt[ibin][ibin2]->GetNbinsX()+1;ibin3++){
                h_R_mpf_px[ibin][ibin2]->SetBinContent(ibin3,f_R_mpf_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_px[ibin][ibin2]->SetBinError(ibin3,f_R_mpf_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_px[ibin][ibin2]->SetBinContent(ibin3,f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_px[ibin][ibin2]->SetBinError(ibin3,f_R_jtpt_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_mpf_data_px[ibin][ibin2]->SetBinContent(ibin3,f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_mpf_data_px[ibin][ibin2]->SetBinError(ibin3,f_R_mpf_data_proj[ibin][ibin2][ibin3]->GetParError(1));

                h_R_jtpt_data_px[ibin][ibin2]->SetBinContent(ibin3,f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParameter(1));
                h_R_jtpt_data_px[ibin][ibin2]->SetBinError(ibin3,f_R_jtpt_data_proj[ibin][ibin2][ibin3]->GetParError(1));
            }

            h_R_mpf_px[ibin][ibin2]->SetMarkerColor(kRed+2); h_R_mpf_px[0][ibin2]->SetLineColor(kRed+2);
            h_R_mpf_data_px[ibin][ibin2]->SetMarkerColor(kRed+2); h_R_mpf_data_px[0][ibin2]->SetLineColor(kRed+2);

            h_R_jtpt_px[ibin][ibin2]->SetMarkerColor(kGreen-2); h_R_jtpt_px[0][ibin2]->SetLineColor(kGreen-2);
            h_R_jtpt_data_px[ibin][ibin2]->SetMarkerColor(kGreen-2); h_R_jtpt_data_px[0][ibin2]->SetLineColor(kGreen-2);
        }
    }    


    TCanvas *c_R_mpf = new TCanvas("c_R_mpf","c_R_mpf",1200,450);
    c_R_mpf->Divide(3,1);
    gStyle->SetOptStat(0);

    TLegend *legenda = new TLegend(0.25,0.6,0.65,0.98);
    legenda ->SetLineColor(kWhite);
    legenda ->AddEntry(h_R_mpf_px[0][3], "MPF ; Pythia8", "lepf");
    if(is_pp) legenda ->AddEntry(h_R_mpf_data_px[0][3], "MPF ; pp data", "lepf");
    else legenda ->AddEntry(h_R_mpf_data_px[0][3], "MPF ; PbPb data (50-100%)", "lepf");

    TLegend *legendb = new TLegend(0.25,0.6,0.65,0.98);
    legendb ->SetLineColor(kWhite);
    legendb ->AddEntry(h_R_jtpt_px[0][3], "p_{T bal} ; Pythia8", "lepf");
    if(is_pp) legendb ->AddEntry(h_R_jtpt_data_px[0][3], "p_{T bal} ; pp data", "lepf");
    else  legendb ->AddEntry(h_R_jtpt_data_px[0][3], "p_{T bal} ; PbPb data (50-100%)", "lepf");

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=3; ibin2<4;ibin2++){
            c_R_mpf->cd(ibin2+1-2);
    	    h_R_mpf_px[ibin][ibin2]->GetYaxis()->SetTitle("Absolute Response");
            h_R_mpf_px[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
            h_R_mpf_px[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_R_mpf_px[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
    	    h_R_mpf_px[ibin][ibin2]->GetXaxis()->SetRangeUser(40.,200.);
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
        }
	}

    TH1D *h_R_jtpt_px_data_MC;
    TH1D *h_R_mpf_px_data_MC;

    TCanvas *c_R_jtpt = new TCanvas("c_R_jtpt","c_R_jtpt",1200,450);
    c_R_jtpt->Divide(3,1);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCbins;ibin++){
        for(int ibin2=3; ibin2<4;ibin2++){
            c_R_jtpt->cd(ibin2+1-2);
            h_R_jtpt_px[ibin][ibin2]->GetYaxis()->SetTitle("Absolute Response");
            h_R_jtpt_px[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_R_jtpt_px[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_R_jtpt_px[ibin][ibin2]->GetXaxis()->SetRangeUser(40.,200.);
            //h_R_jtpt_px[ibin][ibin2]->SetMarkerStyle(4);        
            h_R_jtpt_px[ibin][ibin2]->SetMarkerStyle(20);        
            h_R_jtpt_data_px[ibin][ibin2]->SetMarkerStyle(4);        
            h_R_jtpt_px[ibin][ibin2]->Draw("e0 same");
            h_R_jtpt_data_px[ibin][ibin2]->Draw("e0 same");
            l1->Draw("same"); //l2->Draw("same"); l3->Draw("same");

            if(ibin==0 and ibin2==3) legendb->Draw("same");
        }
    }











}
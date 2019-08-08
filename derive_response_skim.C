#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "TLorentzVector.h"
#include <TF1.h>
#include "TNtuple.h"
#include "TProfile.h"
#include "assert.h"

#include "JetCorrector.h"

using namespace std;

const int nCbins = 4;
float CBins[nCbins+1] = {-1., 20., 60., 100., 200.};

int mycbin=-999, pbin=-999;

const int npthat_bins = 6;
double pthatbins[npthat_bins+1] = {15,30,50,80,120,170,9999};
double xsecs[npthat_bins+1] = {2.452E-04, 6.914E-05, 1.631E-05, 3.067E-06, 7.068E-07, 1.812E-07, 0};
int pthatEntries[npthat_bins] = {1000000, 1000000, 1000000, 1000000, 1000000, 1000000};

int findAlphaBin(float alpha, int nbins_alpha, double* xbins_alpha){

  for(int i=0; i<nbins_alpha; i++){
    if(alpha<xbins_alpha[i+1]) return i;
  }
  return nbins_alpha-1;
}

int findEtaBin(float eta, int nbins_eta, double* xbins_eta){

  for(int i=0; i<nbins_eta; i++){
    if(abs(eta)<xbins_eta[i+1]) return i;
  }
  return nbins_eta-1;
}

int findPtBin(float pt, int nbins_pt, double* xbins_pt){

  for(int i=0; i<nbins_pt; i++){
    if(pt<xbins_pt[i+1]) return i;
  }
  return nbins_pt-1;
}

TLorentzVector findMissEtJets(int nref, float* pt_F, float* rawpt_F, float* eta_F, float* phi_F, float* m_F, float low_corrpt){
  
  double sum_ex = 0.0;
  double sum_ey = 0.0;
  TLorentzVector missEtJets(0,0,0,0);

  for(int i=0; i<nref; i++){
	  if (pt_F[i]<12.) continue; 
	  if (rawpt_F[i]<low_corrpt && rawpt_F[i]>6){
		  sum_ex += rawpt_F[i]*cos(phi_F[i]);
		  sum_ey += rawpt_F[i]*sin(phi_F[i]);  
	  }
	  //else {
	  if (rawpt_F[i]>low_corrpt){
		  sum_ex += pt_F[i]*cos(phi_F[i]);
		  sum_ey += pt_F[i]*sin(phi_F[i]);  
	  }
  }
  double met = sqrt( sum_ex*sum_ex + sum_ey*sum_ey );
  missEtJets.SetPxPyPzE( -sum_ex, -sum_ey, 0.0, met);

  return missEtJets;
}

void derive_response_skim(bool is_pp=1, bool is_data=1){

    vector<string> Files;
    if(is_pp){
        Files.push_back("Spring18_ppRef5TeV_V1_MC_L2Relative_AK4PF.txt");
        Files.push_back("Spring18_ppRef5TeV_V1_DATA_L2Residual_AK4PF.txt");
    }
    else{
        Files.push_back("Autumn18_HI_V4_MC_L2Relative_AK4PF.txt");
        Files.push_back("Autumn18_HI_V2_DATA_L2Residual_AK4PF.txt");
    }

    string corr_file;
    if(is_pp) corr_file = "Spring18_ppRef5TeV_V1_MC_L2Relative_AK4PF.txt";
    else corr_file = "Autumn18_HI_V1_MC_L2Relative_AK4PF.txt";

    JetCorrector L1(corr_file);
    JetCorrector JEC(Files);

    float lowest_corrpt;
    if(is_pp) lowest_corrpt = 15.;
    else lowest_corrpt = 20.;

    double xbins_alpha[] = {0.,0.05,0.1,0.2,0.3,0.4,0.5};
    const int nbins_alpha = sizeof(xbins_alpha)/sizeof(double)-1;

	double xbins_pt[] = {40.0,50.0,70.0,95.0,120.0,150.0,200.0,300.0,500.0}; 
    //double xbins_pt[] = {40.0,45.0,50.0,55.,60.,65.,70.0,75.,80.,85.,90.,95.0,120.0,150.0,200.0,300.0,500.0}; 
	const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

    double xbins_eta[] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 4.013, 5.191};
    const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;

    //define histograms
    TH1D *h_vz = new TH1D("h_vz", "h_vz", 40.,-20.,20.);
    TH1D *h_pthat = new TH1D("h_pthat", "h_pthat", 40.,0.,400.);
    TH1D *h_hibin = new TH1D("h_hibin", "h_hibin", 200,0.,200.);
    h_vz->Sumw2();
    h_pthat->Sumw2();
    h_hibin->Sumw2();

    TH1D *h_jetpt[nCbins];
    TH1D *h_jeteta[nCbins];
    TH1D *h_jetphi[nCbins];

    TH1D *h_dphi[nCbins];

    TH1D *h_phopt[nCbins];
    TH1D *h_phoeta[nCbins];
    TH1D *h_phophi[nCbins];

    TH1D *h_metpt[nCbins];
    TH1D *h_meteta[nCbins];
    TH1D *h_metphi[nCbins];

    TH2D *h_R_jtpt[nCbins][nbins_alpha];
    TH2D *h_R_mpf[nCbins][nbins_alpha];

    for(int ibin=0; ibin<nCbins; ibin++){

        h_jetpt[ibin] = new TH1D(Form("h_jetpt_cent%d", ibin),"",nbins_pt, xbins_pt);
        h_jeteta[ibin] = new TH1D(Form("h_jeteta_cent%d",ibin),"",nbins_eta, xbins_eta);
        h_jetphi[ibin] = new TH1D(Form("h_jetphi_cent%d",ibin),"",30,-3.142,3.142);
        h_jetpt[ibin]->Sumw2();
        h_jeteta[ibin]->Sumw2();
        h_jetphi[ibin]->Sumw2();

        h_metpt[ibin] = new TH1D(Form("h_metpt_cent%d", ibin),"",nbins_pt, xbins_pt);
        h_meteta[ibin] = new TH1D(Form("h_meteta_cent%d",ibin),"",nbins_eta, xbins_eta);
        h_metphi[ibin] = new TH1D(Form("h_metphi_cent%d",ibin),"",30,-3.142,3.142);
        h_metpt[ibin]->Sumw2();
        h_meteta[ibin]->Sumw2();
        h_metphi[ibin]->Sumw2();

        h_phopt[ibin] = new TH1D(Form("h_phopt_cent%d", ibin),"",nbins_pt, xbins_pt);
        h_phoeta[ibin] = new TH1D(Form("h_phoeta_cent%d",ibin),"",nbins_eta, xbins_eta);
        h_phophi[ibin] = new TH1D(Form("h_phophi_cent%d",ibin),"",30,-3.142,3.142);
        h_phopt[ibin]->Sumw2();
        h_phoeta[ibin]->Sumw2();
        h_phophi[ibin]->Sumw2();

        h_dphi[ibin] = new TH1D(Form("h_dphi_cent%d", ibin),"",100,0.,3.14);
        h_dphi[ibin] -> Sumw2();

	    for(int ibin2=0; ibin2<nbins_alpha; ibin2++){
	      h_R_jtpt[ibin][ibin2] = new TH2D(Form("h_R_jtpt_cent%d_alpha%d",ibin,ibin2),"",nbins_pt, xbins_pt,30,0.,3.);
	      h_R_jtpt[ibin][ibin2]->Sumw2();

          h_R_mpf[ibin][ibin2] = new TH2D(Form("h_R_mpf_cent%d_alpha%d",ibin,ibin2),"",nbins_pt, xbins_pt,30,0.,3.);
          h_R_mpf[ibin][ibin2]->Sumw2();
	    }
	}

    TFile *my_file;

    if(is_data){
      if(is_pp) my_file = TFile::Open("HighEGJet_Run2017G-17Nov2017-v2_FOREST_9410_pho_v1_part0.root");
      else my_file = TFile::Open("HIHardProbes_HIRun2018A-04Apr2019-v1_FOREST_1033patch1_pho40.root");       
    }
    else{
      if(is_pp) my_file = TFile::Open("QCDPhoton_TuneCP5_5p02TeV_pythia8_FOREST_9410_v4.root");
      else my_file = TFile::Open("PbPbMC2018_L3Res_Jul14.root");
    }

    TTree *mixing_tree = (TTree*) my_file->Get("pj");

    int MAXJETS = 1000;
    UInt_t run=-999, lumi = -999;
    Float_t vz = -999, hihf = -999, pthat = -999; 
    Int_t hiBin = -999;
    //Int_t nref=0, nPho=0;

    float pthat_weight=0.;

    //CUTS
    const double eta_min = 0.7;
    const double eta_max = 1.3;
    const double sec_eta_max = 5.;
    const double phoEt_min = 40.;
    const double jtpt_min = 12.;
/*
    const float hovere_max = 0.006778;
    const float see_min = 0.;
    const float see_max = 0.009555;
    const float iso_max = -0.010780;
*/
    const float hovere_max_pp = 0.009732;
    const float see_min_pp = 0.;
    const float see_max_pp = 0.009905;
    const float iso_max_pp = -0.014755;

    const float hovere_max_PbPb = 0.119947;
    const float see_min_PbPb = 0.;
    const float see_max_PbPb = 0.010392;
    const float iso_max_PbPb = 2.099277;

    Float_t jtpt_arr[MAXJETS], jteta_arr[MAXJETS], jtphi_arr[MAXJETS], jtm_arr[MAXJETS], rawpt_arr[MAXJETS];

    //vector<float> *phoEt=0,*phoEta=0,*phoPhi=0;
    vector<float> *pfcIso3pTgt2p0subUE=0, *pfpIso3subUE=0, *pfnIso3subUE=0;
    vector<float> *phoEta=0, *phoPhi=0, *phoEt=0, *pho_swissCrx=0, *pho_seedTime=0, *phoSigmaIEtaIEta_2012=0, *phoHoverE=0, *pho_ecalClusterIsoR3=0, *pho_hcalRechitIsoR3=0, *pho_trackIsoR3PtCut20=0;
    vector<float> *jtpt=0, *rawpt=0, *jteta=0, *jtphi=0, *jtm=0, *refpt = 0, *trackmax = 0; 
    Int_t n_evt, nref, nPFpart, nPho;
    Int_t HLT_HIGEDPhoton20_v1, HLT_HIGEDPhoton20_v1_Prescl, HLT_HIGEDPhoton30_v1, HLT_HIGEDPhoton30_v1_Prescl, HLT_HIGEDPhoton40_v1, HLT_HIGEDPhoton40_v1_Prescl;

    Float_t corrpt[MAXJETS];

    double alpha;

	mixing_tree->SetBranchAddress("vz",&vz);
	mixing_tree->SetBranchAddress("hiBin",&hiBin);
    mixing_tree->SetBranchAddress("nref",&nref);
    mixing_tree->SetBranchAddress("nPho",&nPho);
    mixing_tree->SetBranchAddress("weight",&pthat_weight);
    if(!is_data) mixing_tree->SetBranchAddress("pthat",&pthat);

    mixing_tree->SetBranchAddress("jtpt",&jtpt);
    mixing_tree->SetBranchAddress("rawpt",&rawpt);
    mixing_tree->SetBranchAddress("jteta",&jteta);
    mixing_tree->SetBranchAddress("jtphi",&jtphi);

    mixing_tree->SetBranchAddress("phoEt",&phoEt);
    mixing_tree->SetBranchAddress("phoEta",&phoEta);
    mixing_tree->SetBranchAddress("phoPhi",&phoPhi);
    mixing_tree->SetBranchAddress("pho_swissCrx",&pho_swissCrx);
    mixing_tree->SetBranchAddress("pho_seedTime",&pho_seedTime);
    mixing_tree->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);
    mixing_tree->SetBranchAddress("phoHoverE",&phoHoverE);
/*
    mixing_tree->SetBranchAddress("pfcIso3pTgt2p0subUE",&pfcIso3pTgt2p0subUE);
    mixing_tree->SetBranchAddress("pfpIso3subUE",&pfpIso3subUE);
    mixing_tree->SetBranchAddress("pfnIso3subUE",&pfnIso3subUE);
*/
    mixing_tree->SetBranchAddress("pho_trackIsoR3PtCut20",&pho_trackIsoR3PtCut20);
    mixing_tree->SetBranchAddress("pho_hcalRechitIsoR3",&pho_hcalRechitIsoR3);
    mixing_tree->SetBranchAddress("pho_ecalClusterIsoR3",&pho_ecalClusterIsoR3);

	n_evt = mixing_tree->GetEntriesFast();
    //n_evt = 3904;
    cout<<n_evt<<endl;
    int ntrigevt = 0;
	
    for(int evi = 0; evi < n_evt; evi++) {
        //cout<<"here0"<<endl;
        if(evi && evi%100000==0) cout << "evi: "<< evi << " of " << n_evt << endl;

        mixing_tree->GetEntry(evi);

        if(is_data) pthat_weight = 1.;
        if(is_pp) hiBin=1;

        h_vz->Fill(vz, pthat_weight);
        h_pthat->Fill(pthat, pthat_weight);
        h_hibin->Fill(hiBin, pthat_weight);

		for (int cbin = 0; cbin < nCbins; cbin++){ 
			if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){
			  mycbin = cbin; 
			}
		}        

        ///photon loop
        double highest_pho_Et = -999.;
        int highest_pho_Et_ind = -1;
        int n_high_pho=0;

        for(int pho = 0; pho < nPho ; pho++) {

            //photon cuts
            if(phoEt->at(pho) < phoEt_min || fabs(phoEta->at(pho)) < eta_min || fabs(phoEta->at(pho)) > eta_max) continue;
            if(is_pp && phoHoverE->at(pho) > hovere_max_pp) continue;
            else if(!is_pp && phoHoverE->at(pho) > hovere_max_PbPb) continue; 

            n_high_pho++;

            if(phoEt->at(pho) > highest_pho_Et){
                highest_pho_Et_ind = pho;
                highest_pho_Et = phoEt->at(pho);
            }
        }//photon loop

        /* require leading photon */
        if (highest_pho_Et_ind < 0.) continue;

        if(is_pp){
            if (phoSigmaIEtaIEta_2012->at(highest_pho_Et_ind) > see_max_pp || phoSigmaIEtaIEta_2012->at(highest_pho_Et_ind) < see_min_pp) continue;
            if ((pho_ecalClusterIsoR3->at(highest_pho_Et_ind)+pho_hcalRechitIsoR3->at(highest_pho_Et_ind)+pho_trackIsoR3PtCut20->at(highest_pho_Et_ind)) > iso_max_pp) continue;
        }
        else{
            if (phoSigmaIEtaIEta_2012->at(highest_pho_Et_ind) > see_max_PbPb || phoSigmaIEtaIEta_2012->at(highest_pho_Et_ind) < see_min_PbPb) continue;
            if ((pho_ecalClusterIsoR3->at(highest_pho_Et_ind)+pho_hcalRechitIsoR3->at(highest_pho_Et_ind)+pho_trackIsoR3PtCut20->at(highest_pho_Et_ind)) > iso_max_PbPb) continue;
        }

        if(highest_pho_Et_ind == -1) continue;
        ntrigevt++;
        
        h_phopt[mycbin]->Fill(phoEt->at(highest_pho_Et_ind), pthat_weight);
        h_phoeta[mycbin]->Fill(fabs(phoEta->at(highest_pho_Et_ind)), pthat_weight);
        h_phophi[mycbin]->Fill(phoPhi->at(highest_pho_Et_ind), pthat_weight);

        //jet loop
        double highest_corrpt = -999.;
        int highest_corrpt_ind = -1;

        for(int jet = 0; jet < nref ; jet++) {
            //MC truth correctons only applicable till 15(20) GeV for PbPb(pp)
            if(rawpt->at(jet) < lowest_corrpt) {
                corrpt[jet] = rawpt->at(jet);
            }
            else{
                if(!is_data){
                    L1.SetJetPT(rawpt->at(jet));
                    L1.SetJetEta(jteta->at(jet));
                    L1.SetJetPhi(jtphi->at(jet));
                    corrpt[jet] = L1.GetCorrectedPT();
                    //corrpt[jet] = rawpt->at(jet);
                }
                else{
                    JEC.SetJetPT(rawpt->at(jet));
                    JEC.SetJetEta(jteta->at(jet));
                    JEC.SetJetPhi(jtphi->at(jet));
                    corrpt[jet] = JEC.GetCorrectedPT();
                }
            }

            jtpt_arr[jet] = corrpt[jet];
            rawpt_arr[jet] = rawpt->at(jet);
            jteta_arr[jet] = jteta->at(jet);
            jtphi_arr[jet] = jtphi->at(jet);

            //continue if jet below 12 GeV and outisde of eta 1.3
            if(corrpt[jet] < 12. || fabs(jteta->at(jet)) < eta_min || fabs(jteta->at(jet)) > eta_max) continue;        

            double dphi3 = fabs(jtphi->at(jet) - phoPhi->at(highest_pho_Et_ind));

            if (dphi3>(2*3.14159)) dphi3-=(2*3.14159);
            if (dphi3>3.14159) dphi3=2*3.14159-dphi3;

            if(dphi3 > 2.7){
                if(corrpt[jet] > highest_corrpt){
                    highest_corrpt_ind = jet;
                    highest_corrpt = corrpt[jet];
                }
            }
        }

        if(highest_corrpt_ind < 0.) continue;
        
        h_dphi[mycbin]->Fill(fabs(phoPhi->at(highest_pho_Et_ind)-jtphi_arr[highest_corrpt_ind]), pthat_weight);

        //find 2nd leading jet
        double sec_highest_jetpt = 0.;
        int sec_highest_jetpt_ind = -1;

        for(int jet = 0; jet < nref ; jet++) {

            //continue if this is leading jet
            if(jet == highest_corrpt_ind || corrpt[jet] < 5. || fabs(jteta->at(jet)) > 5.) continue;        
            
            double dphi3 = fabs(jtphi->at(jet) - phoPhi->at(highest_pho_Et_ind));

            if(dphi3>(2*3.14159)) dphi3-=(2*3.14159);
            if(dphi3>3.14159) dphi3=2*3.14159-dphi3;

            //continue if in the same hemisphere as photon
            if(dphi3 < 1.57) continue;

            if(corrpt[jet] > sec_highest_jetpt){
                sec_highest_jetpt_ind = jet;
                sec_highest_jetpt = corrpt[jet];
            }
        }

        TLorentzVector missEtJets = findMissEtJets(nref, jtpt_arr, rawpt_arr, jteta_arr, jtphi_arr, jtm_arr, lowest_corrpt);
        h_metpt[mycbin]->Fill(missEtJets.Et(), pthat_weight);
        h_meteta[mycbin]->Fill(fabs(missEtJets.Eta()), pthat_weight);
        h_metphi[mycbin]->Fill(missEtJets.Phi(), pthat_weight);

        alpha = sec_highest_jetpt / phoEt->at(highest_pho_Et_ind);

        int alphabin = findAlphaBin(alpha, nbins_alpha, xbins_alpha);

        h_jetpt[mycbin]->Fill(corrpt[highest_corrpt_ind], pthat_weight);
        h_jeteta[mycbin]->Fill(fabs(jteta->at(highest_corrpt_ind)), pthat_weight);
        h_jetphi[mycbin]->Fill(jtphi->at(highest_corrpt_ind), pthat_weight);

        //Fill responses
        //if(corrpt[highest_corrpt_ind] < phoEt->at(highest_pho_Et_ind) && phoEt->at(highest_pho_Et_ind) < 80.) cout<<corrpt[highest_corrpt_ind]<<"  "<<phoEt->at(highest_pho_Et_ind)<<endl;
        //pt bal method
        h_R_jtpt[mycbin][alphabin]->Fill(phoEt->at(highest_pho_Et_ind), corrpt[highest_corrpt_ind] / phoEt->at(highest_pho_Et_ind), pthat_weight);

        //MPF method
        double dphi2 = fabs(missEtJets.Phi() - phoPhi->at(highest_pho_Et_ind));

        if(dphi2>(2*3.14159)) dphi2-=(2*3.14159);
        if (dphi2>3.14159) dphi2=2*3.14159-dphi2;

        h_R_mpf[mycbin][alphabin]->Fill(phoEt->at(highest_pho_Et_ind), (1. + (missEtJets.Et()*cos(dphi2) / phoEt->at(highest_pho_Et_ind))), pthat_weight);

    }//event loop

    cout<<ntrigevt<<endl;
	cout<<"writing"<<endl; 
    TFile *output_file;

    if(is_data){
      if(is_pp) output_file = new TFile("ppdata2017_etaouter_L3Res_Aug9_histos.root", "RECREATE");
      else output_file = new TFile("PbPbdata2018_L3Res_Aug9_histos.root", "RECREATE");       
    }
    else{
      if(is_pp) output_file = new TFile("ppMC2017_etaouter_L3Res_Aug9_histos.root", "RECREATE");
      else output_file = new TFile("PbPbMC2018_L3Res_Aug9_histos.root", "RECREATE");
    }
	output_file->cd();

    h_vz->Write();
    h_pthat->Write();
    h_hibin->Write(); 

 	for(int ibin=0; ibin<nCbins; ibin++){
        h_jetpt[ibin]->Write();
        h_jeteta[ibin]->Write();
        h_jetphi[ibin]->Write();      

        h_metpt[ibin]->Write();
        h_meteta[ibin]->Write();
        h_metphi[ibin]->Write(); 

        h_phopt[ibin]->Write();
        h_phoeta[ibin]->Write();
        h_phophi[ibin]->Write(); 

        h_dphi[ibin]->Write();

        for(int ibin2=0; ibin2<nbins_alpha; ibin2++){
            h_R_jtpt[ibin][ibin2]->Write(); 
            h_R_mpf[ibin][ibin2]->Write(); 
        }
	}
    
	output_file->Close();

	cout<<"done"<<endl;

}

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
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

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

using namespace std;

enum enum_dataset_types {e_Data2015, e_Data_pp, e_HydJet, e_Pythia, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2018","Data_pp","HydJet","Pythia"};

int dataset_type_code = -999;

//arg 1 = which data set, arg 2 = output file number
void make_photonntuples(bool doCrab=0, int jobID=0, int endfile = 199, int dataset_type_code = 3, int output_file_num = 1)
{

    bool is_data = false;

    if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

    cout << "dataset code: " << dataset_type_code << endl;

    //-----------------------------
    // Set JFF-dependent corrections
    //-----------------------------

    float reco_eta, reco_phi, reco_pt, pfPt_temp, pfEta_temp, pfPhi_temp, pfId_temp, pfVsPt_temp;

    double corrected_pt, residual_corrected_pt, r;
    bool do_PbPb=1;

    int radius = 4;

    if(dataset_type_code== 1 || dataset_type_code== 3) do_PbPb = 0;

    string in_file_name;

    if(doCrab){
        in_file_name = Form("job_input_file_list_%d.txt",jobID);
    }
    else if(is_data&&!do_PbPb){
        in_file_name = "pp5TeV_HighPtPD_Apr2016.txt";
    }else if(is_data&&do_PbPb){
        in_file_name = "PbPbData_2015_lxplusskim_MartaAnalysis.txt";
    }else if(dataset_type_code > 10){
        in_file_name = "ppMC_Pythia6_5TeV_partonaxis.txt";
        //in_file_name = "Pythia8_ak4Calo_5TeV.txt";
    }else if(dataset_type_code > 1&&dataset_type_code <11){
        //in_file_name = "pythiaHydjet_2p76TeV_forest.txt";
        in_file_name = "Pythia6_Hydjet_PbPbMix_PurdueList.txt";
    }else{
        cerr<<"need to set up to run on that sample..."<<endl;
    }

    cout << "trying a filelist named "<< in_file_name << endl;

    //MC
    TString output_file_base = "./";

    output_file_base +=dataset_type_strs[dataset_type_code];

    TString output_file_extension = "";   
    //output_file_extension += output_file_num;   
    output_file_extension += ".root";
    TFile *output_file = new TFile((TString) (output_file_base+output_file_extension), "RECREATE");
    TTree *mixing_tree = new TTree("mixing_tree", "");

    const int MAXPARTICLES = 100000;

    Int_t HBHENoiseFilterResultRun2Loose, eventSelection, pprimaryVertexFilter, phfCoincFilter3, pclusterCompatibilityFilter, pBeamScrapingFilter;
    Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_Prescl, HLT_Jet100_Prescl;
    Int_t pPAcollisionEventSelectionPA = -999;

    vector<float> calo_jteta, calo_jtphi, calo_jtpt, calo_corrpt, calo_jtm, calo_rawpt;
    vector<float> leading_phoEt, leading_phoEta, leading_phoPhi;
    vector<float> leading_pho_swissCrx, leading_pho_seedTime, leading_phoSigmaIEtaIEta_2012, leading_phoHoverE, leading_pho_ecalClusterIsoR4, leading_pho_hcalRechitIsoR4, leading_pho_trackIsoR4PtCut20;
    vector<float> calo_refpt, calo_refeta, calo_refphi;
    vector<int> calo_parton_flavor;
    vector<int> trkChg; //Run2
    vector<float> trkEta, trkPhi, trkPt;
    vector<int> sube, chg, pdg;
    vector<float> pt, phi, eta;

    Int_t hiBin = -999;
    Float_t pthat = -999;
    Float_t vz = -999, hiHF= -999.;//, sumpt[15];

    Int_t HLT_HIGEDPhoton20_v1=0, HLT_HIGEDPhoton40_Cent30_100_v1=0, HLT_HIGEDPhoton40_Cent50_100_v1=0; 
    Int_t HLT_HIIslandPhoton40_Eta1p5_v1=0, HLT_HIGEDPhoton40_v1=0, HLT_HIGEDPhoton30_v1=0;
    Int_t HLT_HIGEDPhoton20_v1_Prescl=0, HLT_HIGEDPhoton40_Cent30_100_v1_Prescl=0, HLT_HIGEDPhoton40_Cent50_100_v1_Prescl=0; 
    Int_t HLT_HIIslandPhoton40_Eta1p5_v1_Prescl=0, HLT_HIGEDPhoton40_v1_Prescl=0, HLT_HIGEDPhoton30_v1_Prescl=0;

    Int_t HLT_HISinglePhoton10_Eta1p5ForPPRef_v8, HLT_HISinglePhoton10_Eta1p5ForPPRef_v8_Prescl, HLT_HISinglePhoton15_Eta1p5ForPPRef_v8, HLT_HISinglePhoton15_Eta1p5ForPPRef_v8_Prescl;
    Int_t HLT_HISinglePhoton20_Eta1p5ForPPRef_v8, HLT_HISinglePhoton20_Eta1p5ForPPRef_v8_Prescl, HLT_HISinglePhoton30_Eta1p5ForPPRef_v8, HLT_HISinglePhoton30_Eta1p5ForPPRef_v8_Prescl;
    Int_t HLT_HISinglePhoton40_Eta1p5ForPPRef_v8, HLT_HISinglePhoton40_Eta1p5ForPPRef_v8_Prescl;

    mixing_tree->Branch("vz", &vz);
    if(do_PbPb){
        mixing_tree->Branch("hiBin", &hiBin);
    }
    if(!is_data){
        mixing_tree->Branch("pthat", &pthat);
    }

    if(do_PbPb){
        mixing_tree->Branch("HLT_HIGEDPhoton20_v1", &HLT_HIGEDPhoton20_v1);
        mixing_tree->Branch("HLT_HIGEDPhoton20_v1_Prescl", &HLT_HIGEDPhoton20_v1_Prescl);
        mixing_tree->Branch("HLT_HIGEDPhoton30_v1", &HLT_HIGEDPhoton30_v1);
        mixing_tree->Branch("HLT_HIGEDPhoton30_v1_Prescl", &HLT_HIGEDPhoton30_v1_Prescl);
        mixing_tree->Branch("HLT_HIGEDPhoton40_v1", &HLT_HIGEDPhoton40_v1);
        mixing_tree->Branch("HLT_HIGEDPhoton40_v1_Prescl", &HLT_HIGEDPhoton40_v1_Prescl);
        mixing_tree->Branch("HLT_HIGEDPhoton40_Cent50_100_v1", &HLT_HIGEDPhoton40_Cent50_100_v1);
        mixing_tree->Branch("HLT_HIGEDPhoton40_Cent50_100_v1_Prescl", &HLT_HIGEDPhoton40_Cent50_100_v1_Prescl);
        mixing_tree->Branch("HLT_HIGEDPhoton40_Cent30_100_v1", &HLT_HIGEDPhoton40_Cent30_100_v1);
        mixing_tree->Branch("HLT_HIGEDPhoton40_Cent30_100_v1_Prescl", &HLT_HIGEDPhoton40_Cent30_100_v1_Prescl);
        mixing_tree->Branch("HLT_HIIslandPhoton40_Eta1p5_v1", &HLT_HIIslandPhoton40_Eta1p5_v1);
        mixing_tree->Branch("HLT_HIIslandPhoton40_Eta1p5_v1_Prescl", &HLT_HIIslandPhoton40_Eta1p5_v1_Prescl);

        mixing_tree->Branch("pprimaryVertexFilter", &pprimaryVertexFilter);
        mixing_tree->Branch("collisionEventSelectionAOD", &eventSelection);
        mixing_tree->Branch("phfCoincFilter3Th3", &phfCoincFilter3);
        mixing_tree->Branch("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
        mixing_tree->Branch("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose);
        mixing_tree->Branch("pBeamScrapingFilter", &pBeamScrapingFilter);
    }
    else{
        mixing_tree->Branch("HLT_HISinglePhoton10_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton10_Eta1p5ForPPRef_v8);                        
        mixing_tree->Branch("HLT_HISinglePhoton10_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton10_Eta1p5ForPPRef_v8_Prescl);                        
        mixing_tree->Branch("HLT_HISinglePhoton15_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton15_Eta1p5ForPPRef_v8);                        
        mixing_tree->Branch("HLT_HISinglePhoton15_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton15_Eta1p5ForPPRef_v8_Prescl);                        
        mixing_tree->Branch("HLT_HISinglePhoton20_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton20_Eta1p5ForPPRef_v8);                        
        mixing_tree->Branch("HLT_HISinglePhoton20_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton20_Eta1p5ForPPRef_v8_Prescl);                        
        mixing_tree->Branch("HLT_HISinglePhoton30_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton30_Eta1p5ForPPRef_v8);                        
        mixing_tree->Branch("HLT_HISinglePhoton30_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton30_Eta1p5ForPPRef_v8_Prescl);                        
        mixing_tree->Branch("HLT_HISinglePhoton40_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton40_Eta1p5ForPPRef_v8);                        
        mixing_tree->Branch("HLT_HISinglePhoton40_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton40_Eta1p5ForPPRef_v8_Prescl);                        

        mixing_tree->Branch("pPAprimaryVertexFilter",&pprimaryVertexFilter);
        mixing_tree->Branch("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
        mixing_tree->Branch("pBeamScrapingFilter",&pBeamScrapingFilter);
    } 

    mixing_tree->Branch("calo_jtpt", &calo_jtpt);
    mixing_tree->Branch("calo_rawpt", &calo_rawpt);
    mixing_tree->Branch("calo_jteta", &calo_jteta);
    mixing_tree->Branch("calo_jtphi", &calo_jtphi);
    mixing_tree->Branch("calp_jtm",&calo_jtm);

    if(!is_data){
        mixing_tree->Branch("calo_refpt",&calo_refpt);
    }

    mixing_tree->Branch("phoEt",&leading_phoEt);
    mixing_tree->Branch("phoEta",&leading_phoEta);
    mixing_tree->Branch("phoPhi",&leading_phoPhi);
    mixing_tree->Branch("pho_swissCrx",&leading_pho_swissCrx);
    mixing_tree->Branch("pho_seedTime",&leading_pho_seedTime);
    mixing_tree->Branch("phoSigmaIEtaIEta_2012",&leading_phoSigmaIEtaIEta_2012);
    mixing_tree->Branch("phoHoverE",&leading_phoHoverE);
    mixing_tree->Branch("pho_ecalClusterIsoR4",&leading_pho_ecalClusterIsoR4);
    mixing_tree->Branch("pho_hcalRechitIsoR4",&leading_pho_hcalRechitIsoR4);
    mixing_tree->Branch("pho_trackIsoR4PtCut20",&leading_pho_trackIsoR4PtCut20);

    std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
    if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
    std::string filename;
    int ifile=0;

    const int MAXJETS = 500;

    Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS], t_calo_rawpt[MAXJETS], t_calo_jtm[MAXJETS];
    Float_t t_calo_refpt[MAXJETS], t_calo_refeta[MAXJETS], t_calo_refphi[MAXJETS];

    Int_t n_evt, nTrk, mult, calo_nref, nPFpart, nCSpart, nPho;
    
    while(instr>>filename && ifile<endfile){
        filename.erase(std::remove(filename.begin(), filename.end(), '"'), filename.end());
        filename.erase(std::remove(filename.begin(), filename.end(), ','), filename.end());
        filename.erase(std::remove(filename.begin(), filename.end(), '['), filename.end());
        filename.erase(std::remove(filename.begin(), filename.end(), ']'), filename.end());
        cout<<"File name is "<< filename <<endl;
        ifile++;

        //TFile *my_file = TFile::Open(filename.c_str());

        //if(!my_file){
            int pos = filename.find_first_of('s');
            string reducedfn = filename.substr(pos-1);
            string xrdPrefix = "root://cmsxrootd.fnal.gov//";
            cout << "local file not detected. Trying " << xrdPrefix+reducedfn << endl;
            TFile *my_file = TFile::Open((xrdPrefix+reducedfn).c_str());
                        //TFile::Open((xrdPrefix+reducedfn).c_str());
                //if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }    
        //}
      

        //TFile *my_file = TFile::Open("root://cmsxrootd.fnal.gov///store/group/phys_heavyions/stepobr/pp2017/HighEGJet/crab_20190705_002317/190704_222359/0001/HiForestAOD_1472.root"); 
        //TFile *my_file = TFile::Open("root://cmsxrootd.fnal.gov///store/group/phys_heavyions/dhangal/pp2017MCForest_5TeV_AllQCDPhoton_pthat170/QCDPhoton_pThat-170_TuneCP5_5p02TeV_pythia8/crab_pp2017MCForest_5TeV_AllQCDPhoton_pthat170/190703_173407/0000/HiForestAOD_44.root"); 


        if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }  

        if(my_file->IsZombie()) { 
            cout << "Is zombie" << std::endl;
            continue;
        }   

        TTree *akCs4PFtree; 
        TTree *skimanalysis_tree; 
        TTree *hlt_tree; 
        TTree *photon_tree; 
        TTree *event_tree;
        TTree *pfcand_tree;

        event_tree = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
        if(!do_PbPb && is_data)hlt_tree = (TTree*) my_file->Get("hltanalysisReco/HltTree");
        else hlt_tree = (TTree*) my_file->Get("hltanalysis/HltTree");
        akCs4PFtree = (TTree*) my_file->Get("ak4PFJetAnalyzer/t");
        skimanalysis_tree = (TTree*) my_file->Get("skimanalysis/HltTree");
        if(do_PbPb) photon_tree = (TTree*) my_file->Get("ggHiNtuplizer/EventTree");
        else photon_tree = (TTree*) my_file->Get("ggHiNtuplizerGED/EventTree");

        UInt_t run=-999, lumi = -999;

        //vector<float> *phoEt=0,*phoEta=0,*phoPhi=0;
        vector<float> *phoEta=0, *phoPhi=0, *phoEt=0, *pho_swissCrx=0, *pho_seedTime=0, *phoSigmaIEtaIEta_2012=0, *phoHoverE=0, *pho_ecalClusterIsoR4=0, *pho_hcalRechitIsoR4=0, *pho_trackIsoR4PtCut20=0;
        vector<float> *pfcand_pt=0, *pfcand_eta=0, *pfcand_phi=0;
        vector<int> *pfcand_ID = 0;

        double alpha;

        event_tree->SetBranchAddress("run",&run);
        event_tree->SetBranchAddress("lumi",&lumi);
        event_tree->SetBranchAddress("vz",&vz);
        if(do_PbPb) {
            event_tree->SetBranchAddress("hiBin",&hiBin);
            event_tree->SetBranchAddress("hiHF",&hiHF);
        }
        if(!is_data) event_tree->SetBranchAddress("pthat",&pthat);

        if(do_PbPb){
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton20_v1",&HLT_HIGEDPhoton20_v1);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton20_v1_Prescl",&HLT_HIGEDPhoton20_v1_Prescl);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton30_v1",&HLT_HIGEDPhoton30_v1);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton30_v1_Prescl",&HLT_HIGEDPhoton30_v1_Prescl);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton40_v1",&HLT_HIGEDPhoton40_v1);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton40_v1_Prescl",&HLT_HIGEDPhoton40_v1_Prescl);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton40_Cent50_100_v1",&HLT_HIGEDPhoton40_Cent50_100_v1);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton40_Cent50_100_v1_Prescl",&HLT_HIGEDPhoton40_Cent50_100_v1_Prescl);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton40_Cent30_100_v1",&HLT_HIGEDPhoton40_Cent30_100_v1);            
          hlt_tree->SetBranchAddress("HLT_HIGEDPhoton40_Cent30_100_v1_Prescl",&HLT_HIGEDPhoton40_Cent30_100_v1_Prescl);            
          hlt_tree->SetBranchAddress("HLT_HIIslandPhoton40_Eta1p5_v1",&HLT_HIIslandPhoton40_Eta1p5_v1);            
          hlt_tree->SetBranchAddress("HLT_HIIslandPhoton40_Eta1p5_v1_Prescl",&HLT_HIIslandPhoton40_Eta1p5_v1_Prescl);            

          skimanalysis_tree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
          skimanalysis_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
          skimanalysis_tree->SetBranchAddress("collisionEventSelectionAOD",&eventSelection);
          skimanalysis_tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
          skimanalysis_tree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
          skimanalysis_tree->SetBranchAddress("phfCoincFilter3Th3",&phfCoincFilter3);
        }
        else{
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton10_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton10_Eta1p5ForPPRef_v8);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton10_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton10_Eta1p5ForPPRef_v8_Prescl);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton15_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton15_Eta1p5ForPPRef_v8);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton15_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton15_Eta1p5ForPPRef_v8_Prescl);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton20_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton20_Eta1p5ForPPRef_v8);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton20_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton20_Eta1p5ForPPRef_v8_Prescl);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton30_Eta1p5ForPPRef_v8);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton30_Eta1p5ForPPRef_v8_Prescl);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton40_Eta1p5ForPPRef_v8",&HLT_HISinglePhoton40_Eta1p5ForPPRef_v8);                        
          hlt_tree->SetBranchAddress("HLT_HISinglePhoton40_Eta1p5ForPPRef_v8_Prescl",&HLT_HISinglePhoton40_Eta1p5ForPPRef_v8_Prescl);                        

          skimanalysis_tree->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryVertexFilter);
          skimanalysis_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
          skimanalysis_tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
        } 

        akCs4PFtree->SetBranchAddress("nref",&calo_nref);
        akCs4PFtree->SetBranchAddress("jtpt",t_calo_jtpt);
        akCs4PFtree->SetBranchAddress("rawpt",t_calo_rawpt);
        akCs4PFtree->SetBranchAddress("jteta",t_calo_jteta);
        akCs4PFtree->SetBranchAddress("jtphi",t_calo_jtphi);
        akCs4PFtree->SetBranchAddress("jtm",t_calo_jtm);
        if(!is_data) akCs4PFtree->SetBranchAddress("refpt",t_calo_refpt);

        photon_tree->SetBranchAddress("nPho",&nPho);
        photon_tree->SetBranchAddress("phoEt",&phoEt);
        photon_tree->SetBranchAddress("phoEta",&phoEta);
        photon_tree->SetBranchAddress("phoPhi",&phoPhi);
        photon_tree->SetBranchAddress("pho_swissCrx",&pho_swissCrx);
        photon_tree->SetBranchAddress("pho_seedTime",&pho_seedTime);
        photon_tree->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);
        photon_tree->SetBranchAddress("phoHoverE",&phoHoverE);
        photon_tree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
        photon_tree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
        photon_tree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);

        n_evt = event_tree->GetEntriesFast();

        //n_evt=10;
        cout << "Entries: "<< n_evt << endl;
        for(int evi = 0; evi < n_evt; evi++) {

            if(evi && evi%1000==0) cout << "evi: "<< evi << " of " << n_evt << endl;

            event_tree->GetEntry(evi);
            skimanalysis_tree->GetEntry(evi);
            hlt_tree->GetEntry(evi);
            akCs4PFtree->GetEntry(evi);
            photon_tree->GetEntry(evi);

            if(fabs(vz) >= 15.) continue;

            bool found_photon = false;

            //photon loop
            for(int pho = 0; pho < (int) phoEt->size() ; pho++) {
                if(phoEt->at(pho) < 40. || fabs(phoEta->at(pho)) > 1.3) continue;                 
                //else if(abs(pho_seedTime->at(pho)) > 3. || pho_swissCrx->at(pho) >0.9 || phoSigmaIEtaIEta_2012->at(pho) < 0.002 || phoSigmaIEtaIEta_2012->at(pho) > 0.01) continue;
                //else if (phoHoverE->at(pho) > 0.1) continue;
                //else if ((pho_ecalClusterIsoR4->at(pho) + pho_hcalRechitIsoR4->at(pho) + pho_trackIsoR4PtCut20->at(pho)) > 1) continue;

                found_photon = true;

                leading_phoEt.push_back(phoEt->at(pho));
                leading_phoEta.push_back(phoEta->at(pho));
                leading_phoPhi.push_back(phoPhi->at(pho));
                leading_phoHoverE.push_back(phoHoverE->at(pho));
                leading_pho_seedTime.push_back(pho_seedTime->at(pho));
                leading_pho_swissCrx.push_back(pho_swissCrx->at(pho));
                leading_pho_hcalRechitIsoR4.push_back(pho_hcalRechitIsoR4->at(pho));
                leading_pho_ecalClusterIsoR4.push_back(pho_ecalClusterIsoR4->at(pho));
                leading_phoSigmaIEtaIEta_2012.push_back(phoSigmaIEtaIEta_2012->at(pho));
                leading_pho_trackIsoR4PtCut20.push_back(pho_trackIsoR4PtCut20->at(pho));
            }

            if(found_photon==0) continue; 

            //start calo jet loop
            for(int j4i = 0; j4i < calo_nref ; j4i++) {

                if(fabs(t_calo_jteta[j4i]) > 5. || t_calo_rawpt[j4i] < 5.) continue;

                calo_jteta.push_back(t_calo_jteta[j4i]);
                calo_jtphi.push_back(t_calo_jtphi[j4i]);
                calo_jtpt.push_back(t_calo_jtpt[j4i]);
                calo_rawpt.push_back(t_calo_rawpt[j4i]);
                calo_jtm.push_back(t_calo_jtm[j4i]);
                if(!is_data){ 
                    calo_refpt.push_back(t_calo_refpt[j4i]);
                }
            } /// calo jet loop

            ///// Fill it
            mixing_tree->Fill();

            calo_jteta.clear();
            calo_jtm.clear();
            calo_jtphi.clear();
            calo_jtpt.clear();
            calo_rawpt.clear();

            leading_phoEt.clear();
            leading_phoEta.clear();
            leading_phoPhi.clear();

            leading_pho_swissCrx.clear();
            leading_pho_seedTime.clear();
            leading_pho_ecalClusterIsoR4.clear();
            leading_pho_hcalRechitIsoR4.clear();
            leading_pho_trackIsoR4PtCut20.clear();
            leading_phoSigmaIEtaIEta_2012.clear();
            leading_phoHoverE.clear();

            if(!is_data){
                calo_refpt.clear();
            }

        }  ///event loop
        my_file->Close();

    }//file loop

    cout<<"writing"<<endl;

    output_file->cd();
    mixing_tree->Write();
    output_file->Close();

    cout<<"done"<<endl;

}

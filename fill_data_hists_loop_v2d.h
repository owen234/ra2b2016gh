//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 22 16:30:51 2016 by ROOT version 6.05/02
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef fill_data_hists_loop_v2d_h
#define fill_data_hists_loop_v2d_h

#include "TLorentzVector.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1F.h"

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class fill_data_hists_loop_v2d {
public :

   //----------
   void set_bi() ;

   void translate_global_bin( int gbi, int& tbi_nj, int& tbi_nb, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) ;
   void translate_global_bin_nbsum( int gbi_nbsum, int& tbi_nj, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) ;
   int  ht_bi_from_htmht( int bi_htmht ) ;
   int  mht_bi_from_htmht( int bi_htmht ) ;

   void set_bin_labels( TH1F* hp ) ;
   void set_bin_labels_div_by_nb( TH1F* hp ) ;
   void set_bin_labels_mhtc_plot( TH1F* hp ) ;

   char samplename[1000] ;

   //---------

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          RunNum;
   UInt_t          LumiBlockNum;
   ULong64_t       EvtNum;
   Bool_t          BadChargedCandidateFilter;
   Bool_t          BadPFMuonFilter;
   Int_t           BTags;
   Double_t        CaloMET;
   Double_t        CaloMETPhi;
   Int_t           CSCTightHaloFilter;
   Double_t        DeltaPhi1;
   Double_t        DeltaPhi2;
   Double_t        DeltaPhi3;
   Double_t        DeltaPhi4;
   Int_t           EcalDeadCellTriggerPrimitiveFilter;
   Int_t           eeBadScFilter;
   Int_t           globalTightHalo2016Filter;
   Int_t           HBHEIsoNoiseFilter;
   Int_t           HBHENoiseFilter;
   Double_t        HT;
   Bool_t          JetID;
   vector<TLorentzVector> *Jets;
   vector<int>     *Jets_hadronFlavor;
   vector<double>  *Jets_muonEnergyFraction;
   vector<int>     *Jets_partonFlavor;
   vector<double>  *Jets_qgAxis2;
   vector<double>  *Jets_qgLikelihood;
   vector<int>     *Jets_qgMult;
   vector<double>  *Jets_qgPtD;
   Double_t        MET;
   Double_t        METPhi;
   Double_t        MHT;
   Double_t        MHTPhi;
   Int_t           NJets;
   Int_t           NVtx;
   Double_t        PFCaloMETRatio;
   vector<string>  *TriggerNames;
   vector<int>     *TriggerPass;
   Char_t          passBadMuFilter;
   Char_t          noMuonJet;
   Double_t        minDeltaPhi;
   Char_t          inLDP;
   Char_t          hasHadTau;
   Char_t          pass_pfmet100_trig;
   Char_t          pass_pfmetnomu100_trig;
   Char_t          pass_pfmet110_trig;
   Char_t          pass_pfmetnomu110_trig;
   Char_t          pass_pfmet120_trig;
   Char_t          pass_pfmetnomu120_trig;
   Char_t          pass_ht300_met100_trig;
   Char_t          pass_ht800_trig;
   Char_t          pass_ht900_trig;

   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_LumiBlockNum;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_BadChargedCandidateFilter;   //!
   TBranch        *b_BadPFMuonFilter;   //!
   TBranch        *b_BTags;   //!
   TBranch        *b_CaloMET;   //!
   TBranch        *b_CaloMETPhi;   //!
   TBranch        *b_CSCTightHaloFilter;   //!
   TBranch        *b_DeltaPhi1;   //!
   TBranch        *b_DeltaPhi2;   //!
   TBranch        *b_DeltaPhi3;   //!
   TBranch        *b_DeltaPhi4;   //!
   TBranch        *b_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_eeBadScFilter;   //!
   TBranch        *b_globalTightHalo2016Filter;   //!
   TBranch        *b_HBHEIsoNoiseFilter;   //!
   TBranch        *b_HBHENoiseFilter;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_JetID;   //!
   TBranch        *b_Jets;   //!
   TBranch        *b_Jets_hadronFlavor;   //!
   TBranch        *b_Jets_muonEnergyFraction;   //!
   TBranch        *b_Jets_partonFlavor;   //!
   TBranch        *b_Jets_qgAxis2;   //!
   TBranch        *b_Jets_qgLikelihood;   //!
   TBranch        *b_Jets_qgMult;   //!
   TBranch        *b_Jets_qgPtD;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTPhi;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_PFCaloMETRatio;   //!
   TBranch        *b_TriggerNames;   //!
   TBranch        *b_TriggerPass;   //!
   TBranch        *b_passBadMuFilter;   //!
   TBranch        *b_noMuonJet;   //!
   TBranch        *b_minDeltaPhi;   //!
   TBranch        *b_inLDP;   //!
   TBranch        *b_hasHadTau;   //!
   TBranch        *b_pass_pfmet100_trig;   //!
   TBranch        *b_pass_pfmetnomu100_trig;   //!
   TBranch        *b_pass_pfmet110_trig;   //!
   TBranch        *b_pass_pfmetnomu110_trig;   //!
   TBranch        *b_pass_pfmet120_trig;   //!
   TBranch        *b_pass_pfmetnomu120_trig;   //!
   TBranch        *b_pass_ht300_met100_trig;   //!
   TBranch        *b_pass_ht800_trig;   //!
   TBranch        *b_pass_ht900_trig;   //!

   fill_data_hists_loop_v2d(TTree *tree=0, const char* samplename_arg="" );
   virtual ~fill_data_hists_loop_v2d();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool verb=false, int nloop=-1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fill_data_hists_loop_v2d_cxx
fill_data_hists_loop_v2d::fill_data_hists_loop_v2d(TTree *tree, const char* samplename_arg ) : fChain(0) 
{
   sprintf( samplename, "%s", samplename_arg ) ;
   printf("\n\n Creating loop for sample %s\n\n", samplename ) ;

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
     //----------
      chain->Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016B-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016C-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016D-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016E-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016F-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016G-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_signal/tree_MET_2016B-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_signal/tree_MET_2016C-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_signal/tree_MET_2016D-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_signal/tree_MET_2016E-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_signal/tree_MET_2016F-slimskim.root/tree");
      chain->Add("fnal-prod-v10-skims-slimmed/tree_signal/tree_MET_2016G-slimskim.root/tree");
     //----------
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

fill_data_hists_loop_v2d::~fill_data_hists_loop_v2d()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fill_data_hists_loop_v2d::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fill_data_hists_loop_v2d::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fill_data_hists_loop_v2d::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Jets = 0;
   Jets_hadronFlavor = 0;
   Jets_muonEnergyFraction = 0;
   Jets_partonFlavor = 0;
   Jets_qgAxis2 = 0;
   Jets_qgLikelihood = 0;
   Jets_qgMult = 0;
   Jets_qgPtD = 0;
   TriggerNames = 0;
   TriggerPass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("BadChargedCandidateFilter", &BadChargedCandidateFilter, &b_BadChargedCandidateFilter);
   fChain->SetBranchAddress("BadPFMuonFilter", &BadPFMuonFilter, &b_BadPFMuonFilter);
   fChain->SetBranchAddress("BTags", &BTags, &b_BTags);
   fChain->SetBranchAddress("CaloMET", &CaloMET, &b_CaloMET);
   fChain->SetBranchAddress("CaloMETPhi", &CaloMETPhi, &b_CaloMETPhi);
   fChain->SetBranchAddress("CSCTightHaloFilter", &CSCTightHaloFilter, &b_CSCTightHaloFilter);
   fChain->SetBranchAddress("DeltaPhi1", &DeltaPhi1, &b_DeltaPhi1);
   fChain->SetBranchAddress("DeltaPhi2", &DeltaPhi2, &b_DeltaPhi2);
   fChain->SetBranchAddress("DeltaPhi3", &DeltaPhi3, &b_DeltaPhi3);
   fChain->SetBranchAddress("DeltaPhi4", &DeltaPhi4, &b_DeltaPhi4);
   fChain->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, &b_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("eeBadScFilter", &eeBadScFilter, &b_eeBadScFilter);
   fChain->SetBranchAddress("globalTightHalo2016Filter", &globalTightHalo2016Filter, &b_globalTightHalo2016Filter);
   fChain->SetBranchAddress("HBHEIsoNoiseFilter", &HBHEIsoNoiseFilter, &b_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("HBHENoiseFilter", &HBHENoiseFilter, &b_HBHENoiseFilter);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("JetID", &JetID, &b_JetID);
   fChain->SetBranchAddress("Jets", &Jets, &b_Jets);
   fChain->SetBranchAddress("Jets_hadronFlavor", &Jets_hadronFlavor, &b_Jets_hadronFlavor);
   fChain->SetBranchAddress("Jets_muonEnergyFraction", &Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction);
   fChain->SetBranchAddress("Jets_partonFlavor", &Jets_partonFlavor, &b_Jets_partonFlavor);
   fChain->SetBranchAddress("Jets_qgAxis2", &Jets_qgAxis2, &b_Jets_qgAxis2);
   fChain->SetBranchAddress("Jets_qgLikelihood", &Jets_qgLikelihood, &b_Jets_qgLikelihood);
   fChain->SetBranchAddress("Jets_qgMult", &Jets_qgMult, &b_Jets_qgMult);
   fChain->SetBranchAddress("Jets_qgPtD", &Jets_qgPtD, &b_Jets_qgPtD);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTPhi", &MHTPhi, &b_MHTPhi);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("PFCaloMETRatio", &PFCaloMETRatio, &b_PFCaloMETRatio);
   fChain->SetBranchAddress("TriggerNames", &TriggerNames, &b_TriggerNames);
   fChain->SetBranchAddress("TriggerPass", &TriggerPass, &b_TriggerPass);
   fChain->SetBranchAddress("passBadMuFilter", &passBadMuFilter, &b_passBadMuFilter);
   fChain->SetBranchAddress("noMuonJet", &noMuonJet, &b_noMuonJet);
   fChain->SetBranchAddress("minDeltaPhi", &minDeltaPhi, &b_minDeltaPhi);
   fChain->SetBranchAddress("inLDP", &inLDP, &b_inLDP);
   fChain->SetBranchAddress("hasHadTau", &hasHadTau, &b_hasHadTau);
   fChain->SetBranchAddress("pass_pfmet100_trig", &pass_pfmet100_trig, &b_pass_pfmet100_trig);
   fChain->SetBranchAddress("pass_pfmetnomu100_trig", &pass_pfmetnomu100_trig, &b_pass_pfmetnomu100_trig);
   fChain->SetBranchAddress("pass_pfmet110_trig", &pass_pfmet110_trig, &b_pass_pfmet110_trig);
   fChain->SetBranchAddress("pass_pfmetnomu110_trig", &pass_pfmetnomu110_trig, &b_pass_pfmetnomu110_trig);
   fChain->SetBranchAddress("pass_pfmet120_trig", &pass_pfmet120_trig, &b_pass_pfmet120_trig);
   fChain->SetBranchAddress("pass_pfmetnomu120_trig", &pass_pfmetnomu120_trig, &b_pass_pfmetnomu120_trig);
   fChain->SetBranchAddress("pass_ht300_met100_trig", &pass_ht300_met100_trig, &b_pass_ht300_met100_trig);
   fChain->SetBranchAddress("pass_ht800_trig", &pass_ht800_trig, &b_pass_ht800_trig);
   fChain->SetBranchAddress("pass_ht900_trig", &pass_ht900_trig, &b_pass_ht900_trig);
   Notify();
}

Bool_t fill_data_hists_loop_v2d::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fill_data_hists_loop_v2d::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fill_data_hists_loop_v2d::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   if (entry > 0 ) return 1 ;
   return 1;
}
#endif // #ifdef fill_data_hists_loop_v2d_cxx

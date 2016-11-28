#ifndef data_turnon1_c
#define data_turnon1_c

#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include <string>
#include "binning.h"
#include "histio.c"
#include "TLorentzVector.h"


   void add_overflow( TH1F* hp ) ;

  //---------

   void data_turnon1(  ) {

      setup_bins();

      TChain ch_ht("tree") ;
      TChain ch_met("tree") ;

      TCanvas* can = new TCanvas( "can_data_turnon1", "Data trigger efficiency", 900, 800 ) ;

     //-----------
   // ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016B-slimskim.root") ;
   // ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016C-slimskim.root") ;
   // ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016D-slimskim.root") ;
   // ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016B-slimskim.root") ;
   // ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016C-slimskim.root") ;
   // ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016D-slimskim.root") ;
     //-----------
      ch_ht.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_JetHT_2016B-slimskim.root") ;
      ch_ht.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_JetHT_2016C-slimskim.root") ;
      ch_ht.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_JetHT_2016D-slimskim.root") ;
      ch_ht.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_JetHT_2016E-slimskim.root") ;
      ch_ht.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_JetHT_2016F-slimskim.root") ;
      ch_ht.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_JetHT_2016G-slimskim.root") ;
      ch_met.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016B-slimskim.root") ;
      ch_met.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016C-slimskim.root") ;
      ch_met.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016D-slimskim.root") ;
      ch_met.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016E-slimskim.root") ;
      ch_met.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016F-slimskim.root") ;
      ch_met.Add("fnal-prod-v10-skims-slimmed/tree_LDP/tree_MET_2016G-slimskim.root") ;
     //-----------


      int nbin(10) ;
      double xbins[11] = {0., 200., 220, 240, 260, 280, 300, 340, 400, 500, 700 } ;

      TH1F* denom_h_mht_nj[10];
      TH1F* numer_h_mht_nj[10];
      TH1F* h_eff_nj[10];

      for ( int nj = 0; nj < nb_nj; nj++)
      {
         TString nj_str         ; nj_str         .Form("%d",nj+1);

         denom_h_mht_nj[nj] = new TH1F( "h_mht_nj"+nj_str+"_denom", "MHT, Nj"+nj_str+", HT800" ,  nbin, xbins ) ;
         numer_h_mht_nj[nj] = new TH1F( "h_mht_nj"+nj_str+"_numer", "MHT, Nj"+nj_str+", MET100",  nbin, xbins ) ;
         h_eff_nj[nj]       = new TH1F( "h_eff_nj"+nj_str         , "Eff, Nj"+nj_str           , nbin, xbins ) ;
      }//nj

         int NJets=0 ;
         double HT=0, MHT=0 ;
         double CaloMET = 0, MET = 0, PFCaloMETRatio = 0;
         std::vector < TLorentzVector > *Jets = 0;
         Int_t eeBadScFilter = 0, HBHENoiseFilter = 0, EcalDeadCellTriggerPrimitiveFilter = 0, NVtx = 0, HBHEIsoNoiseFilter = 0 , BTags;
         Bool_t JetID = 0, pass_ht800_trig = 0, pass_pfmet100_trig = 0, pass_pfmetnomu100_trig = 0, pass_pfmet110_trig = 0, pass_pfmetnomu110_trig = 0, pass_pfmet120_trig = 0, pass_pfmetnomu120_trig = 0;



         ch_ht.SetBranchAddress("HT",&HT);
         ch_ht.SetBranchAddress("NJets",&NJets);
         ch_ht.SetBranchAddress("BTags",&BTags);
         ch_ht.SetBranchAddress("MHT",&MHT);
         ch_ht.SetBranchAddress("PFCaloMETRatio",&PFCaloMETRatio);
         ch_ht.SetBranchAddress("MET",&MET);
         ch_ht.SetBranchAddress("CaloMET",&CaloMET);
         ch_ht.SetBranchAddress("eeBadScFilter", &eeBadScFilter);
         ch_ht.SetBranchAddress("HBHENoiseFilter", &HBHENoiseFilter);
         ch_ht.SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter);
         ch_ht.SetBranchAddress("NVtx", &NVtx);
         ch_ht.SetBranchAddress("JetID",&JetID);
         ch_ht.SetBranchAddress("HBHEIsoNoiseFilter", &HBHEIsoNoiseFilter);
         ch_ht.SetBranchAddress("pass_ht800_trig", &pass_ht800_trig);
         ch_ht.SetBranchAddress("pass_pfmet100_trig", &pass_pfmet100_trig);
         ch_ht.SetBranchAddress("pass_pfmetnomu100_trig", &pass_pfmetnomu100_trig);
         ch_ht.SetBranchAddress("pass_pfmet110_trig", &pass_pfmet110_trig);
         ch_ht.SetBranchAddress("pass_pfmetnomu110_trig", &pass_pfmetnomu110_trig);
         ch_ht.SetBranchAddress("pass_pfmet120_trig", &pass_pfmet120_trig);
         ch_ht.SetBranchAddress("pass_pfmetnomu120_trig", &pass_pfmetnomu120_trig);
   
	 unsigned int nentreis = ch_ht.GetEntries();
	 std::cout << "number of entries: " << nentreis << std::endl;
      for ( unsigned int i = 0; i < nentreis; i++ )
      {

         ch_ht.GetEntry(i);

      	 if ( !HBHEIsoNoiseFilter                              ) continue;
         if ( !HBHENoiseFilter                                 ) continue;
         if ( !eeBadScFilter                                   ) continue;
         if ( !EcalDeadCellTriggerPrimitiveFilter              ) continue;
         if ( !JetID                                           ) continue;
         if ( !pass_ht800_trig                                 ) continue;
         if ( NVtx <= 0                                        ) continue;
         if ( HT <= 900.                                       ) continue;
         if ( !(MET>100 && CaloMET>80 && PFCaloMETRatio<5)     ) continue;


         for ( int nj = 0; nj < nb_nj; nj++)
         {
            if ( !(NJets > bin_edges_nj[nj] && NJets < bin_edges_nj[nj+1] ) ) continue;
	 
            denom_h_mht_nj[nj]->Fill(MHT);

	    //applying pfmet filters

            if ( !(pass_pfmet100_trig || pass_pfmetnomu100_trig || pass_pfmet110_trig || pass_pfmetnomu110_trig || pass_pfmet120_trig || pass_pfmetnomu120_trig ) ) continue;

            numer_h_mht_nj[nj]->Fill(MHT);

         }//nj
      }//i
      for ( int nj = 0; nj < nb_nj; nj++)
      {

         denom_h_mht_nj[nj]->Draw("SAME");
	 numer_h_mht_nj[nj]->Draw("SAME");

         add_overflow( denom_h_mht_nj[nj] ) ;
         add_overflow( numer_h_mht_nj[nj] ) ;

         for ( int bi=1; bi<=nbin; bi++ ) {

            float eff_val(0.) ;
            float eff_err(0.) ;
            float denom, numer ;

            denom = denom_h_mht_nj[nj] -> GetBinContent( bi ) ;
            numer = numer_h_mht_nj[nj] -> GetBinContent( bi ) ;
            if (denom>0) {
               eff_val = numer/denom ;
               eff_err = sqrt( eff_val*(1.-eff_val)/denom ) ;
            }
            h_eff_nj[nj] -> SetBinContent( bi, eff_val ) ;
            h_eff_nj[nj] -> SetBinError( bi, eff_err ) ;


         } // bi
 
         gStyle -> SetOptStat(0) ;

         h_eff_nj[nj] -> SetMaximum( 1.05 ) ;
         h_eff_nj[nj] -> SetMinimum( 0.4 ) ;

         h_eff_nj[nj] -> SetMarkerStyle(20+nj) ;
         h_eff_nj[nj] -> SetMarkerColor(2+nj) ;
      }//nj

      h_eff_nj[0] -> Draw() ;
      gPad -> SetGridx(1) ;
      gPad -> SetGridy(1) ;
      for ( int nj = 1; nj < nb_nj; nj++) {
         h_eff_nj[nj] -> Draw("same") ;
      }
      can -> SaveAs( "outputfiles/data-turnon1.pdf" ) ;


      gSystem -> Exec( "mkdir -p outputfiles" ) ;

      saveHist("outputfiles/data-turnon.root", "h*") ;


   } // data_turnon1

 //=========================================================

   void add_overflow( TH1F* hp ) {

      float last_val = hp -> GetBinContent( hp -> GetNbinsX() ) ;
      float overflow_val = hp -> GetBinContent( hp -> GetNbinsX() + 1 ) ;

      hp -> SetBinContent( hp -> GetNbinsX(), last_val + overflow_val ) ;

   } // add_overflow


 //=========================================================

#endif

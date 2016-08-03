#include "TChain.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include <string>

#ifndef binning_h
#include "binning.h"
#endif

#ifndef histio_c
#include "histio.c"
#endif

   void add_overflow( TH1F* hp ) ;

  //---------

   void data_turnon1(  ) {

      setup_bins();

      TChain ch_ht("tree") ;
      TChain ch_met("tree") ;

     //-----------
      ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016B-slimskim.root") ;
      ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016C-slimskim.root") ;
      ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016D-slimskim.root") ;
      ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016B-slimskim.root") ;
      ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016C-slimskim.root") ;
      ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016D-slimskim.root") ;
     //-----------


      int nbin(10) ;
      double xbins[11] = {0., 200., 220, 240, 260, 280, 300, 340, 400, 500, 700 } ;

      TH1F* denom_h_mht_nj[10];
      TH1F* numer_h_mht_nj[10];
      TH1F* h_eff_nj[10];

      for ( int nj = 0; nj < nb_nj; nj++)
      {

         TString nj_str         ; nj_str         .Form("%d",nj+1); 
         TString nj_cut_low_str ; nj_cut_low_str .Form("%f",bin_edges_nj[nj]);
         TString nj_cut_high_str; nj_cut_high_str.Form("%f",bin_edges_nj[nj+1]);

         denom_h_mht_nj[nj] = new TH1F( "h_mht_nj"+nj_str+"_denom", "MHT, Nj"+nj_str+", HT800" ,  nbin, xbins ) ;
         numer_h_mht_nj[nj] = new TH1F( "h_mht_nj"+nj_str+"_numer", "MHT, Nj"+nj_str+", MET100",  nbin, xbins ) ;
         h_eff_nj[nj]       = new TH1F( "h_eff_nj"+nj_str         , "Eff, Nj"+nj_str           , nbin, xbins ) ;

         ch_ht.Draw("MHT>>h_mht_nj"+nj_str+"_denom","HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID &&  HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig && NJets >= "+ nj_cut_low_str + " && NJets < " + nj_cut_high_str) ;
         ch_met.Draw("MHT>>h_mht_nj"+nj_str+"_numer", "HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig && (pass_pfmet100_trig||pass_pfmetnomu100_trig) && NJets >= "+ nj_cut_low_str + " && NJets < " + nj_cut_high_str ) ;

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
         h_eff_nj[nj] -> SetMinimum( 0.5 ) ;

         h_eff_nj[nj] -> SetMarkerStyle(20+nj) ;
         h_eff_nj[nj] -> SetMarkerColor(2+nj) ;
      }//nj

      h_eff_nj[0] -> Draw() ;
      gPad -> SetGridx(1) ;
      gPad -> SetGridy(1) ;
      h_eff_nj[1] -> Draw("same") ;
      h_eff_nj[2] -> Draw("same") ;
      h_eff_nj[3] -> Draw("same") ;


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



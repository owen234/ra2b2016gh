//
//  Organizing principles of binning.
//
//    For each dimension
//      bin index 0 means below bins used in analysis (e.g. HT<500 is HT bin 0).
//      bin index 1 is the first bin used
//      This also goes for btagging, so btag bin 1 is nb=0.
//      In the bin edges arrays, it goes like this
//        bin_edge_x[0] = 500. ; // low edge of bin 1.
//        bin_edge_x[1] = 750. ; // high edge of bin 1, low edge of bin 2.
//        ...
//

#ifndef fill_data_hists_loop_v2d_cxx
#define fill_data_hists_loop_v2d_cxx
#include "fill_data_hists_loop_v2d.h"
#include <TH2.h>
#include "TH1F.h"
#include "TH2F.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TMath.h"
#include <iostream>
#include <vector>

#include "lumi_taken.h"
#include "binning.h"
#include "histio.c"

using namespace std ;


void fill_data_hists_loop_v2d::Loop( bool verb, int nloop )
{

   //bool blind(true);
   bool blind(false);

   setup_bins() ;

   gDirectory -> Delete( "h*" ) ;

   if (fChain == 0) return;

   TH1F* h_hdp_nb[10], *h_ldp_nb[10];

   TH1::SetDefaultSumw2(); // to make sure all TH1 histograms have Sumw2 enabled

   for ( int nb = 0; nb < nb_nb; nb++)
   {

      TString nb_str         ; nb_str         .Form("%d",nb);
      h_hdp_nb[nb] = new TH1F( "h_hdp_nb"+nb_str, "HDP events, Nb="+nb_str, no_bin_bjet_w_exclusion_w_mhtc[nb], 0.5, no_bin_bjet_w_exclusion_w_mhtc[nb] + 0.5 ) ;
      set_bin_labels_nb( h_hdp_nb[nb], nb) ;

      h_ldp_nb[nb] = new TH1F( "h_ldp_nb"+nb_str, "ldp events, Nb="+nb_str, no_bin_bjet_w_exclusion_w_mhtc[nb], 0.5, no_bin_bjet_w_exclusion_w_mhtc[nb] + 0.5 ) ;
      set_bin_labels_nb( h_ldp_nb[nb], nb ) ;

   }//nb


   TH1F* h_hdp = new TH1F( "h_hdp", "HDP events", nb_global_w_exclusion_w_mhtc, 0.5, nb_global_w_exclusion_w_mhtc + 0.5 ) ;
   set_bin_labels_div_by_nb( h_hdp ) ;
   TH1F* h_nbsum_hdp = new TH1F( "h_nbsum_hdp", "HDP events", max_no_bin_bjet_w_exclusion_w_mhtc, 0.5, max_no_bin_bjet_w_exclusion_w_mhtc + 0.5 ) ;
   set_bin_labels( h_nbsum_hdp ) ;

   int no_bin_mhtc = nb_global_w_exclusion_w_mhtc - nb_global_after_exclusion ;  

   TH1F* h_mhtc_hdp = new TH1F( "h_mhtc_hdp", "HDP events, MHTC", no_bin_mhtc, 0.5, no_bin_mhtc + 0.5 ) ;
   set_bin_labels_mhtc_plot( h_mhtc_hdp ) ;

   TH1F* h_mhtc_ldp = new TH1F( "h_mhtc_ldp", "ldp events, MHTC", no_bin_mhtc, 0.5, no_bin_mhtc+0.5 ) ;
   set_bin_labels_mhtc_plot( h_mhtc_ldp ) ;

   int no_bin_mhtc_nbsum = max_no_bin_bjet_w_exclusion_w_mhtc - max_no_bin_bjet;


   TH1F* h_ldp = new TH1F( "h_ldp", "ldp events", nb_global_w_exclusion_w_mhtc, 0.5, nb_global_w_exclusion_w_mhtc + 0.5 ) ;
   set_bin_labels_div_by_nb( h_ldp ) ;
   TH1F* h_nbsum_ldp = new TH1F( "h_nbsum_ldp", "LDP events", max_no_bin_bjet_w_exclusion_w_mhtc, 0.5, max_no_bin_bjet_w_exclusion_w_mhtc + 0.5 ) ;
   set_bin_labels( h_nbsum_ldp ) ;


   TH1F* h_mht_all = new TH1F( "h_mht_all", "MHT, all events", 100, 0., 1000. ) ;
   TH1F* h_mht_badmu = new TH1F( "h_mht_badmu", "MHT, badmu", 100, 0., 1000. ) ;
   TH1F* h_mht_allrejected = new TH1F( "h_mht_allrejected", "MHT, all rejected", 100, 0., 1000. ) ;

   TH1F* h_met_over_calomet_all = new TH1F( "h_met_over_calomet_all", "MET/CaloMET, all events", 100., 0., 10. ) ;
   TH1F* h_met_over_calomet_badmu = new TH1F( "h_met_over_calomet_badmu", "MET/CaloMET, badmu", 100., 0., 10. ) ;
   TH1F* h_met_over_calomet_allrejected = new TH1F( "h_met_over_calomet_allrejected", "MET/CaloMET, all rejected", 100., 0., 10. ) ;

   TH1F* h_mdp_all = new TH1F( "h_mdp_all", "Min Delta phi, 4 leading jets", 64, 0., 3.2 ) ;

   TH1F* h_mht = new TH1F( "h_mht", "MHT", 80, 150., 350. ) ;
   TH2F* h_mht_vs_met = new TH2F( "h_mht_vs_met", "MHT vs MET", 100, 0., 400., 100, 0., 400. ) ;
   TH2F* h_mht_vs_calomet = new TH2F( "h_mht_vs_calomet", "MHT vs Calo MET", 100, 0., 400., 100, 0., 400. ) ;
   TH2F* h_met_vs_calomet = new TH2F( "h_met_vs_calomet", "MET vs Calo MET", 100, 0., 400., 100, 0., 400. ) ;
   TH1F* h_mht_met_gt_mht_minus_100 = new TH1F( "h_mht_met_gt_mht_minus_100", "MHT, MET>(MHT-100)", 80, 150., 350. ) ;
   TH1F* h_mht_met_gt_mht_minus_100_calomet80 = new TH1F( "h_mht_met_gt_mht_minus_100_calomet80", "MHT, MET>(MHT-100), CaloMET>80", 80, 150., 350. ) ;
   TH1F* h_mht_calomet80_pfmet100 = new TH1F( "h_mht_calomet80_pfmet100", "MHT, CaloMET>80, pfmet>100", 80, 150., 350. ) ;

   TH1F* h_calomet = new TH1F( "h_calomet", "Calomet", 100, 0., 400. ) ;
   TH1F* h_met = new TH1F( "h_met", "MET", 100, 0., 400. ) ;

   Long64_t nentries = fChain->GetEntries();

   printf("\n\n") ;
   printf("  Looping over sample: %s\n", samplename ) ;
   printf("  Number of entries: %lld\n\n", nentries ) ;

  //FILE* ofp = fopen( "debug.txt", "w" ) ;

   Long64_t loopmax = nentries ;
   if ( nloop > 0 ) loopmax = nloop ;

   TStopwatch sw ;
   sw.Start() ;
   int time(0) ;
   float projected_remaining(999999.) ;
   for (Long64_t jentry=0; jentry<loopmax;jentry++) 
   {

      Long64_t ievt = jentry ;

      if ( ievt%1000 == 0 ) { // timer printing stuff
         int thistime = sw.RealTime() ;
         sw.Continue() ;
         if ( thistime < 2 ) {
            printf("   %10llu out of %10llu  (%6.1f%%) \r", ievt, nentries, 100.*ievt/(1.*nentries) ) ;
         } else {
            if ( thistime > time ) projected_remaining = (1.*thistime)/(1.*ievt)*(nentries-ievt) ;
            if ( projected_remaining < 100 ) {
               printf("   %10llu out of %10llu  (%6.1f%%)    seconds remaining %4.0f                       \r", ievt, nentries, 100.*ievt/(1.*nentries), projected_remaining ) ;
            } else if ( projected_remaining < 3600 ) {
               printf("   %10llu out of %10llu  (%6.1f%%)    time remaining     %2d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries),
                    TMath::Nint(projected_remaining)/60, TMath::Nint(projected_remaining)%60 ) ;
            } else {
               printf("   %10llu out of %10llu  (%6.1f%%)    time remaining  %2d:%02d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries),
                    TMath::Nint(projected_remaining)/3600, (TMath::Nint(projected_remaining)%3600)/60, TMath::Nint(projected_remaining)%60 ) ;
            }
         }
         fflush(stdout) ;
         time = thistime ;
      } // timer printing stuff

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      fChain->GetEntry(jentry);

      bi_htmht         = find_htmht_bin_no(HT, MHT);
      if ( bi_htmht < 0 ) continue ;


      double hw = 1.;

     //-- take out the trash
      bool badMuonEvent(false) ;
      bool isjunk(false) ;
      for ( unsigned long ji=0; ji<Jets->size(); ji++ ) {
         if ( Jets->at(ji).Pt() < 200 ) continue ;
         if ( Jets_muonEnergyFraction->at(ji) < 0.5 ) continue ;
         double dPhi = Jets->at(ji).Phi() - METPhi ;
         if ( dPhi >  3.1415926 ) dPhi = dPhi - 2*3.14159 ;
         if ( dPhi < -3.1415926 ) dPhi = dPhi + 2*3.14159 ;
         if ( fabs( dPhi ) > 3.1415926 - 0.40 ) { 
            badMuonEvent = true ;
            break ;
         }
      } // ji
      if ( badMuonEvent ) {
       //printf("\n\n *** Bad Muon event: %lld\n", jentry ) ;
       //printf("    MET = %7.1f  ,  MHT = %7.1f ,  CaloMET = %7.1f ,   METPhi = %6.3f,   MHTPhi = %6.3f\n", MET, MHT, CaloMET, METPhi, MHTPhi ) ;
       //for ( unsigned long ji=0; ji<Jets->size(); ji++ ) {
       //   if ( Jets->at(ji).Pt() < 30 ) break ;
       //   printf(" jet %2lu :  Pt = %7.1f,  phi = %6.3f,  eta = %6.3f,  Muon fr = %5.3f\n",
       //        ji, Jets->at(ji).Pt(), Jets->at(ji).Phi(), Jets->at(ji).Eta(), Jets_muonEnergyFraction->at(ji) ) ;
       //   fflush(stdout) ;
       //} // ji
      } // badMuonEvent?


      if ( CaloMET <= 0 )  isjunk = true ;
      if ( eeBadScFilter != 1 )  isjunk = true ;
      ////////////if ( CSCTightHaloFilter != 1 )  isjunk = true ; ///// Don't apply this!  July 19.
      if ( HBHEIsoNoiseFilter != 1 )  isjunk = true ;
      if ( HBHENoiseFilter != 1 )  isjunk = true ;
      if ( EcalDeadCellTriggerPrimitiveFilter != 1 )  isjunk = true ; ///// new for v9
      if ( NVtx <= 0 )  isjunk = true ;  // **** new for v9
      if ( ! JetID ) isjunk = true ;  // **** new for v9
      ///////if ( badMuonEvent ) isjunk = true ;
      ///////if ( MET/CaloMET > 5 )  isjunk = true ;
      if ( ! noMuonJet ) isjunk = true ;  // **** new for v9, should be redundant with above
      if ( PFCaloMETRatio >= 5 ) isjunk = true ;  // **** new for v9, should be redundant with above
      if ( globalTightHalo2016Filter!=1 ) isjunk = true ;  // **** new for v9
      if ( ! BadChargedCandidateFilter ) isjunk = true ;  // **** new for v9
      if ( ! BadPFMuonFilter ) isjunk = true ;  // **** new for v9

      /////////if ( !(pass_pfmet100_trig||pass_pfmetnomu100_trig) ) isjunk = true ;
      if ( !(pass_pfmet100_trig||pass_pfmetnomu100_trig || pass_pfmet110_trig||pass_pfmetnomu110_trig || pass_pfmet120_trig||pass_pfmetnomu120_trig) ) isjunk = true ;


  /// if ( RunNum==274344 && LumiBlockNum==615 && EvtNum==953880058 ) {
  ///    printf("\n\n ******* The missing event.\n\n") ;
  ///    printf("   Is junk is %s\n", isjunk?"T":"F") ;
  ///    printf("   CaloMET: %f\n", CaloMET) ;
  ///    printf("   eeBadScFilter  %d\n", eeBadScFilter ) ;
  ///    printf("   CSCTightHaloFilter  %d\n", CSCTightHaloFilter ) ;
  ///    printf("   HBHEIsoNoiseFilter  %d\n", HBHEIsoNoiseFilter ) ;
  ///    printf("\n\n") ;
  /// }



      h_mht_all -> Fill( MHT, hw ) ;
      if ( CaloMET > 0 ) h_met_over_calomet_all -> Fill( MET/CaloMET, hw ) ;

      if ( badMuonEvent ) {
         h_mht_badmu -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) h_met_over_calomet_badmu -> Fill( MET/CaloMET, hw ) ;
      }


      if ( isjunk ) {
         h_mht_allrejected -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) h_met_over_calomet_allrejected -> Fill( MET/CaloMET, hw ) ;
      }

      if ( (RunNum==275936 && LumiBlockNum==430) ) isjunk = true ; //*** bad lumi section.  July 19th.

      if ( isjunk ) continue ;

    //-- new baseline

      bi_nj            = find_njet_bin_no (NJets);
      if ( bi_nj < 0    ) continue ;
      bi_nb            = find_nbjet_bin_no(BTags);
      if ( bi_nb < 0    ) continue ;


      h_mht -> Fill( MHT, hw ) ;
      h_mht_vs_met -> Fill( MET, MHT, hw ) ;
      h_mht_vs_calomet -> Fill( CaloMET, MHT, hw ) ;
      h_met_vs_calomet -> Fill( CaloMET, MET, hw ) ;
      if ( MET>(MHT-100) ) h_mht_met_gt_mht_minus_100 -> Fill( MHT, hw ) ;
      if ( MET>(MHT-100) && CaloMET>80 ) h_mht_met_gt_mht_minus_100_calomet80 -> Fill( MHT, hw ) ;
      if ( CaloMET>80 && MET>100 ) h_mht_calomet80_pfmet100 -> Fill( MHT, hw ) ;
      h_met -> Fill( MET, hw ) ;
      h_calomet -> Fill( CaloMET, hw ) ;


      bool in_ldp(true) ;
      if ( DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3 ) in_ldp = false ;
      char tag_ldp[10] ;
      if ( in_ldp ) { sprintf( tag_ldp, "LDP" ) ; } else { sprintf( tag_ldp, "HDP" ) ; }

      double minDeltaPhi = 99. ;
      if ( DeltaPhi1 < minDeltaPhi ) minDeltaPhi = DeltaPhi1 ;
      if ( DeltaPhi2 < minDeltaPhi ) minDeltaPhi = DeltaPhi2 ;
      if ( DeltaPhi3 < minDeltaPhi ) minDeltaPhi = DeltaPhi3 ;
      if ( DeltaPhi4 < minDeltaPhi ) minDeltaPhi = DeltaPhi4 ;

      h_mdp_all -> Fill( minDeltaPhi, hw ) ;

      translate_htmht_bin_to_ht_and_mht_bins (bi_htmht, bi_ht, bi_mht);
      bi_global        = global_bin_with_mhtc ( bi_nj, bi_nb, bi_htmht);
      bi_nbsum_global  = global_bin_with_mhtc_sum_nb ( bi_nj, bi_htmht);
      bi_nbsum_global_nb  = global_bin_with_mhtc_nb ( bi_nj, bi_nb, bi_htmht);

      if ( bi_global < 0 ) continue; // remove excluded bins



      int bi_mhtc_plot       = global_bin_only_mhtc ( bi_nj, bi_nb, bi_htmht ) ;
 
      if ( !in_ldp ) {
         if ( !(blind && bi_mht>1) ) {
            h_hdp -> Fill( bi_global, hw ) ;
            //if ( bi_global == 30 ) fprintf( ofp, "%20u %20u %20lld\n", RunNum, LumiBlockNum, EvtNum ) ;	    
            h_nbsum_hdp -> Fill( bi_nbsum_global, hw ) ;

            for ( int nb = 0; nb < nb_nb; nb++)
               if ( BTags > bin_edges_nb[nb] && BTags < bin_edges_nb[nb+1] ) h_hdp_nb[nb] -> Fill( bi_nbsum_global_nb, hw ) ;
	 }//(blind && bi_mht>1)?

         if ( bi_mht == 1 ) h_mhtc_hdp -> Fill( bi_mhtc_plot , hw ) ;
      } else {
            h_ldp -> Fill( bi_global, hw ) ;
            h_nbsum_ldp -> Fill( bi_nbsum_global, hw ) ;

         for ( int nb = 0; nb < nb_nb; nb++)
            if ( BTags > bin_edges_nb[nb] && BTags < bin_edges_nb[nb+1] ) h_ldp_nb[nb] -> Fill( bi_nbsum_global_nb, hw ) ;

         if ( bi_mht == 1 ) h_mhtc_ldp -> Fill( bi_mhtc_plot , hw ) ;
      }


      if ( verb ) {
         printf("\n\n ===== RunNum %6d, LumiBlockNum %5d, EvtNum %9llu\n", RunNum, LumiBlockNum, EvtNum ) ;	      
         printf("  Number of bins:\n") ;
         printf("     QCD:     %d HT,  %d MHT,  %2d HTMHT,   %d Njet,  %d Nb,   %3d total\n", nb_ht[1], nb_mht, nb_htmht, nb_nj, nb_nb, nb_global ) ;
         printf("  Essential event vars:\n" ) ;
         printf("    NJets = %2d,  bin %d\n", NJets, bi_nj ) ;
         printf("    Nb    = %2d,  bin %d\n", BTags, bi_nb ) ;
         printf("    MHT   = %7.1f,  bin %2d\n", MHT, bi_mht ) ;
         printf("    HT    = %7.1f,  bin %2d\n", HT, bi_ht ) ;
         printf("      htmht bin %2d\n", bi_htmht ) ;
         printf("      global bin %3d\n", bi_global ) ;
         printf("    Dphi1 = %6.3f,  Dphi1 = %6.3f, Dphi1 = %6.3f, Dphi1 = %6.3f,  %s\n", DeltaPhi1, DeltaPhi2, DeltaPhi3, DeltaPhi4, tag_ldp ) ;
      } // verb?



   } // jentry

   //fclose( ofp ) ;


   printf("\n\n\n Done.\n\n\n") ;

   char histfile[10000] ;
   if ( strlen( samplename ) > 0 ) {
sprintf( histfile, "outputfiles/hists-data-v2d-%s.root", samplename ) ;
   } else {
      sprintf( histfile, "outputfiles/hists-data-v2d.root" ) ;
   }
   saveHist( histfile, "h*" ) ;

} // Loop

//=====================================================================================================

   
void fill_data_hists_loop_v2d::set_bin_labels_nb( TH1F* hp, int tbi_nb) {
   //tbi_nb starts from zero
   if ( hp == 0x0 ) return ;

   int control_bi = 0, print_bi = 0;
   int bi = 0;
   for ( int tbi_nj = 1; tbi_nj <= nb_nj;    tbi_nj++){
      for ( int tbi_htmht = 1; tbi_htmht <= nb_htmht; tbi_htmht++) {
      int tbi_ht, tbi_mht ;
      translate_htmht_bin_to_ht_and_mht_bins (tbi_htmht, tbi_ht, tbi_mht ) ;

      if ( is_this_bin_excluded (tbi_nj-1, tbi_nb, tbi_ht-1, tbi_mht-1 ) ) continue;
      bi++;
      char binlabel[100] ;
      //---------
      // sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
      //---------
      int print_bi_htmht = tbi_htmht - nb_ht[1] ;
      int print_bi_mht = tbi_mht - 1 ;
      if ( tbi_mht == 1 ) {
         control_bi++;
         sprintf( binlabel, "Nj%d-MHTC-HT%d (C%1d)   C%2d", tbi_nj, tbi_ht, tbi_ht, control_bi ) ;
      } else {
            print_bi++;
            sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, print_bi_mht, tbi_ht, print_bi_htmht, print_bi ) ;
      }
      //---------
      hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bin_htmht
   } // bin_nj
   hp -> GetXaxis() -> LabelsOption( "v" ) ;

}//set_bin_labels_nb

//=====================================================================================================


void fill_data_hists_loop_v2d::set_bin_labels( TH1F* hp ) {

   if ( hp == 0x0 ) return ;

   int control_bi = 0, print_bi = 0;
   int bi = 0;
   for ( int tbi_nj = 1; tbi_nj <= nb_nj;    tbi_nj++){
      for ( int tbi_htmht = 1; tbi_htmht <= nb_htmht; tbi_htmht++) {
      int tbi_ht, tbi_mht ;
      translate_htmht_bin_to_ht_and_mht_bins (tbi_htmht, tbi_ht, tbi_mht ) ;

      if ( is_this_bin_excluded_nbsum (tbi_nj-1, tbi_ht-1, tbi_mht-1 ) ) continue;
      bi++;
      char binlabel[100] ;
      //---------
      // sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
      //---------
      int print_bi_htmht = tbi_htmht - nb_ht[1] ;
      int print_bi_mht = tbi_mht - 1 ;
      if ( tbi_mht == 1 ) {
         control_bi++;
         sprintf( binlabel, "Nj%d-MHTC-HT%d (C%1d)   C%2d", tbi_nj, tbi_ht, tbi_ht, control_bi ) ;
      } else {
            print_bi++;
            sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, print_bi_mht, tbi_ht, print_bi_htmht, print_bi ) ;
      }
      //---------
      hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bin_htmht
   } // bin_nj
   hp -> GetXaxis() -> LabelsOption( "v" ) ;

} // set_bin_labels

//===================================

   void fill_data_hists_loop_v2d::set_bin_labels_div_by_nb( TH1F* hp ) {

      if ( hp == 0x0 ) return ;
      int control_bi = 0, print_bi = 0;
      for ( int bi=1; bi<=nb_global_w_exclusion_w_mhtc; bi++ ) {
         int tbi_nj; int tbi_nb; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_qcd_bin_to_nj_nb_ht_mht ( bi, tbi_nj,  tbi_nb, tbi_ht, tbi_mht);
         translate_ht_and_mht_bin_to_htmht_bins( tbi_ht, tbi_mht, tbi_htmht);
	 if ( is_this_bin_excluded(tbi_nj-1,  tbi_nb-1, tbi_ht-1, tbi_mht-1) ) continue;

         char binlabel[100] ;
         //---------
         // sprintf( binlabel, "NJets%d-MHT%d-HT%d (%2d)   %2d", tbi_nj-1, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
         //---------
         int print_bi_htmht = tbi_htmht - nb_ht[1] ;
         int print_bi_mht = tbi_mht - 1 ;
         int print_bi_nb = tbi_nb -1 ;
         if ( tbi_mht == 1 ) {
            control_bi++;
            sprintf( binlabel, "Nj%d-Nb%d-MHTC-HT%d (C%1d)   C%2d  %3d", tbi_nj, print_bi_nb, tbi_ht, tbi_ht, control_bi, bi ) ;
         } else {
               print_bi++;
               sprintf( binlabel, "Nj%d-Nb%d-MHT%d-HT%d (%2d)   S%3d  %3d", tbi_nj, print_bi_nb, print_bi_mht, tbi_ht, print_bi_htmht, print_bi, bi ) ;
         }
         //---------
         hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bi

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_div_by_nb

   //===================================

   void fill_data_hists_loop_v2d::set_bin_labels_mhtc_plot( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      int bi_plot(0) ;
      for ( int tbi_nj=1; tbi_nj<=nb_nj; tbi_nj++ ) {
         for ( int tbi_nb=1; tbi_nb<=nb_nb; tbi_nb++ ) {
            for ( int tbi_ht=1; tbi_ht<=(nb_ht[1]); tbi_ht++ ) {
               if ( is_this_bin_excluded(tbi_nj-1, tbi_nb-1, tbi_ht-1, 0) ) continue;
               bi_plot ++ ;
               char binlabel[100] ;
               sprintf( binlabel, "Nj%d-Nb%d-MHTC-HT%d  C%2d", tbi_nj, tbi_nb-1, tbi_ht, bi_plot ) ;
               hp -> GetXaxis() -> SetBinLabel( bi_plot, binlabel ) ;
            } // tbi_ht
         } // tbi_nb
      } // tbi_nj

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_mhtc_plot

   //===================================

#endif

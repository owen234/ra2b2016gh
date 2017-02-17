#ifndef fill_hists_loop_v2d_c
#define fill_hists_loop_v2d_c

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

#define fill_hists_loop_v2d_cxx
#include "fill_hists_loop_v2d.h"
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


int bi_nj, bi_nb, bi_ht, bi_mht, bi_htmht;

void fill_hists_loop_v2d::Loop( bool verb, int nloop )
{

   setup_bins() ;



   bool islostlep(false) ;
   bool ishadtau(false) ;
   bool isqcd(false) ;
   TString ts_samplename( samplename ) ;

   if ( ts_samplename.Contains("lostlep") ) {
      printf("\n\n This is a lost lepton run through.\n\n") ;
      islostlep = true ;
   }
   if ( ts_samplename.Contains("hadtau") ) {
      printf("\n\n This is a had tau run through.\n\n") ;
      ishadtau = true ;
   }
   if ( ts_samplename.Contains("qcd") ) {
      printf("\n\n This is a qcd run through.\n\n") ;
      isqcd = true ;
   }

   gDirectory -> Delete( "h*" ) ;

   if ( isqcd ) loadHist( "outputfiles/data-turnon.root" ) ;

   TH1F* h_turnon[10] ;
   if ( isqcd ) {

      for ( int nj = 1; nj <= nb_nj; nj++)
      {

         TString nj_str         ; nj_str         .Form("%d",nj);

         h_turnon[nj] = (TH1F*) gDirectory -> FindObject( "h_eff_nj"+nj_str ) ;
         if ( h_turnon[nj] == 0x0 ) { printf("\n\n *** missing nj%d turnon.\n\n", nj) ; return ; }

      }

   }



   if (fChain == 0) return;

   TH1F* h_hdp_nb[10], *h_ldp_nb[10];

   TH1::SetDefaultSumw2(); // to make sure all TH1 histograms have Sumw2 enabled

   for ( int nb = 0; nb < nb_nb; nb++)
   {

      TString nb_str         ; nb_str         .Form("%d",nb);
      h_hdp_nb[nb] = new TH1F( "h_hdp_nb"+nb_str, "HDP events, Nb="+nb_str, nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ;
      set_bin_labels( h_hdp_nb[nb] ) ;

      h_ldp_nb[nb] = new TH1F( "h_ldp_nb"+nb_str, "ldp events, Nb="+nb_str, nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ;
      set_bin_labels( h_ldp_nb[nb] ) ;


   }


   TH1F* h_hdp = new TH1F( "h_hdp", "HDP events", nb_global, 0.5, nb_global + 0.5 ) ;
   set_bin_labels_div_by_nb( h_hdp ) ;
   TH1F* h_nbsum_hdp = new TH1F( "h_nbsum_hdp", "HDP events", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ;
   set_bin_labels( h_nbsum_hdp ) ;

   TH1F* h_mhtc_hdp = new TH1F( "h_mhtc_hdp", "HDP events, MHTC", nb_nj*nb_nb*nb_ht[1], 0.5, nb_nj*nb_nb*nb_ht[1]+0.5 ) ;
   set_bin_labels_mhtc_plot( h_mhtc_hdp ) ;
   TH1F* h_mhtc_nbsum_hdp = new TH1F( "h_mhtc_nbsum_hdp", "HDP events, MHTC, Nbsum", nb_nj*nb_ht[1], 0.5, nb_nj*nb_ht[1]+0.5 ) ;
   set_bin_labels_mhtc_nbsum_plot( h_mhtc_nbsum_hdp ) ;

   TH1F* h_ldp = new TH1F( "h_ldp", "ldp events", nb_global, 0.5, nb_global + 0.5 ) ;
   set_bin_labels_div_by_nb( h_ldp ) ;
   TH1F* h_nbsum_ldp = new TH1F( "h_nbsum_ldp", "LDP events", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ;
   set_bin_labels( h_nbsum_ldp ) ;
   TH1F* h_mhtc_ldp = new TH1F( "h_mhtc_ldp", "ldp events, MHTC", nb_nj*nb_nb*nb_ht[1], 0.5, nb_nj*nb_nb*nb_ht[1]+0.5 ) ;
   set_bin_labels_mhtc_plot( h_mhtc_ldp ) ;
   TH1F* h_mhtc_nbsum_ldp = new TH1F( "h_mhtc_nbsum_ldp", "ldp events, MHTC, Nbsum", nb_nj*nb_ht[1], 0.5, nb_nj*nb_ht[1]+0.5 ) ;
   set_bin_labels_mhtc_nbsum_plot( h_mhtc_nbsum_ldp ) ;

   TH1F* h_jet_eta_badmu = new TH1F( "h_jet_eta_badmu", "Jet eta, bad muon", 55, -5.5, 5.5 ) ;

   TH1F* h_mht_all = new TH1F( "h_mht_all", "MHT, all events", 100, 0., 1000. ) ;
   TH1F* h_mht_badmu = new TH1F( "h_mht_badmu", "MHT, badmu", 100, 0., 1000. ) ;
   TH1F* h_mht_met100_calomet80_rejected = new TH1F( "h_mht_met100_calomet80_rejected", "MHT, MET<100 or CaloMET<80", 100, 0., 1000. ) ;
   TH1F* h_mht_filters = new TH1F( "h_mht_filters", "MHT, all rejected", 100, 0., 1000. ) ;
   TH1F* h_mht_allrejected = new TH1F( "h_mht_allrejected", "MHT, all rejected", 100, 0., 1000. ) ;

   TH1F* h_met_over_calomet_all = new TH1F( "h_met_over_calomet_all", "MET/CaloMET, all events", 100., 0., 10. ) ;
   TH1F* h_met_over_calomet_badmu = new TH1F( "h_met_over_calomet_badmu", "MET/CaloMET, badmu", 100., 0., 10. ) ;
   TH1F* h_met_over_calomet_met100_calomet80_rejected = new TH1F( "h_met_over_calomet_met100_calomet80_rejected", "MET/CaloMET, MET<100 or CaloMET<80", 100., 0., 10. ) ;
   TH1F* h_met_over_calomet_filters = new TH1F( "h_met_over_calomet_filters", "MET/CaloMET, all rejected", 100., 0., 10. ) ;
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
   TH1F* h_calomet_met_gt100 = new TH1F( "h_calomet_met_gt100", "Calomet, MET>100", 100, 0., 400. ) ;
   TH1F* h_met_calomet_gt80 = new TH1F( "h_met_calomet_gt80", "MET, Calomet>80", 100, 0., 400. ) ;

   TH1F* h_ht_after_cleaning = new TH1F( "h_ht_after_cleaning", "HT, after all QCD junk cleaning", 50, 0., 5000. ) ;

   TH1F* h_ldp_weight[1000] ;
   for ( int sbi=1; sbi<=nb_global; sbi++ ) {
      char hname[100] ;
      sprintf( hname, "h_ldp_weight_bin%03d", sbi ) ;
      h_ldp_weight[sbi] = new TH1F( hname, hname, 100, 0., 2. ) ;
   }

   TH1F* h_ht_ldp_mhtc = new TH1F( "h_ht_ldp_mhtc", "HT, LDP, MHTC", 50, 0., 5000. ) ;
   TH1F* h_ht_ldp_mht1 = new TH1F( "h_ht_ldp_mht1", "HT, LDP, MHT1", 50, 0., 5000. ) ;
   TH1F* h_ht_ldp_mht2 = new TH1F( "h_ht_ldp_mht2", "HT, LDP, MHT2", 50, 0., 5000. ) ;
   TH1F* h_ht_ldp_mht3 = new TH1F( "h_ht_ldp_mht3", "HT, LDP, MHT3", 50, 0., 5000. ) ;
   TH1F* h_ht_ldp_mht4 = new TH1F( "h_ht_ldp_mht4", "HT, LDP, MHT4", 50, 0., 5000. ) ;

   Long64_t nentries = fChain->GetEntries();

   printf("\n\n") ;
   printf("  Looping over sample: %s\n", samplename ) ;
   printf("  Number of entries: %lld\n\n", nentries ) ;


   Long64_t loopmax = nentries ;
   if ( nloop > 0 ) loopmax = nloop ;

   Long64_t nbytes = 0, nb = 0;

   TStopwatch sw ;
   sw.Start() ;
   int time(0) ;
   float projected_remaining(999999.) ;

   for (Long64_t jentry=0; jentry<loopmax;jentry++) {

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
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      set_bi() ;


      //////////////double hw = Weight * lumi_ ;
      double hw = Weight * puWeight * lumi_ ; // add pu weight

      if ( islostlep && hasHadTau ) continue ;
      if ( ishadtau && !hasHadTau ) continue ;

     //-- take out the trash
      bool badMuonEvent(false) ;
      bool isjunk(false) ;
      double badmu_jet_eta(-9.) ;
      for ( unsigned long ji=0; ji<Jets->size(); ji++ ) {
         if ( Jets->at(ji).Pt() < 200 ) continue ;
         if ( Jets_muonEnergyFraction->at(ji) < 0.5 ) continue ;
         double dPhi = Jets->at(ji).Phi() - METPhi ;
         if ( dPhi >  3.1415926 ) dPhi = dPhi - 2*3.14159 ;
         if ( dPhi < -3.1415926 ) dPhi = dPhi + 2*3.14159 ;
         if ( fabs( dPhi ) > 3.1415926 - 0.40 ) { 
            badMuonEvent = true ;
            badmu_jet_eta = Jets->at(ji).Eta() ;
            break ;
         }
      } // ji
      if ( badMuonEvent ) { h_jet_eta_badmu -> Fill( badmu_jet_eta, hw ) ; }
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



      h_mht_all -> Fill( MHT, hw ) ;
      if ( CaloMET > 0 ) h_met_over_calomet_all -> Fill( MET/CaloMET, hw ) ;

      if ( CaloMET <= 0 )  isjunk = true ;
      if ( eeBadScFilter < 1 )  isjunk = true ;
      ////////////if ( CSCTightHaloFilter < 1 )  isjunk = true ; //*** Don't apply this!  July 19.
      if ( HBHEIsoNoiseFilter < 1 )  isjunk = true ;
      if ( HBHENoiseFilter < 1 )  isjunk = true ;
      if ( EcalDeadCellTriggerPrimitiveFilter < 1 )  isjunk = true ; // **** new for v9
      if ( NVtx <= 0 )  isjunk = true ;  // **** new for v9
      if ( ! JetID ) isjunk = true ;  // **** new for v9


      if ( isjunk ) {
         h_mht_filters -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) h_met_over_calomet_filters -> Fill( MET/CaloMET, hw ) ;
      }

      if ( MET<100 || CaloMET<80 ) {
         h_mht_met100_calomet80_rejected -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) {
            h_met_over_calomet_met100_calomet80_rejected -> Fill( MET/CaloMET, hw ) ;
         }
      }

      if ( badMuonEvent ) {
         h_mht_badmu -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) h_met_over_calomet_badmu -> Fill( MET/CaloMET, hw ) ;
      }


      if ( badMuonEvent ) isjunk = true ; // --- Feb 7, 2017 : add this back in due to bug(?) in noMuonJet
      ///////if ( MET/CaloMET > 5 )  isjunk = true ;
      if ( ! noMuonJet ) isjunk = true ;  // **** new for v9, should be redundant with above
      if ( PFCaloMETRatio > 5 ) isjunk = true ;  // **** new for v9, should be redundant with above

      if ( isjunk ) {
         h_mht_allrejected -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) h_met_over_calomet_allrejected -> Fill( MET/CaloMET, hw ) ;
      }


      if ( isjunk ) continue ;

    //-- new baseline
      if ( NJets < bin_edges_nj[0] ) continue ;
      if ( HT < bin_edges_ht[1][0] ) continue ;
      ////////if ( bi_mht <= 0 ) continue ;



     //--- if this is QCD, correct for the trigger turnon
      if ( isqcd ) {

         int turnon_bin(0) ;
         float trig_eff(1.) ;
     /// if ( bi_nj>=1 && bi_nj<= nb_nj ) {
     ///    if ( MHT < (h_turnon[bi_nj] -> GetXaxis() -> GetXmax()) ) {
     ///       turnon_bin = h_turnon[bi_nj] -> GetXaxis() -> FindBin( MHT ) ;
     ///    } else {
     ///       turnon_bin = h_turnon[bi_nj] -> GetNbinsX() ;
     ///    }
     /// }
     /// trig_eff = h_turnon[bi_nj] -> GetBinContent( turnon_bin ) ;
         hw = hw * trig_eff ;

         /////////////////if ( ! (pass_pfmet100_trig || pass_pfmetnomu100_trig) ) continue ;

      }




      h_mht -> Fill( MHT, hw ) ;
      h_mht_vs_met -> Fill( MET, MHT, hw ) ;
      h_mht_vs_calomet -> Fill( CaloMET, MHT, hw ) ;
      h_met_vs_calomet -> Fill( CaloMET, MET, hw ) ;
      if ( MET>(MHT-100) ) h_mht_met_gt_mht_minus_100 -> Fill( MHT, hw ) ;
      if ( MET>(MHT-100) && CaloMET>80 ) h_mht_met_gt_mht_minus_100_calomet80 -> Fill( MHT, hw ) ;
      if ( CaloMET>80 && MET>100 ) h_mht_calomet80_pfmet100 -> Fill( MHT, hw ) ;
      h_met -> Fill( MET, hw ) ;
      h_calomet -> Fill( CaloMET, hw ) ;
      if ( CaloMET>80 ) h_met_calomet_gt80 -> Fill( MET, hw ) ;
      if ( MET>100 ) h_calomet_met_gt100 -> Fill( CaloMET, hw ) ;

     //******** adding these new cleanup cuts
      if ( MET<100 ) continue ;
      if ( CaloMET<80 ) continue ;
     //********


      h_ht_after_cleaning -> Fill( HT, hw ) ;


      ////////********//////////if ( MHT < 250 ) continue ;

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

      int bi_mhtc_plot = (bi_nj-1)*nb_nb*(nb_ht[1]) + (bi_nb-1)*(nb_ht[1]) + bi_ht ;
      int bi_mhtc_nbsum_plot = (bi_nj-1)*(nb_ht[1]) + bi_ht ;
      if ( !in_ldp ) {
         h_hdp -> Fill( bi_global, hw ) ;
         h_nbsum_hdp -> Fill( bi_nbsum_global, hw ) ;

         for ( int nb = 0; nb < nb_nb; nb++)
            if ( BTags > bin_edges_nb[nb] && BTags < bin_edges_nb[nb+1] ) h_hdp_nb[nb] -> Fill( bi_nbsum_global, hw ) ;

         if ( bi_mht == 1 ) h_mhtc_hdp -> Fill( bi_mhtc_plot , hw ) ;
         if ( bi_mht == 1 ) h_mhtc_nbsum_hdp -> Fill( bi_mhtc_nbsum_plot , hw ) ;
      } else {
         h_ldp -> Fill( bi_global, hw ) ;
         if ( bi_global > 0 && bi_global < nb_global ) h_ldp_weight[bi_global] -> Fill( hw ) ;
         h_nbsum_ldp -> Fill( bi_nbsum_global, hw ) ;

         for ( int nb = 0; nb < nb_nb; nb++)
            if ( BTags > bin_edges_nb[nb] && BTags < bin_edges_nb[nb+1] ) h_ldp_nb[nb] -> Fill( bi_nbsum_global, hw ) ;

         if ( bi_mht == 1 ) h_mhtc_ldp -> Fill( bi_mhtc_plot , hw ) ;
         if ( bi_mht == 1 ) h_mhtc_nbsum_ldp -> Fill( bi_mhtc_nbsum_plot , hw ) ;
         if ( bi_mht == 1 && HT>300 ) h_ht_ldp_mhtc -> Fill( HT, hw ) ;
         if ( bi_mht == 2 && HT>300 ) h_ht_ldp_mht1 -> Fill( HT, hw ) ;
         if ( bi_mht == 3 && HT>350 ) h_ht_ldp_mht2 -> Fill( HT, hw ) ;
         if ( bi_mht == 4 && HT>500 ) h_ht_ldp_mht3 -> Fill( HT, hw ) ;
         if ( bi_mht == 5 && HT>750 ) h_ht_ldp_mht4 -> Fill( HT, hw ) ;
      }


      if ( verb ) {
         printf("\n\n ===== RunNum %6d, LumiBlockNum %5d, EvtNum %9llu,  Weight %g\n", RunNum, LumiBlockNum, EvtNum, Weight ) ;
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


   TH1F* h_max_ldp_weight = new TH1F( "h_max_ldp_weight", "Max ldp weight", nb_global, 0.5, nb_global + 0.5 ) ;

   for ( int sbi=1; sbi<=nb_global; sbi++ ) {
      float max_weight(0.) ;
      for ( int hbi=1; hbi<=h_ldp_weight[sbi]->GetNbinsX(); hbi++ ) {
         float c = h_ldp_weight[sbi] -> GetBinContent( hbi ) ;
         if ( c > 0 ) {
            max_weight = h_ldp_weight[sbi] -> GetXaxis() -> GetBinCenter( hbi ) ;
         }
      } // hbi
      h_max_ldp_weight -> SetBinContent( sbi, max_weight ) ;
      h_max_ldp_weight -> GetXaxis() -> SetBinLabel( sbi, h_ldp -> GetXaxis() -> GetBinLabel( sbi ) ) ;
   } // sbi

   h_max_ldp_weight -> SetFillColor(11) ;
   h_max_ldp_weight -> GetXaxis() -> LabelsOption( "v" ) ;



   printf("\n\n\n Done.\n\n\n") ;

   char histfile[10000] ;
   if ( strlen( samplename ) > 0 ) {
      sprintf( histfile, "outputfiles/hists-v2d-%s.root", samplename ) ;
   } else {
      sprintf( histfile, "outputfiles/hists-v2d.root" ) ;
   }
   saveHist( histfile, "h*" ) ;

} // Loop

//=====================================================================================================

   void fill_hists_loop_v2d::set_bi( ) {

      bi_ht = 0 ; bi_mht = 0 ; bi_nj = 0 ; bi_nb = 0 ;
      for ( int bi = 1; bi <= nb_nj ; bi++ ) { if ( NJets  >= bin_edges_nj[bi-1]  && NJets < bin_edges_nj[bi]  ) { bi_nj  = bi ; break ; } }
      for ( int bi = 1; bi <= nb_nb ; bi++ ) { if ( BTags  >= bin_edges_nb[bi-1]  && BTags < bin_edges_nb[bi]  ) { bi_nb  = bi ; break ; } }
      for ( int bi = 1; bi <= nb_mht; bi++ ) { if ( MHT >= bin_edges_mht[bi-1] && MHT < bin_edges_mht[bi] ) { bi_mht = bi ; break ; } }
      if ( bi_mht>=1 && bi_mht<=nb_mht ) {
         for ( int bi = 1; bi <= nb_ht[bi_mht] ; bi++ ) { if ( HT  >= bin_edges_ht[bi_mht][bi-1]  && HT < bin_edges_ht[bi_mht][bi]  ) { bi_ht  = bi ; break ; } }
      }

      bi_htmht = 0 ;
      if ( bi_ht >=1 && bi_mht >=1 ) {
         bi_htmht = bi_ht ;
         if ( bi_mht >= 2 ) bi_htmht += 3 ;
         if ( bi_mht >= 3 ) bi_htmht += 3 ;
         if ( bi_mht >= 4 ) bi_htmht += 3 ;
         if ( bi_mht >= 5 ) bi_htmht += 2 ;
      }

      bi_global = 0 ;
      bi_nbsum_global = 0 ;
      if ( bi_nj >= 1 && bi_nb >= 1 && bi_htmht >= 1 ) {
         bi_global = (bi_nj-1) * nb_nb * nb_htmht  +  (bi_nb-1) * nb_htmht   + bi_htmht ;
         bi_nbsum_global = (bi_nj-1) * nb_htmht  +  bi_htmht ;
      }

   } // set_bi

//=====================================================================================================

   void fill_hists_loop_v2d::translate_global_bin( int gbi, int& tbi_nj, int& tbi_nb, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) {

      if ( gbi <= 0 ) { printf("\n\n *** translate_global_bin : illegal global bin index %d\n\n", gbi ) ; gSystem -> Exit(-1) ; }

      tbi_nj = (gbi-1) / ( nb_nb*nb_htmht) + 1 ;

      tbi_nb =  ( ((gbi-1)/nb_htmht) + 1 )  % nb_nb  ;
      if ( tbi_nb == 0 ) tbi_nb = nb_nb ;

      tbi_htmht = gbi % nb_htmht ;
      if ( tbi_htmht == 0 ) tbi_htmht = nb_htmht ;

      tbi_ht = ht_bi_from_htmht( tbi_htmht ) ;
      tbi_mht = mht_bi_from_htmht( tbi_htmht ) ;

      //printf("  global bin %3d :  bi_nj = %2d ,  bi_nb = %2d ,  bi_htmht = %2d ,  bi_ht = %2d ,  bi_mht = %d\n",
      //    gbi, tbi_nj, tbi_nb, tbi_htmht, tbi_ht, tbi_mht ) ;

   } // translate_global_bin

//=====================================================================================================

   void fill_hists_loop_v2d::translate_global_bin_nbsum( int gbi_nbsum, int& tbi_nj, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) {

      if ( gbi_nbsum <= 0 ) { printf("\n\n *** translate_global_bin_nbsum : illegal global bin index %d\n\n", gbi_nbsum ) ; gSystem -> Exit(-1) ; }

      tbi_nj = (gbi_nbsum-1) / (nb_htmht) + 1 ;

      tbi_htmht = gbi_nbsum % nb_htmht ;
      if ( tbi_htmht == 0 ) tbi_htmht = nb_htmht ;

      tbi_ht = ht_bi_from_htmht( tbi_htmht ) ;
      tbi_mht = mht_bi_from_htmht( tbi_htmht ) ;

      //printf("  qcd nbsum global bin %3d :  bi_nj = %2d ,  bi_htmht = %2d ,  bi_ht = %2d ,  bi_mht = %d\n",
      //    gbi_nbsum, tbi_nj, tbi_htmht, tbi_ht, tbi_mht ) ;

   } // translate_global_bin_nbsum

//=====================================================================================================

   int  fill_hists_loop_v2d::ht_bi_from_htmht( int abi_htmht ) {

      if ( abi_htmht <= 0 ) { printf("\n\n *** ht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      if ( abi_htmht > nb_htmht ) { printf("\n\n *** ht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      int rv(0) ;
      int sum1(0) ;
      int sum2(0) ;
      for ( int mbi=1; mbi<=nb_mht; mbi++ ) {
         sum2 += nb_ht[mbi] ;
         if ( abi_htmht > sum1 && abi_htmht <= sum2 ) {
            rv = abi_htmht - sum1 ;
            return rv ;
         }
         sum1 = sum2 ;
      }
      return rv ;

   } // ht_bi_from_htmht

   //===================================

   int  fill_hists_loop_v2d::mht_bi_from_htmht( int abi_htmht ) {

      if ( abi_htmht <= 0 ) { printf("\n\n *** mht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      if ( abi_htmht > nb_htmht ) { printf("\n\n *** mht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      int rv(0) ;
      int sum1(0) ;
      int sum2(0) ;
      for ( int mbi=1; mbi<=nb_mht; mbi++ ) {
         sum2 += nb_ht[mbi] ;
         if ( abi_htmht > sum1 && abi_htmht <= sum2 ) {
            rv = mbi ;
            return rv ;
         }
         sum1 = sum2 ;
      }
      return rv ;

   } // ht_bi_from_htmht

   //===================================

   void fill_hists_loop_v2d::set_bin_labels( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      for ( int bi=1; bi<=nb_global/nb_nb; bi++ ) {
         int tbi_nj; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin_nbsum( bi, tbi_nj, tbi_htmht, tbi_ht, tbi_mht ) ;
         char binlabel[100] ;
         //---------
         // sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
         //---------
         int print_bi_htmht = tbi_htmht - nb_ht[1] ;
         int print_bi_mht = tbi_mht - 1 ;
         int print_bi = bi - nb_ht[1]*tbi_nj ;
         int control_bi = (tbi_nj-1)*nb_ht[1] + tbi_ht ;
         if ( tbi_mht == 1 ) {
            sprintf( binlabel, "NJets%d-MHTC-HT%d (C%1d)   C%2d", tbi_nj-1, tbi_ht, tbi_ht, control_bi ) ;
         } else {
            if ( tbi_mht<=3 ) {
               sprintf( binlabel, "NJets%d-MHT%d-HT%d (%2d)   %2d", tbi_nj-1, print_bi_mht, tbi_ht, print_bi_htmht, print_bi ) ;
            } else {
               sprintf( binlabel, "NJets%d-MHT%d-HT%d (%2d)   %2d", tbi_nj-1, print_bi_mht, tbi_ht+1, print_bi_htmht, print_bi ) ;
            }
         }
         //---------
         hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bi

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels

   //===================================

   void fill_hists_loop_v2d::set_bin_labels_div_by_nb( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      for ( int bi=1; bi<=nb_global; bi++ ) {
         int tbi_nj; int tbi_nb; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin( bi, tbi_nj, tbi_nb, tbi_htmht, tbi_ht, tbi_mht ) ;
         char binlabel[100] ;
         //---------
         // sprintf( binlabel, "NJets%d-MHT%d-HT%d (%2d)   %2d", tbi_nj-1, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
         //---------
         int print_bi_htmht = tbi_htmht - nb_ht[1] ;
         int print_bi_mht = tbi_mht - 1 ;
         int print_bi_nb = tbi_nb -1 ;
         int print_bi = (tbi_nj-1)*(nb_htmht-nb_ht[1])*(nb_nb) + (tbi_nb-1)*(nb_htmht-nb_ht[1]) + tbi_htmht -nb_ht[1] ;
         //int control_bi = (tbi_nj-1)*nb_ht[1] + tbi_ht ;
         int control_bi = (tbi_nj-1)*(nb_nb)*(nb_ht[1]) + (tbi_nb-1)*(nb_ht[1]) + tbi_ht ;
         if ( tbi_mht == 1 ) {
            sprintf( binlabel, "NJets%d_BTags%d-MHTC-HT%d (C%1d)   C%2d  %3d", tbi_nj-1, print_bi_nb, tbi_ht, tbi_ht, control_bi, bi ) ;
         } else {
            if ( tbi_mht<=3 ) {
               sprintf( binlabel, "NJets%d_BTags%d-MHT%d-HT%d (%2d)   S%3d  %3d", tbi_nj-1, print_bi_nb, print_bi_mht, tbi_ht, print_bi_htmht, print_bi, bi ) ;
            } else {
               sprintf( binlabel, "NJets%d_BTags%d-MHT%d-HT%d (%2d)   S%3d  %3d", tbi_nj-1, print_bi_nb, print_bi_mht, tbi_ht+1, print_bi_htmht, print_bi, bi ) ;
            }
         }
         //---------
         hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bi

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_div_by_nb

   //===================================

   void fill_hists_loop_v2d::set_bin_labels_mhtc_plot( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      int bi_plot(0) ;
      for ( int tbi_nj=1; tbi_nj<=nb_nj; tbi_nj++ ) {
         for ( int tbi_nb=1; tbi_nb<=nb_nb; tbi_nb++ ) {
            for ( int tbi_ht=1; tbi_ht<=(nb_ht[1]); tbi_ht++ ) {
               bi_plot ++ ;
               char binlabel[100] ;
               sprintf( binlabel, "NJets%d_BTags%d-MHTC-HT%d  C%2d", tbi_nj-1, tbi_nb-1, tbi_ht, bi_plot ) ;
               hp -> GetXaxis() -> SetBinLabel( bi_plot, binlabel ) ;
            } // tbi_ht
         } // tbi_nb
      } // tbi_nj

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_mhtc_plot

   //===================================

   void fill_hists_loop_v2d::set_bin_labels_mhtc_nbsum_plot( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      int bi_plot(0) ;
      for ( int tbi_nj=1; tbi_nj<=nb_nj; tbi_nj++ ) {
            for ( int tbi_ht=1; tbi_ht<=(nb_ht[1]); tbi_ht++ ) {
               bi_plot ++ ;
               char binlabel[100] ;
               sprintf( binlabel, "NJets%d-MHTC-HT%d  C%2d", tbi_nj-1, tbi_ht, bi_plot ) ;
               hp -> GetXaxis() -> SetBinLabel( bi_plot, binlabel ) ;
            } // tbi_ht
      } // tbi_nj

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_mhtc_nbsum_plot

   //===================================

#endif

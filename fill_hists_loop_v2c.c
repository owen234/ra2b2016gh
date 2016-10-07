
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

#define fill_hists_loop_v2c_cxx
#include "fill_hists_loop_v2c.h"
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
#include "histio.c"
#include <iostream>
#include <vector>
#include "lumi_taken.h"

using namespace std ;


void fill_hists_loop_v2c::Loop( bool verb, int nloop )
{

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

   TH1F* h_turnon[5] ;
   if ( isqcd ) {
      h_turnon[1] = (TH1F*) gDirectory -> FindObject( "h_eff_nj1" ) ;
      h_turnon[2] = (TH1F*) gDirectory -> FindObject( "h_eff_nj2" ) ;
      h_turnon[3] = (TH1F*) gDirectory -> FindObject( "h_eff_nj3" ) ;
      h_turnon[4] = (TH1F*) gDirectory -> FindObject( "h_eff_nj4" ) ;
      if ( h_turnon[1] == 0x0 ) { printf("\n\n *** missing nj1 turnon.\n\n") ; return ; }
      if ( h_turnon[2] == 0x0 ) { printf("\n\n *** missing nj2 turnon.\n\n") ; return ; }
      if ( h_turnon[3] == 0x0 ) { printf("\n\n *** missing nj3 turnon.\n\n") ; return ; }
      if ( h_turnon[4] == 0x0 ) { printf("\n\n *** missing nj4 turnon.\n\n") ; return ; }
   }

   setup_bins() ;

   if (fChain == 0) return;


   TH1F* h_hdp = new TH1F( "h_hdp", "HDP events", nb_global, 0.5, nb_global + 0.5 ) ; h_hdp -> Sumw2() ;
   set_bin_labels_div_by_nb( h_hdp ) ;
   TH1F* h_nbsum_hdp = new TH1F( "h_nbsum_hdp", "HDP events", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nbsum_hdp -> Sumw2() ;
   set_bin_labels( h_nbsum_hdp ) ;
   TH1F* h_nb0_hdp = new TH1F( "h_nb0_hdp", "HDP events, Nb=0", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb0_hdp -> Sumw2() ;
   set_bin_labels( h_nb0_hdp ) ;
   TH1F* h_nb1_hdp = new TH1F( "h_nb1_hdp", "HDP events, Nb=1", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb1_hdp -> Sumw2() ;
   set_bin_labels( h_nb1_hdp ) ;
   TH1F* h_nb2_hdp = new TH1F( "h_nb2_hdp", "HDP events, Nb=2", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb2_hdp -> Sumw2() ;
   set_bin_labels( h_nb2_hdp ) ;
   TH1F* h_nb3_hdp = new TH1F( "h_nb3_hdp", "HDP events, Nb=3", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb3_hdp -> Sumw2() ;
   set_bin_labels( h_nb3_hdp ) ;
   TH1F* h_mhtc_hdp = new TH1F( "h_mhtc_hdp", "HDP events, MHTC", nb_nj*nb_nb*nb_ht[1], 0.5, nb_nj*nb_nb*nb_ht[1]+0.5 ) ; h_mhtc_hdp -> Sumw2() ;
   set_bin_labels_mhtc_plot( h_mhtc_hdp ) ;
   TH1F* h_mhtc_nbsum_hdp = new TH1F( "h_mhtc_nbsum_hdp", "HDP events, MHTC, Nbsum", nb_nj*nb_ht[1], 0.5, nb_nj*nb_ht[1]+0.5 ) ; h_mhtc_nbsum_hdp -> Sumw2() ;
   set_bin_labels_mhtc_nbsum_plot( h_mhtc_nbsum_hdp ) ;

   TH1F* h_ldp = new TH1F( "h_ldp", "ldp events", nb_global, 0.5, nb_global + 0.5 ) ; h_ldp -> Sumw2() ;
   set_bin_labels_div_by_nb( h_ldp ) ;
   TH1F* h_nbsum_ldp = new TH1F( "h_nbsum_ldp", "LDP events", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nbsum_ldp -> Sumw2() ;
   set_bin_labels( h_nbsum_ldp ) ;
   TH1F* h_nb0_ldp = new TH1F( "h_nb0_ldp", "ldp events, Nb=0", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb0_ldp -> Sumw2() ;
   set_bin_labels( h_nb0_ldp ) ;
   TH1F* h_nb1_ldp = new TH1F( "h_nb1_ldp", "ldp events, Nb=1", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb1_ldp -> Sumw2() ;
   set_bin_labels( h_nb1_ldp ) ;
   TH1F* h_nb2_ldp = new TH1F( "h_nb2_ldp", "ldp events, Nb=2", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb2_ldp -> Sumw2() ;
   set_bin_labels( h_nb2_ldp ) ;
   TH1F* h_nb3_ldp = new TH1F( "h_nb3_ldp", "ldp events, Nb=3", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb3_ldp -> Sumw2() ;
   set_bin_labels( h_nb3_ldp ) ;
   TH1F* h_mhtc_ldp = new TH1F( "h_mhtc_ldp", "ldp events, MHTC", nb_nj*nb_nb*nb_ht[1], 0.5, nb_nj*nb_nb*nb_ht[1]+0.5 ) ; h_mhtc_ldp -> Sumw2() ;
   set_bin_labels_mhtc_plot( h_mhtc_ldp ) ;
   TH1F* h_mhtc_nbsum_ldp = new TH1F( "h_mhtc_nbsum_ldp", "ldp events, MHTC, Nbsum", nb_nj*nb_ht[1], 0.5, nb_nj*nb_ht[1]+0.5 ) ; h_mhtc_nbsum_ldp -> Sumw2() ;
   set_bin_labels_mhtc_nbsum_plot( h_mhtc_nbsum_ldp ) ;

   TH1F* h_jet_eta_badmu = new TH1F( "h_jet_eta_badmu", "Jet eta, bad muon", 55, -5.5, 5.5 ) ;

   TH1F* h_mht_all = new TH1F( "h_mht_all", "MHT, all events", 100, 0., 1000. ) ; h_mht_all -> Sumw2() ;
   TH1F* h_mht_badmu = new TH1F( "h_mht_badmu", "MHT, badmu", 100, 0., 1000. ) ; h_mht_badmu -> Sumw2() ;
   TH1F* h_mht_met100_calomet80_rejected = new TH1F( "h_mht_met100_calomet80_rejected", "MHT, MET<100 or CaloMET<80", 100, 0., 1000. ) ; h_mht_met100_calomet80_rejected -> Sumw2() ;
   TH1F* h_mht_filters = new TH1F( "h_mht_filters", "MHT, all rejected", 100, 0., 1000. ) ; h_mht_filters -> Sumw2() ;
   TH1F* h_mht_allrejected = new TH1F( "h_mht_allrejected", "MHT, all rejected", 100, 0., 1000. ) ; h_mht_allrejected -> Sumw2() ;

   TH1F* h_met_over_calomet_all = new TH1F( "h_met_over_calomet_all", "MET/CaloMET, all events", 100., 0., 10. ) ; h_met_over_calomet_all -> Sumw2() ;
   TH1F* h_met_over_calomet_badmu = new TH1F( "h_met_over_calomet_badmu", "MET/CaloMET, badmu", 100., 0., 10. ) ; h_met_over_calomet_badmu -> Sumw2() ;
   TH1F* h_met_over_calomet_met100_calomet80_rejected = new TH1F( "h_met_over_calomet_met100_calomet80_rejected", "MET/CaloMET, MET<100 or CaloMET<80", 100., 0., 10. ) ; h_met_over_calomet_met100_calomet80_rejected -> Sumw2() ;
   TH1F* h_met_over_calomet_filters = new TH1F( "h_met_over_calomet_filters", "MET/CaloMET, all rejected", 100., 0., 10. ) ; h_met_over_calomet_filters -> Sumw2() ;
   TH1F* h_met_over_calomet_allrejected = new TH1F( "h_met_over_calomet_allrejected", "MET/CaloMET, all rejected", 100., 0., 10. ) ; h_met_over_calomet_allrejected -> Sumw2() ;

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

   TH1F* h_ldp_weight[209] ;
   for ( int sbi=1; sbi<=208; sbi++ ) {
      char hname[100] ;
      sprintf( hname, "h_ldp_weight_bin%03d", sbi ) ;
      h_ldp_weight[sbi] = new TH1F( hname, hname, 100, 0., 2. ) ;
   }



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


      double hw = Weight * lumi_ ;

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
      if ( CSCTightHaloFilter < 1 )  isjunk = true ;  // branch types are different in 80X samples
      if ( HBHEIsoNoiseFilter < 1 )  isjunk = true ;  // branch types are different in 80X samples
      if ( HBHENoiseFilter < 1 )  isjunk = true ;  // branch types are different in 80X samples

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


      if ( badMuonEvent ) isjunk = true ;
      if ( MET/CaloMET > 5 )  isjunk = true ;

      if ( isjunk ) {
         h_mht_allrejected -> Fill( MHT, hw ) ;
         if ( CaloMET > 0 ) h_met_over_calomet_allrejected -> Fill( MET/CaloMET, hw ) ;
      }


      if ( isjunk ) continue ;

    //-- new baseline
      if ( NJets < 3 ) continue ;
      if ( HT < 300 ) continue ;
      ////////if ( bi_mht <= 0 ) continue ;



     //--- if this is QCD, correct for the trigger turnon
      if ( isqcd ) {
         int turnon_bin(0) ;
         float trig_eff(1.) ;
         if ( bi_nj>=1 && bi_nj<=4 ) {
            if ( MHT < (h_turnon[bi_nj] -> GetXaxis() -> GetXmax()) ) {
               turnon_bin = h_turnon[bi_nj] -> GetXaxis() -> FindBin( MHT ) ;
            } else {
               turnon_bin = h_turnon[bi_nj] -> GetNbinsX() ;
            }
         }
         trig_eff = h_turnon[bi_nj] -> GetBinContent( turnon_bin ) ;
         hw = hw * trig_eff ;
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
         if ( BTags == 0 ) h_nb0_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 1 ) h_nb1_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 2 ) h_nb2_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags >= 3 ) h_nb3_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( bi_mht == 1 ) h_mhtc_hdp -> Fill( bi_mhtc_plot , hw ) ;
         if ( bi_mht == 1 ) h_mhtc_nbsum_hdp -> Fill( bi_mhtc_nbsum_plot , hw ) ;
      } else {
         h_ldp -> Fill( bi_global, hw ) ;
         if ( bi_global > 0 && bi_global < 208 ) h_ldp_weight[bi_global] -> Fill( hw ) ;
         h_nbsum_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 0 ) h_nb0_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 1 ) h_nb1_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 2 ) h_nb2_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags >= 3 ) h_nb3_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( bi_mht == 1 ) h_mhtc_ldp -> Fill( bi_mhtc_plot , hw ) ;
         if ( bi_mht == 1 ) h_mhtc_nbsum_ldp -> Fill( bi_mhtc_nbsum_plot , hw ) ;
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


   TH1F* h_max_ldp_weight = new TH1F( "h_max_ldp_weight", "Max ldp weight", 208, 0.5, 208.5 ) ;

   for ( int sbi=1; sbi<=208; sbi++ ) {
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
      sprintf( histfile, "outputfiles/hists-v2c-%s.root", samplename ) ;
   } else {
      sprintf( histfile, "outputfiles/hists-v2c.root" ) ;
   }
   saveHist( histfile, "h*" ) ;

} // Loop

//=====================================================================================================

   void fill_hists_loop_v2c::setup_bins() {

      int bi ;

      bi = 0 ;
      bin_edges_nj[bi] = 2.5 ; bi++ ;
      bin_edges_nj[bi] = 4.5 ; bi++ ;
      bin_edges_nj[bi] = 6.5 ; bi++ ;
      bin_edges_nj[bi] = 8.5 ; bi++ ;
      bin_edges_nj[bi] = 99.5 ;
      nb_nj = bi ;

      bi = 0 ;
      bin_edges_nb[bi] = -0.5 ; bi++ ;
      bin_edges_nb[bi] =  0.5 ; bi++ ;
      bin_edges_nb[bi] =  1.5 ; bi++ ;
      bin_edges_nb[bi] =  2.5 ; bi++ ;
      bin_edges_nb[bi] =  99.5 ;
      nb_nb = bi ;

      bi = 0 ;
      bin_edges_mht[bi] =   250. ; bi++ ;
      ////////////////////////bin_edges_mht[bi] =   200. ; bi++ ;
      bin_edges_mht[bi] =   300. ; bi++ ;
      bin_edges_mht[bi] =   350. ; bi++ ;
      bin_edges_mht[bi] =   500. ; bi++ ;
      bin_edges_mht[bi] =   750. ; bi++ ;
      bin_edges_mht[bi] = 20000. ;
      nb_mht = bi ;

      int mbi(1) ;
      nb_htmht = 0 ;
      //--- MHT bin 0
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   300. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 1
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   300. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 2
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   350. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 3
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 4
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   750. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1500. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;



      nb_global = nb_htmht * nb_nb * nb_nj ;

      printf("\n\n") ;
      printf("  nj bins : %d\n", nb_nj ) ; for ( int i=1; i<=nb_nj; i++ ) { printf("    nj bin %2d : [%.1f,%.1f]\n", i, bin_edges_nj[i-1], bin_edges_nj[i] ) ; }
      printf("  nb bins : %d\n", nb_nb ) ; for ( int i=1; i<=nb_nb; i++ ) { printf("    nb bin %2d : [%.1f,%.1f]\n", i, bin_edges_nb[i-1], bin_edges_nb[i] ) ; }
      printf("  mht bins : %d\n", nb_mht ) ; for ( int i=1; i<=nb_mht; i++ ) { printf("    mht bin %2d : [%.1f,%.1f]\n", i, bin_edges_mht[i-1], bin_edges_mht[i] ) ; }
      for ( int mbi=1; mbi<=nb_mht; mbi++ ) {
         printf("  ht bins, mht%d : %d\n", mbi, nb_ht[mbi] ) ; for ( int i=1; i<=nb_ht[mbi]; i++ ) { printf("    ht bin %2d : [%.1f,%.1f]\n", i, bin_edges_ht[mbi][i-1], bin_edges_ht[mbi][i] ) ; }
      }
      printf("\n\n") ;

      TH1F* h_nj_bin_index = new TH1F( "h_nj_bin_index", "Njet bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_nb_bin_index = new TH1F( "h_nb_bin_index", "Nb bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_htmht_bin_index = new TH1F( "h_htmht_bin_index", "HTMHT bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_ht_bin_index = new TH1F( "h_ht_bin_index", "HT bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_mht_bin_index = new TH1F( "h_mht_bin_index", "MHT bin index", nb_global, 0.5, nb_global+0.5 ) ;

      printf("\n\n") ;
      for ( int bi=1; bi<=nb_global; bi++ ) {
         int tbi_nj; int tbi_nb; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin( bi, tbi_nj, tbi_nb, tbi_htmht, tbi_ht, tbi_mht ) ;
         h_nj_bin_index -> SetBinContent( bi, tbi_nj ) ;
         h_nb_bin_index -> SetBinContent( bi, tbi_nb ) ;
         h_htmht_bin_index -> SetBinContent( bi, tbi_htmht ) ;
         h_ht_bin_index -> SetBinContent( bi, tbi_ht ) ;
         h_mht_bin_index -> SetBinContent( bi, tbi_mht ) ;
      } // bi


      TH1F* h_nbsum_nj_bin_index = new TH1F( "h_nbsum_nj_bin_index", "Njet bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;
      TH1F* h_nbsum_htmht_bin_index = new TH1F( "h_nbsum_htmht_bin_index", "HTMHT bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;
      TH1F* h_nbsum_ht_bin_index = new TH1F( "h_nbsum_ht_bin_index", "HT bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;
      TH1F* h_nbsum_mht_bin_index = new TH1F( "h_nbsum_mht_bin_index", "MHT bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;

      printf("\n\n") ;
      for ( int bi=1; bi<=nb_global/nb_nb; bi++ ) {
         int tbi_nj; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin_nbsum( bi, tbi_nj, tbi_htmht, tbi_ht, tbi_mht ) ;
         h_nbsum_nj_bin_index -> SetBinContent( bi, tbi_nj ) ;
         h_nbsum_htmht_bin_index -> SetBinContent( bi, tbi_htmht ) ;
         h_nbsum_ht_bin_index -> SetBinContent( bi, tbi_ht ) ;
         h_nbsum_mht_bin_index -> SetBinContent( bi, tbi_mht ) ;
      } // bi

      printf("\n\n") ;

   } // setup_bins

//=====================================================================================================

   void fill_hists_loop_v2c::set_bi( ) {

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

   void fill_hists_loop_v2c::translate_global_bin( int gbi, int& tbi_nj, int& tbi_nb, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) {

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

   void fill_hists_loop_v2c::translate_global_bin_nbsum( int gbi_nbsum, int& tbi_nj, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) {

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

   int  fill_hists_loop_v2c::ht_bi_from_htmht( int abi_htmht ) {

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

   int  fill_hists_loop_v2c::mht_bi_from_htmht( int abi_htmht ) {

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

   void fill_hists_loop_v2c::set_bin_labels( TH1F* hp ) {

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
            sprintf( binlabel, "Nj%d-MHTC-HT%d (C%1d)   C%2d", tbi_nj, tbi_ht, tbi_ht, control_bi ) ;
         } else {
            sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, print_bi_mht, tbi_ht, print_bi_htmht, print_bi ) ;
         }
         //---------
         hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bi

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels

   //===================================

   void fill_hists_loop_v2c::set_bin_labels_div_by_nb( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      for ( int bi=1; bi<=nb_global; bi++ ) {
         int tbi_nj; int tbi_nb; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin( bi, tbi_nj, tbi_nb, tbi_htmht, tbi_ht, tbi_mht ) ;
         char binlabel[100] ;
         //---------
         // sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
         //---------
         int print_bi_htmht = tbi_htmht - nb_ht[1] ;
         int print_bi_mht = tbi_mht - 1 ;
         int print_bi_nb = tbi_nb -1 ;
         int print_bi = (tbi_nj-1)*(nb_htmht-nb_ht[1])*(nb_nb) + (tbi_nb-1)*(nb_htmht-nb_ht[1]) + tbi_htmht -nb_ht[1] ;
         //int control_bi = (tbi_nj-1)*nb_ht[1] + tbi_ht ;
         int control_bi = (tbi_nj-1)*(nb_nb)*(nb_ht[1]) + (tbi_nb-1)*(nb_ht[1]) + tbi_ht ;
         if ( tbi_mht == 1 ) {
            sprintf( binlabel, "Nj%d-Nb%d-MHTC-HT%d (C%1d)   C%2d  %3d", tbi_nj, print_bi_nb, tbi_ht, tbi_ht, control_bi, bi ) ;
         } else {
            sprintf( binlabel, "Nj%d-Nb%d-MHT%d-HT%d (%2d)   S%3d  %3d", tbi_nj, print_bi_nb, print_bi_mht, tbi_ht, print_bi_htmht, print_bi, bi ) ;
         }
         //---------
         hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bi

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_div_by_nb

   //===================================

   void fill_hists_loop_v2c::set_bin_labels_mhtc_plot( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      int bi_plot(0) ;
      for ( int tbi_nj=1; tbi_nj<=nb_nj; tbi_nj++ ) {
         for ( int tbi_nb=1; tbi_nb<=nb_nb; tbi_nb++ ) {
            for ( int tbi_ht=1; tbi_ht<=(nb_ht[1]); tbi_ht++ ) {
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

   void fill_hists_loop_v2c::set_bin_labels_mhtc_nbsum_plot( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      int bi_plot(0) ;
      for ( int tbi_nj=1; tbi_nj<=nb_nj; tbi_nj++ ) {
            for ( int tbi_ht=1; tbi_ht<=(nb_ht[1]); tbi_ht++ ) {
               bi_plot ++ ;
               char binlabel[100] ;
               sprintf( binlabel, "Nj%d-MHTC-HT%d  C%2d", tbi_nj, tbi_ht, bi_plot ) ;
               hp -> GetXaxis() -> SetBinLabel( bi_plot, binlabel ) ;
            } // tbi_ht
      } // tbi_nj

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels_mhtc_nbsum_plot

   //===================================



#ifndef make_lostlep_input_files2_c 
#define make_lostlep_input_files2_c

#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "binning.h"
#include "get_hist.h"


   void make_lostlep_input_files2(
                                 const char* ldp_input_root_file = "non-qcd-inputs-fall16c/LLPrediction_QCDLDP_notCombined.root",
                                 const char* hdp_input_root_file = "non-qcd-inputs-fall16c/LLPrediction_QCDHDP_notCombined.root",
                                 const char* output_text_file = "outputfiles/combine-input-lostlep.txt",
                                 const char* nbsum_text_file = "outputfiles/nbsum-input-lostlep.txt"
                               ) {

      setup_bins();

      //bool verb(false) ;
      bool verb(true) ;

      gDirectory -> Delete( "h*" ) ;

      TFile* tf_ldp = new TFile( ldp_input_root_file, "read" ) ;
      if ( tf_ldp == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", ldp_input_root_file ) ; return ; }
      if ( !(tf_ldp -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", ldp_input_root_file ) ; return ; }

      TFile* tf_hdp = new TFile( hdp_input_root_file, "read" ) ;
      if ( tf_hdp == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", hdp_input_root_file ) ; return ; }
      if ( !(tf_hdp -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", hdp_input_root_file ) ; return ; }

      tf_ldp -> cd( "Prediction_data" ) ;
      tf_hdp -> cd( "Prediction_data" ) ;

      printf("\n") ;
      tf_ldp -> ls() ;
      printf("\n") ;

      FILE* ofp_combine ;
      if ( (ofp_combine = fopen( output_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening combine output file: %s\n\n", output_text_file ) ;
         return ;
      }

      FILE* ofp_nbsum ;
      if ( (ofp_nbsum = fopen( nbsum_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum output file: %s\n\n", nbsum_text_file ) ;
         return ;
      }

      FILE* ofp_nbsum_stat_syst ;
      TString systfile_nbsum( nbsum_text_file ) ;
      systfile_nbsum.ReplaceAll("input","stat-syst") ;
      if ( (ofp_nbsum_stat_syst = fopen( systfile_nbsum.Data(), "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum stat-syst output file: %s\n\n", systfile_nbsum.Data() ) ;
         return ;
      }

      FILE* ofp_combine_stat_syst ;
      TString systfile_combine( output_text_file ) ;
      systfile_combine.ReplaceAll("input","stat-syst") ;
      if ( (ofp_combine_stat_syst = fopen( systfile_combine.Data(), "w" ))==NULL ) {
         printf( "\n\n *** Problem opening combine stat-syst output file: %s\n\n", systfile_combine.Data() ) ;
         return ;
      }

      TH1F* h_ldp = (TH1F*) tf_ldp -> Get( "Prediction_data/totalPred_LL" ) ;
      if ( h_ldp == 0x0 ) { printf("\n\n *** Missing totalPred_LL for LDP\n\n") ; return ; }

      TH1F* h_hdp = (TH1F*) tf_hdp -> Get( "Prediction_data/totalPred_LL" ) ;
      if ( h_hdp == 0x0 ) { printf("\n\n *** Missing totalPred_LL for HDP\n\n") ; return ; }

      TH1* h_systerr_lowdphi[100] ;
      TH1* h_systerr_highdphi[100] ;
      char systerr_name[100][100] ;
      int n_systerr(0) ;

      {
         int si(0) ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredIsoTrackSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredIsoTrackSysUp_LL" ) ;
         sprintf( systerr_name[si], "IsoTrk" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMTWSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMTWSysUp_LL" ) ;
         sprintf( systerr_name[si], "MTW" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredPuritySysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredPuritySysUp_LL" ) ;
         sprintf( systerr_name[si], "Purity" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredSingleLepPuritySysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredSingleLepPuritySysUp_LL" ) ;
         sprintf( systerr_name[si], "SingleLepPurity" ) ;
         si++ ;

         ///////////h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredDiLepFoundSysUp_LL" ) ;
         ///////////h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredDiLepFoundSysUp_LL" ) ;
         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredDiLepSRSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredDiLepSRSysUp_LL" ) ;
         sprintf( systerr_name[si], "DiLepFound" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuIsoSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuIsoSysUp_LL" ) ;
         sprintf( systerr_name[si], "MuIso" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuRecoSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuRecoSysUp_LL" ) ;
         sprintf( systerr_name[si], "MuReco" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuAccSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuAccSysUp_LL" ) ;
         sprintf( systerr_name[si], "MuAcc" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuAccQsquareSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuAccQsquareSysUp_LL" ) ;
         sprintf( systerr_name[si], "MuAccQ2" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecIsoSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecIsoSysUp_LL" ) ;
         sprintf( systerr_name[si], "ElecIso" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecRecoSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecRecoSysUp_LL" ) ;
         sprintf( systerr_name[si], "ElecReco" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecAccSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecAccSysUp_LL" ) ;
         sprintf( systerr_name[si], "ElecAcc" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecAccQsquareSysUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecAccQsquareSysUp_LL" ) ;
         sprintf( systerr_name[si], "ElecAccQ2" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredNonClosureUp_LL" ) ;
         h_systerr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredNonClosureUp_LL" ) ;
         sprintf( systerr_name[si], "NonClosure" ) ;
         si++ ;


         n_systerr = si ;
      }




      TH1* h_staterr_lowdphi[100] ;
      TH1* h_staterr_highdphi[100] ;
      char staterr_name[100][100] ;
      int n_staterr(0) ;

      {
         int si(0) ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredIsoTrackStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredIsoTrackStatUp_LL" ) ;
         sprintf( staterr_name[si], "IsoTrk" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMTWStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMTWStatUp_LL" ) ;
         sprintf( staterr_name[si], "MTW" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredPurityStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredPurityStatUp_LL" ) ;
         sprintf( staterr_name[si], "Purity" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredSingleLepPurityStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredSingleLepPurityStatUp_LL" ) ;
         sprintf( staterr_name[si], "SingleLepPurity" ) ;
         si++ ;

         //////////////////////h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredDiLepFoundStatUp_LL" ) ;
         //////////////////////h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredDiLepFoundStatUp_LL" ) ;
         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredDiLepSRStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredDiLepSRStatUp_LL" ) ;
         sprintf( staterr_name[si], "DiLepFound" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuIsoStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuIsoStatUp_LL" ) ;
         sprintf( staterr_name[si], "MuIso" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuRecoStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuRecoStatUp_LL" ) ;
         sprintf( staterr_name[si], "MuReco" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredMuAccStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredMuAccStatUp_LL" ) ;
         sprintf( staterr_name[si], "MuAcc" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecIsoStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecIsoStatUp_LL" ) ;
         sprintf( staterr_name[si], "ElecIso" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecRecoStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecRecoStatUp_LL" ) ;
         sprintf( staterr_name[si], "ElecReco" ) ;
         si++ ;

         h_staterr_lowdphi[si]  = get_hist( tf_ldp, "Prediction_data/totalPredElecAccStatUp_LL" ) ;
         h_staterr_highdphi[si]  = get_hist( tf_hdp, "Prediction_data/totalPredElecAccStatUp_LL" ) ;
         sprintf( staterr_name[si], "ElecAcc" ) ;
         si++ ;

         n_staterr = si ;
      }










      int bi_hist(0) ;
      int bi_control(0) ;
      int bi_search(0) ;
      int count_all(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               count_all++ ;

               int bi_hist = global_bin_with_mhtc( bi_nj, bi_nb, bi_htmht ) ;

               double ldp_val(0.) ;
               double ldp_hist_err(0.) ;

               double hdp_val(0.) ;
               double hdp_hist_err(0.) ;

               if ( bi_hist > 0 ) {
                  ldp_val = h_ldp -> GetBinContent( bi_hist ) ;
                  ldp_hist_err = h_ldp -> GetBinError( bi_hist ) ;

                  hdp_val = h_hdp -> GetBinContent( bi_hist ) ;
                  hdp_hist_err = h_hdp -> GetBinError( bi_hist ) ;
               }



               double total_syst_err2_lowdphi(0.) ;
               double total_syst_err2_highdphi(0.) ;
               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_lowdphi(0.) ;
                  double syst_highdphi(0.) ;
                  if ( bi_hist > 0 ) {
                     syst_lowdphi = h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                     syst_highdphi = h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  }
                  double syst_lowdphi_events(0.) ;
                  double syst_highdphi_events(0.) ;
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     syst_lowdphi_events  = syst_lowdphi  * ldp_val ;
                     syst_highdphi_events = syst_highdphi * hdp_val ;
                  }
                  total_syst_err2_lowdphi  += pow( syst_lowdphi_events , 2. ) ;
                  total_syst_err2_highdphi += pow( syst_highdphi_events, 2. ) ;
               } // si

               double total_stat_err2_lowdphi(0.) ;
               double total_stat_err2_highdphi(0.) ;
               for ( int si=0; si<n_staterr; si++ ) {
                  double stat_lowdphi(0.) ;
                  double stat_highdphi(0.) ;
                  if ( bi_hist > 0 ) {
                     stat_lowdphi = h_staterr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                     stat_highdphi = h_staterr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  }
                  double stat_lowdphi_events(0.) ;
                  double stat_highdphi_events(0.) ;
                  if ( stat_lowdphi > 0 && stat_highdphi > 0 ) {
                     stat_lowdphi_events  = stat_lowdphi  * ldp_val ;
                     stat_highdphi_events = stat_highdphi * hdp_val ;
                  }
                  total_stat_err2_lowdphi  += pow( stat_lowdphi_events , 2. ) ;
                  total_stat_err2_highdphi += pow( stat_highdphi_events, 2. ) ;
                  printf(" LDP stat err : %12s %6.3f  %8.1f\n", staterr_name[si], stat_lowdphi, stat_lowdphi_events ) ;
               } // si

               printf(" LDP total stat err :  %8.1f   %8.1f\n", sqrt(total_stat_err2_lowdphi), ldp_hist_err  ) ;
               ///////////////double combined_stat_err_lowdphi = sqrt( pow( total_stat_err2_lowdphi, 2. ) + pow( ldp_hist_err, 2. ) ) ; // BUG
               double combined_stat_err_lowdphi = sqrt( total_stat_err2_lowdphi + pow( ldp_hist_err, 2. ) ) ;
               double total_syst_err_lowdphi = sqrt( total_syst_err2_lowdphi ) ;

               ///////////double combined_stat_err_highdphi = sqrt( pow( total_stat_err2_highdphi, 2. ) + pow( hdp_hist_err, 2. ) ) ; // BUG
               double combined_stat_err_highdphi = sqrt( total_stat_err2_highdphi + pow( hdp_hist_err, 2. ) ) ;
               double total_syst_err_highdphi = sqrt( total_syst_err2_highdphi ) ;





               TString hist_bin_label( h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               int bi_ht = 0, bi_mht = 0;

               if ( bi_htmht == 1 ) { bi_ht = 1; bi_mht = 1; }
               if ( bi_htmht == 2 ) { bi_ht = 2; bi_mht = 1; }
               if ( bi_htmht == 3 ) { bi_ht = 3; bi_mht = 1; }

               if ( bi_htmht == 4 ) { bi_ht = 1; bi_mht = 2; }
               if ( bi_htmht == 5 ) { bi_ht = 2; bi_mht = 2; }
               if ( bi_htmht == 6 ) { bi_ht = 3; bi_mht = 2; }

               if ( bi_htmht == 7 ) { bi_ht = 1; bi_mht = 3; }
               if ( bi_htmht == 8 ) { bi_ht = 2; bi_mht = 3; }
               if ( bi_htmht == 9 ) { bi_ht = 3; bi_mht = 3; }

               if ( bi_htmht ==10 ) { bi_ht = 2; bi_mht = 4; }
               if ( bi_htmht ==11 ) { bi_ht = 3; bi_mht = 4; }

               if ( bi_htmht ==12 ) { bi_ht = 2; bi_mht = 5; }
               if ( bi_htmht ==13 ) { bi_ht = 3; bi_mht = 5; }

               if ( bi_mht == 1 ) {
                  bi_control ++ ;
               } else {
                  bi_search ++ ;
               }

               char mhtchar[10] ;
               if ( bi_mht == 1 ) {
                  sprintf( mhtchar, "C" ) ;
               } else {
                  sprintf( mhtchar, "%d", bi_mht-1 ) ;
               }

               char label[1000] ;
               sprintf( label, " %3d  %s %3d  Nj%d-Nb%d-MHT%s-HT%d",
                   count_all, (bi_mht==1)?"C":"S", (bi_mht==1)?bi_control:bi_search,
                   bi_nj, bi_nb-1, mhtchar, bi_ht ) ;

               //// printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;

               double ldp_stat_over_sqrtn(0.8), ldp_syst_over_n(0.15) ;
               if ( ldp_val > 0 ) {
                  ldp_stat_over_sqrtn = combined_stat_err_lowdphi / sqrt( ldp_val ) ;
                  ldp_syst_over_n = total_syst_err_lowdphi / ldp_val ;
               }
               double hdp_stat_over_sqrtn(0.8), hdp_syst_over_n(0.15) ;
               if ( hdp_val > 0 ) {
                  hdp_stat_over_sqrtn = combined_stat_err_highdphi / sqrt( hdp_val ) ;
                  hdp_syst_over_n = total_syst_err_highdphi / hdp_val ;
               }

               printf(               "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, combined_stat_err_lowdphi, total_syst_err_lowdphi,   hdp_val, combined_stat_err_highdphi, total_syst_err_highdphi ) ;

               fprintf( ofp_combine, "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, combined_stat_err_lowdphi, total_syst_err_lowdphi,   hdp_val, combined_stat_err_highdphi, total_syst_err_highdphi ) ;

               if ( bi_htmht > 3 ) {
                  char this_label[100] ;
                  sprintf( this_label, "%3d S-Nj%d-Nb%d-MHT%d-HT%d (S%d)", bi_search, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht-3 ) ;
                  fprintf( ofp_combine_stat_syst, "   %-30s   %6.3f  %6.3f       %6.3f  %6.3f\n",
                      this_label,
                      ldp_stat_over_sqrtn, ldp_syst_over_n,   hdp_stat_over_sqrtn, hdp_syst_over_n ) ;
               } else {
                  char this_label[100] ;
                  sprintf( this_label, "%3d C-Nj%d-Nb%d-MHT%d-HT%d (C%d)", bi_control, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht ) ;
                  fprintf( ofp_combine_stat_syst, "   %-30s   %6.3f  %6.3f       %6.3f  %6.3f\n",
                      this_label,
                      ldp_stat_over_sqrtn, ldp_syst_over_n,   hdp_stat_over_sqrtn, hdp_syst_over_n ) ;
               }

            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp_combine ) ;
      printf("\n\n Wrote %s\n\n", output_text_file ) ;

     //---------

      for ( int bi_ht=1; bi_ht<=nBinsHT; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            float nbsum_lowdphi_val(0.) ;
            float nbsum_lowdphi_err2(0.) ;

            float nbsum_highdphi_val(0.) ;
            float nbsum_highdphi_err2(0.) ;

            double total_syst_lowdphi_events[100] ;
            double total_syst_highdphi_events[100] ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_lowdphi_events[si]  = 0 ;
               total_syst_highdphi_events[si] = 0 ;
            }
            double total_stat2_lowdphi_events[100] ;
            double total_stat2_highdphi_events[100] ;
            for ( int si=0; si<n_staterr; si++ ) {
               total_stat2_lowdphi_events[si]  = 0 ;
               total_stat2_highdphi_events[si] = 0 ;
            }

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               int bi_hist = global_bin_with_mhtc( bi_nj, bi_nb, bi_ht ) ;

               double lowdphi_val(0.) ;
               double lowdphi_err(0.) ;
               if ( bi_hist > 0 ) {
                  lowdphi_val = h_ldp -> GetBinContent( bi_hist ) ;
                  lowdphi_err = h_ldp -> GetBinError( bi_hist ) ;
               }

               nbsum_lowdphi_val += lowdphi_val ;
               nbsum_lowdphi_err2 += pow( lowdphi_err, 2. ) ;

               double highdphi_val = h_hdp -> GetBinContent( bi_hist ) ;
               double highdphi_err = h_hdp -> GetBinError( bi_hist ) ;

               nbsum_highdphi_val += highdphi_val ;
               nbsum_highdphi_err2 += pow( highdphi_err, 2. ) ;

               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_lowdphi(0.) ;
                  double syst_highdphi(0.) ;
                  if ( bi_hist > 0 ) {
                     syst_lowdphi = h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                     syst_highdphi = h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  }
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst_lowdphi_events[si]  += syst_lowdphi  * lowdphi_val ;
                     total_syst_highdphi_events[si] += syst_highdphi * highdphi_val ;
                  }
               } // si

               for ( int si=0; si<n_staterr; si++ ) {
                  double stat_lowdphi(0.) ;
                  double stat_highdphi(0.) ;
                  if ( bi_hist > 0 ) {
                     stat_lowdphi = h_staterr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                     stat_highdphi = h_staterr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  }
                  if ( stat_lowdphi > 0 && stat_highdphi > 0 ) {
                     total_stat2_lowdphi_events[si]  += pow( stat_lowdphi  * lowdphi_val, 2 ) ;
                     total_stat2_highdphi_events[si] += pow( stat_highdphi * highdphi_val, 2 ) ;
                  }
               } // si
            } // bi_nb

            double total_syst_err2_lowdphi(0.) ;
            double total_syst_err2_highdphi(0.) ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_err2_lowdphi += pow( total_syst_lowdphi_events[si], 2.)   ;
               total_syst_err2_highdphi += pow( total_syst_highdphi_events[si], 2.)  ;
               if (verb) printf(" syst %2d : %25s :   %7.1f  %7.1f\n", si, systerr_name[si], total_syst_lowdphi_events[si], total_syst_highdphi_events[si] ) ;
            }
            double total_stat_err2_lowdphi(0.) ;
            double total_stat_err2_highdphi(0.) ;
            for ( int si=0; si<n_staterr; si++ ) {
               total_stat_err2_lowdphi +=  total_stat2_lowdphi_events[si]  ;
               total_stat_err2_highdphi +=  total_stat2_highdphi_events[si]  ;
               if (verb) printf(" stat %2d : %25s :   %7.1f  %7.1f\n", si, staterr_name[si], sqrt(total_stat2_lowdphi_events[si]), sqrt(total_stat2_highdphi_events[si]) ) ;
            }
            double total_syst_err_lowdphi = sqrt( total_syst_err2_lowdphi ) ;
            double total_syst_err_highdphi = sqrt( total_syst_err2_highdphi ) ;
            double total_stat_err_lowdphi = sqrt( nbsum_lowdphi_err2 + total_stat_err2_lowdphi ) ;
            double total_stat_err_highdphi = sqrt( nbsum_highdphi_err2 + total_stat_err2_highdphi ) ;

            double ldp_syst_over_n(0.5) ;
            double ldp_stat_over_sqrtn(0.5) ;
            if ( nbsum_lowdphi_val > 0 ) {
               ldp_syst_over_n =  total_syst_err_lowdphi / nbsum_lowdphi_val  ;
               ldp_stat_over_sqrtn = total_stat_err_lowdphi / sqrt(nbsum_lowdphi_val) ;
            }
            double hdp_syst_over_n(0.5) ;
            double hdp_stat_over_sqrtn(0.5) ;
            if ( nbsum_highdphi_val > 0 ) {
               hdp_syst_over_n =  total_syst_err_highdphi / nbsum_highdphi_val  ;
               hdp_stat_over_sqrtn = total_stat_err_highdphi / sqrt(nbsum_highdphi_val) ;
            }

            printf( "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f (%6.3f)     %8.1f +/- %5.1f +/- %5.1f (%6.3f)\n", bi_nj, bi_ht,
                  nbsum_lowdphi_val, total_stat_err_lowdphi, total_syst_err_lowdphi,  ldp_syst_over_n,
                  nbsum_highdphi_val, total_stat_err_highdphi, total_syst_err_highdphi,  hdp_syst_over_n ) ;
            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f    %8.1f +/- %5.1f +/- %5.1f\n", bi_nj, bi_ht,
                  nbsum_lowdphi_val, total_stat_err_lowdphi, total_syst_err_lowdphi,
                  nbsum_highdphi_val, total_stat_err_highdphi, total_syst_err_highdphi ) ;
            fprintf( ofp_nbsum_stat_syst, "   Nj%d-HT%d    %6.3f  %6.3f       %6.3f  %6.3f\n", bi_nj, bi_ht, 
                ldp_stat_over_sqrtn, ldp_syst_over_n,   hdp_stat_over_sqrtn, hdp_syst_over_n ) ;

         } // bi_nj
      } // bi_ht

      fclose( ofp_nbsum ) ;
      fclose( ofp_nbsum_stat_syst ) ;
      fclose( ofp_combine_stat_syst ) ;
      printf("\n\n Wrote %s\n\n", nbsum_text_file ) ;
      printf("\n\n Wrote %s\n\n", systfile_nbsum.Data() ) ;
      printf("\n\n Wrote %s\n\n", systfile_combine.Data() ) ;



   } // make_lostlep_input_files2


#endif

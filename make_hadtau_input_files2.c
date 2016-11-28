#ifndef make_hadtau_input_files2_c
#define make_hadtau_input_files2_c

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include "binning.h"
#include "get_hist.h"

   bool  hadtau_gbi_array_ready(false) ;
   int   hadtau_gbi_with_mhtc[6][5][14] ;
   int hadtau_global_bin_with_mhtc( int nji, int nbi, int htmhti ) ;

   //------

   void make_hadtau_input_files2( const char* input_root_file  = "non-qcd-inputs-fall16a/HadTauEstimation_data_formatted.root",
                                  const char* output_text_file = "outputfiles/combine-input-hadtau.txt",
                                  const char* nbsum_text_file  = "outputfiles/nbsum-input-hadtau.txt"
                               ) {

      setup_bins();

      gDirectory -> Delete( "h*" ) ;

      TFile* tf_input = new TFile( input_root_file, "read" ) ;
      if ( tf_input == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }
      if ( !(tf_input -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }

      printf("\n") ;
      tf_input -> ls() ;
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


      TH1* h_pred_lowdphi = get_hist( tf_input, "QCDBin_LowDphi_nominal" ) ;

      TH1* h_pred_highdphi = get_hist( tf_input, "QCDBin_HiDphi_nominal" ) ;




      TH1* h_systerr_lowdphi[100] ;
      TH1* h_systerr_highdphi[100] ;
      char systerr_name[100][100] ;
      int n_systerr(0) ;

      {
         int si(0) ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_closureUncertainty" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_closureUncertainty" ) ;
         sprintf( systerr_name[si], "Closure" ) ;
         si++ ;

     //----------

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_BMistagUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_BMistagUp" ) ;
         sprintf( systerr_name[si], "Btag" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuRecoSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuRecoSysUp" ) ;
         sprintf( systerr_name[si], "MuRecoSys" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuIsoSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuIsoSysUp" ) ;
         sprintf( systerr_name[si], "MuIsoSys" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuRecoIsoUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuRecoIsoUp" ) ;
         sprintf( systerr_name[si], "MuRecoIso" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_JECSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_JECSysUp" ) ;
         sprintf( systerr_name[si], "JEC" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MTSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MTSysUp" ) ;
         sprintf( systerr_name[si], "MT" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MtEffStat" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MtEffStat" ) ;
         sprintf( systerr_name[si], "MtEffStat" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_IsoTrkVetoEffUncertaintySys" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_IsoTrkVetoEffUncertaintySys" ) ;
         sprintf( systerr_name[si], "IsoTrkVetoSys" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_AccStat" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_AccStat" ) ;
         sprintf( systerr_name[si], "AccStat" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_AccSysPDFUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_AccSysPDFUp" ) ;
         sprintf( systerr_name[si], "AccSysPDF" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_AccSysScaleUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_AccSysScaleUp" ) ;
         sprintf( systerr_name[si], "AccSysScale" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuFromTauStat" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuFromTauStat" ) ;
         sprintf( systerr_name[si], "MuFromTauStat" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_DileptonUncertainty" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_DileptonUncertainty" ) ;
         sprintf( systerr_name[si], "Dilepton" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_TrigEffUncertainty" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_TrigEffUncertainty" ) ;
         sprintf( systerr_name[si], "TrigEff" ) ;
         si++ ;

         n_systerr = si ;
      }




      int bi_hist(0) ;
      int bi_control(0) ;
      int bi_search(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               int bi_hist = hadtau_global_bin_with_mhtc( bi_nj, bi_nb, bi_htmht ) ;
               if ( bi_hist > 0 ) {
                  int hbl_nj, hbl_nb, hbl_mht, hbl_htmht ;
                  sscanf( h_pred_lowdphi->GetXaxis()->GetBinLabel( bi_hist ), "NJets%d_BTags%d_MHT%d_HT%d", &hbl_nj, &hbl_nb, &hbl_mht, &hbl_htmht ) ;
                  if ( bi_nj != (hbl_nj+1) ) { printf("\n\n*** Inconsistent Njets.  %d != %d\n", bi_nj, hbl_nj+1) ; gSystem -> Exit(-1) ; }
                  if ( bi_nb != (hbl_nb+1) ) { printf("\n\n*** Inconsistent Nb. %d != %d\n", bi_nb, hbl_nb+1 ) ; gSystem -> Exit(-1) ; }
               }

               double ldp_val(0.) ;
               double ldp_hist_err(0.) ;
               if ( bi_hist > 0 ) {
                  ldp_val = h_pred_lowdphi -> GetBinContent( bi_hist ) ;
                  ldp_hist_err = h_pred_lowdphi -> GetBinError( bi_hist ) ;
               }
               double ldp_err = 1. ;
               if ( ldp_val > 0 ) ldp_err = ldp_hist_err ;

               double hdp_val(0.) ;
               double hdp_hist_err(0.) ;
               if ( bi_hist > 0 ) {
                  hdp_val = h_pred_highdphi -> GetBinContent( bi_hist ) ;
                  hdp_hist_err = h_pred_highdphi -> GetBinError( bi_hist ) ;
               }
               double hdp_err = 1. ;
               if ( hdp_val > 0 ) hdp_err = hdp_hist_err ;


               double total_syst2_lowdphi_events(0.) ;
               double total_syst2_highdphi_events(0.) ;
               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_lowdphi(0.) ;
                  double syst_highdphi(0.) ;
                  if ( bi_hist > 0 ) {
                     syst_lowdphi = fabs( h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                     syst_highdphi = fabs( h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                  }
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst2_lowdphi_events  += pow( syst_lowdphi  * ldp_val, 2. ) ;
                     total_syst2_highdphi_events += pow( syst_highdphi * hdp_val, 2. ) ;
                  }
               } // si

               double ldp_syst = sqrt( total_syst2_lowdphi_events ) ;
               double hdp_syst = sqrt( total_syst2_highdphi_events ) ;

               TString hist_bin_label ;
               if ( bi_hist > 0 ) { hist_bin_label = h_pred_lowdphi -> GetXaxis() -> GetBinLabel( bi_hist ) ; }

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
                   bi_hist, (bi_mht==1)?"C":"S", (bi_mht==1)?bi_control:bi_search,
                   bi_nj, bi_nb-1, mhtchar, bi_ht ) ;

               //// printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;

               double ldp_stat_over_sqrtn(0.8), ldp_syst_over_n(0.15) ;
               if ( ldp_val > 0 ) {
                  ldp_stat_over_sqrtn = ldp_err / sqrt( ldp_val ) ;
                  ldp_syst_over_n = ldp_syst / ldp_val ;
                  if ( ldp_stat_over_sqrtn > 1.5 ) ldp_stat_over_sqrtn = 1.0 ;
               }
               double hdp_stat_over_sqrtn(0.8), hdp_syst_over_n(0.15) ;
               if ( hdp_val > 0 ) {
                  hdp_stat_over_sqrtn = hdp_err / sqrt( hdp_val ) ;
                  hdp_syst_over_n = hdp_syst / hdp_val ;
                  if ( hdp_stat_over_sqrtn > 1.5 ) hdp_stat_over_sqrtn = 1.0 ;
               }

               printf(               "%s      %8.1f +/- %5.1f +/- %5.1f          %8.1f +/- %5.1f +/- %5.1f\n",
                   label,
                   ldp_val, ldp_err, ldp_syst,
                   hdp_val, hdp_err, hdp_syst ) ;

               fprintf( ofp_combine, "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_err, ldp_syst,   hdp_val, hdp_err, hdp_syst ) ;

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

     //-----

      for ( int bi_ht=1; bi_ht<=nBinsHT; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            float ldp_nbsum_val(0.) ;
            float hdp_nbsum_val(0.) ;

            float ldp_nbsum_stat_err2(0.) ;
            float hdp_nbsum_stat_err2(0.) ;

            double total_syst_lowdphi_events[100] ;
            double total_syst_highdphi_events[100] ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_lowdphi_events[si]  = 0 ;
               total_syst_highdphi_events[si] = 0 ;
            }

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               bi_hist = hadtau_global_bin_with_mhtc( bi_nj, bi_nb, bi_ht ) ;

               if ( bi_hist > 0 ) {
                  int hbl_nj, hbl_nb, hbl_mht, hbl_htmht ;
                  sscanf( h_pred_lowdphi->GetXaxis()->GetBinLabel( bi_hist ), "NJets%d_BTags%d_MHT%d_HT%d", &hbl_nj, &hbl_nb, &hbl_mht, &hbl_htmht ) ;
                  if ( bi_nj != (hbl_nj+1) ) { printf("\n\n*** Inconsistent Njets.  %d != %d\n", bi_nj, hbl_nj+1) ; gSystem -> Exit(-1) ; }
                  if ( bi_nb != (hbl_nb+1) ) { printf("\n\n*** Inconsistent Nb. %d != %d\n", bi_nb, hbl_nb+1 ) ; gSystem -> Exit(-1) ; }
               }

               double ldp_val(0.) ;
               double ldp_hist_err(0.) ;
               if ( bi_hist > 0 ) {
                  ldp_val = h_pred_lowdphi -> GetBinContent( bi_hist ) ;
                  ldp_hist_err = h_pred_lowdphi -> GetBinError( bi_hist ) ;
               }
               double ldp_err = 1. ;
               if ( ldp_val > 0 ) ldp_err = ldp_hist_err ;

               double hdp_val(0.) ;
               double hdp_hist_err (0.) ;
               if ( bi_hist > 0 ) {
                  hdp_val = h_pred_highdphi -> GetBinContent( bi_hist ) ;
                  hdp_hist_err = h_pred_highdphi -> GetBinError( bi_hist ) ;
               }
               double hdp_err = 1. ;
               if ( hdp_val > 0 ) hdp_err = hdp_hist_err ;


               //////////////double ldp_syst = 0.10 * ldp_val ; //*** flat 10% syst for now.
               //////////////double hdp_syst = 0.10 * hdp_val ; //*** flat 10% syst for now.

               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_lowdphi(0.) ;
                  double syst_highdphi(0.) ;
                  if ( bi_hist > 0 ) {
                     syst_lowdphi = fabs( h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                     syst_highdphi = fabs( h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                  }
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst_lowdphi_events[si]  += syst_lowdphi  * ldp_val ;
                     total_syst_highdphi_events[si] += syst_highdphi * hdp_val ;
                  }
               } // si

               ldp_nbsum_val += ldp_val ;
               hdp_nbsum_val += hdp_val ;

               ldp_nbsum_stat_err2 += pow( ldp_err, 2. ) ;
               hdp_nbsum_stat_err2 += pow( hdp_err, 2. ) ;

               TString hist_bin_label ;
               if ( bi_hist > 0 ) { hist_bin_label = h_pred_lowdphi -> GetXaxis() -> GetBinLabel( bi_hist ) ; }

               char label[1000] ;
               sprintf( label, " %3d  Nj%d-Nb%d-MHTC-HT%d", bi_hist, bi_nj, bi_nb-1, bi_ht ) ;

            } // bi_nb

            ///printf("    Systematics:\n") ;
            double total_syst_err2_lowdphi(0.) ;
            double total_syst_err2_highdphi(0.) ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_err2_lowdphi += pow( total_syst_lowdphi_events[si], 2.)   ;
               total_syst_err2_highdphi += pow( total_syst_highdphi_events[si], 2.)  ;
               ///printf("  %2d : %25s :   %7.1f  %7.1f\n", si, systerr_name[si], total_syst_lowdphi_events[si], total_syst_highdphi_events[si] ) ;
            }
            double total_syst_err_lowdphi = sqrt( total_syst_err2_lowdphi ) ;
            double total_syst_err_highdphi = sqrt( total_syst_err2_highdphi ) ;


            float ldp_nbsum_stat_err = sqrt( ldp_nbsum_stat_err2 ) ;
            float hdp_nbsum_stat_err = sqrt( hdp_nbsum_stat_err2 ) ;

            double ldp_syst_over_n(0.) ;
            double ldp_stat_over_sqrtn(0.) ;
            if ( ldp_nbsum_val > 0 ) ldp_syst_over_n = total_syst_err_lowdphi / ldp_nbsum_val  ;
            if ( ldp_nbsum_val > 0 ) ldp_stat_over_sqrtn = ldp_nbsum_stat_err / sqrt(ldp_nbsum_val)  ;

            double hdp_syst_over_n(0.) ;
            double hdp_stat_over_sqrtn(0.) ;
            if ( hdp_nbsum_val > 0 ) hdp_syst_over_n = total_syst_err_highdphi / hdp_nbsum_val  ;
            if ( hdp_nbsum_val > 0 ) hdp_stat_over_sqrtn = hdp_nbsum_stat_err / sqrt(hdp_nbsum_val)  ;

            printf(          "   Nj%d-HT%d      %8.1f +/- %5.1f +/- %5.1f  (%6.3f)        %8.1f +/- %5.1f +/- %5.1f   (%6.3f)\n",
                bi_nj, bi_ht,
                ldp_nbsum_val, ldp_nbsum_stat_err, total_syst_err_lowdphi,    ldp_syst_over_n,
                hdp_nbsum_val, hdp_nbsum_stat_err, total_syst_err_highdphi,   hdp_syst_over_n ) ;

            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f    %8.1f +/- %5.1f +/- %5.1f\n", bi_nj, bi_ht,
                  ldp_nbsum_val, ldp_nbsum_stat_err, total_syst_err_lowdphi,
                  hdp_nbsum_val, hdp_nbsum_stat_err, total_syst_err_highdphi   ) ;
            fprintf( ofp_nbsum_stat_syst, "   Nj%d-HT%d    %6.3f  %6.3f       %6.3f  %6.3f\n", bi_nj, bi_ht, 
                  ldp_stat_over_sqrtn, ldp_syst_over_n,    hdp_stat_over_sqrtn, hdp_syst_over_n ) ;

         } // bi_nj
      } // bi_ht

      fclose( ofp_nbsum ) ;
      fclose( ofp_nbsum_stat_syst ) ;
      printf("\n\n Wrote %s\n\n", nbsum_text_file ) ;
      printf("\n\n Wrote %s\n\n", systfile_nbsum.Data() ) ;
      printf("\n\n Wrote %s\n\n", systfile_combine.Data() ) ;



   } // make_hadtau_input_files2

//=========================================

int hadtau_global_bin_with_mhtc ( int arg_nji, int arg_nbi, int arg_htmhti ) {

   if ( !hadtau_gbi_array_ready ) {
      int gbi(0) ;
      int gbi_no_mhtc(0) ;      
      for ( int nji=1; nji<=nb_nj; nji++ ) {
         int nbmax=nb_nb ;
         if ( nji==1 ) nbmax -- ;
         for ( int nbi=1; nbi<=nbmax; nbi ++ ) {
            for ( int htmhti=1; htmhti<=nb_htmht; htmhti++ ) {
               int hti, mhti ;
               translate_htmht_bin_to_ht_and_mht_bins( htmhti, hti, mhti ) ;
               bool excluded = is_this_bin_excluded( nji-1, nbi-1, hti-1, mhti-1 ) ;
               if ( mhti==1 ) {
                  gbi++ ;
               } else {
                  if ( !excluded ) gbi++ ;
               }
               if ( !excluded ) {
                  hadtau_gbi_with_mhtc[nji][nbi][htmhti] = gbi ;
               } else {
                  hadtau_gbi_with_mhtc[nji][nbi][htmhti] = -1 ;
               } // else if !excluded
            } //htmhti
         } // nbi
      } // mhi
      hadtau_gbi_array_ready = true ;
   }

   return hadtau_gbi_with_mhtc[arg_nji][arg_nbi][arg_htmhti] ;

} // hadtau_global_bin_with_mhtc

//=========================================



#endif

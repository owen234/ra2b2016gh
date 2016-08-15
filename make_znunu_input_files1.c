#ifndef make_znunu_input_files1_c
#define make_znunu_input_files1_c

#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "get_hist.h"
#include "binning.h"

   void make_znunu_input_files1( const char* ldp_input_root_file = "non-qcd-inputs-topup2/ZinvHistos_ldp.root",
                                 const char* hdp_input_root_file = "non-qcd-inputs-topup2/ZinvHistos_hdp.root",
                                 const char* output_text_file = "outputfiles/combine-input-znunu.txt",
                                 const char* nbsum_text_file  = "outputfiles/nbsum-input-znunu.txt"
                               ) {
      setup_bins();


      gDirectory -> Delete( "h*" ) ;

      TFile* tf_ldp = new TFile( ldp_input_root_file, "read" ) ;
      if ( tf_ldp == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", ldp_input_root_file ) ; return ; }
      if ( !(tf_ldp -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", ldp_input_root_file ) ; return ; }

      TFile* tf_hdp = new TFile( hdp_input_root_file, "read" ) ;
      if ( tf_hdp == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", hdp_input_root_file ) ; return ; }
      if ( !(tf_hdp -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", hdp_input_root_file ) ; return ; }

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


      TH1F* h_ldp_nonzero = (TH1F*) tf_ldp -> Get( "ZinvBGpred" ) ;
      if ( h_ldp_nonzero == 0x0 ) { printf("\n\n *** Missing ZinvBGpred\n\n") ; return ; }

      TH1F* h_hdp_nonzero = (TH1F*) tf_hdp -> Get( "ZinvBGpred" ) ;
      if ( h_hdp_nonzero == 0x0 ) { printf("\n\n *** Missing ZinvBGpred\n\n") ; return ; }



      TH1F* h_ldp_zero = (TH1F*) tf_ldp -> Get( "ZinvBG0EVpred" ) ;
      if ( h_ldp_zero == 0x0 ) { printf("\n\n *** Missing ZinvBG0EVpred\n\n") ; return ; }

      TH1F* h_hdp_zero = (TH1F*) tf_hdp -> Get( "ZinvBG0EVpred" ) ;
      if ( h_hdp_zero == 0x0 ) { printf("\n\n *** Missing ZinvBG0EVpred\n\n") ; return ; }

      TH1F* h_ldp = (TH1F*) h_ldp_nonzero -> Clone( "h_ldp" ) ;
      h_ldp -> Add( h_ldp_zero ) ;

      TH1F* h_hdp = (TH1F*) h_hdp_nonzero -> Clone( "h_hdp" ) ;
      h_hdp -> Add( h_hdp_zero ) ;

      TH1* h_systerr_lowdphi[100] ;
      TH1* h_systerr_highdphi[100] ;
      char systerr_name[100][100] ;
      int n_systerr(0) ;

      {
         int si(0) ;

         TH1* h_syst_ldp_nonzero ;
         TH1* h_syst_ldp_zero ;
         TH1* h_syst_ldp ;

         TH1* h_syst_hdp_nonzero ;
         TH1* h_syst_hdp_zero ;
         TH1* h_syst_hdp ;

        //--- Note: These are in events.

         h_syst_ldp_nonzero = get_hist( tf_ldp, "ZinvBGsysUp" ) ;
         h_syst_ldp_zero = get_hist( tf_ldp, "ZinvBG0EVsysUp" ) ;
         h_syst_ldp = (TH1*) h_syst_ldp_nonzero -> Clone( "ZinvBGsysUp_combined" ) ;
         h_syst_ldp -> Add( h_syst_ldp_zero ) ;

         h_systerr_lowdphi[si] = h_syst_ldp ;

         h_syst_hdp_nonzero = get_hist( tf_hdp, "ZinvBGsysUp" ) ;
         h_syst_hdp_zero = get_hist( tf_hdp, "ZinvBG0EVsysUp" ) ;
         h_syst_hdp = (TH1*) h_syst_hdp_nonzero -> Clone( "ZinvBGsysUp_combined" ) ;
         h_syst_hdp -> Add( h_syst_hdp_zero ) ;

         h_systerr_highdphi[si] = h_syst_hdp ;


         sprintf( systerr_name[si], "Total" ) ;
         si++ ;

         n_systerr = si ;
      }







      int bi_hist(0) ;
      int bi_control(0) ;
      int bi_search(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               bi_hist ++ ;

               double ldp_val = h_ldp -> GetBinContent( bi_hist ) ;
               double ldp_hist_err = h_ldp -> GetBinError( bi_hist ) ;

               double hdp_val = h_hdp -> GetBinContent( bi_hist ) ;
               double hdp_hist_err = h_hdp -> GetBinError( bi_hist ) ;

               double total_syst_lowdphi_events[100] ;
               double total_syst_highdphi_events[100] ;
               for ( int si=0; si<n_systerr; si++ ) {
                  total_syst_lowdphi_events[si]  = 0 ;
                  total_syst_highdphi_events[si] = 0 ;
               }
               for ( int si=0; si<n_systerr; si++ ) {
                  ////////////////double syst_lowdphi = h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  ////////////////double syst_highdphi = h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  double syst_lowdphi = h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) ;
                  double syst_highdphi = h_systerr_highdphi[si] -> GetBinContent( bi_hist ) ;
                  ///////printf(" DEBUG1 syst hist content:  ldp %9.3f  hdp %9.3f\n", 
                      ///////h_systerr_lowdphi[si] -> GetBinContent( bi_hist ),
                      ///////h_systerr_highdphi[si] -> GetBinContent( bi_hist ) ) ;
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst_lowdphi_events[si]  += syst_lowdphi  ;
                     total_syst_highdphi_events[si] += syst_highdphi  ;
                  }
               } // si
               double total_syst_err2_lowdphi(0.) ;
               double total_syst_err2_highdphi(0.) ;
               for ( int si=0; si<n_systerr; si++ ) {
                  total_syst_err2_lowdphi += pow( total_syst_lowdphi_events[si], 2.)   ;
                  total_syst_err2_highdphi += pow( total_syst_highdphi_events[si], 2.)  ;
                  //////printf("  %2d : %25s :   %7.1f  %7.1f\n", si, systerr_name[si], total_syst_lowdphi_events[si], total_syst_highdphi_events[si] ) ;
               }
               double total_syst_err_lowdphi = sqrt( total_syst_err2_lowdphi ) ;
               double total_syst_err_highdphi = sqrt( total_syst_err2_highdphi ) ;
               //////////printf("  DEBUG2 syst ldp %9.3f   syst hdp %9.3f\n", total_syst_err_lowdphi, total_syst_err_highdphi ) ;

               TString hist_bin_label( h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               int bi_ht, bi_mht ;

               if ( bi_htmht == 1 ) { bi_ht = 1; bi_mht = 1; }
               if ( bi_htmht == 2 ) { bi_ht = 2; bi_mht = 1; }
               if ( bi_htmht == 3 ) { bi_ht = 3; bi_mht = 1; }

               if ( bi_htmht == 4 ) { bi_ht = 1; bi_mht = 2; }
               if ( bi_htmht == 5 ) { bi_ht = 2; bi_mht = 2; }
               if ( bi_htmht == 6 ) { bi_ht = 3; bi_mht = 2; }

               if ( bi_htmht == 7 ) { bi_ht = 1; bi_mht = 3; }
               if ( bi_htmht == 8 ) { bi_ht = 2; bi_mht = 3; }
               if ( bi_htmht == 9 ) { bi_ht = 3; bi_mht = 3; }

               /////////////////if ( bi_htmht ==10 ) { bi_ht = 1; bi_mht = 4; }
               /////////////////if ( bi_htmht ==11 ) { bi_ht = 2; bi_mht = 4; }

               /////////////////if ( bi_htmht ==12 ) { bi_ht = 1; bi_mht = 5; }
               /////////////////if ( bi_htmht ==13 ) { bi_ht = 2; bi_mht = 5; }

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

               double ldp_stat_over_sqrtn(0.8), ldp_syst_over_n(0.20) ;
               if ( ldp_val > 0 ) {
                  ldp_stat_over_sqrtn = ldp_hist_err / sqrt( ldp_val ) ;
                  if ( total_syst_err_lowdphi > 0 ) ldp_syst_over_n = total_syst_err_lowdphi / ldp_val ;
               }
               double hdp_stat_over_sqrtn(0.8), hdp_syst_over_n(0.20) ;
               if ( hdp_val > 0 ) {
                  hdp_stat_over_sqrtn = hdp_hist_err / sqrt( hdp_val ) ;
                  if ( total_syst_err_highdphi ) hdp_syst_over_n = total_syst_err_highdphi / hdp_val ;
               }

               printf(               "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_hist_err, total_syst_err_lowdphi,   hdp_val, hdp_hist_err, total_syst_err_highdphi ) ;

               fprintf( ofp_combine, "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_hist_err, total_syst_err_lowdphi,   hdp_val, hdp_hist_err, total_syst_err_highdphi ) ;

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

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               bi_hist = (bi_nj-1)*(nb_nb)*(nb_htmht) + (bi_nb-1)*(nb_htmht) + bi_ht ;

               double ldp_val = h_ldp -> GetBinContent( bi_hist ) ;
               double ldp_hist_err = h_ldp -> GetBinError( bi_hist ) ;

               double hdp_val = h_hdp -> GetBinContent( bi_hist ) ;
               double hdp_hist_err = h_hdp -> GetBinError( bi_hist ) ;


               nbsum_lowdphi_val += ldp_val ;
               nbsum_highdphi_val += hdp_val ;

               nbsum_lowdphi_err2 += pow( ldp_hist_err, 2. ) ;
               nbsum_highdphi_err2 += pow( hdp_hist_err, 2. ) ;

               for ( int si=0; si<n_systerr; si++ ) {
                  ////////double syst_lowdphi = h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  ////////double syst_highdphi = h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ;
                  double syst_lowdphi = h_systerr_lowdphi[si] -> GetBinContent( bi_hist )  ;
                  double syst_highdphi = h_systerr_highdphi[si] -> GetBinContent( bi_hist ) ;
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst_lowdphi_events[si]  += syst_lowdphi  ;
                     total_syst_highdphi_events[si] += syst_highdphi  ;
                  }
               } // si


               TString hist_bin_label( h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               char label[1000] ;
               sprintf( label, " %3d  Nj%d-Nb%d-MHTC-HT%d", bi_hist, bi_nj, bi_nb-1, bi_ht ) ;

       ////    printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;


            } // bi_nb

            double total_syst_err2_lowdphi(0.) ;
            double total_syst_err2_highdphi(0.) ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_err2_lowdphi += pow( total_syst_lowdphi_events[si], 2.)   ;
               total_syst_err2_highdphi += pow( total_syst_highdphi_events[si], 2.)  ;
               //printf("  %2d : %25s :   %7.1f  %7.1f\n", si, systerr_name[si], total_syst_lowdphi_events[si], total_syst_highdphi_events[si] ) ;
            }
            double total_syst_err_lowdphi = sqrt( total_syst_err2_lowdphi ) ;
            double total_syst_err_highdphi = sqrt( total_syst_err2_highdphi ) ;

            double ldp_syst_over_n(0.2) ;
            double ldp_stat_over_sqrtn(0.8) ;
            if ( nbsum_lowdphi_val > 0 ) {
               ldp_syst_over_n =  total_syst_err_lowdphi / nbsum_lowdphi_val  ;
               ldp_stat_over_sqrtn = 1. ;
            }
            double hdp_syst_over_n(0.2) ;
            double hdp_stat_over_sqrtn(0.8) ;
            if ( nbsum_highdphi_val > 0 ) {
               hdp_syst_over_n =  total_syst_err_highdphi / nbsum_highdphi_val  ;
               hdp_stat_over_sqrtn = 1. ;
            }
            printf( "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f  (%6.3f)      %8.1f +/- %5.1f +/- %5.1f  (%6.3f)\n", bi_nj, bi_ht,
                  nbsum_lowdphi_val, sqrt(nbsum_lowdphi_err2), total_syst_err_lowdphi,  ldp_syst_over_n,
                  nbsum_highdphi_val, sqrt(nbsum_highdphi_err2), total_syst_err_highdphi,  hdp_syst_over_n ) ;
            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f    %8.1f +/- %5.1f +/- %5.1f\n", bi_nj, bi_ht,
                  nbsum_lowdphi_val, sqrt(nbsum_lowdphi_err2), total_syst_err_lowdphi,
                  nbsum_highdphi_val, sqrt(nbsum_highdphi_err2), total_syst_err_highdphi ) ;
            fprintf( ofp_nbsum_stat_syst, "   Nj%d-HT%d    %6.3f  %6.3f      %6.3f  %6.3f\n", bi_nj, bi_ht,
                   ldp_stat_over_sqrtn, ldp_syst_over_n,      hdp_stat_over_sqrtn, hdp_syst_over_n ) ;

         } // bi_nj
      } // bi_ht

      fclose( ofp_nbsum ) ;
      fclose( ofp_nbsum_stat_syst ) ;
      printf("\n\n Wrote %s\n\n", nbsum_text_file ) ;
      printf("\n\n Wrote %s\n\n", systfile_nbsum.Data() ) ;
      printf("\n\n Wrote %s\n\n", systfile_combine.Data() ) ;



   } // make_znunu_input_files1

#endif

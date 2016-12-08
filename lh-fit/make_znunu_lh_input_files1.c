#ifndef make_znunu_lh_input_files1_c
#define make_znunu_lh_input_files1_c

#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "../get_hist.h"
#include "../binning.h"

   void make_znunu_lh_input_files1( const char* ldp_input_root_file = "../non-qcd-inputs-fall16b/ZinvHistos_ldp.root",
                                 const char* hdp_input_root_file = "../non-qcd-inputs-fall16b/ZinvHistos_hdp.root",
                                 const char* output_text_file = "outputfiles/combine-input-znunu.txt",
                                 const char* nbsum_ldp_text_file  = "outputfiles/nbsum-ldp-input-znunu.txt",
                                 const char* nbsum_hdp_text_file  = "outputfiles/nbsum-hdp-input-znunu.txt"
                               ) {
      setup_bins();

      gSystem -> Exec( "mkdir -p outputfiles" ) ;

      double trig_eff[6] = {0., 0.982, 0.985, 0.995, 1.00, 1.00 } ; // inclusive MET trig eff, index 1 = MHTC.
                                                                    // see email from Aditee on Nov 3, 2016.


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


      FILE* ofp_ldp_nbsum ;
      if ( (ofp_ldp_nbsum = fopen( nbsum_ldp_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum LDP output file: %s\n\n", nbsum_ldp_text_file ) ;
         return ;
      }

      FILE* ofp_hdp_nbsum ;
      if ( (ofp_hdp_nbsum = fopen( nbsum_hdp_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum HDP output file: %s\n\n", nbsum_hdp_text_file ) ;
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

      ////////////printf("\n\n *** Debug dump of LDP hist:\n") ;
      ////////////h_ldp -> Print("all") ;

      TH1F* h_ldp_total_syst = get_hist( tf_ldp, "ZinvBGsysUp" ) ;
      TH1F* h_hdp_total_syst = get_hist( tf_hdp, "ZinvBGsysUp" ) ;

      TH1* h_systerr_ldp[100] ;
      TH1* h_systerr_hdp[100] ;
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


         h_systerr_ldp[si] = get_hist( tf_ldp, "hgJstat" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hgJstat" ) ;
         sprintf( systerr_name[si], "hgJstat" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvDYstat" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvDYstat" ) ;
         sprintf( systerr_name[si], "hzvvDYstat" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvDYMCstat" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvDYMCstat" ) ;
         sprintf( systerr_name[si], "hzvvDYMCstat" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hgJZgRerr" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hgJZgRerr" ) ;
         sprintf( systerr_name[si], "hgJZgRerr" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvgJPurErr" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvgJPurErr" ) ;
         sprintf( systerr_name[si], "hzvvgJPurErr" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvZgDRerrUp" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvZgDRerrUp" ) ;
         sprintf( systerr_name[si], "hzvvZgDRerrUp" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvScaleErr" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvScaleErr" ) ;
         sprintf( systerr_name[si], "hzvvScaleErr" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvDYsysNjUp" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvDYsysNjUp" ) ;
         sprintf( systerr_name[si], "hzvvDYsysNjUp" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvDYsysKin" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvDYsysKin" ) ;
         sprintf( systerr_name[si], "hzvvDYsysKin" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvDYsysPur" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvDYsysPur" ) ;
         sprintf( systerr_name[si], "hzvvDYsysPur" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvgJEtrgErr" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvgJEtrgErr" ) ;
         sprintf( systerr_name[si], "hzvvgJEtrgErr" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hgJSFerr" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hgJSFerr" ) ;
         sprintf( systerr_name[si], "hgJSFerr" ) ;
         si++ ;
         n_systerr = si ;

         h_systerr_ldp[si] = get_hist( tf_ldp, "hzvvgJFdirErr" ) ;
         h_systerr_hdp[si] = get_hist( tf_hdp, "hzvvgJFdirErr" ) ;
         sprintf( systerr_name[si], "hzvvgJFdirErr" ) ;
         si++ ;
         n_systerr = si ;


      }



      printf("\n\n\n") ;
      printf(" ===== List of syst hists:\n") ;
      for ( int si=0; si<n_systerr; si++ ) {
         printf("  %2d : %25s\n", si, systerr_name[si] ) ;
      } // si
      printf("\n\n\n") ;





      int bi_hist(0) ;
      int bi_control(0) ;
      int bi_search(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               bi_hist = global_bin_with_mhtc( bi_nj, bi_nb, bi_htmht ) ;

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

               double total_syst_ldp_events[100] ;
               double total_syst_hdp_events[100] ;
               for ( int si=0; si<n_systerr; si++ ) {
                  total_syst_ldp_events[si]  = 0 ;
                  total_syst_hdp_events[si] = 0 ;
               }
               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_ldp(0.) ;
                  double syst_hdp(0.) ;
                  if ( bi_hist > 0 ) {
                     syst_ldp = h_systerr_ldp[si] -> GetBinContent( bi_hist ) ;
                     syst_hdp = h_systerr_hdp[si] -> GetBinContent( bi_hist ) ;
                  }
                  if ( syst_ldp > 0 && syst_hdp > 0 ) {
                     ///////total_syst_ldp_events[si]  += syst_ldp  ;
                     ///////total_syst_hdp_events[si] += syst_hdp  ;
                     total_syst_ldp_events[si]  += (syst_ldp-1)*ldp_val  ;
                     total_syst_hdp_events[si] += (syst_hdp-1)*hdp_val  ;
                  }
               } // si
               double total_syst_err2_ldp(0.) ;
               double total_syst_err2_hdp(0.) ;
               for ( int si=0; si<n_systerr; si++ ) {
                  total_syst_err2_ldp += pow( total_syst_ldp_events[si], 2.)   ;
                  total_syst_err2_hdp += pow( total_syst_hdp_events[si], 2.)  ;
                  printf("  %2d : %25s :   %7.1f  %7.1f\n", si, systerr_name[si], total_syst_ldp_events[si], total_syst_hdp_events[si] ) ;
               }
               double total_syst_err_ldp = sqrt( total_syst_err2_ldp ) ;
               double total_syst_err_hdp = sqrt( total_syst_err2_hdp ) ;
               printf("  DEBUG2      syst ldp %9.3f   syst hdp %9.3f\n", total_syst_err_ldp, total_syst_err_hdp ) ;
               printf("  DEBUG3 hist syst ldp %9.3f        hdp %9.3f\n", 
                    h_ldp_total_syst->GetBinContent( bi_hist ), h_hdp_total_syst->GetBinContent( bi_hist ) ) ;
               printf("\n") ;

               TString hist_bin_label ;
               if ( bi_hist > 0 ) { hist_bin_label = h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ; }

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

               double ldp_stat_over_sqrtn(0.8), ldp_syst_over_n(0.20) ;
               if ( ldp_val > 0 ) {
                  ldp_stat_over_sqrtn = ldp_hist_err / sqrt( ldp_val ) ;
                  if ( total_syst_err_ldp > 0 ) ldp_syst_over_n = total_syst_err_ldp / ldp_val ;
               }
               double hdp_stat_over_sqrtn(0.8), hdp_syst_over_n(0.20) ;
               if ( hdp_val > 0 ) {
                  hdp_stat_over_sqrtn = hdp_hist_err / sqrt( hdp_val ) ;
                  if ( total_syst_err_hdp ) hdp_syst_over_n = total_syst_err_hdp / hdp_val ;
               }

               ldp_val                 = ldp_val                 * trig_eff[bi_mht] ;
               ldp_hist_err            = ldp_hist_err            * trig_eff[bi_mht] ;
               total_syst_err_ldp  = total_syst_err_ldp  * trig_eff[bi_mht] ;
               hdp_val                 = hdp_val                 * trig_eff[bi_mht] ;
               hdp_hist_err            = hdp_hist_err            * trig_eff[bi_mht] ;
               total_syst_err_hdp = total_syst_err_hdp * trig_eff[bi_mht] ;


               printf(               "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_hist_err, total_syst_err_ldp,   hdp_val, hdp_hist_err, total_syst_err_hdp ) ;

               fprintf( ofp_combine, "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_hist_err, total_syst_err_ldp,   hdp_val, hdp_hist_err, total_syst_err_hdp ) ;



            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp_combine ) ;
      printf("\n\n Wrote %s\n\n", output_text_file ) ;











     //-----

      fprintf( ofp_ldp_nbsum, "  znunu ldp nsyst=%d\n", n_systerr ) ;
      fprintf( ofp_hdp_nbsum, "  znunu hdp nsyst=%d\n", n_systerr ) ;

      fprintf( ofp_ldp_nbsum, "  bin label    value   stat  " ) ;
      fprintf( ofp_hdp_nbsum, "  bin label    value   stat  " ) ;
      fprintf( ofp_ldp_nbsum, " %15s ", "total-syst" ) ;
      fprintf( ofp_hdp_nbsum, " %15s ", "total-syst" ) ;
      for ( int si=0; si<n_systerr; si++ ) {
         fprintf( ofp_ldp_nbsum, " %15s ", systerr_name[si] ) ;
         fprintf( ofp_hdp_nbsum, " %15s ", systerr_name[si] ) ;
      } // si
      fprintf( ofp_ldp_nbsum, "\n") ;
      fprintf( ofp_hdp_nbsum, "\n") ;

      for ( int bi_ht=1; bi_ht<=nBinsHT; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            if ( bi_ht==1 && bi_nj>3 ) continue ;


            float nbsum_ldp_val(0.) ;
            float nbsum_ldp_err2(0.) ;

            float nbsum_hdp_val(0.) ;
            float nbsum_hdp_err2(0.) ;

            double total_syst_ldp_events[100] ;
            double total_syst_hdp_events[100] ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_ldp_events[si]  = 0 ;
               total_syst_hdp_events[si] = 0 ;
            }

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               bi_hist = global_bin_with_mhtc( bi_nj, bi_nb, bi_ht ) ;

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


               nbsum_ldp_val += ldp_val ;
               nbsum_hdp_val += hdp_val ;

               nbsum_ldp_err2 += pow( ldp_hist_err, 2. ) ;
               nbsum_hdp_err2 += pow( hdp_hist_err, 2. ) ;

               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_ldp(0.) ;
                  double syst_hdp(0.) ;
                  if ( bi_hist > 0 ) {
                     syst_ldp = ((h_systerr_ldp[si] -> GetBinContent( bi_hist )) - 1.)*ldp_val  ;
                     syst_hdp = ((h_systerr_hdp[si] -> GetBinContent( bi_hist )) - 1.)*hdp_val ;
                  }
                  if ( syst_ldp > 0 && syst_hdp > 0 ) {
                     total_syst_ldp_events[si]  += syst_ldp  ;
                     total_syst_hdp_events[si] += syst_hdp  ;
                  }
               } // si


               TString hist_bin_label ;
               if ( bi_hist > 0 ) { hist_bin_label =  h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ; }

               char label[1000] ;
               sprintf( label, " %3d  Nj%d-Nb%d-MHTC-HT%d", bi_hist, bi_nj, bi_nb-1, bi_ht ) ;

       ////    printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;


            } // bi_nb

            nbsum_ldp_val       = nbsum_ldp_val        * trig_eff[1] ;
            nbsum_ldp_err2      = nbsum_ldp_err2       * pow(trig_eff[1],2) ;
            nbsum_hdp_val      = nbsum_hdp_val       * trig_eff[1] ;
            nbsum_hdp_err2     = nbsum_hdp_err2      * pow(trig_eff[1],2) ;

            fprintf( ofp_ldp_nbsum, "   Nj%d-HT%d   %8.1f   %5.1f", bi_nj, bi_ht, nbsum_ldp_val, sqrt(nbsum_ldp_err2) ) ;
            fprintf( ofp_hdp_nbsum, "   Nj%d-HT%d   %8.1f   %5.1f", bi_nj, bi_ht, nbsum_hdp_val, sqrt(nbsum_hdp_err2) ) ;

            double total_syst_err2_ldp(0.) ;
            double total_syst_err2_hdp(0.) ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_err2_ldp += pow( total_syst_ldp_events[si], 2.)   ;
               total_syst_err2_hdp += pow( total_syst_hdp_events[si], 2.)  ;
            }
            double total_syst_err_ldp = sqrt( total_syst_err2_ldp ) ;
            double total_syst_err_hdp = sqrt( total_syst_err2_hdp ) ;
            double rel_ldp_total_err(0.) ;
            if ( nbsum_ldp_val > 0 ) rel_ldp_total_err = total_syst_err_ldp / nbsum_ldp_val ;
            double rel_hdp_total_err(0.) ;
            if ( nbsum_hdp_val > 0 ) rel_hdp_total_err = total_syst_err_hdp / nbsum_hdp_val ;
            fprintf( ofp_ldp_nbsum, "  %5.1f (%7.4f)", total_syst_err_ldp, rel_ldp_total_err ) ;
            fprintf( ofp_hdp_nbsum, "  %5.1f (%7.4f)", total_syst_err_hdp, rel_hdp_total_err ) ;

            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_ldp_events[si] = total_syst_ldp_events[si] * trig_eff[1] ;
               total_syst_hdp_events[si] = total_syst_hdp_events[si] * trig_eff[1] ;
               double rel_ldp_err(0.) ;
               if ( nbsum_ldp_val > 0 ) { rel_ldp_err = total_syst_ldp_events[si] / nbsum_ldp_val ; }
               double rel_hdp_err(0.) ;
               if ( nbsum_hdp_val > 0 ) { rel_hdp_err = total_syst_hdp_events[si] / nbsum_hdp_val ; }
               fprintf( ofp_ldp_nbsum, "  %5.1f (%7.4f)", total_syst_ldp_events[si], rel_ldp_err ) ;
               fprintf( ofp_hdp_nbsum, "  %5.1f (%7.4f)", total_syst_hdp_events[si], rel_hdp_err ) ;
            }

            fprintf( ofp_ldp_nbsum, "\n" ) ;
            fprintf( ofp_hdp_nbsum, "\n" ) ;

         } // bi_nj
      } // bi_ht

      fclose( ofp_ldp_nbsum ) ;
      fclose( ofp_hdp_nbsum ) ;
      printf("\n\n Wrote %s\n\n", nbsum_ldp_text_file ) ;



   } // make_znunu_lh_input_files1
#endif

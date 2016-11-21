#ifndef make_qcdmc_input_files1_c
#define make_qcdmc_input_files1_c

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"

#include "binning.h"

   void make_qcdmc_input_files1( const char* input_root_file = "outputfiles/hists-v2d-qcd.root",
                                 const char* nbsum_text_file = "outputfiles/nbsum-input-qcd.txt",
                                 const char* output_hist_file = "outputfiles/modelfit-input-qcdmc.root"
                               ) {

      setup_bins();
      gDirectory -> Delete( "h*" ) ;

      TFile* tf = new TFile( input_root_file, "read" ) ;
      if ( tf == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }
      if ( !(tf -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }

//      printf("\n") ;
//      tf -> ls() ;
//      printf("\n") ;

      FILE* ofp_nbsum ;
      if ( (ofp_nbsum = fopen( nbsum_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum output file: %s\n\n", nbsum_text_file ) ;
         return ;
      }


      FILE* ofp_pars ;
      if ( (ofp_pars = fopen( "outputfiles/kqcd-parameters-from-qcdmc.txt", "w" ))==NULL ) {
         printf("\n\n *** Problem opening outputfiles/kqcd-parameters-from-qcdmc.txt\n\n") ; gSystem -> Exit(-1) ;
      }

      TH1F* h_ldp = (TH1F*) tf -> Get( "h_ldp" ) ;
      if ( h_ldp == 0x0 ) { printf("\n\n *** Missing h_ldp\n\n") ; return ; }

      TH1F* h_hdp = (TH1F*) tf -> Get( "h_hdp" ) ;
      if ( h_hdp == 0x0 ) { printf("\n\n *** Missing h_hdp\n\n") ; return ; }

//      int nb_nj(5) ;
//      int nb_nb(4) ;
//      int nb_htmht(13) ;

      int bi_hist(0) ;

      TH1F* h_ratio = new TH1F( "h_ratio", "H/L ratio", nBinsHT*nb_nj, 0.5, nBinsHT*nb_nj + 0.5 ) ;

      int bi_ratio_hist(0) ;

      for ( int bi_ht=1; bi_ht<=3; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            double ldp_nbsum_val(0.) ;
            double ldp_nbsum_err2(0.) ;

            double hdp_nbsum_val(0.) ;
            double hdp_nbsum_err2(0.) ;

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               bi_hist = (bi_nj-1)*(nb_nb)*(nb_htmht) + (bi_nb-1)*(nb_htmht) + bi_ht ;

               double ldp_val = h_ldp -> GetBinContent( bi_hist ) ;
               double ldp_hist_err = h_ldp -> GetBinError( bi_hist ) ;

               double hdp_val = h_hdp -> GetBinContent( bi_hist ) ;
               double hdp_hist_err = h_hdp -> GetBinError( bi_hist ) ;

               ldp_nbsum_val += ldp_val ;
               hdp_nbsum_val += hdp_val ;

               ldp_nbsum_err2 += pow( ldp_hist_err, 2. ) ;
               hdp_nbsum_err2 += pow( hdp_hist_err, 2. ) ;

               TString hist_bin_label( h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               char label[1000] ;
               sprintf( label, " %3d  Nj%d-Nb%d-MHTC-HT%d", bi_hist, bi_nj, bi_nb-1, bi_ht ) ;

       ////    printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;




            } // bi_nb

            double ldp_nbsum_err = sqrt( ldp_nbsum_err2 ) ;
            double hdp_nbsum_err = sqrt( hdp_nbsum_err2 ) ;

            double ratio_val(0.) ;
            double ratio_err(0.) ;

            if ( ldp_nbsum_val > 0. ) {
               ratio_val = hdp_nbsum_val / ldp_nbsum_val ;
               if ( hdp_nbsum_val > 0. ) {
                  ratio_err = ratio_val * sqrt( pow( ldp_nbsum_err/ldp_nbsum_val, 2. ) + pow( hdp_nbsum_err/hdp_nbsum_val, 2. ) ) ;
               }
            }

            printf( "   Nj%d-HT%d   %8.1f +/- %5.1f     %8.1f +/- %5.1f      R(H/L) = %6.3f +/- %6.3f\n", 
                bi_nj, bi_ht, ldp_nbsum_val, ldp_nbsum_err,  hdp_nbsum_val, hdp_nbsum_err, ratio_val, ratio_err ) ;

            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f +/- %5.1f     %8.1f +/- %5.1f      R(H/L) = %6.3f +/- %6.3f\n", 
                bi_nj, bi_ht, ldp_nbsum_val, ldp_nbsum_err,  hdp_nbsum_val, hdp_nbsum_err, ratio_val, ratio_err ) ;

            char label[100] ;
            sprintf( label, "Nj%d-HT%d", bi_nj, bi_ht ) ;

            bool excluded = is_this_bin_excluded( bi_nj-1, 0, bi_ht-1, 0 ) ;
            if ( !excluded ) {
               fprintf( ofp_pars, "R_qcd_ldp_Nj%d_HT%d %6.4f +/- %6.4f\n", bi_nj, bi_ht, ratio_val, ratio_err ) ;
            }

            bi_ratio_hist ++ ;
            if ( !(bi_ht==1 && bi_nj>nb_nj-2) ) {
               h_ratio -> SetBinContent( bi_ratio_hist, ratio_val ) ;
               h_ratio -> SetBinError( bi_ratio_hist, ratio_err ) ;
            } else {
               h_ratio -> SetBinContent( bi_ratio_hist, -9 ) ;
               h_ratio -> SetBinError( bi_ratio_hist, 0. ) ;
            }
            h_ratio -> GetXaxis() -> SetBinLabel( bi_ratio_hist, label ) ;


         } // bi_nj
      } // bi_ht


      h_ratio -> GetXaxis() -> LabelsOption( "v" ) ;

      gStyle -> SetPadBottomMargin(0.20) ;
      gStyle -> SetOptStat(0) ;
      h_ratio -> SetMarkerStyle(20) ;
      h_ratio -> SetMinimum(-0.1) ;
      h_ratio -> Draw() ;
      gPad -> SetGridy(1) ;


      fclose( ofp_nbsum ) ;
      printf("\n\n Wrote %s\n\n", nbsum_text_file ) ;

      fclose( ofp_pars ) ;

      TFile rf( output_hist_file, "recreate" ) ;
      h_ratio -> Write() ;
      rf.Close() ;
      printf("\n\n Saved ratio histogram in %s\n\n", output_hist_file ) ;


   } // make_qcdmc_input_files1


#endif

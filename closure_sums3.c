#ifndef closure_sums3_c
#define closure_sums3_c

#include <stdio.h>
#include "TH1F.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include "binning.h"
#include "read_pars.h"
#include "histio.c"
#include <string.h>
#include "num_to_str.h"


   int sum_index_array[1000] ;
   int sum_nbins ;
   char sum_name[100] ;

   float qcd_ldp_val[1000] ;
   float qcd_ldp_err[1000] ;
   float qcd_hdp_val[1000] ;
   float qcd_hdp_err[1000] ;
   char  qcd_label[1000][100] ;

   float model_sum_val ;
   float model_sum_err ;
   float model_sum_err_fit ;
   float model_sum_err_syst ;
   float model_sum_err_mc ;
   float qcd_sum_val ;
   float qcd_sum_err ;

   bool skip_high_weight_search_bins ;

   bool include_qcdmc_correction ;

   TH1F* h_ratio_qcdmc_minus_model ;

   void htmht_bin_to_ht_and_mht_bins( int bi_htmht, int& bi_ht, int& bi_mht ) ;
   void read_qcd( const char* qcd_file ) ;
   void calc_sum( bool verb ) ;
   void print_search_bin_table() ;


  //----

   void closure_sums3(
          bool verb = true,
          bool arg_include_qcdmc_correction = false,
          const char* model_pars_file = "outputfiles/model-pars-qcdmc3.txt",
          const char* model_ratio_hist_file = "outputfiles/model-ratio-hist1.root",
          const char* qcd_file  = "outputfiles/qcdmc-counts.txt"
                     ) {

      include_qcdmc_correction = arg_include_qcdmc_correction ;
      setup_bins();
      gDirectory -> Delete( "h*" ) ;

      loadHist( model_ratio_hist_file ) ;
      h_ratio_qcdmc_minus_model = (TH1F*) gDirectory -> FindObject( "h_ratio_qcdmc_minus_model" ) ;
      if ( h_ratio_qcdmc_minus_model == 0x0 ) { printf("\n\n *** missing h_ratio_qcdmc_minus_model.\n\n" ) ; return ; }

      read_pars( model_pars_file ) ;
      read_qcd( qcd_file ) ;

      //skip_high_weight_search_bins = true ;
      skip_high_weight_search_bins = false ;

      print_search_bin_table() ;




      TH1F* h_closure_ht_model_fit   = new TH1F( "h_closure_ht_model_fit", "HT bins, model", nBinsHT, 0.5, nBinsHT + 0.5 ) ;
      TH1F* h_closure_ht_model_fit_and_syst   = new TH1F( "h_closure_ht_model_fit_and_syst", "HT bins, model", nBinsHT, 0.5, nBinsHT + 0.5 ) ;
      TH1F* h_closure_ht_model_total = new TH1F( "h_closure_ht_model_total", "HT bins, model", nBinsHT, 0.5, nBinsHT + 0.5 ) ;
      TH1F* h_closure_ht_qcdmc = new TH1F( "h_closure_ht_qcdmc", "HT bins, QCD MC", nBinsHT, 0.5, nBinsHT + 0.5 ) ;

      for ( int bin_ht = nBinsHT; bin_ht > 0; bin_ht--)
      {
         sprintf( sum_name, "Search HT%d",bin_ht ) ;
         sum_nbins = 0 ;
         int bi_hist(0);
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;
                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != bin_ht ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;

         calc_sum( verb ) ;


         int hbi(bin_ht) ;
         h_closure_ht_model_fit -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_ht_model_fit_and_syst -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_ht_model_total -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_ht_model_fit -> SetBinError( hbi, model_sum_err_fit ) ;
         h_closure_ht_model_fit_and_syst -> SetBinError( hbi, sqrt( pow(model_sum_err_fit,2.)+pow(model_sum_err_syst,2.)) ) ;
         h_closure_ht_model_total -> SetBinError( hbi, model_sum_err ) ;
         h_closure_ht_model_total -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;
         h_closure_ht_qcdmc -> SetBinContent( hbi, qcd_sum_val ) ;
         h_closure_ht_qcdmc -> SetBinError( hbi, qcd_sum_err ) ;
         h_closure_ht_qcdmc -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;



      }// bin_ht









     //---- 4 Njet bins
      TH1F* h_closure_njet_model_fit   = new TH1F( "h_closure_njet_model_fit", "njet bins, model", nb_nj, 0.5, nb_nj + 0.5 ) ;
      TH1F* h_closure_njet_model_fit_and_syst   = new TH1F( "h_closure_njet_model_fit_and_syst", "njet bins, model", nb_nj, 0.5, nb_nj + 0.5 ) ;
      TH1F* h_closure_njet_model_total = new TH1F( "h_closure_njet_model_total", "njet bins, model", nb_nj, 0.5, nb_nj + 0.5 ) ;
      TH1F* h_closure_njet_qcdmc = new TH1F( "h_closure_njet_qcdmc", "njet bins, QCD MC", nb_nj, 0.5, nb_nj + 0.5 ) ;


      for ( int bin_nj = 1; bin_nj <= nb_nj; bin_nj++)
      {
         sprintf( sum_name, "Search Nj%d",bin_nj ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;
                  
                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_mht == 1 ) continue ;
                  if ( bi_nj != bin_nj ) continue;
                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
	 } //bi_nj
         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(bin_nj) ;
         h_closure_njet_model_fit -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_njet_model_fit_and_syst -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_njet_model_total -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_njet_model_fit -> SetBinError( hbi, model_sum_err_fit ) ;
         h_closure_njet_model_fit_and_syst -> SetBinError( hbi, sqrt( pow(model_sum_err_fit,2.)+pow(model_sum_err_syst,2.)) ) ;
         h_closure_njet_model_total -> SetBinError( hbi, model_sum_err ) ;
         h_closure_njet_model_total -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;
         h_closure_njet_qcdmc -> SetBinContent( hbi, qcd_sum_val ) ;
         h_closure_njet_qcdmc -> SetBinError( hbi, qcd_sum_err ) ;
         h_closure_njet_qcdmc -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;

      }//bi_nj


      for ( int bin_nj = 1; bin_nj <= nb_nj; bin_nj++)
      {
         sprintf( sum_name, "Search HT3, Nj%d",bin_nj ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;
  
                  if ( bi_nj != bin_nj ) continue;
                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj
         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }//bi_nj





      TH1F* h_closure_nb_model_fit   = new TH1F( "h_closure_nb_model_fit", "nb bins, model", nb_nb, 0.5, nb_nb + 0.5 ) ;
      TH1F* h_closure_nb_model_fit_and_syst   = new TH1F( "h_closure_nb_model_fit_and_syst", "nb bins, model", nb_nb, 0.5, nb_nb + 0.5 ) ;
      TH1F* h_closure_nb_model_total = new TH1F( "h_closure_nb_model_total", "nb bins, model", nb_nb, 0.5, nb_nb + 0.5 ) ;
      TH1F* h_closure_nb_qcdmc = new TH1F( "h_closure_nb_qcdmc", "nb bins, QCD MC", nb_nb, 0.5, nb_nb + 0.5 ) ;

      for ( int bin_nb = 1; bin_nb <= nb_nb; bin_nb++)
      {
         sprintf( sum_name, "Search Nb%d", bin_nb-1 ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;
   
                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;
   
                  if ( bi_nb != bin_nb ) continue;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            }//bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(bin_nb) ;
         h_closure_nb_model_fit -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_nb_model_fit_and_syst -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_nb_model_total -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_nb_model_fit -> SetBinError( hbi, model_sum_err_fit ) ;
         h_closure_nb_model_fit_and_syst -> SetBinError( hbi, sqrt( pow(model_sum_err_fit,2.)+pow(model_sum_err_syst,2.)) ) ;
         h_closure_nb_model_total -> SetBinError( hbi, model_sum_err ) ;
         h_closure_nb_model_total -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;
         h_closure_nb_qcdmc -> SetBinContent( hbi, qcd_sum_val ) ;
         h_closure_nb_qcdmc -> SetBinError( hbi, qcd_sum_err ) ;
         h_closure_nb_qcdmc -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;

      } // bin_nb




      for ( int bin_nb = 1; bin_nb <= nb_nb; bin_nb++)
      {
         sprintf( sum_name, "Search HT3, Nb%d", bin_nb-1) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nb != bin_nb ) continue;
                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } //bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }//bin_nb


     //--- 4 MHT bins

      TH1F* h_closure_mht_model_fit   = new TH1F( "h_closure_mht_model_fit", "mht bins, model", nb_mht-1, 0.5, nb_mht - 0.5 ) ;
      TH1F* h_closure_mht_model_fit_and_syst   = new TH1F( "h_closure_mht_model_fit_and_syst", "mht bins, model", nb_mht-1, 0.5, nb_mht - 0.5 ) ;
      TH1F* h_closure_mht_model_total = new TH1F( "h_closure_mht_model_total", "mht bins, model", nb_mht-1, 0.5, nb_mht - 0.5 ) ;
      TH1F* h_closure_mht_qcdmc = new TH1F( "h_closure_mht_qcdmc", "mht bins, QCD MC", nb_mht-1, 0.5, nb_mht - 0.5 ) ;

      for ( int bin_mht = 1; bin_mht < nb_mht; bin_mht++)
      {


         sprintf( sum_name, "Search MHT%d",bin_mht);
 

         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_mht != bin_mht+1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(bin_mht) ;
         h_closure_mht_model_fit -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_mht_model_fit_and_syst -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_mht_model_total -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_mht_model_fit -> SetBinError( hbi, model_sum_err_fit ) ;
         h_closure_mht_model_fit_and_syst -> SetBinError( hbi, sqrt( pow(model_sum_err_fit,2.)+pow(model_sum_err_syst,2.)) ) ;
         h_closure_mht_model_total -> SetBinError( hbi, model_sum_err ) ;
         h_closure_mht_model_total -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;
         h_closure_mht_qcdmc -> SetBinContent( hbi, qcd_sum_val ) ;
         h_closure_mht_qcdmc -> SetBinError( hbi, qcd_sum_err ) ;
         h_closure_mht_qcdmc -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;


      }





      for ( int bin_mht = 1; bin_mht < nb_mht; bin_mht++)
      {
         sprintf( sum_name, "Search HT3, MHT%d",bin_mht ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht != bin_mht+1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }// bin_mht

      for ( int bin_nj=1; bin_nj<=nb_nj; bin_nj++ ) 
      {
         sprintf( sum_name, "Search HT2, Nj%d", bin_nj ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 2 ) continue ;
                  if ( bi_nj != bin_nj ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }//bin_nj

      for ( int bin_nj=1; bin_nj<=nb_nj; bin_nj++ )
      {
         if ( bin_nj > nb_nb - 2 ) continue;
         sprintf( sum_name, "Search HT1, Nj%d", bin_nj ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 1 ) continue ;
                  if ( bi_nj != bin_nj ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }//bin_nj




     //--- boxes
      TString htmht_str = num_to_str(nb_htmht-nBinsHT);
      TH1F* h_closure_10boxes_model_fit   = new TH1F( "h_closure_10boxes_model_fit", htmht_str+" MHT-HT boxes, model", nb_htmht-nBinsHT, 0.5, nb_htmht-nBinsHT+0.5 ) ;
      TH1F* h_closure_10boxes_model_fit_and_syst   = new TH1F( "h_closure_10boxes_model_fit_and_syst", htmht_str+
		                                               " MHT-HT boxes, model", nb_htmht-nBinsHT, 0.5, nb_htmht-nBinsHT+0.5 ) ;
      TH1F* h_closure_10boxes_model_total = new TH1F( "h_closure_10boxes_model_total", htmht_str+" MHT-HT boxes, model", nb_htmht-nBinsHT, 0.5, nb_htmht-nBinsHT+0.5 ) ;
      TH1F* h_closure_10boxes_qcdmc = new TH1F( "h_closure_10boxes_qcdmc", htmht_str+" MHT-HT boxes, QCD MC", nb_htmht-nBinsHT, 0.5, nb_htmht-nBinsHT+0.5 ) ;

      for ( int bin_htmht=3+1; bin_htmht<=nb_htmht; bin_htmht++ )
      {
         sprintf( sum_name, "Search box%d", bin_htmht-3 ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != bin_htmht ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(bin_htmht-3) ;
         h_closure_10boxes_model_fit -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_10boxes_model_fit_and_syst -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_10boxes_model_total -> SetBinContent( hbi, model_sum_val ) ;
         h_closure_10boxes_model_fit -> SetBinError( hbi, model_sum_err_fit ) ;
         h_closure_10boxes_model_fit_and_syst -> SetBinError( hbi, sqrt( pow(model_sum_err_fit,2.)+pow(model_sum_err_syst,2.)) ) ;
         h_closure_10boxes_model_total -> SetBinError( hbi, model_sum_err ) ;
         h_closure_10boxes_model_total -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;
         h_closure_10boxes_qcdmc -> SetBinContent( hbi, qcd_sum_val ) ;
         h_closure_10boxes_qcdmc -> SetBinError( hbi, qcd_sum_err ) ;
         h_closure_10boxes_qcdmc -> GetXaxis() -> SetBinLabel( hbi, sum_name ) ;

      }


      if ( include_qcdmc_correction ) {
         saveHist( "outputfiles/closure-sums3-wcor.root", "h*" ) ;
      } else {
         saveHist( "outputfiles/closure-sums3.root", "h*" ) ;
      }


   } // closure_sums3

  //=======================================================================================

   void read_qcd( const char* qcd_file ) {

      ifstream ifs_qcd_ldp ;
      ifs_qcd_ldp.open( qcd_file ) ;
      if ( !ifs_qcd_ldp.good() ) { printf("\n\n *** Problem opening %s\n\n", qcd_file ) ; return ; }


      int bi_hist(0) ;
      int bi_control(0) ;
      int bi_search(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               TString line ;
               line.ReadLine( ifs_qcd_ldp ) ;
               if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

               bi_hist ++ ;
               if ( bi_htmht <= 3 ) {
                  bi_control ++ ;
               } else {
                  bi_search ++ ;
               }

               int bi_ht, bi_mht ;

               htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;


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

               int rbi1, rbi2 ;
               char rsc[10] ;
               char rlabel[100] ;
               float ldp_val, ldp_err, hdp_val, hdp_err ;
               sscanf( line.Data(), " %d %s %d %s  %f %f %f %f", &rbi1, rsc, &rbi2, rlabel, &ldp_val, &ldp_err, &hdp_val, &hdp_err ) ;
//               if ( rbi1 != bi_hist ) { printf("\n\n *** inconsistency!  %d : %s\n\n", bi_hist, line.Data() ) ; return ; }

               qcd_ldp_val[bi_hist] = ldp_val ;
               qcd_ldp_err[bi_hist] = ldp_err ;
               qcd_hdp_val[bi_hist] = hdp_val ;
               qcd_hdp_err[bi_hist] = hdp_err ;
               sprintf( qcd_label[bi_hist], "%s", label ) ;
               printf( " %30s   HT ind %d,  Nj ind %d,  MHT ind %d, Nb ind %d    %8.1f  %8.1f   %8.1f  %8.1f\n", label,
                  bi_ht, bi_nj, bi_mht, bi_nb,
                  qcd_ldp_val[bi_hist], qcd_ldp_err[bi_hist],
                  qcd_hdp_val[bi_hist], qcd_hdp_err[bi_hist]
                  ) ;

            } // bi_htmht
         } // bi_nb
      } // bi_nj

   } // read_qcd

  //=======================================================================================

   void calc_sum( bool verb ) {

         model_sum_val = 0. ;
         qcd_sum_val = 0. ;
         float model_sum_correction_val(0.) ;
         float model_sum_correction_err2(0.) ;
         float qcd_sum_err2 = 0. ;
         bool skip_this_search_bin[1000] ;

         for ( int i=1; i<=nb_global_w_exclusion_w_mhtc; i++ ) { skip_this_search_bin[i] = false ; }

         for ( int sbi=0; sbi<sum_nbins; sbi++ ) {

            int hbi = sum_index_array[sbi] ;
            if ( skip_high_weight_search_bins && qcd_hdp_err[hbi] > 10. ) {
               skip_this_search_bin[hbi] = true ;
               printf("  Skipping %3d  %30s  %8.1f +/- %8.1f\n", hbi, qcd_label[hbi], qcd_hdp_val[hbi], qcd_hdp_err[hbi] ) ;
               continue ;
            }



           //--- skip HT1 at high Njet.
            int bi_nj, bi_nb, bi_ht, bi_mht;
            translate_qcd_bin_to_nj_nb_ht_mht (hbi, bi_nj, bi_nb, bi_ht, bi_mht);
            if ( bi_ht == 1 && bi_nj > nb_nj-2 ) continue ;

            float model_ratio(0.) ;
            model_ratio =  par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_ht_mht[bi_ht][bi_mht] * par_val_nb[bi_nb] ;
            float ratio_correction_val(0.) ;
            float ratio_correction_err(0.) ;
            {
               ratio_correction_val = h_ratio_qcdmc_minus_model -> GetBinContent( translate_qcd_bin_to_search_bin(hbi) ) ;
               ratio_correction_err = h_ratio_qcdmc_minus_model -> GetBinError( translate_qcd_bin_to_search_bin(hbi) ) ;
               if ( verb ) {
                  printf("  Correcting QCD bin %3d (%30s) with search bin %d (%30s) : correction = %6.4f +/- %6.4f\n",
                      hbi, qcd_label[hbi], translate_qcd_bin_to_search_bin(hbi), h_ratio_qcdmc_minus_model -> GetXaxis() -> GetBinLabel( translate_qcd_bin_to_search_bin(hbi) ),
                      ratio_correction_val, ratio_correction_err ) ;
               }
            }

            float model_val = qcd_ldp_val[hbi] * model_ratio ;
            float model_correction_val = qcd_ldp_val[hbi] * ratio_correction_val ;

            model_sum_val += model_val ;
            model_sum_correction_val += model_correction_val ;
            model_sum_correction_err2 += pow( qcd_ldp_val[hbi] * ratio_correction_err, 2. ) ;
            qcd_sum_val += qcd_hdp_val[hbi] ;
            qcd_sum_err2 += pow( qcd_hdp_err[hbi], 2. ) ;

            if (verb) printf("  %3d  %30s   model = %8.1f  (correction %8.1f) ,  QCD = %8.1f +/- %8.1f\n",
              hbi, qcd_label[hbi], model_val, model_correction_val, qcd_hdp_val[hbi], qcd_hdp_err[hbi] ) ;

         } // sbi
         float model_sum_correction_err = sqrt( model_sum_correction_err2 ) ;

         qcd_sum_err = sqrt( qcd_sum_err2 ) ;



         float partial_sum[10][100]{} ;

         for ( int sbi=0; sbi<sum_nbins; sbi++ ) {

            int hbi = sum_index_array[sbi] ;

            int bi_nj, bi_nb, bi_ht, bi_mht;
            translate_qcd_bin_to_nj_nb_ht_mht (hbi, bi_nj, bi_nb, bi_ht, bi_mht);

            if ( skip_this_search_bin[hbi] ) continue ;

           //--- skip HT1 at high Njet.
            if ( bi_ht == 1 && bi_nj > nb_nj-2 ) continue ;

            float model_ratio(0.) ;

	    model_ratio =  par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_ht_mht[bi_ht][bi_mht] * par_val_nb[bi_nb] ;
	    float model_val = qcd_ldp_val[hbi] * model_ratio ;

          //---


            if ( par_val_ht  [bi_ht] > 0. )
               partial_sum[1][bi_ht] += model_val / par_val_ht[bi_ht] ;

            if ( par_val_njet[bi_nj] > 0. )
               partial_sum[2][bi_nj] += model_val / par_val_njet[bi_nj] ;
	    for ( int bin_htmht = 1; bin_htmht <= nb_htmht; bin_htmht++)
	    {
               int bin_ht = 0, bin_mht = 0;
               htmht_bin_to_ht_and_mht_bins(bin_htmht, bin_ht,bin_mht);
               if ( bin_ht == bi_ht && bin_mht == bi_mht && par_val_ht_mht[bi_ht][bi_mht] > 0. )
                  partial_sum[3][bin_htmht] += model_val / par_val_ht_mht[bi_ht][bi_mht] ;

	    }

            if ( par_val_nb[bi_nb] > 0. )
               partial_sum[4][bi_nb] += model_val / par_val_nb[bi_nb] ;

        }//sbi
         float par_err_fit[10][100] ;
         float par_err_syst[10][100] ;
         char  par_name[10][36][100] ; 

         int npar(0) ;
         { 
	 for ( int bin_ht = 1; bin_ht <= nBinsHT; bin_ht++)
	 {
            par_err_fit[1][bin_ht] = par_err_ht_fit[bin_ht] ; par_err_syst[1][bin_ht] = par_err_ht_syst[bin_ht] ; sprintf( par_name[1][bin_ht], "HT%d", bin_ht ) ;

         }//bin_ht

         for ( int bin_nj = 1; bin_nj <= nb_nj; bin_nj++)
	 {
            par_err_fit[2][bin_nj] = par_err_njet_fit[bin_nj] ; par_err_syst[2][bin_nj] = par_err_njet_syst[bin_nj] ; sprintf( par_name[2][bin_nj], "Njet%d", bin_nj ) ;
	 }//bin_nj

	 for ( int bin_htmht = 1; bin_htmht <= nb_htmht; bin_htmht++)
         {

	    int bin_ht, bin_mht;
            htmht_bin_to_ht_and_mht_bins (bin_htmht, bin_ht, bin_mht);

	    char mht_str[10] ; sprintf(mht_str, "%d", bin_mht - 1);

	    if (  bin_mht == 1 ) sprintf(mht_str,"C");

            if ( bin_ht == 1 ) sprintf( par_name[3][bin_htmht], "MHT%s-HTL",mht_str);
            if ( bin_ht == 2 ) sprintf( par_name[3][bin_htmht], "MHT%s-HTM",mht_str);
            if ( bin_ht == 3 ) sprintf( par_name[3][bin_htmht], "MHT%s-HTH",mht_str);

            par_err_fit[3][bin_htmht] = 0. ; par_err_syst[3][bin_htmht] = par_err_ht_mht[bin_ht][bin_mht] ;


	 }//bin_htmht



	 for ( int bin_nb = 1; bin_nb <= nb_nb; bin_nb++)
	 {
            par_err_fit[4][bin_nb] = 0. ; par_err_syst[4][bin_nb] = par_err_nb[bin_nb] ; sprintf( par_name[4][bin_nb], "Nb%d", bin_nb-1 ) ;
	 }//bin_nb

         }

         float model_sum_err2_fit(0.), model_sum_err2_syst(0.) ;

         for ( int bg = 1; bg < 5; bg++)
	 {
            char bg_name[10];
	    int nb_bg = 0;
            if ( bg == 1 ) {nb_bg = nBinsHT ; strcpy(bg_name,"HT");}
            if ( bg == 2 ) {nb_bg = nb_nj   ; strcpy(bg_name,"NJets");}
            if ( bg == 3 ) {nb_bg = nb_htmht; strcpy(bg_name,"HT_MHT");}
            if ( bg == 4 ) {nb_bg = nb_nb; strcpy(bg_name,"NbJets");}

	 for ( int bin_bg = 1; bin_bg <= nb_bg; bin_bg++)
	 {

            float par_err = sqrt( pow( par_err_fit[bg][bin_bg], 2.) + pow( par_err_syst[bg][bin_bg], 2. ) ) ;
            if (verb) {


            printf("  %5s : partial sum %8.1f, err %8.4f (%8.4f (fit) +/- %8.4f (syst)), contrib %8.1f\n",
                  par_name[bg][bin_bg], partial_sum[bg][bin_bg], par_err, par_err_fit[bg][bin_bg], par_err_syst[bg][bin_bg],
                  partial_sum[bg][bin_bg] * par_err ) ;
            } //verb
            model_sum_err2_fit += pow( partial_sum[bg][bin_bg] * par_err_fit[bg][bin_bg], 2. ) ;
            model_sum_err2_syst += pow( partial_sum[bg][bin_bg] * par_err_syst[bg][bin_bg], 2. ) ;

	 }//bin_bg
         }//bg



         model_sum_err_fit = sqrt( model_sum_err2_fit ) ;
         model_sum_err_mc = sqrt( model_sum_correction_err2 ) ;
         model_sum_err_syst = sqrt( model_sum_err2_syst ) ;

         model_sum_err = sqrt( model_sum_err2_fit + model_sum_correction_err2 + model_sum_err2_syst ) ;

         if ( include_qcdmc_correction ) {
            model_sum_val += model_sum_correction_val ;
         }

         printf("  Sums:  model = %8.1f +/- %8.1f (%8.5f (fit), %8.5f (syst))     QCD = %8.1f +/- %8.1f\n",
            model_sum_val, model_sum_err, model_sum_err_fit, model_sum_err_syst,
            qcd_sum_val, qcd_sum_err ) ;

     } // calc_sum

  //=======================================================================================

   void print_search_bin_table() {

      printf("\n\n") ;

      TH1F* h_ratio_all = (TH1F*) gDirectory -> FindObject( "h_ratio_all" ) ;
      if ( h_ratio_all == 0x0 ) { printf("\n\n *** Missing h_ratio_all hist.\n\n") ; gSystem->Exit(-1) ; }

      int hbi = 0;
      for ( int i = 1; i <= nb_global_after_exclusion; i++ )
      {

	 hbi++;
         int bi_nj, bi_nb, bi_ht, bi_mht;
         translate_search_bin_to_nj_nb_ht_mht (hbi, bi_nj, bi_nb, bi_ht, bi_mht);

         float model_ratio(0.) ;
         model_ratio =  par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_ht_mht[bi_ht][bi_mht] * par_val_nb[bi_nb] ;

         float model_ratio_from_hist = h_ratio_all -> GetBinContent( hbi ) ;//amin

         if ( bi_ht == 1 && bi_nj > nb_nj-2 ) continue ;

         if ( fabs(model_ratio-model_ratio_from_hist) > 0.00001 ) { printf("\n *** inconsistent model ratios in Njet bin = %d, Nbjet = bin %d, HT bin = %d and MHT bin = %d.\n\n",
			bi_nj, bi_nb, bi_ht, bi_mht ) ; }

         float ratio_correction_val(0.) ;
         float ratio_correction_err(0.) ;
         {
            ratio_correction_val = h_ratio_qcdmc_minus_model -> GetBinContent( hbi ) ;
            ratio_correction_err = h_ratio_qcdmc_minus_model -> GetBinError( hbi ) ;
         }
         float model_val = qcd_ldp_val[hbi] * model_ratio ;

         float model_correction_val = qcd_ldp_val[hbi] * ratio_correction_val ;
         float model_correction_err = qcd_ldp_val[hbi] * ratio_correction_err ;
	 printf("  %3d   %30s :  Model %6.1f +/- %6.1f,  cor = %9.3f * %7.4f = %6.1f +/- %6.1f,  total = %6.1f +/- %6.1f,   HDP QCD MC  %6.1f +/- %6.1f  ",
             hbi, qcd_label[hbi],
             model_val, 0.,
             qcd_ldp_val[hbi], ratio_correction_val, model_correction_val, model_correction_err,
             model_val+model_correction_val, 0.,
             qcd_hdp_val[hbi], qcd_hdp_err[hbi] ) ;

         if ( fabs( model_val+model_correction_val - qcd_hdp_val[hbi] ) > 0.01 && include_qcdmc_correction ) {
            printf(" ***\n") ;
         } else {
            printf("\n") ;
         }

      } // hbi

      printf("\n\n") ;

   } // print_search_bin_table


#endif

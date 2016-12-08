


#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "../binning.h"
#include "../get_hist.h"
#include "../histio.c"

#include <fstream>

   int sum_index_array[1000] ;
   int sum_nbins ;

   float binsum_val(0.) ;
   float binsum_err(0.) ;
   float binsum_err_nocorrelations(0.) ;

   int n_syst_cols(23) ;

   float syst_rel_err[100] ;
   float mcc_rel_err ;
   float nqcd_val ;
   float nldp_err ;


   ifstream ifs_combine_table ;
   FILE* ofp_text ;

  //------------

   void read_line( int bin_number ) ;
   void calc_sum( bool verb ) ;

  //------------

   void binsums_qcd1(   const char* input_combine_file = "outputfiles/combine-input-all.txt",
                        const char* output_text_file = "outputfiles/binsums-qcd.txt",
                        const char* output_root_file = "outputfiles/binsums-qcd.root"
                    ) {

      setup_bins();

      //bool verb(false) ;
      bool verb(true) ;

      gDirectory -> Delete( "h*" ) ;

      ifs_combine_table.open( input_combine_file ) ;
      if ( !ifs_combine_table.good() ) { printf("\n\n *** Bad input combine file: %s\n\n", input_combine_file ) ; gSystem -> Exit(-1) ; }

      if ( (ofp_text = fopen( output_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening binsums output file: %s\n\n", output_text_file ) ;
         return ;
      }

      int all_bins_summed ;
      TH1F* hp ;


    //============= Njet bins

      hp = new TH1F( "h_njet_proj", "Njets bins", nb_nj, 0.5, nb_nj+0.5 ) ;

      all_bins_summed = 0 ;

      for ( int bin_nj = 1; bin_nj <= nb_nj; bin_nj++) {

         sum_nbins = 0 ;
         int bi_hist(0) ;

         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_nj != bin_nj ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } //bi_nj

         all_bins_summed += sum_nbins ;

         printf("\n\n Calculating search Nj%d\n", bin_nj ) ;
         calc_sum( verb ) ;

         char sum_label[100] ;
         sprintf( sum_label, "bin_Nj%d_proj", bin_nj ) ;

         fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
         for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
         fprintf( ofp_text, "\n" ) ;

         printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

         int hb = bin_nj ;
         hp -> SetBinContent( hb, binsum_val ) ;
         hp -> SetBinError( hb, binsum_err ) ;
         hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

      }//bin_nj

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, Njet: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }





    //============= HT bins

      hp = new TH1F( "h_ht_proj", "HT bins", nBinsHT, 0.5, nBinsHT+0.5 ) ;

      all_bins_summed = 0 ;

      for ( int bin_ht = 1; bin_ht <= nBinsHT; bin_ht++) {

         sum_nbins = 0 ;
         int bi_hist(0) ;

         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_ht != bin_ht ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } //bi_nj

         all_bins_summed += sum_nbins ;

         printf("\n\n Calculating search HT%d\n", bin_ht ) ;
         calc_sum( verb ) ;

         char sum_label[100] ;
         sprintf( sum_label, "bin_HT%d_proj", bin_ht ) ;

         fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
         for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
         fprintf( ofp_text, "\n" ) ;

         printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

         int hb = bin_ht ;
         hp -> SetBinContent( hb, binsum_val ) ;
         hp -> SetBinError( hb, binsum_err ) ;
         hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

      }//bin_ht

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, HT: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }






    //============= Nb bins

      hp = new TH1F( "h_nb_proj", "Nb bins", nb_nb, 0.5, nb_nb+0.5 ) ;

      all_bins_summed = 0 ;

      for ( int bin_nb = 1; bin_nb <= nb_nb; bin_nb++) {

         sum_nbins = 0 ;
         int bi_hist(0) ;

         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_nb != bin_nb ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } //bi_nj

         all_bins_summed += sum_nbins ;

         printf("\n\n Calculating search Nb%d\n", bin_nb-1 ) ;
         calc_sum( verb ) ;

         char sum_label[100] ;
         sprintf( sum_label, "bin_Nb%d_proj", bin_nb-1 ) ;

         fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
         for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
         fprintf( ofp_text, "\n" ) ;

         printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

         int hb = bin_nb ;
         hp -> SetBinContent( hb, binsum_val ) ;
         hp -> SetBinError( hb, binsum_err ) ;
         hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

      }//bin_nb

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, Nb: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }






    //============= MHT bins

      hp = new TH1F( "h_mht_proj", "MHT bins", nb_mht-1, 0.5, nb_mht-1+0.5 ) ;

      all_bins_summed = 0 ;

      for ( int bin_mb = 2; bin_mb <= nb_mht; bin_mb++) {

         sum_nbins = 0 ;
         int bi_hist(0) ;

         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_mht != bin_mb ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } //bi_nj

         all_bins_summed += sum_nbins ;

         printf("\n\n Calculating search MHT%d\n", bin_mb-1 ) ;
         calc_sum( verb ) ;

         char sum_label[100] ;
         sprintf( sum_label, "bin_MHT%d_proj", bin_mb-1 ) ;

         fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
         for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
         fprintf( ofp_text, "\n" ) ;

         printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

         int hb = bin_mb-1 ;
         hp -> SetBinContent( hb, binsum_val ) ;
         hp -> SetBinError( hb, binsum_err ) ;
         hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

      }//bin_mb

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, MHT: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }





    //============= HT-MHT bins

      hp = new TH1F( "h_htmhtboxes_proj", "10 HT-MHT boxes", nb_htmht-3, 0.5, nb_htmht-3+0.5 ) ;

      all_bins_summed = 0 ;

      for ( int bin_htmht = 4; bin_htmht <= nb_htmht; bin_htmht++) {

         sum_nbins = 0 ;
         int bi_hist(0) ;

         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_htmht != bin_htmht ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } //bi_nj

         all_bins_summed += sum_nbins ;

         printf("\n\n Calculating search HT-MHT box %d\n", bin_htmht-3 ) ;
         calc_sum( verb ) ;

         char sum_label[100] ;
         sprintf( sum_label, "bin_HTMHTbox%02d_proj", bin_htmht-3 ) ;

         fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
         for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
         fprintf( ofp_text, "\n" ) ;

         printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

         int hb = bin_htmht-3 ;
         hp -> SetBinContent( hb, binsum_val ) ;
         hp -> SetBinError( hb, binsum_err ) ;
         hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

      }//bin_htmht

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, HT-MHT boxes: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }






    //============= All Njet + HT-MHT bin combinations


      all_bins_summed = 0 ;


      for ( int bin_htmht = 4; bin_htmht <= nb_htmht; bin_htmht++) {

         char hname[100] ;
         char htitle[100] ;
         sprintf( hname, "h_njhtmht%02d_proj", bin_htmht-3 ) ;
         sprintf( htitle, "Njets, HTMHT box %d", bin_htmht-3 ) ;
         hp = new TH1F( hname, htitle, nb_nj, 0.5, nb_nj+0.5 ) ;

         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            sum_nbins = 0 ;

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_htmht != bin_htmht ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb

            all_bins_summed += sum_nbins ;

            printf("\n\n Calculating search HT-MHT box %d\n", bin_htmht-3 ) ;
            calc_sum( verb ) ;

            char sum_label[100] ;
            sprintf( sum_label, "bin_HTMHTbox%02dnj%d_proj", bin_htmht-3, bi_nj ) ;

            if ( sum_nbins == 0 ) continue ;

            fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
            for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
            fprintf( ofp_text, "\n" ) ;

            printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

            int hb = bi_nj ;
            hp -> SetBinContent( hb, binsum_val ) ;
            hp -> SetBinError( hb, binsum_err ) ;
            hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

         } //bi_nj

      }//bin_htmht

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, Njet + HT-MHT box combinations: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }



    //============= All Nb + HT-MHT bin combinations


      all_bins_summed = 0 ;


      for ( int bin_htmht = 4; bin_htmht <= nb_htmht; bin_htmht++) {

         char hname[100] ;
         char htitle[100] ;
         sprintf( hname, "h_nbhtmht%02d_proj", bin_htmht-3 ) ;
         sprintf( htitle, "Nb, HTMHT box %d", bin_htmht-3 ) ;
         hp = new TH1F( hname, htitle, nb_nb, 0.5, nb_nb+0.5 ) ;


         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            sum_nbins = 0 ;

            for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
               for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  if (is_this_bin_excluded(bi_nj-1,bi_nb-1,bi_htmht-1)) continue;

                  int bi_ht, bi_mht ;

                  htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

                  int search_bi = global_search_bin( bi_nj, bi_nb, bi_htmht ) ;
                  if ( search_bi <= 0 ) { printf( "\n\n *** Illegal search bin index: bi_nj=%d, bi_nb=%d, bi_htmht=%d\n\n", bi_nj, bi_nb, bi_htmht ) ; gSystem->Exit(-1) ; }

                  if ( bi_mht == 1 ) continue ; // don't include MHTC (bi_mht==1)
                  if ( bi_htmht != bin_htmht ) continue;
                  sum_index_array[sum_nbins] = search_bi ;
                  sum_nbins++ ;

               } // bi_htmht
            } //bi_nj

            all_bins_summed += sum_nbins ;

            printf("\n\n Calculating search HT-MHT box %d\n", bin_htmht-3 ) ;
            calc_sum( verb ) ;

            char sum_label[100] ;
            sprintf( sum_label, "bin_HTMHTbox%02dnb%d_proj", bin_htmht-3, bi_nb-1 ) ;

            if ( sum_nbins == 0 ) continue ;

            fprintf( ofp_text, " %20s  %8.1f +/- %8.1f  %3d ", sum_label, binsum_val, binsum_err, sum_nbins ) ;
            for ( int i=0; i<sum_nbins; i++ ) { fprintf( ofp_text, " %3d ", sum_index_array[i] ) ; }
            fprintf( ofp_text, "\n" ) ;

            printf( "%20s : %8.1f +/- %8.1f (without correlations %8.1f)\n", sum_label, binsum_val, binsum_err, binsum_err_nocorrelations ) ;

            int hb = bi_nb ;
            hp -> SetBinContent( hb, binsum_val ) ;
            hp -> SetBinError( hb, binsum_err ) ;
            hp -> GetXaxis() -> SetBinLabel( hb, sum_label ) ;

         } // bi_nb

      }//bin_htmht

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

      printf(" All bins summed, Nb + HT-MHT box combinations: %3d\n", all_bins_summed ) ;
      if ( all_bins_summed != 174 ) { printf("\n\n *** Wrong total.\n\n") ; gSystem -> Exit(-1) ; }








      printf("\n\n Closing output text file: %s\n\n", output_text_file ) ;
      fclose( ofp_text ) ;


      saveHist( output_root_file, "h*" ) ;


   } // binsums_qcd1

  //===============================================================================================================

   void read_line( int bin_number ) {

      ifs_combine_table.seekg(0) ;

      TString line ;
      line.ReadLine( ifs_combine_table ) ; // column headers line

      while ( ifs_combine_table.good() ) {

         TString token ;

       //----
         token.ReadToken( ifs_combine_table ) ;
         if ( !ifs_combine_table.good() ) break ;
         int this_bin ;
         sscanf( token.Data(), "%d", &this_bin ) ;
         if ( this_bin <= 0 ) { printf("\n\n *** illegal bin number : %d %s\n\n", this_bin, token.Data() ) ; gSystem -> Exit(-1) ; }

       //----
         token.ReadToken( ifs_combine_table ) ;
         char this_bin_label[100] ;
         sprintf( this_bin_label, "%s", token.Data() ) ;

       //----
         token.ReadToken( ifs_combine_table ) ;
         int this_n_ldp ;
         sscanf( token.Data(), "%d", &this_n_ldp ) ;

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_nonqcd_val ;
         sscanf( token.Data(), "%f", &this_nonqcd_val ) ;

       //----
         token.ReadToken( ifs_combine_table ) ; // +/-

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_nonqcd_err ;
         sscanf( token.Data(), "%f", &this_nonqcd_err ) ;

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_rqcd_val ;
         sscanf( token.Data(), "%f", &this_rqcd_val ) ;

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_nqcd_val ;
         sscanf( token.Data(), "%f", &this_nqcd_val ) ;

       //----
         float this_syst_rel_err[100] ;
         for ( int si=0; si<n_syst_cols; si++ ) {
            token.ReadToken( ifs_combine_table ) ;
            this_syst_rel_err[si] = 0. ;
            if ( token.Contains("-") ) continue ;
            float par ;
            sscanf( token.Data(), "%f", &par ) ;
            this_syst_rel_err[si] = par - 1. ;
            if ( this_syst_rel_err[si] < 0 ) { printf("\n\n *** Illegal systematics parameter: %s\n\n", token.Data() ) ; gSystem -> Exit(-1) ; }
         } // si

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_mcc_rel_err(0.) ;
         if ( !(token.Contains("-")) ) {
            float par ;
            sscanf( token.Data(), "%f", &par ) ;
            this_mcc_rel_err = par - 1. ;
            if ( this_mcc_rel_err < 0 ) { printf("\n\n *** Illegal systematics parameter: %s\n\n", token.Data() ) ; gSystem -> Exit(-1) ; }
         }

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_rqcd_val2 ;
         sscanf( token.Data(), "%f", &this_rqcd_val2 ) ;
         if ( this_rqcd_val2 != this_rqcd_val ) { printf("\n\n *** Read inconsistent Rqcd vals.  bin %d  %.4f != %.4f\n\n", this_bin, this_rqcd_val, this_rqcd_val2 ) ;  gSystem -> Exit(-1) ; }

       //----
         token.ReadToken( ifs_combine_table ) ; // +/-

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_rqcd_err ;
         sscanf( token.Data(), "%f", &this_rqcd_err ) ;

       //----
         token.ReadToken( ifs_combine_table ) ; // ,

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_nqcd_val2 ;
         sscanf( token.Data(), "%f", &this_nqcd_val2 ) ;
         if ( this_nqcd_val2 != this_nqcd_val ) { printf("\n\n *** Read inconsistent nqcd vals.  bin %d  %.4f != %.4f\n\n", this_bin, this_nqcd_val, this_nqcd_val2 ) ;  gSystem -> Exit(-1) ; }

       //----
         token.ReadToken( ifs_combine_table ) ; // +/-

       //----
         token.ReadToken( ifs_combine_table ) ;
         float this_nqcd_err_from_file ;
         sscanf( token.Data(), "%f", &this_nqcd_err_from_file ) ;

       //----
         token.ReadToken( ifs_combine_table ) ; // (
       //----
         token.ReadToken( ifs_combine_table ) ; // rel error in %
       //----
         token.ReadToken( ifs_combine_table ) ; // %)

         ///////printf("  %3d : %35s  %8d\n", this_bin, this_bin_label, this_n_ldp ) ;

         if ( this_bin == bin_number ) {
            for ( int si=0; si<n_syst_cols; si++ ) { syst_rel_err[si] = this_syst_rel_err[si] ; }
            nqcd_val = this_nqcd_val ;
            if ( this_n_ldp == 0 ) this_n_ldp = 1. ;
            nldp_err = this_rqcd_val * sqrt( this_n_ldp + pow( this_nonqcd_err, 2. ) ) ;
            mcc_rel_err = this_mcc_rel_err ;
            //////printf("  %4d : %s : %8.1f\n", this_bin, this_bin_label, nqcd_val ) ;
            return ;
         }


      }

      printf("\n\n *** Did not find bin number %d\n\n", bin_number ) ;
      gSystem -> Exit(-1) ;


   } // read_line

  //===============================================================================================================


   void calc_sum( bool verb ) {

      binsum_val = 0. ;
      binsum_err = 0. ;

      float syst_err_sum[100] ;
      for ( int i=0; i<100; i++ ) { syst_err_sum[i] = 0. ; }

      float sum_nldp_err2(0.) ;
      float sum_mcc_err2(0.) ;

      float sum_err2_nocorrelations(0.) ;

      for ( int sbi=0; sbi<sum_nbins; sbi++ ) {

         read_line( sum_index_array[sbi] ) ;

         binsum_val += nqcd_val ;

         float this_bin_syst_err2(0.) ;
         for ( int si=0; si<n_syst_cols; si++ ) {
            syst_err_sum[si] += nqcd_val * syst_rel_err[si] ;
            this_bin_syst_err2 += pow( nqcd_val * syst_rel_err[si], 2. ) ;
         }
         sum_nldp_err2 += pow( nldp_err, 2. ) ;
         sum_mcc_err2 += pow( nqcd_val * mcc_rel_err, 2. ) ;

         float this_bin_err = sqrt( this_bin_syst_err2 + pow( nldp_err, 2. ) + pow( nqcd_val * mcc_rel_err, 2. ) ) ;

         sum_err2_nocorrelations += pow( this_bin_err, 2. ) ;

         if ( verb ) printf("  calc_sum %3d :  %9.2f +/- %6.2f  (%5.1f, %5.1f, %5.1f)\n",
            sum_index_array[sbi], nqcd_val, this_bin_err, sqrt( this_bin_syst_err2 ), nldp_err, mcc_rel_err ) ;

      } // sbi

      float total_syst_err2(0.) ;
      for ( int si=0; si<n_syst_cols; si++ ) {
         total_syst_err2 += pow( syst_err_sum[si], 2. ) ;
      } // si

      float total_err = sqrt( total_syst_err2 + sum_nldp_err2 + sum_mcc_err2 ) ;
      binsum_err = total_err ;

      binsum_err_nocorrelations = sqrt( sum_err2_nocorrelations ) ;


   } // calc_sum


  //===============================================================================================================









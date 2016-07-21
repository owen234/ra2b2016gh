

#include "TH1F.h"
#include "TSystem.h"
#include <fstream>



      int sum_index_array[209] ;
      int sum_nbins ;
      char sum_name[100] ;

      float par_val_ht[5] ;
      float par_err_ht_fit[5] ;
      float par_err_ht_syst[5] ;

      float par_val_njet[5] ;
      float par_err_njet_fit[5] ;
      float par_err_njet_syst[5] ;

      float par_val_mht_hth[6] ;
      float par_err_mht_hth[6] ;

      float par_val_mht_htm[6] ;
      float par_err_mht_htm[6] ;

      float par_val_mht_htl[6] ;
      float par_err_mht_htl[6] ;

      float par_val_nb[5] ;
      float par_err_nb[5] ;


      float qcd_ldp_val[209] ;
      float qcd_ldp_err[209] ;
      float qcd_hdp_val[209] ;
      float qcd_hdp_err[209] ;
      char  qcd_label[209][100] ;

      int   ht_ind[209] ;
      int   nj_ind[209] ;
      int   mht_ind[209] ;
      int   nb_ind[209] ;

      int nb_nj(4) ;
      int nb_nb(4) ;
      int nb_htmht(13) ;
      int nb_ht(3) ;

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

   void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2 ) ;
   void set_ht_and_mht_ind_from_htmht_ind( int bi_htmht, int& bi_ht, int& bi_mht ) ;
   void read_pars( const char* model_pars_file ) ;
   void read_qcd( const char* qcd_file ) ;
   void calc_sum( bool verb ) ;
   int  translate_208_to_160( int bi_208 ) ;
   void print_160bin_table() ;


#include "histio.c"

  //----

   void closure_sums3(
          bool verb = false,
          bool arg_include_qcdmc_correction = false,
          const char* model_pars_file = "model-pars-qcdmc3.txt",
          const char* model_ratio_hist_file = "outputfiles/model-ratio-hist1.root",
          const char* qcd_file  = "outputfiles/qcdmc-counts.txt"
       ) {

      include_qcdmc_correction = arg_include_qcdmc_correction ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( model_ratio_hist_file ) ;
      h_ratio_qcdmc_minus_model = (TH1F*) gDirectory -> FindObject( "h_ratio_qcdmc_minus_model" ) ;
      if ( h_ratio_qcdmc_minus_model == 0x0 ) { printf("\n\n *** missing h_ratio_qcdmc_minus_model.\n\n" ) ; return ; }

      read_pars( model_pars_file ) ;
      read_qcd( qcd_file ) ;

      //skip_high_weight_search_bins = true ;
      skip_high_weight_search_bins = false ;

      print_160bin_table() ;




     //---- 3 HT bins

      TH1F* h_closure_ht_model_fit   = new TH1F( "h_closure_ht_model_fit", "HT bins, model", 3, 0.5, 3.5 ) ;
      TH1F* h_closure_ht_model_fit_and_syst   = new TH1F( "h_closure_ht_model_fit_and_syst", "HT bins, model", 3, 0.5, 3.5 ) ;
      TH1F* h_closure_ht_model_total = new TH1F( "h_closure_ht_model_total", "HT bins, model", 3, 0.5, 3.5 ) ;
      TH1F* h_closure_ht_qcdmc = new TH1F( "h_closure_ht_qcdmc", "HT bins, QCD MC", 3, 0.5, 3.5 ) ;

      {
         sprintf( sum_name, "Search HT3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(3) ;
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


      }


      {
         sprintf( sum_name, "Search HT2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(2) ;
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

      }


      {
         sprintf( sum_name, "Search HT1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(1) ;
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

      }













     //---- 4 Njet bins

      TH1F* h_closure_njet_model_fit   = new TH1F( "h_closure_njet_model_fit", "njet bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_njet_model_fit_and_syst   = new TH1F( "h_closure_njet_model_fit_and_syst", "njet bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_njet_model_total = new TH1F( "h_closure_njet_model_total", "njet bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_njet_qcdmc = new TH1F( "h_closure_njet_qcdmc", "njet bins, QCD MC", 4, 0.5, 4.5 ) ;


      {
         sprintf( sum_name, "Search Nj1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nj != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(1) ;
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

      }

      {
         sprintf( sum_name, "Search Nj2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nj != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(2) ;
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

      }

      {
         sprintf( sum_name, "Search Nj3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nj != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(3) ;
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

      }

      {
         sprintf( sum_name, "Search Nj4" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nj != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(4) ;
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

      }





      {
         sprintf( sum_name, "Search HT3, Nj1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nj != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }

      {
         sprintf( sum_name, "Search HT3, Nj2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nj != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }

      {
         sprintf( sum_name, "Search HT3, Nj3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nj != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }

      {
         sprintf( sum_name, "Search HT3, Nj4" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nj != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }













     //---- 4 Nb bins

      TH1F* h_closure_nb_model_fit   = new TH1F( "h_closure_nb_model_fit", "nb bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_nb_model_fit_and_syst   = new TH1F( "h_closure_nb_model_fit_and_syst", "nb bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_nb_model_total = new TH1F( "h_closure_nb_model_total", "nb bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_nb_qcdmc = new TH1F( "h_closure_nb_qcdmc", "nb bins, QCD MC", 4, 0.5, 4.5 ) ;

      {
         sprintf( sum_name, "Search Nb0" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nb != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(1) ;
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

      }


      {
         sprintf( sum_name, "Search Nb1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nb != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(2) ;
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

      }


      {
         sprintf( sum_name, "Search Nb2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nb != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(3) ;
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

      }


      {
         sprintf( sum_name, "Search Nb3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_nb != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(4) ;
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

      }









      {
         sprintf( sum_name, "Search HT3, Nb0" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nb != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }



      {
         sprintf( sum_name, "Search HT3, Nb1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nb != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }



      {
         sprintf( sum_name, "Search HT3, Nb2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nb != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }



      {
         sprintf( sum_name, "Search HT3, Nb3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_nb != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }



     //--- 4 MHT bins

      TH1F* h_closure_mht_model_fit   = new TH1F( "h_closure_mht_model_fit", "mht bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_mht_model_fit_and_syst   = new TH1F( "h_closure_mht_model_fit_and_syst", "mht bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_mht_model_total = new TH1F( "h_closure_mht_model_total", "mht bins, model", 4, 0.5, 4.5 ) ;
      TH1F* h_closure_mht_qcdmc = new TH1F( "h_closure_mht_qcdmc", "mht bins, QCD MC", 4, 0.5, 4.5 ) ;

      {
         sprintf( sum_name, "Search MHT1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_mht != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(1) ;
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

      {
         sprintf( sum_name, "Search MHT2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_mht != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(2) ;
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

      {
         sprintf( sum_name, "Search MHT3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_mht != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(3) ;
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

      {
         sprintf( sum_name, "Search MHT4" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_mht != 5 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(4) ;
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










      {
         sprintf( sum_name, "Search HT3, MHT1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }


      {
         sprintf( sum_name, "Search HT3, MHT2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }


      {
         sprintf( sum_name, "Search HT3, MHT3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }



      {
         sprintf( sum_name, "Search HT3, MHT4" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 3 ) continue ;
                  if ( bi_mht != 5 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }


      {
         sprintf( sum_name, "Search HT2, Nj1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 2 ) continue ;
                  if ( bi_nj != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }

      {
         sprintf( sum_name, "Search HT2, Nj2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 2 ) continue ;
                  if ( bi_nj != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }


      {
         sprintf( sum_name, "Search HT2, Nj3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 2 ) continue ;
                  if ( bi_nj != 3 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }


      {
         sprintf( sum_name, "Search HT2, Nj4" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 2 ) continue ;
                  if ( bi_nj != 4 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }

      {
         sprintf( sum_name, "Search HT1, Nj1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 1 ) continue ;
                  if ( bi_nj != 1 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }

      {
         sprintf( sum_name, "Search HT1, Nj2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_ht != 1 ) continue ;
                  if ( bi_nj != 2 ) continue ;
                  if ( bi_mht == 1 ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

      }









     //--- boxes

      TH1F* h_closure_10boxes_model_fit   = new TH1F( "h_closure_10boxes_model_fit", "10 MHT-HT boxes, model", 10, 0.5, 10.5 ) ;
      TH1F* h_closure_10boxes_model_fit_and_syst   = new TH1F( "h_closure_10boxes_model_fit_and_syst", "10 MHT-HT boxes, model", 10, 0.5, 10.5 ) ;
      TH1F* h_closure_10boxes_model_total = new TH1F( "h_closure_10boxes_model_total", "10 MHT-HT boxes, model", 10, 0.5, 10.5 ) ;
      TH1F* h_closure_10boxes_qcdmc = new TH1F( "h_closure_10boxes_qcdmc", "10 MHT-HT boxes, QCD MC", 10, 0.5, 10.5 ) ;

      {
         sprintf( sum_name, "Search box1" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+1) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(1) ;
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


      {
         sprintf( sum_name, "Search box2" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+2) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(2) ;
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


      {
         sprintf( sum_name, "Search box3" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+3) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(3) ;
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

      {
         sprintf( sum_name, "Search box4" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+4) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(4) ;
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

      {
         sprintf( sum_name, "Search box5" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+5) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(5) ;
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

      {
         sprintf( sum_name, "Search box6" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+6) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(6) ;
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


      {
         sprintf( sum_name, "Search box7" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+7) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(7) ;
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


      {
         sprintf( sum_name, "Search box8" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+8) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(8) ;
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



      {
         sprintf( sum_name, "Search box9" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+9) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(9) ;
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

      {
         sprintf( sum_name, "Search box10" ) ;
         sum_nbins = 0 ;
         int bi_hist(0) ;
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
               for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

                  bi_hist ++ ;
                  int bi_ht, bi_mht ;

                  set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

                  if ( bi_htmht != (3+10) ) continue ;

                  sum_index_array[sum_nbins] = bi_hist ;
                  sum_nbins++ ;

               } // bi_htmht
            } // bi_nb
         } // bi_nj

         printf("\n\n Calculating %s\n", sum_name ) ;
         calc_sum( verb ) ;

         int hbi(10) ;
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

  /// for ( int bi=1; bi<=208; bi++ ) {
  ///    translate_208_to_160( bi ) ;
  /// } // bi

   } // closure_sums3

  //=======================================================================================

   void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2 ) {

      val = 0. ;
      err1 = 0. ;
      err2 = 0. ;

      ifs.seekg(0) ;

      TString line ;
      while ( ifs.good() ) {
         line.ReadLine( ifs ) ;
         char line_parname[100] ;
         float line_val, line_err1, line_err2 ;
         sscanf( line.Data(), "%s  %f %f %f", line_parname, &line_val, &line_err1, &line_err2 ) ;
         if ( strcmp( line_parname, pname ) == 0 ) {
            val = line_val ;
            err1 = line_err1 ;
            err2 = line_err2 ;
            return ;
         }
      }

      printf("\n\n *** get_par : Failed to find parameter %s\n\n", pname ) ;
      gSystem -> Exit(-1) ;

   } // get_par

  //=======================================================================================

   void set_ht_and_mht_ind_from_htmht_ind( int bi_htmht, int& bi_ht, int& bi_mht ) {

      if ( bi_htmht < 1 || bi_htmht > 13 ) {
         printf("\n\n wtf???\n\n") ;
         gSystem -> Exit(-1) ;
      }

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

   } // set_ht_and_mht_ind_from_htmht_ind

  //=======================================================================================

      void read_pars( const char* model_pars_file ) {

         ifstream ifs_model_pars ;
         ifs_model_pars.open( model_pars_file ) ;
         if ( !ifs_model_pars.good() ) { printf("\n\n *** Problem opening %s\n\n", model_pars_file ) ; gSystem->Exit(-1) ; }

         float val, err1, err2 ;
         char pname[100] ;

       //---
         sprintf( pname, "Kqcd_HT1" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_ht[1] = val ;
         par_err_ht_fit[1] = err1 ;
         par_err_ht_syst[1] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Kqcd_HT2" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_ht[2] = val ;
         par_err_ht_fit[2] = err1 ;
         par_err_ht_syst[2] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Kqcd_HT3" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_ht[3] = val ;
         par_err_ht_fit[3] = err1 ;
         par_err_ht_syst[3] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

       //---
         sprintf( pname, "Sqcd_njet1" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_njet[1] = val ;
         par_err_njet_fit[1] = err1 ;
         par_err_njet_syst[1] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_njet2" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_njet[2] = val ;
         par_err_njet_fit[2] = err1 ;
         par_err_njet_syst[2] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_njet3" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_njet[3] = val ;
         par_err_njet_fit[3] = err1 ;
         par_err_njet_syst[3] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_njet4" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_njet[4] = val ;
         par_err_njet_fit[4] = err1 ;
         par_err_njet_syst[4] = val*err2 ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;



       //---
         sprintf( pname, "Sqcd_mhtc_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_hth[1] = val ;
         par_err_mht_hth[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht1_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_hth[2] = val ;
         par_err_mht_hth[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht2_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_hth[3] = val ;
         par_err_mht_hth[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht3_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_hth[4] = val ;
         par_err_mht_hth[4] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht4_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_hth[5] = val ;
         par_err_mht_hth[5] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

       //---
         sprintf( pname, "Sqcd_mhtc_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htm[1] = val ;
         par_err_mht_htm[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht1_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htm[2] = val ;
         par_err_mht_htm[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht2_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htm[3] = val ;
         par_err_mht_htm[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht3_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htm[4] = val ;
         par_err_mht_htm[4] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht4_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htm[5] = val ;
         par_err_mht_htm[5] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

       //---
         sprintf( pname, "Sqcd_mhtc_htl" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htl[1] = val ;
         par_err_mht_htl[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht1_htl" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htl[2] = val ;
         par_err_mht_htl[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht2_htl" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_mht_htl[3] = val ;
         par_err_mht_htl[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;





       //---
         sprintf( pname, "Sqcd_nb0" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_nb[1] = val ;
         par_err_nb[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_nb1" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_nb[2] = val ;
         par_err_nb[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_nb2" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_nb[3] = val ;
         par_err_nb[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_nb3" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2 ) ;
         par_val_nb[4] = val ;
         par_err_nb[4] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

      } // read_pars

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

               bi_hist ++ ;
               if ( bi_htmht <= 3 ) {
                  bi_control ++ ;
               } else {
                  bi_search ++ ;
               }

               int bi_ht, bi_mht ;

               set_ht_and_mht_ind_from_htmht_ind( bi_htmht, bi_ht, bi_mht ) ;

               ht_ind[bi_hist]  = bi_ht ;
               nj_ind[bi_hist]  = bi_nj ;
               mht_ind[bi_hist] = bi_mht ;
               nb_ind[bi_hist]  = bi_nb ;

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

               TString line ;
               line.ReadLine( ifs_qcd_ldp ) ;
               int rbi1, rbi2 ;
               char rsc[10] ;
               char rlabel[100] ;
               float ldp_val, ldp_err, hdp_val, hdp_err ;
               sscanf( line.Data(), " %d %s %d %s  %f %f %f %f", &rbi1, rsc, &rbi2, rlabel, &ldp_val, &ldp_err, &hdp_val, &hdp_err ) ;
               if ( rbi1 != bi_hist ) { printf("\n\n *** inconsistency!  %d : %s\n\n", bi_hist, line.Data() ) ; return ; }

               qcd_ldp_val[bi_hist] = ldp_val ;
               qcd_ldp_err[bi_hist] = ldp_err ;
               qcd_hdp_val[bi_hist] = hdp_val ;
               qcd_hdp_err[bi_hist] = hdp_err ;
               sprintf( qcd_label[bi_hist], "%s", label ) ;

               printf( " %30s   HT ind %d,  Nj ind %d,  MHT ind %d, Nb ind %d    %8.1f  %8.1f   %8.1f  %8.1f\n", label,
                  ht_ind[bi_hist], nj_ind[bi_hist], mht_ind[bi_hist], nb_ind[bi_hist],
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
         bool skip_this_search_bin[209] ;

         for ( int i=1; i<=208; i++ ) { skip_this_search_bin[i] = false ; }

         for ( int sbi=0; sbi<sum_nbins; sbi++ ) {

            int hbi = sum_index_array[sbi] ;

            if ( skip_high_weight_search_bins && qcd_hdp_err[hbi] > 10. ) {
               skip_this_search_bin[hbi] = true ;
               printf("  Skipping %3d  %30s  %8.1f +/- %8.1f\n", hbi, qcd_label[hbi], qcd_hdp_val[hbi], qcd_hdp_err[hbi] ) ;
               continue ;
            }

           //--- skip HT1 at high Njet.
            if ( ht_ind[hbi] == 1 && nj_ind[hbi] > 2 ) continue ;


            float model_ratio(0.) ;
            if ( ht_ind[hbi]==3 ) {
               model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_hth[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
            }
            if ( ht_ind[hbi]==2 ) {
               model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_htm[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
            }
            if ( ht_ind[hbi]==1 ) {
               model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_htl[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
            }

            float ratio_correction_val(0.) ;
            float ratio_correction_err(0.) ;
            {
               int bin_160 = translate_208_to_160( hbi ) ;
               ratio_correction_val = h_ratio_qcdmc_minus_model -> GetBinContent( bin_160 ) ;
               ratio_correction_err = h_ratio_qcdmc_minus_model -> GetBinError( bin_160 ) ;
               if ( verb ) {
                  printf("  Correcting 208 bin %3d (%30s) with 160 bin %d (%30s) : correction = %6.4f +/- %6.4f\n",
                      hbi, qcd_label[hbi], bin_160, h_ratio_qcdmc_minus_model -> GetXaxis() -> GetBinLabel( bin_160 ),
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








         float partial_sum[36] ;
         int pi(0) ;
         for ( int i=0; i<36; i++ ) { partial_sum[i] = 0. ; }

         for ( int sbi=0; sbi<sum_nbins; sbi++ ) {

            int hbi = sum_index_array[sbi] ;

            if ( skip_this_search_bin[hbi] ) continue ;

           //--- skip HT1 at high Njet.
            if ( ht_ind[hbi] == 1 && nj_ind[hbi] > 2 ) continue ;

            float model_ratio(0.) ;
            if ( ht_ind[hbi] == 3 ) {
               model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_hth[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
            }
            if ( ht_ind[hbi] == 2 ) {
               model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_htm[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
            }
            if ( ht_ind[hbi] == 1 ) {
               model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_htl[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
            }
            float model_val = qcd_ldp_val[hbi] * model_ratio ;

          //---

            pi = 0 ;
            if ( ht_ind[hbi] == 1 && par_val_ht[1] > 0. ) {
               partial_sum[pi] += model_val / par_val_ht[1] ;
            }

            pi++ ;
            if ( ht_ind[hbi] == 2 && par_val_ht[2] > 0. ) {
               partial_sum[pi] += model_val / par_val_ht[2] ;
            }

            pi++ ;
            if ( ht_ind[hbi] == 3 && par_val_ht[3] > 0. ) {
               partial_sum[pi] += model_val / par_val_ht[3] ;
            }

          //---

            pi++ ;
            if ( nj_ind[hbi] == 1 && par_val_njet[1] > 0. ) {
               partial_sum[pi] += model_val / par_val_njet[1] ;
            }

            pi++ ;
            if ( nj_ind[hbi] == 2 && par_val_njet[2] > 0. ) {
               partial_sum[pi] += model_val / par_val_njet[2] ;
            }

            pi++ ;
            if ( nj_ind[hbi] == 3 && par_val_njet[3] > 0. ) {
               partial_sum[pi] += model_val / par_val_njet[3] ;
            }

            pi++ ;
            if ( nj_ind[hbi] == 4 && par_val_njet[4] > 0. ) {
               partial_sum[pi] += model_val / par_val_njet[4] ;
            }

          //---

            pi++ ;
            if ( mht_ind[hbi] == 1 && ht_ind[hbi] == 3 && par_val_mht_hth[1] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_hth[1] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 2 && ht_ind[hbi] == 3 && par_val_mht_hth[2] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_hth[2] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 3 && ht_ind[hbi] == 3 && par_val_mht_hth[3] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_hth[3] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 4 && ht_ind[hbi] == 3 && par_val_mht_hth[4] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_hth[4] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 5 && ht_ind[hbi] == 3 && par_val_mht_hth[5] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_hth[5] ;
            }

          //---

            pi++ ;
            if ( mht_ind[hbi] == 1 && ht_ind[hbi] == 2 && par_val_mht_htm[1] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htm[1] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 2 && ht_ind[hbi] == 2 && par_val_mht_htm[2] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htm[2] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 3 && ht_ind[hbi] == 2 && par_val_mht_htm[3] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htm[3] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 4 && ht_ind[hbi] == 2 && par_val_mht_htm[4] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htm[4] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 5 && ht_ind[hbi] == 2 && par_val_mht_htm[5] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htm[5] ;
            }

          //---

            pi++ ;
            if ( mht_ind[hbi] == 1 && ht_ind[hbi] == 1 && par_val_mht_htl[1] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htl[1] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 2 && ht_ind[hbi] == 1 && par_val_mht_htl[2] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htl[2] ;
            }

            pi++ ;
            if ( mht_ind[hbi] == 3 && ht_ind[hbi] == 1 && par_val_mht_htl[3] > 0. ) {
               partial_sum[pi] += model_val / par_val_mht_htl[3] ;
            }

          //---

            pi++ ;
            if ( nb_ind[hbi] == 1 && par_val_nb[1] > 0. ) {
               partial_sum[pi] += model_val / par_val_nb[1] ;
            }

            pi++ ;
            if ( nb_ind[hbi] == 2 && par_val_nb[2] > 0. ) {
               partial_sum[pi] += model_val / par_val_nb[2] ;
            }

            pi++ ;
            if ( nb_ind[hbi] == 3 && par_val_nb[3] > 0. ) {
               partial_sum[pi] += model_val / par_val_nb[3] ;
            }

            pi++ ;
            if ( nb_ind[hbi] == 4 && par_val_nb[4] > 0. ) {
               partial_sum[pi] += model_val / par_val_nb[4] ;
            }

         } // sbi

         float par_err_fit[36] ;
         float par_err_syst[36] ;
         char  par_name[36][100] ;

         int npar(0) ;
         { int pi(0) ;
            par_err_fit[pi] = par_err_ht_fit[1] ; par_err_syst[pi] = par_err_ht_syst[1] ; sprintf( par_name[pi], "HT1" ) ; pi++ ;
            par_err_fit[pi] = par_err_ht_fit[2] ; par_err_syst[pi] = par_err_ht_syst[2] ; sprintf( par_name[pi], "HT2" ) ; pi++ ;
            par_err_fit[pi] = par_err_ht_fit[3] ; par_err_syst[pi] = par_err_ht_syst[3] ; sprintf( par_name[pi], "HT3" ) ; pi++ ;
            par_err_fit[pi] = par_err_njet_fit[1] ; par_err_syst[pi] = par_err_njet_syst[1] ; sprintf( par_name[pi], "Njet1" ) ; pi++ ;
            par_err_fit[pi] = par_err_njet_fit[2] ; par_err_syst[pi] = par_err_njet_syst[2] ; sprintf( par_name[pi], "Njet2" ) ; pi++ ;
            par_err_fit[pi] = par_err_njet_fit[3] ; par_err_syst[pi] = par_err_njet_syst[3] ; sprintf( par_name[pi], "Njet3" ) ; pi++ ;
            par_err_fit[pi] = par_err_njet_fit[4] ; par_err_syst[pi] = par_err_njet_syst[4] ; sprintf( par_name[pi], "Njet4" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_hth[1] ; sprintf( par_name[pi], "MHTC-HTH" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_hth[2] ; sprintf( par_name[pi], "MHT1-HTH" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_hth[3] ; sprintf( par_name[pi], "MHT2-HTH" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_hth[4] ; sprintf( par_name[pi], "MHT3-HTH" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_hth[5] ; sprintf( par_name[pi], "MHT4-HTH" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htm[1] ; sprintf( par_name[pi], "MHTC-HTM" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htm[2] ; sprintf( par_name[pi], "MHT1-HTM" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htm[3] ; sprintf( par_name[pi], "MHT2-HTM" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htm[4] ; sprintf( par_name[pi], "MHT3-HTM" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htm[5] ; sprintf( par_name[pi], "MHT4-HTM" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htl[1] ; sprintf( par_name[pi], "MHTC-HTL" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htl[2] ; sprintf( par_name[pi], "MHT1-HTL" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_mht_htl[3] ; sprintf( par_name[pi], "MHT2-HTL" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_nb[1] ; sprintf( par_name[pi], "Nb0" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_nb[2] ; sprintf( par_name[pi], "Nb1" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_nb[3] ; sprintf( par_name[pi], "Nb2" ) ; pi++ ;
            par_err_fit[pi] = 0. ; par_err_syst[pi] = par_err_nb[4] ; sprintf( par_name[pi], "Nb3" ) ; pi++ ;
            npar = pi ;
         }

         float model_sum_err2_fit(0.), model_sum_err2_syst(0.) ;
         for ( int pi=0; pi<npar; pi++ ) {
            float par_err = sqrt( pow( par_err_fit[pi], 2.) + pow( par_err_syst[pi], 2. ) ) ;
            if (verb) {
               printf("  %2d %5s : partial sum %8.1f, err %8.4f (%8.4f (fit) +/- %8.4f (syst)), contrib %8.1f\n",
                  pi, par_name[pi], partial_sum[pi], par_err, par_err_fit[pi], par_err_syst[pi],
                  partial_sum[pi] * par_err ) ;
            }
            model_sum_err2_fit += pow( partial_sum[pi] * par_err_fit[pi], 2. ) ;
            model_sum_err2_syst += pow( partial_sum[pi] * par_err_syst[pi], 2. ) ;
         } // bi

         model_sum_err_fit = sqrt( model_sum_err2_fit ) ;
         model_sum_err_mc = sqrt( model_sum_correction_err2 ) ;
         model_sum_err_syst = sqrt( model_sum_err2_syst ) ;
         model_sum_err = sqrt( model_sum_err2_fit + model_sum_correction_err2 + model_sum_err2_syst ) ;

         if ( include_qcdmc_correction ) {
            model_sum_val += model_sum_correction_val ;
         }

         printf("  Sums:  model = %8.1f +/- %8.1f (%8.1f (fit), %8.1f (syst))     QCD = %8.1f +/- %8.1f\n",
            model_sum_val, model_sum_err, model_sum_err_fit, model_sum_err_syst,
            qcd_sum_val, qcd_sum_err ) ;

     } // calc_sum

  //=======================================================================================

   int  translate_208_to_160( int bi_208 ) {

      if ( bi_208 < 1 || bi_208 > 208 ) return -1 ;
      int bi_nj = (bi_208-1)/(4*13) + 1 ;
      int bi_nb = ((bi_208-1)%(4*13))/13 + 1 ;
      int bi_htmht13 = (bi_208-1)%13+1 ;

      int bi_160 = (bi_nj-1)*(4*10) + (bi_nb-1)*(10) + bi_htmht13-3 ;
      if ( bi_htmht13 < 4 ) bi_160 = -1 ;

   // if ( bi_htmht13 <= 3 ) {
   //    printf("  208 bin : %3d  ,  Nj %d,  Nb %d,  HTMHT %d\n", bi_208, bi_nj, bi_nb, bi_htmht13 ) ;
   // } else {
   //    printf("  208 bin : %3d  ,  Nj %d,  Nb %d,  HTMHT %d,  160 bin %3d\n", bi_208, bi_nj, bi_nb, bi_htmht13, bi_160 ) ;
   // }

      return bi_160 ;

   } // translate_208_to_160

  //=======================================================================================

   void print_160bin_table() {

      printf("\n\n") ;

      TH1F* h_ratio_all = (TH1F*) gDirectory -> FindObject( "h_ratio_all" ) ;
      if ( h_ratio_all == 0x0 ) { printf("\n\n *** Missing h_ratio_all hist.\n\n") ; gSystem->Exit(-1) ; }

      for ( int hbi=1; hbi<=208; hbi++ ) {

         int bi_160 = translate_208_to_160( hbi ) ;
         if ( bi_160 < 1 ) continue ;

        //--- skip HT1 at high Njet.
         if ( ht_ind[hbi] == 1 && nj_ind[hbi] > 2 ) continue ;


         float model_ratio(0.) ;
         if ( ht_ind[hbi]==3 ) {
            model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_hth[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
         }
         if ( ht_ind[hbi]==2 ) {
            model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_htm[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
         }
         if ( ht_ind[hbi]==1 ) {
            model_ratio =  par_val_ht[ht_ind[hbi]] * par_val_njet[nj_ind[hbi]] * par_val_mht_htl[mht_ind[hbi]] * par_val_nb[nb_ind[hbi]] ;
         }

         float model_ratio_from_hist = h_ratio_all -> GetBinContent( bi_160 ) ;

    ///  printf("  %3d   %30s :  arrays Ht%d Nj%d MHT%d Nb%d : vals HT %6.4f Nj %6.4f MHT %6.4f Nb %6.4f : Model ratio, calc %7.5f  (from hist %7.5f  %30s)\n",
    ///    bi_160, qcd_label[hbi],
    ///    ht_ind[hbi], nj_ind[hbi], mht_ind[hbi], nb_ind[hbi],
    ///    par_val_ht[ht_ind[hbi]], par_val_njet[nj_ind[hbi]], par_val_mht_hth[mht_ind[hbi]], par_val_nb[nb_ind[hbi]],
    ///    model_ratio, model_ratio_from_hist, h_ratio_all -> GetXaxis() -> GetBinLabel( bi_160 ) ) ;

         if ( fabs(model_ratio-model_ratio_from_hist) > 0.00001 ) { printf("\n *** inconsistent model ratios.\n\n" ) ; }

         float ratio_correction_val(0.) ;
         float ratio_correction_err(0.) ;
         {
            ratio_correction_val = h_ratio_qcdmc_minus_model -> GetBinContent( bi_160 ) ;
            ratio_correction_err = h_ratio_qcdmc_minus_model -> GetBinError( bi_160 ) ;
         }
         float model_val = qcd_ldp_val[hbi] * model_ratio ;

         float model_correction_val = qcd_ldp_val[hbi] * ratio_correction_val ;
         float model_correction_err = qcd_ldp_val[hbi] * ratio_correction_err ;

         printf("  %3d   %30s :  Model %6.1f +/- %6.1f,  cor = %9.3f * %7.4f = %6.1f +/- %6.1f,  total = %6.1f +/- %6.1f,   HDP QCD MC  %6.1f +/- %6.1f  ",
             bi_160, qcd_label[hbi],
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

   } // print_160bin_table


  //=======================================================================================







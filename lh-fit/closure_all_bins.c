
#include "../binning.h"
#include "../histio.c"
#include "../get_hist.h"
#include "TSystem.h"

#include <fstream>

   //-- bi_* variables all start counting from 1, not 0.

   bool decode_bin_label( const char* bin_label, int& bi_nj, int& bi_nb, int& bi_mht, int& bi_ht ) ;
   void get_smht_par( ifstream &ifs, const char* pname, double &val, double &err ) ;

  //--------

   void closure_all_bins( const char* qcdmc_hist_file = "../outputfiles/hists-v2d-qcd.root",
                          const char* fit_pars_file = "../outputfiles/kqcd-parameters-from-qcdmc.txt",
                          //const char* fit_pars_file = "outputfiles/lhfit-results-ws-lhfit-test3/kqcd-parameter-fit-results.txt",
                          const char* mht_pars_file = "../outputfiles/model-pars-data3.txt"
   ) {


      setup_bins() ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( qcdmc_hist_file ) ;

      TH1F* h_ldp = get_hist( "h_ldp" ) ;
      TH1F* h_hdp = get_hist( "h_hdp" ) ;



      int search_bin_count(0) ;

      for ( int bi=1; bi<=(h_ldp->GetNbinsX()); bi++ ) {

         char bin_label[100] ;
         sscanf( h_ldp->GetXaxis()->GetBinLabel(bi) , "%s", bin_label ) ;

         int bi_nj, bi_nb, bi_mht, bi_ht ;
         bool excluded = decode_bin_label( bin_label, bi_nj, bi_nb, bi_mht, bi_ht ) ;

         if ( excluded ) continue ;
         if ( bi_mht==1 ) continue ;

         search_bin_count ++ ;

         printf("  %3d %3d : %s : %d, %d, %d, %d %s\n", search_bin_count, bi, bin_label, bi_nj, bi_nb, bi_mht, bi_ht, (excluded?"excluded":"ok") ) ;

      } // bi




    //----

      printf("\n\n Fit parameters:\n") ;

      double fit_rqcd_val[10][10] ;
      double fit_rqcd_err[10][10] ;
      for ( int nj=1; nj<=nb_nj; nj++ ) {
         for ( int ht=1; ht<=nb_ht[1]; ht++ ) {
            fit_rqcd_val[nj][ht] = 0. ;
            fit_rqcd_err[nj][ht] = 0. ;
         }
      }

      ifstream ifs_fit_pars ;
      ifs_fit_pars.open( fit_pars_file ) ;
      if ( !ifs_fit_pars.good() ) { printf("\n\n *** bad fit pars file: %s\n\n\n", fit_pars_file ) ; gSystem -> Exit(-1) ; }
      while ( ifs_fit_pars.good() ) {
         TString line ;
         line.ReadLine( ifs_fit_pars ) ;
         if ( !ifs_fit_pars.good() ) break ;
         int bi_nj, bi_ht ;
         double val, err ;
         sscanf( line.Data(), " R_qcd_ldp_Nj%d_HT%d %lf +/- %lf", &bi_nj, &bi_ht, &val, &err ) ;
         printf( "%s :  bi_nj=%d, bi_ht=%d, val = %6.4f, err = %6.4f\n", line.Data(), bi_nj, bi_ht, val, err ) ;
         fit_rqcd_val[bi_nj][bi_ht] = val ;
         fit_rqcd_err[bi_nj][bi_ht] = err ;
      }
      ifs_fit_pars.close() ;

      printf("\n\n") ;


    //----

      printf("\n\n MHT parameters:\n") ;

      double smht_val[10][10] ;
      double smht_err[10][10] ;
      for ( int mht=1; mht<=nb_mht; mht++ ) {
         for ( int ht=1; ht<=nb_ht[1]; ht++ ) {
            smht_val[mht][ht] = 0. ;
            smht_err[mht][ht] = 0. ;
         } // ht
      } // mit

      ifstream ifs_mht_pars ;
      ifs_mht_pars.open( mht_pars_file ) ;
      if ( !ifs_mht_pars.good() ) { printf("\n\n *** Bad MHT pars file : %s\n\n\n", mht_pars_file ) ; gSystem -> Exit(-1) ; }

      get_smht_par( ifs_mht_pars, "Sqcd_mht1_hth", smht_val[2][3], smht_err[2][3] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht2_hth", smht_val[3][3], smht_err[3][3] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht3_hth", smht_val[4][3], smht_err[4][3] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht4_hth", smht_val[5][3], smht_err[5][3] ) ;

      get_smht_par( ifs_mht_pars, "Sqcd_mht1_htm", smht_val[2][2], smht_err[2][2] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht2_htm", smht_val[3][2], smht_err[3][2] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht3_htm", smht_val[4][2], smht_err[4][2] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht4_htm", smht_val[5][2], smht_err[5][2] ) ;

      get_smht_par( ifs_mht_pars, "Sqcd_mht1_htl", smht_val[2][1], smht_err[2][1] ) ;
      get_smht_par( ifs_mht_pars, "Sqcd_mht2_htl", smht_val[3][1], smht_err[3][1] ) ;

      printf("\n\n") ;


    //----




      TH1F* h_ldp_search = new TH1F( "h_ldp_search", "LDP, search bins", search_bin_count, 0.5, search_bin_count + 0.5 ) ;
      TH1F* h_hdp_search = new TH1F( "h_hdp_search", "HDP, search bins", search_bin_count, 0.5, search_bin_count + 0.5 ) ;
      TH1F* h_rqcd = new TH1F( "h_rqcd", "Rqcd", search_bin_count, 0.5, search_bin_count + 0.5 ) ;
      TH1F* h_hdp_model_search = new TH1F( "h_hdp_model_search", "HDP, model prediction, search bins", search_bin_count, 0.5, search_bin_count + 0.5 ) ;

      int sbi(0) ;
      for ( int bi=1; bi<=(h_ldp->GetNbinsX()); bi++ ) {

         char bin_label[100] ;
         sscanf( h_ldp->GetXaxis()->GetBinLabel(bi) , "%s", bin_label ) ;

         int bi_nj, bi_nb, bi_mht, bi_ht ;
         bool excluded = decode_bin_label( bin_label, bi_nj, bi_nb, bi_mht, bi_ht ) ;

         if ( excluded ) continue ;
         if ( bi_mht==1 ) continue ;

         sbi ++ ;

         char new_bin_label[100] ;
         sprintf( new_bin_label, "%s %3d", bin_label, sbi ) ;



         double rqcd_val = fit_rqcd_val[bi_nj][bi_ht] * smht_val[bi_mht][bi_ht] ;
         double rqcd_err(0.) ;

         if ( fit_rqcd_val[bi_nj][bi_ht] > 0 && smht_val[bi_mht][bi_ht] > 0 ) {
            rqcd_err = rqcd_val * sqrt( pow( fit_rqcd_err[bi_nj][bi_ht]/fit_rqcd_val[bi_nj][bi_ht], 2 )
                                      + pow( smht_err[bi_mht][bi_ht]/smht_val[bi_mht][bi_ht], 2 ) ) ;
         }

         printf(" %3d %s : Rqcd = (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) = %6.4f +/- %6.4f\n",
           sbi, bin_label, fit_rqcd_val[bi_nj][bi_ht], fit_rqcd_err[bi_nj][bi_ht],
                           smht_val[bi_mht][bi_ht], smht_err[bi_mht][bi_ht],
                           rqcd_val, rqcd_err ) ;

         h_ldp_search -> SetBinContent( sbi, h_ldp->GetBinContent( bi ) ) ;
         h_ldp_search -> SetBinError( sbi, h_ldp->GetBinError( bi ) ) ;

         h_hdp_search -> SetBinContent( sbi, h_hdp->GetBinContent( bi ) ) ;
         h_hdp_search -> SetBinError( sbi, h_hdp->GetBinError( bi ) ) ;

         h_rqcd -> SetBinContent( sbi, rqcd_val ) ;
         h_rqcd -> SetBinError( sbi, rqcd_err ) ;

         double ldp_val = h_ldp->GetBinContent( bi ) ;
         double ldp_err = h_ldp->GetBinError( bi ) ;

         double model_hdp_val = ldp_val * rqcd_val ;
         double model_hdp_err(0.) ;
         if ( ldp_val > 0 && rqcd_val > 0 ) {
            model_hdp_err = model_hdp_val * sqrt( pow( ldp_err/ldp_val, 2. ) + pow( rqcd_err/rqcd_val, 2. ) ) ;
         } else {
            model_hdp_err = (1.*rqcd_val) * sqrt( pow( ldp_err/1., 2. ) + pow( rqcd_err/rqcd_val, 2. ) ) ;
            printf(" %3d %s : ldp_val = %f\n", sbi, bin_label, ldp_val ) ;
         }

         h_hdp_model_search -> SetBinContent( sbi, model_hdp_val ) ;
         h_hdp_model_search -> SetBinError( sbi, model_hdp_err ) ;

         h_ldp_search -> GetXaxis() -> SetBinLabel( sbi, new_bin_label ) ;
         h_hdp_search -> GetXaxis() -> SetBinLabel( sbi, new_bin_label ) ;
         h_rqcd -> GetXaxis() -> SetBinLabel( sbi, new_bin_label ) ;
         h_hdp_model_search -> GetXaxis() -> SetBinLabel( sbi, new_bin_label ) ;



      } // bi

      TH1F* h_hdp_model_search_copy = (TH1F*) h_hdp_model_search->Clone( "h_hdp_model_search_copy" ) ;
      for ( int i=1; i<=search_bin_count; i++ ) { h_hdp_model_search_copy -> SetBinError( i, 0.0000001 ) ; }
      h_hdp_model_search_copy -> SetLineColor(4) ;

      h_ldp_search -> GetXaxis() -> LabelsOption( "v" ) ;
      h_hdp_search -> GetXaxis() -> LabelsOption( "v" ) ;
      h_rqcd -> GetXaxis() -> LabelsOption( "v" ) ;
      h_hdp_model_search -> GetXaxis() -> LabelsOption( "v" ) ;

      h_hdp_model_search -> SetFillColor( kRed-10 ) ;

      h_hdp_search -> SetMarkerStyle(20) ;

      h_hdp_search -> Draw() ;
      h_hdp_model_search -> Draw("E2 same") ;
      h_hdp_search -> Draw( "same" ) ;
      h_hdp_model_search_copy -> Draw("e same") ;

      gPad -> SetGridy(1) ;



   } // closure_all_bins

  //=======================================================================================================


    //-- bi_* variables all start counting from 1, not 0.

   bool decode_bin_label( const char* bin_label, int& bi_nj, int& bi_nb, int& bi_mht, int& bi_ht ) {
      char mhtchar[10] ;
      sscanf( bin_label, "NJets%d_BTags%d-MHT%1s-HT%d", &bi_nj, &bi_nb, mhtchar, &bi_ht ) ;
      if ( strcmp( mhtchar, "C" ) == 0 ) {
         bi_mht = 0 ;
      } else {
         sscanf( mhtchar, "%d", &bi_mht ) ;
      }
      bool excluded = is_this_bin_excluded( bi_nj, bi_nb, bi_ht-1, bi_mht ) ;
      //-- bi_* variables all start counting from 1, not 0.
      bi_nj ++ ;
      bi_nb ++ ;
      bi_mht ++ ;
      return excluded ;
   } // decode_bin_label

  //=======================================================================================================

   void get_smht_par( ifstream &ifs, const char* pname, double &val, double &err ) {

      ifs.seekg(0) ;
      while ( ifs.good() ) {
         TString line ;
         line.ReadLine( ifs ) ;
         if ( !ifs.good() ) break ;
         char line_pname[100] ;
         double line_val, line_stat_err, line_syst_err ;
         sscanf( line.Data(), "%s %lf %lf %lf", line_pname, &line_val, &line_stat_err, &line_syst_err ) ;
         if ( strcmp( line_pname, pname ) == 0 ) {
            val = line_val ;
            err = sqrt( line_stat_err*line_stat_err + (line_val*line_syst_err*line_val*line_syst_err) ) ;
            printf("  %s  val = %6.4f, err = %6.4f\n", pname, val, err ) ;
            return ;
         }
      }
      printf("\n\n *** get_smht_par : did not find parameter %s\n\n", pname ) ;
      gSystem -> Exit(-1) ;

   } // get_smht_par

  //=======================================================================================================











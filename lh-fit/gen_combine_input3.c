#ifndef gen_combine_input3_c

#include "TH1F.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "THStack.h"
#include "TPad.h"
#include "TLegend.h"

#include "../histio.c"
#include "../binning.h"
#include "../read_pars.h"
#include "read_np_matrix.h"


//--------

void gen_combine_input3(
   const char* model_pars_file          = "../outputfiles/model-pars-data3.txt",
   const char* model_pars_file2         = "outputfiles/lhfit-results-ws-lhfit3/kqcd-parameter-fit-results.txt",
   const char* np_matrix_file           = "outputfiles/nuisance-parameter-matrix.txt",
   const char* mc_minus_model_hist_file = "outputfiles/model-ratio-hist4.root",
   const char* data_file                = "../outputfiles/combine-input-data.txt",
   const char* lostlep_file             = "outputfiles/combine-input-lostlep.txt",
   const char* hadtau_file              = "outputfiles/combine-input-hadtau.txt",
   const char* znunu_file               = "outputfiles/combine-input-znunu.txt", 
   const char* combine_output_file      = "outputfiles/combine-input-all.txt",
   const char* qcdratio_output_file     = "outputfiles/qcd-ratios.txt",
   const char* hist_output_file         = "outputfiles/gci-output.root"
                       ) {

   gDirectory -> Delete( "h*" ) ;
   setup_bins();
   read_pars( model_pars_file ) ;
   read_pars2( model_pars_file2 ) ;
   read_np_matrix( np_matrix_file ) ;


   gStyle -> SetOptStat(0) ;
   gStyle -> SetPadGridY(1) ;

   bool print_dashes( true ) ;

   TString line ;

   TFile* tf_mc_minus_model = new TFile( mc_minus_model_hist_file, "READ" ) ;
   if ( tf_mc_minus_model == 0x0 ) { printf("\n\n *** bad mc_minus_model_hist_file file %s\n\n", mc_minus_model_hist_file ) ; return ; }
   if ( ! (tf_mc_minus_model->IsOpen()) ) { printf("\n\n *** bad mc_minus_model_hist_file file %s\n\n", mc_minus_model_hist_file ) ; return ; }
   TH1F* h_ratio_qcdmc_minus_model = (TH1F*) tf_mc_minus_model -> Get( "h_ratio_qcdmc_minus_model" ) ;
   if ( h_ratio_qcdmc_minus_model == 0x0 ) { printf("\n\n *** Can't find h_ratio_qcdmc_minus_model in %s\n\n", mc_minus_model_hist_file ) ; return ; }

   TH1F* h_hdp = new TH1F("h_hdp ","hdp prediction using other backgrounds' estimations", nb_global_after_exclusion, 0.5, ( (double) nb_global_after_exclusion ) + 0.5 ) ;

   ifstream ifs_data ;
   ifs_data.open( data_file ) ;
   if ( !ifs_data.good() ) { printf("\n\n *** Problem opening data file: %s\n\n", data_file ) ; return ; }


   ifstream ifs_lostlep ;
   ifs_lostlep.open( lostlep_file ) ;
   if ( !ifs_lostlep.good() ) { printf("\n\n *** Problem opening lostlep file: %s\n\n", lostlep_file ) ; return ; }

   ifstream ifs_hadtau ;
   ifs_hadtau.open( hadtau_file ) ;
   if ( !ifs_hadtau.good() ) { printf("\n\n *** Problem opening hadtau file: %s\n\n", hadtau_file ) ; return ; }

   ifstream ifs_znunu ;
   ifs_znunu.open( znunu_file ) ;
   if ( !ifs_znunu.good() ) { printf("\n\n *** Problem opening znunu file: %s\n\n", znunu_file ) ; return ; }

   FILE* ofp_combine ;
   if ( (ofp_combine = fopen( combine_output_file, "w" ))==NULL ) {
      printf( "\n\n *** Problem opening combine output file: %s\n\n", combine_output_file ) ;
      return ;
   }

    FILE* ofp_qcdratio ;
   if ( (ofp_qcdratio = fopen( qcdratio_output_file, "w" ))==NULL ) {
      printf( "\n\n *** Problem opening qcdratio output file: %s\n\n", qcdratio_output_file ) ;
      return ;
   }


   FILE* ofp_nonqcd_hdp ;
   if ( (ofp_nonqcd_hdp = fopen( "outputfiles/nonqcd-hdp-summary.txt", "w"))==NULL ) {
      printf( "\n\n *** Problem opening outputfiles/nonqcd-hdp-summary.txt output file.\n\n") ;
      gSystem -> Exit(-1) ;
   }

   FILE* ofp_debug ;
   if ( (ofp_debug = fopen( "outputfiles/gci3-debug.txt", "w"))==NULL ) {
      printf( "\n\n *** Problem opening outputfiles/gci3-debug.txt output file.\n\n") ;
      gSystem -> Exit(-1) ;
   }


   TH1F* h_ratio_all = new TH1F( "h_ratio_all", "QCD H/L ratio", nb_global_after_exclusion, 0.5, ( (double) nb_global_after_exclusion ) + 0.5 ) ;
   TH1F* h_ratio_nj[10];
   for ( int bi_njet=1; bi_njet<=nb_nj; bi_njet++ )
   {
      TString nj_str;
      nj_str.Form("%d",bi_njet);
      h_ratio_nj[bi_njet] = new TH1F( "h_ratio_nj"+nj_str, "QCD H/L ratio, njets"+nj_str, no_bin_njet[bi_njet-1], 0.5, ( (double) no_bin_njet[bi_njet-1] ) + 0.5 ) ;
   }

   TH1F* h_ratio_ht[10];

   for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ )
   {
      TString ht_str;
      ht_str.Form("%d",bi_ht);
      h_ratio_ht[bi_ht] = new TH1F( "h_ratio_ht"+ht_str, "QCD H/L ratio, HT"+ht_str, no_bin_ht[bi_ht-1], 0.5, ( (double) no_bin_ht[bi_ht-1] ) + 0.5 ) ;
   }//bi_ht


   TH1F* h_ldp_all_bg[10]; // 0: lostlep  1: hadtau, 2: znunu and 3: data
   TH1F* h_ldp_ht_bg[10][10];


   for ( int bg = 0; bg < 4; bg++)
   {
      TString bg_str;
      if ( bg == 0 ) bg_str = "lostlep";
      if ( bg == 1 ) bg_str = "hadtau";
      if ( bg == 2 ) bg_str = "znunu";
      if ( bg == 3 ) bg_str = "data";

      h_ldp_all_bg[bg] = new TH1F( "h_ldp_all_"+bg_str, "LDP, "+bg_str+", all", nb_global_after_exclusion, 0.5, ( (double) nb_global_after_exclusion )+0.5 ) ;

      for ( int bi_ht =1; bi_ht<=nBinsHT; bi_ht++ )
      {
         TString ht_str;
         ht_str.Form("%d",bi_ht);
         h_ldp_ht_bg[bi_ht][bg] = new TH1F( "h_ldp_ht"+ht_str+"_"+bg_str, "LDP, "+bg_str+", ht"+ht_str, no_bin_ht[bi_ht-1], 0.5, (double) no_bin_ht[bi_ht-1] + 0.5 ) ;
      }//bi_ht
   }//bg


      fprintf( ofp_combine, "          Search bin            N_ldp   Nnonqcd +/- err    Rqcd  Rqcd*Nqldp ");

      //////////////////for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ ) { fprintf( ofp_combine, "Kht%d   ", bi_ht); }

      //////////////////fprintf( ofp_combine, " ");
      //////////////////for ( int bi_nj = 1; bi_nj<=nb_nj; bi_nj++ ) {
      //////////////////   if ( bi_nj != njet_bin_to_be_fixed_in_qcd_model_fit+1 ) { fprintf( ofp_combine, "Snj%d   ", bi_nj); }
      //////////////////}

      fprintf( ofp_combine, " Kj1h1  Kj1h2  Kj1h3  Kj2h1  Kj2h2  Kj2h3  Kj3h1  Kj3h2  Kj3h3  Kj4h2  Kj4h3  Kj5h2  Kj5h3   " ) ;

      fprintf( ofp_combine, " Sh1m1  Sh1m2   Sh2m1  Sh2m2  Sh2m3  Sh2m4  Sh3m1  Sh3m2  Sh3m3  Sh3m4      MCC      Rqcd +/- err        Model Nqcd_hdp\n") ;

      fprintf( ofp_nonqcd_hdp, "    Label           |          lostlep            |            hadtau           |             Znunu           ||  Total non-QCD     |       QCD BG     || frac  | err ratio\n" ) ;
      fprintf( ofp_nonqcd_hdp, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" ) ;
      int bi_hist(0) ;

      int    nobs_ldp[10][10][100] ;
      int    nobs_hdp[10][10][100] ;

      int ht_plot_bi[10]{};
      int nj_plot_bi[10]{};

   // //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // printf("\n\n ######### debug ##########\n") ;
   // for ( int htmhti=4; htmhti<=13; htmhti++ ) {
   //    int hti, mhti ;
   //    htmht_bin_to_ht_and_mht_bins ( htmhti, hti, mhti );
   //    printf( " htmhti=%3d  par_val_ht_mht[%d][%d] = %7.4f\n", htmhti, hti, mhti, par_val_ht_mht[hti][mhti] ) ;
   // } // htmhti
   // printf(" ######### debug ##########\n\n\n") ;
   // //++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               int bi_ht = 0, bi_mht = 0 ;

               htmht_bin_to_ht_and_mht_bins ( bi_htmht, bi_ht, bi_mht );

               if ( is_this_bin_excluded(bi_nj - 1, bi_nb -1, bi_ht - 1, bi_mht - 1 ) ) {
                       //Because those bins we want to remove still exist in background text files. If we change input text files, we can remove these lines
                       line.ReadLine( ifs_data ) ; 
                       line.ReadLine( ifs_lostlep ) ;
                       line.ReadLine( ifs_hadtau ) ;
                       line.ReadLine( ifs_znunu ) ;
                       continue;
               }

               int global_bi, region_bi ;
               char region_tag[5] ;
               char label[100] ;
               int r_nobs_ldp, r_nobs_hdp ;
               float nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst ;
               float nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ;

               if ( bi_htmht > 3 ) {
                  ht_plot_bi[bi_ht]++;
               }
              //-------
               line.ReadLine( ifs_data ) ;
               sscanf( line.Data(), "%d %s %d %s  %d  %d", &global_bi, region_tag, &region_bi, label, &r_nobs_ldp, &r_nobs_hdp ) ;

               region_bi = gbi_search_bins[bi_nj][bi_nb][bi_htmht]; /// Because new binning is not compatible with text files we have right now.
               printf( "  Data    : %3d %s %3d  %s  Nldp = %5d ,                        Nhdp = %5d\n",
                   global_bi, region_tag, region_bi, label, r_nobs_ldp, r_nobs_hdp ) ;
               nobs_ldp[bi_nj][bi_nb][bi_htmht] = r_nobs_ldp ;
               nobs_hdp[bi_nj][bi_nb][bi_htmht] = r_nobs_hdp ;

               if ( bi_htmht > 3 ) {
                  h_ldp_all_bg[3] -> SetBinContent( region_bi, r_nobs_ldp ) ;
                  h_ldp_all_bg[3] -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                  h_ldp_ht_bg[bi_ht][3] -> SetBinContent( ht_plot_bi[bi_ht], r_nobs_ldp ) ;
                  h_ldp_ht_bg[bi_ht][3] -> GetXaxis() -> SetBinLabel( ht_plot_bi[bi_ht], label ) ;
               }

               double total_bg_ldp_val = 0. ;
               double total_bg_ldp_err2 = 0. ;

               double total_bg_hdp_val = 0. ;
               double total_bg_hdp_err2 = 0. ;

               if ( bi_htmht > 3 ) fprintf( ofp_nonqcd_hdp, " %s  ", label ) ;
               if ( bi_htmht > 3 ) fprintf( ofp_debug, " %s  ", label ) ;


              //-------
               int hdp_value = r_nobs_hdp;

               for ( int bg = 0; bg < 3; bg++) {

                  if ( bg == 0 ) line.ReadLine( ifs_lostlep ) ;
                  if ( bg == 1 ) line.ReadLine( ifs_hadtau ) ;
                  if ( bg == 2 ) line.ReadLine( ifs_znunu ) ;


                  char bg_str[10];
                  if ( bg == 0 ) strcpy(bg_str, "lostlep");
                  if ( bg == 1 ) strcpy(bg_str, "hadtau");
                  if ( bg == 2 ) strcpy(bg_str, "znunu");

                  sscanf( line.Data(), "%d %s %d %s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                       &global_bi, region_tag, &region_bi, label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;

                  region_bi = gbi_search_bins[bi_nj][bi_nb][bi_htmht]; /// Because new binning is not compatible with text files we have right now.


                  printf( " %8s : %3d %s %3d  %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                      bg_str, global_bi, region_tag, region_bi, label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

                  total_bg_ldp_val += nbg_ldp_val ;
                  total_bg_hdp_val += nbg_hdp_val ;
                  total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
                  total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

                  if ( bi_htmht > 3 ) {
                     h_ldp_all_bg[bg] -> SetBinContent( region_bi, nbg_ldp_val ) ;
                     h_ldp_all_bg[bg] -> SetBinError( region_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                     h_ldp_all_bg[bg] -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                     h_ldp_ht_bg[bi_ht][bg] -> SetBinContent( ht_plot_bi[bi_ht], nbg_ldp_val ) ;
                     h_ldp_ht_bg[bi_ht][bg] -> SetBinError( ht_plot_bi[bi_ht], sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                     h_ldp_ht_bg[bi_ht][bg] -> GetXaxis() -> SetBinLabel( ht_plot_bi[bi_ht], label ) ;

                     hdp_value -= nbg_hdp_val; 
                     h_hdp -> SetBinContent(region_bi, hdp_value);
                     h_hdp -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                     h_hdp -> GetXaxis() -> LabelsOption( "v" ) ;

                     fprintf( ofp_nonqcd_hdp, " | %6.1f +/- %5.1f +/- %5.1f ", nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

                  }//if bi_htmht

               }//bg


               double total_bg_ldp_err = sqrt( total_bg_ldp_err2 ) ;
               double total_bg_hdp_err = sqrt( total_bg_hdp_err2 ) ;

               if ( bi_htmht > 3 ) fprintf( ofp_nonqcd_hdp, " || %7.1f +/- %5.1f ", total_bg_hdp_val, total_bg_hdp_err ) ;

               printf( "total BG  : %3d %s %3d  %s  Nldp = %7.1f +/- %5.1f           ,  Nhdp = %7.1f +/- %5.1f          \n",
                   global_bi, region_tag, region_bi, label, total_bg_ldp_val, total_bg_ldp_err,   total_bg_hdp_val, total_bg_hdp_err ) ;

               double qcd_ldp_val = r_nobs_ldp - total_bg_ldp_val ;
               double qcd_ldp_err = sqrt( r_nobs_ldp + total_bg_ldp_err*total_bg_ldp_err ) ;

               if ( bi_htmht <= 3 ) {
                  double qcd_hdp_val = r_nobs_hdp - total_bg_hdp_val ;
                  double qcd_hdp_err = sqrt( r_nobs_hdp + total_bg_hdp_err*total_bg_hdp_err ) ;

                  printf( "QCD       : %3d %s %3d  %s  Nldp = %7.1f +/- %5.1f           ,  Nhdp = %7.1f +/- %5.1f          \n",
                      global_bi, region_tag, region_bi, label, qcd_ldp_val, qcd_ldp_err,   qcd_hdp_val, qcd_hdp_err ) ;
               }

               printf("\n") ;


               if ( bi_htmht > 3 ) {

                  char combine_label[100] ;

                  sprintf( combine_label, "%3d NJets%d_BTags%d_MHT%d_HT%d", region_bi, bi_nj-1, bi_nb-1, bi_mht-2, bi_ht-1 ) ;

                  char owen_label[100] ;
                  sprintf( owen_label, "%3d Nj%d-Nb%d-MHT%d-HT%d (%d)", region_bi, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht-3 ) ;

                  float Rqcd_val(0.) ;

                  /////////////////////Rqcd_val = par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_ht_mht[bi_ht][bi_mht] ;
                  Rqcd_val = par_val_ht_njet[bi_ht][bi_nj] * par_val_ht_mht[bi_ht][bi_mht] ;

                  printf("    debug2:  %s Rqcd_val = par_val_ht_njet[%d][%d] * par_val_ht_mht[%d][%d] = %7.4f * %7.4f = %7.4f\n",
                      owen_label, bi_ht, bi_nj, bi_ht, bi_mht, par_val_ht_njet[bi_ht][bi_nj], par_val_ht_mht[bi_ht][bi_mht], Rqcd_val ) ;

                  /////////////////if ( bi_htmht==4 && bi_nj>2 ) Rqcd_val = 0 ; // shut things off in the first ht bin in the highest 2 njets bins.
                  /////////////////if ( bi_htmht==7 && bi_nj>2 ) Rqcd_val = 0 ; // shut things off in the first ht bin in the highest 2 njets bins.

                  bi_hist++ ;
                  float mc_minus_model_val = h_ratio_qcdmc_minus_model -> GetBinContent( bi_hist ) ;
                  float mc_minus_model_err = h_ratio_qcdmc_minus_model -> GetBinError( bi_hist ) ;

                  //printf("    debug3:  %s , hist bin label %s h_ratio_qcdmc_minus_model = %8.4f +/- %8.4f\n", owen_label, 
                  //h_ratio_qcdmc_minus_model -> GetXaxis()->GetBinLabel( bi_hist), mc_minus_model_val, mc_minus_model_err ) ;		  

                  if ( mc_minus_model_val < 0. ) mc_minus_model_val = 0. ;
                  if ( mc_minus_model_val > 0.5 ) {
                     mc_minus_model_val = 0.5 ;
                     mc_minus_model_err = 0.5 ;
                  }
                  if ( mc_minus_model_err > 0.5 ) {
                     mc_minus_model_err = 0.5 ;
                  }

                  float deltar_over_r(0.) ;
                  if ( Rqcd_val > 0 ) deltar_over_r = mc_minus_model_val / Rqcd_val ;
                  float correction_val = 0.5 * mc_minus_model_val ;
                  float correction_err(0.) ;
                  if ( mc_minus_model_err > mc_minus_model_val ) {
                     correction_err = 0.5 * mc_minus_model_err ;
                  } else {
                     correction_err = 0.5 * mc_minus_model_val ;
                  }
                  float corrected_Rqcd_val = Rqcd_val + correction_val ;
                  float correction_rel_err(0.) ;
                  if ( corrected_Rqcd_val > 0 ) correction_rel_err = correction_val / corrected_Rqcd_val ; // do it this way to avoid huge rel error.

                  float relative_error_scale_factor(1.) ;
                  if ( correction_val > 0 && corrected_Rqcd_val > 0. ) relative_error_scale_factor = Rqcd_val / corrected_Rqcd_val ;


                  Rqcd_val = corrected_Rqcd_val ;



                  float nqcd_hdp_val = Rqcd_val * qcd_ldp_val ;

                  ////////////////std::vector<float> rel_err_ht (10,1.); // creating a [10] vector with inital value 1.
                  ////////////////rel_err_ht[bi_ht] += par_rel_err_ht[bi_ht] ;

                  ////////////////std::vector <float> rel_err_nj(10,1.);// creating a [10] vector with inital value 1.
                  ////////////////if ( bi_nj != njet_bin_to_be_fixed_in_qcd_model_fit+1 ) rel_err_nj[bi_nj] += par_rel_err_njet[bi_nj] ;

                  float rel_err_ht_njet[10][10] ;
                  for ( int lbi_ht=1; lbi_ht<=nBinsHT; lbi_ht++ ) {
                     for ( int lbi_nj=1; lbi_nj<=nb_nj; lbi_nj++ ) {
                        rel_err_ht_njet[lbi_ht][lbi_nj] = 1. ;
                     } // lbi_nj
                  } // lbi_ht
                  rel_err_ht_njet[bi_ht][bi_nj] += par_rel_err_ht_njet[bi_ht][bi_nj] * relative_error_scale_factor ;



                  std::vector<std::vector<float>> rel_err_ht_mht(10,std::vector<float>(10,1.)) ; //creating a [10][10] vector with inital value 1.
                  rel_err_ht_mht[bi_ht][bi_mht] += par_rel_err_ht_mht[bi_ht][bi_mht] * relative_error_scale_factor ;

                  std::vector <float> rel_err_nb(10,1.);
                  rel_err_nb[bi_nb-1] += par_rel_err_nb[bi_nb] * relative_error_scale_factor ;

                  float rel_err_mcclosure(1.) ;
                  rel_err_mcclosure += correction_rel_err ;


                  float Rqcd_rel_err2(0.) ;

                  //////////////// for ( int lbi_ht = 1; lbi_ht<=nBinsHT; lbi_ht++ ) {
                  ////////////////    Rqcd_rel_err2 += pow( rel_err_ht[lbi_ht]-1, 2. ) ;
                  //////////////// } // lbi_ht

                  //////////////// for ( int lbi_nj = 1; lbi_nj<=nb_nj; lbi_nj++ ) {
                  ////////////////    if ( lbi_nj != njet_bin_to_be_fixed_in_qcd_model_fit+1 ) {
                  ////////////////       Rqcd_rel_err2 += pow( rel_err_nj[lbi_nj]-1, 2. ) ;
                  ////////////////    }
                  //////////////// } // lbi_nj


                  for ( int npi=0; npi<number_of_np; npi++ ) {
                     Rqcd_rel_err2 += pow( (np_mat[ np_mat_row_index_nj_ht[bi_nj][bi_ht] ][npi])*relative_error_scale_factor, 2. ) ;
                     fprintf( ofp_debug, " %5.3f ", (np_mat[ np_mat_row_index_nj_ht[bi_nj][bi_ht] ][npi])*relative_error_scale_factor ) ;
                  } // npi

                  for ( int lbi_ht = 1; lbi_ht<=nBinsHT; lbi_ht++ ) {
                     ////////////////for ( int lbi_mht=1; lbi_mht<nb_mht; lbi_mht++ )  // *** bug
                     for ( int lbi_mht=2; lbi_mht<=nb_mht; lbi_mht++ ) {
                        //////////if ( lbi_ht == 1 && lbi_mht > 2 ) continue; //*** bug
                        if ( lbi_ht == 1 && lbi_mht > 3 ) continue;
                        Rqcd_rel_err2 += pow( rel_err_ht_mht[lbi_ht][lbi_mht]-1, 2. ) ;
                        fprintf( ofp_debug, "  %6.4f ", rel_err_ht_mht[lbi_ht][lbi_mht]-1 ) ;
                     } // lbi_mht
                  } // lbi_ht

                  for ( int lbi_nb =1; lbi_nb<=nb_nj; lbi_nb++ ) {
                     Rqcd_rel_err2 += pow( rel_err_nb[lbi_nb-1]-1, 2. ) ;
                        fprintf( ofp_debug, "  %6.4f ", rel_err_nb[lbi_nb-1]-1) ;
                  } // lbi_nb

                  Rqcd_rel_err2 += pow( rel_err_mcclosure-1, 2. ) ;
                  fprintf( ofp_debug, "  %6.4f  |  %6.4f\n", rel_err_mcclosure-1, sqrt( Rqcd_rel_err2 ) ) ;
                  float Rqcd_err = Rqcd_val * sqrt( Rqcd_rel_err2 ) ;

                  fprintf( ofp_qcdratio, " %3d %28s  Kht%d_nj%d * Smht%d = (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) = %7.4f +/- %7.4f\n",
                      region_bi, owen_label, bi_ht, bi_nj, (bi_mht-1),
                      par_val_ht_njet[bi_ht][bi_nj], par_err_ht_njet[bi_ht][bi_nj],
                      par_val_ht_mht[bi_ht][bi_mht], par_err_ht_mht[bi_ht][bi_mht],
                      Rqcd_val, Rqcd_err ) ;

                  float nqcd_hdp_err(0.) ;
                  if ( Rqcd_val > 0 && qcd_ldp_val > 0  ) nqcd_hdp_err = nqcd_hdp_val * sqrt( pow( Rqcd_err/Rqcd_val, 2. ) + pow( qcd_ldp_err/qcd_ldp_val, 2.) ) ;


                  fprintf( ofp_combine, " %28s %6d  %8.2f +/- %6.2f ",
                     combine_label, nobs_ldp[bi_nj][bi_nb][bi_htmht], total_bg_ldp_val, total_bg_ldp_err ) ;

                  fprintf( ofp_combine, " %6.4f %7.2f  ", Rqcd_val, nqcd_hdp_val ) ;

                  if ( print_dashes ) {

                     //////////////////for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ )
                     //////////////////   if ( rel_err_ht[bi_ht] == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_ht[bi_ht] ) ; }

                     //////////////////fprintf( ofp_combine, " " ) ;
                     //////////////////for ( int bi_nj =1; bi_nj<=nb_nj; bi_nj++ ) {
                     //////////////////   if ( bi_nj != njet_bin_to_be_fixed_in_qcd_model_fit+1 ) {
                     //////////////////      if ( rel_err_nj[bi_nj] == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nj[bi_nj] ) ; }
                     //////////////////   }
                     //////////////////}

                     for ( int npi=0; npi<number_of_np; npi++ ) {
                        fprintf( ofp_combine, " %6.3f", 1.+(np_mat[ np_mat_row_index_nj_ht[bi_nj][bi_ht] ][npi])*relative_error_scale_factor ) ;
                     } // npi

                     fprintf( ofp_combine, "  " ) ;
                     for ( int lbi_ht =1; lbi_ht<=nBinsHT; lbi_ht++ ) {
                        //////////for ( int lbi_mht=1; lbi_mht<nb_mht; lbi_mht++ )  //*** bug
                        for ( int lbi_mht=2; lbi_mht<=nb_mht; lbi_mht++ ) {
                           ////////if ( lbi_ht == 1 && lbi_mht > 2 ) continue;
                           if ( lbi_ht == 1 && lbi_mht > 3 ) continue; //*** bug
                           if ( rel_err_ht_mht[lbi_ht][lbi_mht] == 1. ) { fprintf( ofp_combine, "   -   " ) ; }
                           else { fprintf( ofp_combine, " %6.3f", rel_err_ht_mht[lbi_ht][lbi_mht] ) ; }
                        }
                        if ( lbi_ht != nBinsHT ) fprintf( ofp_combine, " " ) ;
                     }
                     fprintf( ofp_combine, " " ) ;

                     fprintf( ofp_combine, " " ) ;
                     if ( rel_err_mcclosure == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_mcclosure ) ; }

                  } else {

                     ///////////////////// for ( int bi_ht =1; bi_ht<=nBinsHT; bi_ht++ ) { fprintf( ofp_combine, " %6.3f", rel_err_ht[bi_ht] ) ; }
                     ///////////////////// for ( int bi_nj =1; bi_nj<=nb_nj  ; bi_nb++ ) {
                     /////////////////////    if ( bi_nj != njet_bin_to_be_fixed_in_qcd_model_fit+1) {
                     /////////////////////       fprintf( ofp_combine, " %6.3f", rel_err_nj[bi_nj] ) ;
                     /////////////////////    }
                     ///////////////////// }

                     for ( int npi=0; npi<number_of_np; npi++ ) {
                        fprintf( ofp_combine, " %6.3f", 1.+(np_mat[ np_mat_row_index_nj_ht[bi_nj][bi_ht] ][npi])*relative_error_scale_factor ) ;
                     } // npi
                     fprintf( ofp_combine, " " ) ;

                     for ( int bi_ht =1; bi_ht<=nBinsHT; bi_ht++ ) {
                        ///////////////for ( int bi_mht=1; bi_mht<nb_mht; bi_mht++ )  //*** bug
                        for ( int bi_mht=2; bi_mht<=nb_mht; bi_mht++ ) {
                           ////////////if ( bi_ht == 1 && bi_mht > 2 ) continue; //*** bug
                           if ( bi_ht == 1 && bi_mht > 3 ) continue;
                           fprintf( ofp_combine, " %6.3f", rel_err_ht_mht[bi_ht][bi_mht] ) ;
                        }
                     }
                     fprintf( ofp_combine, " %6.3f   ", rel_err_mcclosure ) ;

                  }//else print_dashes

                  fprintf( ofp_combine, "    " ) ;
                  fprintf( ofp_combine, " %6.4f +/- %6.4f , ", Rqcd_val, Rqcd_err ) ;
                  fprintf( ofp_combine, " %7.2f +/- %5.2f   ", nqcd_hdp_val, nqcd_hdp_err ) ;

                  fprintf( ofp_nonqcd_hdp, " | %5.1f +/- %5.1f ", nqcd_hdp_val, nqcd_hdp_err ) ;

                  double qcd_frac(0.) ;
                  if ( (total_bg_hdp_val+nqcd_hdp_val) > 0 && nqcd_hdp_val>0 ) {
                     qcd_frac = nqcd_hdp_val/(total_bg_hdp_val+nqcd_hdp_val) ;
                  }
                  double qcd_err_ratio(0.) ;
                  if ( total_bg_hdp_err > 0 ) {
                     qcd_err_ratio = nqcd_hdp_err / total_bg_hdp_err ;
                  }
                  fprintf( ofp_nonqcd_hdp, " || %5.1f | %5.2f |\n", 100*qcd_frac, qcd_err_ratio ) ;

                  float overall_rel_err(0.) ;
                  if ( nqcd_hdp_val > 0.1 ) {
                     overall_rel_err = nqcd_hdp_err/nqcd_hdp_val ;
                  }
                  fprintf( ofp_combine, " ( %5.0f %%) ", 100. * overall_rel_err ) ;


                  h_ratio_all -> SetBinContent( region_bi, Rqcd_val ) ;
                  h_ratio_all -> SetBinError( region_bi, Rqcd_err ) ;
                  h_ratio_all -> GetXaxis() -> SetBinLabel( region_bi, owen_label ) ;
                  nj_plot_bi[bi_nj]++;
                  h_ratio_nj[bi_nj] -> SetBinContent( nj_plot_bi[bi_nj], Rqcd_val ) ;
                  h_ratio_nj[bi_nj] -> SetBinError( nj_plot_bi[bi_nj], Rqcd_err ) ;
                  h_ratio_nj[bi_nj] -> GetXaxis() -> SetBinLabel( nj_plot_bi[bi_nj], owen_label ) ;

                  h_ratio_ht[bi_ht] -> SetBinContent( ht_plot_bi[bi_ht], Rqcd_val ) ;
                  h_ratio_ht[bi_ht] -> SetBinError( ht_plot_bi[bi_ht], Rqcd_val ) ;
                  h_ratio_ht[bi_ht] -> GetXaxis() -> SetBinLabel( ht_plot_bi[bi_ht], owen_label ) ;


                  fprintf( ofp_combine, "\n" ) ;
               }//bi_htmht

            } // bi_htmht
            fprintf( ofp_nonqcd_hdp, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" ) ;
         } // bi_nb
         fprintf( ofp_nonqcd_hdp, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" ) ;
      } // bi_nj

      fclose( ofp_combine ) ;
      printf("\n\n Wrote combine input file : %s\n\n", combine_output_file ) ;

      fclose( ofp_qcdratio ) ;
      printf("\n\n Wrote qcd ratio file : %s\n\n", qcdratio_output_file ) ;

      fclose( ofp_nonqcd_hdp ) ;

      fclose( ofp_debug ) ;


      h_ratio_all -> SetMarkerStyle(20) ;
      h_ratio_all -> GetXaxis() -> LabelsOption( "v" ) ;

      for ( int bi_nj =1; bi_nj<=nb_nj; bi_nj++ )
      {
         h_ratio_nj[bi_nj] -> SetMarkerStyle(20) ;
         h_ratio_nj[bi_nj] -> GetXaxis() -> LabelsOption( "v" ) ;
      }
      for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ )
      {
         h_ratio_ht[bi_ht] -> SetMarkerStyle(20) ;
         h_ratio_ht[bi_ht] -> GetXaxis() -> LabelsOption( "v" ) ;
      }

      for ( int bg = 0; bg < 4; bg++)
      {
         h_ldp_all_bg[bg] -> GetXaxis() -> LabelsOption( "v" ) ;
         for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ )
            h_ldp_ht_bg[bi_ht][bg] -> GetXaxis() -> LabelsOption( "v" ) ;
      }


      saveHist( hist_output_file, "h*" ) ;



    //----

      gStyle -> SetPadBottomMargin(0.35) ;

      TCanvas* can_all = new TCanvas( "can_all", "All bins", 1900, 500 ) ;
      THStack* hstack_ldp_all = new THStack( "hstack_ldp_all", "hstack_ldp_all" ) ;
      TLegend* legend = new TLegend( 0.80, 0.80, 0.95, 0.95 ) ;

      for ( int bg = 0; bg < 4; bg++)
         if ( bg != 3 ) hstack_ldp_all -> Add( h_ldp_all_bg[bg] ) ;

      h_ldp_all_bg[0] -> SetFillColor( kBlue-10 ) ;
      h_ldp_all_bg[1] -> SetFillColor( kCyan-10 ) ;
      h_ldp_all_bg[2] -> SetFillColor( kGreen-7 ) ;
      h_ldp_all_bg[3] -> SetLineWidth(3) ;

      legend -> AddEntry( h_ldp_all_bg[2], "znunu" ) ;
      legend -> AddEntry( h_ldp_all_bg[1], "hadtau" ) ;
      legend -> AddEntry( h_ldp_all_bg[0], "lostlep" ) ;
      legend -> AddEntry( h_ldp_all_bg[3], "data" ) ;


      h_ldp_all_bg[3] -> Draw() ;
      hstack_ldp_all -> Draw( "hist same" ) ;
      hstack_ldp_all -> Draw( "same" ) ;
      h_ldp_all_bg[3] -> Draw( "axis same" ) ;
      h_ldp_all_bg[3] -> Draw( "axig same" ) ;
      h_ldp_all_bg[3] -> Draw( "same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;


    //----

      TCanvas* can_ht;
      THStack* hstack_ldp_ht[10];

      for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ )
      {
         TString ht_str;
         ht_str.Form("%d",bi_ht);	   

         hstack_ldp_ht[bi_ht] = new THStack( "hstack_ldp_ht"+ht_str, "hstack_ldp_ht"+ht_str ) ;

         can_ht = new TCanvas( "can_ht"+ht_str, "HT"+ht_str+" bins", 1200, 500 ) ;

         for ( int bg = 0; bg < 4; bg++) {
            if ( bg != 3 ) {
               hstack_ldp_ht[bi_ht] -> Add( h_ldp_ht_bg[bi_ht][bg] ) ;
            }
         }

         h_ldp_ht_bg[bi_ht][0] -> SetFillColor( kBlue-10 ) ;
         h_ldp_ht_bg[bi_ht][1] -> SetFillColor( kCyan-10 ) ;
         h_ldp_ht_bg[bi_ht][2] -> SetFillColor( kGreen-7 ) ;
         h_ldp_ht_bg[bi_ht][3] -> SetLineWidth(3) ;

         h_ldp_ht_bg[bi_ht][3] -> Draw() ;
         hstack_ldp_ht[bi_ht]  -> Draw( "hist same" ) ;
         hstack_ldp_ht[bi_ht]  -> Draw( "same" ) ;
         h_ldp_ht_bg[bi_ht][3] -> Draw( "axis same" ) ;
         h_ldp_ht_bg[bi_ht][3] -> Draw( "axig same" ) ;
         h_ldp_ht_bg[bi_ht][3] -> Draw( "same" ) ;

         legend -> Draw() ;

         gPad -> SetGridy(1) ;
         gPad -> SetLogy(1) ;
      }//bi_ht

    //----


   } // gen_combine_input3


  //=======================================================================================

#endif


#include "TH1F.h"
#include "TStyle.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>

#include "histio.c"

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

      float par_rel_err_ht[5] ;
      float par_rel_err_njet[5] ;
      float par_rel_err_mht_hth[6] ;
      float par_rel_err_mht_htm[6] ;
      float par_rel_err_mht_htl[6] ;
      float par_rel_err_nb[5] ;

  //--------
   void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2, float& rel_err ) ;
   void read_pars( const char* model_pars_file ) ;
  //--------

   void gen_combine_input2(
         const char* model_pars_file = "model-pars-data4.txt",
         const char* mc_minus_model_hist_file = "outputfiles/model-ratio-hist1.root",
         const char* data_file    = "outputfiles/combine-input-data.txt",
         const char* lostlep_file = "outputfiles/combine-input-lostlep.txt",
         const char* hadtau_file  = "outputfiles/combine-input-hadtau.txt",
         const char* znunu_file   = "outputfiles/combine-input-znunu.txt",
         const char* combine_output_file = "outputfiles/combine-input-all.txt",
         const char* qcdratio_output_file = "outputfiles/qcd-ratios.txt",
         const char* hist_output_file = "outputfiles/gci-output.root"
                         ) {

      read_pars( model_pars_file ) ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      //bool print_dashes( true ) ;
      bool print_dashes( false ) ;

      TString line ;

      ////////ifstream ifs_modelfit ;
      ////////ifs_modelfit.open( modelfit_file ) ;
      ////////if ( !ifs_modelfit.good() ) { printf("\n\n *** Problem opening modelfit file: %s\n\n", modelfit_file ) ; return ; }

      TFile* tf_mc_minus_model = new TFile( mc_minus_model_hist_file, "READ" ) ;
      if ( tf_mc_minus_model == 0x0 ) { printf("\n\n *** bad mc_minus_model_hist_file file %s\n\n", mc_minus_model_hist_file ) ; return ; }
      if ( ! (tf_mc_minus_model->IsOpen()) ) { printf("\n\n *** bad mc_minus_model_hist_file file %s\n\n", mc_minus_model_hist_file ) ; return ; }
      TH1F* h_ratio_qcdmc_minus_model = (TH1F*) tf_mc_minus_model -> Get( "h_ratio_qcdmc_minus_model" ) ;
      if ( h_ratio_qcdmc_minus_model == 0x0 ) { printf("\n\n *** Can't find h_ratio_qcdmc_minus_model in %s\n\n", mc_minus_model_hist_file ) ; return ; }

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

  ////double kqcd_ht_val[4] ;
  ////double kqcd_ht_err[4] ;
  ////double kqcd_ht_rel_err[4] ;

  ////double sqcd_nj_val[5] ;
  ////double sqcd_nj_err[5] ;
  ////double sqcd_nj_rel_err[5] ;

  ////double sqcd_mht_val[5] ;
  ////double sqcd_mht_err[5] ;
  ////double sqcd_mht_rel_err[5] ;

  ////double sqcd_nb_val[5] ;
  ////double sqcd_nb_err[5] ;
  ////double sqcd_nb_rel_err[5] ;

      int nb_ht(3) ;
      int nb_nj(4) ;
      int nb_nb(4) ;
      int nb_htmht(13) ;


  ////for ( int bi_ht=1; bi_ht<=nb_ht; bi_ht++ ) {

  ////   char label[100] ;
  ////   float val, err, rel_err ;

  ////   line.ReadLine( ifs_modelfit ) ;
  ////   sscanf( line.Data(), "%s %f +/- %f", label, &val, &err ) ;
  ////   printf("   kqcd, HT%d : %s  %7.5f +/- %7.5f\n", bi_ht, label, val, err ) ;
  ////   rel_err = 1. ;
  ////   if ( val != 0 ) rel_err = err / val ;
  ////   kqcd_ht_val[bi_ht] = val ;
  ////   kqcd_ht_err[bi_ht] = err ;
  ////   kqcd_ht_rel_err[bi_ht] = rel_err ;

  ////}

  ////sqcd_nj_val[1] = 1. ; sqcd_nj_err[1] = 0. ; sqcd_nj_rel_err[1] = 0. ;

  ////for ( int bi_nj=2; bi_nj<=nb_nj; bi_nj++ ) {

  ////   char label[100] ;
  ////   float val, err, rel_err ;

  ////   line.ReadLine( ifs_modelfit ) ;
  ////   sscanf( line.Data(), "%s %f +/- %f", label, &val, &err ) ;
  ////   printf("   sqcd, Nj%d : %s  %7.5f +/- %7.5f\n", bi_nj, label, val, err ) ;
  ////   rel_err = 1. ;
  ////   if ( val != 0 ) rel_err = err / val ;
  ////   sqcd_nj_val[bi_nj] = val ;
  ////   sqcd_nj_err[bi_nj] = err ;
  ////   sqcd_nj_rel_err[bi_nj] = rel_err ;

  ////}

  ////sqcd_mht_val[1] = 0.65 ; sqcd_mht_err[1] = 0.195 ; sqcd_mht_rel_err[1] = 0.30 ;
  ////sqcd_mht_val[2] = 0.50 ; sqcd_mht_err[2] = 0.300 ; sqcd_mht_rel_err[2] = 0.60 ;
  ////sqcd_mht_val[3] = 0.20 ; sqcd_mht_err[3] = 0.240 ; sqcd_mht_rel_err[3] = 1.20 ;
  ////sqcd_mht_val[4] = 0.10 ; sqcd_mht_err[4] = 0.300 ; sqcd_mht_rel_err[4] = 3.00 ;

  ////sqcd_nb_val[1] = 1.00 ; sqcd_nb_err[1] = 0.100 ; sqcd_nb_rel_err[1] = 0.10 ;
  ////sqcd_nb_val[2] = 1.00 ; sqcd_nb_err[2] = 0.200 ; sqcd_nb_rel_err[2] = 0.20 ;
  ////sqcd_nb_val[3] = 1.00 ; sqcd_nb_err[3] = 0.400 ; sqcd_nb_rel_err[3] = 0.40 ;
  ////sqcd_nb_val[4] = 1.00 ; sqcd_nb_err[4] = 0.800 ; sqcd_nb_rel_err[4] = 0.80 ;




      TH1F* h_ratio_all = new TH1F( "h_ratio_all", "QCD H/L ratio", 160, 0.5, 160+0.5 ) ;

      TH1F* h_ratio_nj1 = new TH1F( "h_ratio_nj1", "QCD H/L ratio, njets1", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ratio_nj2 = new TH1F( "h_ratio_nj2", "QCD H/L ratio, njets2", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ratio_nj3 = new TH1F( "h_ratio_nj3", "QCD H/L ratio, njets3", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ratio_nj4 = new TH1F( "h_ratio_nj4", "QCD H/L ratio, njets4", 40, 0.5, 40+0.5 ) ;

      TH1F* h_ratio_ht1 = new TH1F( "h_ratio_ht1", "QCD H/L ratio, HT1", 32, 0.5, 32+0.5 ) ;
      TH1F* h_ratio_ht2 = new TH1F( "h_ratio_ht2", "QCD H/L ratio, HT2", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ratio_ht3 = new TH1F( "h_ratio_ht3", "QCD H/L ratio, HT3", 64, 0.5, 64+0.5 ) ;


      TH1F* h_ldp_all_lostlep = new TH1F( "h_ldp_all_lostlep", "LDP, lostlep, all", 160, 0.5, 160+0.5 ) ;
      TH1F* h_ldp_all_hadtau = new TH1F( "h_ldp_all_hadtau", "LDP, hadtau, all", 160, 0.5, 160+0.5 ) ;
      TH1F* h_ldp_all_znunu = new TH1F( "h_ldp_all_znunu", "LDP, znunu, all", 160, 0.5, 160+0.5 ) ;
      TH1F* h_ldp_all_data = new TH1F( "h_ldp_all_data", "LDP, data, all", 160, 0.5, 160+0.5 ) ;

      TH1F* h_ldp_nj1_lostlep = new TH1F( "h_ldp_nj1_lostlep", "LDP, lostlep, nj1", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj1_hadtau = new TH1F( "h_ldp_nj1_hadtau", "LDP, hadtau, nj1", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj1_znunu = new TH1F( "h_ldp_nj1_znunu", "LDP, znunu, nj1", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj1_data = new TH1F( "h_ldp_nj1_data", "LDP, data, nj1", 40, 0.5, 40+0.5 ) ;

      TH1F* h_ldp_nj2_lostlep = new TH1F( "h_ldp_nj2_lostlep", "LDP, lostlep, nj2", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj2_hadtau = new TH1F( "h_ldp_nj2_hadtau", "LDP, hadtau, nj2", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj2_znunu = new TH1F( "h_ldp_nj2_znunu", "LDP, znunu, nj2", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj2_data = new TH1F( "h_ldp_nj2_data", "LDP, data, nj2", 40, 0.5, 40+0.5 ) ;

      TH1F* h_ldp_nj3_lostlep = new TH1F( "h_ldp_nj3_lostlep", "LDP, lostlep, nj3", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj3_hadtau = new TH1F( "h_ldp_nj3_hadtau", "LDP, hadtau, nj3", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj3_znunu = new TH1F( "h_ldp_nj3_znunu", "LDP, znunu, nj3", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj3_data = new TH1F( "h_ldp_nj3_data", "LDP, data, nj3", 40, 0.5, 40+0.5 ) ;

      TH1F* h_ldp_nj4_lostlep = new TH1F( "h_ldp_nj4_lostlep", "LDP, lostlep, nj4", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj4_hadtau = new TH1F( "h_ldp_nj4_hadtau", "LDP, hadtau, nj4", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj4_znunu = new TH1F( "h_ldp_nj4_znunu", "LDP, znunu, nj4", 40, 0.5, 40+0.5 ) ;
      TH1F* h_ldp_nj4_data = new TH1F( "h_ldp_nj4_data", "LDP, data, nj4", 40, 0.5, 40+0.5 ) ;


      TH1F* h_ldp_ht1_lostlep = new TH1F( "h_ldp_ht1_lostlep", "LDP, lostlep, ht1", 32, 0.5, 32+0.5 ) ;
      TH1F* h_ldp_ht1_hadtau = new TH1F( "h_ldp_ht1_hadtau", "LDP, hadtau, ht1", 32, 0.5, 32+0.5 ) ;
      TH1F* h_ldp_ht1_znunu = new TH1F( "h_ldp_ht1_znunu", "LDP, znunu, ht1", 32, 0.5, 32+0.5 ) ;
      TH1F* h_ldp_ht1_data = new TH1F( "h_ldp_ht1_data", "LDP, data, ht1", 32, 0.5, 32+0.5 ) ;

      TH1F* h_ldp_ht2_lostlep = new TH1F( "h_ldp_ht2_lostlep", "LDP, lostlep, ht2", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ldp_ht2_hadtau = new TH1F( "h_ldp_ht2_hadtau", "LDP, hadtau, ht2", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ldp_ht2_znunu = new TH1F( "h_ldp_ht2_znunu", "LDP, znunu, ht2", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ldp_ht2_data = new TH1F( "h_ldp_ht2_data", "LDP, data, ht2", 64, 0.5, 64+0.5 ) ;

      TH1F* h_ldp_ht3_lostlep = new TH1F( "h_ldp_ht3_lostlep", "LDP, lostlep, ht3", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ldp_ht3_hadtau = new TH1F( "h_ldp_ht3_hadtau", "LDP, hadtau, ht3", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ldp_ht3_znunu = new TH1F( "h_ldp_ht3_znunu", "LDP, znunu, ht3", 64, 0.5, 64+0.5 ) ;
      TH1F* h_ldp_ht3_data = new TH1F( "h_ldp_ht3_data", "LDP, data, ht3", 64, 0.5, 64+0.5 ) ;


      int    nobs_ldp[5][5][14] ;
      double nonqcd_val_ldp[5][5][14] ;
      double nonqcd_err_ldp[5][5][14] ;

      int    nobs_hdp[5][5][14] ;
      double nonqcd_val_hdp[5][5][14] ;
      double nonqcd_err_hdp[5][5][14] ;

      //fprintf( ofp_combine, "          Search bin            N_ldp   Nnonqcd +/- err    Rqcd  Rqcd*Nqldp Kht1   Kht2   Kht3    Snj2   Snj3   Snj4    Smht1  Smht2  Smht3  Smht4   Snb0   Snb1   Snb3   Snb3     Rqcd +/- err        Model Nqcd_hdp\n") ;
      fprintf( ofp_combine, "          Search bin            N_ldp   Nnonqcd +/- err    Rqcd  Rqcd*Nqldp Kht1   Kht2   Kht3    Snj2   Snj3   Snj4    Sh1m1  Sh1m2   Sh2m1  Sh2m2  Sh2m3  Sh2m4  Sh3m1  Sh3m2  Sh3m3  Sh3m4      MCC      Rqcd +/- err        Model Nqcd_hdp\n") ;

      int bi_hist(0) ;

      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

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

               if ( bi_htmht ==10 ) { bi_ht = 2; bi_mht = 4; }
               if ( bi_htmht ==11 ) { bi_ht = 3; bi_mht = 4; }

               if ( bi_htmht ==12 ) { bi_ht = 2; bi_mht = 5; }
               if ( bi_htmht ==13 ) { bi_ht = 3; bi_mht = 5; }

               int global_bi, region_bi ;
               char region_tag[5] ;
               char label[100] ;
               int r_nobs_ldp, r_nobs_hdp ;
               float nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst ;
               float nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ;

               int ht3_plot_bi(-1) ;
               int ht2_plot_bi(-1) ;
               int ht1_plot_bi(-1) ;
               if ( bi_htmht > 3 ) {
                  ht3_plot_bi = (bi_nj-1)*(nb_nb)*4 + (bi_nb-1)*4 + bi_mht-1 ;
                  ht2_plot_bi = (bi_nj-1)*(nb_nb)*4 + (bi_nb-1)*4 + bi_mht-1 ;
                  ht1_plot_bi = (bi_nj-1)*(nb_nb)*2 + (bi_nb-1)*2 + bi_mht-1 ;
               }

              //-------
               line.ReadLine( ifs_data ) ;
               sscanf( line.Data(), "%d %s %d %s  %d  %d", &global_bi, region_tag, &region_bi, label, &r_nobs_ldp, &r_nobs_hdp ) ;
               printf( "  Data    : %3d %s %3d  %s  Nldp = %5d ,                        Nhdp = %5d\n",
                   global_bi, region_tag, region_bi, label, r_nobs_ldp, r_nobs_hdp ) ;

               nobs_ldp[bi_nj][bi_nb][bi_htmht] = r_nobs_ldp ;
               nobs_hdp[bi_nj][bi_nb][bi_htmht] = r_nobs_hdp ;

               if ( bi_htmht > 3 ) {
                  h_ldp_all_data -> SetBinContent( region_bi, r_nobs_ldp ) ;
                  h_ldp_all_data -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_data -> SetBinContent( ht3_plot_bi, r_nobs_ldp ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_data -> GetXaxis() -> SetBinLabel( ht3_plot_bi, label ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_data -> SetBinContent( ht2_plot_bi, r_nobs_ldp ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_data -> GetXaxis() -> SetBinLabel( ht2_plot_bi, label ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_data -> SetBinContent( ht1_plot_bi, r_nobs_ldp ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_data -> GetXaxis() -> SetBinLabel( ht1_plot_bi, label ) ;
               }


               double total_bg_ldp_val = 0. ;
               double total_bg_ldp_err2 = 0. ;

               double total_bg_hdp_val = 0. ;
               double total_bg_hdp_err2 = 0. ;

              //-------
               line.ReadLine( ifs_lostlep ) ;
               sscanf( line.Data(), "%d %s %d %s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                    &global_bi, region_tag, &region_bi, label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;
               printf( "  Lostlep : %3d %s %3d  %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                   global_bi, region_tag, region_bi, label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

               total_bg_ldp_val += nbg_ldp_val ;
               total_bg_hdp_val += nbg_hdp_val ;
               total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
               total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

               if ( bi_htmht > 3 ) {
                  h_ldp_all_lostlep -> SetBinContent( region_bi, nbg_ldp_val ) ;
                  h_ldp_all_lostlep -> SetBinError( region_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  h_ldp_all_lostlep -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_lostlep -> SetBinContent( ht3_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_lostlep -> SetBinError( ht3_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_lostlep -> GetXaxis() -> SetBinLabel( ht3_plot_bi, label ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_lostlep -> SetBinContent( ht2_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_lostlep -> SetBinError( ht2_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_lostlep -> GetXaxis() -> SetBinLabel( ht2_plot_bi, label ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_lostlep -> SetBinContent( ht1_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_lostlep -> SetBinError( ht1_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_lostlep -> GetXaxis() -> SetBinLabel( ht1_plot_bi, label ) ;
               }


              //-------
               line.ReadLine( ifs_hadtau ) ;
               sscanf( line.Data(), "%d %s %d %s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                    &global_bi, region_tag, &region_bi, label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;
               printf( "  hadtau  : %3d %s %3d  %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                   global_bi, region_tag, region_bi, label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

               total_bg_ldp_val += nbg_ldp_val ;
               total_bg_hdp_val += nbg_hdp_val ;
               total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
               total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

               if ( bi_htmht > 3 ) {
                  h_ldp_all_hadtau -> SetBinContent( region_bi, nbg_ldp_val ) ;
                  h_ldp_all_hadtau -> SetBinError( region_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  h_ldp_all_hadtau -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_hadtau -> SetBinContent( ht3_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_hadtau -> SetBinError( ht3_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_hadtau -> GetXaxis() -> SetBinLabel( ht3_plot_bi, label ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_hadtau -> SetBinContent( ht2_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_hadtau -> SetBinError( ht2_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_hadtau -> GetXaxis() -> SetBinLabel( ht2_plot_bi, label ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_hadtau -> SetBinContent( ht1_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_hadtau -> SetBinError( ht1_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_hadtau -> GetXaxis() -> SetBinLabel( ht1_plot_bi, label ) ;
               }

              //-------
               line.ReadLine( ifs_znunu ) ;
               sscanf( line.Data(), "%d %s %d %s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                    &global_bi, region_tag, &region_bi, label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;
               printf( "  znunu   : %3d %s %3d  %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                   global_bi, region_tag, region_bi, label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

               total_bg_ldp_val += nbg_ldp_val ;
               total_bg_hdp_val += nbg_hdp_val ;
               total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
               total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

               if ( bi_htmht > 3 ) {
                  h_ldp_all_znunu -> SetBinContent( region_bi, nbg_ldp_val ) ;
                  h_ldp_all_znunu -> SetBinError( region_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  h_ldp_all_znunu -> GetXaxis() -> SetBinLabel( region_bi, label ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_znunu -> SetBinContent( ht3_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_znunu -> SetBinError( ht3_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 3 ) h_ldp_ht3_znunu -> GetXaxis() -> SetBinLabel( ht3_plot_bi, label ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_znunu -> SetBinContent( ht2_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_znunu -> SetBinError( ht2_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 2 ) h_ldp_ht2_znunu -> GetXaxis() -> SetBinLabel( ht2_plot_bi, label ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_znunu -> SetBinContent( ht1_plot_bi, nbg_ldp_val ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_znunu -> SetBinError( ht1_plot_bi, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
                  if ( bi_ht == 1 ) h_ldp_ht1_znunu -> GetXaxis() -> SetBinLabel( ht1_plot_bi, label ) ;
               }


               double total_bg_ldp_err = sqrt( total_bg_ldp_err2 ) ;
               double total_bg_hdp_err = sqrt( total_bg_hdp_err2 ) ;

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
                  sprintf( combine_label, "%3d NJets%d-BTags%d-MHT%d-HT%d", region_bi, bi_nj-1, bi_nb-1, bi_mht-2, bi_htmht-4 ) ;

                  char owen_label[100] ;
                  sprintf( owen_label, "%3d Nj%d-Nb%d-MHT%d-HT%d (%d)", region_bi, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht-3 ) ;

                  ////////////float Rqcd_val = kqcd_ht_val[bi_ht] * sqcd_nj_val[bi_nj] * sqcd_mht_val[bi_mht-1] ;
                  float Rqcd_val(0.) ;
                  if ( bi_ht == 1 ) {
                     Rqcd_val = par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_mht_htl[bi_mht-1] ;
                  } else if ( bi_ht == 2 ) {
                     Rqcd_val = par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_mht_htm[bi_mht-1] ;
                  } else if ( bi_ht == 3 ) {
                     Rqcd_val = par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_mht_hth[bi_mht-1] ;
                  }

                  if ( bi_htmht==4 && bi_nj>2 ) Rqcd_val = 0 ; // shut things off in the first ht bin in the highest 2 njets bins.
                  if ( bi_htmht==7 && bi_nj>2 ) Rqcd_val = 0 ; // shut things off in the first ht bin in the highest 2 njets bins.

                  bi_hist++ ;
                  float mc_minus_model_val = h_ratio_qcdmc_minus_model -> GetBinContent( bi_hist ) ;
                  float mc_minus_model_err = h_ratio_qcdmc_minus_model -> GetBinError( bi_hist ) ;
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

                  Rqcd_val = corrected_Rqcd_val ;



                  ////////printf("  check:  %3d %35s  Model R = %6.3f ,  MC-model = %6.3f +/- %6.3f     corrected R = %6.3f, correction = %6.3f +/- %6.3f (%6.3f)\n", bi_hist, h_ratio_qcdmc_minus_model->GetXaxis()->GetBinLabel( bi_hist ),
                      ////////Rqcd_val, mc_minus_model_val, mc_minus_model_err,
                      ////////corrected_Rqcd_val, correction_val, correction_err, correction_rel_err ) ;






                  float nqcd_hdp_val = Rqcd_val * qcd_ldp_val ;

                  float rel_err_ht1(1.), rel_err_ht2(1.), rel_err_ht3(1.) ;
                  ///////if ( bi_ht == 1 ) rel_err_ht1 += kqcd_ht_rel_err[1] ;
                  ///////if ( bi_ht == 2 ) rel_err_ht2 += kqcd_ht_rel_err[2] ;
                  ///////if ( bi_ht == 3 ) rel_err_ht3 += kqcd_ht_rel_err[3] ;
                  if ( bi_ht == 1 ) rel_err_ht1 += par_rel_err_ht[1] ;
                  if ( bi_ht == 2 ) rel_err_ht2 += par_rel_err_ht[2] ;
                  if ( bi_ht == 3 ) rel_err_ht3 += par_rel_err_ht[3] ;

                  float rel_err_nj2(1.), rel_err_nj3(1.), rel_err_nj4(1.) ;
                  //////if ( bi_nj == 2 ) rel_err_nj2 += sqcd_nj_rel_err[2] ;
                  //////if ( bi_nj == 3 ) rel_err_nj3 += sqcd_nj_rel_err[3] ;
                  //////if ( bi_nj == 4 ) rel_err_nj4 += sqcd_nj_rel_err[4] ;
                  if ( bi_nj == 2 ) rel_err_nj2 += par_rel_err_njet[2] ;
                  if ( bi_nj == 3 ) rel_err_nj3 += par_rel_err_njet[3] ;
                  if ( bi_nj == 4 ) rel_err_nj4 += par_rel_err_njet[4] ;

                  ///////float rel_err_mht1(1.), rel_err_mht2(1.), rel_err_mht3(1.), rel_err_mht4(1.) ;
                  ///////if ( bi_mht == 2 ) rel_err_mht1 += sqcd_mht_rel_err[1] ;
                  ///////if ( bi_mht == 3 ) rel_err_mht2 += sqcd_mht_rel_err[2] ;
                  ///////if ( bi_mht == 4 ) rel_err_mht3 += sqcd_mht_rel_err[3] ;
                  ///////if ( bi_mht == 5 ) rel_err_mht4 += sqcd_mht_rel_err[4] ;
                  float rel_err_htl_mht1(1.), rel_err_htl_mht2(1.) ;
                  if ( bi_ht == 1 && bi_mht == 2 ) rel_err_htl_mht1 += par_rel_err_mht_htl[2] ;
                  if ( bi_ht == 1 && bi_mht == 3 ) rel_err_htl_mht2 += par_rel_err_mht_htl[3] ;
                  //--
                  float rel_err_htm_mht1(1.), rel_err_htm_mht2(1.), rel_err_htm_mht3(1.), rel_err_htm_mht4(1.) ;
                  if ( bi_ht == 2 && bi_mht == 2 ) rel_err_htm_mht1 += par_rel_err_mht_htm[2] ;
                  if ( bi_ht == 2 && bi_mht == 3 ) rel_err_htm_mht2 += par_rel_err_mht_htm[3] ;
                  if ( bi_ht == 2 && bi_mht == 4 ) rel_err_htm_mht3 += par_rel_err_mht_htm[4] ;
                  if ( bi_ht == 2 && bi_mht == 5 ) rel_err_htm_mht4 += par_rel_err_mht_htm[5] ;
                  //--
                  float rel_err_hth_mht1(1.), rel_err_hth_mht2(1.), rel_err_hth_mht3(1.), rel_err_hth_mht4(1.) ;
                  if ( bi_ht == 3 && bi_mht == 2 ) rel_err_hth_mht1 += par_rel_err_mht_hth[2] ;
                  if ( bi_ht == 3 && bi_mht == 3 ) rel_err_hth_mht2 += par_rel_err_mht_hth[3] ;
                  if ( bi_ht == 3 && bi_mht == 4 ) rel_err_hth_mht3 += par_rel_err_mht_hth[4] ;
                  if ( bi_ht == 3 && bi_mht == 5 ) rel_err_hth_mht4 += par_rel_err_mht_hth[5] ;

                  float rel_err_nb0(1.), rel_err_nb1(1.), rel_err_nb2(1.), rel_err_nb3(1.) ;
                  ///////if ( bi_nb == 1 ) rel_err_nb0 += sqcd_nb_rel_err[1] ;
                  ///////if ( bi_nb == 2 ) rel_err_nb1 += sqcd_nb_rel_err[2] ;
                  ///////if ( bi_nb == 3 ) rel_err_nb2 += sqcd_nb_rel_err[3] ;
                  ///////if ( bi_nb == 4 ) rel_err_nb3 += sqcd_nb_rel_err[4] ;
                  if ( bi_nb == 1 ) rel_err_nb0 += par_rel_err_nb[1] ;
                  if ( bi_nb == 2 ) rel_err_nb1 += par_rel_err_nb[2] ;
                  if ( bi_nb == 3 ) rel_err_nb2 += par_rel_err_nb[3] ;
                  if ( bi_nb == 4 ) rel_err_nb3 += par_rel_err_nb[4] ;

                  float rel_err_mcclosure(1.) ;
                  rel_err_mcclosure += correction_rel_err ;


                  float Rqcd_rel_err2(0.) ;
                  Rqcd_rel_err2 += pow( rel_err_ht1-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_ht2-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_ht3-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nj2-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nj3-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nj4-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_htl_mht1-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_htl_mht2-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_htm_mht1-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_htm_mht2-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_htm_mht3-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_htm_mht4-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_hth_mht1-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_hth_mht2-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_hth_mht3-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_hth_mht4-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nb0-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nb1-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nb2-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_nb3-1, 2. ) ;
                  Rqcd_rel_err2 += pow( rel_err_mcclosure-1, 2. ) ;
                  float Rqcd_err = Rqcd_val * sqrt( Rqcd_rel_err2 ) ;

                  if ( bi_ht == 1 ) {
                     fprintf( ofp_qcdratio, " %3d %28s  Kht%d * Snj%d * Smht%d = (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) = %7.4f +/- %7.4f\n",
                         region_bi, owen_label, bi_ht, bi_nj, (bi_mht-1),
                         par_val_ht[bi_ht], sqrt( pow( par_err_ht_fit[bi_ht], 2.) + pow( par_val_ht[bi_ht]*par_err_ht_syst[bi_ht], 2.) ),
                         par_val_njet[bi_nj], sqrt( pow( par_err_njet_fit[bi_nj], 2.) + pow( par_val_njet[bi_nj]*par_err_njet_syst[bi_nj], 2.) ),
                         par_val_mht_htl[bi_mht-1], par_err_mht_htl[bi_mht-1],
                         Rqcd_val, Rqcd_err ) ;
                  } else if ( bi_ht == 2 ) {
                     fprintf( ofp_qcdratio, " %3d %28s  Kht%d * Snj%d * Smht%d = (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) = %7.4f +/- %7.4f\n",
                         region_bi, owen_label, bi_ht, bi_nj, (bi_mht-1),
                         par_val_ht[bi_ht], sqrt( pow( par_err_ht_fit[bi_ht], 2.) + pow( par_val_ht[bi_ht]*par_err_ht_syst[bi_ht], 2.) ),
                         par_val_njet[bi_nj], sqrt( pow( par_err_njet_fit[bi_nj], 2.) + pow( par_val_njet[bi_nj]*par_err_njet_syst[bi_nj], 2.) ),
                         par_val_mht_htm[bi_mht-1], par_err_mht_htm[bi_mht-1],
                         Rqcd_val, Rqcd_err ) ;
                  } else if ( bi_ht == 3 ) {
                     fprintf( ofp_qcdratio, " %3d %28s  Kht%d * Snj%d * Smht%d = (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) * (%6.4f +/- %6.4f) = %7.4f +/- %7.4f\n",
                         region_bi, owen_label, bi_ht, bi_nj, (bi_mht-1),
                         par_val_ht[bi_ht], sqrt( pow( par_err_ht_fit[bi_ht], 2.) + pow( par_val_ht[bi_ht]*par_err_ht_syst[bi_ht], 2.) ),
                         par_val_njet[bi_nj], sqrt( pow( par_err_njet_fit[bi_nj], 2.) + pow( par_val_njet[bi_nj]*par_err_njet_syst[bi_nj], 2.) ),
                         par_val_mht_hth[bi_mht-1], par_err_mht_hth[bi_mht-1],
                         Rqcd_val, Rqcd_err ) ;
                  }

                  float nqcd_hdp_err(0.) ;
                  if ( Rqcd_val > 0 && qcd_ldp_val > 0  ) nqcd_hdp_err = nqcd_hdp_val * sqrt( pow( Rqcd_err/Rqcd_val, 2. ) + pow( qcd_ldp_err/qcd_ldp_val, 2.) ) ;


                  fprintf( ofp_combine, " %28s %6d  %8.2f +/- %6.2f ",
                     combine_label, nobs_ldp[bi_nj][bi_nb][bi_htmht], total_bg_ldp_val, total_bg_ldp_err ) ;

                  fprintf( ofp_combine, " %6.4f %7.2f ", Rqcd_val, nqcd_hdp_val ) ;

                  if ( print_dashes ) {

                     if ( rel_err_ht1 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_ht1 ) ; }
                     if ( rel_err_ht2 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_ht2 ) ; }
                     if ( rel_err_ht3 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_ht3 ) ; }
                     fprintf( ofp_combine, " " ) ;

                     if ( rel_err_nj2 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nj2 ) ; }
                     if ( rel_err_nj3 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nj3 ) ; }
                     if ( rel_err_nj4 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nj4 ) ; }
                     fprintf( ofp_combine, " " ) ;

                     if ( rel_err_htl_mht1 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_htl_mht1 ) ; }
                     if ( rel_err_htl_mht2 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_htl_mht2 ) ; }
                     fprintf( ofp_combine, " " ) ;

                     if ( rel_err_htm_mht1 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_htm_mht1 ) ; }
                     if ( rel_err_htm_mht2 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_htm_mht2 ) ; }
                     if ( rel_err_htm_mht3 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_htm_mht3 ) ; }
                     if ( rel_err_htm_mht4 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_htm_mht4 ) ; }
                     fprintf( ofp_combine, " " ) ;

                     if ( rel_err_hth_mht1 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_hth_mht1 ) ; }
                     if ( rel_err_hth_mht2 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_hth_mht2 ) ; }
                     if ( rel_err_hth_mht3 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_hth_mht3 ) ; }
                     if ( rel_err_hth_mht4 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_hth_mht4 ) ; }
                     fprintf( ofp_combine, " " ) ;

                     //////if ( rel_err_nb0 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nb0 ) ; }
                     //////if ( rel_err_nb1 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nb1 ) ; }
                     //////if ( rel_err_nb2 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nb2 ) ; }
                     //////if ( rel_err_nb3 == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_nb3 ) ; }
                     //////fprintf( ofp_combine, " " ) ;

                     fprintf( ofp_combine, " " ) ;
                     if ( rel_err_mcclosure == 1. ) { fprintf( ofp_combine, "   -   " ) ; } else { fprintf( ofp_combine, " %6.3f", rel_err_mcclosure ) ; }

                  } else {

                     fprintf( ofp_combine, " %6.3f %6.3f %6.3f ", rel_err_ht1, rel_err_ht2, rel_err_ht3 ) ;
                     fprintf( ofp_combine, " %6.3f %6.3f %6.3f ", rel_err_nj2, rel_err_nj3, rel_err_nj4 ) ;
                     fprintf( ofp_combine, " %6.3f %6.3f ", rel_err_htl_mht1, rel_err_htl_mht2 ) ;
                     fprintf( ofp_combine, " %6.3f %6.3f %6.3f %6.3f ", rel_err_htm_mht1, rel_err_htm_mht2, rel_err_htm_mht3, rel_err_htm_mht4 ) ;
                     fprintf( ofp_combine, " %6.3f %6.3f %6.3f %6.3f ", rel_err_hth_mht1, rel_err_hth_mht2, rel_err_hth_mht3, rel_err_hth_mht4 ) ;
                     ///////fprintf( ofp_combine, " %6.3f %6.3f %6.3f %6.3f ", rel_err_nb0, rel_err_nb1, rel_err_nb2, rel_err_nb3 ) ;
                     fprintf( ofp_combine, " %6.3f   ", rel_err_mcclosure ) ;

                  }
                  fprintf( ofp_combine, "  " ) ;
                  fprintf( ofp_combine, " %6.4f +/- %6.4f , ", Rqcd_val, Rqcd_err ) ;
                  fprintf( ofp_combine, " %7.2f +/- %5.2f   ", nqcd_hdp_val, nqcd_hdp_err ) ;




                  h_ratio_all -> SetBinContent( region_bi, Rqcd_val ) ;
                  h_ratio_all -> SetBinError( region_bi, Rqcd_err ) ;
                  h_ratio_all -> GetXaxis() -> SetBinLabel( region_bi, owen_label ) ;

                  int ht1_plot_bi = (bi_nj-1)*(nb_nb)*2 + (bi_nb-1)*2 + bi_mht-1 ;
                  int ht2_plot_bi = (bi_nj-1)*(nb_nb)*4 + (bi_nb-1)*4 + bi_mht-1 ;
                  int ht3_plot_bi = (bi_nj-1)*(nb_nb)*4 + (bi_nb-1)*4 + bi_mht-1 ;
                  int nj_plot_bi = (bi_nb-1)*10 + bi_htmht-3 ;

                  if ( bi_nj==1 ) h_ratio_nj1 -> SetBinContent( nj_plot_bi, Rqcd_val ) ;
                  if ( bi_nj==1 ) h_ratio_nj1 -> SetBinError( nj_plot_bi, Rqcd_err ) ;
                  if ( bi_nj==1 ) h_ratio_nj1 -> GetXaxis() -> SetBinLabel( nj_plot_bi, owen_label ) ;

                  if ( bi_nj==2 ) h_ratio_nj2 -> SetBinContent( nj_plot_bi, Rqcd_val ) ;
                  if ( bi_nj==2 ) h_ratio_nj2 -> SetBinError( nj_plot_bi, Rqcd_err ) ;
                  if ( bi_nj==2 ) h_ratio_nj2 -> GetXaxis() -> SetBinLabel( nj_plot_bi, owen_label ) ;

                  if ( bi_nj==3 ) h_ratio_nj3 -> SetBinContent( nj_plot_bi, Rqcd_val ) ;
                  if ( bi_nj==3 ) h_ratio_nj3 -> SetBinError( nj_plot_bi, Rqcd_err ) ;
                  if ( bi_nj==3 ) h_ratio_nj3 -> GetXaxis() -> SetBinLabel( nj_plot_bi, owen_label ) ;

                  if ( bi_nj==4 ) h_ratio_nj4 -> SetBinContent( nj_plot_bi, Rqcd_val ) ;
                  if ( bi_nj==4 ) h_ratio_nj4 -> SetBinError( nj_plot_bi, Rqcd_err ) ;
                  if ( bi_nj==4 ) h_ratio_nj4 -> GetXaxis() -> SetBinLabel( nj_plot_bi, owen_label ) ;

                  if ( bi_ht == 1 ) h_ratio_ht1 -> SetBinContent( ht1_plot_bi, Rqcd_val ) ;
                  if ( bi_ht == 1 ) h_ratio_ht1 -> SetBinError( ht1_plot_bi, Rqcd_val ) ;
                  if ( bi_ht == 1 ) h_ratio_ht1 -> GetXaxis() -> SetBinLabel( ht1_plot_bi, owen_label ) ;

                  if ( bi_ht == 2 ) h_ratio_ht2 -> SetBinContent( ht2_plot_bi, Rqcd_val ) ;
                  if ( bi_ht == 2 ) h_ratio_ht2 -> SetBinError( ht2_plot_bi, Rqcd_val ) ;
                  if ( bi_ht == 2 ) h_ratio_ht2 -> GetXaxis() -> SetBinLabel( ht2_plot_bi, owen_label ) ;

                  if ( bi_ht == 3 ) h_ratio_ht3 -> SetBinContent( ht3_plot_bi, Rqcd_val ) ;
                  if ( bi_ht == 3 ) h_ratio_ht3 -> SetBinError( ht3_plot_bi, Rqcd_val ) ;
                  if ( bi_ht == 3 ) h_ratio_ht3 -> GetXaxis() -> SetBinLabel( ht3_plot_bi, owen_label ) ;


                  fprintf( ofp_combine, "\n" ) ;
               }

            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp_combine ) ;
      printf("\n\n Wrote combine input file : %s\n\n", combine_output_file ) ;

      fclose( ofp_qcdratio ) ;
      printf("\n\n Wrote qcd ratio file : %s\n\n", qcdratio_output_file ) ;

      h_ratio_all -> SetMarkerStyle(20) ;
      h_ratio_nj1 -> SetMarkerStyle(20) ;
      h_ratio_nj2 -> SetMarkerStyle(20) ;
      h_ratio_nj3 -> SetMarkerStyle(20) ;
      h_ratio_nj4 -> SetMarkerStyle(20) ;
      h_ratio_ht1 -> SetMarkerStyle(20) ;
      h_ratio_ht2 -> SetMarkerStyle(20) ;
      h_ratio_ht3 -> SetMarkerStyle(20) ;

      h_ratio_all -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj1 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj2 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj3 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj4 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_ht1 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_ht2 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_ht3 -> GetXaxis() -> LabelsOption( "v" ) ;

      saveHist( hist_output_file,"h*") ;



    //----

      gStyle -> SetPadBottomMargin(0.35) ;

      TCanvas* can_all = new TCanvas( "can_all", "All 160 bins", 1900, 500 ) ;

      h_ldp_all_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_all_hadtau -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_all_znunu -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_all_data -> GetXaxis() -> LabelsOption( "v" ) ;

      h_ldp_all_lostlep -> SetFillColor( kBlue-10 ) ;
      h_ldp_all_hadtau  -> SetFillColor( kCyan-10 ) ;
      h_ldp_all_znunu   -> SetFillColor( kGreen-7 ) ;
      h_ldp_all_data -> SetLineWidth(3) ;

      TLegend* legend = new TLegend( 0.80, 0.80, 0.95, 0.95 ) ;
      legend -> AddEntry( h_ldp_all_znunu, "znunu" ) ;
      legend -> AddEntry( h_ldp_all_hadtau, "hadtau" ) ;
      legend -> AddEntry( h_ldp_all_lostlep, "lostlep" ) ;
      legend -> AddEntry( h_ldp_all_data, "data" ) ;

      THStack* hstack_ldp_all = new THStack( "hstack_ldp_all", "hstack_ldp_all" ) ;
      hstack_ldp_all -> Add( h_ldp_all_lostlep ) ;
      hstack_ldp_all -> Add( h_ldp_all_hadtau ) ;
      hstack_ldp_all -> Add( h_ldp_all_znunu ) ;

      h_ldp_all_data -> Draw() ;
      hstack_ldp_all -> Draw( "hist same" ) ;
      hstack_ldp_all -> Draw( "same" ) ;
      h_ldp_all_data -> Draw( "axis same" ) ;
      h_ldp_all_data -> Draw( "axig same" ) ;
      h_ldp_all_data -> Draw( "same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;


    //----

      TCanvas* can_ht3 = new TCanvas( "can_ht3", "HT3 bins", 1200, 500 ) ;

      h_ldp_ht3_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht3_hadtau -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht3_znunu -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht3_data -> GetXaxis() -> LabelsOption( "v" ) ;

      h_ldp_ht3_lostlep -> SetFillColor( kBlue-10 ) ;
      h_ldp_ht3_hadtau  -> SetFillColor( kCyan-10 ) ;
      h_ldp_ht3_znunu   -> SetFillColor( kGreen-7 ) ;
      h_ldp_ht3_data -> SetLineWidth(3) ;

      THStack* hstack_ldp_ht3 = new THStack( "hstack_ldp_ht3", "hstack_ldp_ht3" ) ;
      hstack_ldp_ht3 -> Add( h_ldp_ht3_lostlep ) ;
      hstack_ldp_ht3 -> Add( h_ldp_ht3_hadtau ) ;
      hstack_ldp_ht3 -> Add( h_ldp_ht3_znunu ) ;

      h_ldp_ht3_data -> Draw() ;
      hstack_ldp_ht3 -> Draw( "hist same" ) ;
      hstack_ldp_ht3 -> Draw( "same" ) ;
      h_ldp_ht3_data -> Draw( "axis same" ) ;
      h_ldp_ht3_data -> Draw( "axig same" ) ;
      h_ldp_ht3_data -> Draw( "same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;


    //----

      TCanvas* can_ht2 = new TCanvas( "can_ht2", "ht2 bins", 1200, 500 ) ;

      h_ldp_ht2_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht2_hadtau -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht2_znunu -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht2_data -> GetXaxis() -> LabelsOption( "v" ) ;

      h_ldp_ht2_lostlep -> SetFillColor( kBlue-10 ) ;
      h_ldp_ht2_hadtau  -> SetFillColor( kCyan-10 ) ;
      h_ldp_ht2_znunu   -> SetFillColor( kGreen-7 ) ;
      h_ldp_ht2_data -> SetLineWidth(3) ;

      THStack* hstack_ldp_ht2 = new THStack( "hstack_ldp_ht2", "hstack_ldp_ht2" ) ;
      hstack_ldp_ht2 -> Add( h_ldp_ht2_lostlep ) ;
      hstack_ldp_ht2 -> Add( h_ldp_ht2_hadtau ) ;
      hstack_ldp_ht2 -> Add( h_ldp_ht2_znunu ) ;

      h_ldp_ht2_data -> Draw() ;
      hstack_ldp_ht2 -> Draw( "hist same" ) ;
      hstack_ldp_ht2 -> Draw( "same" ) ;
      h_ldp_ht2_data -> Draw( "axis same" ) ;
      h_ldp_ht2_data -> Draw( "axig same" ) ;
      h_ldp_ht2_data -> Draw( "same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;



    //----

      TCanvas* can_ht1 = new TCanvas( "can_ht1", "ht1 bins", 1200, 500 ) ;

      h_ldp_ht1_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht1_hadtau -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht1_znunu -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_ht1_data -> GetXaxis() -> LabelsOption( "v" ) ;

      h_ldp_ht1_lostlep -> SetFillColor( kBlue-10 ) ;
      h_ldp_ht1_hadtau  -> SetFillColor( kCyan-10 ) ;
      h_ldp_ht1_znunu   -> SetFillColor( kGreen-7 ) ;
      h_ldp_ht1_data -> SetLineWidth(3) ;

      THStack* hstack_ldp_ht1 = new THStack( "hstack_ldp_ht1", "hstack_ldp_ht1" ) ;
      hstack_ldp_ht1 -> Add( h_ldp_ht1_lostlep ) ;
      hstack_ldp_ht1 -> Add( h_ldp_ht1_hadtau ) ;
      hstack_ldp_ht1 -> Add( h_ldp_ht1_znunu ) ;

      h_ldp_ht1_data -> Draw() ;
      hstack_ldp_ht1 -> Draw( "hist same" ) ;
      hstack_ldp_ht1 -> Draw( "same" ) ;
      h_ldp_ht1_data -> Draw( "axis same" ) ;
      h_ldp_ht1_data -> Draw( "axig same" ) ;
      h_ldp_ht1_data -> Draw( "same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;




   } // gen_combine_input2


  //=======================================================================================

      void read_pars( const char* model_pars_file ) {

         ifstream ifs_model_pars ;
         ifs_model_pars.open( model_pars_file ) ;
         if ( !ifs_model_pars.good() ) { printf("\n\n *** Problem opening %s\n\n", model_pars_file ) ; gSystem->Exit(-1) ; }

         float val, err1, err2, rel_err ;
         char pname[100] ;

       //---
         sprintf( pname, "Kqcd_HT1" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_ht[1] = val ;
         par_err_ht_fit[1] = err1 ;
         par_err_ht_syst[1] = val*err2 ;
         par_rel_err_ht[1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;


         sprintf( pname, "Kqcd_HT2" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_ht[2] = val ;
         par_err_ht_fit[2] = err1 ;
         par_err_ht_syst[2] = val*err2 ;
         par_rel_err_ht[2] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Kqcd_HT3" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_ht[3] = val ;
         par_err_ht_fit[3] = err1 ;
         par_err_ht_syst[3] = val*err2 ;
         par_rel_err_ht[3] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

       //---
         sprintf( pname, "Sqcd_njet1" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_njet[1] = val ;
         par_err_njet_fit[1] = err1 ;
         par_err_njet_syst[1] = val*err2 ;
         par_rel_err_njet[1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_njet2" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_njet[2] = val ;
         par_err_njet_fit[2] = err1 ;
         par_err_njet_syst[2] = val*err2 ;
         par_rel_err_njet[2] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_njet3" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_njet[3] = val ;
         par_err_njet_fit[3] = err1 ;
         par_err_njet_syst[3] = val*err2 ;
         par_rel_err_njet[3] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_njet4" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_njet[4] = val ;
         par_err_njet_fit[4] = err1 ;
         par_err_njet_syst[4] = val*err2 ;
         par_rel_err_njet[4] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;



       //---
         sprintf( pname, "Sqcd_mhtc_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_hth[1] = val ;
         par_err_mht_hth[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_hth[1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht1_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_hth[2] = val ;
         par_err_mht_hth[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_hth[2] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht2_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_hth[3] = val ;
         par_err_mht_hth[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_hth[3] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht3_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_hth[4] = val ;
         par_err_mht_hth[4] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_hth[4] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht4_hth" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_hth[5] = val ;
         par_err_mht_hth[5] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_hth[5] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

       //---
         sprintf( pname, "Sqcd_mhtc_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htm[1] = val ;
         par_err_mht_htm[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htm[1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht1_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htm[2] = val ;
         par_err_mht_htm[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htm[2] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht2_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htm[3] = val ;
         par_err_mht_htm[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htm[3] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht3_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htm[4] = val ;
         par_err_mht_htm[4] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htm[4] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht4_htm" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htm[5] = val ;
         par_err_mht_htm[5] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htm[5] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

       //---
         sprintf( pname, "Sqcd_mhtc_htl" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htl[1] = val ;
         par_err_mht_htl[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htl[1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht1_htl" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htl[2] = val ;
         par_err_mht_htl[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htl[2] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_mht2_htl" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_mht_htl[3] = val ;
         par_err_mht_htl[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_mht_htl[3] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;





       //---
         sprintf( pname, "Sqcd_nb0" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_nb[1] = val ;
         par_err_nb[1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_nb[1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_nb1" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_nb[2] = val ;
         par_err_nb[2] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_nb[2] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_nb2" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_nb[3] = val ;
         par_err_nb[3] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_nb[3] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

         sprintf( pname, "Sqcd_nb3" ) ;
         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_nb[4] = val ;
         par_err_nb[4] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_nb[4] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;

      } // read_pars

  //=======================================================================================

   void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2, float& rel_err ) {

      val = 0. ;
      err1 = 0. ;
      err2 = 0. ;
      rel_err = 0. ;

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
            if ( val > 0 ) {
               float err_total = sqrt( err1*err1 + val*err2 * val*err2 ) ;
               rel_err = err_total / val ;
            }
            return ;
         }
      }

      printf("\n\n *** get_par : Failed to find parameter %s\n\n", pname ) ;
      gSystem -> Exit(-1) ;

   } // get_par

  //=======================================================================================




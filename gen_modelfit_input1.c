#ifndef gen_modelfit_input1_c
#define gen_modelfit_input1_c

#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPad.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <fstream>

#include "binning.h"

   void gen_modelfit_input(
         const char* data_file    = "outputfiles/nbsum-input-data.txt",
         const char* lostlep_file = "outputfiles/nbsum-input-lostlep.txt",
         const char* hadtau_file  = "outputfiles/nbsum-input-hadtau.txt",
         const char* znunu_file   = "outputfiles/nbsum-input-znunu.txt",
         const char* output_hist_file = "outputfiles/modelfit-input-data.root"
                         ) {


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

      int    nobs_ldp[10][10] ;
      double nonqcd_val_ldp[10][10] ;
      double nonqcd_err_ldp[10][10] ;

      int    nobs_hdp[10][10] ;
      double nonqcd_val_hdp[10][10] ;
      double nonqcd_err_hdp[10][10] ;

      TH1F* h_ratio = new TH1F( "h_ratio", "H/L ratio", nBinsHT * nb_nj, 0.5, nBinsHT * nb_n + 0.5 ) ;

      TH1F* h_ldp_lostlep = new TH1F( "h_ldp_lostlep", "ldp, lostlep", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;
      TH1F* h_ldp_hadtau = new TH1F( "h_ldp_hadtau", "ldp, hadtau", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;
      TH1F* h_ldp_znunu = new TH1F( "h_ldp_znunu", "ldp, znunu", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;
      TH1F* h_ldp_data = new TH1F( "h_ldp_data", "ldp, data", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;


      TH1F* h_hdp_lostlep = new TH1F( "h_hdp_lostlep", "hdp, lostlep", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;
      TH1F* h_hdp_hadtau = new TH1F( "h_hdp_hadtau", "hdp, hadtau", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;
      TH1F* h_hdp_znunu = new TH1F( "h_hdp_znunu", "hdp, znunu", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;
      TH1F* h_hdp_data = new TH1F( "h_hdp_data", "hdp, data", nBinsHT * nb_n, 0.5, nBinsHT * nb_n + 0.5 ) ;


      int bi_hist(0) ;

      for ( int bi_ht=1; bi_ht<=nBinsHT; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

               bi_hist ++ ;

               TString line ;
               int global_bi, region_bi ;
               char region_tag[5] ;
               char label[100] ;
               int r_nobs_ldp, r_nobs_hdp ;
               float nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst ;
               float nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ;


              //-------
               line.ReadLine( ifs_data ) ;
               sscanf( line.Data(), "%s  %d  %d", label, &r_nobs_ldp, &r_nobs_hdp ) ;
               printf( "  Data    :   %s  Nldp = %5d ,                        Nhdp = %5d\n",
                    label, r_nobs_ldp, r_nobs_hdp ) ;

               nobs_ldp[bi_ht][bi_nj] = r_nobs_ldp ;
               nobs_hdp[bi_ht][bi_nj] = r_nobs_hdp ;

               h_ldp_data -> SetBinContent( bi_hist, r_nobs_ldp ) ;
               h_ldp_data -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

               h_hdp_data -> SetBinContent( bi_hist, r_nobs_hdp ) ;
               h_hdp_data -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;


               double total_bg_ldp_val = 0. ;
               double total_bg_ldp_err2 = 0. ;

               double total_bg_hdp_val = 0. ;
               double total_bg_hdp_err2 = 0. ;

              //-------
               line.ReadLine( ifs_lostlep ) ;
               sscanf( line.Data(), "%s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                    label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;
               printf( "  Lostlep :   %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                    label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

               total_bg_ldp_val += nbg_ldp_val ;
               total_bg_hdp_val += nbg_hdp_val ;
               total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
               total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

               h_ldp_lostlep -> SetBinContent( bi_hist, nbg_ldp_val ) ;
               h_ldp_lostlep -> SetBinError( bi_hist, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
               h_ldp_lostlep -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

               h_hdp_lostlep -> SetBinContent( bi_hist, nbg_hdp_val ) ;
               h_hdp_lostlep -> SetBinError( bi_hist, sqrt( pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ) ) ;
               h_hdp_lostlep -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

              //-------
               line.ReadLine( ifs_hadtau ) ;
               sscanf( line.Data(), "%s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                    label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;
               printf( "  Hadtau  :   %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                    label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

               total_bg_ldp_val += nbg_ldp_val ;
               total_bg_hdp_val += nbg_hdp_val ;
               total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
               total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

               h_ldp_hadtau -> SetBinContent( bi_hist, nbg_ldp_val ) ;
               h_ldp_hadtau -> SetBinError( bi_hist, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
               h_ldp_hadtau -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

               h_hdp_hadtau -> SetBinContent( bi_hist, nbg_hdp_val ) ;
               h_hdp_hadtau -> SetBinError( bi_hist, sqrt( pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ) ) ;
               h_hdp_hadtau -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

              //-------
               line.ReadLine( ifs_znunu ) ;
               sscanf( line.Data(), "%s  %f +/- %f +/- %f    %f +/- %f +/- %f",
                    label,  &nbg_ldp_val, &nbg_ldp_stat, &nbg_ldp_syst,  &nbg_hdp_val, &nbg_hdp_stat, &nbg_hdp_syst ) ;
               printf( "  Znunu   :   %s  Nldp = %7.1f +/- %5.1f +/- %5.1f ,  Nhdp = %7.1f +/- %5.1f +/- %5.1f\n",
                    label, nbg_ldp_val, nbg_ldp_stat, nbg_ldp_syst,  nbg_hdp_val, nbg_hdp_stat, nbg_hdp_syst ) ;

               total_bg_ldp_val += nbg_ldp_val ;
               total_bg_hdp_val += nbg_hdp_val ;
               total_bg_ldp_err2 += pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ;
               total_bg_hdp_err2 += pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ;

               h_ldp_znunu -> SetBinContent( bi_hist, nbg_ldp_val ) ;
               h_ldp_znunu -> SetBinError( bi_hist, sqrt( pow( nbg_ldp_stat, 2. ) + pow( nbg_ldp_syst, 2. ) ) ) ;
               h_ldp_znunu -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

               h_hdp_znunu -> SetBinContent( bi_hist, nbg_hdp_val ) ;
               h_hdp_znunu -> SetBinError( bi_hist, sqrt( pow( nbg_hdp_stat, 2. ) + pow( nbg_hdp_syst, 2. ) ) ) ;
               h_hdp_znunu -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;


               double total_bg_ldp_err = sqrt( total_bg_ldp_err2 ) ;
               double total_bg_hdp_err = sqrt( total_bg_hdp_err2 ) ;

               printf( " total BG :   %s  Nldp = %7.1f +/- %5.1f           ,  Nhdp = %7.1f +/- %5.1f          \n",
                   label, total_bg_ldp_val, total_bg_ldp_err,   total_bg_hdp_val, total_bg_hdp_err ) ;

               double qcd_ldp_val = r_nobs_ldp - total_bg_ldp_val ;
               double qcd_ldp_err = sqrt( r_nobs_ldp + total_bg_ldp_err*total_bg_ldp_err ) ;

               double qcd_hdp_val = r_nobs_hdp - total_bg_hdp_val ;
               double qcd_hdp_err = sqrt( r_nobs_hdp + total_bg_hdp_err*total_bg_hdp_err ) ;
               double qcd_hdp_rel_err(0.) ;
               if ( qcd_hdp_val > 0 ) qcd_hdp_rel_err = qcd_hdp_err / qcd_hdp_val ;

               double ratio_val(0.) ;
               double ratio_err(0.) ;
               double ratio_rel_err(0.) ;
               if ( qcd_ldp_val != 0 ) {
                  ratio_val = qcd_hdp_val / qcd_ldp_val ;
                  if ( qcd_hdp_val != 0 ) {
                     ratio_err = fabs(ratio_val) * sqrt( pow( qcd_ldp_err/qcd_ldp_val, 2.) + pow( qcd_hdp_err/qcd_hdp_val, 2. ) ) ;
                     ratio_rel_err = ratio_err / ratio_val ;
                  }
               }

               printf( " QCD      :   %s  Nldp = %7.1f +/- %5.1f           ,  Nhdp = %7.1f +/- %5.1f  (%6.3f)      R(H/L) = %6.3f +/- %6.3f  (%6.3f)\n",
                   label, qcd_ldp_val, qcd_ldp_err,   qcd_hdp_val, qcd_hdp_err, qcd_hdp_rel_err,    ratio_val, ratio_err, ratio_rel_err ) ;

               printf("\n") ;

               if ( !(bi_ht==1 && bi_nj>nb_nj-2) ) {  // skip top two njets bins for lowest HT
                  h_ratio -> SetBinContent( bi_hist, ratio_val ) ;
                  h_ratio -> SetBinError( bi_hist, ratio_err ) ;
               } else {
                  h_ratio -> SetBinContent( bi_hist, -9 ) ;
                  h_ratio -> SetBinError( bi_hist, 0. ) ;
               }
               h_ratio -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;


         } // bi_nj
      } // bi_ht




      char fname[1000] ;

      h_ldp_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_hadtau -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_znunu -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_data -> GetXaxis() -> LabelsOption( "v" ) ;

      h_hdp_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
      h_hdp_hadtau -> GetXaxis() -> LabelsOption( "v" ) ;
      h_hdp_znunu -> GetXaxis() -> LabelsOption( "v" ) ;
      h_hdp_data -> GetXaxis() -> LabelsOption( "v" ) ;

      h_ldp_lostlep -> SetFillColor( kBlue-10 ) ;
      h_ldp_hadtau  -> SetFillColor( kCyan-10 ) ;
      h_ldp_znunu   -> SetFillColor( kGreen-7 ) ;
      h_ldp_data -> SetLineWidth(3) ;

      h_hdp_lostlep -> SetFillColor( kBlue-10 ) ;
      h_hdp_hadtau  -> SetFillColor( kCyan-10 ) ;
      h_hdp_znunu   -> SetFillColor( kGreen-7 ) ;
      h_hdp_data -> SetLineWidth(3) ;

      TLegend* legend = new TLegend( 0.80, 0.80, 0.95, 0.95 ) ;
      legend -> AddEntry( h_ldp_znunu, "znunu" ) ;
      legend -> AddEntry( h_ldp_hadtau, "hadtau" ) ;
      legend -> AddEntry( h_ldp_lostlep, "lostlep" ) ;
      legend -> AddEntry( h_ldp_data, "data" ) ;


      h_ratio -> GetXaxis() -> LabelsOption( "v" ) ;

      gStyle -> SetPadBottomMargin(0.20) ;
      gStyle -> SetOptStat(0) ;

     //-------

      TCanvas* can_ratio = new TCanvas( "can_ratio", "Ratio", 900, 800 ) ;

      h_ratio -> SetMarkerStyle(20) ;
      h_ratio -> SetMinimum(-0.1) ;
      h_ratio -> Draw() ;
      gPad -> SetGridy(1) ;


     //-------

      TCanvas* can_ldp = new TCanvas( "can_ldp", "LDP", 900, 800 ) ;

      THStack* hstack_ldp = new THStack( "hstack_ldp", "hstack_ldp" ) ;
      hstack_ldp -> Add( h_ldp_lostlep ) ;
      hstack_ldp -> Add( h_ldp_hadtau ) ;
      hstack_ldp -> Add( h_ldp_znunu ) ;

      h_ldp_data -> Draw() ;
      hstack_ldp -> Draw( "hist same" ) ;
      hstack_ldp -> Draw( "same" ) ;
      h_ldp_data -> Draw( "axis same" ) ;
      h_ldp_data -> Draw( "axig same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;

      can_ldp -> Update() ; can_ldp -> Draw() ;
      can_ldp -> SaveAs( "outputfiles/modelfit-input-ldp-liny-full.pdf" ) ;

      h_ldp_data -> SetMaximum( 1.1 * (h_ldp_data -> GetBinContent(1)) ) ;
      can_ldp -> Update() ; can_ldp -> Draw() ;
      can_ldp -> SaveAs( "outputfiles/modelfit-input-ldp-liny-zoom1.pdf" ) ;

      h_ldp_data -> SetMaximum( 1.1 * (h_ldp_data -> GetBinContent(11)) ) ;
      can_ldp -> Update() ; can_ldp -> Draw() ;
      can_ldp -> SaveAs( "outputfiles/modelfit-input-ldp-liny-zoom2.pdf" ) ;

      h_ldp_data -> SetMaximum( 1.1 * (h_ldp_data -> GetBinContent(2)) ) ;
      can_ldp -> Update() ; can_ldp -> Draw() ;
      can_ldp -> SaveAs( "outputfiles/modelfit-input-ldp-liny-zoom3.pdf" ) ;

      h_ldp_data -> SetMaximum( 3 * (h_ldp_data -> GetBinContent(5)) ) ;
      h_ldp_data -> SetMinimum(0.5) ;
      gPad -> SetLogy(1) ;
      can_ldp -> Update() ; can_ldp -> Draw() ;
      can_ldp -> SaveAs( "outputfiles/modelfit-input-ldp-logy-full.pdf" ) ;



     //-------

      TCanvas* can_hdp = new TCanvas( "can_hdp", "hdp", 900, 800 ) ;

      THStack* hstack_hdp = new THStack( "hstack_hdp", "hstack_hdp" ) ;
      hstack_hdp -> Add( h_hdp_lostlep ) ;
      hstack_hdp -> Add( h_hdp_hadtau ) ;
      hstack_hdp -> Add( h_hdp_znunu ) ;

      h_hdp_data -> Draw() ;
      hstack_hdp -> Draw( "hist same" ) ;
      hstack_hdp -> Draw( "same" ) ;
      h_hdp_data -> Draw( "axis same" ) ;
      h_hdp_data -> Draw( "axig same" ) ;

      legend -> Draw() ;

      gPad -> SetGridy(1) ;

      can_hdp -> Update() ; can_hdp -> Draw() ;
      can_hdp -> SaveAs( "outputfiles/modelfit-input-hdp-liny-full.pdf" ) ;

      h_hdp_data -> SetMaximum( 1.1 * (h_hdp_data -> GetBinContent(5)) ) ;
      can_hdp -> Update() ; can_hdp -> Draw() ;
      can_hdp -> SaveAs( "outputfiles/modelfit-input-hdp-liny-zoom1.pdf" ) ;

      h_hdp_data -> SetMaximum( 1.1 * (h_hdp_data -> GetBinContent(9)) ) ;
      can_hdp -> Update() ; can_hdp -> Draw() ;
      can_hdp -> SaveAs( "outputfiles/modelfit-input-hdp-liny-zoom2.pdf" ) ;

      h_hdp_data -> SetMaximum( 3 * (h_hdp_data -> GetBinContent(1)) ) ;
      h_hdp_data -> SetMinimum(0.5) ;
      gPad -> SetLogy(1) ;
      can_hdp -> Update() ; can_hdp -> Draw() ;
      can_hdp -> SaveAs( "outputfiles/modelfit-input-hdp-logy-full.pdf" ) ;





      TFile rf( output_hist_file, "recreate" ) ;
      h_ratio -> Write() ;
      rf.Close() ;
      printf("\n\n Saved ratio histogram in %s\n\n", output_hist_file ) ;



   } // gen_combine_input




#endif

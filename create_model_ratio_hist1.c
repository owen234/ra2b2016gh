#ifndef create_model_ratio_hist1_c
#define create_model_ratio_hist1_c

#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>
#include "binning.h"
#include "histio.c"
#include "read_pars.h"

   TH1F* get_hist( const char* hname );

   void create_model_ratio_hist1( const char* model_pars_file = "outputfiles/model-pars-qcdmc3.txt",
                                  const char* qcd_ratio_file = "outputfiles/qcdmc-ratio-v3.root" ) {
      setup_bins(); 
      gDirectory -> Delete( "h*" ) ;

      loadHist( qcd_ratio_file, "qcdmc" ) ;

      read_pars( model_pars_file ) ;

      TH1F* h_ratio_all = new TH1F( "h_ratio_all", "QCD model H/L ratio", (nb_htmht-3) * nb_nb * nb_nj, 0.5, (nb_htmht-3) * nb_nb * nb_nj + 0.5 ) ;

      TH1F* h_max_ldp_weight_search_bins = get_hist( "h_max_ldp_weight_search_bins_qcdmc" ) ;
      TH1F* h_ldp_search_bins = get_hist( "h_ldp_search_bins_qcdmc" ) ;
      TH1F* h_hdp_search_bins = get_hist( "h_hdp_search_bins_qcdmc" ) ;
      TH1F* h_ratio_qcdmc = get_hist( "h_ratio_qcdmc" ) ;

      int bi_hist(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=4; bi_htmht<=nb_htmht; bi_htmht++ ) {

               bi_hist++ ;

               int bi_ht, bi_mht ;
               htmht_bin_to_ht_and_mht_bins( bi_htmht, bi_ht, bi_mht ) ;

               char label[100] ;
               sprintf( label, " %3d Nj%d-Nb%d-MHT%d-HT%d (%d)", bi_hist, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht-3 ) ;

               double model_ratio_val = 0;
               double model_ratio_err = 0;

                  model_ratio_val = par_val_ht[bi_ht] * par_val_njet[bi_nj] * par_val_ht_mht[bi_ht][bi_mht] * par_val_nb[bi_nb] ;
                  model_ratio_err = model_ratio_val * sqrt(
                         pow( par_err_ht_fit[bi_ht]/par_val_ht[bi_ht], 2. )
                      +  pow( par_err_ht_syst[bi_ht]/par_val_ht[bi_ht], 2. )
                      +  pow( par_err_njet_fit[bi_nj]/par_val_njet[bi_nj], 2. )
                      +  pow( par_err_njet_syst[bi_nj]/par_val_njet[bi_nj], 2. )
                      +  pow( par_err_ht_mht[bi_ht][bi_mht]/par_val_ht_mht[bi_ht][bi_mht], 2. )
                      +  pow( par_err_nb[bi_nb]/par_val_nb[bi_nb], 2. )
                    ) ;
                  printf("  %s : Nj %6.4f Nb %6.4f MHT %6.4f HT %6.4f  model ratio = %6.4f +/- %6.4f\n", label,
                    par_val_njet[bi_nj], par_val_nb[bi_nb], par_val_ht_mht[bi_ht][bi_mht], par_val_ht[bi_ht], model_ratio_val, model_ratio_err  ) ;

               h_ratio_all -> GetXaxis() -> SetBinLabel( bi_hist, label ) ;

               if ( bi_nj>(nb_nj-2) && bi_ht==1 ) continue ; // skip top two njets bins for lowest HT.

               h_ratio_all -> SetBinContent( bi_hist, model_ratio_val ) ;
               h_ratio_all -> SetBinError( bi_hist, model_ratio_err ) ;

            } // bi_htmht
         } // bi_nb
      } // bi_nj

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin(0.30) ;

      h_ratio_all -> SetMarkerStyle( 22 ) ;
      h_ratio_all -> SetMarkerColor( 2 ) ;

      h_ratio_all -> GetXaxis() -> LabelsOption("v") ;
      h_ratio_all -> Draw() ;
      gPad -> SetGridy(1) ;


     //---------------

      TH1F* h_ratio_qcdmc_minus_model = new TH1F( "h_ratio_qcdmc_minus_model", "QCD H/L ratio difference (QCD MC - model)", (nb_htmht-3) * nb_nb * nb_nj, 0.5, (nb_htmht-3) * nb_nb * nb_nj + 0.5 ) ;

      printf("\n\n") ;
      for ( int bi=1; bi<=(nb_htmht-3) * nb_nb * nb_nj; bi++ ) { //loop over search bins
         float model_val = h_ratio_all -> GetBinContent( bi ) ;
         float qcdmc_val = h_ratio_qcdmc -> GetBinContent( bi ) ;
         float ldp_val = h_ldp_search_bins -> GetBinContent( bi ) ;
         float hdp_val = h_hdp_search_bins -> GetBinContent( bi ) ;
         float max_ldp_weight = h_max_ldp_weight_search_bins -> GetBinContent( bi ) ;
         char label[100] ;
         sprintf( label, "%s", h_ratio_all -> GetXaxis() -> GetBinLabel( bi ) ) ;
         float diff_val(0.) ;
         float diff_err(0.) ;
         printf(" debug1 : model bin label = %s , qcdmc bin label = %s\n", h_ratio_all -> GetXaxis() -> GetBinLabel( bi ), h_ratio_qcdmc -> GetXaxis() -> GetBinLabel( bi ) ) ;
         if ( hdp_val > 0 ) {
            diff_val = qcdmc_val - model_val ;
            diff_err = diff_val ;
            printf("  %40s : LDP %7.1f  HDP %7.1f   max LDP weight %5.3f, diff err = %5.3f\n", label, ldp_val, hdp_val, max_ldp_weight, diff_err ) ;
         } else {
            diff_val = 0. ;
            if ( ldp_val > 0 ) {
               diff_err = max_ldp_weight / ldp_val ;
               printf("  %40s : LDP %7.1f  HDP %7.1f   max LDP weight %5.3f,  zero HDP H/L err = %5.3f\n", label, ldp_val, hdp_val, max_ldp_weight, diff_err ) ;
            } else {
               //diff_err = 0.5 ;
               //diff_err = 0.2;
               diff_err = 0.0;
               printf("  %40s : LDP %7.1f  HDP %7.1f   max LDP weight %5.3f,  *** both zero\n", label, ldp_val, hdp_val, max_ldp_weight ) ;
            }
         }
         h_ratio_qcdmc_minus_model -> SetBinContent( bi, diff_val ) ;
         h_ratio_qcdmc_minus_model -> SetBinError( bi, diff_err ) ;
         h_ratio_qcdmc_minus_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
      } // bi
      printf("\n\n") ;

      h_ratio_qcdmc_minus_model -> GetXaxis() -> LabelsOption( "v" ) ;


      saveHist("outputfiles/model-ratio-hist1.root", "h*" ) ;

   } // create_model_ratio_hist1


#ifndef get_hist_
#define get_hist_

   TH1F* get_hist( const char* hname ) {
      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** Missing histogram : %s\n\n", hname ) ;
         gDirectory -> ls() ;
         gSystem -> Exit( -1 ) ;
      }
      return hp ;
   } // get_hist

#endif
//===============================================================================
#endif

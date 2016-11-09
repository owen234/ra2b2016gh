#include "TSystem.h"
#include "histio.c"
#include <fstream>

#include "binning.h"
#include "read_pars.h"
#include "num_to_str.h"

  //--------

   void calc_model_ratios_v4(
              const char* model_pars_file = "outputfiles/model-pars-qcdmc3.txt"
           ) {
      setup_bins();
      read_pars( model_pars_file ) ;

      par_val_nb[5] = par_val_nb[1] ; // val for sum same as nb0
      par_err_nb[5] = par_err_nb[1] ; // val for sum same as nb0

      gDirectory -> Delete( "h*" ) ;


      TH1F* h_model_ratio[4][6] ;
      TH1F* h_model_dr   [4][6] ;

      for ( int bi_ht = 1; bi_ht <= nBinsHT; bi_ht++)
      {
         TString ht_str, ht_str2;
         if ( bi_ht == 1 ) { ht_str = "htl"; ht_str2 = ", HT low";    }
         if ( bi_ht == 2 ) { ht_str = "htm"; ht_str2 = ", HT medium"; }
         if ( bi_ht == 3 ) { ht_str = "hth"; ht_str2 = ", HT high";   }

         for ( int bi_nb = 1; bi_nb <= nb_nb; bi_nb++)
	 {

         h_model_ratio[bi_ht][bi_nb] = new TH1F( "h_model_ratio_"+ht_str+"_nb"+num_to_str(bi_nb-1), "Model, H/L ratio, Nb"       +num_to_str(bi_nb-1)+ht_str2, 
                                               nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

         h_model_dr   [bi_ht][bi_nb] = new TH1F( "h_model_dr_"   +ht_str+"_nb"+num_to_str(bi_nb-1), "Model, double H/L ratio, Nb"+num_to_str(bi_nb-1)+ht_str2,
                                               nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

         }//ni_nb
         h_model_ratio[bi_ht][nb_nb+1] = new TH1F( "h_model_ratio_"+ht_str+"_nbsum", "Model, H/L ratio, Nbsum"       +ht_str2, nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
         h_model_dr   [bi_ht][nb_nb+1] = new TH1F( "h_model_dr_"   +ht_str+"_nbsum", "Model, double H/L ratio, Nbsum"+ht_str2, nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

      }//bi_ht

      int bi_search(0) ;

      for ( int nji = 1; nji <= nb_nj; nji++ ) {
         for ( int nbi = 1; nbi <= (nb_nb+1) ; nbi++ ) {
            for ( int mbi = 1; mbi <= nb_mht; mbi++ ) {
               char label[100] ;
               if ( mbi==1 ) {
                  sprintf( label, "Nj%d-MHTC", nji ) ;
               } else {
                  sprintf( label, "Nj%d-MHT%d", nji, mbi-1 ) ;
               }


               int htb_low(0) ;
               if ( mbi<4 ) { htb_low = 1 ; } else { htb_low = 2 ; }
               for ( int hti = htb_low; hti<=nBinsHT; hti++ ) {

                  if ( is_this_bin_excluded(nji, nbi, hti, mbi) ) continue;
   
                  if ( mbi > 1 && nbi <= nb_nb ) bi_search ++ ;

                  //if ( hti == 1 && nji > 2 ) continue ;

                  float par_val_mht(0.) ;
                  float par_err_mht(0.) ;
                  TH1F* h_ratio(0x0) ;
                  TH1F* h_dr(0x0) ;
                  int hbi(0) ;

                  par_val_mht = par_val_ht_mht[hti][mbi] ;
                  par_err_mht = par_err_ht_mht[hti][mbi] ;
                  h_ratio     = h_model_ratio [hti][nbi] ;
                  h_dr        = h_model_dr    [hti][nbi] ;
                  hbi = (nji-1)*nb_mht + mbi ;

                  float model_ratio_val(0.) ;
                  float model_ratio_err2(0.) ;

                  model_ratio_val =  par_val_ht[hti] * par_val_njet[nji] * par_val_mht * par_val_nb[nbi] ;


                  if ( par_val_ht[hti]   > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_ht_fit[hti]    / par_val_ht[hti]  , 2. ) ;
                  if ( par_val_ht[hti]   > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_ht_syst[hti]   / par_val_ht[hti]  , 2. ) ;

                  if ( par_val_njet[nji] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_njet_fit[nji]  / par_val_njet[nji], 2. ) ;
                  if ( par_val_njet[nji] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_njet_syst[nji] / par_val_njet[nji], 2. ) ;

                  if ( par_val_nb[nbi]   > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_nb[nbi]        / par_val_nb[nbi]  , 2. ) ;

                  if ( par_val_mht       > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_mht            / par_val_mht      , 2. ) ;


                  float model_ratio_err = sqrt( model_ratio_err2 ) ;


                  float model_dr_val(0.) ;
                  float model_dr_err(0.) ;

                  model_dr_val = par_val_mht ;
                  model_dr_err = par_err_mht ;



                  h_ratio -> SetBinContent( hbi, model_ratio_val ) ;
                  h_ratio -> SetBinError  ( hbi, model_ratio_err ) ;
                  h_ratio -> GetXaxis() -> SetBinLabel( hbi, label ) ;

                  if ( mbi > 1 ) {
                     h_dr -> SetBinContent( hbi, model_dr_val ) ;
                     h_dr -> SetBinError( hbi, model_dr_err ) ;
                     h_dr -> GetXaxis() -> SetBinLabel( hbi, label ) ;
                  }

                  int print_bi(-1) ;
                  if ( mbi > 1 && nbi <= nb_nb ) {
                     print_bi = bi_search ;
                  } else {
                     print_bi = -1 ;
                  }

                   printf("  %3d :  Nj%d Nb%d MHT%d HT%d :  HT %7.4f * Nj %7.4f * MHT %7.4f * Nb %7.4f  = %7.4f +/- %7.4f\n",
                          print_bi, nji, nbi-1, mbi, hti,
                          par_val_ht[hti], par_val_njet[nji], par_val_mht, par_val_nb[nbi],
                          model_ratio_val, model_ratio_err ) ;


               } // hti

            } // mbi
         } // nbi
      } // nji


      printf("\n\n Saving histograms in : outputfiles/calc-model-ratios-v4.root\n") ;
      saveHist( "outputfiles/calc-model-ratios-v4.root", "h*" ) ;


   } // calc_model_ratios_v4

  //=======================================================================================

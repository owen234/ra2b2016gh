

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

      float par_val_nb[6] ;
      float par_err_nb[6] ;

      int nb_ht(3) ;
      int nb_nj(4) ;
      int nb_nb(4) ;
      int nb_htmht(13) ;
      int nb_mht(5) ;


  //--------

   void read_pars( const char* model_pars_file ) ;
   void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2 ) ;

  //--------

   void calc_model_ratios_v4(
              const char* model_pars_file = "model-pars-qcdmc3.txt"
           ) {

      read_pars( model_pars_file ) ;

      par_val_nb[5] = par_val_nb[1] ; // val for sum same as nb0
      par_err_nb[5] = par_err_nb[1] ; // val for sum same as nb0

      gDirectory -> Delete( "h*" ) ;


      TH1F* h_model_ratio_hth[6] ;
      h_model_ratio_hth[1] = new TH1F( "h_model_ratio_hth_nb0", "Model, H/L ratio, Nb0, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_hth[2] = new TH1F( "h_model_ratio_hth_nb1", "Model, H/L ratio, Nb1, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_hth[3] = new TH1F( "h_model_ratio_hth_nb2", "Model, H/L ratio, Nb2, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_hth[4] = new TH1F( "h_model_ratio_hth_nb3", "Model, H/L ratio, Nb3, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_hth[5] = new TH1F( "h_model_ratio_hth_nbsum", "Model, H/L ratio, Nbsum, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

      TH1F* h_model_ratio_htm[6] ;
      h_model_ratio_htm[1] = new TH1F( "h_model_ratio_htm_nb0", "Model, H/L ratio, Nb0, HT medium", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_htm[2] = new TH1F( "h_model_ratio_htm_nb1", "Model, H/L ratio, Nb1, HT medium", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_htm[3] = new TH1F( "h_model_ratio_htm_nb2", "Model, H/L ratio, Nb2, HT medium", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_htm[4] = new TH1F( "h_model_ratio_htm_nb3", "Model, H/L ratio, Nb3, HT medium", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_ratio_htm[5] = new TH1F( "h_model_ratio_htm_nbsum", "Model, H/L ratio, Nbsum, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

      TH1F* h_model_ratio_htl[6] ;
      h_model_ratio_htl[1] = new TH1F( "h_model_ratio_htl_nb0", "Model, H/L ratio, Nb0, HT low", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_ratio_htl[2] = new TH1F( "h_model_ratio_htl_nb1", "Model, H/L ratio, Nb1, HT low", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_ratio_htl[3] = new TH1F( "h_model_ratio_htl_nb2", "Model, H/L ratio, Nb2, HT low", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_ratio_htl[4] = new TH1F( "h_model_ratio_htl_nb3", "Model, H/L ratio, Nb3, HT low", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_ratio_htl[5] = new TH1F( "h_model_ratio_htl_nbsum", "Model, H/L ratio, Nbsum, HT high", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;




      TH1F* h_model_dr_hth[6] ;
      h_model_dr_hth[1] = new TH1F( "h_model_dr_hth_nb0", "Model, double H/L ratio, Nb0, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_hth[2] = new TH1F( "h_model_dr_hth_nb1", "Model, double H/L ratio, Nb1, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_hth[3] = new TH1F( "h_model_dr_hth_nb2", "Model, double H/L ratio, Nb2, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_hth[4] = new TH1F( "h_model_dr_hth_nb3", "Model, double H/L ratio, Nb3, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_hth[5] = new TH1F( "h_model_dr_hth_nbsum", "Model, double H/L ratio, Nbsum, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

      TH1F* h_model_dr_htm[6] ;
      h_model_dr_htm[1] = new TH1F( "h_model_dr_htm_nb0", "Model, double H/L ratio, Nb0, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_htm[2] = new TH1F( "h_model_dr_htm_nb1", "Model, double H/L ratio, Nb1, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_htm[3] = new TH1F( "h_model_dr_htm_nb2", "Model, double H/L ratio, Nb2, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_htm[4] = new TH1F( "h_model_dr_htm_nb3", "Model, double H/L ratio, Nb3, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;
      h_model_dr_htm[5] = new TH1F( "h_model_dr_htm_nbsum", "Model, double H/L ratio, Nbsum, HT high", nb_nj*nb_mht, 0.5, nb_nj*nb_mht+0.5 ) ;

      TH1F* h_model_dr_htl[6] ;
      h_model_dr_htl[1] = new TH1F( "h_model_dr_htl_nb0", "Model, double H/L ratio, Nb0, HT high", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_dr_htl[2] = new TH1F( "h_model_dr_htl_nb1", "Model, double H/L ratio, Nb1, HT high", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_dr_htl[3] = new TH1F( "h_model_dr_htl_nb2", "Model, double H/L ratio, Nb2, HT high", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_dr_htl[4] = new TH1F( "h_model_dr_htl_nb3", "Model, double H/L ratio, Nb3, HT high", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;
      h_model_dr_htl[5] = new TH1F( "h_model_dr_htl_nbsum", "Model, double H/L ratio, Nbsum, HT high", nb_nj*(nb_mht-2), 0.5, nb_nj*(nb_mht-2)+0.5 ) ;



      int bi_160(0) ;

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
               for ( int hti = htb_low; hti<=nb_ht; hti++ ) {

                  if ( mbi > 1 && nbi <= nb_nb ) bi_160 ++ ;

                  //if ( hti == 1 && nji > 2 ) continue ;

                  float par_val_mht(0.) ;
                  float par_err_mht(0.) ;
                  TH1F* h_ratio(0x0) ;
                  TH1F* h_dr(0x0) ;
                  int hbi(0) ;
                  if ( hti==3 ) {
                     par_val_mht = par_val_mht_hth[mbi] ;
                     par_err_mht = par_err_mht_hth[mbi] ;
                     h_ratio = h_model_ratio_hth[nbi] ;
                     h_dr = h_model_dr_hth[nbi] ;
                     hbi = (nji-1)*nb_mht + mbi ;
                  }
                  if ( hti==2 ) {
                     par_val_mht = par_val_mht_htm[mbi] ;
                     par_err_mht = par_err_mht_htm[mbi] ;
                     h_ratio = h_model_ratio_htm[nbi] ;
                     h_dr = h_model_dr_htm[nbi] ;
                     hbi = (nji-1)*nb_mht + mbi ;
                  }
                  if ( hti==1 ) {
                     par_val_mht = par_val_mht_htl[mbi] ;
                     par_err_mht = par_err_mht_htl[mbi] ;
                     h_ratio = h_model_ratio_htl[nbi] ;
                     h_dr = h_model_dr_htl[nbi] ;
                     hbi = (nji-1)*(nb_mht-2) + mbi ;
                  }

                  float model_ratio_val(0.) ;
                  float model_ratio_err2(0.) ;

                  model_ratio_val =  par_val_ht[hti] * par_val_njet[nji] * par_val_mht * par_val_nb[nbi] ;


                  if ( par_val_ht[hti] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_ht_fit[hti] / par_val_ht[hti], 2. ) ;
                  if ( par_val_ht[hti] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_ht_syst[hti] / par_val_ht[hti], 2. ) ;

                  if ( par_val_njet[nji] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_njet_fit[nji] / par_val_njet[nji], 2. ) ;
                  if ( par_val_njet[nji] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_njet_syst[nji] / par_val_njet[nji], 2. ) ;

                  if ( par_val_nb[nbi] > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_nb[nbi] / par_val_nb[nbi], 2. ) ;

                  if ( par_val_mht > 0 )   model_ratio_err2 += pow( model_ratio_val * par_err_mht / par_val_mht, 2. ) ;


                  float model_ratio_err = sqrt( model_ratio_err2 ) ;


                  float model_dr_val(0.) ;
                  float model_dr_err(0.) ;

                  model_dr_val = par_val_mht ;
                  model_dr_err = par_err_mht ;



                  h_ratio -> SetBinContent( hbi, model_ratio_val ) ;
                  h_ratio -> SetBinError( hbi, model_ratio_err ) ;
                  h_ratio -> GetXaxis() -> SetBinLabel( hbi, label ) ;

                  if ( mbi > 1 ) {
                     h_dr -> SetBinContent( hbi, model_dr_val ) ;
                     h_dr -> SetBinError( hbi, model_dr_err ) ;
                     h_dr -> GetXaxis() -> SetBinLabel( hbi, label ) ;
                  }

                  int print_bi(-1) ;
                  if ( mbi > 1 && nbi <= nb_nb ) {
                     print_bi = bi_160 ;
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





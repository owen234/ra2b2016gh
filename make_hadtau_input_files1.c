
   TH1* get_hist( TFile* tf, const char* hname ) ;

   void make_hadtau_input_files1( const char* input_root_file  = "non-qcd-inputs-topup2/ARElog58_12.9ifb_UpdatedSysAndMuSys_HadTauEstimation_data_formatted_New.root",
                                  const char* output_text_file = "outputfiles/combine-input-hadtau.txt",
                                  const char* nbsum_text_file  = "outputfiles/nbsum-input-hadtau.txt"
                               ) {

      gDirectory -> Delete( "h*" ) ;

      TFile* tf_input = new TFile( input_root_file, "read" ) ;
      if ( tf_input == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }
      if ( !(tf_input -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }

      printf("\n") ;
      tf_input -> ls() ;
      printf("\n") ;

      FILE* ofp_combine ;
      if ( (ofp_combine = fopen( output_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening combine output file: %s\n\n", output_text_file ) ;
         return ;
      }

      FILE* ofp_nbsum ;
      if ( (ofp_nbsum = fopen( nbsum_text_file, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum output file: %s\n\n", nbsum_text_file ) ;
         return ;
      }

      FILE* ofp_nbsum_stat_syst ;
      TString systfile_nbsum( nbsum_text_file ) ;
      systfile_nbsum.ReplaceAll("input","stat-syst") ;
      if ( (ofp_nbsum_stat_syst = fopen( systfile_nbsum.Data(), "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum stat-syst output file: %s\n\n", systfile_nbsum.Data() ) ;
         return ;
      }

      FILE* ofp_combine_stat_syst ;
      TString systfile_combine( output_text_file ) ;
      systfile_combine.ReplaceAll("input","stat-syst") ;
      if ( (ofp_combine_stat_syst = fopen( systfile_combine.Data(), "w" ))==NULL ) {
         printf( "\n\n *** Problem opening combine stat-syst output file: %s\n\n", systfile_combine.Data() ) ;
         return ;
      }


      ///////TH1F* h_ldp = (TH1F*) tf -> Get( "QCDBin_LowDphi_nominal_fullstatuncertainty" ) ;
      TH1* h_pred_lowdphi = get_hist( tf_input, "QCDBin_LowDphi_nominal_fullstatuncertainty" ) ;

      /////////TH1F* h_hdp = (TH1F*) tf -> Get( "QCDBin_HiDphi_nominal_fullstatuncertainty" ) ;
      TH1* h_pred_highdphi = get_hist( tf_input, "QCDBin_HiDphi_nominal_fullstatuncertainty" ) ;





      TH1* h_systerr_lowdphi[100] ;
      TH1* h_systerr_highdphi[100] ;
      char systerr_name[100][100] ;
      int n_systerr(0) ;

      {
         int si(0) ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_closureUncertainty" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_closureUncertainty" ) ;
         sprintf( systerr_name[si], "Closure" ) ;
         si++ ;

     //----------

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_BMistagUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_BMistagUp" ) ;
         sprintf( systerr_name[si], "Btag" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuRecoSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuRecoSysUp" ) ;
         sprintf( systerr_name[si], "MuRecoSys" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuIsoSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuIsoSysUp" ) ;
         sprintf( systerr_name[si], "MuIsoSys" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuRecoIsoUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuRecoIsoUp" ) ;
         sprintf( systerr_name[si], "MuRecoIso" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_JECSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_JECSysUp" ) ;
         sprintf( systerr_name[si], "JEC" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MTSysUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MTSysUp" ) ;
         sprintf( systerr_name[si], "MT" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MtEffStat" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MtEffStat" ) ;
         sprintf( systerr_name[si], "MtEffStat" ) ;
         si++ ;

    ///  h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_IsoTrkVetoEffUncertaintySys" ) ;
    ///  h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_IsoTrkVetoEffUncertaintySys" ) ;
    ///  sprintf( systerr_name[si], "IsoTrkVetoSys" ) ;
    ///  si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_AccStat" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_AccStat" ) ;
         sprintf( systerr_name[si], "AccStat" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_AccSysPDFUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_AccSysPDFUp" ) ;
         sprintf( systerr_name[si], "AccSysPDF" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_AccSysScaleUp" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_AccSysScaleUp" ) ;
         sprintf( systerr_name[si], "AccSysScale" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_MuFromTauStat" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_MuFromTauStat" ) ;
         sprintf( systerr_name[si], "MuFromTauStat" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_DileptonUncertainty" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_DileptonUncertainty" ) ;
         sprintf( systerr_name[si], "Dilepton" ) ;
         si++ ;

         h_systerr_lowdphi[si]  = get_hist( tf_input, "QCDBin_LowDphi_TrigEffUncertainty" ) ;
         h_systerr_highdphi[si] = get_hist( tf_input, "QCDBin_HiDphi_TrigEffUncertainty" ) ;
         sprintf( systerr_name[si], "TrigEff" ) ;
         si++ ;

         n_systerr = si ;
      }



      int nb_nj(4) ;
      int nb_nb(4) ;
      int nb_htmht(13) ;
      int nb_ht(3) ;

      int bi_hist(0) ;
      int bi_control(0) ;
      int bi_search(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               bi_hist ++ ;

               double ldp_val = h_pred_lowdphi -> GetBinContent( bi_hist ) ;
               double ldp_hist_err = h_pred_lowdphi -> GetBinError( bi_hist ) ;
               double ldp_err = 1. ;
               if ( ldp_val > 0 ) ldp_err = ldp_hist_err ;

               double hdp_val = h_pred_highdphi -> GetBinContent( bi_hist ) ;
               double hdp_hist_err = h_pred_highdphi -> GetBinError( bi_hist ) ;
               double hdp_err = 1. ;
               if ( hdp_val > 0 ) hdp_err = hdp_hist_err ;


               double total_syst2_lowdphi_events(0.) ;
               double total_syst2_highdphi_events(0.) ;
               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_lowdphi = fabs( h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                  double syst_highdphi = fabs( h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst2_lowdphi_events  += pow( syst_lowdphi  * ldp_val, 2. ) ;
                     total_syst2_highdphi_events += pow( syst_highdphi * hdp_val, 2. ) ;
                  }
               } // si

               double ldp_syst = sqrt( total_syst2_lowdphi_events ) ;
               double hdp_syst = sqrt( total_syst2_highdphi_events ) ;

               TString hist_bin_label( h_pred_lowdphi -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

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

               if ( bi_mht == 1 ) {
                  bi_control ++ ;
               } else {
                  bi_search ++ ;
               }

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

               //// printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;

               double ldp_stat_over_sqrtn(0.8), ldp_syst_over_n(0.15) ;
               if ( ldp_val > 0 ) {
                  ldp_stat_over_sqrtn = ldp_err / sqrt( ldp_val ) ;
                  ldp_syst_over_n = ldp_syst / ldp_val ;
                  if ( ldp_stat_over_sqrtn > 1.5 ) ldp_stat_over_sqrtn = 1.0 ;
               }
               double hdp_stat_over_sqrtn(0.8), hdp_syst_over_n(0.15) ;
               if ( hdp_val > 0 ) {
                  hdp_stat_over_sqrtn = hdp_err / sqrt( hdp_val ) ;
                  hdp_syst_over_n = hdp_syst / hdp_val ;
                  if ( hdp_stat_over_sqrtn > 1.5 ) hdp_stat_over_sqrtn = 1.0 ;
               }

               printf(               "%s      %8.1f +/- %5.1f +/- %5.1f          %8.1f +/- %5.1f +/- %5.1f\n",
                   label,
                   ldp_val, ldp_err, ldp_syst,
                   hdp_val, hdp_err, hdp_syst ) ;

               fprintf( ofp_combine, "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_err, ldp_syst,   hdp_val, hdp_err, hdp_syst ) ;

               if ( bi_htmht > 3 ) {
                  char this_label[100] ;
                  sprintf( this_label, "%3d S-Nj%d-Nb%d-MHT%d-HT%d (S%d)", bi_search, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht-3 ) ;
                  fprintf( ofp_combine_stat_syst, "   %-30s   %6.3f  %6.3f       %6.3f  %6.3f\n",
                      this_label,
                      ldp_stat_over_sqrtn, ldp_syst_over_n,   hdp_stat_over_sqrtn, hdp_syst_over_n ) ;
               } else {
                  char this_label[100] ;
                  sprintf( this_label, "%3d C-Nj%d-Nb%d-MHT%d-HT%d (C%d)", bi_control, bi_nj, bi_nb-1, bi_mht-1, bi_ht, bi_htmht ) ;
                  fprintf( ofp_combine_stat_syst, "   %-30s   %6.3f  %6.3f       %6.3f  %6.3f\n",
                      this_label,
                      ldp_stat_over_sqrtn, ldp_syst_over_n,   hdp_stat_over_sqrtn, hdp_syst_over_n ) ;
               }

            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp_combine ) ;
      printf("\n\n Wrote %s\n\n", output_text_file ) ;

     //-----

      for ( int bi_ht=1; bi_ht<=nb_ht; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            float ldp_nbsum_val(0.) ;
            float hdp_nbsum_val(0.) ;

            float ldp_nbsum_stat_err2(0.) ;
            float hdp_nbsum_stat_err2(0.) ;

            double total_syst_lowdphi_events[100] ;
            double total_syst_highdphi_events[100] ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_lowdphi_events[si]  = 0 ;
               total_syst_highdphi_events[si] = 0 ;
            }

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               bi_hist = (bi_nj-1)*(nb_nb)*(nb_htmht) + (bi_nb-1)*(nb_htmht) + bi_ht ;

               double ldp_val = h_pred_lowdphi -> GetBinContent( bi_hist ) ;
               double ldp_hist_err = h_pred_lowdphi -> GetBinError( bi_hist ) ;
               double ldp_err = 1. ;
               if ( ldp_val > 0 ) ldp_err = ldp_hist_err ;

               double hdp_val = h_pred_highdphi -> GetBinContent( bi_hist ) ;
               double hdp_hist_err = h_pred_highdphi -> GetBinError( bi_hist ) ;
               double hdp_err = 1. ;
               if ( hdp_val > 0 ) hdp_err = hdp_hist_err ;


               //////////////double ldp_syst = 0.10 * ldp_val ; //*** flat 10% syst for now.
               //////////////double hdp_syst = 0.10 * hdp_val ; //*** flat 10% syst for now.

               for ( int si=0; si<n_systerr; si++ ) {
                  double syst_lowdphi = fabs( h_systerr_lowdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                  double syst_highdphi = fabs( h_systerr_highdphi[si] -> GetBinContent( bi_hist ) - 1. ) ;
                  if ( syst_lowdphi > 0 && syst_highdphi > 0 ) {
                     total_syst_lowdphi_events[si]  += syst_lowdphi  * ldp_val ;
                     total_syst_highdphi_events[si] += syst_highdphi * hdp_val ;
                  }
               } // si

               ldp_nbsum_val += ldp_val ;
               hdp_nbsum_val += hdp_val ;

               ldp_nbsum_stat_err2 += pow( ldp_err, 2. ) ;
               hdp_nbsum_stat_err2 += pow( hdp_err, 2. ) ;

               TString hist_bin_label( h_pred_lowdphi -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               char label[1000] ;
               sprintf( label, " %3d  Nj%d-Nb%d-MHTC-HT%d", bi_hist, bi_nj, bi_nb-1, bi_ht ) ;

            } // bi_nb

            ///printf("    Systematics:\n") ;
            double total_syst_err2_lowdphi(0.) ;
            double total_syst_err2_highdphi(0.) ;
            for ( int si=0; si<n_systerr; si++ ) {
               total_syst_err2_lowdphi += pow( total_syst_lowdphi_events[si], 2.)   ;
               total_syst_err2_highdphi += pow( total_syst_highdphi_events[si], 2.)  ;
               ///printf("  %2d : %25s :   %7.1f  %7.1f\n", si, systerr_name[si], total_syst_lowdphi_events[si], total_syst_highdphi_events[si] ) ;
            }
            double total_syst_err_lowdphi = sqrt( total_syst_err2_lowdphi ) ;
            double total_syst_err_highdphi = sqrt( total_syst_err2_highdphi ) ;


            float ldp_nbsum_stat_err = sqrt( ldp_nbsum_stat_err2 ) ;
            float hdp_nbsum_stat_err = sqrt( hdp_nbsum_stat_err2 ) ;

            double ldp_syst_over_n(0.) ;
            double ldp_stat_over_sqrtn(0.) ;
            if ( ldp_nbsum_val > 0 ) ldp_syst_over_n = total_syst_err_lowdphi / ldp_nbsum_val  ;
            if ( ldp_nbsum_val > 0 ) ldp_stat_over_sqrtn = ldp_nbsum_stat_err / sqrt(ldp_nbsum_val)  ;

            double hdp_syst_over_n(0.) ;
            double hdp_stat_over_sqrtn(0.) ;
            if ( hdp_nbsum_val > 0 ) hdp_syst_over_n = total_syst_err_highdphi / hdp_nbsum_val  ;
            if ( hdp_nbsum_val > 0 ) hdp_stat_over_sqrtn = hdp_nbsum_stat_err / sqrt(hdp_nbsum_val)  ;

            printf(          "   Nj%d-HT%d      %8.1f +/- %5.1f +/- %5.1f  (%6.3f)        %8.1f +/- %5.1f +/- %5.1f   (%6.3f)\n",
                bi_nj, bi_ht,
                ldp_nbsum_val, ldp_nbsum_stat_err, total_syst_err_lowdphi,    ldp_syst_over_n,
                hdp_nbsum_val, hdp_nbsum_stat_err, total_syst_err_highdphi,   hdp_syst_over_n ) ;

            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f    %8.1f +/- %5.1f +/- %5.1f\n", bi_nj, bi_ht,
                  ldp_nbsum_val, ldp_nbsum_stat_err, total_syst_err_lowdphi,
                  hdp_nbsum_val, hdp_nbsum_stat_err, total_syst_err_highdphi   ) ;
            fprintf( ofp_nbsum_stat_syst, "   Nj%d-HT%d    %6.3f  %6.3f       %6.3f  %6.3f\n", bi_nj, bi_ht, 
                  ldp_stat_over_sqrtn, ldp_syst_over_n,    hdp_stat_over_sqrtn, hdp_syst_over_n ) ;

         } // bi_nj
      } // bi_ht

      fclose( ofp_nbsum ) ;
      fclose( ofp_nbsum_stat_syst ) ;
      printf("\n\n Wrote %s\n\n", nbsum_text_file ) ;
      printf("\n\n Wrote %s\n\n", systfile_nbsum.Data() ) ;
      printf("\n\n Wrote %s\n\n", systfile_combine.Data() ) ;



   } // make_hadtau_input_files1




  //========================================================================================================

   TH1* get_hist( TFile* tf, const char* hname ) {
      TH1* hp   = (TH1*) tf -> Get( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** Can't find hist %s.\n\n", hname ) ;
         gSystem -> Exit(-1) ;
      }
      return hp ;
   } // get_hist

  //========================================================================================================





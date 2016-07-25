

   void make_signalmc_input_files1( const char* input_root_file = "outputfiles/hists-v2d-T1bbbbH.root",
                                    const char* output_text_file = "outputfiles/combine-input-T1bbbbH.txt",
                                    const char* nbsum_text_file = "outputfiles/nbsum-input-T1bbbbH.txt"
                               ) {

      gDirectory -> Delete( "h*" ) ;

      TFile* tf = new TFile( input_root_file, "read" ) ;
      if ( tf == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }
      if ( !(tf -> IsOpen() ) ) { printf("\n\n *** Bad input file: %s\n\n", input_root_file ) ; return ; }

      printf("\n") ;
      tf -> ls() ;
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

      TH1F* h_ldp = (TH1F*) tf -> Get( "h_ldp" ) ;
      if ( h_ldp == 0x0 ) { printf("\n\n *** Missing h_ldp\n\n") ; return ; }

      TH1F* h_hdp = (TH1F*) tf -> Get( "h_hdp" ) ;
      if ( h_hdp == 0x0 ) { printf("\n\n *** Missing h_hdp\n\n") ; return ; }

      int nb_nj(4) ;
      int nb_nb(4) ;
      int nb_htmht(13) ;
      int nb_ht(3) ;

      int bi_hist(0) ;


      for ( int bi_ht=1; bi_ht<=nb_ht; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            double ldp_nbsum_val(0.) ;
            double ldp_nbsum_err2(0.) ;

            double hdp_nbsum_val(0.) ;
            double hdp_nbsum_err2(0.) ;

            for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {

               bi_hist = (bi_nj-1)*(nb_nb)*(nb_htmht) + (bi_nb-1)*(nb_htmht) + bi_ht ;

               double ldp_val = h_ldp -> GetBinContent( bi_hist ) ;

               double hdp_val = h_hdp -> GetBinContent( bi_hist ) ;

               ldp_nbsum_val += ldp_val ;
               hdp_nbsum_val += hdp_val ;

               TString hist_bin_label( h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               char label[1000] ;
               sprintf( label, " %3d  Nj%d-Nb%d-MHTC-HT%d", bi_hist, bi_nj, bi_nb-1, bi_ht ) ;

       ////    printf("  label : %s   ,  hist label %s\n", label, hist_bin_label.Data() ) ;

            } // bi_nb

            printf( "   Nj%d-HT%d   %8.1f      %8.1f\n", 
                bi_nj, bi_ht, ldp_nbsum_val, hdp_nbsum_val ) ;

            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f      %8.1f\n", 
                bi_nj, bi_ht, ldp_nbsum_val, hdp_nbsum_val ) ;


         } // bi_nj
      } // bi_ht


      fclose( ofp_nbsum ) ;
      printf("\n\n Wrote %s\n\n", nbsum_text_file ) ;





     //----------
      bi_hist = 0 ;
      int bi_control(0) ;
      int bi_search(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {

               bi_hist ++ ;

               double ldp_val = h_ldp -> GetBinContent( bi_hist ) ;

               double hdp_val = h_hdp -> GetBinContent( bi_hist ) ;


               TString hist_bin_label( h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

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


               printf(               "%s      %8.1f          %8.1f\n", label,    ldp_val,    hdp_val ) ;

               fprintf( ofp_combine, "%s      %8.1f          %8.1f\n", label,    ldp_val,    hdp_val ) ;


            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp_combine ) ;
      printf("\n\n Wrote %s\n\n", output_text_file ) ;

     //---------


   } // make_signalmc_input_files1




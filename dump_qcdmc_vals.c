


#include "histio.c"


   void dump_qcdmc_vals(
         const char* infile = "outputfiles/hists-v2d-qcd.root",
         const char* outfile = "outputfiles/qcdmc-counts.txt"
         ) {

      gDirectory -> Delete( "h*" ) ;

      loadHist( infile ) ;

      TH1F* h_ldp = (TH1F*) gDirectory -> FindObject( "h_ldp" ) ;
      if ( h_ldp == 0x0 ) { printf("\n\n *** can't find h_ldp in %s\n\n", infile ) ; return ; }

      TH1F* h_hdp = (TH1F*) gDirectory -> FindObject( "h_hdp" ) ;
      if ( h_hdp == 0x0 ) { printf("\n\n *** can't find h_hdp in %s\n\n", infile ) ; return ; }

      FILE* ofp ;
      if ( (ofp = fopen( outfile, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening output file: %s\n\n", outfile ) ;
         return ;
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
               if ( bi_htmht <= 3 ) {
                  bi_control ++ ;
               } else {
                  bi_search ++ ;
               }

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

               printf( " %30s    %12.3f  %12.3f      %12.3f  %12.3f\n", label,
                  h_ldp -> GetBinContent( bi_hist ),
                  h_ldp -> GetBinError( bi_hist ),
                  h_hdp -> GetBinContent( bi_hist ),
                  h_hdp -> GetBinError( bi_hist )
                  ) ;

               fprintf( ofp,  " %30s    %12.3f  %12.3f      %12.3f  %12.3f\n", label,
                  h_ldp -> GetBinContent( bi_hist ),
                  h_ldp -> GetBinError( bi_hist ),
                  h_hdp -> GetBinContent( bi_hist ),
                  h_hdp -> GetBinError( bi_hist )
                  ) ;



            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp ) ;
      printf("\n\n Wrote %s\n\n", outfile ) ;


   } // dump_qcdmc_vals





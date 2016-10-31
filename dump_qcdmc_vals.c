#ifndef dump_qcdmc_vals_c
#define dump_qcdmc_vals_c

#include "TSystem.h"

#include "histio.c"
#include "binning.h"

   void dump_qcdmc_vals(
         const char* infile = "outputfiles/hists-v2d-qcd.root",
         const char* outfile = "outputfiles/qcdmc-counts.txt"
         ) {

      setup_bins();
      
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

               int bi_ht = 0 , bi_mht = 0 ;

               htmht_bin_to_ht_and_mht_bins ( bi_htmht, bi_ht, bi_mht );

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

#endif

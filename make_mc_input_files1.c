#ifndef make_mc_input_files1_c
#define make_mc_input_files1_c

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"

#include <fstream>

#include "binning.h"

#include "histio.c"

   TH1F* get_hist( const char* hname ) ;

   void make_mc_input_files1(
                              const char* sample_name = "lostlep",
                              const char* iodir  = "outputfiles/"
                             ) {

      setup_bins();
      gDirectory -> Delete( "h*" ) ;

      char fname[10000] ;
      TString line ;

      sprintf( fname, "%s/hists-v2d-%s.root", iodir, sample_name ) ;
      loadHist( fname ) ;

      printf("\n") ;
      gDirectory -> ls() ;
      printf("\n") ;

      sprintf( fname, "%s/mc-combine-input-%s.txt", iodir, sample_name ) ;
      FILE* ofp_combine ;
      if ( (ofp_combine = fopen( fname, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening combine output file: %s\n\n", fname ) ;
         return ;
      }

      sprintf( fname, "%s/mc-nbsum-input-%s.txt", iodir, sample_name ) ;
      FILE* ofp_nbsum ;
      if ( (ofp_nbsum = fopen( fname, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening nbsum output file: %s\n\n", fname ) ;
         return ;
      }

     ////--- Have to wait until the other BG teams give new inputs with Njets=2 for this part

///   ifstream ifs_combine_stat_syst ;
///   sprintf( fname, "%s/combine-stat-syst-%s.txt", iodir, sample_name ) ;
///   ifs_combine_stat_syst.open( fname ) ;
///   if ( !ifs_combine_stat_syst.good() ) { printf("\n\n *** Can't open %s\n\n", fname ) ; return ; }

///   ifstream ifs_nbsum_stat_syst ;
///   sprintf( fname, "%s/nbsum-stat-syst-%s.txt", iodir, sample_name ) ;
///   ifs_nbsum_stat_syst.open( fname ) ;
///   if ( !ifs_nbsum_stat_syst.good() ) { printf("\n\n *** Can't open %s\n\n", fname ) ; return ; }


      TH1F* h_ldp = get_hist( "h_ldp" ) ;
      TH1F* h_hdp = get_hist( "h_hdp" ) ;

      int bi_hist(0) ;
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


               float ldp_stat_over_sqrtn = 0.7 ; // guess until input below is ready including Njets=2.
               float ldp_syst_over_n = 0.15 ; // guess until input below is ready including Njets=2.
               float hdp_stat_over_sqrtn = 0.7 ; // guess until input below is ready including Njets=2.
               float hdp_syst_over_n = 0.15 ;  // guess until input below is ready including Njets=2.

      /////    line.ReadLine( ifs_combine_stat_syst ) ;
      /////    int cssbi ;
      /////    char csslabel1[100] ;
      /////    char csslabel2[100] ;
      /////    float ldp_stat_over_sqrtn(1.), ldp_syst_over_n(1.),   hdp_stat_over_sqrtn(1.), hdp_syst_over_n(1.) ;
      /////    sscanf( line.Data(), " %d %s %s  %f %f  %f %f", &cssbi, csslabel1, csslabel2,
      /////      &ldp_stat_over_sqrtn, &ldp_syst_over_n, &hdp_stat_over_sqrtn, &hdp_syst_over_n ) ;
      /////    if ( bi_htmht <=3 && cssbi != bi_control ) {
      /////       printf("\n\n *1* Inconsistency in input file bins:  %d != %d, %s is not %s\n\n",
      /////          cssbi, bi_control, label, csslabel1 ) ;
      /////       return ;
      /////    }
      /////    if ( bi_htmht >3 && cssbi != bi_search ) {
      /////       printf("\n\n *2* Inconsistency in input file bins:  %d != %d, %s is not %s\n\n",
      /////          cssbi, bi_search, label, csslabel1 ) ;
      /////       return ;
      /////    }


               double ldp_stat = ldp_stat_over_sqrtn * sqrt( ldp_val ) ;
               double ldp_syst = ldp_syst_over_n * ldp_val  ;

               double hdp_stat = hdp_stat_over_sqrtn * sqrt( hdp_val ) ;
               double hdp_syst = hdp_syst_over_n * hdp_val  ;

               printf(               "%s      %8.1f +/- %5.1f +/- %5.1f          %8.1f +/- %5.1f +/- %5.1f\n",
                   label,
                   ldp_val, ldp_stat, ldp_syst,
                   hdp_val, hdp_stat, hdp_syst ) ;

               fprintf( ofp_combine, "%s      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                   label,    ldp_val, ldp_stat, ldp_syst,   hdp_val, hdp_stat, hdp_syst ) ;



            } // bi_htmht
         } // bi_nb
      } // bi_nj

      fclose( ofp_combine ) ;

     //-----

      printf("\n\n nb_ht[1] = %d\n\n", nb_ht[1] ) ;
      for ( int bi_ht=1; bi_ht<=nb_ht[1]; bi_ht++ ) {
         for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {

            float ldp_nbsum_val(0.) ;
            float hdp_nbsum_val(0.) ;

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

             float ldp_stat_over_sqrtn = 0.7 ; // guess until input below is ready including Njets=2.
             float ldp_syst_over_n = 0.15 ; // guess until input below is ready including Njets=2.
             float hdp_stat_over_sqrtn = 0.7 ; // guess until input below is ready including Njets=2.
             float hdp_syst_over_n = 0.15 ;  // guess until input below is ready including Njets=2.

    //////  line.ReadLine( ifs_nbsum_stat_syst ) ;
    //////  char csslabel1[100] ;
    //////  float ldp_stat_over_sqrtn(1.), ldp_syst_over_n(1.),   hdp_stat_over_sqrtn(1.), hdp_syst_over_n(1.) ;
    //////  sscanf( line.Data(), " %s  %f %f  %f %f", csslabel1,
    //////    &ldp_stat_over_sqrtn, &ldp_syst_over_n, &hdp_stat_over_sqrtn, &hdp_syst_over_n ) ;

            double ldp_nbsum_stat_err = ldp_stat_over_sqrtn * sqrt( ldp_nbsum_val ) ;
            double ldp_nbsum_syst_err = ldp_syst_over_n * ldp_nbsum_val  ;

            double hdp_nbsum_stat_err = hdp_stat_over_sqrtn * sqrt( hdp_nbsum_val ) ;
            double hdp_nbsum_syst_err = hdp_syst_over_n * hdp_nbsum_val  ;


            printf(          "   Nj%d-HT%d      %8.1f +/- %5.1f +/- %5.1f         %8.1f +/- %5.1f +/- %5.1f\n",
                bi_nj, bi_ht,
                ldp_nbsum_val, ldp_nbsum_stat_err, ldp_nbsum_syst_err,
                hdp_nbsum_val, hdp_nbsum_stat_err, hdp_nbsum_syst_err  ) ;

            fprintf( ofp_nbsum, "   Nj%d-HT%d   %8.1f +/- %5.1f +/- %5.1f    %8.1f +/- %5.1f +/- %5.1f\n", bi_nj, bi_ht,
                  ldp_nbsum_val, ldp_nbsum_stat_err, ldp_nbsum_syst_err,
                  hdp_nbsum_val, hdp_nbsum_stat_err, hdp_nbsum_syst_err ) ;

         } // bi_nj
      } // bi_ht

      fclose( ofp_nbsum ) ;



   } // make_hadtau_input_files1



//===============================================================================

   TH1F* get_hist( const char* hname ) {
      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** Missing histogram : %s\n\n", hname ) ;
         gDirectory -> ls() ;
         gSystem -> Exit( -1 ) ;
      }
      return hp ;
   } // get_hist

//===============================================================================

#endif


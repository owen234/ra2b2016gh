
#include "TString.h"
#include "TMath.h"

#include <fstream>


   void make_fakedata_input_files1(
                const char* input_prefix = "outputfiles/mc-nbsum-input"
          ) {

      char fname[1000] ;

      ifstream ifs_hadtau ;
      ifstream ifs_lostlep ;
      ifstream ifs_znunu ;
      ifstream ifs_qcd ;

      sprintf( fname, "%s-hadtau.txt", input_prefix ) ;
      ifs_hadtau.open( fname ) ;
      if ( !ifs_hadtau.good() ) { printf("\n\n *** Bad input file: %s\n\n", fname ) ; return ; }

      sprintf( fname, "%s-lostlep.txt", input_prefix ) ;
      ifs_lostlep.open( fname ) ;
      if ( !ifs_lostlep.good() ) { printf("\n\n *** Bad input file: %s\n\n", fname ) ; return ; }

      sprintf( fname, "%s-znunu.txt", input_prefix ) ;
      ifs_znunu.open( fname ) ;
      if ( !ifs_znunu.good() ) { printf("\n\n *** Bad input file: %s\n\n", fname ) ; return ; }

      sprintf( fname, "%s-qcd.txt", input_prefix ) ;
      ifs_qcd.open( fname ) ;
      if ( !ifs_qcd.good() ) { printf("\n\n *** Bad input file: %s\n\n", fname ) ; return ; }


      sprintf( fname, "%s-fakedata.txt", input_prefix ) ;
      FILE* ofp ;
      if ( ( ofp = fopen( fname, "w" ) ) == NULL ) {
         printf( "\n\n *** Problem opening output file: %s\n\n", fname ) ; return ;
      }

      TString line ;

      while ( ifs_hadtau.good() ) {

         char label[100] ;
         float ldp_val, ldp_stat, ldp_syst ;
         float hdp_val, hdp_stat, hdp_syst ;

         float ldp_sum(0.) ;
         float hdp_sum(0.) ;

         line.ReadLine( ifs_hadtau ) ;
         if ( !ifs_hadtau.good() ) break ;
         sscanf( line.Data(), "%s %f +/- %f +/- %f  %f +/- %f +/- %f",
             label, &ldp_val, &ldp_stat, &ldp_syst,  &hdp_val, &hdp_stat, &hdp_syst ) ;

         printf("  %s  hadtau   LDP=%8.1f  HDP=%8.1f\n", label, ldp_val, hdp_val ) ;
         ldp_sum += ldp_val ;
         hdp_sum += hdp_val ;


         line.ReadLine( ifs_lostlep ) ;
         sscanf( line.Data(), "%s %f +/- %f +/- %f  %f +/- %f +/- %f",
             label, &ldp_val, &ldp_stat, &ldp_syst,  &hdp_val, &hdp_stat, &hdp_syst ) ;

         printf("  %s  lostlep  LDP=%8.1f  HDP=%8.1f\n", label, ldp_val, hdp_val ) ;
         ldp_sum += ldp_val ;
         hdp_sum += hdp_val ;


         line.ReadLine( ifs_znunu ) ;
         sscanf( line.Data(), "%s %f +/- %f +/- %f  %f +/- %f +/- %f",
             label, &ldp_val, &ldp_stat, &ldp_syst,  &hdp_val, &hdp_stat, &hdp_syst ) ;

         printf("  %s  znunu    LDP=%8.1f  HDP=%8.1f\n", label, ldp_val, hdp_val ) ;
         ldp_sum += ldp_val ;
         hdp_sum += hdp_val ;


         line.ReadLine( ifs_qcd ) ;
         sscanf( line.Data(), "%s %f +/- %f +/- %f  %f +/- %f +/- %f",
             label, &ldp_val, &ldp_stat, &ldp_syst,  &hdp_val, &hdp_stat, &hdp_syst ) ;

         printf("  %s  QCD      LDP=%8.1f  HDP=%8.1f\n", label, ldp_val, hdp_val ) ;
         ldp_sum += ldp_val ;
         hdp_sum += hdp_val ;

         int fakedata_ldp = TMath::Nint( ldp_sum ) ;
         int fakedata_hdp = TMath::Nint( hdp_sum ) ;
         printf("  %s  fakedata LDP=%8d  HDP=%8d\n\n", label, fakedata_ldp, fakedata_hdp ) ;

         fprintf( ofp, "%s    %8d    %8d\n", label, fakedata_ldp, fakedata_hdp ) ;

      }

      fclose( ofp ) ;


   } // make_fakedata_input_files1



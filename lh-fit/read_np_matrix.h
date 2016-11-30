
#include "TSystem.h"
#include "TString.h"


   int number_of_np ;
   int np_mat_row_index_nj_ht[20][20] ;
   double np_mat[20][20] ;

   void read_np_matrix( const char* infile ) {

      ifstream ifs ;
      ifs.open( infile ) ;
      if ( !ifs.good() ) { printf("\n\n *** read_np_matrix : Bad input file: %s\n\n", infile ) ; gSystem -> Exit(-1) ; }

      TString token ;

      token.ReadToken( ifs ) ;
      sscanf( token.Data(), "%d", &number_of_np ) ;
      printf("  Number of nuisance parameters: %d\n", number_of_np ) ;

      for ( int npi=0; npi<number_of_np; npi++ ) {
         token.ReadToken( ifs ) ;
         int bi_nj, bi_ht ;
         sscanf( token.Data(), "Nj%d_HT%d", &bi_nj, &bi_ht ) ;
         if ( bi_nj <= 0 ) { printf("\n\n *** read_np_matrix : illegal Nj index %d\n\n", bi_nj ) ; gSystem -> Exit(-1) ; }
         if ( bi_ht <= 0 ) { printf("\n\n *** read_np_matrix : illegal HT index %d\n\n", bi_ht ) ; gSystem -> Exit(-1) ; }
         np_mat_row_index_nj_ht[bi_nj][bi_ht] = npi ;
         printf("  %s  : bi_nj=%d, bi_ht=%d,  nuisance parameter index %2d\n", token.Data(), bi_nj, bi_ht, npi ) ;
      } // npi

      for ( int rowi=0; rowi<number_of_np; rowi++ ) {
         token.ReadToken( ifs ) ;
         int bi_nj, bi_ht ;
         sscanf( token.Data(), "Nj%d_HT%d", &bi_nj, &bi_ht ) ;
         if ( bi_nj <= 0 ) { printf("\n\n *** read_np_matrix : illegal Nj index %d\n\n", bi_nj ) ; gSystem -> Exit(-1) ; }
         if ( bi_ht <= 0 ) { printf("\n\n *** read_np_matrix : illegal HT index %d\n\n", bi_ht ) ; gSystem -> Exit(-1) ; }
         if ( rowi != np_mat_row_index_nj_ht[bi_nj][bi_ht] ) { printf("\n\n *** read_np_matrix : inconsistency, row %d\n", rowi ) ; gSystem -> Exit(-1) ; }
         printf( "  %s  : ", token.Data() ) ;
         for ( int coli=0; coli<number_of_np; coli++ ) {
            token.ReadToken( ifs ) ;
            double rel_err ;
            sscanf( token.Data(), "%lf", &rel_err ) ;
            np_mat[rowi][coli] = rel_err ;
            printf( " %6.4f ", np_mat[rowi][coli] ) ;
         }
         printf("\n") ;
         token.ReadToken( ifs ) ; // read the extra column
      } // rowi



   } // read_np_matrix


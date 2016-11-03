#ifndef read_pars_h 
#define read_pars_h

#include <fstream>
#include "TSystem.h"

#include "binning.h"

float par_val_ht       [10] ;
float par_err_ht_fit   [10] ;
float par_err_ht_syst  [10] ;

float par_val_njet     [10] ;
float par_err_njet_fit [10] ;
float par_err_njet_syst[10] ;

float par_val_ht_mht[10][10] ; // first index for ht and second one for mht 
float par_err_ht_mht[10][10] ; // first index for ht and second one for mht 

float par_val_nb       [10] ;
float par_err_nb       [10] ;

float par_rel_err_ht   [10] ;
float par_rel_err_njet [10] ;
float par_rel_err_ht_mht[10][10] ;
float par_rel_err_nb   [10] ;


void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2, float& rel_err );

void read_pars( const char* model_pars_file ) {

   ifstream ifs_model_pars ;
   ifs_model_pars.open( model_pars_file ) ;
   if ( !ifs_model_pars.good() ) { printf("\n\n *** Problem opening %s\n\n", model_pars_file ) ; gSystem->Exit(-1) ; }

   float val, err1, err2, rel_err ;
   char pname[100] ;

   //---

   for ( int bi_ht = 1; bi_ht<=nBinsHT; bi_ht++ )
   {
      sprintf( pname, "Kqcd_HT%d", bi_ht ) ;
      get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
      par_val_ht     [bi_ht] = val ;
      par_err_ht_fit [bi_ht] = err1 ;
      par_err_ht_syst[bi_ht] = val*err2 ;
      par_rel_err_ht [bi_ht] = rel_err ;
      printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;
   }//bi_ht

       //---

   for ( int bi_njet=1; bi_njet<=nb_nj; bi_njet++ )
   {
      sprintf( pname, "Sqcd_njet%d", bi_njet ) ;
      get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
      par_val_njet     [bi_njet] = val ;
      par_err_njet_fit [bi_njet] = err1 ;
      par_err_njet_syst[bi_njet] = val*err2 ;
      par_rel_err_njet [bi_njet] = rel_err ;
      printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;
   }//bi_njet
       //---

   for ( int bi_ht = nBinsHT; bi_ht>=1; bi_ht-- )
      for ( int bi_mht=0; bi_mht<nb_mht; bi_mht++ )
      {
         char ht_level [10], mhtc[10] = "c";
         if ( bi_ht == 1 ) strcpy(ht_level, "htl");
         if ( bi_ht == 2 ) strcpy(ht_level, "htm");
         if ( bi_ht == 3 ) strcpy(ht_level, "hth");

         if ( bi_ht == 1 && bi_mht > 2 ) continue;

         if ( bi_mht == 0 ) sprintf( pname, "Sqcd_mht%s_%s",mhtc  , ht_level ) ;
         else               sprintf( pname, "Sqcd_mht%d_%s",bi_mht, ht_level ) ;

         get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
         par_val_ht_mht    [bi_ht][bi_mht+1] = val ;
         par_err_ht_mht    [bi_ht][bi_mht+1] = sqrt( err1*err1 + val*val*err2*err2 ) ;
         par_rel_err_ht_mht[bi_ht][bi_mht+1] = rel_err ;
         printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;
      }//bi_mht and bi_ht


   ///

   for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ )
   {
      sprintf( pname, "Sqcd_nb%d", bi_nb-1 ) ;
      get_par( ifs_model_pars, pname, val, err1, err2, rel_err ) ;
      par_val_nb    [bi_nb] = val ;
      par_err_nb    [bi_nb] = sqrt ( err1*err1 + val*val*err2*err2 ) ;
      par_rel_err_nb[bi_nb] = rel_err ;
      printf("   Read %s : %.4f %.4f %.4f\n", pname, val, err1, err2 ) ;
   }//bi_nb

   } // read_pars

  //=======================================================================================

   void get_par( ifstream& ifs, const char* pname, float& val, float& err1, float& err2, float& rel_err ) {

      val = 0. ;
      err1 = 0. ;
      err2 = 0. ;
      rel_err = 0. ;

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
            if ( val > 0 ) {
               float err_total = sqrt( err1*err1 + val*err2 * val*err2 ) ;
               rel_err = err_total / val ;
            }
            return ;
         }
      }

      printf("\n\n *** get_par : Failed to find parameter %s\n\n", pname ) ;
      gSystem -> Exit(-1) ;
   } // get_par
#endif
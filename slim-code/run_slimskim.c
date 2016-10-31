
#include "TSystem.h"
#include "TString.h"
#include <fstream>

#include "doSkimSlim.c"

   void run_slimskim( const char* thedir = "/data/strange2/owen/fnal-prod-v9-skims/tree_LDP", const char* wildcard="*.root", bool blind=false ) {

      char command[10000] ;

      sprintf( command, "ls -1 %s/%s > rootfiles.txt", thedir, wildcard ) ;
      gSystem -> Exec( command ) ;

      char outdir[10000] ;
      sprintf( outdir, "%s/slim", thedir ) ;


      printf("\n\n") ;

      ifstream ifs ;
      ifs.open( "rootfiles.txt" ) ;
      int fi(0) ;
      while ( ifs.good() ) {
         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( !ifs.good() ) break ;
         doSkimSlim( ts.Data(), outdir, true, true, blind ) ;
         fi ++ ;
      }

      printf("\n\n") ;



   } // run_slimskim

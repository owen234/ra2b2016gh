


#include "do_ttbar_filter_skim.c"


   void run_ttbar_filter_skims( const char* io_dir = "/data/strange2/owen/fnal-prod-v9-skims" ) {

      char infile[10000] ;
      char outfile[10000] ;

      sprintf( infile,  "%s/tree_signal/slim/tree_TTJets_SingleLeptFromT-slimskim.root", io_dir ) ;
      sprintf( outfile, "%s/tree_signal/slim/tree_TTJets_SingleLeptFromT_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 1, infile, outfile ) ;

      sprintf( infile,  "%s/tree_signal/slim/tree_TTJets_SingleLeptFromTbar-slimskim.root", io_dir ) ;
      sprintf( outfile, "%s/tree_signal/slim/tree_TTJets_SingleLeptFromTbar_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 1, infile, outfile ) ;

      sprintf( infile,  "%s/tree_signal/slim/tree_TTJets_DiLept-slimskim.root", io_dir );
      sprintf( outfile, "%s/tree_signal/slim/tree_TTJets_DiLept_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 2, infile, outfile ) ;

      sprintf( infile,  "%s/tree_signal/slim/tree_TTJets-slimskim.root", io_dir ) ;
      sprintf( outfile, "%s/tree_signal/slim/tree_TTJets_hadonly_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 3, infile, outfile ) ;



      sprintf( infile,  "%s/tree_LDP/slim/tree_TTJets_SingleLeptFromT-slimskim.root", io_dir ) ;
      sprintf( outfile, "%s/tree_LDP/slim/tree_TTJets_SingleLeptFromT_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 1, infile, outfile ) ;

      sprintf( infile,  "%s/tree_LDP/slim/tree_TTJets_SingleLeptFromTbar-slimskim.root", io_dir ) ;
      sprintf( outfile, "%s/tree_LDP/slim/tree_TTJets_SingleLeptFromTbar_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 1, infile, outfile ) ;

      sprintf( infile,  "%s/tree_LDP/slim/tree_TTJets_DiLept-slimskim.root", io_dir );
      sprintf( outfile, "%s/tree_LDP/slim/tree_TTJets_DiLept_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 2, infile, outfile ) ;

      sprintf( infile,  "%s/tree_LDP/slim/tree_TTJets-slimskim.root", io_dir ) ;
      sprintf( outfile, "%s/tree_LDP/slim/tree_TTJets_hadonly_ght_lt600-slimskim.root", io_dir ) ;
      do_ttbar_filter_skim( 3, infile, outfile ) ;

   }

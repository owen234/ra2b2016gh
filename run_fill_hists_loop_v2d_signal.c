
#include "TChain.h"

#include "fill_hists_loop_v2d.c"

   void run_fill_hists_loop_v2d_signal( const char* indir = "fnal-prod-v9-skims-slimmed" ) {

      char fpat[10000] ;
      char sample_name[100] ;
      TChain* ch ;
      int n_added ;
      fill_hists_loop_v2d* fhl ;


   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "T1bbbbH" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_T1bbbb_1500_100-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_T1bbbb_1500_100-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;


   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "T1bbbbC" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_T1bbbb_1000_900-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_T1bbbb_1000_900-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;


   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "T1qqqqH" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_T1qqqq_1400_100-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_T1qqqq_1400_100-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;


   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "T1qqqqC" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_T1qqqq_1000_800-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_T1qqqq_1000_800-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;


   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "T1ttttH" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_T1tttt_1500_100-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_T1tttt_1500_100-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;

   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "T1ttttC" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_T1tttt_1200_800-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_T1tttt_1200_800-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;


   //----------------------



   } // run_fill_hists_loop_v2d_signal



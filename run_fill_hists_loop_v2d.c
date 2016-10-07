
#include "TChain.h"

#include "fill_hists_loop_v2d.c"

   void run_fill_hists_loop_v2d( const char* indir = "fnal-prod-v9-skims-slimmed" ) {

      char fpat[10000] ;
      char sample_name[100] ;
      TChain* ch ;
      int n_added ;
      fill_hists_loop_v2d* fhl ;

   //----------------------

///   ch = new TChain("tree") ;
///   sprintf( sample_name, "qcd" ) ;

///   sprintf( fpat, "%s/tree_LDP/tree_QCD*.root", indir ) ;
///   n_added = ch->Add(fpat) ;
///   printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

///   sprintf( fpat, "%s/tree_signal/tree_QCD*.root", indir ) ;
///   n_added = ch->Add(fpat) ;
///   printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

///   fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
///   fhl -> Loop() ;


   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "znunu" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_ZJetsToNuNu*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_ZJetsToNuNu*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;

   //----------------------

      ch = new TChain("tree") ;
      sprintf( sample_name, "lostlep" ) ;

      sprintf( fpat, "%s/tree_LDP/tree_WJetsToLNu*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_LDP/tree_TTJets_DiLept_ght_lt600-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_LDP/tree_TTJets_HT-*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_LDP/tree_TTJets_SingleLept*_ght_lt600*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_LDP/tree_TTJets_hadonly_ght_lt600-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;


      sprintf( fpat, "%s/tree_signal/tree_WJetsToLNu*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_TTJets_DiLept_ght_lt600-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_TTJets_HT-*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_TTJets_SingleLept*_ght_lt600*.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;

      sprintf( fpat, "%s/tree_signal/tree_TTJets_hadonly_ght_lt600-slimskim.root", indir ) ;
      n_added = ch->Add(fpat) ;
      printf("  Added %d files matching %s to sample %s\n", n_added, fpat, sample_name ) ;


      sprintf( sample_name, "lostlep" ) ;
      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;

      sprintf( sample_name, "hadtau" ) ;
      fhl = new fill_hists_loop_v2d( ch, sample_name ) ;
      fhl -> Loop() ;

   //----------------------




   } // run_fill_hists_loop_v2d



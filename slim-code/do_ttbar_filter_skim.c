
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TClassTable.h"
#include "TMath.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

#include <iostream>

#include <vector>
using std::vector ;

  void do_ttbar_filter_skim( int sample_type=-1,
                             const char* infile_name = "/Users/owen/work/cms/ra2b-2015/tree-maker/fnal-prod-v2-skims/tree_signal/tree_TTJets.root",
                             const char* outfile_name = "/Users/owen/work/cms/ra2b-2015/tree-maker/fnal-prod-v2-skims/tree_signal/tree_TTJets_hadonly_ght_lt600.root" ) {

      if ( !TClassTable::GetDict("TLorentzVector")) {
         printf("\n\n Loading library for TLorentzVector.\n\n") ;
         gSystem -> Load( "$ROOTSYS/lib/libPhysics.so" ) ;
      }

      if ( sample_type < 1 ) {
         printf("\n\n  Give a sample type from this menu:\n") ;
         printf("    1 :  Single lepton\n") ;
         printf("    2 :  Double lepton\n") ;
         printf("    3 :  Inclusive ttbar\n") ;
         printf("    4 :  HT bins\n") ;
         printf("\n\n") ;
         return ;
      }

      TFile* infile = new TFile( infile_name ) ;
      if ( ! (infile->IsOpen()) ) return ;

      TTree* in_tree = (TTree*) infile->Get("tree") ;

      Long64_t nentries = in_tree -> GetEntries() ;

      printf("\n\n Number of entries: %llu\n\n", nentries ) ;

      in_tree -> SetBranchStatus("*",1) ; // enable all branches.



     //--- Vars needed to decide whether or not to save the event.
      vector<TLorentzVector> *GenEls(0x0);
      TBranch* b_GenEls ;
      in_tree -> SetBranchAddress( "GenEls", &GenEls, &b_GenEls ) ;

      vector<TLorentzVector> *GenMus(0x0);
      TBranch* b_GenMus ;
      in_tree -> SetBranchAddress( "GenMus", &GenMus, &b_GenMus ) ;

      vector<TLorentzVector> *GenTaus(0x0);
      TBranch* b_GenTaus ;
      in_tree -> SetBranchAddress( "GenTaus", &GenTaus, &b_GenTaus ) ;

   // Double_t GenHT ;
   // TBranch* b_GenHT ;
   // in_tree -> SetBranchAddress( "GenHT", &GenHT, &b_GenHT ) ;

      Double_t madHT ;
      TBranch* b_madHT ;
      in_tree -> SetBranchAddress( "madHT", &madHT, &b_madHT ) ;




     //--- Open output file
      printf("\n\n Output file: %s\n\n", outfile_name ) ;
      TFile* outfile = new TFile( outfile_name, "recreate" ) ;
      TTree* outReducedTree = in_tree->CloneTree(0) ;



      int nrej_genlep(0) ;
      int nrej_genht(0) ;
      int nkeep(0) ;

     //--- Loop over the events.
      TStopwatch sw ;
      sw.Start() ;
      int time(0) ;
      float projected_remaining(999999.) ;
      for ( Long64_t ievt=0; ievt<nentries; ievt++ ) {

         if ( ievt%1000 == 0 ) {
            int thistime = sw.RealTime() ;
            sw.Continue() ;
            if ( thistime < 2 ) {
               printf("   %10llu out of %10llu  (%6.1f%%) \r", ievt, nentries, 100.*ievt/(1.*nentries) ) ;
            } else {
               if ( thistime > time ) projected_remaining = (1.*thistime)/(1.*ievt)*(nentries-ievt) ;
               if ( projected_remaining < 100 ) {
                  printf("   %10llu out of %10llu  (%6.1f%%)    seconds remaining %4.0f                       \r", ievt, nentries, 100.*ievt/(1.*nentries), projected_remaining ) ;
               } else if ( projected_remaining < 3600 ) {
                  printf("   %10llu out of %10llu  (%6.1f%%)    time remaining     %2d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries),
                       TMath::Nint(projected_remaining)/60, TMath::Nint(projected_remaining)%60 ) ;
               } else {
                  printf("   %10llu out of %10llu  (%6.1f%%)    time remaining  %2d:%02d:%02d   \r", ievt, nentries, 100.*ievt/(1.*nentries),
                       TMath::Nint(projected_remaining)/3600, (TMath::Nint(projected_remaining)%3600)/60, TMath::Nint(projected_remaining)%60 ) ;
               }
            }
            cout << flush ;
            time = thistime ;
         }

         in_tree -> GetEntry(ievt) ;


            if ( sample_type < 4 ) {
               if ( madHT >= 600 ) {
                  nrej_genht++ ;
                  continue ;
               }
            } else {
               if ( madHT < 600 ) {
                  nrej_genht++ ;
                  continue ;
               }
            }
            if ( sample_type == 3 ) {
               //printf("  GenEls, GenMus, GenTaus: %lu, %lu, %lu\n", GenEls->size(), GenMus->size(), GenTaus->size() ) ;
               if ((GenEls->size() + GenMus->size() + GenTaus->size()) > 0) {
                  nrej_genlep++ ;
                  continue ;
               }
            }


         outReducedTree->Fill() ;
         nkeep++ ;

      } // ievt.

      printf("\n\n\n Done.\n\n\n") ;

      printf("\n\n Output file:  %s\n\n\n", outfile_name ) ;

      if ( time > 3600 ) {
         printf( "   Total time:  %2d:%02d:%02d \n\n\n", time/3600, (time%3600)/60, time%60 ) ;
      } else if ( time > 100 ) {
         printf( "   Total time:     %02d:%02d \n\n\n", time/60, time%60 ) ;
      } else {
         printf( "   Total time:     %d seconds \n\n\n", time ) ;
      }

      printf("\n\n") ;
      printf("  Events kept: %d\n", nkeep ) ;
      printf("  Events rejected with genHT cut: %d\n", nrej_genht ) ;
      printf("  Events rejected with gen leptons cut: %d\n", nrej_genlep ) ;
      printf("\n\n") ;

      outReducedTree->AutoSave() ;

      delete outfile ;

  }


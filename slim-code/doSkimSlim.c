
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
#pragma link C++ class std::vector<bool>+;
#endif

#include <iostream>

#include <vector>
using std::vector ;


  void doSkimSlim( const char* infile_name, const char* outdir = "", bool doSkim = true, bool doSlim = true  ) {

      printf("\n\n") ;
      printf("  doSkimSlim : input file : %s\n", infile_name ) ;
      printf("  doSkimSlim : output dir : %s\n", outdir ) ;
      if ( doSkim ) printf("  skimming\n" ) ;
      if ( doSlim ) printf("  slimming\n" ) ;

      if ( !TClassTable::GetDict("TLorentzVector")) {
         printf("\n\n Loading library for TLorentzVector.\n\n") ;
         gSystem -> Load( "$ROOTSYS/lib/libPhysics.so" ) ;
      }

      TFile* infile = new TFile( infile_name ) ;
      if ( ! (infile->IsOpen()) ) return ;

      TTree* inReducedTree = (TTree*) infile->Get("tree") ;

      Long64_t nentries = inReducedTree -> GetEntries() ;

      printf("\n\n Number of entries: %llu\n\n", nentries ) ;

      if ( doSlim ) {

         inReducedTree -> SetBranchStatus("*",0) ; // disable all branches.

        //--- Now enable only the ones we want to save in the output file.
        //    Add anything you care about here.

         inReducedTree -> SetBranchStatus("RunNum",1) ;
         inReducedTree -> SetBranchStatus("LumiBlockNum",1) ;
         inReducedTree -> SetBranchStatus("EvtNum",1) ;
         inReducedTree -> SetBranchStatus("Weight",1) ;

         inReducedTree -> SetBranchStatus("BTags",1) ;
         inReducedTree -> SetBranchStatus("CaloMET",1) ;
         inReducedTree -> SetBranchStatus("CaloMETPhi",1) ;
         inReducedTree -> SetBranchStatus("CrossSection",1) ;
         inReducedTree -> SetBranchStatus("CSCTightHaloFilter",1) ;
         inReducedTree -> SetBranchStatus("DeltaPhi1",1) ;
         inReducedTree -> SetBranchStatus("DeltaPhi2",1) ;
         inReducedTree -> SetBranchStatus("DeltaPhi3",1) ;
         inReducedTree -> SetBranchStatus("DeltaPhi4",1) ;
         inReducedTree -> SetBranchStatus("eeBadScFilter",1) ;
         inReducedTree -> SetBranchStatus("GenMET",1) ;
         inReducedTree -> SetBranchStatus("GenMETPhi",1) ;
         inReducedTree -> SetBranchStatus("GenMHT",1) ;
         inReducedTree -> SetBranchStatus("GenMHTPhi",1) ;
         inReducedTree -> SetBranchStatus("HBHEIsoNoiseFilter",1) ;
         inReducedTree -> SetBranchStatus("HBHENoiseFilter",1) ;
         inReducedTree -> SetBranchStatus("HT",1) ;
         inReducedTree -> SetBranchStatus("MET",1) ;
         inReducedTree -> SetBranchStatus("METPhi",1) ;
         inReducedTree -> SetBranchStatus("MHT",1) ;
         inReducedTree -> SetBranchStatus("MHTPhi",1) ;
         inReducedTree -> SetBranchStatus("NJets",1) ;
         inReducedTree -> SetBranchStatus("NVtx",1) ;
         inReducedTree -> SetBranchStatus("SusyLSPMass",1) ;
         inReducedTree -> SetBranchStatus("SusyMotherMass",1) ;
         inReducedTree -> SetBranchStatus("TriggerNames",1) ;
         inReducedTree -> SetBranchStatus("TriggerPass",1) ;
         inReducedTree -> SetBranchStatus("Jets",1) ;
         inReducedTree -> SetBranchStatus("Jets_muonEnergyFraction",1) ;
         inReducedTree -> SetBranchStatus("Jets_hadronFlavor",1) ;
         inReducedTree -> SetBranchStatus("Jets_partonFlavor",1) ;
         inReducedTree -> SetBranchStatus("Jets_qgLikelihood",1) ;
         inReducedTree -> SetBranchStatus("Jets_qgAxis2",1) ;
         inReducedTree -> SetBranchStatus("Jets_qgMult",1) ;
         inReducedTree -> SetBranchStatus("Jets_qgPtD",1) ;

         inReducedTree -> SetBranchStatus("GenEls",1) ;
         inReducedTree -> SetBranchStatus("GenMus",1) ;
         inReducedTree -> SetBranchStatus("GenTaus",1) ;
         inReducedTree -> SetBranchStatus("GenTau_GenTauHad",1) ;
         inReducedTree -> SetBranchStatus("GenHT",1) ;
         inReducedTree -> SetBranchStatus("madHT",1) ;
         
         inReducedTree -> SetBranchStatus("JetID",1) ;
         inReducedTree -> SetBranchStatus("EcalDeadCellTriggerPrimitiveFilter",1) ;

        //-- new for fnal-prod-v9
         inReducedTree -> SetBranchStatus("globalTightHalo2016Filter",1) ;
         inReducedTree -> SetBranchStatus("BadChargedCandidateFilter",1) ;
         inReducedTree -> SetBranchStatus("BadPFMuonFilter",1) ;
         inReducedTree -> SetBranchStatus("PFCaloMETRatio",1) ;
         inReducedTree -> SetBranchStatus("noMuonJet",1) ;

         inReducedTree -> SetBranchStatus("puWeight",1) ;
         
      } else {

         inReducedTree -> SetBranchStatus("*",1) ; // enable all branches.

      }








     //--- Vars needed to decide whether or not to save the event.
     //    Also vars needed to construct new vars to be saved.


      Int_t NJets ;
      Double_t  HT, MHT ;
      inReducedTree -> SetBranchAddress( "NJets", &NJets ) ;
      inReducedTree -> SetBranchAddress( "HT", &HT ) ;
      inReducedTree -> SetBranchAddress( "MHT", &MHT ) ;

      Double_t DeltaPhi1, DeltaPhi2, DeltaPhi3, DeltaPhi4 ;
      inReducedTree -> SetBranchAddress( "DeltaPhi1", &DeltaPhi1 ) ;
      inReducedTree -> SetBranchAddress( "DeltaPhi2", &DeltaPhi2 ) ;
      inReducedTree -> SetBranchAddress( "DeltaPhi3", &DeltaPhi3 ) ;
      inReducedTree -> SetBranchAddress( "DeltaPhi4", &DeltaPhi4 ) ;

      Double_t METPhi ;
      inReducedTree -> SetBranchAddress( "METPhi", &METPhi ) ;

      std::vector<string>  *TriggerNames = 0 ;
      inReducedTree -> SetBranchAddress( "TriggerNames", &TriggerNames ) ;

      std::vector<int>     *TriggerPass = 0 ;
      inReducedTree -> SetBranchAddress( "TriggerPass", &TriggerPass ) ;




      vector<TLorentzVector> *Jets(0x0);
      TBranch* b_Jets ;
      inReducedTree -> SetBranchAddress( "Jets", &Jets, &b_Jets ) ;

      vector<double>  *Jets_muonEnergyFraction(0x0);
      TBranch* b_Jets_muonEnergyFraction ;
      inReducedTree -> SetBranchAddress( "Jets_muonEnergyFraction", &Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction ) ;

      vector<bool>    *GenTau_GenTauHad ;
      TBranch* b_GenTau_GenTauHad ;
      GenTau_GenTauHad = 0 ;
      inReducedTree -> SetBranchAddress( "GenTau_GenTauHad", &GenTau_GenTauHad, &b_GenTau_GenTauHad ) ;

      vector<TLorentzVector> *GenTaus(0x0);
      TBranch* b_GenTaus ;
      inReducedTree -> SetBranchAddress( "GenTaus", &GenTaus, &b_GenTaus ) ;







     //--- Open output file
      TString outfile_name ;
      if ( strlen(outdir) > 0 ) {
         char command[10000] ;
         printf("\n Making output directory : %s\n", outdir ) ;
         sprintf( command, "mkdir -p %s", outdir ) ;
         gSystem -> Exec( command ) ;
         sprintf( command, "basename %s .root", infile_name ) ;
         TString filename_only = gSystem -> GetFromPipe( command ) ;
         if ( doSlim ) {
            outfile_name = TString( outdir ) + TString("/") + filename_only + TString("-slimskim.root") ;
         } else {
            outfile_name = TString( outdir ) + TString("/") + filename_only + TString("-skim.root") ;
         }
      } else {
         outfile_name = TString( infile_name ) ;
         if ( doSlim ) {
            if ( outfile_name.Contains("-skim.root") ) {
               outfile_name.ReplaceAll( "-skim.root", "-slimskim.root" ) ;
            } else {
               if ( doSkim ) {
                  outfile_name.ReplaceAll( ".root", "-slimskim.root" ) ;
               } else {
                  outfile_name.ReplaceAll( ".root", "-slim.root" ) ;
               }
            }
         } else {
            outfile_name.ReplaceAll( ".root", "-skim.root" ) ;
         }
         if ( outfile_name.CompareTo( infile_name ) == 0 ) {
            printf("\n\n *** Input and output file names same.  Input doesn't contain .root in name?\n") ;
            printf("    input: %s \n", infile_name ) ;
            printf("    output: %s \n", outfile_name.Data() ) ;
            return ;
         }
      }
      printf("\n\n Output file: %s\n\n", outfile_name.Data() ) ;
      char command[10000] ;
      sprintf( command, "ls %s >& /dev/null", outfile_name.Data() ) ;
      int returnstat = gSystem->Exec( command ) ;
      if ( returnstat == 0 ) {
         char mvfile[10000] ;
         sprintf( mvfile, "%s-old", outfile_name.Data() ) ;
         printf("\n\n *** Output file already exists.  Moving it to %s\n\n", mvfile ) ;
         sprintf( command, "mv %s %s", outfile_name.Data(), mvfile ) ;
         gSystem->Exec( command ) ;
      }
      TFile* outfile = new TFile( outfile_name, "recreate" ) ;
      TTree* outReducedTree = inReducedTree->CloneTree(0) ;




      //--- New variables to be added to output tree.

      bool passBadMuFilter ;
      outReducedTree -> Branch( "passBadMuFilter", &passBadMuFilter, "passBadMuFilter/B" ) ;

      bool noMuonJet ; // Kevin's version: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction#Run2ProductionV9
      outReducedTree -> Branch( "noMuonJet", &noMuonJet, "noMuonJet/B" ) ;

      double minDeltaPhi ;
      outReducedTree -> Branch( "minDeltaPhi", &minDeltaPhi, "minDeltaPhi/D" ) ;

      bool inLDP ;
      outReducedTree -> Branch( "inLDP", &inLDP, "inLDP/B" ) ;

      bool hasHadTau ;
      outReducedTree -> Branch( "hasHadTau", &hasHadTau, "hasHadTau/B" ) ;



      bool pass_pfmet100_trig ;
      outReducedTree -> Branch( "pass_pfmet100_trig", &pass_pfmet100_trig, "pass_pfmet100_trig/B" ) ;

      bool pass_pfmetnomu100_trig ;
      outReducedTree -> Branch( "pass_pfmetnomu100_trig", &pass_pfmetnomu100_trig, "pass_pfmetnomu100_trig/B" ) ;

      bool pass_ht300_met100_trig ;
      outReducedTree -> Branch( "pass_ht300_met100_trig", &pass_ht300_met100_trig, "pass_ht300_met100_trig/B" ) ;

      bool pass_ht800_trig ;
      outReducedTree -> Branch( "pass_ht800_trig", &pass_ht800_trig, "pass_ht800_trig/B" ) ;

      bool pass_ht900_trig ;
      outReducedTree -> Branch( "pass_ht900_trig", &pass_ht900_trig, "pass_ht900_trig/B" ) ;









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

         GenTau_GenTauHad = 0 ;
         ///inReducedTree -> LoadTree(ievt) ;
         inReducedTree -> GetEntry(ievt) ;

         if ( doSkim ) {

            if ( NJets < 2 ) continue ;  //** changed from 2 to 3
            if ( MHT < 200. ) continue ; //*** changed from 250 to 200.
            if ( HT < 300. ) continue ;

         }



         
         minDeltaPhi = 99. ;
         if ( DeltaPhi1 < minDeltaPhi ) minDeltaPhi = DeltaPhi1 ;
         if ( DeltaPhi2 < minDeltaPhi ) minDeltaPhi = DeltaPhi2 ;
         if ( DeltaPhi3 < minDeltaPhi ) minDeltaPhi = DeltaPhi3 ;
         if ( DeltaPhi4 < minDeltaPhi ) minDeltaPhi = DeltaPhi4 ;







         passBadMuFilter = true ;
         for ( unsigned long ji=0; ji<Jets->size(); ji++ ) {
            if ( Jets->at(ji).Pt() < 200 ) continue ;
            if ( Jets_muonEnergyFraction->at(ji) < 0.5 ) continue ;
            double dPhi = Jets->at(ji).Phi() - METPhi ;
            if ( dPhi >  3.1415926 ) dPhi = dPhi - 2*3.14159 ;
            if ( dPhi < -3.1415926 ) dPhi = dPhi + 2*3.14159 ;
            if ( fabs( dPhi ) > 3.1415926 - 0.40 ) {
               passBadMuFilter = false ;
               break ;
            }
         } // ji

        //--- Kevin's version : https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction#Run2ProductionV9
         noMuonJet = true;
         for(unsigned j = 0; j < Jets->size(); ++j){
            double dPhi = Jets->at(j).Phi() - METPhi ;
            if ( dPhi >  3.1415926 ) dPhi = dPhi - 2*3.14159 ;
            if ( dPhi < -3.1415926 ) dPhi = dPhi + 2*3.14159 ;
            //if(Jets->at(j).Pt() > 200 && Jets_muonEnergyFraction->at(j) > 0.5 && DeltaPhi(Jets->at(j).Phi(),METPhi) > (TMath::Pi() - 0.4))
            if(Jets->at(j).Pt() > 200 && Jets_muonEnergyFraction->at(j) > 0.5 && dPhi > (TMath::Pi() - 0.4)){
               noMuonJet = false;
               break;
            }
         } 






         inLDP = true ;
         if ( DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3 ) inLDP = false ;

         hasHadTau = false ;
         if ( GenTau_GenTauHad != 0x0 ) {
            for ( unsigned long ti=0; ti<GenTau_GenTauHad->size(); ti++ ) {
               //printf( " gentau loop: %lu  %s\n", ti, GenTau_GenTauHad->at(ti)?"T":"F" ) ; fflush(stdout) ;
               if ( GenTau_GenTauHad->at(ti) ) { hasHadTau = true ; }
            } // ti
         }

         pass_pfmet100_trig = false ;
         pass_pfmetnomu100_trig = false ;
         pass_ht300_met100_trig = false ;
         pass_ht800_trig = false ;
         pass_ht900_trig = false ;
         for ( unsigned long ti=0; ti<TriggerPass->size(); ti++ ) {
            TString tname( TriggerNames->at(ti).c_str() ) ;
            int tstatus = TriggerPass->at(ti) ;
            if ( tname.Contains( "HLT_PFMET100_PFMHT100_IDTight_v" ) && tstatus > 0 ) pass_pfmet100_trig = true ;
            if ( tname.Contains( "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v" ) && tstatus > 0 ) pass_pfmetnomu100_trig = true ;
            if ( tname.Contains( "HLT_PFHT300_PFMET100_v" ) && tstatus > 0 ) pass_ht300_met100_trig = true ;
            if ( tname.Contains( "HLT_PFHT800_v" ) && tstatus > 0 ) pass_ht800_trig = true ;
            if ( tname.Contains( "HLT_PFHT900_v" ) && tstatus > 0 ) pass_ht900_trig = true ;
            //////printf("  %3lu : %2d %s\n", ti, tstatus, tname.Data() ) ;
         } // ti
         


         outReducedTree->Fill() ;

      } // ievt.

      printf("\n\n\n Done.\n\n\n") ;

      printf("\n\n Output file:  %s\n\n\n", outfile_name.Data() ) ;

      if ( time > 3600 ) {
         printf( "   Total time:  %2d:%02d:%02d \n\n\n", time/3600, (time%3600)/60, time%60 ) ;
      } else if ( time > 100 ) {
         printf( "   Total time:     %02d:%02d \n\n\n", time/60, time%60 ) ;
      } else {
         printf( "   Total time:     %d seconds \n\n\n", time ) ;
      }

      outReducedTree->AutoSave() ;

      delete outfile ;

  }


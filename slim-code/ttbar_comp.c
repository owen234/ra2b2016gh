
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDirectory.h"


   void ttbar_comp( const char* tree_dir = "/data/strange2/owen/fnal-prod-v9-skims/tree_signal/slim" ) {

      gDirectory -> Delete( "h*" ) ;

      TChain* ch_incl = new TChain( "tree" ) ;
      TChain* ch_ht = new TChain( "tree" ) ;
      TChain* ch_sl = new TChain( "tree" ) ;
      TChain* ch_dl = new TChain( "tree" ) ;

      TChain* ch_combo = new TChain( "tree" ) ;

      char fpat[10000] ;

      int nadded(0) ;

      //-----

      sprintf( fpat, "%s/tree_TTJets-slimskim.root", tree_dir ) ;
      nadded = ch_incl -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_HT*-slimskim.root", tree_dir ) ;
      nadded = ch_ht -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_SingleLeptFromT-slimskim.root", tree_dir ) ;
      nadded = ch_sl -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_SingleLeptFromTbar-slimskim.root", tree_dir ) ;
      nadded = ch_sl -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_DiLept-slimskim.root", tree_dir ) ;
      nadded = ch_dl -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      //-----

      sprintf( fpat, "%s/tree_TTJets_hadonly_ght_lt600-slimskim.root", tree_dir ) ;
      nadded = ch_combo -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_HT*-slimskim.root", tree_dir ) ;
      nadded = ch_combo -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_SingleLeptFromT_ght_lt600-slimskim.root", tree_dir ) ;
      nadded = ch_combo -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_SingleLeptFromTbar_ght_lt600-slimskim.root", tree_dir ) ;
      nadded = ch_combo -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      sprintf( fpat, "%s/tree_TTJets_DiLept_ght_lt600-slimskim.root", tree_dir ) ;
      nadded = ch_combo -> Add( fpat ) ;
      printf(" added %d with fpat %s\n", nadded, fpat ) ;

      //-----

////  sprintf( fpat, "%s/tree_TTJets_combo-slimskim.root", tree_dir ) ;
////  nadded = ch_combo -> Add( fpat ) ;
////  printf(" added %d with fpat %s\n", nadded, fpat ) ;

      //-----

      int bins = 80 ;
      double xl = 0 ;
      double xh = 4000 ;

      TH1F* h_incl = new TH1F( "h_incl", "incl", bins, xl, xh ) ;
      TH1F* h_incl_ght_gt600 = new TH1F( "h_incl_ght_gt600", "incl, gen HT>600", bins, xl, xh ) ;
      TH1F* h_incl_ght_lt600 = new TH1F( "h_incl_ght_lt600", "incl, gen HT<600", bins, xl, xh ) ;
      TH1F* h_incl_nolep = new TH1F( "h_incl_nolep", "incl, no leptons", bins, xl, xh ) ;
      TH1F* h_incl_nolep_ght_gt600 = new TH1F( "h_incl_nolep_ght_gt600", "incl, no leptons, gen HT>600", bins, xl, xh ) ;
      TH1F* h_incl_nolep_ght_lt600 = new TH1F( "h_incl_nolep_ght_lt600", "incl, no leptons, gen HT<600", bins, xl, xh ) ;

      TH1F* h_ht = new TH1F( "h_ht", "HT bins", bins, xl, xh ) ;

      TH1F* h_sl = new TH1F( "h_sl", "SL", bins, xl, xh ) ;
      TH1F* h_sl_ght_gt600 = new TH1F( "h_sl_ght_gt600", "SL, gen HT>600", bins, xl, xh ) ;
      TH1F* h_sl_ght_lt600 = new TH1F( "h_sl_ght_lt600", "SL, gen HT<600", bins, xl, xh ) ;

      TH1F* h_dl = new TH1F( "h_dl", "DL", bins, xl, xh ) ;
      TH1F* h_dl_ght_gt600 = new TH1F( "h_dl_ght_gt600", "DL, gen HT>600", bins, xl, xh ) ;
      TH1F* h_dl_ght_lt600 = new TH1F( "h_dl_ght_lt600", "DL, gen HT<600", bins, xl, xh ) ;

      TH1F* h_combo = new TH1F( "h_combo", "Combo", bins, xl, xh ) ;

      h_combo -> Sumw2() ;

      h_incl -> Sumw2() ;
      h_incl_ght_gt600 -> Sumw2() ;
      h_incl_ght_lt600 -> Sumw2() ;
      h_incl_nolep -> Sumw2() ;
      h_incl_nolep_ght_gt600 -> Sumw2() ;
      h_incl_nolep_ght_lt600 -> Sumw2() ;

      h_ht -> Sumw2() ;

      h_sl -> Sumw2() ;
      h_sl_ght_gt600 -> Sumw2() ;
      h_sl_ght_lt600 -> Sumw2() ;

      h_dl -> Sumw2() ;
      h_dl_ght_gt600 -> Sumw2() ;
      h_dl_ght_lt600 -> Sumw2() ;


      char cuts[10000] ;

      TCanvas* can = new TCanvas( "can_ttbar_comp", "TTbar comp", 1100, 900 ) ;


      //sprintf( cuts, "( (@GenEls->size()==0 && @GenMus->size()==0 && @GenTaus->size()==0) && (@Muons.size()+@Electrons.size()) == 0 && (isoElectronTracks+isoMuonTracks+isoPionTracks)==0 && (DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3) && JetID>0)*Weight" ) ;
      sprintf( cuts, "( (@GenEls->size()==0 && @GenMus->size()==0 && @GenTaus->size()==0) )*Weight" ) ;
      ch_incl -> Draw( "HT>>h_incl_nolep", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT>600 && (@GenEls->size()==0 && @GenMus->size()==0 && @GenTaus->size()==0) )*Weight" ) ;
      ch_incl -> Draw( "HT>>h_incl_nolep_ght_gt600", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT<600 && (@GenEls->size()==0 && @GenMus->size()==0 && @GenTaus->size()==0) )*Weight" ) ;
      ch_incl -> Draw( "HT>>h_incl_nolep_ght_lt600", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT>600 )*Weight" ) ;
      ch_incl -> Draw( "HT>>h_incl_ght_gt600", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT<600 )*Weight" ) ;
      ch_incl -> Draw( "HT>>h_incl_ght_lt600", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "Weight" ) ;
      ch_incl -> Draw( "HT>>h_incl", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "Weight" ) ;
      ch_combo -> Draw( "HT>>h_combo", cuts ) ;
      can -> Update() ; can -> Draw() ;




      sprintf( cuts, "Weight" ) ;
      ch_ht -> Draw( "HT>>h_ht", cuts ) ;
      can -> Update() ; can -> Draw() ;



      sprintf( cuts, "Weight" ) ;
      ch_sl -> Draw( "HT>>h_sl", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT>600)*Weight" ) ;
      ch_sl -> Draw( "HT>>h_sl_ght_gt600", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT<600)*Weight" ) ;
      ch_sl -> Draw( "HT>>h_sl_ght_lt600", cuts ) ;
      can -> Update() ; can -> Draw() ;



      sprintf( cuts, "Weight" ) ;
      ch_dl -> Draw( "HT>>h_dl", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT>600)*Weight" ) ;
      ch_dl -> Draw( "HT>>h_dl_ght_gt600", cuts ) ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "( madHT<600)*Weight" ) ;
      ch_dl -> Draw( "HT>>h_dl_ght_lt600", cuts ) ;
      can -> Update() ; can -> Draw() ;


      TH1F* h_had_sl_dl = (TH1F*) h_incl_nolep -> Clone( "h_had_sl_dl" ) ;
      h_had_sl_dl -> Add( h_sl ) ;
      h_had_sl_dl -> Add( h_dl ) ;

      TH1F* h_had_sl_dl_ght_lt600 = (TH1F*) h_incl_nolep_ght_lt600 -> Clone( "h_had_sl_dl_ght_lt600" ) ;
      h_had_sl_dl_ght_lt600 -> Add( h_sl_ght_lt600 ) ;
      h_had_sl_dl_ght_lt600 -> Add( h_dl_ght_lt600 ) ;

      TH1F* h_stitch = (TH1F*) h_incl_nolep_ght_lt600 -> Clone( "h_stitch" ) ;
      h_stitch -> Add( h_sl_ght_lt600 ) ;
      h_stitch -> Add( h_dl_ght_lt600 ) ;
      h_stitch -> Add( h_ht ) ;


      h_stitch -> SetLineColor(4) ;
      h_incl_nolep_ght_lt600 -> SetLineColor( 6 ) ;
      h_sl_ght_lt600 -> SetLineColor(2) ;
      h_dl_ght_lt600 -> SetLineColor(3) ;
      h_ht -> SetLineColor( 7 ) ;

      h_incl -> Draw() ;
      h_stitch -> Draw( "same" ) ;

      h_sl_ght_lt600 -> Draw("same") ;
      h_dl_ght_lt600 -> Draw("same") ;
      h_incl_nolep_ght_lt600 -> Draw("same") ;
      h_ht -> Draw("same") ;

      gPad -> SetLogy(1) ;


   } // ttbar_comp






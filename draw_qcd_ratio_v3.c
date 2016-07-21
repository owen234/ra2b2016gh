
#include "TDirectory.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"


#include "histio.c"

   TH1F* get_hist( const char* hname ) ;
   void  draw_boundaries() ;
   void  draw_boundaries_nj( int nhtb ) ;

   const int nb_nj(4) ;
   const int nb_nb(4) ;
   const int nb_mht(5) ;
   const int nb_htmht(13) ;

   //---------

   void draw_qcd_ratio_v3( const char* infile = "outputfiles/hists-v2d-qcd.root", const char* outputdir = "outputfiles/" ) {

      TLine* line0 = new TLine() ;
      line0 -> SetLineColor(4) ;
      TString tstring ;
      TLine* rline = new TLine() ;
      rline -> SetLineColor(2) ;

      double ratio_max = 1.0 ;

      char command[10000] ;
      sprintf( command, "mkdir -p %s", outputdir ) ;
      gSystem -> Exec( command ) ;


      gStyle -> SetOptStat(0) ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( infile ) ;

      char hname[1000] ;
      char label[100] ;

      TH1F* h_ldp = get_hist( "h_ldp" ) ;
      TH1F* h_hdp = get_hist( "h_hdp" ) ;
      TH1F* h_max_ldp_weight = get_hist( "h_max_ldp_weight" ) ;

      TH1F* h_ratio = new TH1F( "h_ratio", "QCD H/L ratio", 160, 0.5, 160.5 ) ;
      TH1F* h_ratio_nb0 = new TH1F( "h_ratio_nb0", "QCD H/L ratio, Nb0", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nb1 = new TH1F( "h_ratio_nb1", "QCD H/L ratio, Nb1", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nb2 = new TH1F( "h_ratio_nb2", "QCD H/L ratio, Nb2", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nb3 = new TH1F( "h_ratio_nb3", "QCD H/L ratio, Nb3", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nj1 = new TH1F( "h_ratio_nj1", "QCD H/L ratio, Nj1", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nj2 = new TH1F( "h_ratio_nj2", "QCD H/L ratio, Nj2", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nj3 = new TH1F( "h_ratio_nj3", "QCD H/L ratio, Nj3", 40, 0.5, 40.5 ) ;
      TH1F* h_ratio_nj4 = new TH1F( "h_ratio_nj4", "QCD H/L ratio, Nj4", 40, 0.5, 40.5 ) ;
      TH1F* h_max_ldp_weight_160bins = new TH1F( "h_max_ldp_weight_160bins", "max LDP weight", 160, 0.5, 160.5 ) ;
      TH1F* h_ldp_160bins = new TH1F( "h_ldp_160bins", "LDP counts, 160 bins", 160, 0.5, 160.5 ) ;
      TH1F* h_hdp_160bins = new TH1F( "h_hdp_160bins", "HDP counts, 160 bins", 160, 0.5, 160.5 ) ;

      int bi_hist(0) ;
      int bi_search_hist(0) ;
      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht; bi_htmht++ ) {
               bi_hist++ ;
               sprintf( label, "%s", h_ldp -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;
               //printf( " %3d : Nj%d Nb%d HTMHT%d : %30s\n", bi_hist, bi_nj, bi_nb-1, bi_htmht, label ) ;
               if ( bi_htmht > 3 ) {
                  bi_search_hist ++ ;
                  double ldp_val = h_ldp -> GetBinContent( bi_hist ) ;
                  double ldp_err = h_ldp -> GetBinError( bi_hist ) ;
                  double hdp_val = h_hdp -> GetBinContent( bi_hist ) ;
                  double hdp_err = h_hdp -> GetBinError( bi_hist ) ;
                  double ratio_val(0.) ;
                  double ratio_err(0.) ;
                  if ( ldp_val > 0 && hdp_val > 0 ) {
                     ratio_val = hdp_val / ldp_val ;
                     ratio_err = ratio_val * sqrt( pow( ldp_err/ldp_val, 2. ) + pow( hdp_err/hdp_val, 2. ) ) ;
                  }
                  h_ratio -> SetBinContent( bi_search_hist, ratio_val ) ;
                  h_ratio -> SetBinError( bi_search_hist, ratio_err ) ;
                  h_ratio -> GetXaxis() -> SetBinLabel( bi_search_hist, label ) ;
                  h_max_ldp_weight_160bins -> SetBinContent( bi_search_hist, h_max_ldp_weight->GetBinContent( bi_hist ) ) ;
                  h_max_ldp_weight_160bins -> GetXaxis() -> SetBinLabel( bi_search_hist, label ) ;
                  h_ldp_160bins -> SetBinContent( bi_search_hist, ldp_val ) ;
                  h_ldp_160bins -> SetBinError( bi_search_hist, ldp_err ) ;
                  h_ldp_160bins -> GetXaxis() -> SetBinLabel( bi_search_hist, label ) ;
                  h_hdp_160bins -> SetBinContent( bi_search_hist, hdp_val ) ;
                  h_hdp_160bins -> SetBinError( bi_search_hist, hdp_err ) ;
                  h_hdp_160bins -> GetXaxis() -> SetBinLabel( bi_search_hist, label ) ;
                  printf( "  search %3d : %30s : R= %6.4f +/- %6.4f\n", bi_search_hist, label, ratio_val, ratio_err ) ;
                  TH1F* hp_nb(0x0) ;
                  if ( bi_nb==1 ) hp_nb = h_ratio_nb0 ;
                  if ( bi_nb==2 ) hp_nb = h_ratio_nb1 ;
                  if ( bi_nb==3 ) hp_nb = h_ratio_nb2 ;
                  if ( bi_nb==4 ) hp_nb = h_ratio_nb3 ;
                  int bi_nb_hist = (bi_nj-1)*10 + bi_htmht-3 ;
                  hp_nb -> SetBinContent( bi_nb_hist, ratio_val ) ;
                  hp_nb -> SetBinError( bi_nb_hist, ratio_err ) ;
                  hp_nb -> GetXaxis() -> SetBinLabel( bi_nb_hist, label ) ;
                  TH1F* hp_nj(0x0) ;
                  if ( bi_nj==1 ) hp_nj = h_ratio_nj1 ;
                  if ( bi_nj==2 ) hp_nj = h_ratio_nj2 ;
                  if ( bi_nj==3 ) hp_nj = h_ratio_nj3 ;
                  if ( bi_nj==4 ) hp_nj = h_ratio_nj4 ;
                  int bi_nj_hist = (bi_nb-1)*10 + bi_htmht-3 ;
                  hp_nj -> SetBinContent( bi_nj_hist, ratio_val ) ;
                  hp_nj -> SetBinError( bi_nj_hist, ratio_err ) ;
                  hp_nj -> GetXaxis() -> SetBinLabel( bi_nj_hist, label ) ;
               } // not control bin.
            } // bi_htmht
         } // bi_nb
      } // bi_nj

      h_ratio -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio -> SetMarkerStyle(20) ;

      h_max_ldp_weight_160bins -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ldp_160bins -> GetXaxis() -> LabelsOption( "v" ) ;
      h_hdp_160bins -> GetXaxis() -> LabelsOption( "v" ) ;



      h_ratio_nb0 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nb0 -> SetMarkerStyle(20) ;
      h_ratio_nb1 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nb1 -> SetMarkerStyle(20) ;
      h_ratio_nb2 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nb2 -> SetMarkerStyle(20) ;
      h_ratio_nb3 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nb3 -> SetMarkerStyle(20) ;

      h_ratio_nj1 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj1 -> SetMarkerStyle(20) ;
      h_ratio_nj2 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj2 -> SetMarkerStyle(20) ;
      h_ratio_nj3 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj3 -> SetMarkerStyle(20) ;
      h_ratio_nj4 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_nj4 -> SetMarkerStyle(20) ;

      h_ratio -> Draw() ;
      gPad -> SetGridy(1) ;

      saveHist("outputfiles/qcdmc-ratio-v3.root","h*") ;

   } // draw_qcd_ratio_v3

//===============================================================================

   TH1F* get_hist( const char* hname ) {
      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** Missing histogram : %s\n\n", hname ) ;
         gDirectory -> ls() ;
         gSystem -> Exit( -1 ) ;
      }
      return hp ;
   } // get_hist

//===============================================================================

   void draw_boundaries() {

      TLine* line1 = new TLine() ;
      line1 -> SetLineColor(2) ;
      line1 -> SetLineStyle(2) ;

      TLine* line2 = new TLine() ;
      line2 -> SetLineColor(4) ;
      line2 -> SetLineStyle(3) ;

      double ymin = gPad -> GetY1() ;
      double ymax = gPad -> GetUymax() ;
      printf("  draw_boundaries : %s, ymin = %.1f, ymax = %.1f\n", gPad->GetTitle(), ymin, ymax ) ;

      for ( int nji=1; nji<=nb_nj; nji++ ) {
         float x =  0.5+nji*nb_htmht ;
         line1 -> DrawLine( x, ymin, x, ymax ) ;
         x = (nji-1)*nb_htmht + 3.5 ;
         line2 -> DrawLine( x, ymin, x, ymax ) ;
         x = (nji-1)*nb_htmht + 6.5 ;
         line2 -> DrawLine( x, ymin, x, ymax ) ;
         x = (nji-1)*nb_htmht + 9.5 ;
         line2 -> DrawLine( x, ymin, x, ymax ) ;
         x = (nji-1)*nb_htmht + 11.5 ;
         line2 -> DrawLine( x, ymin, x, ymax ) ;
      } // nji

   } // draw_boundaries

//===============================================================================


   void draw_boundaries_nj( int nhtb ) {

      TLine* line1 = new TLine() ;
      line1 -> SetLineColor(2) ;
      line1 -> SetLineStyle(2) ;

      double ymin = gPad -> GetY1() ;
      double ymax = gPad -> GetUymax() ;
      printf("  draw_boundaries_nj : %s, ymin = %.1f, ymax = %.1f\n", gPad->GetTitle(), ymin, ymax ) ;

      for ( int nji=1; nji<=nb_nj; nji++ ) {
         float x =  0.5+nji*nhtb ;
         line1 -> DrawLine( x, ymin, x, ymax ) ;
      } // nji

   } // draw_boundaries_nj

//===============================================================================





#ifndef draw_closure_sums1_c
#define draw_closure_sums1_c

#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "histio.c"
#include "get_hist.h"

TH1F* get_hist( const char* hname ) ;

void draw_closure_sums1( const char* plotname = "njet", const char* infile = "outputfiles/closure-sums3.root" ) {

      char hname[100] ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadLeftMargin(0.13) ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( infile ) ;

      gDirectory -> ls() ;

      sprintf( hname, "h_closure_%s_model_fit", plotname ) ;
      TH1F* h_model_fit = get_hist( hname ) ;
      TString tsn( hname ) ;
      tsn.ReplaceAll( "fit", "zero_err" ) ;
      TH1F* h_model_zero_err = (TH1F*) h_model_fit -> Clone( tsn.Data() ) ;
      for ( int bi=1; bi<= h_model_zero_err->GetNbinsX(); bi++ ) {
         h_model_zero_err -> SetBinError( bi, 0.0000001 ) ;
      }

      sprintf( hname, "h_closure_%s_model_fit_and_syst", plotname ) ;
      TH1F* h_model_fit_and_syst = get_hist( hname ) ;

      sprintf( hname, "h_closure_%s_model_total", plotname ) ;
      TH1F* h_model_total = get_hist( hname ) ;

      sprintf( hname, "h_closure_%s_qcdmc", plotname ) ;
      TH1F* h_model_qcdmc = get_hist( hname ) ;

      h_model_qcdmc -> SetMarkerStyle(20) ;
      h_model_qcdmc -> SetMarkerSize(1.5) ;
      h_model_qcdmc -> SetLineWidth(2) ;

      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_dcs" ) ;
      if ( can == 0x0 ) {
         can = new TCanvas( "can_dcs", "", 700, 800 ) ;
      }
      can -> Clear() ;
      if ( strcmp( plotname, "10boxes" ) == 0 ) {
         gPad -> SetBottomMargin( 0.20 ) ;
         h_model_total -> GetXaxis() -> LabelsOption("v") ;
      } else {
         gPad -> SetBottomMargin( 0.10 ) ;
      }

      TLegend* legend(0x0) ;
      if ( strcmp( plotname, "ht" ) == 0 ) {
         legend = new TLegend( 0.60, 0.15, 0.98, 0.40 ) ;
      } else {
         legend = new TLegend( 0.60, 0.70, 0.98, 0.95 ) ;
      }
      legend -> AddEntry( h_model_fit, "Model, Fit" ) ;
      legend -> AddEntry( h_model_fit_and_syst, "Model, Fit + MHT" ) ;
      legend -> AddEntry( h_model_total, "Model, Fit + MHT + MC" ) ;
      legend -> AddEntry( h_model_qcdmc, "QCD Monte Carlo" ) ;


      h_model_fit -> SetFillColor( kRed-7 ) ;
      h_model_fit_and_syst -> SetFillColor( kRed-9 ) ;
      h_model_total -> SetFillColor( kRed-10 ) ;
      h_model_zero_err -> SetLineColor( 4 ) ;

      h_model_total -> SetYTitle( "QCD BG events" ) ;
      h_model_total -> SetTitleSize( 0.05, "y" ) ;
      if ( strcmp( plotname, "njet" ) == 0 ) { h_model_total -> SetTitle( "Njet bins" ) ; }
      if ( strcmp( plotname, "nb" ) == 0 ) { h_model_total -> SetTitle( "Nb bins" ) ; }
      if ( strcmp( plotname, "mht" ) == 0 ) { h_model_total -> SetTitle( "MHT bins" ) ; }
      if ( strcmp( plotname, "ht" ) == 0 ) { h_model_total -> SetTitle( "HT bins" ) ; }
      if ( strcmp( plotname, "10boxes" ) == 0 ) { h_model_total -> SetTitle( "MHT-HT bins" ) ; }

      h_model_total -> Draw( "e2" ) ;
      h_model_fit_and_syst -> Draw( "e2 same" ) ;
      h_model_fit -> Draw("e2 same" ) ;
      h_model_zero_err -> Draw( "same" ) ;
      gPad -> SetGridy(1) ;

      h_model_qcdmc -> Draw( "same" ) ;
      h_model_total -> Draw( "axig same" ) ;
      h_model_total -> Draw( "axis same" ) ;

      legend -> Draw() ;

      char fname[1000] ;
      sprintf( fname, "outputfiles/closure-sum-%s.pdf", plotname ) ;
      can -> SaveAs( fname ) ;

} // draw_closure_sums1

#endif

#include "TDirectory.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TStyle.h"

#include "get_hist.h"
#include "histio.c"

  //--------

   void draw_mcexpect1(
         const char* infile    = "outputfiles/calc-mc-ratios.root",
         const char* iodir     = "outputfiles/mcexpect-plots"
      ) {

     char command[10000] ;
     sprintf( command, "mkdir -p %s", iodir ) ;
     gSystem -> Exec( command ) ;


     char fname[10000] ;

     gStyle -> SetPadBottomMargin( 0.25 ) ;
     gStyle -> SetOptStat(0) ;

     gDirectory -> Delete( "h*" ) ;

     loadHist( infile ) ;


     printf("\n") ;
     gDirectory -> ls( ) ;
     printf("\n") ;



   //-----

     TH1F* h_ldp_lostlep = get_hist( "h_ldp_expected_error_htgrouping_lostlep" ) ;
     TH1F* h_ldp_hadtau  = get_hist( "h_ldp_expected_error_htgrouping_hadtau" ) ;
     TH1F* h_ldp_znunu   = get_hist( "h_ldp_expected_error_htgrouping_znunu" ) ;
     TH1F* h_ldp_qcd     = get_hist( "h_ldp_expected_error_htgrouping_qcd" ) ;

     TH1F* h_ldp_mcsum = (TH1F*) h_ldp_lostlep -> Clone( "h_ldp_mcsum" ) ;
     h_ldp_mcsum -> Add( h_ldp_hadtau ) ;
     h_ldp_mcsum -> Add( h_ldp_znunu ) ;
     h_ldp_mcsum -> Add( h_ldp_qcd ) ;

     h_ldp_lostlep -> SetFillColor( kBlue-10 ) ;
     h_ldp_hadtau  -> SetFillColor( kCyan-10 ) ;
     h_ldp_znunu   -> SetFillColor( kGreen-7 ) ;
     h_ldp_qcd     -> SetFillColor( kRed-9 ) ;


     THStack* hstack_ldp = new THStack( "hstack_ldp", "hstack_ldp" ) ;
     hstack_ldp -> Add( h_ldp_qcd ) ;
     hstack_ldp -> Add( h_ldp_lostlep ) ;
     hstack_ldp -> Add( h_ldp_hadtau ) ;
     hstack_ldp -> Add( h_ldp_znunu ) ;


   //-----

     TH1F* h_hdp_lostlep = get_hist( "h_hdp_expected_error_htgrouping_lostlep" ) ;
     TH1F* h_hdp_hadtau  = get_hist( "h_hdp_expected_error_htgrouping_hadtau" ) ;
     TH1F* h_hdp_znunu   = get_hist( "h_hdp_expected_error_htgrouping_znunu" ) ;
     TH1F* h_hdp_qcd     = get_hist( "h_hdp_expected_error_htgrouping_qcd" ) ;

     TH1F* h_hdp_mcsum = (TH1F*) h_hdp_lostlep -> Clone( "h_hdp_mcsum" ) ;
     h_hdp_mcsum -> Add( h_hdp_hadtau ) ;
     h_hdp_mcsum -> Add( h_hdp_znunu ) ;
     h_hdp_mcsum -> Add( h_hdp_qcd ) ;

     h_hdp_lostlep -> SetFillColor( kBlue-10 ) ;
     h_hdp_hadtau  -> SetFillColor( kCyan-10 ) ;
     h_hdp_znunu   -> SetFillColor( kGreen-7 ) ;
     h_hdp_qcd     -> SetFillColor( kRed-9 ) ;

     THStack* hstack_hdp = new THStack( "hstack_hdp", "hstack_hdp" ) ;
     hstack_hdp -> Add( h_hdp_qcd ) ;
     hstack_hdp -> Add( h_hdp_lostlep ) ;
     hstack_hdp -> Add( h_hdp_hadtau ) ;
     hstack_hdp -> Add( h_hdp_znunu ) ;



/////-----

     TCanvas* c1 = new TCanvas("c1","LDP",900,500) ;
     hstack_ldp -> Draw( "hist" ) ;
     hstack_ldp -> Draw( "same" ) ;
     h_ldp_mcsum -> Draw( "axig same" ) ;
     gPad -> SetGridy(1) ;

     c1 -> Update() ; c1 -> Draw() ;
     sprintf( fname, "%s/mcexpect-ldp-liny-full.pdf", iodir ) ;
     c1 -> SaveAs( fname ) ;

//   h_ldp_data -> SetMaximum( 1.5 * ( h_ldp_data -> GetBinContent( 14 ) ) ) ;
//   c1 -> Update() ; c1 -> Draw() ;
//   sprintf( fname, "%s/datamc-ldp-liny-zoom1.pdf", iodir ) ;
//   c1 -> SaveAs( fname ) ;

//   h_ldp_data -> SetMaximum( 2.0 * ( h_ldp_data -> GetBinContent( 27 ) ) ) ;
//   c1 -> Update() ; c1 -> Draw() ;
//   sprintf( fname, "%s/datamc-ldp-liny-zoom2.pdf", iodir ) ;
//   c1 -> SaveAs( fname ) ;





     TCanvas* c2 = new TCanvas("c2","LDP",900,500) ;

     h_hdp_mcsum -> Draw( "hist" ) ;
     hstack_hdp -> Draw( "hist same" ) ;
     hstack_hdp -> Draw( "same" ) ;
     h_hdp_mcsum -> Draw( "axig same" ) ;
     gPad -> SetGridy(1) ;

     c2 -> Update() ; c2 -> Draw() ;
     sprintf( fname, "%s/mcexpect-hdp-liny-full.pdf", iodir ) ;
     c2 -> SaveAs( fname ) ;

     h_hdp_mcsum -> SetMaximum( 1.25 * ( h_hdp_mcsum -> GetBinContent( 5 ) ) ) ;
     c2 -> Update() ; c2 -> Draw() ;
     sprintf( fname, "%s/mcexpect-hdp-liny-zoom1.pdf", iodir ) ;
     c2 -> SaveAs( fname ) ;

     h_hdp_mcsum -> SetMaximum( 2.0 * ( h_hdp_mcsum -> GetBinContent( 9 ) ) ) ;
     c2 -> Update() ; c2 -> Draw() ;
     sprintf( fname, "%s/mcexpect-hdp-liny-zoom2.pdf", iodir ) ;
     c2 -> SaveAs( fname ) ;


     h_hdp_qcd -> SetMarkerStyle(20) ;

       h_hdp_qcd -> SetMaximum(65 ) ;
       h_hdp_qcd -> SetMinimum(-15 ) ;

     //h_hdp_qcd -> SetMaximum(300 ) ;
     //h_hdp_qcd -> SetMinimum(-50 ) ;

     h_hdp_qcd -> Draw( "hist" ) ;
     h_hdp_qcd -> Draw( "same" ) ;
     h_hdp_qcd -> Draw( "axig same" ) ;
     c2 -> Update() ; c2 -> Draw() ;
     sprintf( fname, "%s/mcexpect-hdp-qcd.pdf", iodir ) ;
     c2 -> SaveAs( fname ) ;


   } // draw_mcexpect1

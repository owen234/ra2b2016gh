
#include "../histio.c"

   TH1F* get_hist( const char* hname, bool push_zero_off = true ) ;

   void draw_ratios_mc_and_model_v5( int draw_htbin = 3,
                                     double drmax = 1.5, double rmax = 0.14,
                                     const char* plottype = "nbsum",
                                     bool show_model = true,
                                     const char* mcfile = "outputfiles/qcd-ratio-nbsum.root",
                                     const char* modelfile = "outputfiles/calc-model-ratios-v5.root",
                                     const char* plotdir = "outputfiles/model-mc-comp-v5/" ) {


      gDirectory -> Delete( "h*" ) ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin( 0.30 ) ;
      gStyle -> SetPadLeftMargin(0.15) ;

      char hname[100] ;
      char fname[10000] ;
      char command[10000] ;

      sprintf( command, "mkdir -p %s", plotdir ) ;
      gSystem -> Exec( command ) ;


      printf("\n\n") ;
      loadHist( mcfile ) ;
      printf("\n") ;
      loadHist( modelfile ) ;
      printf("\n") ;

      gDirectory -> ls() ;

      printf("\n\n") ;

      char ht_bin_str[100] ;
      if ( draw_htbin == 1 ) sprintf( ht_bin_str, "htl" ) ;
      if ( draw_htbin == 2 ) sprintf( ht_bin_str, "htm" ) ;
      if ( draw_htbin == 3 ) sprintf( ht_bin_str, "hth" ) ;


      sprintf( hname, "h_ratio_%s", ht_bin_str ) ;
      TH1F* h_mc_ratio = get_hist( hname ) ;
      sprintf( hname, "h_model_ratio_%s_%s", ht_bin_str, plottype ) ;
      TH1F* h_model_ratio = get_hist( hname ) ;

      h_model_ratio -> SetLineColor(2) ;
      TH1F* h_model_ratio_copy1 = (TH1F*) h_model_ratio -> Clone( "h_model_ratio_copy1" ) ;
      h_model_ratio -> SetFillColor( kRed-10 ) ;


      //////////if ( draw_htbin == 1 ) h_mc_ratio -> GetXaxis() -> SetRange(1,6) ;
      if ( draw_htbin == 1 ) h_mc_ratio -> GetXaxis() -> SetRange(1,9) ;


    //---
      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_draw_ratios_mc_and_model_v5" ) ;
      if ( can == 0x0 ) { can = new TCanvas( "can_draw_ratios_mc_and_model_v5", "", 900, 700 ) ; }
      can -> Clear() ;


      h_mc_ratio -> SetMaximum( rmax ) ;

      h_mc_ratio -> SetYTitle( "R^{QCD}" ) ;
      h_mc_ratio -> SetTitleSize( 0.06, "y" ) ;

      h_mc_ratio -> Draw("e0") ;
      if ( show_model ) {
         h_model_ratio  -> Draw("E2 same") ;
         h_model_ratio_copy1  -> Draw("hist same") ;
         h_mc_ratio -> Draw("same") ;
      }
      gPad -> SetGridy(1) ;
         h_mc_ratio -> Draw("e0 same") ;
      h_mc_ratio -> Draw("axig same") ;

      sprintf( fname, "%s/plot-ratio-%s-%s.pdf", plotdir, ht_bin_str, plottype ) ;
      can -> SaveAs( fname ) ;


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      sprintf( hname, "h_dr_%s", ht_bin_str ) ;
      TH1F* h_mc_dr = get_hist( hname ) ;
      sprintf( hname, "h_model_dr_%s_%s", ht_bin_str, plottype ) ;
      TH1F* h_model_dr = get_hist( hname, false ) ;

      h_model_dr -> SetLineColor(2) ;
      TH1F* h_model_dr_copy1 = (TH1F*) h_model_dr -> Clone( "h_model_dr_copy1" ) ;
      h_model_dr -> SetFillColor( kRed-10 ) ;


      for ( int bi=1; bi<=h_mc_dr->GetNbinsX(); bi++ ) { if ( h_mc_dr->GetBinContent( bi ) <= 0 ) h_mc_dr->SetBinContent(bi, -9 ) ; }

      ////////////if ( draw_htbin == 1 )  h_mc_dr -> GetXaxis() -> SetRange(1,6) ;
      if ( draw_htbin == 1 )  h_mc_dr -> GetXaxis() -> SetRange(1,9) ;


    //---
      TCanvas* candr = (TCanvas*) gDirectory -> FindObject( "candr_draw_drs_mc_and_model_v5" ) ;
      if ( candr == 0x0 ) { candr = new TCanvas( "candr_draw_drs_mc_and_model_v5", "HT low", 900, 700 ) ; }
      candr -> Clear() ;

      h_mc_dr -> SetMaximum( drmax ) ;

      h_mc_dr -> SetYTitle( "R^{QCD}(MHT_{k}) / R^{QCD}(MHT_{c})" ) ;
      h_mc_dr -> SetTitleSize( 0.05, "y" ) ;


      h_mc_dr -> Draw("e0") ;
      if ( show_model ) {
         h_model_dr  -> Draw("E2 same") ;
         h_model_dr_copy1  -> Draw("hist same") ;
         h_mc_dr -> Draw("same") ;
         h_mc_dr -> Draw("e0 same") ;
      }
      gPad -> SetGridy(1) ;
      h_mc_dr -> Draw("axig same") ;
           h_mc_dr -> Draw("e0 same") ;
      candr -> Update() ; candr -> Draw() ;

      sprintf( fname, "%s/plot-dr-%s-%s.pdf", plotdir, ht_bin_str, plottype ) ;
      candr -> SaveAs( fname ) ;



    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   } // draw_ratios_mc_and_model_v5

//===============================================================================

   TH1F* get_hist( const char* hname, bool push_zero_off ) {
      printf("  get_hist:  finding histogram %s\n", hname ) ;
      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** Missing histogram : %s\n\n", hname ) ;
         gDirectory -> ls() ;
         gSystem -> Exit( -1 ) ;
      }
      if ( push_zero_off ) {
         for ( int bi=1; bi<= hp->GetNbinsX(); bi++ ) {
            if ( hp -> GetBinContent(bi) == 0 ) {
               hp -> SetBinContent(bi, -9.) ;
               hp -> SetBinError(bi, 0.) ;
            }
         }
      }
      return hp ;
   } // get_hist

//===============================================================================


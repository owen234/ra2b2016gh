#ifndef  draw_model_vs_mc_c
#define  draw_model_vs_mc_c

#include "TGraphErrors.h"
#include "TSystem.h"
#include "TPad.h"
#include "histio.c"
#include "TCanvas.h"
#include "TStyle.h"

#include "get_hist.h"
#include "binning.h"

   void draw_model_vs_mc( const char* model_file = "outputfiles/gci-output.root",
                          const char* qcdmc_file = "outputfiles/qcdmc-ratio-v3.root",
                          const char* output_dir = "outputfiles/model-vs-mc" ) {

      setup_bins();

      gStyle -> SetOptStat(0) ;

      char fname[10000] ;
      char command[10000] ;
      sprintf( command, "mkdir -p %s", output_dir ) ;
      gSystem -> Exec( command ) ;

      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPadBottomMargin(0.35) ;

      loadHist( model_file, "model" ) ;
      loadHist( qcdmc_file, "qcdmc" ) ;

      gDirectory -> ls( "*ratio*" ) ;

      TH1F* h_model = get_hist( "h_ratio_all_model" ) ;
      TH1F* h_qcdmc = get_hist( "h_ratio_qcdmc" ) ;

      TH1F* h_qcdmc_minus_model = (TH1F*) h_qcdmc -> Clone( "h_qcdmc_minus_model" ) ;

      for ( int bi=1; bi<=h_qcdmc_minus_model-> GetNbinsX(); bi++ ) {
         float model_val = h_model -> GetBinContent( bi ) ;
         float qcdmc_val = h_qcdmc -> GetBinContent( bi ) ;
         if ( qcdmc_val > 0 ) {
            h_qcdmc_minus_model -> SetBinContent( bi, qcdmc_val - model_val ) ;
         } else {
            h_qcdmc_minus_model -> SetBinContent( bi, 0. ) ;
         }
         h_qcdmc_minus_model -> SetBinError( bi, 0. ) ;
      }

      h_model -> SetMarkerStyle(0) ;

      TH1F* h_model_noerrs = (TH1F*) h_model -> Clone( "h_model_noerrs" ) ;
      for ( int bi=1; bi<=h_model_noerrs->GetNbinsX(); bi++ ) { h_model_noerrs -> SetBinError( bi, 0.0000001 ) ; }

      h_model -> SetMarkerStyle(0) ;
      h_model -> SetLineColor(4) ;
      h_model -> SetFillColor(kRed-10) ;
      h_model_noerrs -> SetLineColor(4) ;
      h_model_noerrs -> SetLineWidth(2) ;

      TGraphErrors* gr_model(0x0) ;
      {
         double x[1000], y[1000], ex[1000], ey[1000];

	 int nbins = h_model->GetNbinsX();
         for ( int bi=1; bi <= nbins; bi++ ) {
            x [bi] = bi ;
            ex[bi] = 0.5 ;
            y [bi] = h_model -> GetBinContent( bi ) ;
            ey[bi] = h_model -> GetBinError( bi ) ;
         } // bi

	 nbins = h_qcdmc->GetNbinsX();
         for ( int bi=1; bi <= nbins; bi++ ) {
            if ( h_qcdmc->GetBinContent( bi ) == 0. && h_qcdmc->GetBinError( bi) == 0. ) {
               h_qcdmc->SetBinContent( bi, -9. ) ;
            }// h_qcdmc
	 }//bi
	 gr_model = new TGraphErrors( nb_global_after_exclusion, x, y, ex, ey ) ;
      }
      gr_model -> SetFillColor( kRed-10 ) ;


      TCanvas* can1 = new TCanvas( "can1", "Model vs QCD MC", 1250, 700 ) ;

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      h_model_noerrs -> GetXaxis() -> SetRange( 1, no_bin_njet[0] ) ;

      h_model_noerrs -> Draw( "e" ) ;
      gr_model -> Draw( "0 2" ) ;
      h_model_noerrs -> Draw( "e same" ) ;

      h_qcdmc -> Draw( "same e0" ) ;
      h_qcdmc -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---
      int range_min = 0;
      for ( int bi_nj = 1; bi_nj <= nb_nj; bi_nj++)
      {

         if ( bi_nj > 1 ) {range_min += no_bin_njet[bi_nj-2];}
         h_model_noerrs -> GetXaxis() -> SetRange( range_min+1, range_min+no_bin_njet[bi_nj-1] ) ;
         h_model_noerrs -> SetTitle( "QCD Model H/L ratio, Njet "+num_to_str(bi_nj) );
         h_model_noerrs -> SetMaximum( 1.2 ) ;
         h_model_noerrs -> SetMinimum( -0.1 ) ;

         can1 -> Update() ; can1 -> Draw() ;
         sprintf( fname, "%s/plot-nj%d-full.pdf", output_dir, bi_nj ) ;
         can1 -> SaveAs( fname ) ;

         h_model_noerrs -> SetMaximum( 0.25 ) ;
         h_model_noerrs -> SetMinimum( -0.05 ) ;

         can1 -> Update() ; can1 -> Draw() ;
         sprintf( fname, "%s/plot-nj%d-zoom1.pdf", output_dir, bi_nj ) ;
         can1 -> SaveAs( fname ) ;


         h_model_noerrs -> SetMaximum( 0.06 ) ;
         h_model_noerrs -> SetMinimum( -0.01 ) ;

         can1 -> Update() ; can1 -> Draw() ;
         sprintf( fname, "%s/plot-nj%d-zoom2.pdf", output_dir, bi_nj ) ;
         can1 -> SaveAs( fname ) ;

      }//bi_nj
     //---

      h_model_noerrs -> GetXaxis() -> SetRange( 1, nb_global_after_exclusion ) ;

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-allnj-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_noerrs -> SetMaximum( 0.25 ) ;
      h_model_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-allnj-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_noerrs -> SetMaximum( 0.06 ) ;
      h_model_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-allnj-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;





     //---------------------------------------------------------

      std::vector<TH1F*> h_model_nb;        h_model_nb        . resize(nb_nb);
      std::vector<TH1F*> h_model_nb_noerrs; h_model_nb_noerrs . resize(nb_nb);
      std::vector<TH1F*> h_qcdmc_nb;        h_qcdmc_nb        . resize(nb_nb);

      for ( int bi_nb=0; bi_nb<nb_nb; bi_nb++ ) 
      {
         h_model_nb[bi_nb] = new TH1F( "h_model_nb"+num_to_str(bi_nb), "QCD Model H/L ratio, Nb"+num_to_str(bi_nb), no_bin_bjet[bi_nb], 0.5, no_bin_bjet[bi_nb]+0.5 ) ;
         h_model_nb_noerrs[bi_nb] = new TH1F( "h_model_nb"+num_to_str(bi_nb)+"_noerrs", "QCD Model H/L ratio, Nb"+num_to_str(bi_nb),
                                              no_bin_bjet[bi_nb], 0.5, no_bin_bjet[bi_nb]+0.5 ) ;
         h_qcdmc_nb[bi_nb] = new TH1F( "h_qcdmc_nb"+num_to_str(bi_nb), "QCD MC H/L ratio, Nb"+num_to_str(bi_nb), no_bin_bjet[bi_nb], 0.5, no_bin_bjet[bi_nb]+0.5 ) ;

      }//bi_nj


      int bi_hist(0) ;
      std::vector<int> bi_hist_nb; bi_hist_nb.resize(nb_nb);

      for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=nb_nb; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=nb_htmht-nBinsHT; bi_htmht++ ) {
               if ( is_this_bin_excluded(bi_nj-1, bi_nb-1, bi_htmht-1) ) continue;
               bi_hist++ ;

               float model_val = h_model -> GetBinContent( bi_hist ) ;
               float model_err = h_model -> GetBinError( bi_hist ) ;
               float qcdmc_val = h_qcdmc -> GetBinContent( bi_hist ) ;
               float qcdmc_err = h_qcdmc -> GetBinError( bi_hist ) ;

               char label[100] ;
               sprintf( label, "%s", h_model -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               TH1F* hp_model(0x0) ;
               TH1F* hp_model_noerrs(0x0) ;
               TH1F* hp_qcdmc(0x0) ;

               hp_model = h_model_nb[bi_nb-1] ;
               hp_model_noerrs = h_model_nb_noerrs[bi_nb-1] ;
               hp_qcdmc = h_qcdmc_nb[bi_nb-1] ;

               bi_hist_nb[bi_nb-1]++;

               hp_model -> SetBinContent( bi_hist_nb[bi_nb-1], model_val ) ;
               hp_model -> SetBinError( bi_hist_nb[bi_nb-1], model_err ) ;
               hp_model -> GetXaxis() -> SetBinLabel( bi_hist_nb[bi_nb-1], label ) ;

               hp_model_noerrs -> SetBinContent( bi_hist_nb[bi_nb-1], model_val ) ;
               hp_model_noerrs -> SetBinError( bi_hist_nb[bi_nb-1], 0.00000001 ) ;
               hp_model_noerrs -> GetXaxis() -> SetBinLabel( bi_hist_nb[bi_nb-1], label ) ;

               hp_qcdmc -> SetBinContent( bi_hist_nb[bi_nb-1], qcdmc_val ) ;
               hp_qcdmc -> SetBinError( bi_hist_nb[bi_nb-1], qcdmc_err ) ;
               hp_qcdmc -> GetXaxis() -> SetBinLabel( bi_hist_nb[bi_nb-1], label ) ;

            } // bi_htmht
         } // bi_nb
      } // bi_nj



     //--------
      std::vector<TGraphErrors*> gr_model_nb;
      gr_model_nb.resize(nb_nb);

      for ( int bi_nb=0; bi_nb < nb_nb; bi_nb++ )
      {

         {
            double x[100], y[100], ex[100], ey[100] ;
            for ( int bi=1; bi <= h_model_nb[bi_nb]->GetNbinsX(); bi++ ) 
	    {
               x[bi] = bi ;
               ex[bi] = 0.5 ;
               y[bi] = h_model_nb[bi_nb] -> GetBinContent( bi ) ;
               ey[bi] = h_model_nb[bi_nb] -> GetBinError( bi ) ;
	    }//bi

            for ( int bi=1; bi <= h_qcdmc_nb[bi_nb]->GetNbinsX(); bi++ ) 
            {	 
               if ( h_qcdmc_nb[bi_nb]->GetBinContent( bi ) == 0. && h_qcdmc_nb[bi_nb]->GetBinError( bi) == 0. ) {
                  h_qcdmc_nb[bi_nb]->SetBinContent( bi, -9. ) ;
               }
            } // bi
            gr_model_nb[bi_nb] = new TGraphErrors( nb_global_after_exclusion, x, y, ex, ey ) ;
         }



      gr_model_nb      [bi_nb] -> SetFillColor( kRed-10 ) ;
      h_model_nb_noerrs[bi_nb] -> SetLineColor(4) ;
      h_model_nb_noerrs[bi_nb] -> SetLineWidth(2) ;
      h_qcdmc_nb       [bi_nb] -> SetMarkerStyle(20) ;
      h_model_nb_noerrs[bi_nb] -> GetXaxis() -> LabelsOption("v") ;

      h_model_nb_noerrs[bi_nb] -> SetMaximum( 1.2 ) ;
      h_model_nb_noerrs[bi_nb] -> SetMinimum( -0.1 ) ;

      h_model_nb_noerrs[bi_nb] -> GetXaxis() -> SetRange( 1, 40 ) ;

      h_model_nb_noerrs[bi_nb] -> Draw( "e" ) ;
      gr_model_nb      [bi_nb] -> Draw( "0 2" ) ;
      h_model_nb_noerrs[bi_nb] -> Draw( "e same" ) ;

      h_qcdmc_nb       [bi_nb] -> Draw( "same e0" ) ;
      h_qcdmc_nb       [bi_nb] -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---

      h_model_nb_noerrs[bi_nb] -> SetMaximum( 1.2 ) ;
      h_model_nb_noerrs[bi_nb] -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb%d-full.pdf", output_dir, bi_nb) ;
      can1 -> SaveAs( fname ) ;

      h_model_nb_noerrs[bi_nb] -> SetMaximum( 0.25 ) ;
      h_model_nb_noerrs[bi_nb] -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb%d-zoom1.pdf", output_dir, bi_nb ) ;
      can1 -> SaveAs( fname ) ;


      h_model_nb_noerrs[bi_nb] -> SetMaximum( 0.06 ) ;
      h_model_nb_noerrs[bi_nb] -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb%d-zoom2.pdf", output_dir, bi_nb ) ;
      can1 -> SaveAs( fname ) ;

      }//bi_nb


   } // draw_model_vs_mc

#endif

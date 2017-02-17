#include "TCanvas.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"

#include "../CMS_lumi.C"
#include "../get_hist.h"
#include "../histio.c"
#include "../binning.h"
#include "../num_to_str.h"



   void closure_v5( const char* infile = "outputfiles/model-ratio-hist4.root" ) {

  // Canvas size
  int W = 1200;
  int H = 740;
  int H_ref = 740;
  int W_ref = 800;
  float T = 0.10*H_ref;
  float B = 0.06*H_ref;
  float L = 0.16*W_ref;
  float R = 0.04*W_ref;




      bool do_text_binlabels(false) ;
      setup_bins();
      gStyle -> SetOptStat(0) ;
      gStyle -> SetOptTitle(0) ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( infile ) ;

      TH1F* h_ldp_search_bins_qcdmc = get_hist( "h_ldp_search_bins_qcdmc" ) ;
      TH1F* h_hdp_search_bins_qcdmc = get_hist( "h_hdp_search_bins_qcdmc" ) ;
      TH1F* h_ratio_all = get_hist( "h_ratio_all" ) ;
      TH1F* h_ratio_qcdmc_minus_model = get_hist( "h_ratio_qcdmc_minus_model" ) ;

      int nbins = h_ldp_search_bins_qcdmc->GetNbinsX() ;

      //////TH1F* h_hdp_evts = h_hdp_search_bins_qcdmc ;
      TH1F* h_hdp_evts_orig = h_hdp_search_bins_qcdmc ;
      TH1F* h_hdp_model = new TH1F( "h_hdp_model", "HDP, QCD model", nbins, 0.5, nbins+0.5 ) ;
      TH1F* h_hdp_evts_over_model = new TH1F( "h_hdp_evts_over_model", "HDP, QCD MC events / model", nbins, 0.5, nbins+0.5 ) ;
      TH1F* h_hdp_model_over_model = new TH1F( "h_hdp_model_over_model", "HDP, QCD model / model", nbins, 0.5, nbins+0.5 ) ;

      TH1F* h_hdp_evts = new TH1F( "h_hdp_evts","", nb_global_after_exclusion, 0.5, nb_global_after_exclusion + 0.5 ) ;

      for ( int bi=1; bi<=nbins; bi++ ) {

         h_hdp_evts -> SetBinContent( bi, h_hdp_evts_orig -> GetBinContent( bi ) ) ;
         h_hdp_evts -> SetBinError( bi, h_hdp_evts_orig -> GetBinError( bi ) ) ;

         char label[100] ;
         sprintf( label, "%s", h_ldp_search_bins_qcdmc -> GetXaxis() -> GetBinLabel( bi ) ) ;

         float ldp_val = h_ldp_search_bins_qcdmc -> GetBinContent( bi ) ;
         float ldp_err = h_ldp_search_bins_qcdmc -> GetBinError( bi ) ;

         float hdp_val = h_hdp_search_bins_qcdmc -> GetBinContent( bi ) ;
         float hdp_err = h_hdp_search_bins_qcdmc -> GetBinError( bi ) ;

         float ratio_model_val = h_ratio_all -> GetBinContent( bi ) ;
         float ratio_model_err = h_ratio_all -> GetBinError( bi ) ;

         float ratio_diff_val = h_ratio_qcdmc_minus_model -> GetBinContent( bi ) ;
         float ratio_diff_err = h_ratio_qcdmc_minus_model -> GetBinError( bi ) ;




         ///float ratio_val = ratio_model_val + ratio_diff_val ;
         ///float ratio_err = sqrt( pow( ratio_model_err, 2. ) + pow( ratio_diff_err, 2. ) ) ;

         /////float ratio_val = ratio_model_val  ;
         /////float ratio_err = sqrt( pow( ratio_model_err, 2. ) + pow( ratio_diff_err, 2. ) ) ;

         float ratio_val = ratio_model_val ;
         float ratio_err = ratio_model_err ;





         float hdp_pred_val = ldp_val * ratio_val ;
         float hdp_pred_err(0.) ;
         if ( ldp_val > 0 && ratio_val > 0 ) {
            hdp_pred_err = hdp_pred_val * sqrt( pow( ldp_err/ldp_val, 2. ) + pow( ratio_err/ratio_val, 2. ) ) ;
         }
         h_hdp_model -> SetBinContent( bi, hdp_pred_val ) ;
         h_hdp_model -> SetBinError( bi, hdp_pred_err ) ;
         h_hdp_model -> GetXaxis() -> SetBinLabel( bi, label ) ;

         if ( hdp_pred_val > 0 ) {
            h_hdp_evts_over_model -> SetBinContent( bi, (hdp_val/hdp_pred_val) ) ;
            h_hdp_evts_over_model -> SetBinError( bi, (hdp_err/hdp_pred_val) ) ;
            h_hdp_model_over_model -> SetBinContent( bi, 1. ) ;
            h_hdp_model_over_model -> SetBinError( bi, (hdp_pred_err/hdp_pred_val) ) ;
         }
         if ( do_text_binlabels ) {
            h_hdp_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
            h_hdp_model_over_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
            h_hdp_evts_over_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
         }

      } // bi

      h_hdp_model -> GetXaxis() -> LabelsOption( "v" ) ;

      if ( do_text_binlabels ) {
         h_hdp_evts_over_model  -> GetXaxis() -> LabelsOption( "v" ) ;
         h_hdp_model_over_model -> GetXaxis() -> LabelsOption( "v" ) ;
      }


    //-------------------
      //double ymax(3100.) ;
      //double ymax(18100.) ;
      double ymax(68100.) ;
      double ymin(0.002) ;

      h_hdp_evts -> SetMarkerStyle( 20 ) ;


      int width[10][10] = {};

      for ( int bi_nj=0; bi_nj<nb_nj; bi_nj++ ){
         for ( int bi_nb=0; bi_nb<nb_nb; bi_nb++ ){
            for ( int bi_htmht=3; bi_htmht<nb_htmht; bi_htmht++ ){
              if ( !is_this_bin_excluded(bi_nj, bi_nb, bi_htmht) ) width[bi_nj][bi_nb]++;
            }//bi_htmht
         }//bi_nb
      }//bi_nj



      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "canvas" ) ;
      if ( can1 == 0x0 ) { can1 = new TCanvas( "canvas", "QCD closure", 10, 10, W, H ) ; }
      can1 -> Clear() ;

      TCanvas* canvas = can1 ;

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  canvas->Divide(1, 2);

     //--------

////  TPad* pad_top = new TPad( "pad_top", "", 0., 0.35, 1., 1. ) ;
////  pad_top -> SetTopMargin(0.15) ;
////  pad_top -> SetBottomMargin( 0.0 ) ;
////  pad_top -> SetRightMargin(0.03) ;
////  pad_top -> Draw() ;
////  pad_top -> cd() ;


  TPad* canvas_up = (TPad*) canvas->GetListOfPrimitives()->FindObject("canvas_1");
  TPad* canvas_dw = (TPad*) canvas->GetListOfPrimitives()->FindObject("canvas_2");

  TPad* pad_top = canvas_up ;

  //
  // define the size
  double up_height     = 0.8;  // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.30; // please tune so that the smaller canvas size will work in your environment
  double font_size_dw  = 0.1;  // please tune the font size parameter for bottom figure
  double dw_height     = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.04; // KH, added to put the bottom one closer to the top panel

  //
  // set pad size
  canvas_up->SetPad(0., 1 - up_height,    1., 1.00);
  canvas_dw->SetPad(0., 0.,               1., dw_height+dw_height_offset);
  //
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_up->SetTopMargin(0.12);
  canvas_up->SetLeftMargin(0.1);
  //
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_dw->SetBottomMargin(0.35);
  canvas_dw->SetTopMargin(0);
  canvas_dw->SetLeftMargin(0.1);

      h_hdp_model -> SetLineColor( 4 ) ;
      TH1F* h_hdp_model_copy1 = (TH1F*) h_hdp_model -> Clone( "h_hdp_model_copy1" ) ;
      for ( int bi=1; bi<=h_hdp_model_copy1->GetNbinsX(); bi++ ) { h_hdp_model_copy1 -> SetBinError( bi, 0.0000001 ) ; }
      h_hdp_model -> SetFillColor( kRed-10 ) ;

      //double ymax(3900.) ;
      ///////////double ymax(890.) ;
      h_hdp_evts -> SetMaximum(ymax) ;
      h_hdp_evts -> SetMinimum(ymin) ;
      //h_hdp_evts -> SetTitleSize( 0.075, "y" ) ;
      h_hdp_evts -> SetTitleSize( 0.085, "y" ) ;
      h_hdp_evts -> SetTitleOffset( 0.45, "y" ) ;
      h_hdp_evts -> SetYTitle( "Events" ) ;
      h_hdp_evts -> SetLabelSize( 0.06, "y" ) ;
      h_hdp_evts -> SetMarkerSize( 1.2 ) ;
      h_hdp_evts -> GetYaxis() -> SetTickLength( 0.015 ) ;
      h_hdp_evts -> GetXaxis() -> SetTickLength( 0.015 ) ;
      h_hdp_model_copy1 -> GetYaxis() -> SetTickLength( 0.015 ) ;

      h_hdp_evts -> Draw() ;
      h_hdp_model -> Draw( "E2 same" ) ;
      h_hdp_evts -> Draw( "same" ) ;
      h_hdp_model_copy1 -> Draw( "e same" ) ;
      ///////////h_hdp_model_copy1 -> Draw( "e same axig" ) ;
      ///////////h_hdp_model_copy1 -> Draw( "e same axis" ) ;
      h_hdp_model -> Draw( "same axig" ) ;
      h_hdp_model -> Draw( "same axis" ) ;
      gPad -> SetLogy() ;

      TLine* linenj = new TLine() ;
      linenj -> SetLineStyle(2) ;
      TLine* linenb = new TLine() ;
      linenb -> SetLineStyle(3) ;

      int position = 0;

      for ( int bi_nj=0; bi_nj<nb_nj; bi_nj++ ) {
         for ( int bi_nb=0; bi_nb<nb_nb-1; bi_nb++ ) {
            position += width[bi_nj][bi_nb];
            if ( bi_nb <= bi_nj + 1) linenb -> DrawLine( position+0.5, ymin, position+0.5, exp(0.75*(log(ymax)-log(ymin))+log(ymin)) ) ;
         }
         position += width[bi_nj][nb_nb-1];
         linenj -> DrawLine( position+0.5, ymin, position+0.5, ymax ) ;
      }


      TLatex* njlabels = new TLatex() ;
      njlabels -> SetTextAlign( 21 ) ;
      position = 0;
      for ( int bi_nj=0; bi_nj<nb_nj; bi_nj++ )
      {
         int width_nj = 0;
         for ( int bi_nb=0; bi_nb<nb_nb; bi_nb++ ) { width_nj += width[bi_nj][bi_nb]; }
         TString str = num_to_str( bin_edges_nj[bi_nj]+0.5, 0) + " #leq N_{jet} #leq " + num_to_str(bin_edges_nj[bi_nj+1]-0.5, 0);
	 if ( bin_edges_nj[bi_nj]+0.5 == bin_edges_nj[bi_nj+1]-0.5 ) str = "N_{jet} = " + num_to_str(bin_edges_nj[bi_nj+1]-0.5, 0);
         if ( bi_nj == nb_nj - 1 ) str = "N_{jet} #geq " + num_to_str(bin_edges_nj[bi_nj]+0.5, 0);
	 njlabels -> DrawLatex( position + width_nj/2, exp(0.92*(log(ymax)-log(ymin))+log(ymin)), str) ;
         position += width_nj;
      }//bi_nj


      TLatex* nblabels = new TLatex() ;
      nblabels -> SetTextAlign( 21 ) ;

      position = 0;
      for ( int bi_nb=0; bi_nb<nb_nb; bi_nb++ ) { position += width[0][bi_nb];}

      nblabels -> DrawLatex(  position+3+width[1][0]/2, exp(0.81*(log(ymax)-log(ymin))+log(ymin)), "N_{b-jet}" ) ;
      nblabels -> DrawLatex(           3+width[0][0]/2, exp(0.81*(log(ymax)-log(ymin))+log(ymin)), "N_{b-jet}" ) ;

      int position0 = 0;
      for ( int bi_nb=0; bi_nb<nb_nb; bi_nb++ )
      {
         TString str;
         if ( bi_nb == nb_nb - 1 ) str = "#geq " + num_to_str ( bin_edges_nb[bi_nb]+0.5 , 0);
         else str = num_to_str(bin_edges_nb[bi_nb]+0.5, 0);

	 if ( width[1][bi_nb] != 0. ) 
            nblabels -> DrawLatex( position  + width[1][bi_nb]/2, exp(0.75*(log(ymax)-log(ymin))+log(ymin)), str ) ;
         if ( width[0][bi_nb] != 0. ) 
            nblabels -> DrawLatex( position0 + width[0][bi_nb]/2, exp(0.75*(log(ymax)-log(ymin))+log(ymin)), str ) ;

	 
         position  += width[1][bi_nb];
         position0 += width[0][bi_nb];

      }//bi_nb


      TLegend* legend = new TLegend( 0.61, 0.50, 0.87, 0.75 ) ;
      legend -> SetFillColor(0) ;
      legend -> SetLineColor(1) ;
      legend -> SetBorderSize(1) ;
      legend -> SetHeader( "QCD background" ) ;
      legend -> AddEntry( h_hdp_evts, "Direct from simulation" ) ;
      legend -> AddEntry( h_hdp_model, "Treat simulation like data" ) ;
      legend -> Draw() ;


     //--------

      can1 -> cd() ;
      TPad* pad_bot = new TPad( "pad_bot", "", 0., 0., 1., 0.35 ) ;
      pad_bot -> SetTopMargin( 0.02 ) ;
      pad_bot -> SetBottomMargin( 0.40 ) ;
      pad_bot -> SetRightMargin(0.03) ;
      pad_bot -> Draw() ;
      pad_bot -> cd() ;

      ymax = 5.0 ;
      ymin = 0.0 ;


      h_hdp_evts_over_model -> SetMaximum( ymax ) ;
      h_hdp_evts_over_model -> SetMinimum( ymin ) ;
      h_hdp_evts_over_model -> SetMarkerStyle( 20 ) ;
      h_hdp_evts_over_model -> SetMarkerSize( 1.2 ) ;
      h_hdp_evts_over_model -> SetLabelSize( 0.14, "y" ) ;
      if ( !do_text_binlabels ) h_hdp_evts_over_model -> SetLabelSize( 0.12, "x" ) ;
      h_hdp_evts_over_model -> GetYaxis() -> SetNdivisions( 504 ) ;

      //h_hdp_evts_over_model -> SetYTitle( "Direct/Pred." ) ;
      h_hdp_evts_over_model -> SetYTitle( "#frac{Direct}{Prediction}" ) ;
      //h_hdp_evts_over_model -> SetTitleSize( 0.1, "y" ) ;
      h_hdp_evts_over_model -> SetTitleSize( 0.13, "y" ) ;
      h_hdp_evts_over_model -> SetXTitle( "Search region bin number") ;
      h_hdp_evts_over_model -> SetTitleSize( 0.165, "x" ) ;
      h_hdp_evts_over_model -> GetXaxis() -> SetTickLength( 0.065 ) ;
      h_hdp_evts_over_model -> SetTitleOffset( 0.25, "y" ) ;

      h_hdp_model_over_model -> SetFillColor( kRed-10 ) ;

      TLine* line1 = new TLine() ;
      line1 -> SetLineColor(4) ;

      for ( int bi=1; bi<= h_hdp_evts_over_model -> GetNbinsX() ; bi++ ) {
         if ( h_hdp_evts_over_model-> GetBinContent( bi ) <= 0 ) {
            h_hdp_evts_over_model->SetBinContent( bi, -9. ) ;
         }
      }

      h_hdp_evts_over_model -> Draw("e0") ;
      h_hdp_model_over_model -> Draw( "E2 same" ) ;
      line1 -> DrawLine(0,1.,nb_global_after_exclusion,1.) ;
      h_hdp_evts_over_model -> Draw( "same" ) ;
      h_hdp_evts_over_model -> Draw( "axis same" ) ;

      position = 0;
      for ( int bi_nj=0; bi_nj<nb_nj; bi_nj++ ) {
         for ( int bi_nb=0; bi_nb<nb_nb-1; bi_nb++ ) {
            position += width[bi_nj][bi_nb];
            if ( bi_nb <= bi_nj + 1)  linenb -> DrawLine( position+0.5, ymin, position+0.5, ymax ) ;
         }
         position += width[bi_nj][nb_nb-1];
         linenj -> DrawLine( position+0.5, ymin, position+0.5, ymax ) ;
      }




      TString lumiline( lumi_13TeV + " (13 TeV)" ) ;
      lumi_sqrtS = lumiline ;

      writeExtraText = true;
      extraText = "        Simulation" ;

      CMS_lumi( can1, 0, 0 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      can1 -> SaveAs( "outputfiles/closure-all-v5.pdf" ) ;

   } // closure_v5




#include "CMS_lumi.C"

#include "histio.c"

   TH1F* get_hist( const char* hname ) ;

   void closure_v4( const char* infile = "outputfiles/model-ratio-hist1.root" ) {

      bool do_text_binlabels(false) ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetOptTitle(0) ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( infile ) ;

      TH1F* h_ldp_160bins_qcdmc = get_hist( "h_ldp_160bins_qcdmc" ) ;
      TH1F* h_hdp_160bins_qcdmc = get_hist( "h_hdp_160bins_qcdmc" ) ;
      TH1F* h_ratio_all = get_hist( "h_ratio_all" ) ;
      TH1F* h_ratio_qcdmc_minus_model = get_hist( "h_ratio_qcdmc_minus_model" ) ;

      int nbins = h_ldp_160bins_qcdmc->GetNbinsX() ;

      //////TH1F* h_hdp_evts = h_hdp_160bins_qcdmc ;
      TH1F* h_hdp_evts_orig = h_hdp_160bins_qcdmc ;
      TH1F* h_hdp_model = new TH1F( "h_hdp_model", "HDP, QCD model", nbins, 0.5, nbins+0.5 ) ;
      TH1F* h_hdp_evts_over_model = new TH1F( "h_hdp_evts_over_model", "HDP, QCD MC events / model", nbins, 0.5, nbins+0.5 ) ;
      TH1F* h_hdp_model_over_model = new TH1F( "h_hdp_model_over_model", "HDP, QCD model / model", nbins, 0.5, nbins+0.5 ) ;

      TH1F* h_hdp_evts = new TH1F( "h_hdp_evts","", 160, 0.5, 160.5 ) ;

      for ( int bi=1; bi<=nbins; bi++ ) {

         h_hdp_evts -> SetBinContent( bi, h_hdp_evts_orig -> GetBinContent( bi ) ) ;
         h_hdp_evts -> SetBinError( bi, h_hdp_evts_orig -> GetBinError( bi ) ) ;

         char label[100] ;
         sprintf( label, "%s", h_ldp_160bins_qcdmc -> GetXaxis() -> GetBinLabel( bi ) ) ;

         float ldp_val = h_ldp_160bins_qcdmc -> GetBinContent( bi ) ;
         float ldp_err = h_ldp_160bins_qcdmc -> GetBinError( bi ) ;

         float hdp_val = h_hdp_160bins_qcdmc -> GetBinContent( bi ) ;
         float hdp_err = h_hdp_160bins_qcdmc -> GetBinError( bi ) ;

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
         ///h_hdp_model -> GetXaxis() -> SetBinLabel( bi, label ) ;

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
         h_hdp_evts_over_model -> GetXaxis() -> LabelsOption( "v" ) ;
         h_hdp_model_over_model -> GetXaxis() -> LabelsOption( "v" ) ;
      }


    //-------------------

      int nb_ht(3) ;
      int nb_nj(4) ;
      int nb_nb(4) ;
      int nb_htmht(13) ;
      int nb_mht(5) ;

      int nb_htmhts(10) ;


      h_hdp_evts -> SetMarkerStyle( 20 ) ;


      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1_closure_v4" ) ;
      if ( can1 == 0x0 ) { can1 = new TCanvas( "can1_closure_v4", "QCD closure", 1600, 900 ) ; }
      can1 -> Clear() ;


     //--------

      TPad* pad_top = new TPad( "pad_top", "", 0., 0.35, 1., 1. ) ;
      pad_top -> SetTopMargin(0.15) ;
      pad_top -> SetBottomMargin( 0.0 ) ;
      pad_top -> SetRightMargin(0.03) ;
      pad_top -> Draw() ;
      pad_top -> cd() ;


      h_hdp_model -> SetLineColor( 4 ) ;
      TH1F* h_hdp_model_copy1 = (TH1F*) h_hdp_model -> Clone( "h_hdp_model_copy1" ) ;
      for ( int bi=1; bi<=h_hdp_model_copy1->GetNbinsX(); bi++ ) { h_hdp_model_copy1 -> SetBinError( bi, 0.0000001 ) ; }
      h_hdp_model -> SetFillColor( kRed-10 ) ;

      //double ymax(3900.) ;
      ///////////double ymax(890.) ;
      double ymax(3100.) ;
      double ymin(0.002) ;
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
      for ( int i=1; i<=nb_nj; i++ ) {
         linenj -> DrawLine( i*nb_nb*nb_htmhts + 0.5, ymin, i*nb_nb*nb_htmhts + 0.5, ymax ) ;
         for ( int j=1; j<nb_nb; j++ ) {
            linenb -> DrawLine( (i-1)*nb_nb*nb_htmhts + j*nb_htmhts + 0.5, ymin,  (i-1)*nb_nb*nb_htmhts + j*nb_htmhts + 0.5, exp(0.75*(log(ymax)-log(ymin))+log(ymin)) ) ;
         }
      }

      TLatex* nblabels = new TLatex() ;
      nblabels -> SetTextAlign( 21 ) ;
      ///nblabels -> DrawLatex(  5, exp(0.81*(log(ymax)-log(ymin))+log(ymin)), "N_{b-jet}" ) ;
      nblabels -> DrawLatex(  9, exp(0.81*(log(ymax)-log(ymin))+log(ymin)), "N_{b-jet}" ) ;
      nblabels -> DrawLatex(  5, exp(0.75*(log(ymax)-log(ymin))+log(ymin)), "0" ) ;
      nblabels -> DrawLatex( 15, exp(0.75*(log(ymax)-log(ymin))+log(ymin)), "1" ) ;
      nblabels -> DrawLatex( 25, exp(0.75*(log(ymax)-log(ymin))+log(ymin)), "2" ) ;
      nblabels -> DrawLatex( 35, exp(0.75*(log(ymax)-log(ymin))+log(ymin)), "#geq 3" ) ;

      TLatex* njlabels = new TLatex() ;
      njlabels -> SetTextAlign( 21 ) ;
      njlabels -> DrawLatex(  20, exp(0.92*(log(ymax)-log(ymin))+log(ymin)), "3 #leq N_{jet} #leq 4" ) ;
      njlabels -> DrawLatex(  60, exp(0.92*(log(ymax)-log(ymin))+log(ymin)), "5 #leq N_{jet} #leq 6" ) ;
      njlabels -> DrawLatex( 100, exp(0.92*(log(ymax)-log(ymin))+log(ymin)), "7 #leq N_{jet} #leq 8" ) ;
      njlabels -> DrawLatex( 140, exp(0.92*(log(ymax)-log(ymin))+log(ymin)), "N_{jet} #geq 9" ) ;

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
      line1 -> DrawLine(0,1.,160,1.) ;
      h_hdp_evts_over_model -> Draw( "same" ) ;
      h_hdp_evts_over_model -> Draw( "axis same" ) ;

      for ( int i=1; i<=nb_nj; i++ ) {
         linenj -> DrawLine( i*nb_nb*nb_htmhts + 0.5, ymin, i*nb_nb*nb_htmhts + 0.5, ymax ) ;
         for ( int j=1; j<nb_nb; j++ ) {
            linenb -> DrawLine( (i-1)*nb_nb*nb_htmhts + j*nb_htmhts + 0.5, ymin,  (i-1)*nb_nb*nb_htmhts + j*nb_htmhts + 0.5, ymax ) ;
         }
      }

      TString lumiline( "12.9 fb^{-1} (13 TeV)" ) ;
      lumi_sqrtS = lumiline ;

      writeExtraText = true;
      extraText = "        Simulation" ;

      CMS_lumi( can1, 0, 0 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      can1 -> SaveAs( "outputfiles/closure-all-v4.pdf" ) ;

   } // closure_v4


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


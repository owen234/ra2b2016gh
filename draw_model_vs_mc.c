
#include "histio.c"

   TH1F* get_hist( const char* hname ) ;

   //void draw_model_vs_mc( const char* model_file = "outputfiles/model-ratio-hist1.root",
   //                       const char* qcdmc_file = "outputfiles/qcdmc-ratio-v3.root",
   //                       const char* output_dir = "outputfiles/mc-model-vs-mc" ) {

   void draw_model_vs_mc( const char* model_file = "outputfiles/gci-output.root",
                          const char* qcdmc_file = "outputfiles/qcdmc-ratio-v3.root",
                          const char* output_dir = "outputfiles/model-vs-mc" ) {

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
         double x[160], y[160], ex[160], ey[160] ;
         for ( int bi=1; bi<=160; bi++ ) {
            x[bi] = bi ;
            ex[bi] = 0.5 ;
            y[bi] = h_model -> GetBinContent( bi ) ;
            ey[bi] = h_model -> GetBinError( bi ) ;
            if ( h_qcdmc->GetBinContent( bi ) == 0. && h_qcdmc->GetBinError( bi) == 0. ) {
               h_qcdmc->SetBinContent( bi, -9. ) ;
            }
         } // bi
         gr_model = new TGraphErrors( 160, x, y, ex, ey ) ;
      }
      gr_model -> SetFillColor( kRed-10 ) ;


      TCanvas* can1 = new TCanvas( "can1", "Model vs QCD MC", 1250, 700 ) ;

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      h_model_noerrs -> GetXaxis() -> SetRange( 1, 40 ) ;

      h_model_noerrs -> Draw( "e" ) ;
      gr_model -> Draw( "0 2" ) ;
      h_model_noerrs -> Draw( "e same" ) ;

      h_qcdmc -> Draw( "same e0" ) ;
      h_qcdmc -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj1-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_noerrs -> SetMaximum( 0.25 ) ;
      h_model_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj1-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_noerrs -> SetMaximum( 0.06 ) ;
      h_model_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj1-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


     //---

      h_model_noerrs -> GetXaxis() -> SetRange( 41, 80 ) ;

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj2-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_noerrs -> SetMaximum( 0.25 ) ;
      h_model_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj2-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_noerrs -> SetMaximum( 0.06 ) ;
      h_model_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj2-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


     //---

      h_model_noerrs -> GetXaxis() -> SetRange( 81, 120 ) ;

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj3-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_noerrs -> SetMaximum( 0.25 ) ;
      h_model_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj3-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_noerrs -> SetMaximum( 0.06 ) ;
      h_model_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj3-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


     //---

      h_model_noerrs -> GetXaxis() -> SetRange( 121, 160 ) ;

      h_model_noerrs -> SetMaximum( 1.2 ) ;
      h_model_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj4-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_noerrs -> SetMaximum( 0.25 ) ;
      h_model_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj4-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_noerrs -> SetMaximum( 0.06 ) ;
      h_model_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nj4-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


     //---

      h_model_noerrs -> SetLabelSize( 0.02, "x" ) ;

      h_model_noerrs -> GetXaxis() -> SetRange( 1, 160 ) ;

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

      TH1F* h_model_nb0 = new TH1F( "h_model_nb0", "QCD Model H/L ratio, Nb0", 40, 0.5, 40.5 ) ;
      TH1F* h_model_nb1 = new TH1F( "h_model_nb1", "QCD Model H/L ratio, Nb1", 40, 0.5, 40.5 ) ;
      TH1F* h_model_nb2 = new TH1F( "h_model_nb2", "QCD Model H/L ratio, Nb2", 40, 0.5, 40.5 ) ;
      TH1F* h_model_nb3 = new TH1F( "h_model_nb3", "QCD Model H/L ratio, Nb3", 40, 0.5, 40.5 ) ;

      TH1F* h_model_nb0_noerrs = new TH1F( "h_model_nb0_noerrs", "QCD Model H/L ratio, Nb0", 40, 0.5, 40.5 ) ;
      TH1F* h_model_nb1_noerrs = new TH1F( "h_model_nb1_noerrs", "QCD Model H/L ratio, Nb1", 40, 0.5, 40.5 ) ;
      TH1F* h_model_nb2_noerrs = new TH1F( "h_model_nb2_noerrs", "QCD Model H/L ratio, Nb2", 40, 0.5, 40.5 ) ;
      TH1F* h_model_nb3_noerrs = new TH1F( "h_model_nb3_noerrs", "QCD Model H/L ratio, Nb3", 40, 0.5, 40.5 ) ;

      TH1F* h_qcdmc_nb0 = new TH1F( "h_qcdmc_nb0", "QCD MC H/L ratio, Nb0", 40, 0.5, 40.5 ) ;
      TH1F* h_qcdmc_nb1 = new TH1F( "h_qcdmc_nb1", "QCD MC H/L ratio, Nb1", 40, 0.5, 40.5 ) ;
      TH1F* h_qcdmc_nb2 = new TH1F( "h_qcdmc_nb2", "QCD MC H/L ratio, Nb2", 40, 0.5, 40.5 ) ;
      TH1F* h_qcdmc_nb3 = new TH1F( "h_qcdmc_nb3", "QCD MC H/L ratio, Nb3", 40, 0.5, 40.5 ) ;

      int bi_hist(0) ;
      for ( int bi_nj=1; bi_nj<=4; bi_nj++ ) {
         for ( int bi_nb=1; bi_nb<=4; bi_nb++ ) {
            for ( int bi_htmht=1; bi_htmht<=10; bi_htmht++ ) {

               bi_hist++ ;
               int bi_hist_nb = (bi_nj-1)*10 + bi_htmht ;

               float model_val = h_model -> GetBinContent( bi_hist ) ;
               float model_err = h_model -> GetBinError( bi_hist ) ;
               float qcdmc_val = h_qcdmc -> GetBinContent( bi_hist ) ;
               float qcdmc_err = h_qcdmc -> GetBinError( bi_hist ) ;

               char label[100] ;
               sprintf( label, "%s", h_model -> GetXaxis() -> GetBinLabel( bi_hist ) ) ;

               TH1F* hp_model(0x0) ;
               TH1F* hp_model_noerrs(0x0) ;
               TH1F* hp_qcdmc(0x0) ;

               if ( bi_nb==1 ) {
                  hp_model = h_model_nb0 ;
                  hp_model_noerrs = h_model_nb0_noerrs ;
                  hp_qcdmc = h_qcdmc_nb0 ;
               } else if ( bi_nb==2 ) {
                  hp_model = h_model_nb1 ;
                  hp_model_noerrs = h_model_nb1_noerrs ;
                  hp_qcdmc = h_qcdmc_nb1 ;
               } else if ( bi_nb==3 ) {
                  hp_model = h_model_nb2 ;
                  hp_model_noerrs = h_model_nb2_noerrs ;
                  hp_qcdmc = h_qcdmc_nb2 ;
               } else if ( bi_nb==4 ) {
                  hp_model = h_model_nb3 ;
                  hp_model_noerrs = h_model_nb3_noerrs ;
                  hp_qcdmc = h_qcdmc_nb3 ;
               }

               hp_model -> SetBinContent( bi_hist_nb, model_val ) ;
               hp_model -> SetBinError( bi_hist_nb, model_err ) ;
               hp_model -> GetXaxis() -> SetBinLabel( bi_hist_nb, label ) ;

               hp_model_noerrs -> SetBinContent( bi_hist_nb, model_val ) ;
               hp_model_noerrs -> SetBinError( bi_hist_nb, 0.00000001 ) ;
               hp_model_noerrs -> GetXaxis() -> SetBinLabel( bi_hist_nb, label ) ;

               hp_qcdmc -> SetBinContent( bi_hist_nb, qcdmc_val ) ;
               hp_qcdmc -> SetBinError( bi_hist_nb, qcdmc_err ) ;
               hp_qcdmc -> GetXaxis() -> SetBinLabel( bi_hist_nb, label ) ;

            } // bi_htmht
         } // bi_nb
      } // bi_nj



     //--------

      TGraphErrors* gr_model_nb0(0x0) ;
      {
         double x[40], y[40], ex[40], ey[40] ;
         for ( int bi=1; bi<=40; bi++ ) {
            x[bi] = bi ;
            ex[bi] = 0.5 ;
            y[bi] = h_model_nb0 -> GetBinContent( bi ) ;
            ey[bi] = h_model_nb0 -> GetBinError( bi ) ;
            if ( h_qcdmc_nb0->GetBinContent( bi ) == 0. && h_qcdmc_nb0->GetBinError( bi) == 0. ) {
               h_qcdmc_nb0->SetBinContent( bi, -9. ) ;
            }
         } // bi
         gr_model_nb0 = new TGraphErrors( 40, x, y, ex, ey ) ;
      }
      gr_model_nb0 -> SetFillColor( kRed-10 ) ;
      h_model_nb0_noerrs -> SetLineColor(4) ;
      h_model_nb0_noerrs -> SetLineWidth(2) ;
      h_qcdmc_nb0 -> SetMarkerStyle(20) ;
      h_model_nb0_noerrs -> GetXaxis() -> LabelsOption("v") ;

      h_model_nb0_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb0_noerrs -> SetMinimum( -0.1 ) ;

      h_model_nb0_noerrs -> GetXaxis() -> SetRange( 1, 40 ) ;

      h_model_nb0_noerrs -> Draw( "e" ) ;
      gr_model_nb0 -> Draw( "0 2" ) ;
      h_model_nb0_noerrs -> Draw( "e same" ) ;

      h_qcdmc_nb0 -> Draw( "same e0" ) ;
      h_qcdmc_nb0 -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---

      h_model_nb0_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb0_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb0-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_nb0_noerrs -> SetMaximum( 0.25 ) ;
      h_model_nb0_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb0-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_nb0_noerrs -> SetMaximum( 0.06 ) ;
      h_model_nb0_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb0-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;






     //--------

      TGraphErrors* gr_model_nb1(0x0) ;
      {
         double x[40], y[40], ex[40], ey[40] ;
         for ( int bi=1; bi<=40; bi++ ) {
            x[bi] = bi ;
            ex[bi] = 0.5 ;
            y[bi] = h_model_nb1 -> GetBinContent( bi ) ;
            ey[bi] = h_model_nb1 -> GetBinError( bi ) ;
            if ( h_qcdmc_nb1->GetBinContent( bi ) == 0. && h_qcdmc_nb1->GetBinError( bi) == 0. ) {
               h_qcdmc_nb1->SetBinContent( bi, -9. ) ;
            }
         } // bi
         gr_model_nb1 = new TGraphErrors( 40, x, y, ex, ey ) ;
      }
      gr_model_nb1 -> SetFillColor( kRed-10 ) ;
      h_model_nb1_noerrs -> SetLineColor(4) ;
      h_model_nb1_noerrs -> SetLineWidth(2) ;
      h_qcdmc_nb1 -> SetMarkerStyle(20) ;
      h_model_nb1_noerrs -> GetXaxis() -> LabelsOption("v") ;

      h_model_nb1_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb1_noerrs -> SetMinimum( -0.1 ) ;

      h_model_nb1_noerrs -> GetXaxis() -> SetRange( 1, 40 ) ;

      h_model_nb1_noerrs -> Draw( "e" ) ;
      gr_model_nb1 -> Draw( "0 2" ) ;
      h_model_nb1_noerrs -> Draw( "e same" ) ;

      h_qcdmc_nb1 -> Draw( "same e0" ) ;
      h_qcdmc_nb1 -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---

      h_model_nb1_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb1_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb1-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_nb1_noerrs -> SetMaximum( 0.25 ) ;
      h_model_nb1_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb1-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_nb1_noerrs -> SetMaximum( 0.06 ) ;
      h_model_nb1_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb1-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;



     //--------

      TGraphErrors* gr_model_nb2(0x0) ;
      {
         double x[40], y[40], ex[40], ey[40] ;
         for ( int bi=1; bi<=40; bi++ ) {
            x[bi] = bi ;
            ex[bi] = 0.5 ;
            y[bi] = h_model_nb2 -> GetBinContent( bi ) ;
            ey[bi] = h_model_nb2 -> GetBinError( bi ) ;
            if ( h_qcdmc_nb2->GetBinContent( bi ) == 0. && h_qcdmc_nb2->GetBinError( bi) == 0. ) {
               h_qcdmc_nb2->SetBinContent( bi, -9. ) ;
            }
         } // bi
         gr_model_nb2 = new TGraphErrors( 40, x, y, ex, ey ) ;
      }
      gr_model_nb2 -> SetFillColor( kRed-10 ) ;
      h_model_nb2_noerrs -> SetLineColor(4) ;
      h_model_nb2_noerrs -> SetLineWidth(2) ;
      h_qcdmc_nb2 -> SetMarkerStyle(20) ;
      h_model_nb2_noerrs -> GetXaxis() -> LabelsOption("v") ;

      h_model_nb2_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb2_noerrs -> SetMinimum( -0.1 ) ;

      h_model_nb2_noerrs -> GetXaxis() -> SetRange( 1, 40 ) ;

      h_model_nb2_noerrs -> Draw( "e" ) ;
      gr_model_nb2 -> Draw( "0 2" ) ;
      h_model_nb2_noerrs -> Draw( "e same" ) ;

      h_qcdmc_nb2 -> Draw( "same e0" ) ;
      h_qcdmc_nb2 -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---

      h_model_nb2_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb2_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb2-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_nb2_noerrs -> SetMaximum( 0.25 ) ;
      h_model_nb2_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb2-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_nb2_noerrs -> SetMaximum( 0.06 ) ;
      h_model_nb2_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb2-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;



     //--------

      TGraphErrors* gr_model_nb3(0x0) ;
      {
         double x[40], y[40], ex[40], ey[40] ;
         for ( int bi=1; bi<=40; bi++ ) {
            x[bi] = bi ;
            ex[bi] = 0.5 ;
            y[bi] = h_model_nb3 -> GetBinContent( bi ) ;
            ey[bi] = h_model_nb3 -> GetBinError( bi ) ;
            if ( h_qcdmc_nb3->GetBinContent( bi ) == 0. && h_qcdmc_nb3->GetBinError( bi) == 0. ) {
               h_qcdmc_nb3->SetBinContent( bi, -9. ) ;
            }
         } // bi
         gr_model_nb3 = new TGraphErrors( 40, x, y, ex, ey ) ;
      }
      gr_model_nb3 -> SetFillColor( kRed-10 ) ;
      h_model_nb3_noerrs -> SetLineColor(4) ;
      h_model_nb3_noerrs -> SetLineWidth(2) ;
      h_qcdmc_nb3 -> SetMarkerStyle(20) ;
      h_model_nb3_noerrs -> GetXaxis() -> LabelsOption("v") ;

      h_model_nb3_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb3_noerrs -> SetMinimum( -0.1 ) ;

      h_model_nb3_noerrs -> GetXaxis() -> SetRange( 1, 40 ) ;

      h_model_nb3_noerrs -> Draw( "e" ) ;
      gr_model_nb3 -> Draw( "0 2" ) ;
      h_model_nb3_noerrs -> Draw( "e same" ) ;

      h_qcdmc_nb3 -> Draw( "same e0" ) ;
      h_qcdmc_nb3 -> Draw( "axig same" ) ;

      gPad -> SetGridy(1) ;

     //---

      h_model_nb3_noerrs -> SetMaximum( 1.2 ) ;
      h_model_nb3_noerrs -> SetMinimum( -0.1 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb3-full.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;

      h_model_nb3_noerrs -> SetMaximum( 0.25 ) ;
      h_model_nb3_noerrs -> SetMinimum( -0.05 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb3-zoom1.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;


      h_model_nb3_noerrs -> SetMaximum( 0.06 ) ;
      h_model_nb3_noerrs -> SetMinimum( -0.01 ) ;

      can1 -> Update() ; can1 -> Draw() ;
      sprintf( fname, "%s/plot-nb3-zoom2.pdf", output_dir ) ;
      can1 -> SaveAs( fname ) ;







   } // draw_model_vs_mc


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


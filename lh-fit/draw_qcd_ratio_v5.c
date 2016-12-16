
#include "TDirectory.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"

#include "../get_hist.h"

#include "../histio.c"

#include "../binning.h"

   void  draw_boundaries() ;
   void  draw_boundaries_nj( int nhtb ) ;

   /////////const int nb_nj(4) ;
   /////////const int nb_mht(5) ;
   /////////const int nb_htmht(13) ;

   //---------

   void draw_qcd_ratio_v5( const char* plottype = "nbsum",
                            const char* postfix = "",
                            double rzoom1max = 0.25,
                            double drmax = 6,
                            const char* infile = "../outputfiles/hists-v2d-qcd.root",
                            const char* outputdir = "outputfiles/v5-plots/" ) {

      setup_bins() ;

      printf("\n\n Number of Njets bins : %d\n\n", nb_nj ) ;

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
      if ( strlen( postfix ) > 0 ) {
         sprintf( hname, "h_%s_hdp_%s", plottype, postfix ) ;
      } else {
         sprintf( hname, "h_%s_hdp", plottype ) ;
      }
      TH1F* h_hdp = get_hist( hname ) ;
      if ( strlen( postfix ) > 0 ) {
         sprintf( hname, "h_%s_ldp_%s", plottype, postfix ) ;
      } else {
         sprintf( hname, "h_%s_ldp", plottype ) ;
      }
      TH1F* h_ldp = get_hist( hname ) ;

      TH1F* h_ratio = (TH1F*) h_ldp -> Clone( "h_qcd_ratio" ) ;
      tstring = h_ratio -> GetTitle() ;
      tstring.ReplaceAll( "events", "H/L ratio" ) ;
      tstring.ReplaceAll( "LDP", "" ) ;
      tstring.ReplaceAll( "ldp", "" ) ;
      h_ratio -> SetTitle( tstring ) ;

      for ( int bi=1; bi<=h_ldp->GetNbinsX(); bi++ ) {
         float hdp_val = h_hdp -> GetBinContent( bi ) ;
         float hdp_err = h_hdp -> GetBinError( bi ) ;
         float ldp_val = h_ldp -> GetBinContent( bi ) ;
         float ldp_err = h_ldp -> GetBinError( bi ) ;
         float ratio_val = 0. ;
         float ratio_err = 0. ;
         if ( ldp_val > 0. ) {
            ratio_val = hdp_val / ldp_val ;
            float e2 = pow( ratio_val * ldp_err / ldp_val, 2. )  ;
            if ( hdp_val > 0 ) { e2 += pow( ratio_val * hdp_err / hdp_val, 2. ) ; }
            ratio_err = sqrt( e2 ) ;
         }
         h_ratio -> SetBinContent( bi, ratio_val ) ;
         h_ratio -> SetBinError( bi, ratio_err ) ;
      } // bi


      TH1F* h_double_ratio = (TH1F*) h_ldp -> Clone( "h_double_ratio" ) ;
      h_double_ratio -> Reset() ;
      tstring = h_double_ratio -> GetTitle() ;
      tstring.ReplaceAll( "events", "Double H/L ratio" ) ;
      tstring.ReplaceAll( "LDP", "" ) ;
      tstring.ReplaceAll( "ldp", "" ) ;
      h_double_ratio -> SetTitle( tstring ) ;

      TH1F* h_control_ratio_mht0 = new TH1F( "h_control_ratio_mht0", "Control H/L ratio, MHTC [250,300]", 3*nb_nj, 0.5, 3*nb_nj+0.5 ) ;

      TH1F* h_ratio_mht1 = new TH1F( "h_ratio_mht1", "H/L ratio, MHT1 [300,350]", 3*nb_nj, 0.5, 3*nb_nj+0.5 ) ;
      TH1F* h_ratio_mht2 = new TH1F( "h_ratio_mht2", "H/L ratio, MHT2 [350,500]", 3*nb_nj, 0.5, 3*nb_nj+0.5 ) ;
      TH1F* h_ratio_mht3 = new TH1F( "h_ratio_mht3", "H/L ratio, MHT3 [500,750]", 2*nb_nj, 0.5, 2*nb_nj+0.5 ) ;
      TH1F* h_ratio_mht4 = new TH1F( "h_ratio_mht4", "H/L ratio, MHT4 [750,+]", 2*nb_nj, 0.5, 2*nb_nj+0.5 ) ;

      TH1F* h_dr_mht1 = new TH1F( "h_dr_mht1", "double ratio, MHT1 bins", 3*nb_nj, 0.5, 3*nb_nj+0.5 ) ;
      TH1F* h_dr_mht2 = new TH1F( "h_dr_mht2", "double ratio, MHT2 bins", 3*nb_nj, 0.5, 3*nb_nj+0.5 ) ;
      TH1F* h_dr_mht3 = new TH1F( "h_dr_mht3", "double ratio, MHT3 bins", 2*nb_nj, 0.5, 2*nb_nj+0.5 ) ;
      TH1F* h_dr_mht4 = new TH1F( "h_dr_mht4", "double ratio, MHT4 bins", 2*nb_nj, 0.5, 2*nb_nj+0.5 ) ;


      TH1F* h_ratio_htl = new TH1F( "h_ratio_htl", "H/L ratio, HT low", (nb_mht-2)*nb_nj, 0.5, (nb_mht-2)*nb_nj+0.5 ) ;
      TH1F* h_ratio_htm = new TH1F( "h_ratio_htm", "H/L ratio, HT middle", nb_mht*nb_nj, 0.5, nb_mht*nb_nj+0.5 ) ;
      TH1F* h_ratio_hth = new TH1F( "h_ratio_hth", "H/L ratio, HT high", nb_mht*nb_nj, 0.5, nb_mht*nb_nj+0.5 ) ;

      TH1F* h_dr_htl = new TH1F( "h_dr_htl", "H/L double ratio, HT low", (nb_mht-2)*nb_nj, 0.5, (nb_mht-2)*nb_nj+0.5 ) ;
      TH1F* h_dr_htm = new TH1F( "h_dr_htm", "H/L double ratio, HT middle", nb_mht*nb_nj, 0.5, nb_mht*nb_nj+0.5 ) ;
      TH1F* h_dr_hth = new TH1F( "h_dr_hth", "H/L double ratio, HT high", nb_mht*nb_nj, 0.5, nb_mht*nb_nj+0.5 ) ;

      TH1F* h_htm_over_hth_ratio = new TH1F( "h_htm_over_hth_ratio", "middle HT H/L ratio over high HT H/L ratio", nb_mht*nb_nj, 0.5, nb_mht*nb_nj+0.5 ) ;


      for ( int nji=1; nji<=nb_nj; nji++ ) {

         int boff = (nji-1)*nb_htmht ;

         float ht1_val = h_ratio -> GetBinContent( boff+1 ) ;
         float ht1_err = h_ratio -> GetBinError( boff+1 ) ;
         float ht2_val = h_ratio -> GetBinContent( boff+2 ) ;
         float ht2_err = h_ratio -> GetBinError( boff+2 ) ;
         float ht3_val = h_ratio -> GetBinContent( boff+3 ) ;
         float ht3_err = h_ratio -> GetBinError( boff+3 ) ;

         h_control_ratio_mht0 -> SetBinContent( (nji-1)*3 + 1, ht1_val ) ;
         h_control_ratio_mht0 -> SetBinContent( (nji-1)*3 + 2, ht2_val ) ;
         h_control_ratio_mht0 -> SetBinContent( (nji-1)*3 + 3, ht3_val ) ;

         h_control_ratio_mht0 -> SetBinError( (nji-1)*3 + 1, ht1_err ) ;
         h_control_ratio_mht0 -> SetBinError( (nji-1)*3 + 2, ht2_err ) ;
         h_control_ratio_mht0 -> SetBinError( (nji-1)*3 + 3, ht3_err ) ;

         h_control_ratio_mht0 -> GetXaxis() -> SetBinLabel( (nji-1)*3 + 1, h_ratio -> GetXaxis() -> GetBinLabel( boff+1 ) ) ;
         h_control_ratio_mht0 -> GetXaxis() -> SetBinLabel( (nji-1)*3 + 2, h_ratio -> GetXaxis() -> GetBinLabel( boff+2 ) ) ;
         h_control_ratio_mht0 -> GetXaxis() -> SetBinLabel( (nji-1)*3 + 3, h_ratio -> GetXaxis() -> GetBinLabel( boff+3 ) ) ;

         h_ratio_htl -> SetBinContent( (nji-1)*(nb_mht-2) + 1, ht1_val ) ;
         h_ratio_htm -> SetBinContent( (nji-1)*(nb_mht) + 1, ht2_val ) ;
         h_ratio_hth -> SetBinContent( (nji-1)*(nb_mht) + 1, ht3_val ) ;

         h_ratio_htl -> SetBinError( (nji-1)*(nb_mht-2) + 1, ht1_err ) ;
         h_ratio_htm -> SetBinError( (nji-1)*(nb_mht) + 1, ht2_err ) ;
         h_ratio_hth -> SetBinError( (nji-1)*(nb_mht) + 1, ht3_err ) ;

         h_ratio_htl -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht-2) + 1, h_ratio -> GetXaxis() -> GetBinLabel( boff+1 ) ) ;
         h_ratio_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 1, h_ratio -> GetXaxis() -> GetBinLabel( boff+2 ) ) ;
         h_ratio_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 1, h_ratio -> GetXaxis() -> GetBinLabel( boff+3 ) ) ;


      //+++++++++++++++
         float b1_val = h_ratio -> GetBinContent( boff+4 ) ;
         float b1_err = h_ratio -> GetBinError( boff+4 ) ;
         float dr1_val(0.) ;
         float dr1_err(0.) ;
         if ( ht1_val > 0. ) {
            dr1_val = b1_val / ht1_val ;
            if ( b1_val > 0. ) {
               dr1_err = dr1_val * sqrt( pow( b1_err/b1_val, 2. ) + pow( ht1_err/ht1_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+4, dr1_val ) ;
         h_double_ratio -> SetBinError( boff+4, dr1_err ) ;
         h_dr_mht1 -> SetBinContent( (nji-1)*3+1, dr1_val ) ;
         h_dr_mht1 -> SetBinError( (nji-1)*3+1, dr1_err ) ;
         h_dr_mht1 -> GetXaxis() -> SetBinLabel( (nji-1)*3+1, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+4 ) ) ;
         h_ratio_mht1 -> SetBinContent( (nji-1)*3 + 1, b1_val ) ;
         h_ratio_mht1 -> SetBinError( (nji-1)*3 + 1, b1_err ) ;
         h_ratio_mht1 -> GetXaxis() -> SetBinLabel( (nji-1)*3+1, h_ratio -> GetXaxis() -> GetBinLabel( boff+4 ) ) ;
         h_dr_htl -> SetBinContent( (nji-1)*(nb_mht-2) + 2, dr1_val ) ;
         h_dr_htl -> SetBinError( (nji-1)*(nb_mht-2) + 2, dr1_err ) ;
         h_dr_htl -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht-2) + 2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+4 ) ) ;
         h_ratio_htl -> SetBinContent( (nji-1)*(nb_mht-2) + 2, b1_val ) ;
         h_ratio_htl -> SetBinError( (nji-1)*(nb_mht-2) + 2, b1_err ) ;
         h_ratio_htl -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht-2) + 2, h_ratio -> GetXaxis() -> GetBinLabel( boff+4 ) ) ;
       //------
         float b2_val = h_ratio -> GetBinContent( boff+5 ) ;
         float b2_err = h_ratio -> GetBinError( boff+5 ) ;
         float dr2_val(0.) ;
         float dr2_err(0.) ;
         if ( ht2_val > 0. ) {
            dr2_val = b2_val / ht2_val ;
            if ( b2_val > 0. ) {
               dr2_err = dr2_val * sqrt( pow( b2_err/b2_val, 2. ) + pow( ht2_err/ht2_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+5, dr2_val ) ;
         h_double_ratio -> SetBinError( boff+5, dr2_err ) ;
         h_dr_mht1 -> SetBinContent( (nji-1)*3+2, dr2_val ) ;
         h_dr_mht1 -> SetBinError( (nji-1)*3+2, dr2_err ) ;
         h_dr_mht1 -> GetXaxis() -> SetBinLabel( (nji-1)*3+2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+5 ) ) ;
         h_ratio_mht1 -> SetBinContent( (nji-1)*3 + 2, b2_val ) ;
         h_ratio_mht1 -> SetBinError( (nji-1)*3 + 2, b2_err ) ;
         h_ratio_mht1 -> GetXaxis() -> SetBinLabel( (nji-1)*3+2, h_ratio -> GetXaxis() -> GetBinLabel( boff+5 ) ) ;
         h_dr_htm -> SetBinContent( (nji-1)*(nb_mht) + 2, dr2_val ) ;
         h_dr_htm -> SetBinError( (nji-1)*(nb_mht) + 2, dr2_err ) ;
         h_dr_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+5 ) ) ;
         h_ratio_htm -> SetBinContent( (nji-1)*(nb_mht) + 2, b2_val ) ;
         h_ratio_htm -> SetBinError( (nji-1)*(nb_mht) + 2, b2_err ) ;
         h_ratio_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 2, h_ratio -> GetXaxis() -> GetBinLabel( boff+5 ) ) ;
       //------
         float b3_val = h_ratio -> GetBinContent( boff+6 ) ;
         float b3_err = h_ratio -> GetBinError( boff+6 ) ;
         float dr3_val(0.) ;
         float dr3_err(0.) ;
         if ( ht3_val > 0. ) {
            dr3_val = b3_val / ht3_val ;
            if ( b3_val > 0. ) {
               dr3_err = dr3_val * sqrt( pow( b3_err/b3_val, 2. ) + pow( ht3_err/ht3_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+6, dr3_val ) ;
         h_double_ratio -> SetBinError( boff+6, dr3_err ) ;
         h_dr_mht1 -> SetBinContent( (nji-1)*3+3, dr3_val ) ;
         h_dr_mht1 -> SetBinError( (nji-1)*3+3, dr3_err ) ;
         h_dr_mht1 -> GetXaxis() -> SetBinLabel( (nji-1)*3+3, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+6 ) ) ;
         h_ratio_mht1 -> SetBinContent( (nji-1)*3 + 3, b3_val ) ;
         h_ratio_mht1 -> SetBinError( (nji-1)*3 + 3, b3_err ) ;
         h_ratio_mht1 -> GetXaxis() -> SetBinLabel( (nji-1)*3+3, h_ratio -> GetXaxis() -> GetBinLabel( boff+6 ) ) ;
         h_dr_hth -> SetBinContent( (nji-1)*(nb_mht) + 2, dr3_val ) ;
         h_dr_hth -> SetBinError( (nji-1)*(nb_mht) + 2, dr3_err ) ;
         h_dr_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+6 ) ) ;
         h_ratio_hth -> SetBinContent( (nji-1)*(nb_mht) + 2, b3_val ) ;
         h_ratio_hth -> SetBinError( (nji-1)*(nb_mht) + 2, b3_err ) ;
         h_ratio_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 2, h_ratio -> GetXaxis() -> GetBinLabel( boff+6 ) ) ;
       //------


      //+++++++++++++++
         float b4_val = h_ratio -> GetBinContent( boff+7 ) ;
         float b4_err = h_ratio -> GetBinError( boff+7 ) ;
         float dr4_val(0.) ;
         float dr4_err(0.) ;
         if ( ht1_val > 0. ) {
            dr4_val = b4_val / ht1_val ;
            if ( b4_val > 0. ) {
               dr4_err = dr4_val * sqrt( pow( b4_err/b4_val, 2. ) + pow( ht1_err/ht1_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+7, dr4_val ) ;
         h_double_ratio -> SetBinError( boff+7, dr4_err ) ;
         h_dr_mht2 -> SetBinContent( (nji-1)*3+1, dr4_val ) ;
         h_dr_mht2 -> SetBinError( (nji-1)*3+1, dr4_err ) ;
         h_dr_mht2 -> GetXaxis() -> SetBinLabel( (nji-1)*3+1, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+7 ) ) ;
         h_ratio_mht2 -> SetBinContent( (nji-1)*3 + 1, b4_val ) ;
         h_ratio_mht2 -> SetBinError( (nji-1)*3 + 1, b4_err ) ;
         h_ratio_mht2 -> GetXaxis() -> SetBinLabel( (nji-1)*3+1, h_ratio -> GetXaxis() -> GetBinLabel( boff+7 ) ) ;
         h_dr_htl -> SetBinContent( (nji-1)*(nb_mht-2) + 3, dr4_val ) ;
         h_dr_htl -> SetBinError( (nji-1)*(nb_mht-2) + 3, dr4_err ) ;
         h_dr_htl -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht-2) + 3, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+7 ) ) ;
         h_ratio_htl -> SetBinContent( (nji-1)*(nb_mht-2) + 3, b4_val ) ;
         h_ratio_htl -> SetBinError( (nji-1)*(nb_mht-2) + 3, b4_err ) ;
         h_ratio_htl -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht-2) + 3, h_ratio -> GetXaxis() -> GetBinLabel( boff+7 ) ) ;
       //------
         float b5_val = h_ratio -> GetBinContent( boff+8 ) ;
         float b5_err = h_ratio -> GetBinError( boff+8 ) ;
         float dr5_val(0.) ;
         float dr5_err(0.) ;
         if ( ht2_val > 0. ) {
            dr5_val = b5_val / ht2_val ;
            if ( b5_val > 0. ) {
               dr5_err = dr5_val * sqrt( pow( b5_err/b5_val, 2. ) + pow( ht2_err/ht2_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+8, dr5_val ) ;
         h_double_ratio -> SetBinError( boff+8, dr5_err ) ;
         h_dr_mht2 -> SetBinContent( (nji-1)*3+2, dr5_val ) ;
         h_dr_mht2 -> SetBinError( (nji-1)*3+2, dr5_err ) ;
         h_dr_mht2 -> GetXaxis() -> SetBinLabel( (nji-1)*3+2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+8 ) ) ;
         h_ratio_mht2 -> SetBinContent( (nji-1)*3 + 2, b5_val ) ;
         h_ratio_mht2 -> SetBinError( (nji-1)*3 + 2, b5_err ) ;
         h_ratio_mht2 -> GetXaxis() -> SetBinLabel( (nji-1)*3+2, h_ratio -> GetXaxis() -> GetBinLabel( boff+8 ) ) ;
         h_dr_htm -> SetBinContent( (nji-1)*(nb_mht) + 3, dr5_val ) ;
         h_dr_htm -> SetBinError( (nji-1)*(nb_mht) + 3, dr5_err ) ;
         h_dr_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 3, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+8 ) ) ;
         h_ratio_htm -> SetBinContent( (nji-1)*(nb_mht) + 3, b5_val ) ;
         h_ratio_htm -> SetBinError( (nji-1)*(nb_mht) + 3, b5_err ) ;
         h_ratio_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 3, h_ratio -> GetXaxis() -> GetBinLabel( boff+8 ) ) ;
       //------
         float b6_val = h_ratio -> GetBinContent( boff+9 ) ;
         float b6_err = h_ratio -> GetBinError( boff+9 ) ;
         float dr6_val(0.) ;
         float dr6_err(0.) ;
         if ( ht3_val > 0. ) {
            dr6_val = b6_val / ht3_val ;
            if ( b6_val > 0. ) {
               dr6_err = dr6_val * sqrt( pow( b6_err/b6_val, 2. ) + pow( ht3_err/ht3_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+9, dr6_val ) ;
         h_double_ratio -> SetBinError( boff+9, dr6_err ) ;
         h_dr_mht2 -> SetBinContent( (nji-1)*3+3, dr6_val ) ;
         h_dr_mht2 -> SetBinError( (nji-1)*3+3, dr6_err ) ;
         h_dr_mht2 -> GetXaxis() -> SetBinLabel( (nji-1)*3+3, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+9 ) ) ;
         h_ratio_mht2 -> SetBinContent( (nji-1)*3 + 3, b6_val ) ;
         h_ratio_mht2 -> SetBinError( (nji-1)*3 + 3, b6_err ) ;
         h_ratio_mht2 -> GetXaxis() -> SetBinLabel( (nji-1)*3+3, h_ratio -> GetXaxis() -> GetBinLabel( boff+9 ) ) ;
         h_dr_hth -> SetBinContent( (nji-1)*(nb_mht) + 3, dr6_val ) ;
         h_dr_hth -> SetBinError( (nji-1)*(nb_mht) + 3, dr6_err ) ;
         h_dr_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 3, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+9 ) ) ;
         h_ratio_hth -> SetBinContent( (nji-1)*(nb_mht) + 3, b6_val ) ;
         h_ratio_hth -> SetBinError( (nji-1)*(nb_mht) + 3, b6_err ) ;
         h_ratio_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 3, h_ratio -> GetXaxis() -> GetBinLabel( boff+9 ) ) ;
       //------

      //+++++++++++++++
         float b7_val = h_ratio -> GetBinContent( boff+10 ) ;
         float b7_err = h_ratio -> GetBinError( boff+10 ) ;
         float dr7_val(0.) ;
         float dr7_err(0.) ;
         if ( ht2_val > 0. ) {
            dr7_val = b7_val / ht2_val ;
            if ( b7_val > 0. ) {
               dr7_err = dr7_val * sqrt( pow( b7_err/b7_val, 2. ) + pow( ht2_err/ht2_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+10, dr7_val ) ;
         h_double_ratio -> SetBinError( boff+10, dr7_err ) ;
         h_dr_mht3 -> SetBinContent( (nji-1)*2+1, dr7_val ) ;
         h_dr_mht3 -> SetBinError( (nji-1)*2+1, dr7_err ) ;
         h_dr_mht3 -> GetXaxis() -> SetBinLabel( (nji-1)*2+1, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+10 ) ) ;
         h_ratio_mht3 -> SetBinContent( (nji-1)*2 + 1, b7_val ) ;
         h_ratio_mht3 -> SetBinError( (nji-1)*2 + 1, b7_err ) ;
         h_ratio_mht3 -> GetXaxis() -> SetBinLabel( (nji-1)*2+1, h_ratio -> GetXaxis() -> GetBinLabel( boff+10 ) ) ;
         h_dr_htm -> SetBinContent( (nji-1)*(nb_mht) + 4, dr7_val ) ;
         h_dr_htm -> SetBinError( (nji-1)*(nb_mht) + 4, dr7_err ) ;
         h_dr_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 4, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+10 ) ) ;
         h_ratio_htm -> SetBinContent( (nji-1)*(nb_mht) + 4, b7_val ) ;
         h_ratio_htm -> SetBinError( (nji-1)*(nb_mht) + 4, b7_err ) ;
         h_ratio_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 4, h_ratio -> GetXaxis() -> GetBinLabel( boff+10 ) ) ;
       //------
         float b8_val = h_ratio -> GetBinContent( boff+11 ) ;
         float b8_err = h_ratio -> GetBinError( boff+11 ) ;
         float dr8_val(0.) ;
         float dr8_err(0.) ;
         if ( ht3_val > 0. ) {
            dr8_val = b8_val / ht3_val ;
            if ( b8_val > 0. ) {
               dr8_err = dr8_val * sqrt( pow( b8_err/b8_val, 2. ) + pow( ht3_err/ht3_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+11, dr8_val ) ;
         h_double_ratio -> SetBinError( boff+11, dr8_err ) ;
         h_dr_mht3 -> SetBinContent( (nji-1)*2+2, dr8_val ) ;
         h_dr_mht3 -> SetBinError( (nji-1)*2+2, dr8_err ) ;
         h_dr_mht3 -> GetXaxis() -> SetBinLabel( (nji-1)*2+2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+11 ) ) ;
         h_ratio_mht3 -> SetBinContent( (nji-1)*2 + 2, b8_val ) ;
         h_ratio_mht3 -> SetBinError( (nji-1)*2 + 2, b8_err ) ;
         h_ratio_mht3 -> GetXaxis() -> SetBinLabel( (nji-1)*2+2, h_ratio -> GetXaxis() -> GetBinLabel( boff+11 ) ) ;
         h_dr_hth -> SetBinContent( (nji-1)*(nb_mht) + 4, dr8_val ) ;
         h_dr_hth -> SetBinError( (nji-1)*(nb_mht) + 4, dr8_err ) ;
         h_dr_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 4, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+11 ) ) ;
         h_ratio_hth -> SetBinContent( (nji-1)*(nb_mht) + 4, b8_val ) ;
         h_ratio_hth -> SetBinError( (nji-1)*(nb_mht) + 4, b8_err ) ;
         h_ratio_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 4, h_ratio -> GetXaxis() -> GetBinLabel( boff+11 ) ) ;

      //+++++++++++++++
         float b9_val = h_ratio -> GetBinContent( boff+12 ) ;
         float b9_err = h_ratio -> GetBinError( boff+12 ) ;
         float dr9_val(0.) ;
         float dr9_err(0.) ;
         if ( ht2_val > 0. ) {
            dr9_val = b9_val / ht2_val ;
            if ( b9_val > 0. ) {
               dr9_err = dr9_val * sqrt( pow( b9_err/b9_val, 2. ) + pow( ht2_err/ht2_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+12, dr9_val ) ;
         h_double_ratio -> SetBinError( boff+12, dr9_err ) ;
         h_dr_mht4 -> SetBinContent( (nji-1)*2+1, dr9_val ) ;
         h_dr_mht4 -> SetBinError( (nji-1)*2+1, dr9_err ) ;
         h_dr_mht4 -> GetXaxis() -> SetBinLabel( (nji-1)*2+1, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+12 ) ) ;
         h_ratio_mht4 -> SetBinContent( (nji-1)*2 + 1, b9_val ) ;
         h_ratio_mht4 -> SetBinError( (nji-1)*2 + 1, b9_err ) ;
         h_ratio_mht4 -> GetXaxis() -> SetBinLabel( (nji-1)*2+1, h_ratio -> GetXaxis() -> GetBinLabel( boff+12 ) ) ;
         h_dr_htm -> SetBinContent( (nji-1)*(nb_mht) + 5, dr9_val ) ;
         h_dr_htm -> SetBinError( (nji-1)*(nb_mht) + 5, dr9_err ) ;
         h_dr_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 5, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+12 ) ) ;
         h_ratio_htm -> SetBinContent( (nji-1)*(nb_mht) + 5, b9_val ) ;
         h_ratio_htm -> SetBinError( (nji-1)*(nb_mht) + 5, b9_err ) ;
         h_ratio_htm -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 5, h_ratio -> GetXaxis() -> GetBinLabel( boff+12 ) ) ;
       //------
         float b10_val = h_ratio -> GetBinContent( boff+13 ) ;
         float b10_err = h_ratio -> GetBinError( boff+13 ) ;
         float dr10_val(0.) ;
         float dr10_err(0.) ;
         if ( ht3_val > 0. ) {
            dr10_val = b10_val / ht3_val ;
            if ( b10_val > 0. ) {
               dr10_err = dr10_val * sqrt( pow( b10_err/b10_val, 2. ) + pow( ht3_err/ht3_val, 2. ) ) ;
            }
         }
         h_double_ratio -> SetBinContent( boff+13, dr10_val ) ;
         h_double_ratio -> SetBinError( boff+13, dr10_err ) ;
         h_dr_mht4 -> SetBinContent( (nji-1)*2+2, dr10_val ) ;
         h_dr_mht4 -> SetBinError( (nji-1)*2+2, dr10_err ) ;
         h_dr_mht4 -> GetXaxis() -> SetBinLabel( (nji-1)*2+2, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+13 ) ) ;
         h_ratio_mht4 -> SetBinContent( (nji-1)*2 + 2, b10_val ) ;
         h_ratio_mht4 -> SetBinError( (nji-1)*2 + 2, b10_err ) ;
         h_ratio_mht4 -> GetXaxis() -> SetBinLabel( (nji-1)*2+2, h_ratio -> GetXaxis() -> GetBinLabel( boff+13 ) ) ;
         h_dr_hth -> SetBinContent( (nji-1)*(nb_mht) + 5, dr10_val ) ;
         h_dr_hth -> SetBinError( (nji-1)*(nb_mht) + 5, dr10_err ) ;
         h_dr_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 5, h_double_ratio -> GetXaxis() -> GetBinLabel( boff+13 ) ) ;
         h_ratio_hth -> SetBinContent( (nji-1)*(nb_mht) + 5, b10_val ) ;
         h_ratio_hth -> SetBinError( (nji-1)*(nb_mht) + 5, b10_err ) ;
         h_ratio_hth -> GetXaxis() -> SetBinLabel( (nji-1)*(nb_mht) + 5, h_ratio -> GetXaxis() -> GetBinLabel( boff+13 ) ) ;

      } // nji


      printf("\n\n") ;
      for ( int bi=1; bi<=h_ratio_hth-> GetNbinsX(); bi++ ) {
         double htm_val = h_ratio_htm -> GetBinContent( bi ) ;
         double htm_err = h_ratio_htm -> GetBinError( bi ) ;
         double hth_val = h_ratio_hth -> GetBinContent( bi ) ;
         double hth_err = h_ratio_hth -> GetBinError( bi ) ;
         double r_val(0.), r_err(0.) ;
         if ( hth_val > 0 && htm_val > 0 ) {
            r_val = htm_val / hth_val ;
            r_err = r_val * sqrt( pow( htm_err/htm_val, 2 ) + pow( hth_err/hth_val, 2. ) ) ;
         }
         printf("  Med HT H/L ratio : %s : %7.4f +/- %7.4f ,    High HT H/L ratio : %s : %7.4f +/- %7.4f ;  ratio %7.4f +/- %7.4f\n",
            h_ratio_htm -> GetXaxis() -> GetBinLabel(bi), htm_val, htm_err,
            h_ratio_hth -> GetXaxis() -> GetBinLabel(bi), hth_val, hth_err,
            r_val, r_err) ;
         h_htm_over_hth_ratio -> SetBinContent( bi, r_val ) ;
         h_htm_over_hth_ratio -> SetBinError( bi, r_err ) ;
         TString blabel = h_ratio_hth -> GetXaxis() -> GetBinLabel( bi ) ;
         TString blabel2( blabel(0,8) ) ;
         h_htm_over_hth_ratio -> GetXaxis() -> SetBinLabel( bi, blabel2 ) ;
      } // bi
      printf("\n\n") ;
      h_htm_over_hth_ratio -> GetXaxis() -> LabelsOption("v") ;

    //------ Begin: Compute average over Njets for double ratios

      TH1F* h_dr_htl_njave = new TH1F( "h_dr_htl_njave", "H/L double ratio, average over Njets, HT low", nb_mht-2, 0.5, nb_mht-2+0.5 ) ;
      TH1F* h_dr_htm_njave = new TH1F( "h_dr_htm_njave", "H/L double ratio, average over Njets, HT medium", nb_mht, 0.5, nb_mht+0.5 ) ;
      TH1F* h_dr_hth_njave = new TH1F( "h_dr_hth_njave", "H/L double ratio, average over Njets, HT high", nb_mht, 0.5, nb_mht+0.5 ) ;

      {
         printf("\n") ;
         for ( int mbi=2; mbi<=(nb_mht-2); mbi++ ) {
            double sumofw(0.) ;
            double wsum(0.) ;
            for ( int nji=1; nji<=nb_nj; nji++ ) {
               double val = h_dr_htl -> GetBinContent( (nji-1)*(nb_mht-2) + mbi ) ;
               double err = h_dr_htl -> GetBinError(   (nji-1)*(nb_mht-2) + mbi ) ;
               if ( err > 0 ) {
                  double w = 1./(err*err) ;
                  sumofw += w ;
                  wsum += val * w ;
                  printf("   HTL:  %s : nj=%d, mht=%d :  DR = %7.4f +/- %7.4f\n",
                      h_dr_htl -> GetXaxis() -> GetBinLabel( (nji-1)*(nb_mht-2) + mbi ) ,
                      nji, mbi,
                      val, err ) ;
               }
            } // nji
            double ave_val = 0. ;
            double ave_err = 0. ;
            if ( sumofw > 0 ) {
               ave_val = wsum / sumofw ;
               ave_err = 1./sqrt( sumofw ) ;
               printf("   HTL: njave, mht=%d : DR = %7.4f +/- %7.4f\n", mbi, ave_val, ave_err ) ;
            }
            h_dr_htl_njave -> SetBinContent( mbi, ave_val ) ;
            h_dr_htl_njave -> SetBinError( mbi, ave_err ) ;
            char blabel[100] ;
            sprintf( blabel, "LowHT-NjAve-MHT%d", mbi-1 ) ;
            h_dr_htl_njave -> GetXaxis() -> SetBinLabel( mbi, blabel ) ;
         } // mbi.
      }

      {
         printf("\n") ;
         for ( int mbi=2; mbi<=nb_mht; mbi++ ) {
            double sumofw(0.) ;
            double wsum(0.) ;
            for ( int nji=1; nji<=nb_nj; nji++ ) {
               double val = h_dr_htm -> GetBinContent( (nji-1)*(nb_mht) + mbi ) ;
               double err = h_dr_htm -> GetBinError(   (nji-1)*(nb_mht) + mbi ) ;
               if ( err > 0 ) {
                  double w = 1./(err*err) ;
                  sumofw += w ;
                  wsum += val * w ;
                  printf("   HTM:  %s : nj=%d, mht=%d :  DR = %7.4f +/- %7.4f\n",
                      h_dr_htm -> GetXaxis() -> GetBinLabel( (nji-1)*(nb_mht) + mbi ) ,
                      nji, mbi,
                      val, err ) ;
               }
            } // nji
            double ave_val = 0. ;
            double ave_err = 0. ;
            if ( sumofw > 0 ) {
               ave_val = wsum / sumofw ;
               ave_err = 1./sqrt( sumofw ) ;
               printf("   HTM: njave, mht=%d : DR = %7.4f +/- %7.4f\n", mbi, ave_val, ave_err ) ;
            }
            h_dr_htm_njave -> SetBinContent( mbi, ave_val ) ;
            h_dr_htm_njave -> SetBinError( mbi, ave_err ) ;
            char blabel[100] ;
            sprintf( blabel, "MediumHT-NjAve-MHT%d", mbi-1 ) ;
            h_dr_htm_njave -> GetXaxis() -> SetBinLabel( mbi, blabel ) ;
         } // mbi.
      }

      {
         printf("\n") ;
         for ( int mbi=2; mbi<=nb_mht; mbi++ ) {
            double sumofw(0.) ;
            double wsum(0.) ;
            for ( int nji=1; nji<=nb_nj; nji++ ) {
               double val = h_dr_hth -> GetBinContent( (nji-1)*(nb_mht) + mbi ) ;
               double err = h_dr_hth -> GetBinError(   (nji-1)*(nb_mht) + mbi ) ;
               if ( err > 0 ) {
                  double w = 1./(err*err) ;
                  sumofw += w ;
                  wsum += val * w ;
                  printf("   HTH:  %s : nj=%d, mht=%d :  DR = %7.4f +/- %7.4f\n",
                      h_dr_hth -> GetXaxis() -> GetBinLabel( (nji-1)*(nb_mht) + mbi ) ,
                      nji, mbi,
                      val, err ) ;
               }
            } // nji
            double ave_val = 0. ;
            double ave_err = 0. ;
            if ( sumofw > 0 ) {
               ave_val = wsum / sumofw ;
               ave_err = 1./sqrt( sumofw ) ;
               printf("   HTH: njave, mht=%d : DR = %7.4f +/- %7.4f\n", mbi, ave_val, ave_err ) ;
            }
            h_dr_hth_njave -> SetBinContent( mbi, ave_val ) ;
            h_dr_hth_njave -> SetBinError( mbi, ave_err ) ;
            char blabel[100] ;
            sprintf( blabel, "HighHT-NjAve-MHT%d", mbi-1 ) ;
            h_dr_hth_njave -> GetXaxis() -> SetBinLabel( mbi, blabel ) ;
         } // mbi.
      }

    //------ End: Compute average over Njets for double ratios


      h_ratio_htl -> GetXaxis() -> LabelsOption("v") ;
      h_ratio_htm -> GetXaxis() -> LabelsOption("v") ;
      h_ratio_hth -> GetXaxis() -> LabelsOption("v") ;

      h_dr_htl -> GetXaxis() -> LabelsOption("v") ;
      h_dr_htm -> GetXaxis() -> LabelsOption("v") ;
      h_dr_hth -> GetXaxis() -> LabelsOption("v") ;

      h_dr_htl_njave -> GetXaxis() -> LabelsOption("v") ;
      h_dr_htm_njave -> GetXaxis() -> LabelsOption("v") ;
      h_dr_hth_njave -> GetXaxis() -> LabelsOption("v") ;


      gStyle -> SetPadBottomMargin(0.30) ;

      char fname[10000] ;

     //---

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1_draw_qcd_ratio" ) ;
      if ( can1 == 0x0 ) can1 = new TCanvas( "can1_draw_qcd_ratio", "histograms", 900, 600 ) ;
      can1 -> Clear() ;
      can1 -> cd() ;
      h_ldp -> SetLineColor(4) ;
      h_hdp -> SetLineColor(2) ;

      h_ldp -> Draw( "" ) ;
      h_ldp -> Draw( "hist same" ) ;
      h_hdp -> Draw( "same" ) ;
      h_hdp -> Draw( "hist same" ) ;

      can1 -> Update() ; can1 -> Draw() ;

      gPad -> SetLogy(1) ;
      gPad -> SetGridx(0) ;
      gPad -> SetGridy(1) ;

      draw_boundaries() ;

      sprintf( fname, "%s/qcdr-plot-histograms-%s.pdf", outputdir, plottype ) ;
      can1 -> SaveAs( fname ) ;

     //---

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2_draw_qcd_ratio" ) ;
      if ( can2 == 0x0 ) can2 = new TCanvas( "can2_draw_qcd_ratio", "high/low ratio", 900, 600 ) ;
      can2 -> cd() ;

      h_ratio -> SetMaximum(ratio_max) ;
      h_ratio -> SetMinimum(-0.05) ;
      h_ratio -> SetMarkerStyle( 20 ) ;

      h_ratio -> Draw() ;
      rline -> DrawLine( h_ratio -> GetXaxis() -> GetXmin(), 0.,  h_ratio -> GetXaxis() -> GetXmax(), 0. ) ;
      h_ratio -> Draw("same") ;
      can2 -> Update() ; can2 -> Draw() ;


      gPad -> SetGridx(0) ;
      gPad -> SetGridy(1) ;

      draw_boundaries() ;
      line0 -> DrawLine( h_ratio->GetXaxis()->GetXmin(), 0., h_ratio->GetXaxis()->GetXmax(), 0. );

      sprintf( fname, "%s/qcdr-plot-ratio-%s.pdf", outputdir, plottype ) ;
      can2 -> SaveAs( fname ) ;

      h_ratio -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-ratiozoom-%s.pdf", outputdir, plottype ) ;
      can2 -> SaveAs( fname ) ;

     //---

      TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3_draw_qcd_ratio" ) ;
      if ( can3 == 0x0 ) can3 = new TCanvas( "can3_draw_qcd_ratio", "double high/low ratio", 900, 600 ) ;
      can3 -> cd() ;

      h_double_ratio -> SetMaximum(3.00) ;
      h_double_ratio -> SetMinimum(-0.5) ;
      h_double_ratio -> SetMarkerStyle( 20 ) ;

      h_double_ratio -> Draw() ;
      rline -> DrawLine( h_double_ratio -> GetXaxis() -> GetXmin(), 0.,  h_double_ratio -> GetXaxis() -> GetXmax(), 0. ) ;
      h_double_ratio -> Draw("same") ;
      can3 -> Update() ; can3 -> Draw() ;

      gPad -> SetGridx(0) ;
      gPad -> SetGridy(1) ;

      draw_boundaries() ;

      sprintf( fname, "%s/qcdr-plot-dratio-%s.pdf", outputdir, plottype ) ;
      can3 -> SaveAs( fname ) ;

     //---

      h_control_ratio_mht0 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_mht1 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_mht2 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_mht3 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio_mht4 -> GetXaxis() -> LabelsOption( "v" ) ;

      h_dr_mht1 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_dr_mht2 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_dr_mht3 -> GetXaxis() -> LabelsOption( "v" ) ;
      h_dr_mht4 -> GetXaxis() -> LabelsOption( "v" ) ;

      h_dr_mht1 -> SetMarkerStyle(20) ;
      h_dr_mht2 -> SetMarkerStyle(20) ;
      h_dr_mht3 -> SetMarkerStyle(20) ;
      h_dr_mht4 -> SetMarkerStyle(20) ;

      h_dr_mht1 -> SetMinimum(-0.5) ; h_dr_mht1 -> SetMaximum( 3.0) ;
      h_dr_mht2 -> SetMinimum(-0.5) ; h_dr_mht2 -> SetMaximum( 3.0) ;
      h_dr_mht3 -> SetMinimum(-0.5) ; h_dr_mht3 -> SetMaximum( 3.0) ;
      h_dr_mht4 -> SetMinimum(-0.5) ; h_dr_mht4 -> SetMaximum( 3.0) ;

      h_dr_htl -> SetMarkerStyle(20) ;
      h_dr_htm -> SetMarkerStyle(20) ;
      h_dr_hth -> SetMarkerStyle(20) ;

      h_dr_htl -> SetMinimum(-0.5) ; h_dr_htl -> SetMaximum( drmax ) ;
      h_dr_htm -> SetMinimum(-0.5) ; h_dr_htm -> SetMaximum( drmax ) ;
      h_dr_hth -> SetMinimum(-0.5) ; h_dr_hth -> SetMaximum( drmax ) ;

      h_dr_htl_njave -> SetMarkerStyle(20) ;
      h_dr_htm_njave -> SetMarkerStyle(20) ;
      h_dr_hth_njave -> SetMarkerStyle(20) ;

      h_dr_htl_njave -> SetMinimum(-0.2) ; h_dr_htl_njave -> SetMaximum( 2.0) ;
      h_dr_htm_njave -> SetMinimum(-0.2) ; h_dr_htm_njave -> SetMaximum( 2.0) ;
      h_dr_hth_njave -> SetMinimum(-0.2) ; h_dr_hth_njave -> SetMaximum( 2.0) ;






     //---

      TCanvas* can4 = (TCanvas*) gDirectory -> FindObject( "can4_draw_qcd_ratio" ) ;
      if ( can4 == 0x0 ) can4 = new TCanvas( "can4_draw_qcd_ratio", "MHT1: double high/low ratio", 600, 600 ) ;
      can4 -> cd() ;


      h_dr_mht1 -> Draw() ;
      can4 -> Update() ; can4 -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-mht1-dratio-%s.pdf", outputdir, plottype ) ;
      can4 -> SaveAs( fname ) ;



     //---

      TCanvas* can5 = (TCanvas*) gDirectory -> FindObject( "can5_draw_qcd_ratio" ) ;
      if ( can5 == 0x0 ) can5 = new TCanvas( "can5_draw_qcd_ratio", "MHT2: double high/low ratio", 600, 600 ) ;
      can5 -> cd() ;

      h_dr_mht2 -> Draw() ;
      can5 -> Update() ; can5 -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-mht2-dratio-%s.pdf", outputdir, plottype ) ;
      can5 -> SaveAs( fname ) ;



     //---

      TCanvas* can6 = (TCanvas*) gDirectory -> FindObject( "can6_draw_qcd_ratio" ) ;
      if ( can6 == 0x0 ) can6 = new TCanvas( "can6_draw_qcd_ratio", "MHT3: double high/low ratio", 600, 600 ) ;
      can6 -> cd() ;

      h_dr_mht3 -> Draw() ;
      can6 -> Update() ; can6 -> Draw() ;

      draw_boundaries_nj( 2 ) ;
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-mht3-dratio-%s.pdf", outputdir, plottype ) ;
      can6 -> SaveAs( fname ) ;



     //---

      TCanvas* can7 = (TCanvas*) gDirectory -> FindObject( "can7_draw_qcd_ratio" ) ;
      if ( can7 == 0x0 ) can7 = new TCanvas( "can7_draw_qcd_ratio", "MHT4: double high/low ratio", 600, 600 ) ;
      can7 -> cd() ;

      h_dr_mht4 -> Draw() ;
      can7 -> Update() ; can7 -> Draw() ;

      draw_boundaries_nj( 2 ) ;
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-mht4-dratio-%s.pdf", outputdir, plottype ) ;
      can7 -> SaveAs( fname ) ;




     //---

      TCanvas* can8 = (TCanvas*) gDirectory -> FindObject( "can8_draw_qcd_ratio" ) ;
      if ( can8 == 0x0 ) can8 = new TCanvas( "can8_draw_qcd_ratio", "MHT0: control high/low ratio", 600, 600 ) ;
      can8 -> cd() ;

      h_control_ratio_mht0 -> SetMarkerStyle(20) ;
      h_control_ratio_mht0 -> SetMinimum( -0.02 ) ;
      h_control_ratio_mht0 -> SetMaximum(  ratio_max ) ;

      h_control_ratio_mht0 -> Draw() ;
      can8 -> Update() ; can8 -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-mht0-ratio-%s.pdf", outputdir, plottype ) ;
      can8 -> SaveAs( fname ) ;

      h_control_ratio_mht0 -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-mht0-ratiozoom-%s.pdf", outputdir, plottype ) ;
      can8 -> SaveAs( fname ) ;

     //---

      TCanvas* can9 = (TCanvas*) gDirectory -> FindObject( "can9_draw_qcd_ratio" ) ;
      if ( can9 == 0x0 ) can9 = new TCanvas( "can9_draw_qcd_ratio", "MHT1: high/low ratio", 600, 600 ) ;
      can9 -> cd() ;

      h_ratio_mht1 -> SetMarkerStyle(20) ;
      h_ratio_mht1 -> SetMinimum( -0.02 ) ;
      h_ratio_mht1 -> SetMaximum(  ratio_max ) ;

      h_ratio_mht1 -> Draw() ;
      can9 -> Update() ; can9 -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      gPad -> SetGridy(1) ;


      sprintf( fname, "%s/qcdr-plot-mht1-ratio-%s.pdf", outputdir, plottype ) ;
      can9 -> SaveAs( fname ) ;

      h_ratio_mht1 -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-mht1-ratiozoom-%s.pdf", outputdir, plottype ) ;
      can9 -> SaveAs( fname ) ;


     //---

      TCanvas* can10 = (TCanvas*) gDirectory -> FindObject( "can10_draw_qcd_ratio" ) ;
      if ( can10 == 0x0 ) can10 = new TCanvas( "can10_draw_qcd_ratio", "MHT2: high/low ratio", 600, 600 ) ;
      can10 -> cd() ;

      h_ratio_mht2 -> SetMarkerStyle(20) ;
      h_ratio_mht2 -> SetMinimum( -0.02 ) ;
      h_ratio_mht2 -> SetMaximum(  ratio_max ) ;

      h_ratio_mht2 -> Draw() ;
      can10 -> Update() ; can10 -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-mht2-ratio-%s.pdf", outputdir, plottype ) ;
      can10 -> SaveAs( fname ) ;

      h_ratio_mht2 -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-mht2-ratiozoom-%s.pdf", outputdir, plottype ) ;
      can10 -> SaveAs( fname ) ;


     //---

      TCanvas* can11 = (TCanvas*) gDirectory -> FindObject( "can11_draw_qcd_ratio" ) ;
      if ( can11 == 0x0 ) can11 = new TCanvas( "can11_draw_qcd_ratio", "MHT3: high/low ratio", 600, 600 ) ;
      can11 -> cd() ;

      h_ratio_mht3 -> SetMarkerStyle(20) ;
      h_ratio_mht3 -> SetMinimum( -0.02 ) ;
      h_ratio_mht3 -> SetMaximum(  ratio_max ) ;

      h_ratio_mht3 -> Draw() ;
      can11 -> Update() ; can11 -> Draw() ;

      draw_boundaries_nj( 2 ) ;
      gPad -> SetGridy(1) ;


      sprintf( fname, "%s/qcdr-plot-mht3-ratio-%s.pdf", outputdir, plottype ) ;
      can11 -> SaveAs( fname ) ;

      h_ratio_mht3 -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-mht3-ratiozoom-%s.pdf", outputdir, plottype ) ;
      can11 -> SaveAs( fname ) ;

     //---

      TCanvas* can12 = (TCanvas*) gDirectory -> FindObject( "can12_draw_qcd_ratio" ) ;
      if ( can12 == 0x0 ) can12 = new TCanvas( "can12_draw_qcd_ratio", "MHT4: high/low ratio", 600, 600 ) ;
      can12 -> cd() ;

      h_ratio_mht4 -> SetMarkerStyle(20) ;
      h_ratio_mht4 -> SetMinimum( -0.02 ) ;
      h_ratio_mht4 -> SetMaximum(  ratio_max ) ;

      h_ratio_mht4 -> Draw() ;
      can12 -> Update() ; can12 -> Draw() ;

      draw_boundaries_nj( 2 ) ;
      gPad -> SetGridy(1) ;


      sprintf( fname, "%s/qcdr-plot-mht4-ratio-%s.pdf", outputdir, plottype ) ;
      can12 -> SaveAs( fname ) ;

      h_ratio_mht4 -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-mht4-ratiozoom-%s.pdf", outputdir, plottype ) ;
      can12 -> SaveAs( fname ) ;

     //---

      TCanvas* canht = (TCanvas*) gDirectory -> FindObject( "canht_draw_qcd_ratio" ) ;
      if ( canht == 0x0 ) canht = new TCanvas( "canht_draw_qcd_ratio", "MHT4: high/low ratio", 600, 600 ) ;
      canht -> cd() ;

      h_ratio_htl -> SetMarkerStyle(20) ;
      h_ratio_htl -> SetMinimum( -0.02 ) ;
      h_ratio_htl -> SetMaximum(  ratio_max ) ;

      h_ratio_htl -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      line0 -> DrawLine( h_ratio_htl->GetXaxis()->GetXmin(), 0., h_ratio_htl->GetXaxis()->GetXmax(), 0. );
      gPad -> SetGridy(1) ;


      sprintf( fname, "%s/qcdr-plot-htl-ratio-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;

      h_ratio_htl -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-htl-ratiozoom-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;

     //---

      h_ratio_htm -> SetMarkerStyle(20) ;
      h_ratio_htm -> SetMinimum( -0.02 ) ;
      h_ratio_htm -> SetMaximum(  ratio_max ) ;

      h_ratio_htm -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      draw_boundaries_nj( 5 ) ;
      line0 -> DrawLine( h_ratio_htm->GetXaxis()->GetXmin(), 0., h_ratio_htm->GetXaxis()->GetXmax(), 0. );
      gPad -> SetGridy(1) ;


      sprintf( fname, "%s/qcdr-plot-htm-ratio-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;

      h_ratio_htm -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-htm-ratiozoom-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;

     //---

      h_ratio_hth -> SetMarkerStyle(20) ;
      h_ratio_hth -> SetMinimum( -0.02 ) ;
      h_ratio_hth -> SetMaximum(  ratio_max ) ;

      h_ratio_hth -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      draw_boundaries_nj( 5 ) ;
      line0 -> DrawLine( h_ratio_hth->GetXaxis()->GetXmin(), 0., h_ratio_hth->GetXaxis()->GetXmax(), 0. );
      gPad -> SetGridy(1) ;


      sprintf( fname, "%s/qcdr-plot-hth-ratio-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;

      h_ratio_hth -> SetMaximum( rzoom1max ) ;
      sprintf( fname, "%s/qcdr-plot-hth-ratiozoom-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;

     //---

      h_dr_htl -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      draw_boundaries_nj( 3 ) ;
      line0 -> DrawLine( h_dr_htl->GetXaxis()->GetXmin(), 0., h_dr_htl->GetXaxis()->GetXmax(), 0. );
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-htl-dratio-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;


     //---

      h_dr_htm -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      draw_boundaries_nj( 5 ) ;
      line0 -> DrawLine( h_dr_htm->GetXaxis()->GetXmin(), 0., h_dr_htm->GetXaxis()->GetXmax(), 0. );
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-htm-dratio-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;


     //---

      h_dr_hth -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      draw_boundaries_nj( 5 ) ;
      line0 -> DrawLine( h_dr_hth->GetXaxis()->GetXmin(), 0., h_dr_hth->GetXaxis()->GetXmax(), 0. );
      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-hth-dratio-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;


     //---


      h_dr_htl_njave -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-htl-drationjave-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;


     //---

      h_dr_htm_njave -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-htm-drationjave-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;


     //---

      h_dr_hth_njave -> Draw() ;
      canht -> Update() ; canht -> Draw() ;

      gPad -> SetGridy(1) ;

      sprintf( fname, "%s/qcdr-plot-hth-drationjave-%s.pdf", outputdir, plottype ) ;
      canht -> SaveAs( fname ) ;


     //---






      sprintf( fname, "outputfiles/qcd-ratio-%s.root", plottype ) ;
      saveHist( fname, "h*" ) ;

///   printf("\n\n") ;
///   printf("=================================================================================================================================\n") ;
///   for ( int bi=1; bi<=h_ldp->GetNbinsX(); bi++ ) {
///      printf("  %25s : LDP = %9.1f +/- %5.1f ,   HDP = %9.1f +/- %5.1f ,   H/L = %6.3f +/- %6.3f ,   DR = %5.2f +/- %5.2f\n",
///         h_ratio->GetXaxis()->GetBinLabel(bi),
///         h_ldp -> GetBinContent(bi), h_ldp -> GetBinError(bi),
///         h_hdp -> GetBinContent(bi), h_hdp -> GetBinError(bi),
///         h_ratio -> GetBinContent(bi), h_ratio -> GetBinError(bi),
///         h_double_ratio -> GetBinContent(bi), h_double_ratio -> GetBinError(bi)
///         ) ;
///      if ( bi%13 == 0 ) { printf("=============================================================================================================================================\n") ; }
///      if ( (bi-3)%13 == 0 || bi==3 ) { printf("---------------------------------------------------------------------------------------------------------------------------------------------\n") ; }
///   } // bi
///   printf("\n\n") ;

   } // draw_qcd_ratio_v5

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





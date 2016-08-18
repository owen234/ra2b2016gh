#include "TDirectory.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TText.h"
#include "TLine.h"
#include <fstream>
#include "binning.h"
#include <iostream>
#include <string>
#include <sstream>

bool transfer_qcd_parameters(string filein_name, string fileout_name);

void draw_badjet_cat_v3(const char* infile = "outputfiles/syst-2015-v2.root" ) {

      char fname[10000] ;
      char htstr[10];
      gStyle -> SetOptStat(0) ;
      setup_bins();

      TFile* tf = new TFile( infile, "READ" ) ;
      if ( ! tf -> IsOpen() ) { printf( "\n\n *** Bad input file: %s\n\n", infile ) ; return ; }

      FILE * out_file;

      if ( transfer_qcd_parameters("./outputfiles/qcdmc-chi2-fit-model-pars.txt","./outputfiles/model-pars-qcdmc3.txt") )
      {
         out_file = fopen ("./outputfiles/model-pars-qcdmc3.txt","a");// "a" for append
      }
      else 
         out_file = fopen ("./outputfiles/model-pars-qcdmc3.txt","w");// "w" for write
      



      gDirectory -> Delete( "h*" ) ;

      TCanvas* can = new TCanvas( "can_draw_badjet_cat", "", 800, 600 ) ;
      TCanvas* can2 = new TCanvas( "can2_draw_badjet_cat", "", 500, 600 ) ;

      int nbins(5) ;

      float histshift(0.03) ;

      //double fmiss_val = 0.022 ;
      //double fmiss_err = 0.004 ;

      double fmiss_val[6] ;
      double fmiss_err[6] ;

   // fmiss_val[1] = 0.036 ;  fmiss_err[1] = 0.001 ;
   // fmiss_val[2] = 0.028 ;  fmiss_err[2] = 0.002 ;
   // fmiss_val[3] = 0.018 ;  fmiss_err[3] = 0.001 ;
   // fmiss_val[4] = 0.012 ;  fmiss_err[4] = 0.002 ;
   // fmiss_val[5] = 0.011 ;  fmiss_err[5] = 0.002 ;
   

   for ( int ht_level = nBinsHT-1; ht_level >= 0; ht_level--)
   {

      if ( ht_level == 0 ) strcpy (htstr, "htl");
      if ( ht_level == 1 ) strcpy (htstr, "htm");
      if ( ht_level == 2 ) strcpy (htstr, "hth");

      TString ht_str = htstr;

      TH1F* h_rhl_bj_in_dphi = new TH1F( "h_rhl_bj_in_dphi_"+ht_str, "Rhigh/low, bad jet in dphi", nbins, 0.5-histshift, nbins+0.5-histshift ) ;
      TH1F* h_fmissfpass_bj_in_dphi = new TH1F( "h_fmissfpass_bj_in_dphi_"+ht_str, "Fr(bj not in dphi)*fr(pass)", nbins, 0.5+histshift, nbins+0.5+histshift ) ;
      TH1F* h_rprime = new TH1F( "h_rprime_"+ht_str, "Effective Rhigh/low", nbins, 0.5, nbins+0.5 ) ;

 
      nbins = 5 ;
      if ( strcmp( htstr, "htl" ) == 0 ) { nbins = 3 ; }


      for ( int mbi=1; mbi<=nbins; mbi++ ) {

         char hname[100] ;

         if ( mbi==1 ) { sprintf( hname, "h_mdp_%s_mhtc", htstr ) ; } else { sprintf( hname, "h_mdp_%s_mht%d", htstr, mbi-1 ) ; }
         TH1F* hp_all = (TH1F*) tf -> Get( hname ) ;
         if ( hp_all == 0x0 ) { printf("\n\n *** missing %s\n\n", hname ) ; return ; }

         if ( mbi==1 ) { sprintf( hname, "h_mdp_%s_badj_in_dphi_mhtc", htstr ) ; } else { sprintf( hname, "h_mdp_%s_badj_in_dphi_mht%d", htstr, mbi-1 ) ; }
         TH1F* hp_bj_in_dphi = (TH1F*) tf -> Get( hname ) ;
         if ( hp_bj_in_dphi == 0x0 ) { printf("\n\n *** missing %s\n\n", hname ) ; return ; }

         if ( mbi==1 ) { sprintf( hname, "h_mdp_%s_badj_not_in_dphi_mhtc", htstr ) ; } else { sprintf( hname, "h_mdp_%s_badj_not_in_dphi_mht%d", htstr, mbi-1 ) ; }
         TH1F* hp_bj_not_in_dphi = (TH1F*) tf -> Get( hname ) ;
         if ( hp_bj_not_in_dphi == 0x0 ) { printf("\n\n *** missing %s\n\n", hname ) ; return ; }




         if ( mbi==1 ) { sprintf( hname, "h_dphiregion_%s_mhtc", htstr ) ; } else { sprintf( hname, "h_dphiregion_%s_mht%d", htstr, mbi-1 ) ; }
         TH1F* hp_dphiregion_all = (TH1F*) tf -> Get( hname ) ;
         if ( hp_dphiregion_all == 0x0 ) { printf("\n\n *** missing%s\n\n", hname ) ; return ; }

         if ( mbi==1 ) { sprintf( hname, "h_dphiregion_%s_badj_in_dphi_mhtc", htstr ) ; } else { sprintf( hname, "h_dphiregion_%s_badj_in_dphi_mht%d", htstr, mbi-1 ) ; }
         TH1F* hp_dphiregion_bj_in_dphi = (TH1F*) tf -> Get( hname ) ;
         if ( hp_dphiregion_bj_in_dphi == 0x0 ) { printf("\n\n *** missing%s\n\n", hname ) ; return ; }

         if ( mbi==1 ) { sprintf( hname, "h_dphiregion_%s_badj_not_in_dphi_mhtc", htstr ) ; } else { sprintf( hname, "h_dphiregion_%s_badj_not_in_dphi_mht%d", htstr, mbi-1 ) ; }
         TH1F* hp_dphiregion_bj_not_in_dphi = (TH1F*) tf -> Get( hname ) ;
         if ( hp_dphiregion_bj_not_in_dphi == 0x0 ) { printf("\n\n *** missing%s\n\n", hname ) ; return ; }







         double all_val(0.), all_err(0.) ;
         double bj_in_dphi_val(0.), bj_in_dphi_err(0.) ;
         double bj_not_in_dphi_val(0.), bj_not_in_dphi_err(0.) ;

         double all_pass_val(0.), all_pass_err(0.) ;
         double all_fail_val(0.), all_fail_err(0.) ;
         double bj_in_dphi_pass_val(0.), bj_in_dphi_pass_err(0.) ;
         double bj_in_dphi_fail_val(0.), bj_in_dphi_fail_err(0.) ;
         double bj_not_in_dphi_pass_val(0.), bj_not_in_dphi_pass_err(0.) ;
         double bj_not_in_dphi_fail_val(0.), bj_not_in_dphi_fail_err(0.) ;


         all_val = hp_all -> IntegralAndError( 1, hp_all->GetNbinsX(), all_err ) ;
         bj_in_dphi_val = hp_bj_in_dphi -> IntegralAndError( 1, hp_bj_in_dphi->GetNbinsX(), bj_in_dphi_err ) ;
         bj_not_in_dphi_val = hp_bj_not_in_dphi -> IntegralAndError( 1, hp_bj_not_in_dphi->GetNbinsX(), bj_not_in_dphi_err ) ;

         all_pass_val = hp_dphiregion_all -> GetBinContent( 2 ) ;
         all_pass_err = hp_dphiregion_all -> GetBinError( 2 ) ;

         all_fail_val = hp_dphiregion_all -> GetBinContent( 1 ) ;
         all_fail_err = hp_dphiregion_all -> GetBinError( 1 ) ;

         bj_in_dphi_pass_val = hp_dphiregion_bj_in_dphi -> GetBinContent( 2 ) ;
         bj_in_dphi_pass_err = hp_dphiregion_bj_in_dphi -> GetBinError( 2 ) ;

         bj_in_dphi_fail_val = hp_dphiregion_bj_in_dphi -> GetBinContent( 1 ) ;
         bj_in_dphi_fail_err = hp_dphiregion_bj_in_dphi -> GetBinError( 1 ) ;

         bj_not_in_dphi_pass_val = hp_dphiregion_bj_not_in_dphi -> GetBinContent( 2 ) ;
         bj_not_in_dphi_pass_err = hp_dphiregion_bj_not_in_dphi -> GetBinError( 2 ) ;

         bj_not_in_dphi_fail_val = hp_dphiregion_bj_not_in_dphi -> GetBinContent( 1 ) ;
         bj_not_in_dphi_fail_err = hp_dphiregion_bj_not_in_dphi -> GetBinError( 1 ) ;



         double pass_frac_bj_not_in_dphi_val(0.), pass_frac_bj_not_in_dphi_err(0.) ;
         pass_frac_bj_not_in_dphi_val = bj_not_in_dphi_pass_val / ( bj_not_in_dphi_pass_val + bj_not_in_dphi_fail_val ) ;
         if ( bj_not_in_dphi_pass_val>0 && bj_not_in_dphi_fail_val>0 ) {
            pass_frac_bj_not_in_dphi_err = sqrt(
               pow( bj_not_in_dphi_pass_val * bj_not_in_dphi_fail_err , 2. ) /
                  pow( bj_not_in_dphi_pass_val + bj_not_in_dphi_fail_val, 4. )
             + pow( bj_not_in_dphi_fail_val * bj_not_in_dphi_pass_err , 2. ) /
                  pow( bj_not_in_dphi_pass_val + bj_not_in_dphi_fail_val, 4. )
            ) ;
         }

         printf("\n\n") ;
         printf("  MHT bin %d\n", mbi ) ;
         printf(" all = %.2f +/- %.2f,  pass = %.2f +/- %.2f,  fail = %.2f +/- %.2f\n",
            all_val, all_err,  all_pass_val, all_pass_err,  all_fail_val, all_fail_err ) ;
         printf(" bj_in_dphi = %.2f +/- %.2f,  pass = %.2f +/- %.2f,  fail = %.2f +/- %.2f\n",
            bj_in_dphi_val, bj_in_dphi_err,  bj_in_dphi_pass_val, bj_in_dphi_pass_err,  bj_in_dphi_fail_val, bj_in_dphi_fail_err ) ;
         printf(" bj_not_in_dphi = %.2f +/- %.2f,  pass = %.2f +/- %.2f,  fail = %.2f +/- %.2f\n",
            bj_not_in_dphi_val, bj_not_in_dphi_err,  bj_not_in_dphi_pass_val, bj_not_in_dphi_pass_err,  bj_not_in_dphi_fail_val, bj_not_in_dphi_fail_err ) ;
         printf("  pass fraction for bad jet not in dphi:  %5.3f +/- %5.3f\n",
           pass_frac_bj_not_in_dphi_val, pass_frac_bj_not_in_dphi_err ) ;


         bj_not_in_dphi_val = hp_bj_not_in_dphi -> IntegralAndError( 1, hp_bj_not_in_dphi->GetNbinsX(), bj_not_in_dphi_err ) ;

         double frac_bj_not_in_dphi_val(0.), frac_bj_not_in_dphi_err(0.) ;
         if ( all_val > 0 ) {
            frac_bj_not_in_dphi_val = bj_not_in_dphi_val / all_val ;
            frac_bj_not_in_dphi_err = bj_not_in_dphi_err / all_val ;
         }
         fmiss_val[mbi] = frac_bj_not_in_dphi_val ; // *** take it from here now instead of setting by hand at beginning.
         fmiss_err[mbi] = frac_bj_not_in_dphi_err ; // *** take it from here now instead of setting by hand at beginning.

         double r_hl_all_val(0.), r_hl_all_err(0.) ;
         if ( all_fail_val > 0. && all_pass_val > 0 ) {
            r_hl_all_val = all_pass_val / all_fail_val ;
            r_hl_all_err = r_hl_all_val * sqrt( pow(all_pass_err/all_pass_val,2.) + pow(all_fail_err/all_fail_val,2.) ) ;
         }

         double r_hl_bj_in_dphi_val(0.), r_hl_bj_in_dphi_err(0.) ;
         if ( bj_in_dphi_fail_val > 0. && bj_in_dphi_pass_val > 0 ) {
            r_hl_bj_in_dphi_val = bj_in_dphi_pass_val / bj_in_dphi_fail_val ;
            r_hl_bj_in_dphi_err = r_hl_bj_in_dphi_val * sqrt( pow(bj_in_dphi_pass_err/bj_in_dphi_pass_val,2.) + pow(bj_in_dphi_fail_err/bj_in_dphi_fail_val,2.) ) ;
         }

         double r_hl_bj_not_in_dphi_val(0.), r_hl_bj_not_in_dphi_err(0.) ;
         if ( bj_not_in_dphi_fail_val > 0. && bj_not_in_dphi_pass_val > 0 ) {
            r_hl_bj_not_in_dphi_val = bj_not_in_dphi_pass_val / bj_not_in_dphi_fail_val ;
            r_hl_bj_not_in_dphi_err = r_hl_bj_not_in_dphi_val * sqrt( pow(bj_not_in_dphi_pass_err/bj_not_in_dphi_pass_val,2.) + pow(bj_not_in_dphi_fail_err/bj_not_in_dphi_fail_val,2.) ) ;
         }

         h_rhl_bj_in_dphi -> SetBinContent( mbi, r_hl_bj_in_dphi_val ) ;
         h_rhl_bj_in_dphi -> SetBinError( mbi, r_hl_bj_in_dphi_err ) ;

         double fmfp_val = fmiss_val[mbi] * pass_frac_bj_not_in_dphi_val ;
         double fmfp_err = fmfp_val * sqrt( pow( fmiss_err[mbi]/fmiss_val[mbi], 2 ) + pow( pass_frac_bj_not_in_dphi_err/pass_frac_bj_not_in_dphi_val, 2 ) ) ;

         h_fmissfpass_bj_in_dphi -> SetBinContent( mbi, fmfp_val ) ;
         h_fmissfpass_bj_in_dphi -> SetBinError( mbi, fmfp_err ) ;

         double rprime_val = (1.-fmiss_val[mbi])*r_hl_bj_in_dphi_val + fmfp_val ;
         double rprime_err = sqrt(
            pow( (1.-fmiss_val[mbi])*r_hl_bj_in_dphi_err, 2 )
            + pow( (pass_frac_bj_not_in_dphi_val - r_hl_bj_in_dphi_err )*fmiss_err[mbi], 2. )
            + pow( (fmiss_val[mbi]*pass_frac_bj_not_in_dphi_err)  , 2. )
         ) ;

         printf("  MHT bin %d : all = %8.2f +/- %6.2f,   bad jet not in dphi = %8.2f +/- %6.2f,  fraction = %5.3f +/- %5.3f\n",
            mbi, all_val, all_err, bj_not_in_dphi_val, bj_not_in_dphi_err, frac_bj_not_in_dphi_val, frac_bj_not_in_dphi_err ) ;

         printf("  Rhigh/low  :  All = %5.3f +/- %5.3f,   bad jet in dphi = %5.3f +/- %5.3f,   bad jet not in dphi = %5.3f +/- %5.3f\n",
            r_hl_all_val, r_hl_all_err,
            r_hl_bj_in_dphi_val, r_hl_bj_in_dphi_err,
            r_hl_bj_not_in_dphi_val, r_hl_bj_not_in_dphi_err
            ) ;
         printf("  R' (effective high/low ratio) : %5.3f +/- %5.3f\n", rprime_val, rprime_err ) ;

         h_rprime -> SetBinContent( mbi, rprime_val ) ;
         h_rprime -> SetBinError( mbi, rprime_err ) ;

         char binlabel[100] ;
         if ( mbi == 1 ) {
            sprintf( binlabel, "MHTC" ) ;
         } else {
            sprintf( binlabel, "MHT%d", mbi-1 ) ;
         }
         h_rprime -> GetXaxis() -> SetBinLabel( mbi, binlabel ) ;
         h_rprime -> SetLabelSize( 0.07, "x" ) ;



         hp_all -> SetFillColor( kRed - 9 ) ;
         hp_bj_not_in_dphi -> SetFillColor( kBlue - 9 ) ;

         hp_all -> SetMaximum( 1.10 * ( hp_all->GetMaximum()) ) ;
         hp_all -> SetMinimum( -0.05 * ( hp_all->GetMaximum()) ) ;

         sprintf( hname, "%s_c", hp_all -> GetName() ) ;
         TH1F* hp_all_c = (TH1F*) hp_all -> Clone( hname ) ;
         sprintf( hname, "%s_c", hp_bj_not_in_dphi -> GetName() ) ;
         TH1F* hp_bj_not_in_dphi_c = (TH1F*) hp_bj_not_in_dphi -> Clone( hname ) ;

         hp_all -> SetXTitle( "min Delta phi" ) ;
         hp_all -> SetTitleSize( 0.05, "x" ) ;
         hp_all -> SetYTitle( "QCD Events / 0.05" ) ;
         hp_all -> SetTitleSize( 0.05, "y" ) ;

         char ht_title[100] ;
         if ( strcmp( htstr, "hth" ) == 0 ) { sprintf( ht_title, "High HT" ) ; }
         if ( strcmp( htstr, "htm" ) == 0 ) { sprintf( ht_title, "Medium HT" ) ; }
         if ( strcmp( htstr, "htl" ) == 0 ) { sprintf( ht_title, "Low HT" ) ; }

         TText* title_text = new TText() ;
         title_text -> SetTextSize( 0.050 ) ;
         char title[1000] ;
         if ( mbi == 1 ) sprintf( title, "MHT bin C [250,300], %s", ht_title ) ;
         if ( mbi == 2 ) sprintf( title, "MHT bin 1 [300,350], %s", ht_title ) ;
         if ( mbi == 3 ) sprintf( title, "MHT bin 2 [350,500], %s", ht_title ) ;
         if ( mbi == 4 ) sprintf( title, "MHT bin 3 [500,750], %s", ht_title ) ;
         if ( mbi == 5 ) sprintf( title, "MHT bin 4 [>750], %s", ht_title ) ;



         can -> cd() ;
         hp_all -> Draw() ;
         hp_all -> Draw("hist same") ;
         hp_bj_not_in_dphi -> Draw("hist same") ;
         hp_bj_not_in_dphi -> Draw("same") ;
         hp_all -> Draw("same" ) ;
         hp_all -> Draw("same axis" ) ;
         hp_all -> Draw("same axig" ) ;
         TLine* l = new TLine() ;
         l -> SetLineStyle(3) ;
         l -> DrawLine( 0.5, 0, 0.5, hp_all -> GetMaximum() ) ;
         l -> DrawLine( 0.3, 0, 0.3, hp_all -> GetMaximum() ) ;
         title_text -> DrawTextNDC( 0.05, 0.92, title ) ;
         can -> Update() ;

         char pname[100] ;
         sprintf( pname, "inset_pad_%d", mbi ) ;
         gStyle -> SetPadTopMargin(0) ;
         gStyle -> SetPadRightMargin(0) ;
         gStyle -> SetOptTitle(0) ;
         TPad* inset_pad = new TPad( pname, "", 0.25, 0.25, 0.90, 0.90 ) ;
         hp_all_c -> SetMaximum( 0.1 * ( hp_all_c->GetMaximum()) ) ;
         hp_all_c -> SetMinimum( -0.05 * ( hp_all_c->GetMaximum()) ) ;
         inset_pad -> Draw() ;
         inset_pad -> cd() ;
         hp_all_c -> Draw() ;
         hp_all_c -> Draw("hist same") ;
         hp_bj_not_in_dphi_c -> Draw("hist same") ;
         hp_bj_not_in_dphi_c -> Draw("same") ;
         hp_all_c -> Draw("same" ) ;
         hp_all_c -> Draw("same axis" ) ;
         hp_all_c -> Draw("same axig" ) ;
         l -> DrawLine( 0.5, 0, 0.5, hp_all_c -> GetMaximum() ) ;
         l -> DrawLine( 0.3, 0, 0.3, hp_all_c -> GetMaximum() ) ;

         TLegend* legend = new TLegend( 0.20, 0.75, 0.95, 0.95 ) ;
         legend -> AddEntry( hp_all_c, "bad jet in min Delta phi" ) ;
         legend -> AddEntry( hp_bj_not_in_dphi_c, "bad jet not in min Delta phi" ) ;
         legend -> Draw() ;

         TText* text = new TText() ;
         text -> SetTextSize( 0.060 ) ;
         float tx, ty, tdy ;
         char textcs[1000] ;
         tx = 0.28 ;
         ty = 0.65 ;
         tdy = 0.09 ;
         sprintf( textcs, "Frac. bad jet not in Delta phi" ) ;
         text -> DrawTextNDC( tx, ty, textcs ) ;

         ty -= tdy ;
         sprintf( textcs, "        %5.3f +/- %5.3f", frac_bj_not_in_dphi_val, frac_bj_not_in_dphi_err ) ;
         text -> DrawTextNDC( tx, ty, textcs ) ;

   //    sprintf( textcs, "Bad jet in min Delta phi, Rhigh/low: %5.3f +/- %5.3f", r_hl_bj_in_dphi_val, r_hl_bj_in_dphi_err  ) ;
   //    ty -= tdy ;
   //    text -> DrawTextNDC( tx, ty, textcs ) ;

   //    sprintf( textcs, "Bad jet not in min Delta phi, Rhigh/low: %5.3f +/- %5.3f", r_hl_bj_not_in_dphi_val, r_hl_bj_not_in_dphi_err  ) ;
   //    ty -= tdy ;
   //    text -> DrawTextNDC( tx, ty, textcs ) ;

   //    sprintf( textcs, "Bad jet not in min Delta phi, Pass cut frac.: %5.3f +/- %5.3f", pass_frac_bj_not_in_dphi_val, pass_frac_bj_not_in_dphi_err  ) ;
   //    ty -= tdy ;
   //    text -> DrawTextNDC( tx, ty, textcs ) ;

         text -> SetTextSize( 0.065 ) ;
         sprintf( textcs, "R high/low: %5.3f +/- %5.3f", rprime_val, rprime_err  ) ;
         ty -= tdy ;
         ty -= tdy ;
         tx += 0.05 ;
         text -> DrawTextNDC( tx, ty, textcs ) ;



         sprintf( fname, "outputfiles/badjet-cat-mdp-%s-mht%d-liny.pdf", htstr, mbi ) ;
         can -> Update() ; can -> Draw() ;
         can -> SaveAs( fname ) ;

      } // mbi


    //------------


      gStyle -> SetPadTopMargin(0.09) ;
      gStyle -> SetPadBottomMargin(0.10) ;
      gStyle -> SetPadRightMargin(0.05) ;
      gStyle -> SetPadLeftMargin(0.20) ;
      gStyle -> SetOptTitle(0) ;

      TLegend* l2 = new TLegend( 0.50, 0.80, 0.95, 0.95 ) ;
      l2 -> AddEntry( h_rprime, "Overall" ) ;
      l2 -> AddEntry( h_rhl_bj_in_dphi, "Bad jet in Delta phi" ) ;
      l2 -> AddEntry( h_fmissfpass_bj_in_dphi, "Bad jet not in Delta phi" ) ;

      TH1F* h_rprime_c = (TH1F*) h_rprime -> Clone( "h_rprime_c" ) ;
      TH1F* h_fmissfpass_bj_in_dphi_c = (TH1F*) h_fmissfpass_bj_in_dphi -> Clone( "h_fmissfpass_bj_in_dphi_c" ) ;
      TH1F* h_rhl_bj_in_dphi_c = (TH1F*) h_rhl_bj_in_dphi -> Clone( "h_rhl_bj_in_dphi_c" ) ;

      h_rprime -> SetMarkerStyle( 20 ) ;
      h_rprime -> SetLineWidth( 2 ) ;
      h_rprime_c -> SetLineWidth( 2 ) ;
      h_rprime_c -> SetLineStyle( 2 ) ;

      h_fmissfpass_bj_in_dphi -> SetMarkerStyle( 21 ) ;
      h_fmissfpass_bj_in_dphi -> SetLineColor( 4 ) ;
      h_fmissfpass_bj_in_dphi_c -> SetLineColor( 4 ) ;
      h_fmissfpass_bj_in_dphi -> SetMarkerColor( 4 ) ;
      h_fmissfpass_bj_in_dphi -> SetLineWidth( 2 ) ;
      h_fmissfpass_bj_in_dphi_c -> SetLineWidth( 2 ) ;
      h_fmissfpass_bj_in_dphi_c -> SetLineStyle( 2 ) ;

      h_rhl_bj_in_dphi -> SetMarkerStyle( 21 ) ;
      h_rhl_bj_in_dphi -> SetLineColor( 2 ) ;
      h_rhl_bj_in_dphi_c -> SetLineColor( 2 ) ;
      h_rhl_bj_in_dphi -> SetMarkerColor( 2 ) ;
      h_rhl_bj_in_dphi -> SetLineWidth( 2 ) ;
      h_rhl_bj_in_dphi_c -> SetLineWidth( 2 ) ;
      h_rhl_bj_in_dphi_c -> SetLineStyle( 2 ) ;

      h_rprime -> SetYTitle( "Contribution to QCD high/low ratio" ) ;
      h_rprime -> SetTitleSize( 0.05, "y" ) ;
      h_rprime -> SetTitleOffset( 1.5, "y" ) ;
      if ( strcmp( htstr, "htl" ) == 0 ) {
         h_rprime -> SetMaximum(0.20) ;
         h_rprime -> SetTitle( "Low HT") ;
      } else if ( strcmp( htstr, "htm" ) == 0 ) {
         h_rprime -> SetMaximum(0.05) ;
         h_rprime -> SetTitle( "Medium HT") ;
      } else {
         h_rprime -> SetMaximum(0.03) ;
         h_rprime -> SetTitle( "High HT") ;
      }
      gStyle -> SetOptTitle(1) ;

      h_rprime -> Draw("l") ;
      h_fmissfpass_bj_in_dphi -> Draw("l same") ;
      h_fmissfpass_bj_in_dphi_c -> Draw("hist l same") ;
      h_rhl_bj_in_dphi -> Draw("l same") ;
      h_rhl_bj_in_dphi_c -> Draw("hist l same") ;
      h_rprime -> Draw("l same") ;
      h_rprime_c -> Draw("hist l same") ;

      l2 -> Draw() ;

      gPad -> SetGridy(1) ;

      sprintf( fname, "outputfiles/badjet-Rqcd-vs-mht-bin-%s.pdf", htstr ) ;
      can2 -> SaveAs( fname ) ;


    //----------------------------------------
    //   DR and syst

      printf("\n\n\n") ;
      TH1F* h_dr = new TH1F( "h_dr_"+ht_str, "double ratio", nbins-1, 0.5, nbins-1+0.5 ) ;
      TH1F* h_dr_stat = new TH1F( "h_dr_stat_"+ht_str, "double ratio", nbins-1, 0.5, nbins-1+0.5 ) ;
      TH1F* h_dr_syst = new TH1F( "h_dr_syst_"+ht_str, "double ratio", nbins-1, 0.5, nbins-1+0.5 ) ;
      TH1F* h_dr_total = new TH1F( "h_dr_total_"+ht_str, "double ratio", nbins-1, 0.5, nbins-1+0.5 ) ;
      double rc_val = h_rprime -> GetBinContent(1) ;
      double rc_err = h_rprime -> GetBinError(1) ;

         fprintf( out_file, "Sqcd_mhtc_%3s   1.000000     0.00000  0.00\n", htstr) ;

      double store_val = 0, store_rel_err = 0;
      char store_htstr[10];
      for ( int bi=2; bi<=nbins; bi++ ) {
         double r_val = h_rprime -> GetBinContent(bi) ;
         double r_err = h_rprime -> GetBinError(bi) ;
         double bjnidp_r = h_fmissfpass_bj_in_dphi -> GetBinContent( bi ) ;
         double dr_val(0.) ;
         double dr_stat(0.) ;
         double dr_syst(0.) ;
         if ( rc_val > 0 ) {
            dr_val = r_val / rc_val ;
            if ( r_val > 0 ) {
               dr_stat = dr_val * sqrt( pow( rc_err/rc_val, 2. ) + pow( r_err/r_val, 2. ) ) ;
            }
            dr_syst = bjnidp_r / rc_val ;
         }
         double dr_total = sqrt( dr_stat*dr_stat + dr_syst*dr_syst ) ;
         double dr_rel_err(0.) ;
         if ( dr_val > 0 ) dr_rel_err = dr_total / dr_val ;
         printf("  MHT%d :  DR = (%6.4f +/- %6.4f) / ( %6.4f +/- %6.4f) = %6.4f +/- %6.4f (%6.4f, %6.4f), %5.0f%%\n",
            bi-1, r_val, r_err, rc_val, rc_err, dr_val, dr_total, dr_stat, dr_syst, 100*dr_rel_err ) ;

        /// if the variable is equal to zero, equate it to the previous value in the list (if it is not the first variable of its kind)

         if ( dr_val == 0. && bi!=2)
         { 
            fprintf( out_file, "Sqcd_mht%d_%3s   %3.6f     0.00000  %3.2f\n", bi-1, htstr,store_val,store_rel_err) ;
            printf("===================================================================================================================\n");
            printf("Warning: This code automatically set Sqcd_mht%d_%3s = Sqcd_mht%d_%3s, because Sqcd_mht%d_%3s was equal to zero\n", bi-1, htstr,bi-2, store_htstr, bi-1,htstr);
            printf("===================================================================================================================\n");

         }
         else
            fprintf( out_file, "Sqcd_mht%d_%3s   %3.6f     0.00000  %3.2f\n", bi-1, htstr,dr_val,dr_rel_err) ;

         store_val = dr_val;
         store_rel_err = dr_rel_err;
         strcpy(store_htstr,htstr);
         h_dr -> SetBinContent( bi-1, dr_val ) ;
         h_dr -> SetBinError( bi-1, 0.0000001 ) ;
         h_dr_stat -> SetBinContent( bi-1, dr_val ) ;
         h_dr_syst -> SetBinContent( bi-1, dr_val ) ;
         h_dr_total -> SetBinContent( bi-1, dr_val ) ;
         h_dr_stat -> SetBinError( bi-1, dr_stat ) ;
         h_dr_syst -> SetBinError( bi-1, dr_syst ) ;
         h_dr_total -> SetBinError( bi-1, dr_total ) ;
         h_dr_total -> GetXaxis() -> SetBinLabel( bi-1, h_rprime->GetXaxis()->GetBinLabel(bi) ) ;
      } // bi

      h_dr_total -> SetFillColor( kRed-10 ) ;
      h_dr_stat -> SetFillColor( kRed-9 ) ;
      h_dr -> SetLineColor(4) ;
      h_dr -> SetLineWidth(2) ;

      h_dr_total -> Draw( "e2" ) ;
      h_dr_stat -> Draw( "e2 same" ) ;
      h_dr -> Draw("same") ;
      h_dr -> Draw("axis same" ) ;
      h_dr -> Draw("axig same" ) ;
   }//ht_level

   for ( int nb_count = 0; nb_count < nb_nb; nb_count++)
         fprintf( out_file, "Sqcd_nb%d        1.000000     0.00000  0.00\n", nb_count) ;

   } // draw_badjet_cat_v3

bool transfer_qcd_parameters(string filein_name, string fileout_name)
{


    ifstream filein(filein_name);
    ofstream fileout(fileout_name);

    if ( !filein ) { cout << "Warning: file " << filein_name <<" doesn't exist" << endl; return 0; }
    if ( !filein ) { cout << "Warning: cannot open/create file" << fileout_name << endl; return 0; }

    std::string line;
    while( std::getline(filein,line) )
    {

        std::size_t index = line.find("+/-");
        if (index == std::string::npos) { cout << "Warning: The structure of file " << filein_name << "is not as expected"; return 0; }
        line.replace( index, 3, "   ");

        index = line.find("(");
        if (index == std::string::npos) { cout << "Warning: The structure of file " << filein_name << "is not as expected"; return 0; }
        line.replace(line.find('('),line.find(')')-line.find('(')+1,"0.00");

        string var_str, name_str;        
        stringstream convert_temp(line);
        convert_temp >> name_str;
        convert_temp >> var_str;
        stringstream convert(var_str);
        double val;
        if ( !(convert >> val) )  { cout << val << "Warning: The structure of file " << filein_name << "is not as expected"; return 0; }
        fileout << line << std::endl;
    }    

   return 1;
}//transfer_qcd_parameters

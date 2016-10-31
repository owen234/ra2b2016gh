#ifndef modelfit3_c
#define modelfit3_c

#include "TROOT.h"

#include "TText.h"
#include "TLatex.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TMinuit.h"
#include "TMatrixT.h"

#include "histio.c"
#include "binning.h"

#include <iostream>
#include <fstream>

  using std::cout ;
  using std::endl ;


   double data_Rqcd[10][10] ;
   double data_Rqcd_err[10][10] ;

   double fit_Rqcd_HT[10] ;
   double fit_SFqcd_njet[10] ;

   bool only_fit_mht1 ;

   double calc_fit_error( TMinuit* tm, int hbi, int nji, double& simple_model_err ) ;
   void  draw_boundaries( float hmin, float hmax ) ;


   bool first_time ;

  //-----------

   void minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

      int idummy = npar ;
      double fdummy = gin[0] ;
      fdummy = 0. ;
      idummy = iflag ;

      f = 0. ;

      //--- unpack the stupid par vector.
      int parind(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         fit_Rqcd_HT[hbi] = par[parind] ;
         if ( first_time ) { printf( " minuit_fcn : hbi=%d, parind=%d, par value = %g\n", hbi, parind, par[parind] ) ; }
         parind ++ ;
      } // hbi.
      for ( int nji=0; nji<nb_nj; nji++ ) {
         if ( nji == njet_bin_to_be_fixed_in_qcd_model_fit ) {
            fit_SFqcd_njet[nji] = 1.0 ;
         } else {
            fit_SFqcd_njet[nji] = par[parind] ;
            if ( first_time ) { printf( " minuit_fcn : nji=%d, parind=%d, par value = %g\n", nji, parind, par[parind] ) ; }
            parind++ ;
         }
      } // nji.

      if ( first_time ) { printf("\n\n minuit_fcn : first time\n") ; }
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               for ( int nji=0; nji<nb_nj; nji++ ) {
                  if ( hbi==0 && nji>=(nb_nj - 2) ) continue ;  // skip top two njets bins for lowest HT.
                  if ( !( data_Rqcd_err[hbi][nji] > 0. ) ) { continue ; }
                  double delta = data_Rqcd[hbi][nji] - fit_Rqcd_HT[hbi] * fit_SFqcd_njet[nji] ;
                  f += delta*delta / (data_Rqcd_err[hbi][nji] * data_Rqcd_err[hbi][nji] ) ;
                  if ( first_time ) {
                     printf("  minuit_fcn: hbi=%d, nji=%d : data=%g , model = %g * %g = %g\n",
                        hbi, nji, data_Rqcd[hbi][nji],
                        fit_Rqcd_HT[hbi], fit_SFqcd_njet[nji],
                        fit_Rqcd_HT[hbi] * fit_SFqcd_njet[nji] ) ;
                  }
               } // nji.
         } // hbi.

      if ( first_time ) {
         first_time = false ;
         printf("\n\n") ;
      }
   } // minuit_fcn.

  //-----------

  //=============================================================================================================================

   void modelfit3(
                    const char* infile = "outputfiles/modelfit-input-qcdmc.root",
                    const char* outfilebase = "outputfiles/qcdmc-chi2-fit",
                    const char* ratio_histname = "h_ratio",
                    float plot_min = -0.05,
                    float plot_max = 0.3
                    ) {
      setup_bins();
      char fname[1000] ;

      char command[10000] ;
      sprintf( command, "basename %s", infile ) ;
      TString infile_nopath = gSystem -> GetFromPipe( command ) ;

      for ( int i=0; i<10; i++ ) {
         fit_Rqcd_HT[i] = 0. ;
         fit_SFqcd_njet[i] = 0. ;
      }

      first_time = true ;

      gDirectory->Delete("h*") ;

      loadHist( infile ) ;

      char hname[1000] ;
      char htitle[1000] ;

      TFile* tf = new TFile( infile, "READ" ) ;
      if ( tf == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", infile ) ; return ; }
      if ( ! tf -> IsOpen() )  { printf("\n\n *** Bad input file: %s\n\n", infile ) ; return ; }

      TH1F* h_ratio = (TH1F*) tf -> Get( ratio_histname ) ;
      if ( h_ratio == 0x0 ) { printf("\n\n *** missing %s hist in %s.\n\n", hname, infile ) ; return ; }

      int bbi = 0 ;
      int histbin(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         for ( int nji=0; nji<nb_nj; nji++ ) {
            histbin ++ ;
            if ( hbi==0 && nji>=(nb_nj-2) ) continue ; // skip top two njets bins for lowest HT bin.
               char binlabel[100] ;
               sprintf( binlabel, "%s", h_ratio -> GetXaxis() -> GetBinLabel( histbin ) ) ;
               printf( " check %s is nji=%d, hbi=%d\n", binlabel, nji+1, hbi+1 ) ;
               data_Rqcd[hbi][nji] = h_ratio -> GetBinContent( histbin ) ;
               data_Rqcd_err[hbi][nji] = h_ratio -> GetBinError( histbin ) ;
         } // nji
      } // hbi



      TH1F* h_ratio_nb_njet[10][10] ;



      int n_minuit_pars(0) ;
      n_minuit_pars = nBinsHT + nb_nj-1 ;

      TMinuit *myMinuit = new TMinuit( n_minuit_pars ) ; // arg is # of parameters

      myMinuit->SetFCN( minuit_fcn ) ;


      Double_t arglist[10] ;
      Int_t ierflg = 0 ;

      arglist[0] = 1 ;
      myMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //--- do this for chi2 fit.

      int parind(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         sprintf( pname, "Rqcd_HT%d", hbi+1 ) ;
         myMinuit->mnparm( parind, pname, data_Rqcd[hbi][0], 0.03, 0., 2., ierflg ) ;
         parind++ ;
      } // hbi.
      for ( int nji=0; nji<nb_nj; nji++ ) {
         if (nji == njet_bin_to_be_fixed_in_qcd_model_fit ) continue;
         char pname[1000] ;
         sprintf( pname, "SFqcd_njet%d", nji+1 ) ;
         myMinuit->mnparm( parind, pname, 1.0, 0.10, 0., 90., ierflg ) ;
         parind++ ;
      } // nji.

      myMinuit->Migrad() ;
      myMinuit->mncomd("hesse",ierflg) ;


      TH1F* h_model  = new TH1F( "h_model", "Model result", nBinsHT * nb_nj, 0.5, nBinsHT * nb_nj + 0.5 ) ;
      TH1F* h_model2 = new TH1F( "h_model2", "Model result", nBinsHT * nb_nj, 0.5, nBinsHT * nb_nj + 0.5 ) ;
      h_model -> SetLineWidth(2) ;
      h_model -> SetLineColor(kRed+1) ;
      h_model2 -> SetFillColor(kRed-10) ;

      histbin = 0 ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         for ( int nji=0; nji<nb_nj; nji++ ) {
               histbin ++ ;
               if ( hbi==0 && nji>=(nb_nj-2) ) continue ; // skip top two njets bins for lowest HT bin.
               char binlabel[100] ;
               sprintf( binlabel, "%s", h_ratio -> GetXaxis() -> GetBinLabel( histbin ) ) ;
               double model_val = fit_Rqcd_HT[hbi] * fit_SFqcd_njet[nji] ;
               double model_err = 0. ;
               double simple_model_err = 0. ;
               model_err = calc_fit_error( myMinuit, hbi, nji, simple_model_err ) ;
               printf( " %s is nji=%d, hbi=%d    QCD MC ratio = %5.3f +/- %5.3f,   Model = %5.3f +/- %5.3f  (simple %5.3f)\n",
                    binlabel, nji+1, hbi+1,
                    data_Rqcd[hbi][nji], data_Rqcd_err[hbi][nji],
                    model_val, model_err, simple_model_err ) ;
               h_model -> SetBinContent( histbin, model_val ) ;
               h_model2 -> SetBinContent( histbin, model_val ) ;
               h_model2 -> SetBinError( histbin, model_err ) ;
               h_model  -> GetXaxis() -> SetBinLabel( histbin, binlabel ) ;
               h_model2 -> GetXaxis() -> SetBinLabel( histbin, binlabel ) ;
         } // nji
      } // hbi

      h_model  -> GetXaxis() -> LabelsOption("v") ;
      h_model2 -> GetXaxis() -> LabelsOption("v") ;


      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin(0.2) ;
      gStyle -> SetPadLeftMargin(0.15) ;

      TCanvas* can1 = new TCanvas( "can1_modelfit3", "Model fit", 900, 700 ) ;

      h_ratio -> SetMarkerStyle( 20 ) ;

      h_model2 -> SetMinimum( plot_min ) ;
      h_model2 -> SetMaximum( plot_max ) ;

      h_model2 -> SetYTitle("R^{QCD}") ;
      h_model2 -> SetTitleSize( 0.055, "y" ) ;
      h_model2 -> GetXaxis() -> SetLabelSize( 0.05 ) ;

      h_model2 -> Draw("e2") ;
      h_model -> Draw("same") ;
      h_ratio -> Draw("same e0") ;
      h_model2 -> Draw("axis same") ;
      h_model2 -> Draw("axig same") ;

      gPad -> SetGridy(1) ;

      can1 -> Update() ; can1 -> Draw() ;

      sprintf( fname, "%s-modelfit.pdf", outfilebase ) ;
      can1 -> SaveAs( fname ) ;



      sprintf( fname, "%s-model-pars.txt", outfilebase ) ;
      FILE* ofp ;
      if ( (ofp = fopen( fname, "w" ))==NULL ) {
         printf( "\n\n *** Problem opening output file: %s\n\n", fname ) ;
         return ;
      }



      printf("\n\n") ;

      parind = 0 ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "Kqcd_HT%d", hbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         double rel_err(0.) ;
         if ( val != 0 ) { rel_err = err/val ; }
         printf(" %11s  %8.5f +/- %8.5f  (%4.2f)\n", pname, val, err, rel_err ) ;
         fprintf( ofp, " %12s  %8.5f +/- %8.5f  (%4.2f)\n", pname, val, err, rel_err ) ;
         fit_Rqcd_HT[hbi] = val ;
         parind++ ;
      } // hbi.


      for ( int nji=0; nji<nb_nj; nji++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "Sqcd_njet%d", nji+1 ) ;
         double rel_err(0.) ;

         if ( nji == njet_bin_to_be_fixed_in_qcd_model_fit ) 
         {
            val = 1; err = 0; rel_err = 0; 
            printf(" %11s  %6.3f +/- %5.3f  (%4.2f)\n", pname, val, err, rel_err ) ;
            fprintf( ofp, " %12s  %8.5f +/- %8.5f  (%4.2f)\n", pname, val, err, rel_err ) ;
         }
         else
         {
         myMinuit->GetParameter( parind, val, err ) ;
         if ( val != 0 ) { rel_err = err/val ; }
         printf(" %11s  %6.3f +/- %5.3f  (%4.2f)\n", pname, val, err, rel_err ) ;
         fprintf( ofp, " %12s  %8.5f +/- %8.5f  (%4.2f)\n", pname, val, err, rel_err ) ;
         fit_SFqcd_njet[nji] = val ;
         parind++ ;
         }
      } // nji.
      printf("\n\n") ;

      fclose( ofp ) ;
      printf("Saved model pars in %s\n\n", fname ) ;

      return ; //*************


   } // modelfit2d


  //==========================================================================================

   double calc_fit_error( TMinuit* tm, int hbi, int nji, double& simple_model_err ) {

      simple_model_err = 0. ;

      //bool verb = true ;
      bool verb = false ;

      if ( verb ) {
         printf("\n\n ==================================================================================================\n") ;
         printf("     calc_fit_error: computing model error for HT%d, Njet%d\n\n", hbi+1, nji+1 ) ;
      }

      if ( tm == 0x0 ) return -1 ;

      int n_minuit_pars = tm -> GetNumPars() ;
      if ( n_minuit_pars <= 0 ) return -1 ;
      if (verb) printf(" Number of minuit parameters: %d\n", n_minuit_pars ) ;

      Double_t cov_mat[n_minuit_pars][n_minuit_pars] ;
      tm -> mnemat( &cov_mat[0][0], n_minuit_pars ) ;


     //--- Minuit parameter indices, counting from 0
      int ht_pind = hbi ;
      int njet_pind = -1 ;
      if ( nji > 0 ) njet_pind = nBinsHT + nji - 1 ;

      if ( verb ) printf("  Minuit parameter indices:  hbi=%d, pi=%d;     nji=%d, pi=%d\n",
        hbi, ht_pind,  nji, njet_pind ) ;

      double ht_par_val, ht_par_err ;
      tm -> GetParameter( ht_pind, ht_par_val, ht_par_err ) ;
      double njet_par_val, njet_par_err ;
      if ( nji > 0 ) {
         tm -> GetParameter( njet_pind, njet_par_val, njet_par_err ) ;
      } else {
         njet_par_val = 1 ;  njet_par_err = 0. ;
      }

      if ( verb ) printf("   HT par: %5.3f +/- %5.3f\n", ht_par_val, ht_par_err ) ;
      if ( verb ) printf(" Njet par: %5.3f +/- %5.3f\n", njet_par_val, njet_par_err ) ;

      if ( ht_par_val > 0 && njet_par_val > 0 ) {
         simple_model_err = ht_par_val * njet_par_val * sqrt( pow( ht_par_err/ht_par_val, 2. ) + pow( njet_par_err/njet_par_val, 2. ) ) ;
      }



     //--- Calculate the three partial derivatives.
      double df_dhtpar   =  njet_par_val ;
      double df_dnjetpar =  ht_par_val   ;


     //--- Create the partial derivative vector.
      TMatrixT<double> pd_col_vec( n_minuit_pars, 1 ) ;
      TMatrixT<double> pd_row_vec( 1, n_minuit_pars ) ;
      for ( int i=0; i<n_minuit_pars; i++ ) {
         pd_col_vec(i,0) = 0. ;
         pd_row_vec(0,i) = 0. ;
      }
      pd_col_vec( ht_pind, 0 ) = df_dhtpar ;
      pd_row_vec( 0, ht_pind ) = df_dhtpar ;
      if ( nji > 0 ) {
         pd_col_vec( njet_pind, 0 ) = df_dnjetpar ;
         pd_row_vec( 0, njet_pind ) = df_dnjetpar ;
      }
  //  if ( verb ) {
  //     printf("   Partial derivative column vector:\n") ;
  //     pd_col_vec.Print() ;
  //  }
  //  if ( verb ) {
  //     printf("   Partial derivative row vector:\n") ;
  //     pd_row_vec.Print() ;
  //  }

      TMatrixT<double> cov_mat_tm( n_minuit_pars, n_minuit_pars ) ;
      for ( int i=0; i<n_minuit_pars; i++ ) {
         for ( int j=0; j<n_minuit_pars; j++ ) {
            cov_mat_tm( i, j ) = cov_mat[i][j] ;
         } // j
      } // i

//    if ( verb ) {
//       printf("\n\n  ===== Covariance matrix:\n") ;
//       cov_mat_tm.Print() ;
//    }


      TMatrixT<double> cov_times_pd_col( n_minuit_pars, 1 ) ;
      cov_times_pd_col.Mult( cov_mat_tm, pd_col_vec ) ;
 //   if ( verb ) {
 //      printf( "\n\n  ===== Cov mat * pd col vec:\n") ;
 //      cov_times_pd_col.Print() ;
 //   }

      TMatrixT<double> pd_row_times_prod( 1, 1 ) ;
      pd_row_times_prod.Mult( pd_row_vec, cov_times_pd_col ) ;
 //   if ( verb ) {
 //      printf( "\n\n  ===== pd row vec * prod:\n") ;
 //      pd_row_times_prod.Print() ;
 //   }


      double fit_error(0.) ;

      if ( pd_row_times_prod(0,0) > 0 ) {
         fit_error = sqrt( pd_row_times_prod(0,0) ) ;
      }

      if ( verb ) printf("\n Final answer for fit error is %5.3f\n\n", fit_error ) ;

      return fit_error ;

   } // calc_fit_error

  //==========================================================================================

  //--------

   void draw_boundaries( float hmin, float hmax ) {

      TLine* line = new TLine() ;
      line -> SetLineStyle(2) ;

      line -> DrawLine( 11.5, hmin, 11.5, hmax ) ;
      line -> DrawLine( 22.5, hmin, 22.5, hmax ) ;
      line -> DrawLine( 33.5, hmin, 33.5, hmax ) ;
      line -> DrawLine( 44.5, hmin, 44.5, hmax ) ;

      line -> SetLineColor(2) ;
      line -> SetLineStyle(1) ;
      line -> DrawLine( 0.5, 0., 55.5, 0. ) ;

   }

  //--------






#endif

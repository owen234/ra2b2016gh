
#include "TH2F.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "TMatrixT.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TColor.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"

   using namespace RooFit;
   using namespace RooStats;

#include "histio.c"

#include <fstream>

   int npars(13) ;
   float par_val[20] ;
   float par_err[20] ;
   char par_name[20][100] ;
   char par_name_short[20][100] ;

   int par_map_ht_njet[10][10] ;

   int ht_par[20] ;
   int njet_par[20] ;


   RooWorkspace* the_ws ;
   RooAbsPdf* likelihood ;
   RooDataSet* rds ;

   TMatrixD reordered_reoriented_eigen_vector_matrix( npars, npars ) ;
   TVectorD reordered_eigen_vals( npars ) ;

  //---------

   TH1F* get_hist( const char* hname ) ;
   TH2F* make_contour_original_pars_from_lh( int pix, int piy, bool fix_all_model_pars=true ) ;
   TH2F* make_contour_rotated_pars_from_lh( int pix, int piy ) ;

  //---------

   void gen_nuisance_parameter_matrix1( const char* infile = "outputfiles/lhfit-results-ws-lhfit-test3/kqcd-parameter-fit-covmat.txt",
                                        const char* outfile = "outputfiles/nuisance-parameter-matrix.txt",
                                        const char* outfile_simple = "outputfiles/nuisance-parameter-matrix-simple.txt",
                                        const char* wsfile = "outputfiles/ws-lhfit-test3.root" ) {

      gDirectory -> Delete( "h*" ) ;

      const Int_t Number = 3 ;
      Double_t Length[Number] = {0.,0.5, 1.} ;
      Double_t Red[Number] = {0.,1.,1.} ;
      Double_t Green[Number] = {0.,1.,0.} ;
      Double_t Blue[Number] = {1.,1.,0.} ;
      Int_t nb = 90 ;
      TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

      gStyle -> SetOptStat(0) ;

      gStyle->SetPaintTextFormat("7.4f");
      gStyle -> SetPadLeftMargin( 0.15 ) ;
      gStyle -> SetPadRightMargin( 0.15 ) ;

      int cx = 50 ;
      int cy = 50 ;




      FILE* ofp(0x0) ;
      if ( (ofp=fopen( outfile, "w" ) )==NULL ) {
         printf("\n\n *** Problem opening output file: %s\n\n", outfile ) ;
         gSystem -> Exit(-1) ;
      }

      FILE* ofp_simple(0x0) ;
      if ( (ofp_simple=fopen( outfile_simple, "w" ) )==NULL ) {
         printf("\n\n *** Problem opening output file: %s\n\n", outfile_simple ) ;
         gSystem -> Exit(-1) ;
      }




      ifstream ifs ;
      ifs.open( infile ) ;
      if ( !ifs.good() ) { printf("\n\n *** Bad input file: %s\n\n", infile ) ; return ; }


      TString line ;

      printf("\n\n Reading parameters from file %s\n", infile ) ;

      for ( int pi=0; pi<npars; pi++ ) {
         line.ReadLine( ifs ) ;
         sscanf( line.Data(), "%s %f +/- %f", par_name[pi], &(par_val[pi]), &(par_err[pi]) ) ;
         TString pns( par_name[pi] ) ;
         pns.ReplaceAll( "R_qcd_ldp_", "" ) ;
         sprintf( par_name_short[pi], "%s", pns.Data() ) ;
         int bi_nj, bi_ht ;
         sscanf( par_name[pi], "R_qcd_ldp_Nj%d_HT%d", &bi_nj, &bi_ht ) ;
         if ( bi_nj <=0 ) { printf("\n\n *** Illegal Njet index %d : %s\n\n", bi_nj, par_name[pi] ) ; gSystem -> Exit(-1) ; }
         if ( bi_ht <=0 ) { printf("\n\n *** Illegal HT   index %d : %s\n\n", bi_ht, par_name[pi] ) ; gSystem -> Exit(-1) ; }
         par_map_ht_njet[bi_ht][bi_nj] = pi ;
         printf("  %2d : %s , bi_ht=%d, bi_nj=%d,  %6.4f +/- %6.4f\n", pi, par_name[pi], bi_ht, bi_nj, par_val[pi], par_err[pi] ) ;
      } // pi



      line.ReadLine( ifs ) ;

      float cov_mat[20][20] ;
      for ( int ri=0; ri<npars; ri++ ) {
         char label[100], pname[100] ;
         TString token ;
         token.ReadToken( ifs ) ; sprintf( label, "%s", token.Data() ) ;
         token.ReadToken( ifs ) ; sprintf( pname, "%s", token.Data() ) ;
         if ( strcmp( pname, par_name[ri] ) != 0 ) { printf("\n\n *** parameter mismatch: %s %s\n\n", pname, par_name[ri] ) ; return ; }
         for ( int i=0; i<npars; i++ ) {
            token.ReadToken( ifs ) ;
            sscanf( token.Data(), "%f", &(cov_mat[ri][i]) ) ;
         }
      } // ri

      printf(" Done reading parameters from file.\n") ;


      printf("\n\n  Fit parameters:\n") ;
      for ( int pi=0; pi<npars; pi++ ) {
         double rel_err(0.) ;
         if ( par_val[pi] != 0 ) rel_err = par_err[pi] / par_val[pi] ;
         printf("  %d : %12s : %12.5f +/- %12.5f   (%6.4f)\n", pi, par_name[pi], par_val[pi], par_err[pi], rel_err ) ;
      }
      printf("\n\n   Fit covariance matrix:\n") ;
      for ( int ri=0; ri<npars; ri++ ) {
         for ( int ci=0; ci<npars; ci++ ) {
            printf("  %12.8f  ", cov_mat[ri][ci] ) ;
         } // ci
         printf("\n") ;
      } // ri
      printf("\n\n") ;








      TH2F* h_fit_correlation_matrix = new TH2F( "h_fit_correlation_matrix", "Fit correlation matrix", npars, 0.5, npars+0.5,  npars, 0.5, npars+0.5 ) ;
      printf("\n\n   Fit correlation matrix:\n") ;
      for ( int ri=0; ri<npars; ri++ ) {
         for ( int ci=0; ci<npars; ci++ ) {
            float rho = cov_mat[ri][ci] / sqrt( cov_mat[ri][ri] * cov_mat[ci][ci] ) ;
            printf("  %12.8f  ", rho ) ;
            h_fit_correlation_matrix -> SetBinContent( ci+1, npars-ri, rho ) ;
            h_fit_correlation_matrix -> GetXaxis() -> SetBinLabel( ci+1, par_name_short[ci] ) ;
            h_fit_correlation_matrix -> GetYaxis() -> SetBinLabel( npars-ri, par_name_short[ri] ) ;
         } // ci
         printf("\n") ;
      } // ri
      printf("\n\n") ;





      TMatrixDSym cov_mat_tm( npars ) ;
      for ( int ri=0; ri<npars; ri++ ) {
         for ( int ci=0; ci<npars; ci++ ) {
            cov_mat_tm[ri][ci] = cov_mat[ri][ci] ;
         } // ci
      } // ri

      TVectorD eigen_vals( npars ) ;
      TMatrixD eigen_vector_matrix = cov_mat_tm.EigenVectors( eigen_vals ) ;


      printf("\n\n") ;
      for ( int i=0; i<npars; i++ ) {
         printf("   Eigen value %d : %12.8f\n", i, eigen_vals[i] ) ;
      } // i

      printf("\n\n Eigen vector matrix:\n") ;
      for ( int i=0; i<npars; i++ ) {
         for ( int j=0; j<npars; j++ ) {
            printf("  %12.8f  ", eigen_vector_matrix[i][j] ) ;
         } // j
         printf("\n" ) ;
      } // i
      printf("\n\n") ;


      int eigen_vector_index[20] ;
      for ( int i=0; i<npars; i++ ) { eigen_vector_index[i] = -1 ; }

      printf("\n\n  Check of eigen vectors:\n") ;
      for ( int i=0; i<npars; i++ ) {
         printf("  Eigen value %d : %12.8f\n", i, eigen_vals[i] ) ;
         TMatrixT<double> eigen_vector_col(npars,1) ;
         TMatrixT<double> eigen_vector_row(1,npars) ;
         double largest_component_val(0.) ;
         int largest_component_ind(-1) ;
         printf("      Eigen vector elements : ") ;
         for ( int j=0; j<npars; j++ ) {
            eigen_vector_col(j,0) = eigen_vector_matrix[j][i] ;
            eigen_vector_row(0,j) = eigen_vector_matrix[j][i] ;
            if ( fabs( eigen_vector_matrix[j][i] ) > largest_component_val ) {
               if ( eigen_vector_index[j] >= 0 ) {
                  //printf(" already used.\n") ;
               } else {
                  largest_component_val = fabs( eigen_vector_matrix[j][i] ) ;
                  largest_component_ind = j ;
               }
            }
            printf( "  %12.8f  ", eigen_vector_col(j,0) ) ;
         }
         printf("\n") ;
         printf("   Largest component is %s   %12.8f\n", par_name[largest_component_ind], eigen_vector_matrix[largest_component_ind][i] ) ;
         eigen_vector_index[largest_component_ind] = i ;
         TMatrixT<double> cm_times_ev(npars,1) ;
         cm_times_ev.Mult( cov_mat_tm, eigen_vector_col ) ;
         printf("      CM times eigen vector : ") ;
         for ( int j=0; j<npars; j++ ) {
            printf( "  %12.8f  ", cm_times_ev(j,0) ) ;
         }
         printf("\n") ;
         printf("      after div  by EV      : ") ;
         for ( int j=0; j<npars; j++ ) {
            printf( "  %12.8f  ", cm_times_ev(j,0) / eigen_vals[i] ) ;
         }
         printf("\n") ;
         TMatrixT<double> ev_times_cm_times_ev(1,1) ;
         ev_times_cm_times_ev.Mult( eigen_vector_row, cm_times_ev ) ;
         printf("  EVec times CM times EVec (should EVal) :  %12.8f\n", ev_times_cm_times_ev(0,0) ) ;

         printf("\n\n") ;
      } //


      printf("\n\n") ;
      for ( int pi=0; pi<npars; pi++ ) {
         printf(" par %2d : %12s : eigen vector index %d\n", pi, par_name[pi], eigen_vector_index[pi] ) ;
         if ( eigen_vector_index[pi] < 0 || eigen_vector_index[pi] > npars ) { printf("\n\n *** No eigenvector associated with this par.\n\n" ) ; return ; }
      }
      printf("\n\n") ;






      for ( int ci=0; ci<npars; ci++ ) {
         int evi = eigen_vector_index[ci] ;
         reordered_eigen_vals[ci] = eigen_vals[evi] ;
         float sign_factor(1.) ;
         if ( eigen_vector_matrix(ci,evi) < 0 ) sign_factor = -1. ;
         for ( int ri=0; ri<npars; ri++ ) {
            reordered_reoriented_eigen_vector_matrix(ri,ci) = sign_factor * eigen_vector_matrix( ri, evi ) ;
         } // ri
      } // ci

      printf("\n\n Reordered eigen values:\n") ;
      for ( int i=0; i<npars; i++ ) {
         printf("  %12.8f  ", reordered_eigen_vals(i) ) ;
      }



      TH2F* h_eigenvec_matrix = new TH2F( "h_eigenvec_matrix", "Eigenvectors (columns)", npars, 0.5, npars+0.5,  npars, 0.5, npars+0.5 ) ;

      printf("\n") ;
      printf("\n\n Reordered, reoriented, eigen vector matrix:\n") ;
      for ( int i=0; i<npars; i++ ) {
         for ( int j=0; j<npars; j++ ) {
            printf("  %12.8f  ", reordered_reoriented_eigen_vector_matrix[i][j] ) ;
            char xlabel[100] ;
            sprintf( xlabel, "EV-%s", par_name_short[j] ) ;
            h_eigenvec_matrix -> SetBinContent( j+1, npars-i, reordered_reoriented_eigen_vector_matrix[i][j] ) ;
            h_eigenvec_matrix -> GetXaxis() -> SetBinLabel( j+1, xlabel ) ;
            h_eigenvec_matrix -> GetYaxis() -> SetBinLabel( npars-i, par_name_short[i] ) ;
         } // j
         printf("\n" ) ;
      } // i
      printf("\n\n") ;







      TMatrixD rotation_matrix( TMatrixD::kTransposed, eigen_vector_matrix ) ;
      printf("\n\n Rotation matrix:\n") ;
      for ( int ri=0; ri<npars; ri++ ) {
         for ( int ci=0; ci<npars; ci++ ) {
            printf("  %12.8f  ", rotation_matrix(ri,ci) ) ;
         } // ci
         printf("\n") ;
      } // ri



      TMatrixD reordered_reoriented_rotation_matrix( TMatrixD::kTransposed, reordered_reoriented_eigen_vector_matrix ) ;
      printf("\n\n Reordered, reoriented Rotation matrix:\n") ;
      for ( int ri=0; ri<npars; ri++ ) {
         for ( int ci=0; ci<npars; ci++ ) {
            printf("  %12.8f  ", reordered_reoriented_rotation_matrix(ri,ci) ) ;
         } // ci
         printf("\n") ;
      } // ri







      printf("\n\n  Calculation of total relative error using covariance matrix in simple basis.\n") ;
      for ( int pi=0; pi<npars; pi++ ) {
         double err = sqrt( cov_mat[pi][pi] ) / par_val[pi] ;
         printf("   %2d %s  %6.4f\n", pi, par_name[pi], err ) ;
      } // pi












      printf("\n\n  Calculation of total relative error and individual contributions in rotated basis.\n") ;
      printf("          ") ;
      for ( int pi=0; pi<npars; pi++ ) {
         printf("  %s  ", par_name_short[pi] ) ;
      }
      printf("\n") ;
      for ( int rpi=0; rpi<npars; rpi++ ) {
         if ( par_val[rpi] == 0 ) { printf("\n\n *** zero parameter value.\n\n") ; gSystem->Exit(-1) ; }
         float sum_err2(0.) ;
         printf("  %s  ", par_name_short[rpi] ) ;
         for ( int cpi=0; cpi<npars; cpi++ ) {
            int evi = eigen_vector_index[cpi] ;
            double eigen_val = eigen_vals[evi] ;
            double rot_mat_over_par_val = rotation_matrix( evi, rpi ) / par_val[rpi] ;
            double rel_err = sqrt( eigen_val ) * fabs( rot_mat_over_par_val ) ;
            printf(  "  %6.4f   ", rel_err ) ;
            sum_err2 += rel_err * rel_err ;
         } // cpi
         printf("     (%6.4f)\n", sqrt( sum_err2 ) ) ;
      } // rpi

















      printf("\n\n  Calculation of total relative error and individual contributions in rotated basis, using reordered reoriented rotation matrix.\n") ;
      TH2F* h_rel_err_table = new TH2F( "h_rel_err_table",
           "Relative error table, orthogonal parameters basis", npars+2, 0.5, npars+2+0.5,  npars, 0.5, npars+0.5 ) ;

      printf("   %2d     ", npars) ;
      fprintf( ofp, "    %2d    ", npars) ;
      for ( int pi=0; pi<npars; pi++ ) {
         printf("  %s  ", par_name_short[pi] ) ;
         fprintf( ofp, "  %s  ", par_name_short[pi] ) ;
      }
      fprintf( ofp, "\n" ) ;
      printf("\n") ;
      for ( int rpi=0; rpi<npars; rpi++ ) {
         if ( par_val[rpi] == 0 ) { printf("\n\n *** zero parameter value.\n\n") ; gSystem->Exit(-1) ; }
         float sum_err2(0.) ;
         printf("  %s  ", par_name_short[rpi] ) ;
         fprintf( ofp, "  %s  ", par_name_short[rpi] ) ;
         for ( int cpi=0; cpi<npars; cpi++ ) {
            int evi = eigen_vector_index[cpi] ;
            double eigen_val = eigen_vals[evi] ;
            double rot_mat_over_par_val = rotation_matrix( evi, rpi ) / par_val[rpi] ;
            double rel_err = sqrt( eigen_val ) * fabs( rot_mat_over_par_val ) ;
            printf(  "  %6.4f   ", rel_err ) ;
            fprintf( ofp, "  %6.4f   ", rel_err ) ;
            h_rel_err_table -> SetBinContent( cpi+1, npars-rpi, rel_err ) ;
            h_rel_err_table -> GetYaxis() -> SetBinLabel( npars-rpi, par_name_short[rpi] ) ;
            h_rel_err_table -> GetXaxis() -> SetBinLabel( cpi+1, par_name_short[cpi] ) ;
            sum_err2 += rel_err * rel_err ;
         } // cpi
         printf("     (%6.4f)\n", sqrt( sum_err2 ) ) ;
         fprintf( ofp, "     (%6.4f)\n", sqrt( sum_err2 ) ) ;
         h_rel_err_table -> SetBinContent( npars+2, npars-rpi, sqrt( sum_err2 ) ) ;
         h_rel_err_table -> GetXaxis() -> SetBinLabel( npars+2, "Total error" ) ;
      } // rpi






    //---
      TCanvas* can1 = new TCanvas("can1","Combine table, rotated parameters",1300,700 ) ;
      cx += 50 ; cy += 50 ;
      can1 -> SetWindowPosition( cx, cy ) ;

      h_rel_err_table -> SetContour(nb) ;
      h_rel_err_table -> SetMinimum( -0.7 ) ;
      h_rel_err_table -> SetMaximum(  0.7 ) ;

      h_rel_err_table -> Draw("colz") ;
      h_rel_err_table -> Draw("text same") ;














      TH2F* h_rel_err_table_simple_ignore = new TH2F( "h_rel_err_table_simple_ignore",
              "Relative error table, original parameters basis, ignore off-diagonal cov",
              npars+2, 0.5, npars+2+0.5,  npars, 0.5, npars+0.5 ) ;

      printf("\n\n") ;

      printf("\n\n  Calculation ignoring off-diagonal elements of the covariance matrix.\n") ;
      fprintf( ofp_simple, "   %2d     ", npars ) ;
      printf("    %2d    ", npars) ;
      for ( int pi=0; pi<npars; pi++ ) {
         printf("  %s  ", par_name_short[pi] ) ;
         fprintf( ofp_simple, "  %s  ", par_name_short[pi] ) ;
      }
      printf("\n") ;
      fprintf( ofp_simple, "\n") ;

      for ( int rpi=0; rpi<npars; rpi++ ) {
         printf("  %s  ", par_name_short[rpi] ) ;
         fprintf( ofp_simple, "  %s  ", par_name_short[rpi] ) ;
         float total_rel_err(0.) ;
         for ( int cpi=0; cpi<npars; cpi++ ) {
            float rel_err(0.) ;
            if ( cpi==rpi ) {
               rel_err = par_err[cpi]/par_val[cpi] ;
               total_rel_err = rel_err ;
            }
            h_rel_err_table_simple_ignore -> SetBinContent( cpi+1, 13-rpi, rel_err ) ;
            h_rel_err_table_simple_ignore -> GetYaxis() -> SetBinLabel( 13-rpi, par_name_short[rpi] ) ;
            h_rel_err_table_simple_ignore -> GetXaxis() -> SetBinLabel( cpi+1, par_name_short[cpi] ) ;
            printf(  "  %6.4f   ", rel_err ) ;
            fprintf( ofp_simple, "  %6.4f   ", rel_err ) ;
         } // cpi
         h_rel_err_table_simple_ignore -> SetBinContent( npars+2, 13-rpi, total_rel_err ) ;
         h_rel_err_table_simple_ignore -> GetXaxis() -> SetBinLabel( npars+2, "Total error" ) ;
         printf("     (%6.4f)\n", total_rel_err ) ;
         fprintf( ofp_simple, "     (%6.4f)\n", total_rel_err ) ;
      } // rpi

      printf("\n\n") ;



    //---
      TCanvas* can2 = new TCanvas("can2","Combine table, original parameters, ignoring correlations",1300,700 ) ;
      cx += 50 ; cy += 50 ;
      can2 -> SetWindowPosition( cx, cy ) ;

      h_rel_err_table_simple_ignore -> SetContour(nb) ;
      h_rel_err_table_simple_ignore -> SetMinimum( -0.7 ) ;
      h_rel_err_table_simple_ignore -> SetMaximum(  0.7 ) ;

      h_rel_err_table_simple_ignore -> Draw("colz") ;
      h_rel_err_table_simple_ignore -> Draw("text same") ;






    //---
      TCanvas* can3 = new TCanvas("can3","Fit correlation coefficient matrix",1300,700 ) ;
      can3 -> SetWindowPosition( cx, cy ) ;

      h_fit_correlation_matrix -> SetContour(nb) ;
      h_fit_correlation_matrix -> SetMinimum( -1.0 ) ;
      h_fit_correlation_matrix -> SetMaximum(  1.0 ) ;

      h_fit_correlation_matrix -> Draw("colz X+") ;
      h_fit_correlation_matrix -> Draw("text same") ;



    //---
      TCanvas* can4 = new TCanvas("can4","Cov mat Eigen vectors",1300,700 ) ;
      cx += 50 ; cy += 50 ;
      can4 -> SetWindowPosition( cx, cy ) ;


      h_eigenvec_matrix -> SetContour(nb) ;
      h_eigenvec_matrix -> SetMinimum( -0.7 ) ;
      h_eigenvec_matrix -> SetMaximum(  0.7 ) ;

      h_eigenvec_matrix -> Draw("colz") ;
      h_eigenvec_matrix -> Draw("text same") ;


      fclose( ofp_simple ) ;
      printf("\n\n Wrote output file: %s\n", outfile_simple ) ;

      fclose( ofp ) ;
      printf("\n\n Wrote output file: %s\n", outfile ) ;





  return ; //****************************************************************************


   //--- stuff below here is for diagnostics.


     //+++++++++++


      TFile* wstf = new TFile( wsfile ) ;
      the_ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );

      rds = (RooDataSet*) the_ws->obj( "observed_rds" ) ;
      cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      likelihood = the_ws->pdf("likelihood") ;
      if ( likelihood == 0x0 ) { printf("\n\n *** can't find likelihood in workspace.\n\n" ) ; return ; }
      printf("\n\n Likelihood:\n") ;
      likelihood -> Print() ;

      RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), Optimize(0), PrintLevel(0), Hesse(true), Strategy(1) ) ;
      double minNllFloat = fitResult->minNll() ;


      RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
      RooMsgService::instance().getStream(1).removeTopic(Fitting) ;


     //-------

////  TH2F* h_cont_op_p0_vs_p1_from_lh = make_contour_original_pars_from_lh( 0, 1, false ) ;
////  h_cont_op_p0_vs_p1_from_lh -> SetContour(nb) ;
////  h_cont_op_p0_vs_p1_from_lh -> Draw("colz") ;


////  TH2F* h_cont_op_p0_vs_p1_from_lh_fixallmp = make_contour_original_pars_from_lh( 0, 1, true ) ;
////  h_cont_op_p0_vs_p1_from_lh_fixallmp -> SetContour(nb) ;
////  h_cont_op_p0_vs_p1_from_lh_fixallmp -> Draw("colz") ;


////  TH2F* h_cont_rp_p0_vs_p1_from_lh = make_contour_rotated_pars_from_lh( 0, 1 ) ;
////  h_cont_rp_p0_vs_p1_from_lh -> SetContour(nb) ;
////  h_cont_rp_p0_vs_p1_from_lh -> Draw("colz") ;


     //-------


////  TH2F* h_cont_op_p1_vs_p2_from_lh = make_contour_original_pars_from_lh( 1, 2, false ) ;
////  h_cont_op_p1_vs_p2_from_lh -> SetContour(nb) ;
////  h_cont_op_p1_vs_p2_from_lh -> Draw("colz") ;


////  TH2F* h_cont_op_p1_vs_p2_from_lh_fixallmp = make_contour_original_pars_from_lh( 1, 2, true ) ;
////  h_cont_op_p1_vs_p2_from_lh_fixallmp -> SetContour(nb) ;
////  h_cont_op_p1_vs_p2_from_lh_fixallmp -> Draw("colz") ;


////  TH2F* h_cont_rp_p1_vs_p2_from_lh = make_contour_rotated_pars_from_lh( 1, 2 ) ;
////  h_cont_rp_p1_vs_p2_from_lh -> SetContour(nb) ;
////  h_cont_rp_p1_vs_p2_from_lh -> Draw("colz") ;

     //-------


////  TH2F* h_cont_op_p2_vs_p3_from_lh = make_contour_original_pars_from_lh( 2, 3, false ) ;
////  h_cont_op_p2_vs_p3_from_lh -> SetContour(nb) ;
////  h_cont_op_p2_vs_p3_from_lh -> Draw("colz") ;


////  TH2F* h_cont_op_p2_vs_p3_from_lh_fixallmp = make_contour_original_pars_from_lh( 2, 3, true ) ;
////  h_cont_op_p2_vs_p3_from_lh_fixallmp -> SetContour(nb) ;
////  h_cont_op_p2_vs_p3_from_lh_fixallmp -> Draw("colz") ;


////  TH2F* h_cont_rp_p2_vs_p3_from_lh = make_contour_rotated_pars_from_lh( 2, 3 ) ;
////  h_cont_rp_p2_vs_p3_from_lh -> SetContour(nb) ;
////  h_cont_rp_p2_vs_p3_from_lh -> Draw("colz") ;

     //-------


////  TH2F* h_cont_op_p3_vs_p4_from_lh = make_contour_original_pars_from_lh( 3, 4, false ) ;
////  h_cont_op_p3_vs_p4_from_lh -> SetContour(nb) ;
////  h_cont_op_p3_vs_p4_from_lh -> Draw("colz") ;


////  TH2F* h_cont_op_p3_vs_p4_from_lh_fixallmp = make_contour_original_pars_from_lh( 3, 4, true ) ;
////  h_cont_op_p3_vs_p4_from_lh_fixallmp -> SetContour(nb) ;
////  h_cont_op_p3_vs_p4_from_lh_fixallmp -> Draw("colz") ;


////  TH2F* h_cont_rp_p3_vs_p4_from_lh = make_contour_rotated_pars_from_lh( 3, 4 ) ;
////  h_cont_rp_p3_vs_p4_from_lh -> SetContour(nb) ;
////  h_cont_rp_p3_vs_p4_from_lh -> Draw("colz") ;

     //-------


////  TH2F* h_cont_op_p4_vs_p2_from_lh = make_contour_original_pars_from_lh( 2, 4, false ) ;
////  h_cont_op_p4_vs_p2_from_lh -> SetContour(nb) ;
////  h_cont_op_p4_vs_p2_from_lh -> Draw("colz") ;


////  TH2F* h_cont_op_p4_vs_p2_from_lh_fixallmp = make_contour_original_pars_from_lh( 2, 4, true ) ;
////  h_cont_op_p4_vs_p2_from_lh_fixallmp -> SetContour(nb) ;
////  h_cont_op_p4_vs_p2_from_lh_fixallmp -> Draw("colz") ;


////  TH2F* h_cont_rp_p4_vs_p2_from_lh = make_contour_rotated_pars_from_lh( 2, 4 ) ;
////  h_cont_rp_p4_vs_p2_from_lh -> SetContour(nb) ;
////  h_cont_rp_p4_vs_p2_from_lh -> Draw("colz") ;

     //-------


      saveHist("outputfiles/plots.root","h*") ;


   } // gen_nuisance_parameter_matrix1


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

   TH2F* make_contour_original_pars_from_lh( int pix, int piy, bool fix_all_model_pars ) {

      char hname[100] ;
      char htitle[1000] ;

      if ( fix_all_model_pars ) {
         sprintf( hname, "h_cont_op_p%d_vs_p%d_from_lh_fixallmp", piy, pix ) ;
         sprintf( htitle, "nll contour, %s vs %s, all model pars fixed", par_name[piy], par_name[pix] ) ;
      } else {
         sprintf( hname, "h_cont_op_p%d_vs_p%d_from_lh", piy, pix ) ;
         sprintf( htitle, "nll contour, %s vs %s", par_name[piy], par_name[pix] ) ;
      }

      float xl = par_val[pix] - 2.5*par_err[pix] ;
      float xh = par_val[pix] + 2.5*par_err[pix] ;
      float yl = par_val[piy] - 2.5*par_err[piy] ;
      float yh = par_val[piy] + 2.5*par_err[piy] ;
      if ( xl < 0 ) xl = 0. ;
      if ( yl < 0 ) yl = 0. ;

      //int ncp(40) ;
      //int ncp(5) ;
      //int ncp(10) ;
      int ncp(20) ;
      printf("\n\n make_contour_original_pars_from_lh: creating %s , %s\n", hname, htitle ) ;
      TH2F* hp = new TH2F( hname, htitle, ncp, xl, xh, ncp, yl, yh ) ;

      printf(" Scanning on x axis %s  (%6.4f +/- %6.4f) from %6.4f to %6.4f\n", par_name[pix], par_val[pix], par_err[pix], xl, xh ) ;
      printf(" Scanning on y axis %s  (%6.4f +/- %6.4f) from %6.4f to %6.4f\n", par_name[piy], par_val[piy], par_err[piy], yl, yh ) ; fflush( stdout ) ;

      RooRealVar* rv_par[10] ;
      for ( int pi=0; pi<npars; pi++ ) {
         rv_par[pi] = the_ws -> var( par_name[pi] ) ;
         if ( rv_par[pi] == 0x0 ) { printf("\n\n *** can't find %s in ws.\n\n", par_name[pi] ) ; gSystem->Exit(-1) ; }
      }

      double min_nll(0.) ;
      {
         RooFitResult* rfr = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(0) ) ;
         min_nll = rfr->minNll() ;
         delete rfr ;
      }

      float par_startvals[20] ;
      printf("\n\n Starting values at best fit:\n") ;
      for ( int pi=0; pi<npars; pi++ ) {
         par_startvals[pi] = rv_par[pi] -> getVal() ;
         printf("  %12s : %9.4f\n", par_name[pi], par_startvals[pi] ) ;
      } // pi
      printf("\n\n") ;


      for ( int xi=0; xi<ncp; xi++ ) {

         float px = xl + (xh-xl)*(xi+0.5)/ncp ;

         for ( int yi=0; yi<ncp; yi++ ) {

            float py = yl + (yh-yl)*(yi+0.5)/ncp ;

            printf(" grid point %s = %6.4f  %s = %6.4f\n", par_name[pix], px, par_name[piy], py ) ;


            rv_par[pix] -> setVal( px ) ;
            rv_par[piy] -> setVal( py ) ;

           //-------
      //    double nll = likelihood -> getLogVal() ;
           //-------
            if ( fix_all_model_pars ) {
               for ( int pi=0; pi<npars; pi++ ) {
                  RooRealVar* rrv = the_ws -> var( par_name[pi] ) ;
                  rrv -> setConstant( kTRUE ) ;
               } // pi
            } else {
               rv_par[pix] -> setConstant( kTRUE ) ;
               rv_par[piy] -> setConstant( kTRUE ) ;
            }
            RooFitResult* rfr = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(0) ) ;
            double nll = rfr->minNll() ;
            delete rfr ;
           //-------


            //printf("  %s = %6.4f ,  %s = %6.4f ,  chi2 = %f\n", par_name[pix], px, par_name[piy], py, nll ) ;

            hp -> SetBinContent( xi+1, yi+1, 2.*(nll-min_nll) ) ;

         } // yi

      } // xi

      for ( int pi=0; pi<npars; pi++ ) {
         rv_par[pi] -> setVal( par_startvals[pi] ) ;
         rv_par[pi] -> setConstant( kFALSE ) ;
      } // pi

      return hp ;


   } // make_contour_original_pars_from_lh

//===============================================================================

   TH2F* make_contour_rotated_pars_from_lh( int pix, int piy ) {

      char hname[100] ;
      char htitle[1000] ;

      sprintf( hname, "h_cont_rp_p%d_vs_p%d_from_lh", piy, pix ) ;
      sprintf( htitle, "nll contour, rotated %s vs %s", par_name[piy], par_name[pix] ) ;

      float xprime_err = sqrt(reordered_eigen_vals[pix]) ;
      float yprime_err = sqrt(reordered_eigen_vals[piy]) ;

      float xl = - 2.5*xprime_err ;
      float xh = + 2.5*xprime_err ;
      float yl = - 2.5*yprime_err ;
      float yh = + 2.5*yprime_err ;


      //int ncp(40) ;
      //int ncp(5) ;
      //int ncp(10) ;
      int ncp(20) ;
      printf("\n\n make_contour_rotated_pars_from_lh: creating %s , %s\n", hname, htitle ) ;
      TH2F* hp = new TH2F( hname, htitle, ncp, xl, xh, ncp, yl, yh ) ;

      printf(" Scanning on x' axis %s  (+/- %6.4f)\n", par_name[pix], sqrt(reordered_eigen_vals[pix]) ) ;
      printf(" Scanning on y' axis %s  (+/- %6.4f)\n", par_name[piy], sqrt(reordered_eigen_vals[piy]) ) ; fflush( stdout ) ;

      RooRealVar* rv_par[20] ;
      for ( int pi=0; pi<npars; pi++ ) {
         rv_par[pi] = the_ws -> var( par_name[pi] ) ;
         if ( rv_par[pi] == 0x0 ) { printf("\n\n *** can't find %s in ws.\n\n", par_name[pi] ) ; gSystem->Exit(-1) ; }
      }

      double min_nll(0.) ;
      {
         RooFitResult* rfr = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(0) ) ;
         min_nll = rfr->minNll() ;
         delete rfr ;
      }

      float par_startvals[20] ;
      printf("\n\n Starting values at best fit:\n") ;
      for ( int pi=0; pi<npars; pi++ ) {
         par_startvals[pi] = rv_par[pi] -> getVal() ;
         printf("  %12s : %9.4f\n", par_name[pi], par_startvals[pi] ) ;
      } // pi
      printf("\n\n") ;


      for ( int xi=0; xi<ncp; xi++ ) {

         float dxprime = xl + (xh-xl)*(xi+0.5)/ncp ;

         for ( int yi=0; yi<ncp; yi++ ) {

            float dyprime = yl + (yh-yl)*(yi+0.5)/ncp ;
            printf("   dxprime = %6.4f,  dyprime = %6.4f\n", dxprime, dyprime ) ;
      ///   printf("      Unrotated par deltas:  " ) ;
      ///   for ( int pi=0; pi<6; pi++ ) {
      ///      float new_val = dxprime * reordered_reoriented_eigen_vector_matrix[pi][pix] + dyprime * reordered_reoriented_eigen_vector_matrix[pi][piy] ;
      ///      printf("  %9.5f  ", new_val ) ;
      ///   }
            printf("      Unrotated pars:  " ) ;
            for ( int pi=0; pi<npars; pi++ ) {
               float new_val = par_val[pi] + dxprime * reordered_reoriented_eigen_vector_matrix[pi][pix] + dyprime * reordered_reoriented_eigen_vector_matrix[pi][piy] ;
               printf("  %9.5f  ", new_val ) ;
               if ( new_val < 0 ) new_val = 0. ;
               rv_par[pi] -> setVal( new_val ) ;
               rv_par[pi] -> setConstant( kTRUE ) ;
            }
            printf("\n") ;

           //-------
      //    double nll = likelihood -> getLogVal() ;
           //-------
            RooFitResult* rfr = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(0) ) ;
            double nll = rfr->minNll() ;
            delete rfr ;
           //-------

            hp -> SetBinContent( xi+1, yi+1, 2.*(nll-min_nll) ) ;

         } // yi

      } // xi

      for ( int pi=0; pi<npars; pi++ ) {
         rv_par[pi] -> setVal( par_startvals[pi] ) ;
         rv_par[pi] -> setConstant( kFALSE ) ;
      } // pi

      return hp ;


   } // make_contour_rotated_pars_from_lh

//===============================================================================









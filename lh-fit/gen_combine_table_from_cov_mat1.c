
#include "TH2F.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "TMatrixT.h"
#include "TMatrixD.h"
#include "TVectorD.h"
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

   int npars(6) ;
   float par_val[10] ;
   float par_err[10] ;
   char par_name[10][100] ;

   int ht_par[10] ;
   int njet_par[10] ;

   TH1F* h_ratio ;

   RooWorkspace* the_ws ;
   RooAbsPdf* likelihood ;
   RooDataSet* rds ;

   TMatrixD reordered_reoriented_eigen_vector_matrix( 6, 6 ) ;
   TVectorD reordered_eigen_vals( 6 ) ;

  //---------

   TH1F* get_hist( const char* hname ) ;
   TH2F* make_contour_original_pars_from_chi2( int pix, int piy ) ;
   TH2F* make_contour_original_pars_from_lh( int pix, int piy, bool fix_all_model_pars=true ) ;
   TH2F* make_contour_rotated_pars_from_lh( int pix, int piy ) ;

  //---------

   void gen_combine_table_from_cov_mat( const char* infile = "outputfiles/lhfit-results-ws-lhfit-test/kqcd-parameter-fit-covmat.tex",
                                        const char* datarootfile = "../outputfiles/modelfit-input-data.root",
                                        const char* wsfile = "outputfiles/ws-lhfit-test.root" ) {

      gDirectory -> Delete( "h*" ) ;

      loadHist( datarootfile ) ;

      h_ratio = get_hist( "h_ratio" ) ;


      ifstream ifs ;
      ifs.open( infile ) ;
      if ( !ifs.good() ) { printf("\n\n *** Bad input file: %s\n\n", infile ) ; return ; }


      TString line ;

      printf("\n\n Reading parameters from file %s\n", infile ) ;

      for ( int pi=0; pi<npars; pi++ ) {
         line.ReadLine( ifs ) ;
         sscanf( line.Data(), "%s %f +/- %f", par_name[pi], &(par_val[pi]), &(par_err[pi]) ) ;
         if ( strcmp( par_name[pi], "Kqcd_ht1" ) == 0 ) { ht_par[1] = pi ; }
         if ( strcmp( par_name[pi], "Kqcd_ht2" ) == 0 ) { ht_par[2] = pi ; }
         if ( strcmp( par_name[pi], "Kqcd_ht3" ) == 0 ) { ht_par[3] = pi ; }
         if ( strcmp( par_name[pi], "Sqcd_njet2" ) == 0 ) { njet_par[2] = pi ; }
         if ( strcmp( par_name[pi], "Sqcd_njet3" ) == 0 ) { njet_par[3] = pi ; }
         if ( strcmp( par_name[pi], "Sqcd_njet4" ) == 0 ) { njet_par[4] = pi ; }
      } // pi
      par_val[npars] = 1. ;
      par_err[npars] = 0. ;
      sprintf( par_name[npars], "Sqcd_njet1" ) ;
      njet_par[1] = npars ;

      printf("   ht1 par index %d\n", ht_par[1] ) ;
      printf("   ht2 par index %d\n", ht_par[2] ) ;
      printf("   ht3 par index %d\n", ht_par[3] ) ;
      printf(" njet2 par index %d\n", njet_par[2] ) ;
      printf(" njet3 par index %d\n", njet_par[3] ) ;
      printf(" njet4 par index %d\n", njet_par[4] ) ;
      printf(" njet1 par index %d\n", njet_par[1] ) ;


      line.ReadLine( ifs ) ;

      float cov_mat[10][10] ;
      for ( int ri=0; ri<npars; ri++ ) {
         line.ReadLine( ifs ) ;
         char label[100], pname[100] ;
         sscanf( line.Data(), "%s %s %f %f %f %f %f %f", label, pname,
           &(cov_mat[ri][0]), &(cov_mat[ri][1]), &(cov_mat[ri][2]), &(cov_mat[ri][3]), &(cov_mat[ri][4]), &(cov_mat[ri][5]) ) ;
         if ( strcmp( pname, par_name[ri] ) != 0 ) { printf("\n\n *** parameter mismatch: %s %s\n\n", pname, par_name[ri] ) ; return ; }
      } // ri

      printf(" Done reading parameters from file.\n") ;



      printf("\n\n  Fit parameters:\n") ;
      for ( int pi=0; pi<npars; pi++ ) {
         printf("  %d : %12s : %12.5f +/- %12.5f\n", pi, par_name[pi], par_val[pi], par_err[pi] ) ;
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
            h_fit_correlation_matrix -> GetXaxis() -> SetBinLabel( ci+1, par_name[ci] ) ;
            h_fit_correlation_matrix -> GetYaxis() -> SetBinLabel( npars-ri, par_name[ri] ) ;
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

      int eigen_vector_index[10] ;
      for ( int i=0; i<npars; i++ ) { eigen_vector_index[i] = -1 ; }

      printf("\n\n  Check of eigen vectors:\n") ;
      for ( int i=0; i<npars; i++ ) {
         printf("  Eigen value %d : %12.8f\n", i, eigen_vals[i] ) ;
         TMatrixT<double> eigen_vector_col(6,1) ;
         TMatrixT<double> eigen_vector_row(1,6) ;
         double largest_component_val(0.) ;
         int largest_component_ind(-1) ;
         printf("      Eigen vector elements : ") ;
         for ( int j=0; j<6; j++ ) {
            eigen_vector_col(j,0) = eigen_vector_matrix[j][i] ;
            eigen_vector_row(0,j) = eigen_vector_matrix[j][i] ;
            if ( fabs( eigen_vector_matrix[j][i] ) > largest_component_val ) {
               largest_component_val = fabs( eigen_vector_matrix[j][i] ) ;
               largest_component_ind = j ;
            }
            printf( "  %12.8f  ", eigen_vector_col(j,0) ) ;
         }
         printf("\n") ;
         printf("   Largest component is %s\n", par_name[largest_component_ind] ) ;
         eigen_vector_index[largest_component_ind] = i ;
         TMatrixT<double> cm_times_ev(6,1) ;
         cm_times_ev.Mult( cov_mat_tm, eigen_vector_col ) ;
         printf("      CM times eigen vector : ") ;
         for ( int j=0; j<6; j++ ) {
            printf( "  %12.8f  ", cm_times_ev(j,0) ) ;
         }
         printf("\n") ;
         printf("      after div  by EV      : ") ;
         for ( int j=0; j<6; j++ ) {
            printf( "  %12.8f  ", cm_times_ev(j,0) / eigen_vals[i] ) ;
         }
         printf("\n") ;
         TMatrixT<double> ev_times_cm_times_ev(1,1) ;
         ev_times_cm_times_ev.Mult( eigen_vector_row, cm_times_ev ) ;
         printf("  EVec times CM times EVec (should EVal) :  %12.8f\n", ev_times_cm_times_ev(0,0) ) ;

         printf("\n\n") ;
      } //


      printf("\n\n") ;
      for ( int pi=0; pi<6; pi++ ) {
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
            sprintf( xlabel, "EV-%s", par_name[j] ) ;
            h_eigenvec_matrix -> SetBinContent( j+1, npars-i, reordered_reoriented_eigen_vector_matrix[i][j] ) ;
            h_eigenvec_matrix -> GetXaxis() -> SetBinLabel( j+1, xlabel ) ;
            h_eigenvec_matrix -> GetYaxis() -> SetBinLabel( npars-i, par_name[i] ) ;
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
      for ( int nji=2; nji<=4; nji++ ) {
         for ( int hti=1; hti<=3; hti++ ) {
            printf("   Njet%d-HT%d  ", nji, hti ) ;
            float sum_err2(0.) ;
            for ( int rfpi=0; rfpi<npars; rfpi++ ) {
               if ( !( rfpi == ht_par[hti] || rfpi == njet_par[nji] ) ) continue ;
               for ( int cfpi=0; cfpi<npars; cfpi++ ) {
                  if ( !( cfpi == ht_par[hti] || cfpi == njet_par[nji] ) ) continue ;
                  sum_err2 += cov_mat[rfpi][cfpi] / ( par_val[rfpi] * par_val[cfpi] ) ;
               } // cfpi
            } // rfpi
            printf("     (%6.4f)\n", sqrt( sum_err2 ) ) ;
         } // hti
      } // nji




      printf("\n\n  Calculation of total relative error and individual contributions in rotated basis.\n") ;
      printf("                  HT1        HT2       HT3       Nj2       Nj3       Nj4\n") ;

      for ( int nji=2; nji<=4; nji++ ) {
         for ( int hti=1; hti<=3; hti++ ) {
            printf("   Njet%d-HT%d  ", nji, hti ) ;
            double par_ht   = par_val[ ht_par[hti] ] ;
            double par_njet = par_val[ njet_par[nji] ] ;
            if ( par_ht == 0 ) { printf("\n\n *** par val zero.\n\n") ; return ; }
            if ( par_njet == 0 ) { printf("\n\n *** par val zero.\n\n") ; return ; }
            float sum_err2(0.) ;
            for ( int fpi=0; fpi<npars; fpi++ ) {
               int evi = eigen_vector_index[fpi] ;
               double eigen_val = eigen_vals[ evi ] ;
               double rot_mat_over_par_val_sum = rotation_matrix( evi, ht_par[hti] ) / par_ht + rotation_matrix( evi, njet_par[nji] ) / par_njet ;
        ///    printf("  [  sum = %12.8f + %12.8f ] ",
        ///      sqrt(eigen_val) * rotation_matrix( evi, ht_par[hti] ) / par_ht,
        ///      sqrt(eigen_val) * rotation_matrix( evi, njet_par[nji] ) / par_njet ) ;
               double rel_err = sqrt( eigen_val ) * fabs(rot_mat_over_par_val_sum) ;
               printf( "  %6.4f  ", rel_err ) ;
               sum_err2 += rel_err * rel_err ;
            } // fpi
            printf("     (%6.4f)\n", sqrt( sum_err2 ) ) ;
         } // hti
      } // nji





      TH2F* h_rel_err_table = new TH2F( "h_rel_err_table",
           "Relative error table, orthogonal parameters basis", npars+2, 0.5, npars+2+0.5,  12, 0.5, 12+0.5 ) ;

      printf("\n\n  Calculation of total relative error and individual contributions in rotated basis, using reordered reoriented rotation matrix.\n") ;
      printf("                  HT1        HT2       HT3       Nj2       Nj3       Nj4\n") ;

      int rowi(0) ;
      for ( int nji=1; nji<=4; nji++ ) {
         for ( int hti=1; hti<=3; hti++ ) {
            rowi++ ;
            char row_label[100] ;
            sprintf( row_label, "Njet%d-HT%d", nji, hti ) ;
            printf("   %s  ", row_label ) ;
            double par_ht   = par_val[ ht_par[hti] ] ;
            double par_njet = par_val[ njet_par[nji] ] ;
            if ( par_ht == 0 ) { printf("\n\n *** par val zero.\n\n") ; return ; }
            if ( par_njet == 0 ) { printf("\n\n *** par val zero.\n\n") ; return ; }
            float sum_err2(0.) ;
            for ( int fpi=0; fpi<npars; fpi++ ) {
               double eigen_val = reordered_eigen_vals[ fpi ] ;
               double rot_mat_over_par_val_sum(0.) ;
               if ( nji>1 ) {
                  rot_mat_over_par_val_sum = reordered_reoriented_rotation_matrix( fpi, ht_par[hti] ) / par_ht + reordered_reoriented_rotation_matrix( fpi, njet_par[nji] ) / par_njet ;
               } else {
                  rot_mat_over_par_val_sum = reordered_reoriented_rotation_matrix( fpi, ht_par[hti] ) / par_ht  ;
               }
               double rel_err = sqrt( eigen_val ) * rot_mat_over_par_val_sum ;
               printf( " %7.4f  ", rel_err ) ;
               h_rel_err_table -> SetBinContent( fpi+1, 13-rowi, rel_err ) ;
               h_rel_err_table -> GetYaxis() -> SetBinLabel( 13-rowi, row_label ) ;
               h_rel_err_table -> GetXaxis() -> SetBinLabel( fpi+1, par_name[fpi] ) ;
               sum_err2 += rel_err * rel_err ;
            } // fpi
            printf("     (%7.4f)\n", sqrt( sum_err2 ) ) ;
            h_rel_err_table -> SetBinContent( npars+2, 13-rowi, sqrt( sum_err2 ) ) ;
            h_rel_err_table -> GetXaxis() -> SetBinLabel( npars+2, "Total error" ) ;
         } // hti
      } // nji






      TH2F* h_rel_err_table_simple_ignore = new TH2F( "h_rel_err_table_simple_ignore",
              "Relative error table, original parameters basis, ignore off-diagonal cov",
              npars+2, 0.5, npars+2+0.5,  12, 0.5, 12+0.5 ) ;

      printf("\n\n") ;

      printf("\n\n  Calculation ignoring off-diagonal elements of the covariance matrix.\n") ;
      printf("                  HT1        HT2       HT3       Nj2       Nj3       Nj4\n") ;

      rowi = 0 ;
      for ( int nji=1; nji<=4; nji++ ) {
         for ( int hti=1; hti<=3; hti++ ) {
            rowi++ ;
            printf("   Njet%d-HT%d  ", nji, hti ) ;
            char row_label[100] ;
            sprintf( row_label, "Njet%d-HT%d", nji, hti ) ;
            float sum_err2(0.) ;
            for ( int fpi=0; fpi<npars; fpi++ ) {
               double rel_err(0.) ;
               if ( fpi == ht_par[hti] ) { rel_err = par_err[fpi]/par_val[fpi] ; }
               if ( fpi == njet_par[nji] ) { rel_err = par_err[fpi]/par_val[fpi] ; }
               h_rel_err_table_simple_ignore -> SetBinContent( fpi+1, 13-rowi, rel_err ) ;
               h_rel_err_table_simple_ignore -> GetYaxis() -> SetBinLabel( 13-rowi, row_label ) ;
               h_rel_err_table_simple_ignore -> GetXaxis() -> SetBinLabel( fpi+1, par_name[fpi] ) ;
               printf( "  %6.4f  ", rel_err ) ;
               sum_err2 += rel_err * rel_err ;
            } // fpi
            printf("     (%6.4f)\n", sqrt( sum_err2 ) ) ;
            h_rel_err_table_simple_ignore -> SetBinContent( npars+2, 13-rowi, sqrt( sum_err2 ) ) ;
            h_rel_err_table_simple_ignore -> GetXaxis() -> SetBinLabel( npars+2, "Total error" ) ;
         } // hti
      } // nji

      printf("\n\n") ;


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

    //---
      TCanvas* can3 = new TCanvas("can3","Fit correlation coefficient matrix",900,800 ) ;
      can3 -> SetWindowPosition( cx, cy ) ;

      h_fit_correlation_matrix -> SetContour(nb) ;
      h_fit_correlation_matrix -> SetMinimum( -1.0 ) ;
      h_fit_correlation_matrix -> SetMaximum(  1.0 ) ;

      h_fit_correlation_matrix -> Draw("colz X+") ;
      h_fit_correlation_matrix -> Draw("text same") ;

    //---
      TCanvas* can2 = new TCanvas("can2","Combine table, original parameters, ignoring correlations",900,800 ) ;
      cx += 50 ; cy += 50 ;
      can2 -> SetWindowPosition( cx, cy ) ;

      h_rel_err_table_simple_ignore -> SetContour(nb) ;
      h_rel_err_table_simple_ignore -> SetMinimum( -0.7 ) ;
      h_rel_err_table_simple_ignore -> SetMaximum(  0.7 ) ;

      h_rel_err_table_simple_ignore -> Draw("colz") ;
      h_rel_err_table_simple_ignore -> Draw("text same") ;

    //---
      TCanvas* can1 = new TCanvas("can1","Combine table, rotated parameters",900,800 ) ;
      cx += 50 ; cy += 50 ;
      can1 -> SetWindowPosition( cx, cy ) ;

      h_rel_err_table -> SetContour(nb) ;
      h_rel_err_table -> SetMinimum( -0.7 ) ;
      h_rel_err_table -> SetMaximum(  0.7 ) ;

      h_rel_err_table -> Draw("colz") ;
      h_rel_err_table -> Draw("text same") ;


    //---
      TCanvas* can4 = new TCanvas("can4","Cov mat Eigen vectors",900,800 ) ;
      cx += 50 ; cy += 50 ;
      can4 -> SetWindowPosition( cx, cy ) ;


      h_eigenvec_matrix -> SetContour(nb) ;
      h_eigenvec_matrix -> SetMinimum( -0.7 ) ;
      h_eigenvec_matrix -> SetMaximum(  0.7 ) ;

      h_eigenvec_matrix -> Draw("colz") ;
      h_eigenvec_matrix -> Draw("text same") ;









     //+++++++++++

      printf("  h_ratio pointer %p\n", h_ratio ) ; fflush(stdout) ;
      TH2F* h_cont_op_p0_vs_p1_from_chi2 = make_contour_original_pars_from_chi2( 0, 1 ) ;
      TH2F* h_cont_op_p1_vs_p2_from_chi2 = make_contour_original_pars_from_chi2( 1, 2 ) ;
      TH2F* h_cont_op_p2_vs_p3_from_chi2 = make_contour_original_pars_from_chi2( 2, 3 ) ;
      TH2F* h_cont_op_p3_vs_p4_from_chi2 = make_contour_original_pars_from_chi2( 3, 4 ) ;


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

///// TH2F* h_cont_op_p0_vs_p1_from_lh = make_contour_original_pars_from_lh( 0, 1, false ) ;
///// h_cont_op_p0_vs_p1_from_lh -> SetContour(nb) ;
///// h_cont_op_p0_vs_p1_from_lh -> Draw("colz") ;


///// TH2F* h_cont_op_p0_vs_p1_from_lh_fixallmp = make_contour_original_pars_from_lh( 0, 1, true ) ;
///// h_cont_op_p0_vs_p1_from_lh_fixallmp -> SetContour(nb) ;
///// h_cont_op_p0_vs_p1_from_lh_fixallmp -> Draw("colz") ;


///// TH2F* h_cont_rp_p0_vs_p1_from_lh = make_contour_rotated_pars_from_lh( 0, 1 ) ;
///// h_cont_rp_p0_vs_p1_from_lh -> SetContour(nb) ;
///// h_cont_rp_p0_vs_p1_from_lh -> Draw("colz") ;


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


      TH2F* h_cont_op_p4_vs_p2_from_lh = make_contour_original_pars_from_lh( 2, 4, false ) ;
      h_cont_op_p4_vs_p2_from_lh -> SetContour(nb) ;
      h_cont_op_p4_vs_p2_from_lh -> Draw("colz") ;


      TH2F* h_cont_op_p4_vs_p2_from_lh_fixallmp = make_contour_original_pars_from_lh( 2, 4, true ) ;
      h_cont_op_p4_vs_p2_from_lh_fixallmp -> SetContour(nb) ;
      h_cont_op_p4_vs_p2_from_lh_fixallmp -> Draw("colz") ;


      TH2F* h_cont_rp_p4_vs_p2_from_lh = make_contour_rotated_pars_from_lh( 2, 4 ) ;
      h_cont_rp_p4_vs_p2_from_lh -> SetContour(nb) ;
      h_cont_rp_p4_vs_p2_from_lh -> Draw("colz") ;

     //-------


      saveHist("outputfiles/plots.root","h*") ;


   } // gen_combine_table_from_cov_mat


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

   TH2F* make_contour_original_pars_from_chi2( int pix, int piy ) {

      char hname[100] ;
      char htitle[1000] ;

      sprintf( hname, "h_cont_op_p%d_vs_p%d_from_chi2", piy, pix ) ;
      sprintf( htitle, "chi2 contour, %s vs %s", par_name[piy], par_name[pix] ) ;

      float xl = par_val[pix] - 2.5*par_err[pix] ;
      float xh = par_val[pix] + 2.5*par_err[pix] ;
      float yl = par_val[piy] - 2.5*par_err[piy] ;
      float yh = par_val[piy] + 2.5*par_err[piy] ;
      if ( xl < 0 ) xl = 0. ;
      if ( yl < 0 ) yl = 0. ;

      int ncp(40) ;
      //int ncp(5) ;
      printf("\n\n make_contour_original_pars_from_chi2: creating %s , %s\n", hname, htitle ) ;
      TH2F* hp = new TH2F( hname, htitle, ncp, xl, xh, ncp, yl, yh ) ;

      printf(" Scanning on x axis %s  (%6.4f +/- %6.4f) from %6.4f to %6.4f\n", par_name[pix], par_val[pix], par_err[pix], xl, xh ) ;
      printf(" Scanning on y axis %s  (%6.4f +/- %6.4f) from %6.4f to %6.4f\n", par_name[piy], par_val[piy], par_err[piy], yl, yh ) ; fflush( stdout ) ;

      for ( int xi=0; xi<ncp; xi++ ) {

         float px = xl + (xh-xl)*(xi+0.5)/ncp ;

         for ( int yi=0; yi<ncp; yi++ ) {

            float py = yl + (yh-yl)*(yi+0.5)/ncp ;

            ////////printf(" grid point %s = %6.4f  %s = %6.4f\n", par_name[pix], px, par_name[piy], py ) ;

            float chi2(0.) ;
            int ratio_hist_bin(0) ;
            for ( int hti=1; hti<=3; hti++ ) {
               for ( int nji=1; nji<=4; nji++ ) {

                  ratio_hist_bin ++ ;

                  if ( hti==1 && nji>2 ) continue ;

                  float ratio_val = h_ratio -> GetBinContent( ratio_hist_bin ) ;
                  float ratio_err = h_ratio -> GetBinError( ratio_hist_bin ) ;

                  int pi_ht = ht_par[hti] ;
                  int pi_njet = njet_par[nji] ;

                  float kht   = par_val[pi_ht] ;
                  float snjet = par_val[pi_njet] ;

                  if ( pix == pi_ht   ) { kht   = px ; }
                  if ( pix == pi_njet ) { snjet = px ; }
                  if ( piy == pi_ht   ) { kht   = py ; }
                  if ( piy == pi_njet ) { snjet = py ; }

                  float model_val = kht * snjet ;

                  ////////printf("   --- pi_ht = %d, pi_njet = %d\n", pi_ht, pi_njet ) ;
                  //////printf("   --- ht%d njet%d : hist bin %2d (%s) : data = %6.4f +/- %6.4f ,  model = %6.4f * %6.4f = %6.4f,  chi2 contribution %9.4f\n",
                      //////hti, nji, ratio_hist_bin, h_ratio->GetXaxis()->GetBinLabel( ratio_hist_bin ),
                      //////ratio_val, ratio_err, kht, snjet, model_val, pow( (ratio_val - model_val)/ratio_err, 2. ) ) ; fflush( stdout ) ;

                  if ( ratio_err > 0 ) {
                     chi2 += pow( (ratio_val - model_val)/ratio_err, 2. ) ;
                  }


               } // nji
            } // hti

            //////printf("  %s = %6.4f ,  %s = %6.4f ,  chi2 = %8.4f\n", par_name[pix], px, par_name[piy], py, chi2 ) ;

            hp -> SetBinContent( xi+1, yi+1, chi2 ) ;

         } // yi

      } // xi

      return hp ;

   } // make_contour_original_pars_from_chi2

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
      for ( int pi=0; pi<6; pi++ ) {
         rv_par[pi] = the_ws -> var( par_name[pi] ) ;
         if ( rv_par[pi] == 0x0 ) { printf("\n\n *** can't find %s in ws.\n\n", par_name[pi] ) ; gSystem->Exit(-1) ; }
      }

      double min_nll(0.) ;
      {
         RooFitResult* rfr = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(0) ) ;
         min_nll = rfr->minNll() ;
         delete rfr ;
      }

      float par_startvals[10] ;
      printf("\n\n Starting values at best fit:\n") ;
      for ( int pi=0; pi<6; pi++ ) {
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
               for ( int pi=0; pi<6; pi++ ) {
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

      for ( int pi=0; pi<6; pi++ ) {
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

      RooRealVar* rv_par[10] ;
      for ( int pi=0; pi<6; pi++ ) {
         rv_par[pi] = the_ws -> var( par_name[pi] ) ;
         if ( rv_par[pi] == 0x0 ) { printf("\n\n *** can't find %s in ws.\n\n", par_name[pi] ) ; gSystem->Exit(-1) ; }
      }

      double min_nll(0.) ;
      {
         RooFitResult* rfr = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(0) ) ;
         min_nll = rfr->minNll() ;
         delete rfr ;
      }

      float par_startvals[10] ;
      printf("\n\n Starting values at best fit:\n") ;
      for ( int pi=0; pi<6; pi++ ) {
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
            for ( int pi=0; pi<6; pi++ ) {
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

      for ( int pi=0; pi<6; pi++ ) {
         rv_par[pi] -> setVal( par_startvals[pi] ) ;
         rv_par[pi] -> setConstant( kFALSE ) ;
      } // pi

      return hp ;


   } // make_contour_rotated_pars_from_lh

//===============================================================================









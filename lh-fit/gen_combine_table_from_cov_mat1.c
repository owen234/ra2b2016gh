
#include "TH2F.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "TMatrixT.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <fstream>


   void gen_combine_table_from_cov_mat( const char* infile = "outputfiles/lhfit-results-ws-lhfit-test/kqcd-parameter-fit-covmat.tex" ) {

      ifstream ifs ;
      ifs.open( infile ) ;
      if ( !ifs.good() ) { printf("\n\n *** Bad input file: %s\n\n", infile ) ; return ; }

      int npars(6) ;
      float par_val[10] ;
      float par_err[10] ;
      char par_name[10][100] ;

      int ht_par[10] ;
      int njet_par[10] ;

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
      printf("\n\n   Fit correlation matrix:\n") ;
      for ( int ri=0; ri<npars; ri++ ) {
         for ( int ci=0; ci<npars; ci++ ) {
            printf("  %12.8f  ", cov_mat[ri][ci] / sqrt( cov_mat[ri][ri] * cov_mat[ci][ci] ) ) ;
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






      TMatrixD reordered_reoriented_eigen_vector_matrix( npars, npars ) ;
      TVectorD reordered_eigen_vals( npars ) ;

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
      printf("\n") ;
      printf("\n\n Reordered, reoriented, eigen vector matrix:\n") ;
      for ( int i=0; i<npars; i++ ) {
         for ( int j=0; j<npars; j++ ) {
            printf("  %12.8f  ", reordered_reoriented_eigen_vector_matrix[i][j] ) ;
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





      TH2F* h_rel_err_table = new TH2F( "h_rel_err_table", "Relative error table, orthogonal parameters basis", npars+2, 0.5, npars+2+0.5,  9, 0.5, 9+0.5 ) ;

      printf("\n\n  Calculation of total relative error and individual contributions in rotated basis, using reordered reoriented rotation matrix.\n") ;
      printf("                  HT1        HT2       HT3       Nj2       Nj3       Nj4\n") ;

      int rowi(0) ;
      for ( int nji=2; nji<=4; nji++ ) {
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
               double rot_mat_over_par_val_sum = reordered_reoriented_rotation_matrix( fpi, ht_par[hti] ) / par_ht + reordered_reoriented_rotation_matrix( fpi, njet_par[nji] ) / par_njet ;
               double rel_err = sqrt( eigen_val ) * rot_mat_over_par_val_sum ;
               printf( " %7.4f  ", rel_err ) ;
               h_rel_err_table -> SetBinContent( fpi+1, 10-rowi, rel_err ) ;
               h_rel_err_table -> GetYaxis() -> SetBinLabel( 10-rowi, row_label ) ;
               h_rel_err_table -> GetXaxis() -> SetBinLabel( fpi+1, par_name[fpi] ) ;
               sum_err2 += rel_err * rel_err ;
            } // fpi
            printf("     (%7.4f)\n", sqrt( sum_err2 ) ) ;
            h_rel_err_table -> SetBinContent( npars+2, 10-rowi, sqrt( sum_err2 ) ) ;
            h_rel_err_table -> GetXaxis() -> SetBinLabel( npars+2, "Total error" ) ;
         } // hti
      } // nji






      TH2F* h_rel_err_table_simple_ignore = new TH2F( "h_rel_err_table_simple_ignore", "Relative error table, original parameters basis, ignore off-diagonal cov", npars+2, 0.5, npars+2+0.5,  9, 0.5, 9+0.5 ) ;

      printf("\n\n") ;

      printf("\n\n  Calculation ignoring off-diagonal elements of the covariance matrix.\n") ;
      printf("                  HT1        HT2       HT3       Nj2       Nj3       Nj4\n") ;

      rowi = 0 ;
      for ( int nji=2; nji<=4; nji++ ) {
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
               h_rel_err_table_simple_ignore -> SetBinContent( fpi+1, 10-rowi, rel_err ) ;
               h_rel_err_table_simple_ignore -> GetYaxis() -> SetBinLabel( 10-rowi, row_label ) ;
               h_rel_err_table_simple_ignore -> GetXaxis() -> SetBinLabel( fpi+1, par_name[fpi] ) ;
               printf( "  %6.4f  ", rel_err ) ;
               sum_err2 += rel_err * rel_err ;
            } // fpi
            printf("     (%6.4f)\n", sqrt( sum_err2 ) ) ;
            h_rel_err_table_simple_ignore -> SetBinContent( npars+2, 10-rowi, sqrt( sum_err2 ) ) ;
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

    //---
      TCanvas* can1 = new TCanvas("can1","",900,800 ) ;

      h_rel_err_table -> SetContour(nb) ;
      h_rel_err_table -> SetMinimum( -0.7 ) ;
      h_rel_err_table -> SetMaximum(  0.7 ) ;

      h_rel_err_table -> Draw("colz") ;
      h_rel_err_table -> Draw("text same") ;

    //---
      TCanvas* can2 = new TCanvas("can2","",900,800 ) ;

      h_rel_err_table_simple_ignore -> SetContour(nb) ;
      h_rel_err_table_simple_ignore -> SetMinimum( -0.7 ) ;
      h_rel_err_table_simple_ignore -> SetMaximum(  0.7 ) ;

      h_rel_err_table_simple_ignore -> Draw("colz") ;
      h_rel_err_table_simple_ignore -> Draw("text same") ;

   } // gen_combine_table_from_cov_mat





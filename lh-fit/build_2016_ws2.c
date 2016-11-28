
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooUniform.h"
#include "RooStats/ModelConfig.h"

#include "RooProdPdfLogSum.h"
#include "RooPoissonLogEval.h"

#include "../binning.h"

#include <fstream>

   using namespace RooFit ;
   using namespace RooStats ;

   RooArgSet* globalObservables ;
   RooArgSet* allNuisances ;
   RooArgSet* allNuisancePdfs ;
   RooArgSet* pdf_sbIndexList ;
   RooArgSet* allBGmuPars ;

   float get_par_max( float mu ) ;

   RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err, int sbi = -1 ) ;
   RooAbsReal* makeCorrelatedLognormalConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false ) ;

   RooAbsReal* read_bg_line( ifstream& ifs, const char* bin_name ) ;

   RooWorkspace* ws_pointer ;

  //=================================================================================

   void build_2016_ws2(
                            const char* outfile = "outputfiles/ws-lhfit-test2.root",
                            const char* fname_ldp_lostlep = "outputfiles/nbsum-ldp-input-lostlep.txt",
                            const char* fname_ldp_hadtau  = "outputfiles/nbsum-ldp-input-hadtau.txt",
                            const char* fname_ldp_znunu   = "outputfiles/nbsum-ldp-input-znunu.txt",
                            const char* fname_hdp_lostlep = "outputfiles/nbsum-hdp-input-lostlep.txt",
                            const char* fname_hdp_hadtau  = "outputfiles/nbsum-hdp-input-hadtau.txt",
                            const char* fname_hdp_znunu   = "outputfiles/nbsum-hdp-input-znunu.txt",
                            //////const char* fname_data    = "../outputfiles/nbsum-input-fakedata2.txt",
                            const char* fname_data    = "../outputfiles/nbsum-input-data.txt",
                            const char* fname_sigmc   = "../outputfiles/nbsum-input-T1bbbbH.txt"
                          ) {


      setup_bins() ;
      bool no_rounding(true) ;

      char pname[100] ;
      char pname2[100] ;

      bool skip_testfit(false) ;
      bool skip_modelconfig(true) ;


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      printf("\n\n Creating workspace.\n\n") ;

      RooWorkspace workspace("ws") ;
      ws_pointer = &workspace ;

      workspace.autoImportClassCode(true) ;

      globalObservables      = new RooArgSet("globalObservables");
      allNuisances           = new RooArgSet("allNuisances");
      allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      pdf_sbIndexList        = new RooArgSet("pdf_sbIndexList") ;
      allBGmuPars            = new RooArgSet("allBGmuPars") ;
      RooArgSet* observedParametersList = new RooArgSet("observables") ;

      RooArgSet pdflist ;

      sprintf( pname, "sig_strength" ) ;
      RooRealVar* rv_sig_strength = new RooRealVar( pname, pname, 0.0, 0., 10. ) ; // nominal value is zero
      rv_sig_strength -> setConstant(kTRUE) ; // fixed
      rv_sig_strength -> Print() ;
      printf("  %s\n\n", pname ) ;


      RooArgSet sbIndexList("sbIndexList") ;






      printf("\n  QCD Model parameters:\n") ;

      RooAbsReal* rv_qcd_kht[10] ;
      RooAbsReal* rv_qcd_snjet[10] ;

      int n_qcd_kht_pars(3) ;

      rv_qcd_kht[1] = new RooRealVar( "Kqcd_ht1", "Kqcd_ht1", 0.71, 0., 5. ) ;
      rv_qcd_kht[2] = new RooRealVar( "Kqcd_ht2", "Kqcd_ht2", 0.09, 0., 5. ) ;
      rv_qcd_kht[3] = new RooRealVar( "Kqcd_ht3", "Kqcd_ht3", 0.02, 0., 5. ) ;

      printf("\n Njet parameters:\n") ;
      for ( int bi=1; bi<=nb_nj; bi++ ) {
         char pname[100] ;
         sprintf( pname, "Sqcd_njet%d", bi ) ;
         float initial_val(1.0) ;
         if ( bi==1 ) initial_val = 0.67 ;
         if ( bi==2 ) initial_val = 1.00 ;
         if ( bi==3 ) initial_val = 1.45 ;
         if ( bi==4 ) initial_val = 2.85 ;
         if ( bi==5 ) initial_val = 7.96 ;
         /////rv_qcd_snjet[bi] = new RooRealVar( pname, pname, initial_val, 0., 40. ) ;
         rv_qcd_snjet[bi] = new RooRealVar( pname, pname, initial_val, 0., 200. ) ;
         printf(" %s ", pname ) ;
         if ( bi == (njet_bin_to_be_fixed_in_qcd_model_fit+1) ) {
            ((RooRealVar*) rv_qcd_snjet[bi]) -> setConstant( kTRUE ) ;
            printf(" fixed to 1\n") ;
         } else {
            printf("\n") ;
         }
      }



      ifstream ifs_ldp_lostlep ;
      ifstream ifs_ldp_hadtau ;
      ifstream ifs_ldp_znunu ;

      ifstream ifs_hdp_lostlep ;
      ifstream ifs_hdp_hadtau ;
      ifstream ifs_hdp_znunu ;

      ifstream ifs_data ;
      ifstream ifs_sigmc ;


      ifs_ldp_lostlep.open( fname_ldp_lostlep ) ;
      if ( !ifs_ldp_lostlep.good() ) { printf("\n\n *** Bad input lostlep LDP file: %s\n\n", fname_ldp_lostlep ) ; return ; }

      ifs_ldp_hadtau.open( fname_ldp_hadtau ) ;
      if ( !ifs_ldp_hadtau.good() ) { printf("\n\n *** Bad input hadtau LDP file: %s\n\n", fname_ldp_hadtau ) ; return ; }

      ifs_ldp_znunu.open( fname_ldp_znunu ) ;
      if ( !ifs_ldp_znunu.good() ) { printf("\n\n *** Bad input znunu LDP file: %s\n\n", fname_ldp_znunu ) ; return ; }


      ifs_hdp_lostlep.open( fname_hdp_lostlep ) ;
      if ( !ifs_hdp_lostlep.good() ) { printf("\n\n *** Bad input lostlep HDP file: %s\n\n", fname_hdp_lostlep ) ; return ; }

      ifs_hdp_hadtau.open( fname_hdp_hadtau ) ;
      if ( !ifs_hdp_hadtau.good() ) { printf("\n\n *** Bad input hadtau HDP file: %s\n\n", fname_hdp_hadtau ) ; return ; }

      ifs_hdp_znunu.open( fname_hdp_znunu ) ;
      if ( !ifs_hdp_znunu.good() ) { printf("\n\n *** Bad input znunu HDP file: %s\n\n", fname_hdp_znunu ) ; return ; }


      ifs_data.open( fname_data ) ;
      if ( !ifs_data.good() ) { printf("\n\n *** Bad input data file: %s\n\n", fname_data ) ; return ; }

      ifs_sigmc.open( fname_sigmc ) ;
      if ( !ifs_sigmc.good() ) { printf("\n\n *** Bad input sigmc file: %s\n\n", fname_sigmc ) ; return ; }






      printf("\n\n Reading input event counts.\n\n" ) ;

      while ( ifs_data.good() ) {

         printf("\n") ;

         TString ts ;

       //-- data
         ts.ReadLine( ifs_data ) ;
         if ( !ifs_data.good() ) { printf("\n\n  Reached end of data input file.\n\n" ) ; break ; }
         char bin_name[100] ;
         int  data_ldp_val(0), data_hdp_val(0) ;
         sscanf( ts.Data(), "%s %d %d", bin_name, &data_ldp_val, &data_hdp_val ) ;

         int fb_nji(-1), fb_hbi(-1) ;
         sscanf( bin_name, "Nj%d-HT%d", &fb_nji, &fb_hbi ) ;

        //-- skip the garbage bins in the likelihood.
         if ( fb_hbi==1 && fb_nji>(nb_nj-2) ) {
            ts.ReadLine( ifs_sigmc ) ;
            continue ;
         }

    //  //-- fuck the HT1 bins
    //   if ( fb_hbi==1 ) {
    //      ts.ReadLine( ifs_sigmc ) ;
    //      continue ;
    //   }


         printf( " %30s,  Njet index %d,  HT index %d\n", bin_name, fb_nji, fb_hbi ) ;
         printf( "      data    : LDP %5d            ,   HDP %4d\n", data_ldp_val, data_hdp_val ) ;


         fflush(stdout) ;


         sprintf( pname, "R_qcd_ldp_%s", bin_name ) ;
         RooFormulaVar* rv_R_qcd_ldp = new RooFormulaVar( pname, "@0 * @1",
            RooArgSet(
               *(rv_qcd_kht[ fb_hbi ]),
               *(rv_qcd_snjet[ fb_nji ])
                     ) ) ;
         rv_R_qcd_ldp -> Print() ;




         sprintf( pname, "Nldp_%s", bin_name ) ;
         RooRealVar* rv_Nldp = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rv_Nldp -> setVal( data_ldp_val ) ;
         rv_Nldp -> setConstant( kTRUE ) ;
         rv_Nldp -> Print() ;
         observedParametersList -> add( *rv_Nldp ) ;

         sprintf( pname, "Nhdp_%s", bin_name ) ;
         RooRealVar* rv_Nhdp = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rv_Nhdp -> setVal( data_hdp_val ) ;
         rv_Nhdp -> setConstant( kTRUE ) ;
         rv_Nhdp -> Print() ;
         observedParametersList -> add( *rv_Nhdp ) ;




       //-- lostlep
     /// ts.ReadLine( ifs_lostlep ) ;
     /// char lostlep_bin_name[100] ;
     /// float lostlep_ldp_val, lostlep_ldp_err_stat, lostlep_ldp_err_syst,   lostlep_hdp_val, lostlep_hdp_err_stat, lostlep_hdp_err_syst ;
     /// sscanf( ts.Data(), "%s %f +/- %f +/- %f %f +/- %f +/- %f", lostlep_bin_name, &lostlep_ldp_val, &lostlep_ldp_err_stat, &lostlep_ldp_err_syst, &lostlep_hdp_val, &lostlep_hdp_err_stat, &lostlep_hdp_err_syst ) ;
     /// if ( strcmp( bin_name, lostlep_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, lostlep_bin_name ) ; return ; }
     /// printf("      lostlep : LDP %7.1f +/- %5.1f +/- %5.1f,   LDP %6.1f +/- %5.1f +/- %5.1f\n", lostlep_ldp_val, lostlep_ldp_err_stat, lostlep_ldp_err_syst, lostlep_hdp_val, lostlep_hdp_err_stat, lostlep_hdp_err_syst ) ;

     /// RooAbsReal* rv_mu_lostlep_ldp(0x0) ;
     /// sprintf( pname, "mu_lostlep_ldp_statonly_%s", bin_name ) ;
     /// RooAbsReal* rv_mu_lostlep_ldp_statonly = makeLognormalConstraint( pname, lostlep_ldp_val, lostlep_ldp_err_stat ) ;
     /// double lostlep_ldp_syserr_frac(0.001) ;
     /// if ( lostlep_ldp_err_syst > 0. && lostlep_ldp_val > 0. ) { lostlep_ldp_syserr_frac = lostlep_ldp_err_syst / lostlep_ldp_val ; }
     /// sprintf( pname, "mu_lostlep_ldp_systerr_%s", bin_name ) ;
     /// RooAbsReal* rv_mu_lostlep_ldp_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, lostlep_ldp_syserr_frac, "syst_lostlep_ldp" ) ;
     /// sprintf( pname, "mu_lostlep_ldp_%s", bin_name ) ;
     /// rv_mu_lostlep_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_mu_lostlep_ldp_statonly, *rv_mu_lostlep_ldp_systerr ) ) ;


     /// RooAbsReal* rv_mu_lostlep_hdp(0x0) ;
     /// sprintf( pname, "mu_lostlep_hdp_statonly_%s", bin_name ) ;
     /// RooAbsReal* rv_mu_lostlep_hdp_statonly = makeLognormalConstraint( pname, lostlep_hdp_val, lostlep_hdp_err_stat ) ;
     /// double lostlep_hdp_syserr_frac(0.001) ;
     /// if ( lostlep_hdp_err_syst > 0. && lostlep_hdp_val > 0. ) { lostlep_hdp_syserr_frac = lostlep_hdp_err_syst / lostlep_hdp_val ; }
     /// sprintf( pname, "mu_lostlep_hdp_systerr_%s", bin_name ) ;
     /// RooAbsReal* rv_mu_lostlep_hdp_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, lostlep_hdp_syserr_frac, "syst_lostlep_hdp" ) ;
     /// sprintf( pname, "mu_lostlep_hdp_%s", bin_name ) ;
     /// rv_mu_lostlep_hdp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_mu_lostlep_hdp_statonly, *rv_mu_lostlep_hdp_systerr ) ) ;

         RooAbsReal* rv_mu_lostlep_ldp(0x0) ;
         rv_mu_lostlep_ldp = read_bg_line( ifs_ldp_lostlep, bin_name ) ;

         RooAbsReal* rv_mu_lostlep_hdp(0x0) ;
         rv_mu_lostlep_hdp = read_bg_line( ifs_hdp_lostlep, bin_name ) ;



       //-- hadtau
   ////  ts.ReadLine( ifs_hadtau ) ;
   ////  char hadtau_bin_name[100] ;
   ////  float hadtau_ldp_val, hadtau_ldp_err_stat, hadtau_ldp_err_syst,   hadtau_hdp_val, hadtau_hdp_err_stat, hadtau_hdp_err_syst ;
   ////  sscanf( ts.Data(), "%s %f +/- %f +/- %f %f +/- %f +/- %f", hadtau_bin_name, &hadtau_ldp_val, &hadtau_ldp_err_stat, &hadtau_ldp_err_syst, &hadtau_hdp_val, &hadtau_hdp_err_stat, &hadtau_hdp_err_syst ) ;
   ////  if ( strcmp( bin_name, hadtau_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, hadtau_bin_name ) ; return ; }
   ////  printf("      hadtau : LDP %7.1f +/- %5.1f +/- %5.1f,   LDP %6.1f +/- %5.1f +/- %5.1f\n", hadtau_ldp_val, hadtau_ldp_err_stat, hadtau_ldp_err_syst, hadtau_hdp_val, hadtau_hdp_err_stat, hadtau_hdp_err_syst ) ;

   ////  RooAbsReal* rv_mu_hadtau_ldp(0x0) ;
   ////  sprintf( pname, "mu_hadtau_ldp_statonly_%s", bin_name ) ;
   ////  RooAbsReal* rv_mu_hadtau_ldp_statonly = makeLognormalConstraint( pname, hadtau_ldp_val, hadtau_ldp_err_stat ) ;
   ////  double hadtau_ldp_syserr_frac(0.001) ;
   ////  if ( hadtau_ldp_err_syst > 0. && hadtau_ldp_val > 0. ) { hadtau_ldp_syserr_frac = hadtau_ldp_err_syst / hadtau_ldp_val ; }
   ////  sprintf( pname, "mu_hadtau_ldp_systerr_%s", bin_name ) ;
   ////  RooAbsReal* rv_mu_hadtau_ldp_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, hadtau_ldp_syserr_frac, "syst_hadtau_ldp" ) ;
   ////  sprintf( pname, "mu_hadtau_ldp_%s", bin_name ) ;
   ////  rv_mu_hadtau_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_mu_hadtau_ldp_statonly, *rv_mu_hadtau_ldp_systerr ) ) ;

   ////  RooAbsReal* rv_mu_hadtau_hdp(0x0) ;
   ////  sprintf( pname, "mu_hadtau_hdp_statonly_%s", bin_name ) ;
   ////  RooAbsReal* rv_mu_hadtau_hdp_statonly = makeLognormalConstraint( pname, hadtau_hdp_val, hadtau_hdp_err_stat ) ;
   ////  double hadtau_hdp_syserr_frac(0.001) ;
   ////  if ( hadtau_hdp_err_syst > 0. && hadtau_hdp_val > 0. ) { hadtau_hdp_syserr_frac = hadtau_hdp_err_syst / hadtau_hdp_val ; }
   ////  sprintf( pname, "mu_hadtau_hdp_systerr_%s", bin_name ) ;
   ////  RooAbsReal* rv_mu_hadtau_hdp_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, hadtau_hdp_syserr_frac, "syst_hadtau_hdp" ) ;
   ////  sprintf( pname, "mu_hadtau_hdp_%s", bin_name ) ;
   ////  rv_mu_hadtau_hdp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_mu_hadtau_hdp_statonly, *rv_mu_hadtau_hdp_systerr ) ) ;

         RooAbsReal* rv_mu_hadtau_ldp(0x0) ;
         rv_mu_hadtau_ldp = read_bg_line( ifs_ldp_hadtau, bin_name ) ;

         RooAbsReal* rv_mu_hadtau_hdp(0x0) ;
         rv_mu_hadtau_hdp = read_bg_line( ifs_hdp_hadtau, bin_name ) ;




       //-- znunu
   ///   ts.ReadLine( ifs_znunu ) ;
   ///   char znunu_bin_name[100] ;
   ///   float znunu_ldp_val, znunu_ldp_err_stat, znunu_ldp_err_syst,   znunu_hdp_val, znunu_hdp_err_stat, znunu_hdp_err_syst ;
   ///   sscanf( ts.Data(), "%s %f +/- %f +/- %f %f +/- %f +/- %f", znunu_bin_name, &znunu_ldp_val, &znunu_ldp_err_stat, &znunu_ldp_err_syst, &znunu_hdp_val, &znunu_hdp_err_stat, &znunu_hdp_err_syst ) ;
   ///   if ( strcmp( bin_name, znunu_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, znunu_bin_name ) ; return ; }
   ///   printf("      znunu : LDP %7.1f +/- %5.1f +/- %5.1f,   LDP %6.1f +/- %5.1f +/- %5.1f\n", znunu_ldp_val, znunu_ldp_err_stat, znunu_ldp_err_syst, znunu_hdp_val, znunu_hdp_err_stat, znunu_hdp_err_syst ) ;

   ///   RooAbsReal* rv_mu_znunu_ldp(0x0) ;
   ///   sprintf( pname, "mu_znunu_ldp_statonly_%s", bin_name ) ;
   ///   RooAbsReal* rv_mu_znunu_ldp_statonly = makeLognormalConstraint( pname, znunu_ldp_val, znunu_ldp_err_stat ) ;
   ///   double znunu_ldp_syserr_frac(0.001) ;
   ///   if ( znunu_ldp_err_syst > 0. && znunu_ldp_val > 0. ) { znunu_ldp_syserr_frac = znunu_ldp_err_syst / znunu_ldp_val ; }
   ///   sprintf( pname, "mu_znunu_ldp_systerr_%s", bin_name ) ;
   ///   RooAbsReal* rv_mu_znunu_ldp_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, znunu_ldp_syserr_frac, "syst_znunu_ldp" ) ;
   ///   sprintf( pname, "mu_znunu_ldp_%s", bin_name ) ;
   ///   rv_mu_znunu_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_mu_znunu_ldp_statonly, *rv_mu_znunu_ldp_systerr ) ) ;

   ///   RooAbsReal* rv_mu_znunu_hdp(0x0) ;
   ///   sprintf( pname, "mu_znunu_hdp_statonly_%s", bin_name ) ;
   ///   RooAbsReal* rv_mu_znunu_hdp_statonly = makeLognormalConstraint( pname, znunu_hdp_val, znunu_hdp_err_stat ) ;
   ///   double znunu_hdp_syserr_frac(0.001) ;
   ///   if ( znunu_hdp_err_syst > 0. && znunu_hdp_val > 0. ) { znunu_hdp_syserr_frac = znunu_hdp_err_syst / znunu_hdp_val ; }
   ///   sprintf( pname, "mu_znunu_hdp_systerr_%s", bin_name ) ;
   ///   RooAbsReal* rv_mu_znunu_hdp_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, znunu_hdp_syserr_frac, "syst_znunu_hdp" ) ;
   ///   sprintf( pname, "mu_znunu_hdp_%s", bin_name ) ;
   ///   rv_mu_znunu_hdp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_mu_znunu_hdp_statonly, *rv_mu_znunu_hdp_systerr ) ) ;


         RooAbsReal* rv_mu_znunu_ldp(0x0) ;
         rv_mu_znunu_ldp = read_bg_line( ifs_ldp_znunu, bin_name ) ;

         RooAbsReal* rv_mu_znunu_hdp(0x0) ;
         rv_mu_znunu_hdp = read_bg_line( ifs_hdp_znunu, bin_name ) ;







       //-- sigmc
         ts.ReadLine( ifs_sigmc ) ;
         char sigmc_bin_name[100] ;
         float sigmc_ldp_val, sigmc_hdp_val ;
         sscanf( ts.Data(), "%s %f  %f", sigmc_bin_name, &sigmc_ldp_val, &sigmc_hdp_val ) ;
         if ( strcmp( bin_name, sigmc_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, sigmc_bin_name ) ; return ; }
         printf("      sigmc   : LDP %7.1f ,   LDP %6.1f\n", sigmc_ldp_val, sigmc_hdp_val ) ;

         sprintf( pname, "mu_sig0_ldp_%s", bin_name ) ;
         RooRealVar* rv_mu_sig0_ldp = new RooRealVar( pname, pname, sigmc_ldp_val, 0., 1.e6 ) ;
         rv_mu_sig0_ldp -> setConstant( kTRUE ) ;
         rv_mu_sig0_ldp -> Print() ;
         sprintf( pname, "mu_sig_ldp_%s", bin_name ) ;
         RooFormulaVar* rv_mu_sig_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_ldp ) ) ;
         rv_mu_sig_ldp -> Print() ;

         sprintf( pname, "mu_sig0_hdp_%s", bin_name ) ;
         RooRealVar* rv_mu_sig0_hdp = new RooRealVar( pname, pname, sigmc_hdp_val, 0., 1.e6 ) ;
         rv_mu_sig0_hdp -> setConstant( kTRUE ) ;
         rv_mu_sig0_hdp -> Print() ;
         sprintf( pname, "mu_sig_hdp_%s", bin_name ) ;
         RooFormulaVar* rv_mu_sig_hdp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_hdp ) ) ;
         rv_mu_sig_hdp -> Print() ;










         double lostlep_ldp_val = rv_mu_lostlep_ldp -> getVal() ;
         double hadtau_ldp_val  = rv_mu_hadtau_ldp  -> getVal() ;
         double znunu_ldp_val   = rv_mu_znunu_ldp   -> getVal() ;

         double initial_qcd_ldp_val =  data_ldp_val - ( lostlep_ldp_val + hadtau_ldp_val + znunu_ldp_val ) ;
         if ( initial_qcd_ldp_val < 0. ) initial_qcd_ldp_val = 0. ;

         sprintf( pname, "mu_qcd_ldp_%s", bin_name ) ;
         RooRealVar* rv_mu_qcd_ldp = new RooRealVar( pname, pname, initial_qcd_ldp_val, 0., get_par_max( data_ldp_val ) ) ;
         rv_mu_qcd_ldp -> setConstant( kFALSE ) ;
         rv_mu_qcd_ldp -> Print() ;
         allBGmuPars -> add( *rv_mu_qcd_ldp ) ;


         sprintf( pname, "mu_qcd_hdp_%s", bin_name ) ;
         RooFormulaVar* rv_mu_qcd_hdp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_R_qcd_ldp, *rv_mu_qcd_ldp ) ) ;
         rv_mu_qcd_hdp -> Print() ;





         sprintf( pname, "n_ldp_%s", bin_name ) ;
         RooFormulaVar* rv_n_ldp = new RooFormulaVar( pname, "@0 + @1 + @2 + @3 + @4", RooArgSet(
             *rv_mu_qcd_ldp, *rv_mu_lostlep_ldp, *rv_mu_hadtau_ldp, *rv_mu_znunu_ldp, *rv_mu_sig_ldp ) ) ;
         rv_n_ldp -> Print() ;

         sprintf( pname, "n_hdp_%s", bin_name ) ;
         RooFormulaVar* rv_n_hdp = new RooFormulaVar( pname, "@0 + @1 + @2 + @3 + @4", RooArgSet(
             *rv_mu_qcd_hdp, *rv_mu_lostlep_hdp, *rv_mu_hadtau_hdp, *rv_mu_znunu_hdp, *rv_mu_sig_hdp ) ) ;
         rv_n_hdp -> Print() ;





         sprintf( pname, "pdf_ldp_%s", bin_name ) ;
         RooPoissonLogEval* rv_pdf_ldp = new RooPoissonLogEval( pname, pname, *rv_Nldp, *rv_n_ldp, no_rounding ) ;
         rv_pdf_ldp -> Print() ;

         pdflist.add( *rv_pdf_ldp ) ;


         sprintf( pname, "pdf_hdp_%s", bin_name ) ;
         RooPoissonLogEval* rv_pdf_hdp = new RooPoissonLogEval( pname, pname, *rv_Nhdp, *rv_n_hdp, no_rounding ) ;
         rv_pdf_hdp -> Print() ;

         pdflist.add( *rv_pdf_hdp ) ;



      } // reading bins from data file ( ifs_data )



      printf("\n\n Creating and importing dataset into workspace.\n\n") ;

      RooDataSet* dsObserved = new RooDataSet("observed_rds", "observed data values", *observedParametersList ) ;
      dsObserved -> add( *observedParametersList ) ;
      workspace.import( *dsObserved ) ;



      pdflist.add( *allNuisancePdfs ) ;

      printf("\n List of all PDFs\n") ;
      pdflist.Print() ;
      printf("\n") ;

      RooProdPdfLogSum* likelihood = new RooProdPdfLogSum( "likelihood", "likelihood", pdflist ) ;
      likelihood->Print() ;







      if ( !skip_testfit ) {

         printf("\n\n Running a test fit.\n\n") ;

         printf("\n\n =============================================\n\n") ;
         ////likelihood -> fitTo( *dsObserved, PrintLevel(3), Hesse(0), Minos(0) ) ;
         ///likelihood -> fitTo( *dsObserved, Optimize(0), PrintLevel(3), Hesse(0), Minos(0) ) ;
         likelihood -> fitTo( *dsObserved, Optimize(0), PrintLevel(3), Hesse(0), Minos(0) ) ;
         printf("\n\n =============================================\n\n") ;

      }


      if ( !skip_modelconfig ) {

        //-- Set up RooStats models.

         printf("\n\n Setting up S+B model.\n\n") ;

         RooArgSet poi( *rv_sig_strength, "poi" ) ;
         RooUniform signal_prior( "signal_prior", "signal_prior", *rv_sig_strength ) ;

         ModelConfig sbModel ("SbModel");
         sbModel.SetWorkspace( workspace ) ;
         sbModel.SetPdf( *likelihood ) ;
         sbModel.SetParametersOfInterest( poi );
         sbModel.SetPriorPdf(signal_prior);
         sbModel.SetObservables( *observedParametersList );
         sbModel.SetNuisanceParameters( *allNuisances );
         sbModel.SetGlobalObservables( *globalObservables );

         workspace.Print() ;

         printf("\n\n Doing fit for S+B model.\n" ) ; fflush(stdout) ;

         RooAbsReal* pNll = sbModel.GetPdf()->createNLL(*dsObserved);
         RooAbsReal* pProfile = pNll->createProfile(RooArgSet());
         pProfile->getVal();
         RooArgSet* pPoiAndNuisance = new RooArgSet();
         pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
         if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
         printf("\n\n Will save these parameter points that correspond to the fit to data.\n\n") ; fflush(stdout) ;
         pPoiAndNuisance->Print("v");
         sbModel.SetSnapshot(*pPoiAndNuisance);
         workspace.import (sbModel);

         delete pProfile ;
         delete pNll ;
         delete pPoiAndNuisance ;

         printf("\n\n Setting up BG-only model.\n\n") ;

         ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
         bModel.SetName("BModel");
         bModel.SetWorkspace(workspace);

         printf("\n\n Doing fit for BG-only model.\n" ) ; fflush(stdout) ;
         pNll = bModel.GetPdf()->createNLL(*dsObserved);
         pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
         ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.);
         pProfile->getVal();
         pPoiAndNuisance = new RooArgSet();
         pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
         if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
         printf("\n\n Should use these parameter points to generate pseudo data for bkg only.\n\n") ; fflush(stdout) ;
         pPoiAndNuisance->Print("v");
         bModel.SetSnapshot(*pPoiAndNuisance);
         workspace.import (bModel);

         delete pProfile ;
         delete pNll ;
         delete pPoiAndNuisance ;

      } else {

         workspace.import( *likelihood ) ;

      }

      workspace.defineSet( "all_nuisance_pars", *allNuisances, kFALSE ) ;  // for convenience. do not import if missing.  should already be in there.
      workspace.defineSet( "all_nuisance_pdfs", *allNuisancePdfs, kFALSE ) ;  // for convenience. do not import if missing.  should already be in there.
      workspace.defineSet( "all_bg_mu_pars", *allBGmuPars, kFALSE ) ;  // for convenience. do not import if missing.  should already be in there.

      workspace.Print() ;

      printf("\n\n Saving workspace in : %s\n\n", outfile ) ;

      gSystem->Exec(" mkdir -p outputfiles " ) ;

      workspace.writeToFile( outfile ) ;














   } // build_2016_ws2

  //=================================================================================

    float get_par_max( float mu ) {

   /// if ( mu <= 1. ) return 10. ;
   /// if ( mu <= 2. ) return 15. ;
   /// if ( mu <= 50. ) return mu + 4*sqrt(mu) ;
   /// if ( mu <= 100. ) return mu + 3*sqrt(mu) ;
   /// return mu + 2*sqrt(mu) ;

       if ( mu <= 1. ) return 10. ;
       if ( mu <= 2. ) return 15. ;
       if ( mu <= 50. ) return mu + 6*sqrt(mu) ;
       if ( mu <= 100. ) return mu + 5*sqrt(mu) ;
       return mu + 4*sqrt(mu) ;


    } // get_par_max

  //=================================================================================


    RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err, int sbi ) {

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       if ( NP_err <= 0. ) {
          printf(" makeLognormalConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          ws_pointer -> import( *g_mean ) ;
          ws_pointer -> import( *g_sigma ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       char pname[1000] ;
       sprintf( pname, "prim_%s", NP_name ) ;

       printf(" makeLognormalConstraint : creating primary log-normal variable %s\n", pname ) ;
       RooRealVar* np_prim_rrv = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_rrv -> setVal( 0. ) ;
       np_prim_rrv -> setConstant( kFALSE ) ;

       sprintf( pname, "prim_mean_%s", NP_name ) ;
       RooRealVar* np_prim_mean = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_mean->setConstant(kTRUE) ;

       sprintf( pname, "prim_sigma_%s", NP_name ) ;
       RooConstVar* np_prim_sigma = new RooConstVar( pname, pname, 1. ) ;


       char pdfname[1000] ;
       sprintf( pdfname, "pdf_prim_%s", NP_name ) ;
       RooGaussian* np_prim_pdf = new RooGaussian( pdfname, pdfname, *np_prim_rrv, *np_prim_mean, *np_prim_sigma ) ;

       allNuisances -> add( *np_prim_rrv ) ;
       allNuisancePdfs -> add( *np_prim_pdf ) ;
       globalObservables -> add( *np_prim_mean ) ;

       sprintf( pname, "%s_sb_index", np_prim_pdf->GetName() ) ;
       RooConstVar* rv_pdf_sb_index = new RooConstVar( pname, pname, sbi ) ;
       pdf_sbIndexList -> add( *rv_pdf_sb_index ) ;



       //-- compute the log-normal-distributed parameter from the primary parameter.

       //--- This is the new way.  RMS of lognormal is much closer to sigma when sigma is
       //    large, doing it this way.  When sigma/mean is small, they are about the same.
       //    That is, exp(sigma/mean) is close to (sigma/mean + 1).  This one is better when
       //    sigma/mean is not small.  The high-side tail is not as strong.
       //
        RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( ( @1/@0 + 1. ), @2)",
                  RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;
       //------------------------------------------------------------------------------------------



       printf("  makeLognormalConstraint : created log-normal nuisance parameter %s : val = %g\n", NP_name, np_rfv -> getVal() ) ;

       return np_rfv ;


    } // makeLognormalConstraint.


   //==============================================================================================================

    RooAbsReal* makeCorrelatedLognormalConstraint(
            const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign ) {

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* ln_mean  = new RooConstVar( vname, vname, NP_val ) ;

       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* ln_sigma = new RooConstVar( vname, vname, NP_err ) ;


       if ( NP_err <= 0. ) {
          printf("  makeCorrelatedLognormalConstraint: Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          ws_pointer -> import( *ln_mean ) ;
          ws_pointer -> import( *ln_sigma ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       char prim_name[1000] ;
       sprintf( prim_name, "prim_%s", NP_base_name ) ;
       RooRealVar* rrv_np_base_par = (RooRealVar*) allNuisances -> find( prim_name ) ;

       if ( rrv_np_base_par == 0x0 ) {

          printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter - %s\n\n", prim_name ) ;
          rrv_np_base_par = new RooRealVar( prim_name, prim_name, -6.0, 6.0 ) ;
          rrv_np_base_par -> setVal( 0. ) ;
          rrv_np_base_par -> setConstant( kFALSE ) ;
          allNuisances -> add( *rrv_np_base_par ) ;

          char vname[1000] ;
          sprintf( vname, "prim_mean_%s", NP_base_name ) ;
          RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-10.,10. ) ;
          g_mean->setConstant(kTRUE);
          sprintf( vname, "prim_sigma_%s", NP_base_name ) ;
          RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

          char pdfname[100] ;
          sprintf( pdfname, "pdf_prim_%s", NP_base_name ) ;
          printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
          RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;

          allNuisancePdfs -> add( *base_np_pdf ) ;
          globalObservables -> add( *g_mean ) ;

       }


       RooAbsReal* rar(0x0) ;


       char formula[1000] ;

       if ( !changeSign ) {
          sprintf( formula, "@0 * pow( ( @1/@0 + 1.), @2 )" ) ;
       } else {
          sprintf( formula, "@0 * pow( ( @1/@0 + 1.), -1.0 * @2 )" ) ;
       }

       rar = new RooFormulaVar( NP_name, formula, RooArgSet( *ln_mean, *ln_sigma, *rrv_np_base_par ) ) ;

       printf(" makeCorrelatedLognormalConstraint : creating correlated log-normal NP with formula : %s,  %s, val = %g, mean=%g, sigma=%g\n", formula, NP_name, rar->getVal(), NP_val, NP_err ) ;


       return rar ;

    } // makeCorrelatedLognormalConstraint.

   //==============================================================================================================

   RooAbsReal* read_bg_line( ifstream& ifs, const char* bin_name ) {

      bool verb(true) ;

      ifs.seekg(0) ;

      TString ts ;

      ts.ReadLine( ifs ) ;
      char bg_name[100] ;
      char region_name[100] ;
      int nsyst(0) ;
      sscanf( ts.Data(), "%s %s nsyst=%d", bg_name, region_name, &nsyst ) ;

      if (verb) { printf(" read_bg_line : %s %s has %d systs.\n", bg_name, region_name, nsyst ) ; }

      //-- read first 5 words, expecting to end in total-syst
      for ( int i=0; i<5; i++ ) {
         ts.ReadToken( ifs ) ;
         if ( i==4 && strcmp( "total-syst", ts.Data() ) != 0 ) {
            printf(" expected total-syst, got %s\n", ts.Data() ) ;
            gSystem->Exit(-1) ;
         }
      }

      char syst_name[50][100] ;
      for ( int si=0; si<nsyst; si++ ) {
         ts.ReadToken( ifs ) ;
         sprintf( syst_name[si], "%s", ts.Data() ) ;
         ////////printf("   syst %2d : %s\n", si, syst_name[si] ) ;
      }

      bool found_the_bin(false) ;
      while ( ifs.good() ) {
         ts.ReadToken( ifs ) ;
         if ( !ifs.good() ) break ;
         if ( ts.Contains( bin_name ) ) {
            found_the_bin = true ;
            break ;
         }
      }

      if ( !found_the_bin ) {
         printf("\n\n *** read_bg_line: failed to find bin %s\n\n\n", bin_name ) ;
         gSystem -> Exit(-1) ;
      }

      double val(0.), stat_err(0.) ;
      ts.ReadToken( ifs ) ; sscanf( ts.Data(), "%lf", &val ) ;
      ts.ReadToken( ifs ) ; sscanf( ts.Data(), "%lf", &stat_err ) ;

      char formula_string[1000] ;
      RooArgSet* ras_product_terms = new RooArgSet() ;
      sprintf( formula_string, "@0" ) ;

      char pname[1000] ;
      sprintf( pname, "mu_%s_%s_statonly_%s", bg_name, region_name, bin_name ) ;
      RooAbsReal* rv_mu_statonly = makeLognormalConstraint( pname, val, stat_err ) ;
      ras_product_terms -> add( *rv_mu_statonly ) ;

      double total_syst_err_events(0.), total_syst_err_rel(0.) ;
      ts.ReadToDelim( ifs, ')') ;
      sscanf( ts.Data(), "%lf (%lf)", &total_syst_err_events, &total_syst_err_rel ) ;
      if (verb) {
         printf(" read_bg_line : %s %s %s val=%8.1f +/- %5.1f +/- %5.1f (%7.4f)\n",
                bg_name, region_name, bin_name, val, stat_err, total_syst_err_events, total_syst_err_rel ) ;
      }


      for ( int si=0; si<nsyst; si++ ) {
         ts.ReadToDelim( ifs, ')') ;
         double syst_err_evts, syst_err_rel ;
         sscanf( ts.Data(), "%lf (%lf)", &syst_err_evts, &syst_err_rel ) ;
         sprintf( pname, "mu_%s_%s_systerr_%s_%s", bg_name, region_name, syst_name[si], bin_name ) ;
         char systparname[1000] ;
         sprintf( systparname, "syst_%s_%s", bg_name, syst_name[si] ) ;
         RooAbsReal* rv_mu_systerr = makeCorrelatedLognormalConstraint( pname, 1.0, syst_err_rel, systparname ) ;
         ras_product_terms -> add( *rv_mu_systerr ) ;
         char tmpstr[1000] ;
         sprintf( tmpstr, "%s", formula_string ) ;
         sprintf( formula_string, "%s * @%d", tmpstr, si+1 ) ;
      } // si
      sprintf( pname, "mu_%s_%s_%s", bg_name, region_name, bin_name ) ;
      RooAbsReal* rv_mu = new RooFormulaVar( pname, formula_string, *ras_product_terms ) ;
      if ( verb ) { rv_mu -> Print() ; }

      return rv_mu ;


   } // read_bg_line

   //==============================================================================================================



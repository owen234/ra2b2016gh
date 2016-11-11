#include "TDirectory.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TStyle.h"

#include "histio.c"
#include "get_hist.h"

  //--------

   void calc_mc_hlratios1(
         const char* iodir          = "outputfiles/",
         const char* infile_qcd     = "hists-v2d-qcd.root",
         const char* infile_lostlep = "hists-v2c-lostlep.root",
         const char* infile_hadtau  = "hists-v2c-hadtau.root",
         const char* infile_znunu   = "hists-v2c-znunu.root"
      ) {

     char fname[10000] ;

     gStyle -> SetPadBottomMargin( 0.35 ) ;
     gStyle -> SetOptStat(0) ;

     gDirectory -> Delete( "h*" ) ;

     sprintf( fname, "%s/%s", iodir, infile_lostlep ) ;
     loadHist( fname,    "lostlep"    ) ;

     sprintf( fname, "%s/%s", iodir, infile_hadtau ) ;
     loadHist( fname,    "hadtau"    ) ;

     sprintf( fname, "%s/%s", iodir, infile_znunu ) ;
     loadHist( fname,    "znunu"    ) ;

     sprintf( fname, "%s/%s", iodir, infile_qcd ) ;
     loadHist( fname,    "qcd"    ) ;


     printf("\n") ;
     gDirectory -> ls( "h_mhtc_nbsum_ldp*" ) ;
     printf("\n") ;
     gDirectory -> ls( "h_mhtc_nbsum_hdp*" ) ;
     printf("\n") ;

     TString line ;

     ifstream ifs_lostlep ;
     ifs_lostlep.open( "outputfiles/nbsum-stat-syst-lostlep.txt" ) ;
     if ( !ifs_lostlep.good() ) { printf("\n\n *** Can't find outputfiles/nbsum-stat-syst-lostlep.txt\n\n") ; return ; }

     ifstream ifs_hadtau ;
     ifs_hadtau.open( "outputfiles/nbsum-stat-syst-hadtau.txt" ) ;
     if ( !ifs_hadtau.good() ) { printf("\n\n *** Can't find outputfiles/nbsum-stat-syst-hadtau.txt\n\n") ; return ; }

     ifstream ifs_znunu ;
     ifs_znunu.open( "outputfiles/nbsum-stat-syst-znunu.txt" ) ;
     if ( !ifs_znunu.good() ) { printf("\n\n *** Can't find outputfiles/nbsum-stat-syst-znunu.txt\n\n") ; return ; }



     int nb_nj(4) ;
     int nb_ht(3) ;



   //-----

     TH1F* h_ldp_lostlep = get_hist( "h_mhtc_nbsum_ldp_lostlep" ) ;
     TH1F* h_ldp_hadtau  = get_hist( "h_mhtc_nbsum_ldp_hadtau" ) ;
     TH1F* h_ldp_znunu   = get_hist( "h_mhtc_nbsum_ldp_znunu" ) ;
     TH1F* h_ldp_qcd     = get_hist( "h_mhtc_nbsum_ldp_qcd" ) ;

     TH1F* h_ldp_expected_error_lostlep = (TH1F*) h_ldp_lostlep -> Clone( "h_ldp_expected_error_lostlep" ) ;
     TH1F* h_ldp_expected_error_hadtau = (TH1F*) h_ldp_hadtau -> Clone( "h_ldp_expected_error_hadtau" ) ;
     TH1F* h_ldp_expected_error_znunu = (TH1F*) h_ldp_znunu -> Clone( "h_ldp_expected_error_znunu" ) ;
     TH1F* h_ldp_expected_error_qcd = (TH1F*) h_ldp_qcd -> Clone( "h_ldp_expected_error_qcd" ) ;

     TH1F* h_ldp_expected_error_htgrouping_lostlep = new TH1F( "h_ldp_expected_error_htgrouping_lostlep", "LDP, lostlep", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_ldp_expected_error_htgrouping_hadtau = new TH1F( "h_ldp_expected_error_htgrouping_hadtau", "LDP, hadtau", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_ldp_expected_error_htgrouping_znunu = new TH1F( "h_ldp_expected_error_htgrouping_znunu", "LDP, znunu", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_ldp_expected_error_htgrouping_qcd = new TH1F( "h_ldp_expected_error_htgrouping_qcd", "LDP, qcd", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;


   //-----

     TH1F* h_hdp_lostlep = get_hist( "h_mhtc_nbsum_hdp_lostlep" ) ;
     TH1F* h_hdp_hadtau  = get_hist( "h_mhtc_nbsum_hdp_hadtau" ) ;
     TH1F* h_hdp_znunu   = get_hist( "h_mhtc_nbsum_hdp_znunu" ) ;
     TH1F* h_hdp_qcd     = get_hist( "h_mhtc_nbsum_hdp_qcd" ) ;

     TH1F* h_hdp_expected_error_lostlep = (TH1F*) h_hdp_lostlep -> Clone( "h_hdp_expected_error_lostlep" ) ;
     TH1F* h_hdp_expected_error_hadtau = (TH1F*) h_hdp_hadtau -> Clone( "h_hdp_expected_error_hadtau" ) ;
     TH1F* h_hdp_expected_error_znunu = (TH1F*) h_hdp_znunu -> Clone( "h_hdp_expected_error_znunu" ) ;
     TH1F* h_hdp_expected_error_qcd = (TH1F*) h_hdp_qcd -> Clone( "h_hdp_expected_error_qcd" ) ;

     TH1F* h_hdp_expected_error_htgrouping_lostlep = new TH1F( "h_hdp_expected_error_htgrouping_lostlep", "hdp, lostlep", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_hdp_expected_error_htgrouping_hadtau = new TH1F( "h_hdp_expected_error_htgrouping_hadtau", "hdp, hadtau", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_hdp_expected_error_htgrouping_znunu = new TH1F( "h_hdp_expected_error_htgrouping_znunu", "hdp, znunu", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_hdp_expected_error_htgrouping_qcd = new TH1F( "h_hdp_expected_error_htgrouping_qcd", "hdp, qcd", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;

     TH1F* h_qcd_ratio_htgrouping = new TH1F( "h_qcd_ratio_htgrouping", "QCD H/L ratio", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;
     TH1F* h_qcd_ratio_rel_error_htgrouping = new TH1F( "h_qcd_ratio_rel_error_htgrouping", "QCD H/L ratio (err/val)", nb_nj*nb_ht, 0.5, nb_nj*nb_ht+0.5 ) ;


   //-----

     printf("\n\n") ;

     for ( int bi=1; bi<=h_ldp_lostlep->GetNbinsX(); bi++ ) {
        char label[100] ;
        sprintf( label, "%s", h_ldp_lostlep->GetXaxis()->GetBinLabel(bi) ) ;
        printf(" %2d : %s :\n", bi, label ) ;
     } // bi

     for ( int bi_ht=1; bi_ht<=nb_ht; bi_ht++ ) {

        for ( int bi_nj=1; bi_nj<=nb_nj; bi_nj++ ) {


           int bi_hist = (bi_nj-1)*(nb_ht) + bi_ht ;
           char label[100] ;
           char label2[100] ;
           sprintf( label, "%s", h_ldp_lostlep->GetXaxis()->GetBinLabel(bi_hist) ) ;

           int bi_hist_htgrouping = (bi_ht-1)*(nb_nj) + bi_nj ;
           char label_htgrouping[100] ;
           sprintf( label_htgrouping, "HT%d-Nj%d", bi_ht, bi_nj ) ;


           double ldp_lostlep = h_ldp_lostlep -> GetBinContent( bi_hist ) ;
           double ldp_hadtau = h_ldp_hadtau -> GetBinContent( bi_hist ) ;
           double ldp_znunu = h_ldp_znunu -> GetBinContent( bi_hist ) ;
           double ldp_qcd = h_ldp_qcd -> GetBinContent( bi_hist ) ;

           double hdp_lostlep = h_hdp_lostlep -> GetBinContent( bi_hist ) ;
           double hdp_hadtau = h_hdp_hadtau -> GetBinContent( bi_hist ) ;
           double hdp_znunu = h_hdp_znunu -> GetBinContent( bi_hist ) ;
           double hdp_qcd = h_hdp_qcd -> GetBinContent( bi_hist ) ;




           double ldp_stat_over_sqrtn_lostlep(0.), ldp_syst_over_n_lostlep(0.), hdp_stat_over_sqrtn_lostlep(0.), hdp_syst_over_n_lostlep(0.) ;
           line.ReadLine( ifs_lostlep ) ;
           sscanf( line.Data(), "%s %lf %lf  %lf %lf", label2, &ldp_stat_over_sqrtn_lostlep, &ldp_syst_over_n_lostlep,
             &hdp_stat_over_sqrtn_lostlep, &hdp_syst_over_n_lostlep ) ;

           double ldp_stat_over_sqrtn_hadtau(0.), ldp_syst_over_n_hadtau(0.), hdp_stat_over_sqrtn_hadtau(0.), hdp_syst_over_n_hadtau(0.) ;
           line.ReadLine( ifs_hadtau ) ;
           sscanf( line.Data(), "%s %lf %lf  %lf %lf", label2, &ldp_stat_over_sqrtn_hadtau, &ldp_syst_over_n_hadtau,
             &hdp_stat_over_sqrtn_hadtau, &hdp_syst_over_n_hadtau ) ;

           double ldp_stat_over_sqrtn_znunu(0.), ldp_syst_over_n_znunu(0.), hdp_stat_over_sqrtn_znunu(0.), hdp_syst_over_n_znunu(0.) ;
           line.ReadLine( ifs_znunu ) ;
           sscanf( line.Data(), "%s %lf %lf  %lf %lf", label2, &ldp_stat_over_sqrtn_znunu, &ldp_syst_over_n_znunu,
             &hdp_stat_over_sqrtn_znunu, &hdp_syst_over_n_znunu ) ;


           double sum_ldp_allmc  = ldp_lostlep + ldp_hadtau + ldp_znunu + ldp_qcd ;
           double sum_ldp_nonqcd = ldp_lostlep + ldp_hadtau + ldp_znunu ;

           double sum_hdp_allmc  = hdp_lostlep + hdp_hadtau + hdp_znunu + hdp_qcd ;
           double sum_hdp_nonqcd = hdp_lostlep + hdp_hadtau + hdp_znunu ;



           double ldp_lostlep_err(0.) ;
           if ( ldp_lostlep > 0 ) ldp_lostlep_err = sqrt( pow( ldp_stat_over_sqrtn_lostlep, 2.) * ldp_lostlep + pow( ldp_syst_over_n_lostlep * ldp_lostlep, 2.) ) ;
           double hdp_lostlep_err(0.) ;
           if ( hdp_lostlep > 0 ) hdp_lostlep_err = sqrt( pow( hdp_stat_over_sqrtn_lostlep, 2.) * hdp_lostlep + pow( hdp_syst_over_n_lostlep * hdp_lostlep, 2.) ) ;
           printf("  %s  : lostlep :  LDP %6.1f +/- %4.1f +/- %4.1f      HDP %6.1f +/- %4.1f +/- %4.1f\n", label,
              ldp_lostlep, ldp_stat_over_sqrtn_lostlep * sqrt(ldp_lostlep), ldp_syst_over_n_lostlep * ldp_lostlep,
              hdp_lostlep, hdp_stat_over_sqrtn_lostlep * sqrt(hdp_lostlep), hdp_syst_over_n_lostlep * hdp_lostlep ) ;

           h_ldp_expected_error_lostlep -> SetBinError( bi_hist, ldp_lostlep_err ) ;
           h_hdp_expected_error_lostlep -> SetBinError( bi_hist, hdp_lostlep_err ) ;

           h_ldp_expected_error_htgrouping_lostlep -> SetBinContent( bi_hist_htgrouping, ldp_lostlep ) ;
           h_hdp_expected_error_htgrouping_lostlep -> SetBinContent( bi_hist_htgrouping, hdp_lostlep ) ;

           h_ldp_expected_error_htgrouping_lostlep -> SetBinError( bi_hist_htgrouping, ldp_lostlep_err ) ;
           h_hdp_expected_error_htgrouping_lostlep -> SetBinError( bi_hist_htgrouping, hdp_lostlep_err ) ;

           h_ldp_expected_error_htgrouping_lostlep -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;
           h_hdp_expected_error_htgrouping_lostlep -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;



           double ldp_hadtau_err(0.) ;
           if ( ldp_hadtau > 0 ) ldp_hadtau_err = sqrt( pow( ldp_stat_over_sqrtn_hadtau, 2.) * ldp_hadtau + pow( ldp_syst_over_n_hadtau * ldp_hadtau, 2.) ) ;
           double hdp_hadtau_err(0.) ;
           if ( hdp_hadtau > 0 ) hdp_hadtau_err = sqrt( pow( hdp_stat_over_sqrtn_hadtau, 2.) * hdp_hadtau + pow( hdp_syst_over_n_hadtau * hdp_hadtau, 2.) ) ;
           printf("  %s  : hadtau  :  LDP %6.1f +/- %4.1f +/- %4.1f      HDP %6.1f +/- %4.1f +/- %4.1f\n", label,
              ldp_hadtau, ldp_stat_over_sqrtn_hadtau * sqrt(ldp_hadtau), ldp_syst_over_n_hadtau * ldp_hadtau,
              hdp_hadtau, hdp_stat_over_sqrtn_hadtau * sqrt(hdp_hadtau), hdp_syst_over_n_hadtau * hdp_hadtau ) ;

           h_ldp_expected_error_hadtau -> SetBinError( bi_hist, ldp_hadtau_err ) ;
           h_hdp_expected_error_hadtau -> SetBinError( bi_hist, hdp_hadtau_err ) ;

           h_ldp_expected_error_htgrouping_hadtau -> SetBinContent( bi_hist_htgrouping, ldp_hadtau ) ;
           h_hdp_expected_error_htgrouping_hadtau -> SetBinContent( bi_hist_htgrouping, hdp_hadtau ) ;

           h_ldp_expected_error_htgrouping_hadtau -> SetBinError( bi_hist_htgrouping, ldp_hadtau_err ) ;
           h_hdp_expected_error_htgrouping_hadtau -> SetBinError( bi_hist_htgrouping, hdp_hadtau_err ) ;

           h_ldp_expected_error_htgrouping_hadtau -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;
           h_hdp_expected_error_htgrouping_hadtau -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;




           double ldp_znunu_err(0.) ;
           if ( ldp_znunu > 0 ) ldp_znunu_err = sqrt( pow( ldp_stat_over_sqrtn_znunu, 2.) * ldp_znunu + pow( ldp_syst_over_n_znunu * ldp_znunu, 2.) ) ;
           double hdp_znunu_err(0.) ;
           if ( hdp_znunu > 0 ) hdp_znunu_err = sqrt( pow( hdp_stat_over_sqrtn_znunu, 2.) * hdp_znunu + pow( hdp_syst_over_n_znunu * hdp_znunu, 2.) ) ;
           printf("  %s  : znunu   :  LDP %6.1f +/- %4.1f +/- %4.1f      HDP %6.1f +/- %4.1f +/- %4.1f\n", label,
              ldp_znunu, ldp_stat_over_sqrtn_znunu * sqrt(ldp_znunu), ldp_syst_over_n_znunu * ldp_znunu,
              hdp_znunu, hdp_stat_over_sqrtn_znunu * sqrt(hdp_znunu), hdp_syst_over_n_znunu * hdp_znunu ) ;

           h_ldp_expected_error_znunu -> SetBinError( bi_hist, ldp_znunu_err ) ;
           h_hdp_expected_error_znunu -> SetBinError( bi_hist, hdp_znunu_err ) ;

           h_ldp_expected_error_htgrouping_znunu -> SetBinContent( bi_hist_htgrouping, ldp_znunu ) ;
           h_hdp_expected_error_htgrouping_znunu -> SetBinContent( bi_hist_htgrouping, hdp_znunu ) ;

           h_ldp_expected_error_htgrouping_znunu -> SetBinError( bi_hist_htgrouping, ldp_znunu_err ) ;
           h_hdp_expected_error_htgrouping_znunu -> SetBinError( bi_hist_htgrouping, hdp_znunu_err ) ;

           h_ldp_expected_error_htgrouping_znunu -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;
           h_hdp_expected_error_htgrouping_znunu -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;




           double sum_ldp_nonqcd_err = sqrt( pow( ldp_lostlep_err, 2.) + pow( ldp_hadtau_err, 2.) + pow( ldp_znunu_err, 2.) ) ;
           double sum_hdp_nonqcd_err = sqrt( pow( hdp_lostlep_err, 2.) + pow( hdp_hadtau_err, 2.) + pow( hdp_znunu_err, 2.) ) ;


           ///// printf(" %2d : HT%d-Nj%d : %s %s:\n", bi_hist, bi_ht, bi_nj, label, label2 ) ;


           double calc_ldp_qcd_val = sum_ldp_allmc - sum_ldp_nonqcd ;
           double calc_ldp_qcd_err = sqrt( sum_ldp_allmc + pow(sum_ldp_nonqcd_err, 2.) ) ;

           double calc_hdp_qcd_val = sum_hdp_allmc - sum_hdp_nonqcd ;
           double calc_hdp_qcd_err = sqrt( sum_hdp_allmc + pow(sum_hdp_nonqcd_err, 2.) ) ;
           double calc_hdp_qcd_rel_err(0.) ;
           if ( calc_hdp_qcd_val > 0 ) calc_hdp_qcd_rel_err = calc_hdp_qcd_err / calc_hdp_qcd_val ;

           h_ldp_expected_error_qcd -> SetBinError( bi_hist, calc_ldp_qcd_err ) ;
           h_hdp_expected_error_qcd -> SetBinError( bi_hist, calc_hdp_qcd_err ) ;

           h_ldp_expected_error_htgrouping_qcd -> SetBinContent( bi_hist_htgrouping, calc_ldp_qcd_val ) ;
           h_hdp_expected_error_htgrouping_qcd -> SetBinContent( bi_hist_htgrouping, calc_hdp_qcd_val ) ;

           h_ldp_expected_error_htgrouping_qcd -> SetBinError( bi_hist_htgrouping, calc_ldp_qcd_err ) ;
           h_hdp_expected_error_htgrouping_qcd -> SetBinError( bi_hist_htgrouping, calc_hdp_qcd_err ) ;

           h_ldp_expected_error_htgrouping_qcd -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;
           h_hdp_expected_error_htgrouping_qcd -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;




           double ratio_val(0.) ;
           double ratio_err(0.) ;
           double ratio_rel_err(0.) ;
           if ( calc_ldp_qcd_val > 0 && calc_hdp_qcd_val > 0 ) {
              ratio_val = calc_hdp_qcd_val / calc_ldp_qcd_val ;
              ratio_err = ratio_val * sqrt( pow( calc_ldp_qcd_err/calc_ldp_qcd_val, 2. ) + pow( calc_hdp_qcd_err/calc_hdp_qcd_val, 2. ) ) ;
              ratio_rel_err = ratio_err / ratio_val ;
           }

           h_qcd_ratio_htgrouping -> SetBinContent( bi_hist_htgrouping, ratio_val ) ;
           h_qcd_ratio_htgrouping -> SetBinError( bi_hist_htgrouping, ratio_err ) ;

           h_qcd_ratio_htgrouping -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;


           h_qcd_ratio_rel_error_htgrouping -> SetBinContent( bi_hist_htgrouping, ratio_rel_err ) ;
           h_qcd_ratio_rel_error_htgrouping -> GetXaxis() -> SetBinLabel( bi_hist_htgrouping, label_htgrouping ) ;

//    ///  printf(" HT%d-Nj%d :  LDP all = %7.1f , non-QCD = %7.1f , calc QCD = %7.1f +/- %7.1f   ||  ",
//    ///    bi_ht, bi_nj, sum_ldp_allmc, sum_ldp_nonqcd, calc_ldp_qcd_val, calc_ldp_qcd_err ) ;
//    ///  printf(" HDP all = %7.1f , non-QCD = %7.1f , calc QCD = %7.1f +/- %7.1f\n",
//    ///     sum_hdp_allmc, sum_hdp_nonqcd, calc_hdp_qcd_val, calc_hdp_qcd_err ) ;

           printf(" HT%d-Nj%d :  LDP calc QCD = %7.1f +/- %7.1f   ||  ",
             bi_ht, bi_nj, calc_ldp_qcd_val, calc_ldp_qcd_err ) ;
           printf(" HDP calc QCD = %7.1f +/- %7.1f  (%6.3f)     ||  Ratio = %6.3f +/- %6.3f  (%6.3f)\n",
               calc_hdp_qcd_val, calc_hdp_qcd_err, calc_hdp_qcd_rel_err,  ratio_val, ratio_err, ratio_rel_err ) ;

        } // bi_nj

     } // bi_ht

     h_ldp_expected_error_htgrouping_lostlep -> GetXaxis() -> LabelsOption("v") ;
     h_ldp_expected_error_htgrouping_hadtau -> GetXaxis() -> LabelsOption("v") ;
     h_ldp_expected_error_htgrouping_znunu -> GetXaxis() -> LabelsOption("v") ;
     h_ldp_expected_error_htgrouping_qcd -> GetXaxis() -> LabelsOption("v") ;

     h_hdp_expected_error_htgrouping_lostlep -> GetXaxis() -> LabelsOption("v") ;
     h_hdp_expected_error_htgrouping_hadtau -> GetXaxis() -> LabelsOption("v") ;
     h_hdp_expected_error_htgrouping_znunu -> GetXaxis() -> LabelsOption("v") ;
     h_hdp_expected_error_htgrouping_qcd -> GetXaxis() -> LabelsOption("v") ;

     h_qcd_ratio_htgrouping -> GetXaxis() -> LabelsOption("v") ;
     h_qcd_ratio_rel_error_htgrouping -> GetXaxis() -> LabelsOption("v") ;

     sprintf( fname, "%s/calc-mc-ratios.root", iodir ) ;
     saveHist( fname, "h*" ) ;

     printf("\n\n") ;


   } // calc_mc_hlratios1


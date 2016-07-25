#define syst_2015_v2_cxx
#include "syst_2015_v2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TStopwatch.h"
#include "TSystem.h"

#include "histio.c"

double calc_dphi( double phi1, double phi2 ) ;

void syst_2015_v2::Loop( int max_dump, bool verb )
{
   double lumi = 2300 ;

   setup_bins() ;

   printf("\n\n") ;

   TH1F* h_dpt        = new TH1F( "h_dpt", "Pt(rec)-Pt(gen) for worst jet", 60, 0., 600. ) ;
   TH1F* h_dpt_withnu = new TH1F( "h_dpt_withnu", "Pt(rec)-Pt(gen) for worst jet, neutrinos added in to Pt(gen)", 60, 0., 600. ) ;
   h_dpt -> Sumw2() ;
   h_dpt_withnu -> Sumw2() ;

   TH1F* h_worst_ptrec_over_ptgen_all = new TH1F( "h_worst_ptrec_over_ptgen_all", "Pt(rec)/Pt(gen) for worst jet (largest DPt), all", 50, 0., 2.5 ) ;
   TH1F* h_worst_ptrec_over_ptgen_genmet_lt_60 = new TH1F( "h_worst_ptrec_over_ptgen_genmet_lt_60", "Pt(rec)/Pt(gen) for worst jet (largest DPt), GenMET<60", 50, 0., 2.5 ) ;
   TH1F* h_worst_ptrec_over_ptgen_genmet_gt_60 = new TH1F( "h_worst_ptrec_over_ptgen_genmet_gt_60", "Pt(rec)/Pt(gen) for worst jet (largest DPt), GenMET>60", 50, 0., 2.5 ) ;
   h_worst_ptrec_over_ptgen_all -> Sumw2() ;
   h_worst_ptrec_over_ptgen_genmet_lt_60 -> Sumw2() ;
   h_worst_ptrec_over_ptgen_genmet_gt_60 -> Sumw2() ;

  //-----

   TH1F* h_mdp_all = new TH1F( "h_mdp_all", "Min Delta Phi, all baseline", 68, 0., 3.4 ) ; h_mdp_all -> Sumw2() ;
   TH1F* h_mdp_mhtc = new TH1F( "h_mdp_mhtc", "Min Delta Phi, MHTC", 68, 0., 3.4 ) ;       h_mdp_mhtc -> Sumw2() ;
   TH1F* h_mdp_mht1 = new TH1F( "h_mdp_mht1", "Min Delta Phi, MHT1", 68, 0., 3.4 ) ;       h_mdp_mht1 -> Sumw2() ;
   TH1F* h_mdp_mht2 = new TH1F( "h_mdp_mht2", "Min Delta Phi, MHT2", 68, 0., 3.4 ) ;       h_mdp_mht2 -> Sumw2() ;
   TH1F* h_mdp_mht3 = new TH1F( "h_mdp_mht3", "Min Delta Phi, MHT3", 68, 0., 3.4 ) ;       h_mdp_mht3 -> Sumw2() ;
   TH1F* h_mdp_mht4 = new TH1F( "h_mdp_mht4", "Min Delta Phi, MHT4", 68, 0., 3.4 ) ;       h_mdp_mht4 -> Sumw2() ;

   TH1F* h_mdp_all_badj_in_dphi = new TH1F( "h_mdp_all_badj_in_dphi", "Min Delta Phi, all baseline, bad jet in Dphi", 68, 0., 3.4 ) ; h_mdp_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mhtc_badj_in_dphi = new TH1F( "h_mdp_mhtc_badj_in_dphi", "Min Delta Phi, MHTC, bad jet in Dphi", 68, 0., 3.4 ) ;       h_mdp_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht1_badj_in_dphi = new TH1F( "h_mdp_mht1_badj_in_dphi", "Min Delta Phi, MHT1, bad jet in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht2_badj_in_dphi = new TH1F( "h_mdp_mht2_badj_in_dphi", "Min Delta Phi, MHT2, bad jet in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht3_badj_in_dphi = new TH1F( "h_mdp_mht3_badj_in_dphi", "Min Delta Phi, MHT3, bad jet in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht4_badj_in_dphi = new TH1F( "h_mdp_mht4_badj_in_dphi", "Min Delta Phi, MHT4, bad jet in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_mdp_all_badj_not_in_dphi = new TH1F( "h_mdp_all_badj_not_in_dphi", "Min Delta Phi, all baseline, bad jet not in Dphi", 68, 0., 3.4 ) ; h_mdp_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mhtc_badj_not_in_dphi = new TH1F( "h_mdp_mhtc_badj_not_in_dphi", "Min Delta Phi, MHTC, bad jet not in Dphi", 68, 0., 3.4 ) ;       h_mdp_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht1_badj_not_in_dphi = new TH1F( "h_mdp_mht1_badj_not_in_dphi", "Min Delta Phi, MHT1, bad jet not in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht2_badj_not_in_dphi = new TH1F( "h_mdp_mht2_badj_not_in_dphi", "Min Delta Phi, MHT2, bad jet not in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht3_badj_not_in_dphi = new TH1F( "h_mdp_mht3_badj_not_in_dphi", "Min Delta Phi, MHT3, bad jet not in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_mht4_badj_not_in_dphi = new TH1F( "h_mdp_mht4_badj_not_in_dphi", "Min Delta Phi, MHT4, bad jet not in Dphi", 68, 0., 3.4 ) ;       h_mdp_mht4_badj_not_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_all = new TH1F( "h_dphiregion_all", "Dphi region (0=LDP, 1=HDP)", 2, -0.5, 1.5 ) ;            h_dphiregion_all -> Sumw2() ;
   TH1F* h_dphiregion_mhtc = new TH1F( "h_dphiregion_mhtc", "Dphi region (0=LDP, 1=HDP), MHTC", 2, -0.5, 1.5 ) ;    h_dphiregion_mhtc -> Sumw2() ;
   TH1F* h_dphiregion_mht1 = new TH1F( "h_dphiregion_mht1", "Dphi region (0=LDP, 1=HDP), MHT1", 2, -0.5, 1.5 ) ;    h_dphiregion_mht1 -> Sumw2() ;
   TH1F* h_dphiregion_mht2 = new TH1F( "h_dphiregion_mht2", "Dphi region (0=LDP, 1=HDP), MHT2", 2, -0.5, 1.5 ) ;    h_dphiregion_mht2 -> Sumw2() ;
   TH1F* h_dphiregion_mht3 = new TH1F( "h_dphiregion_mht3", "Dphi region (0=LDP, 1=HDP), MHT3", 2, -0.5, 1.5 ) ;    h_dphiregion_mht3 -> Sumw2() ;
   TH1F* h_dphiregion_mht4 = new TH1F( "h_dphiregion_mht4", "Dphi region (0=LDP, 1=HDP), MHT4", 2, -0.5, 1.5 ) ;    h_dphiregion_mht4 -> Sumw2() ;

   TH1F* h_dphiregion_all_badj_in_dphi = new TH1F( "h_dphiregion_all_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet in Dphi", 2, -0.5, 1.5 ) ;          h_dphiregion_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mhtc_badj_in_dphi = new TH1F( "h_dphiregion_mhtc_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht1_badj_in_dphi = new TH1F( "h_dphiregion_mht1_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht2_badj_in_dphi = new TH1F( "h_dphiregion_mht2_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht3_badj_in_dphi = new TH1F( "h_dphiregion_mht3_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht4_badj_in_dphi = new TH1F( "h_dphiregion_mht4_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_all_badj_not_in_dphi = new TH1F( "h_dphiregion_all_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet not in Dphi", 2, -0.5, 1.5 ) ;          h_dphiregion_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mhtc_badj_not_in_dphi = new TH1F( "h_dphiregion_mhtc_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet not in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht1_badj_not_in_dphi = new TH1F( "h_dphiregion_mht1_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet not in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht2_badj_not_in_dphi = new TH1F( "h_dphiregion_mht2_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet not in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht3_badj_not_in_dphi = new TH1F( "h_dphiregion_mht3_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet not in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_mht4_badj_not_in_dphi = new TH1F( "h_dphiregion_mht4_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet not in Dphi", 2, -0.5, 1.5 ) ;  h_dphiregion_mht4_badj_not_in_dphi -> Sumw2() ;


  //-----

   TH1F* h_mdp_hth_all = new TH1F( "h_mdp_hth_all", "Min Delta Phi, all baseline, high HT", 68, 0., 3.4 ) ; h_mdp_hth_all -> Sumw2() ;
   TH1F* h_mdp_hth_mhtc = new TH1F( "h_mdp_hth_mhtc", "Min Delta Phi, MHTC, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mhtc -> Sumw2() ;
   TH1F* h_mdp_hth_mht1 = new TH1F( "h_mdp_hth_mht1", "Min Delta Phi, MHT1, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht1 -> Sumw2() ;
   TH1F* h_mdp_hth_mht2 = new TH1F( "h_mdp_hth_mht2", "Min Delta Phi, MHT2, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht2 -> Sumw2() ;
   TH1F* h_mdp_hth_mht3 = new TH1F( "h_mdp_hth_mht3", "Min Delta Phi, MHT3, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht3 -> Sumw2() ;
   TH1F* h_mdp_hth_mht4 = new TH1F( "h_mdp_hth_mht4", "Min Delta Phi, MHT4, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht4 -> Sumw2() ;

   TH1F* h_mdp_hth_all_badj_in_dphi = new TH1F( "h_mdp_hth_all_badj_in_dphi", "Min Delta Phi, all baseline, bad jet in Dphi, high HT", 68, 0., 3.4 ) ; h_mdp_hth_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mhtc_badj_in_dphi = new TH1F( "h_mdp_hth_mhtc_badj_in_dphi", "Min Delta Phi, MHTC, bad jet in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht1_badj_in_dphi = new TH1F( "h_mdp_hth_mht1_badj_in_dphi", "Min Delta Phi, MHT1, bad jet in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht2_badj_in_dphi = new TH1F( "h_mdp_hth_mht2_badj_in_dphi", "Min Delta Phi, MHT2, bad jet in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht3_badj_in_dphi = new TH1F( "h_mdp_hth_mht3_badj_in_dphi", "Min Delta Phi, MHT3, bad jet in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht4_badj_in_dphi = new TH1F( "h_mdp_hth_mht4_badj_in_dphi", "Min Delta Phi, MHT4, bad jet in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_mdp_hth_all_badj_not_in_dphi = new TH1F( "h_mdp_hth_all_badj_not_in_dphi", "Min Delta Phi, all baseline, bad jet not in Dphi, high HT", 68, 0., 3.4 ) ; h_mdp_hth_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mhtc_badj_not_in_dphi = new TH1F( "h_mdp_hth_mhtc_badj_not_in_dphi", "Min Delta Phi, MHTC, bad jet not in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht1_badj_not_in_dphi = new TH1F( "h_mdp_hth_mht1_badj_not_in_dphi", "Min Delta Phi, MHT1, bad jet not in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht2_badj_not_in_dphi = new TH1F( "h_mdp_hth_mht2_badj_not_in_dphi", "Min Delta Phi, MHT2, bad jet not in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht3_badj_not_in_dphi = new TH1F( "h_mdp_hth_mht3_badj_not_in_dphi", "Min Delta Phi, MHT3, bad jet not in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_hth_mht4_badj_not_in_dphi = new TH1F( "h_mdp_hth_mht4_badj_not_in_dphi", "Min Delta Phi, MHT4, bad jet not in Dphi, high HT", 68, 0., 3.4 ) ;       h_mdp_hth_mht4_badj_not_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_hth_all = new TH1F( "h_dphiregion_hth_all", "Dphi region (0=LDP, 1=HDP), high HT", 2, -0.5, 1.5 ) ;            h_dphiregion_hth_all -> Sumw2() ;
   TH1F* h_dphiregion_hth_mhtc = new TH1F( "h_dphiregion_hth_mhtc", "Dphi region (0=LDP, 1=HDP), MHTC, high HT", 2, -0.5, 1.5 ) ;    h_dphiregion_hth_mhtc -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht1 = new TH1F( "h_dphiregion_hth_mht1", "Dphi region (0=LDP, 1=HDP), MHT1, high HT", 2, -0.5, 1.5 ) ;    h_dphiregion_hth_mht1 -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht2 = new TH1F( "h_dphiregion_hth_mht2", "Dphi region (0=LDP, 1=HDP), MHT2, high HT", 2, -0.5, 1.5 ) ;    h_dphiregion_hth_mht2 -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht3 = new TH1F( "h_dphiregion_hth_mht3", "Dphi region (0=LDP, 1=HDP), MHT3, high HT", 2, -0.5, 1.5 ) ;    h_dphiregion_hth_mht3 -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht4 = new TH1F( "h_dphiregion_hth_mht4", "Dphi region (0=LDP, 1=HDP), MHT4, high HT", 2, -0.5, 1.5 ) ;    h_dphiregion_hth_mht4 -> Sumw2() ;

   TH1F* h_dphiregion_hth_all_badj_in_dphi = new TH1F( "h_dphiregion_hth_all_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet in Dphi, high HT", 2, -0.5, 1.5 ) ;          h_dphiregion_hth_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mhtc_badj_in_dphi = new TH1F( "h_dphiregion_hth_mhtc_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht1_badj_in_dphi = new TH1F( "h_dphiregion_hth_mht1_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht2_badj_in_dphi = new TH1F( "h_dphiregion_hth_mht2_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht3_badj_in_dphi = new TH1F( "h_dphiregion_hth_mht3_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht4_badj_in_dphi = new TH1F( "h_dphiregion_hth_mht4_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_hth_all_badj_not_in_dphi = new TH1F( "h_dphiregion_hth_all_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet not in Dphi, high HT", 2, -0.5, 1.5 ) ;          h_dphiregion_hth_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mhtc_badj_not_in_dphi = new TH1F( "h_dphiregion_hth_mhtc_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet not in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht1_badj_not_in_dphi = new TH1F( "h_dphiregion_hth_mht1_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet not in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht2_badj_not_in_dphi = new TH1F( "h_dphiregion_hth_mht2_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet not in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht3_badj_not_in_dphi = new TH1F( "h_dphiregion_hth_mht3_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet not in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_hth_mht4_badj_not_in_dphi = new TH1F( "h_dphiregion_hth_mht4_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet not in Dphi, high HT", 2, -0.5, 1.5 ) ;  h_dphiregion_hth_mht4_badj_not_in_dphi -> Sumw2() ;



  //-----

   TH1F* h_mdp_htm_all = new TH1F( "h_mdp_htm_all", "Min Delta Phi, all baseline, medium HT", 68, 0., 3.4 ) ; h_mdp_htm_all -> Sumw2() ;
   TH1F* h_mdp_htm_mhtc = new TH1F( "h_mdp_htm_mhtc", "Min Delta Phi, MHTC, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mhtc -> Sumw2() ;
   TH1F* h_mdp_htm_mht1 = new TH1F( "h_mdp_htm_mht1", "Min Delta Phi, MHT1, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht1 -> Sumw2() ;
   TH1F* h_mdp_htm_mht2 = new TH1F( "h_mdp_htm_mht2", "Min Delta Phi, MHT2, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht2 -> Sumw2() ;
   TH1F* h_mdp_htm_mht3 = new TH1F( "h_mdp_htm_mht3", "Min Delta Phi, MHT3, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht3 -> Sumw2() ;
   TH1F* h_mdp_htm_mht4 = new TH1F( "h_mdp_htm_mht4", "Min Delta Phi, MHT4, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht4 -> Sumw2() ;

   TH1F* h_mdp_htm_all_badj_in_dphi = new TH1F( "h_mdp_htm_all_badj_in_dphi", "Min Delta Phi, all baseline, bad jet in Dphi, medium HT", 68, 0., 3.4 ) ; h_mdp_htm_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mhtc_badj_in_dphi = new TH1F( "h_mdp_htm_mhtc_badj_in_dphi", "Min Delta Phi, MHTC, bad jet in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht1_badj_in_dphi = new TH1F( "h_mdp_htm_mht1_badj_in_dphi", "Min Delta Phi, MHT1, bad jet in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht2_badj_in_dphi = new TH1F( "h_mdp_htm_mht2_badj_in_dphi", "Min Delta Phi, MHT2, bad jet in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht3_badj_in_dphi = new TH1F( "h_mdp_htm_mht3_badj_in_dphi", "Min Delta Phi, MHT3, bad jet in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht4_badj_in_dphi = new TH1F( "h_mdp_htm_mht4_badj_in_dphi", "Min Delta Phi, MHT4, bad jet in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_mdp_htm_all_badj_not_in_dphi = new TH1F( "h_mdp_htm_all_badj_not_in_dphi", "Min Delta Phi, all baseline, bad jet not in Dphi, medium HT", 68, 0., 3.4 ) ; h_mdp_htm_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mhtc_badj_not_in_dphi = new TH1F( "h_mdp_htm_mhtc_badj_not_in_dphi", "Min Delta Phi, MHTC, bad jet not in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht1_badj_not_in_dphi = new TH1F( "h_mdp_htm_mht1_badj_not_in_dphi", "Min Delta Phi, MHT1, bad jet not in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht2_badj_not_in_dphi = new TH1F( "h_mdp_htm_mht2_badj_not_in_dphi", "Min Delta Phi, MHT2, bad jet not in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht3_badj_not_in_dphi = new TH1F( "h_mdp_htm_mht3_badj_not_in_dphi", "Min Delta Phi, MHT3, bad jet not in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htm_mht4_badj_not_in_dphi = new TH1F( "h_mdp_htm_mht4_badj_not_in_dphi", "Min Delta Phi, MHT4, bad jet not in Dphi, medium HT", 68, 0., 3.4 ) ;       h_mdp_htm_mht4_badj_not_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_htm_all = new TH1F( "h_dphiregion_htm_all", "Dphi region (0=LDP, 1=HDP), medium HT", 2, -0.5, 1.5 ) ;            h_dphiregion_htm_all -> Sumw2() ;
   TH1F* h_dphiregion_htm_mhtc = new TH1F( "h_dphiregion_htm_mhtc", "Dphi region (0=LDP, 1=HDP), MHTC, medium HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htm_mhtc -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht1 = new TH1F( "h_dphiregion_htm_mht1", "Dphi region (0=LDP, 1=HDP), MHT1, medium HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htm_mht1 -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht2 = new TH1F( "h_dphiregion_htm_mht2", "Dphi region (0=LDP, 1=HDP), MHT2, medium HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htm_mht2 -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht3 = new TH1F( "h_dphiregion_htm_mht3", "Dphi region (0=LDP, 1=HDP), MHT3, medium HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htm_mht3 -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht4 = new TH1F( "h_dphiregion_htm_mht4", "Dphi region (0=LDP, 1=HDP), MHT4, medium HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htm_mht4 -> Sumw2() ;

   TH1F* h_dphiregion_htm_all_badj_in_dphi = new TH1F( "h_dphiregion_htm_all_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet in Dphi, medium HT", 2, -0.5, 1.5 ) ;          h_dphiregion_htm_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mhtc_badj_in_dphi = new TH1F( "h_dphiregion_htm_mhtc_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht1_badj_in_dphi = new TH1F( "h_dphiregion_htm_mht1_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht2_badj_in_dphi = new TH1F( "h_dphiregion_htm_mht2_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht3_badj_in_dphi = new TH1F( "h_dphiregion_htm_mht3_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht4_badj_in_dphi = new TH1F( "h_dphiregion_htm_mht4_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_htm_all_badj_not_in_dphi = new TH1F( "h_dphiregion_htm_all_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet not in Dphi, medium HT", 2, -0.5, 1.5 ) ;          h_dphiregion_htm_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mhtc_badj_not_in_dphi = new TH1F( "h_dphiregion_htm_mhtc_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet not in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht1_badj_not_in_dphi = new TH1F( "h_dphiregion_htm_mht1_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet not in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht2_badj_not_in_dphi = new TH1F( "h_dphiregion_htm_mht2_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet not in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht3_badj_not_in_dphi = new TH1F( "h_dphiregion_htm_mht3_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet not in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htm_mht4_badj_not_in_dphi = new TH1F( "h_dphiregion_htm_mht4_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet not in Dphi, medium HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htm_mht4_badj_not_in_dphi -> Sumw2() ;



  //-----

   TH1F* h_mdp_htl_all = new TH1F( "h_mdp_htl_all", "Min Delta Phi, all baseline, low HT", 68, 0., 3.4 ) ; h_mdp_htl_all -> Sumw2() ;
   TH1F* h_mdp_htl_mhtc = new TH1F( "h_mdp_htl_mhtc", "Min Delta Phi, MHTC, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mhtc -> Sumw2() ;
   TH1F* h_mdp_htl_mht1 = new TH1F( "h_mdp_htl_mht1", "Min Delta Phi, MHT1, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht1 -> Sumw2() ;
   TH1F* h_mdp_htl_mht2 = new TH1F( "h_mdp_htl_mht2", "Min Delta Phi, MHT2, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht2 -> Sumw2() ;
   TH1F* h_mdp_htl_mht3 = new TH1F( "h_mdp_htl_mht3", "Min Delta Phi, MHT3, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht3 -> Sumw2() ;
   TH1F* h_mdp_htl_mht4 = new TH1F( "h_mdp_htl_mht4", "Min Delta Phi, MHT4, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht4 -> Sumw2() ;

   TH1F* h_mdp_htl_all_badj_in_dphi = new TH1F( "h_mdp_htl_all_badj_in_dphi", "Min Delta Phi, all baseline, bad jet in Dphi, low HT", 68, 0., 3.4 ) ; h_mdp_htl_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mhtc_badj_in_dphi = new TH1F( "h_mdp_htl_mhtc_badj_in_dphi", "Min Delta Phi, MHTC, bad jet in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht1_badj_in_dphi = new TH1F( "h_mdp_htl_mht1_badj_in_dphi", "Min Delta Phi, MHT1, bad jet in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht2_badj_in_dphi = new TH1F( "h_mdp_htl_mht2_badj_in_dphi", "Min Delta Phi, MHT2, bad jet in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht3_badj_in_dphi = new TH1F( "h_mdp_htl_mht3_badj_in_dphi", "Min Delta Phi, MHT3, bad jet in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht4_badj_in_dphi = new TH1F( "h_mdp_htl_mht4_badj_in_dphi", "Min Delta Phi, MHT4, bad jet in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_mdp_htl_all_badj_not_in_dphi = new TH1F( "h_mdp_htl_all_badj_not_in_dphi", "Min Delta Phi, all baseline, bad jet not in Dphi, low HT", 68, 0., 3.4 ) ; h_mdp_htl_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mhtc_badj_not_in_dphi = new TH1F( "h_mdp_htl_mhtc_badj_not_in_dphi", "Min Delta Phi, MHTC, bad jet not in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht1_badj_not_in_dphi = new TH1F( "h_mdp_htl_mht1_badj_not_in_dphi", "Min Delta Phi, MHT1, bad jet not in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht2_badj_not_in_dphi = new TH1F( "h_mdp_htl_mht2_badj_not_in_dphi", "Min Delta Phi, MHT2, bad jet not in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht3_badj_not_in_dphi = new TH1F( "h_mdp_htl_mht3_badj_not_in_dphi", "Min Delta Phi, MHT3, bad jet not in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_mdp_htl_mht4_badj_not_in_dphi = new TH1F( "h_mdp_htl_mht4_badj_not_in_dphi", "Min Delta Phi, MHT4, bad jet not in Dphi, low HT", 68, 0., 3.4 ) ;       h_mdp_htl_mht4_badj_not_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_htl_all = new TH1F( "h_dphiregion_htl_all", "Dphi region (0=LDP, 1=HDP), low HT", 2, -0.5, 1.5 ) ;            h_dphiregion_htl_all -> Sumw2() ;
   TH1F* h_dphiregion_htl_mhtc = new TH1F( "h_dphiregion_htl_mhtc", "Dphi region (0=LDP, 1=HDP), MHTC, low HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htl_mhtc -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht1 = new TH1F( "h_dphiregion_htl_mht1", "Dphi region (0=LDP, 1=HDP), MHT1, low HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htl_mht1 -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht2 = new TH1F( "h_dphiregion_htl_mht2", "Dphi region (0=LDP, 1=HDP), MHT2, low HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htl_mht2 -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht3 = new TH1F( "h_dphiregion_htl_mht3", "Dphi region (0=LDP, 1=HDP), MHT3, low HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htl_mht3 -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht4 = new TH1F( "h_dphiregion_htl_mht4", "Dphi region (0=LDP, 1=HDP), MHT4, low HT", 2, -0.5, 1.5 ) ;    h_dphiregion_htl_mht4 -> Sumw2() ;

   TH1F* h_dphiregion_htl_all_badj_in_dphi = new TH1F( "h_dphiregion_htl_all_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet in Dphi, low HT", 2, -0.5, 1.5 ) ;          h_dphiregion_htl_all_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mhtc_badj_in_dphi = new TH1F( "h_dphiregion_htl_mhtc_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mhtc_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht1_badj_in_dphi = new TH1F( "h_dphiregion_htl_mht1_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht1_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht2_badj_in_dphi = new TH1F( "h_dphiregion_htl_mht2_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht2_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht3_badj_in_dphi = new TH1F( "h_dphiregion_htl_mht3_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht3_badj_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht4_badj_in_dphi = new TH1F( "h_dphiregion_htl_mht4_badj_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht4_badj_in_dphi -> Sumw2() ;

   TH1F* h_dphiregion_htl_all_badj_not_in_dphi = new TH1F( "h_dphiregion_htl_all_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), bad jet not in Dphi, low HT", 2, -0.5, 1.5 ) ;          h_dphiregion_htl_all_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mhtc_badj_not_in_dphi = new TH1F( "h_dphiregion_htl_mhtc_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHTC, bad jet not in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mhtc_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht1_badj_not_in_dphi = new TH1F( "h_dphiregion_htl_mht1_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT1, bad jet not in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht1_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht2_badj_not_in_dphi = new TH1F( "h_dphiregion_htl_mht2_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT2, bad jet not in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht2_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht3_badj_not_in_dphi = new TH1F( "h_dphiregion_htl_mht3_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT3, bad jet not in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht3_badj_not_in_dphi -> Sumw2() ;
   TH1F* h_dphiregion_htl_mht4_badj_not_in_dphi = new TH1F( "h_dphiregion_htl_mht4_badj_not_in_dphi", "Dphi region (0=LDP, 1=HDP), MHT4, bad jet not in Dphi, low HT", 2, -0.5, 1.5 ) ;  h_dphiregion_htl_mht4_badj_not_in_dphi -> Sumw2() ;



  //-------

   TH1F* h_hdp = new TH1F( "h_hdp", "HDP events", nb_global, 0.5, nb_global + 0.5 ) ; h_hdp -> Sumw2() ;
   TH1F* h_nbsum_hdp = new TH1F( "h_nbsum_hdp", "HDP events", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nbsum_hdp -> Sumw2() ;
   set_bin_labels( h_nbsum_hdp ) ;
   TH1F* h_nb0_hdp = new TH1F( "h_nb0_hdp", "HDP events, Nb=0", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb0_hdp -> Sumw2() ;
   set_bin_labels( h_nb0_hdp ) ;
   TH1F* h_nb1_hdp = new TH1F( "h_nb1_hdp", "HDP events, Nb=1", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb1_hdp -> Sumw2() ;
   set_bin_labels( h_nb1_hdp ) ;
   TH1F* h_nb2_hdp = new TH1F( "h_nb2_hdp", "HDP events, Nb=2", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb2_hdp -> Sumw2() ;
   set_bin_labels( h_nb2_hdp ) ;
   TH1F* h_nb3_hdp = new TH1F( "h_nb3_hdp", "HDP events, Nb=3", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb3_hdp -> Sumw2() ;
   set_bin_labels( h_nb3_hdp ) ;

   TH1F* h_ldp = new TH1F( "h_ldp", "ldp events", nb_global, 0.5, nb_global + 0.5 ) ; h_ldp -> Sumw2() ;
   TH1F* h_nbsum_ldp = new TH1F( "h_nbsum_ldp", "LDP events", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nbsum_ldp -> Sumw2() ;
   set_bin_labels( h_nbsum_ldp ) ;
   TH1F* h_nb0_ldp = new TH1F( "h_nb0_ldp", "ldp events, Nb=0", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb0_ldp -> Sumw2() ;
   set_bin_labels( h_nb0_ldp ) ;
   TH1F* h_nb1_ldp = new TH1F( "h_nb1_ldp", "ldp events, Nb=1", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb1_ldp -> Sumw2() ;
   set_bin_labels( h_nb1_ldp ) ;
   TH1F* h_nb2_ldp = new TH1F( "h_nb2_ldp", "ldp events, Nb=2", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb2_ldp -> Sumw2() ;
   set_bin_labels( h_nb2_ldp ) ;
   TH1F* h_nb3_ldp = new TH1F( "h_nb3_ldp", "ldp events, Nb=3", nb_global/nb_nb, 0.5, nb_global/nb_nb + 0.5 ) ; h_nb3_ldp -> Sumw2() ;
   set_bin_labels( h_nb3_ldp ) ;

   TH2F* h_mht_vs_ht_nominal = new TH2F( "h_mht_vs_ht_nominal", "MHT vs HT, nominal, Njets>=3", 100, 0., 2000.,  100, 0., 2000. ) ;

  //-------


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t n_dump = 0 ;

   TStopwatch sw ;
   sw.Start() ;
   int time(0) ;
   float projected_remaining(999999.) ;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if ( jentry%1000 == 0 ) { // timer printing stuff
         int thistime = sw.RealTime() ;
         sw.Continue() ;
         if ( thistime < 2 ) {
            printf("   %10llu out of %10llu  (%6.1f%%) \r", jentry, nentries, 100.*jentry/(1.*nentries) ) ;
         } else {
            if ( thistime > time ) projected_remaining = (1.*thistime)/(1.*jentry)*(nentries-jentry) ;
            if ( projected_remaining < 100 ) {
               printf("   %10llu out of %10llu  (%6.1f%%)    seconds remaining %4.0f                       \r", jentry, nentries, 100.*jentry/(1.*nentries), projected_remaining ) ;
            } else if ( projected_remaining < 3600 ) {
               printf("   %10llu out of %10llu  (%6.1f%%)    time remaining     %2d:%02d   \r", jentry, nentries, 100.*jentry/(1.*nentries),
                    TMath::Nint(projected_remaining)/60, TMath::Nint(projected_remaining)%60 ) ;
            } else {
               printf("   %10llu out of %10llu  (%6.1f%%)    time remaining  %2d:%02d:%02d   \r", jentry, nentries, 100.*jentry/(1.*nentries),
                    TMath::Nint(projected_remaining)/3600, (TMath::Nint(projected_remaining)%3600)/60, TMath::Nint(projected_remaining)%60 ) ;
            }
         }
         fflush(stdout) ;
         time = thistime ;
      } // timer printing stuff


      Long64_t ientry = LoadTree(jentry);

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double hw = Weight * lumi ;


     //--- apply new baseline.
      if ( NJets < 3 ) continue ;
      if ( MHT < 250. ) continue ;
      if ( HT < 300. ) continue ;


     //-- take out the trash
      bool badMuonEvent(false) ;
      bool isjunk(false) ;
      for ( unsigned long ji=0; ji<Jets->size(); ji++ ) {
         if ( Jets->at(ji).Pt() < 200 ) continue ;
         if ( Jets_muonEnergyFraction->at(ji) < 0.5 ) continue ;
         double dPhi = Jets->at(ji).Phi() - METPhi ;
         if ( dPhi >  3.1415926 ) dPhi = dPhi - 2*3.14159 ;
         if ( dPhi < -3.1415926 ) dPhi = dPhi + 2*3.14159 ;
         if ( fabs( dPhi ) > 3.1415926 - 0.40 ) { 
            badMuonEvent = true ;
            break ;
         }
      } // ji
      if ( badMuonEvent ) {
         if (verb) printf("\n\n *** Bad Muon event: %lld\n", jentry ) ;
         if (verb) printf("    MET = %7.1f  ,  MHT = %7.1f ,  CaloMET = %7.1f ,   METPhi = %6.3f,   MHTPhi = %6.3f\n", MET, MHT, CaloMET, METPhi, MHTPhi ) ;
         for ( unsigned long ji=0; ji<Jets->size(); ji++ ) {
            if ( Jets->at(ji).Pt() < 30 ) break ;
            if (verb) printf(" jet %2lu :  Pt = %7.1f,  phi = %6.3f,  eta = %6.3f,  Muon fr = %5.3f\n",
                 ji, Jets->at(ji).Pt(), Jets->at(ji).Phi(), Jets->at(ji).Eta(), Jets_muonEnergyFraction->at(ji) ) ;
            if (verb) fflush(stdout) ;
         } // ji
      } // badMuonEvent?

      if ( CaloMET <= 0 )  isjunk = true ;
      if ( eeBadScFilter < 1 )  isjunk = true ;
      if ( CSCTightHaloFilter < 1 )  isjunk = true ;
      if ( HBHEIsoNoiseFilter < 1 )  isjunk = true ;
      if ( HBHENoiseFilter < 1 )  isjunk = true ;
      if ( EcalDeadCellTriggerPrimitiveFilter < 1 )  isjunk = true ; // **** new for v9
      if ( NVtx <= 0 )  isjunk = true ;  // **** new for v9
      if ( ! JetID ) isjunk = true ;  // **** new for v9

      if ( badMuonEvent ) isjunk = true ;
      /////if ( MET/CaloMET > 5 )  isjunk = true ;
      ///////////if ( ! noMuonJet ) isjunk = true ;  // **** new for v9, should be redundant with above
      if ( PFCaloMETRatio > 5 ) isjunk = true ;  // **** new for v9, should be redundant with above

      if ( isjunk ) continue ;
      if ( MET<100 ) continue ;
      if ( CaloMET<80 ) continue ;

      if ( verb ) printf("\n\n\n\n ===== RunNum %6d, LumiBlockNum %5d, EvtNum %9llu\n", RunNum, LumiBlockNum, EvtNum ) ;
      if ( verb ) printf("   Njets %2d,  HT %5.0f,  MHT %5.0f,  MHT phi %6.3f\n", NJets, HT, MHT, MHTPhi ) ;
      if ( verb ) printf("   genHT %5.0f, GenMHT %5.0f, GenMET %5.0f\n", GenHT, GenMHT, GenMET ) ;




      double sum_jetpx(0.) ;
      double sum_jetpy(0.) ;
      int recalc_njets(0) ;
      double recalc_ht(0.) ;
      double recalc_dphi1(10.) ;
      double recalc_dphi2(10.) ;
      double recalc_dphi3(10.) ;
      double recalc_dphi4(10.) ;
      double rjet_ind_dphi1(-1) ;
      double rjet_ind_dphi2(-1) ;
      double rjet_ind_dphi3(-1) ;
      double rjet_ind_dphi4(-1) ;
      for ( unsigned int ji=0; ji<Jets->size(); ji++ ) {

         if ( Jets->at(ji).Pt() < 20 ) break ;
         double dphi = calc_dphi( MHTPhi, Jets->at(ji).Phi() ) ;

         if ( Jets->at(ji).Pt() > 30 ) {
            sum_jetpx += Jets->at(ji).Px() ;
            sum_jetpy += Jets->at(ji).Py() ;
            if ( fabs( Jets->at(ji).Eta() ) <= 2.4 ) {
               recalc_njets++ ;
               recalc_ht += Jets->at(ji).Pt() ;
               if ( recalc_dphi1 > 9 ) {
                  recalc_dphi1 = fabs( dphi ) ;
                  rjet_ind_dphi1 = ji ;
               } else if ( recalc_dphi2 > 9 ) {
                  recalc_dphi2 = fabs( dphi ) ;
                  rjet_ind_dphi2 = ji ;
               } else if ( recalc_dphi3 > 9 ) {
                  recalc_dphi3 = fabs( dphi ) ;
                  rjet_ind_dphi3 = ji ;
               } else if ( recalc_dphi4 > 9 ) {
                  recalc_dphi4 = fabs( dphi ) ;
                  rjet_ind_dphi4 = ji ;
               }
            }
         }

         char tagdphi[3] ;
         char tagpt[3] ;
         if ( Jets->at(ji).Pt() < 30 ) { sprintf( tagpt, "*" ) ; } else { sprintf( tagpt, " " ) ; }
         if ( fabs(dphi) < 0.5 ) { sprintf( tagdphi, "*" ) ; } else { sprintf( tagdphi, " " ) ; }
         char tageta[3] ;
         if ( fabs(Jets->at(ji).Eta()) > 2.4 ) { sprintf( tageta, "*" ) ; } else { sprintf( tageta, " " ) ; }
         if ( verb ) printf("  Jet %3d :  Pt = %5.0f %s,  eta = %6.3f %s,  phi = %6.3f,",
            ji, Jets->at(ji).Pt(), tagpt, Jets->at(ji).Eta(), tageta, Jets->at(ji).Phi() ) ;
         double rjdp(-1.) ;
         if ( verb ) printf("  DeltaPhi = %6.3f ", dphi ) ;
         if ( verb ) printf("%s", tagdphi ) ;

         if ( verb ) printf(" mu mult,E_fr = %d, %.2f ", Jets_muonMultiplicity->at(ji), Jets_muonEnergyFraction->at(ji) ) ;
         if ( verb ) printf(" Neutral E fracs (em,gamma) = %.2f,%.2f", Jets_neutralEmEnergyFraction->at(ji),  Jets_photonEnergyFraction->at(ji) ) ;
         if ( verb ) printf(" Charged E fracs (em,had) = %.2f,%.2f", Jets_chargedEmEnergyFraction->at(ji), Jets_chargedHadronEnergyFraction->at(ji) ) ;


         if ( verb ) printf("\n") ;

      } // ji


      double recalc_mht = sqrt( sum_jetpx*sum_jetpx + sum_jetpy*sum_jetpy ) ;
      double recalc_mht_phi = atan2( -1.*sum_jetpy, -1.*sum_jetpx ) ;

      double min_delta_phi(99.) ;
      if ( recalc_dphi1 < min_delta_phi ) min_delta_phi = recalc_dphi1 ;
      if ( recalc_dphi2 < min_delta_phi ) min_delta_phi = recalc_dphi2 ;
      if ( recalc_dphi3 < min_delta_phi ) min_delta_phi = recalc_dphi3 ;
      if ( recalc_dphi4 < min_delta_phi ) min_delta_phi = recalc_dphi4 ;




      if ( verb ) printf("  MHT phi calc check: recalc = %6.3f,  ntuple = %6.3f\n", recalc_mht_phi, MHTPhi ) ;
      if ( verb ) printf("  MHT calc check:    recalc = %6.1f,  ntuple = %6.1f\n", recalc_mht, MHT ) ;
      if ( verb ) printf("  HT calc check:     recalc = %6.1f,  ntuple = %6.1f\n", recalc_ht, HT ) ;
      if ( verb ) printf("  NJets calc check:  recalc = %2d,  ntuple = %2d\n", recalc_njets, NJets ) ;
      if ( verb ) printf("  Dphi1 calc check:  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi1, DeltaPhi1 ) ;
      if ( verb ) printf("  Dphi2 calc check:  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi2, DeltaPhi2 ) ;
      if ( verb ) printf("  Dphi3 calc check:  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi3, DeltaPhi3 ) ;
      if ( verb && NJets>3 ) printf("  Dphi4 calc check:  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi4, DeltaPhi4 ) ;
      if ( fabs( recalc_mht_phi - MHTPhi ) > 0.01 ) printf("\n\n *** Recalc of MHT phi failed?  recalc = %6.3f,  ntuple = %6.3f\n", recalc_mht_phi, MHTPhi ) ;
      if ( fabs( recalc_mht - MHT ) > 5. ) printf("\n\n *** Recalc of MHT failed?  recalc = %6.1f,  ntuple = %6.1f\n\n", recalc_mht, MHT ) ;
      if ( fabs( recalc_ht - HT ) > 5. ) printf("\n\n *** Recalc of HT failed?  recalc = %6.1f,  ntuple = %6.1f\n\n", recalc_ht, HT ) ;
      if ( recalc_njets != NJets ) printf("\n\n *** Recalc of NJets failed?  recalc = %2d,  ntuple = %2d\n\n", recalc_njets, NJets ) ;
      if ( fabs( recalc_dphi1 - DeltaPhi1 ) > 0.01 ) printf("\n\n *** Recalc of DeltaPhi1 failed?  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi1, DeltaPhi1 ) ;
      if ( fabs( recalc_dphi2 - DeltaPhi2 ) > 0.01 ) printf("\n\n *** Recalc of DeltaPhi2 failed?  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi2, DeltaPhi2 ) ;
      if ( fabs( recalc_dphi3 - DeltaPhi3 ) > 0.01 ) printf("\n\n *** Recalc of DeltaPhi3 failed?  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi3, DeltaPhi3 ) ;
      if ( fabs( recalc_dphi4 - DeltaPhi4 ) > 0.01 ) printf("\n\n *** Recalc of DeltaPhi4 failed?  recalc = %6.3f,  ntuple = %6.3f\n", recalc_dphi4, DeltaPhi4 ) ;






      if ( verb ) printf("\n Gen jets:\n") ;
      for ( unsigned int ji=0; ji<GenJets->size(); ji++ ) {
         if ( GenJets->at(ji).Pt() < 20 ) break ;
         if ( verb ) printf("  Jet %3d :  Pt = %5.0f,  eta = %6.3f,  phi = %6.3f\n",
            ji, GenJets->at(ji).Pt(), GenJets->at(ji).Eta(), GenJets->at(ji).Phi() ) ;
      } // ji

      if ( verb ) printf("\n Gen particles:\n") ;
      for ( unsigned int pi=0; pi<GenParticles->size(); pi++ ) {
         if ( GenParticles->at(pi).Pt() < 20 ) continue ;
         if ( verb ) printf("  Particle %3d :  PDG=%6d, parentPDG=%6d : Pt = %5.0f,  eta = %6.3f,  phi = %6.3f\n",
            pi, GenParticles_PdgId->at(pi), GenParticles_ParentId->at(pi),
            GenParticles->at(pi).Pt(), GenParticles->at(pi).Eta(), GenParticles->at(pi).Phi() ) ;
      } // pi




      if ( verb ) printf("\n Reco - Gen jets matching:\n") ;

      int worst_reco_ind(-1) ;
      double worst_dpt_rec_gen(0.) ;
      double worst_ptrec_over_ptgen(0.) ;

      int worst_reco_withnu_ind(-1) ;
      double worst_dpt_rec_withnu_gen(0.) ;
      double worst_ptrec_over_ptgen_withnu(0.) ;

      for ( unsigned int rji=0; rji<Jets->size(); rji++ ) {

         if ( Jets -> at(rji).Pt() < 25 ) continue ;

         int match_gji ;
         double match_dr = dr_match_rec_to_genjet( rji, match_gji ) ;
         float genpt(0.), geneta(0.), genphi(0.) ;
         float dpt(0.) ;
         float dpt_withnu(0.) ;
         float ptratio(0.) ;
         float ptratio_withnu(0.) ;
         if ( match_gji >= 0 && match_dr < 0.50 ) {
            genpt  = GenJets -> at( match_gji ).Pt() ;
            geneta = GenJets -> at( match_gji ).Eta() ;
            genphi = GenJets -> at( match_gji ).Phi() ;
            dpt = fabs( Jets -> at(rji).Pt() - GenJets -> at( match_gji ).Pt() ) ;
            if ( dpt > worst_dpt_rec_gen ) {
               worst_dpt_rec_gen = dpt ;
               worst_ptrec_over_ptgen = ( Jets -> at(rji).Pt() )/( GenJets -> at( match_gji ).Pt() ) ;
               worst_reco_ind = rji ;
            }
         }

         dpt_withnu = dpt ;

         ptratio = 0. ;
         if ( genpt > 0 ) { ptratio = (Jets->at(rji).Pt()) / genpt ; }

         ptratio_withnu = ptratio ;

         if ( verb ) printf("  Reco Jet %3d, gen Jet %3d, DR=%5.3f :  pflav = %5d, hflav=%5d :  Pt = %5.0f  (%5.0f),  eta = %6.3f  (%6.3f),  phi = %6.3f  (%6.3f),  Pt(rec)/Pt(gen) = %6.3f \n",
            rji, match_gji, match_dr,
            Jets_partonFlavor->at(rji), Jets_hadronFlavor->at(rji),
            Jets->at(rji).Pt(), genpt, Jets->at(rji).Eta(), geneta, Jets->at(rji).Phi(), genphi, ptratio ) ;
         int n_neutrinos_added(0) ;
         if ( match_gji >= 0 && match_dr < 0.5 ) {
            TLorentzVector genjet_nuadded_p4 = GenJets -> at(match_gji) ;
            for ( unsigned int gpi=0; gpi<GenParticles->size(); gpi++ ) {
               if ( !(abs(GenParticles_PdgId->at(gpi))==12 || abs(GenParticles_PdgId->at(gpi))==14 || abs(GenParticles_PdgId->at(gpi))==16) ) continue ;
               float dphi = calc_dphi( GenParticles -> at(gpi).Phi(),  GenJets -> at(match_gji).Phi() ) ;
               float deta = GenParticles -> at(gpi).Eta() - GenJets -> at(match_gji).Eta() ;
               float dr = sqrt( dphi*dphi + deta*deta ) ;
               if ( dr < 0.20 ) {
                     n_neutrinos_added ++ ;
                     genjet_nuadded_p4 += GenParticles->at(gpi) ;
                     if ( verb ) printf("    Matched nu:  gpi=%3d,  DR=%5.3f,   Pt = %5.0f,  eta = %6.3f,  phi = %6.3f\n", 
                          gpi, dr,
                          GenParticles->at(gpi).Pt(), GenParticles->at(gpi).Eta(), GenParticles->at(gpi).Phi()
                          ) ;
               }
            } // gpi
            if ( n_neutrinos_added > 0 ) {
               double ncgj_dr = genjet_nuadded_p4.DeltaR( Jets -> at(rji) ) ;
               ptratio_withnu = 0. ;
               dpt_withnu = fabs( Jets -> at(rji).Pt() - genjet_nuadded_p4.Pt() ) ;
               if ( genjet_nuadded_p4.Pt() > 0 ) { ptratio_withnu = Jets -> at(rji).Pt() / genjet_nuadded_p4.Pt() ; }
               if ( verb ) printf("                              Neutrino corrected gen Jet, DR=%5.3f :   Pt = %5.0f  (%5.0f),  eta = %6.3f  (%6.3f),  phi = %6.3f  (%6.3f),  Pt(rec)/Pt(gen) = %6.3f \n",
                  ncgj_dr, Jets->at(rji).Pt(), genjet_nuadded_p4.Pt(),
                           Jets->at(rji).Eta(), genjet_nuadded_p4.Eta(),
                           Jets->at(rji).Phi(), genjet_nuadded_p4.Phi(), ptratio_withnu ) ;
            }
         } // have a GenJet match?
         if ( dpt_withnu > worst_dpt_rec_withnu_gen ) {
            worst_dpt_rec_withnu_gen = dpt_withnu ;
            worst_ptrec_over_ptgen_withnu = ptratio_withnu ;
            worst_reco_withnu_ind = rji ;
         }
      } // rji


      bool badjet_in_dphi(false) ;
      if ( worst_reco_ind >= 0 ) {
         if ( worst_reco_ind == rjet_ind_dphi1 ) badjet_in_dphi = true ;
         if ( worst_reco_ind == rjet_ind_dphi2 ) badjet_in_dphi = true ;
         if ( worst_reco_ind == rjet_ind_dphi3 ) badjet_in_dphi = true ;
         if ( worst_reco_ind == rjet_ind_dphi4 ) badjet_in_dphi = true ;
      } else {
         badjet_in_dphi = true ;
      }


      if ( verb ) printf("\n") ;
      if ( verb ) printf("  Worst reco         jet ind is %2d : dpt = %7.1f , %s\n", worst_reco_ind, worst_dpt_rec_gen, badjet_in_dphi?"in dphi":"NOT in dphi" ) ;
      if ( verb ) printf("  Worst reco with nu jet ind is %2d : dpt = %7.1f\n\n", worst_reco_withnu_ind, worst_dpt_rec_withnu_gen ) ;





      h_dpt -> Fill( worst_dpt_rec_gen, Weight ) ;
      h_dpt_withnu -> Fill( worst_dpt_rec_withnu_gen, Weight ) ;

      h_worst_ptrec_over_ptgen_all -> Fill( worst_ptrec_over_ptgen, Weight ) ;
      if ( GenMET < 60 ) {
         h_worst_ptrec_over_ptgen_genmet_lt_60 -> Fill( worst_ptrec_over_ptgen, Weight ) ;
      } else {
         h_worst_ptrec_over_ptgen_genmet_gt_60 -> Fill( worst_ptrec_over_ptgen, Weight ) ;
      }




      bool recalc_in_ldp(true) ;
      if ( recalc_dphi1 > 0.5 && recalc_dphi2 > 0.5 && recalc_dphi3 > 0.3 && recalc_dphi4 > 0.3 ) recalc_in_ldp = false ;




      int bi_nj, bi_nb, bi_ht, bi_mht ;
      int bi_htmht ;
      int bi_global, bi_nbsum_global ;

      set_bi( recalc_njets, BTags, recalc_ht, recalc_mht,
              bi_nj, bi_nb, bi_ht, bi_mht,    bi_htmht,   bi_global, bi_nbsum_global ) ;

      if (verb) printf("  Recalculated bin indices:  Nj %d,  Nb %d,  MHT %d,  HT %d,   HTMHT %2d,   %s : global %3d,   global Nbsum %3d\n",
           bi_nj, bi_nb, bi_mht, bi_ht, bi_htmht,
           recalc_in_ldp ? "LDP":"HDP",
           bi_global, bi_nbsum_global ) ;

      bool is_hth(false) ;
      if ( bi_ht == 3 || (bi_mht>3 && bi_ht == 2) ) is_hth = true ;

      bool is_htm(false) ;
      if ( (bi_mht<=3 && bi_ht == 2) || (bi_mht>3 && bi_ht == 1) ) is_htm = true ;

      bool is_htl(false) ;
      if ( (bi_mht<=3 && bi_ht == 1) ) is_htl = true ;


      if ( bi_global < 1 ) continue ; //--- only keep events in search or QCD control.


    //------

      h_mdp_all -> Fill( min_delta_phi, hw ) ;
      if ( bi_mht == 1 ) h_mdp_mhtc -> Fill( min_delta_phi, hw ) ;
      if ( bi_mht == 2 ) h_mdp_mht1 -> Fill( min_delta_phi, hw ) ;
      if ( bi_mht == 3 ) h_mdp_mht2 -> Fill( min_delta_phi, hw ) ;
      if ( bi_mht == 4 ) h_mdp_mht3 -> Fill( min_delta_phi, hw ) ;
      if ( bi_mht == 5 ) h_mdp_mht4 -> Fill( min_delta_phi, hw ) ;

      if ( badjet_in_dphi ) {
         h_mdp_all_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 1 ) h_mdp_mhtc_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 2 ) h_mdp_mht1_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 3 ) h_mdp_mht2_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 4 ) h_mdp_mht3_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 5 ) h_mdp_mht4_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
      } else {
         h_mdp_all_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 1 ) h_mdp_mhtc_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 2 ) h_mdp_mht1_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 3 ) h_mdp_mht2_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 4 ) h_mdp_mht3_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 5 ) h_mdp_mht4_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
      }



      h_dphiregion_all -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      if ( bi_mht == 1 ) h_dphiregion_mhtc -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      if ( bi_mht == 2 ) h_dphiregion_mht1 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      if ( bi_mht == 3 ) h_dphiregion_mht2 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      if ( bi_mht == 4 ) h_dphiregion_mht3 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      if ( bi_mht == 5 ) h_dphiregion_mht4 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;

      if ( badjet_in_dphi ) {
          h_dphiregion_all_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 1 ) h_dphiregion_mhtc_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 2 ) h_dphiregion_mht1_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 3 ) h_dphiregion_mht2_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 4 ) h_dphiregion_mht3_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 5 ) h_dphiregion_mht4_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      } else {
          h_dphiregion_all_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 1 ) h_dphiregion_mhtc_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 2 ) h_dphiregion_mht1_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 3 ) h_dphiregion_mht2_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 4 ) h_dphiregion_mht3_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
          if ( bi_mht == 5 ) h_dphiregion_mht4_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
      }



    //------

      if ( is_hth ) {

         h_mdp_hth_all -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 1 ) h_mdp_hth_mhtc -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 2 ) h_mdp_hth_mht1 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 3 ) h_mdp_hth_mht2 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 4 ) h_mdp_hth_mht3 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 5 ) h_mdp_hth_mht4 -> Fill( min_delta_phi, hw ) ;

         if ( badjet_in_dphi ) {
            h_mdp_hth_all_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 1 ) h_mdp_hth_mhtc_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 2 ) h_mdp_hth_mht1_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 3 ) h_mdp_hth_mht2_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 4 ) h_mdp_hth_mht3_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 5 ) h_mdp_hth_mht4_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         } else {
            h_mdp_hth_all_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 1 ) h_mdp_hth_mhtc_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 2 ) h_mdp_hth_mht1_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 3 ) h_mdp_hth_mht2_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 4 ) h_mdp_hth_mht3_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 5 ) h_mdp_hth_mht4_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         }



         h_dphiregion_hth_all -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 1 ) h_dphiregion_hth_mhtc -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 2 ) h_dphiregion_hth_mht1 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 3 ) h_dphiregion_hth_mht2 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 4 ) h_dphiregion_hth_mht3 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 5 ) h_dphiregion_hth_mht4 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;

         if ( badjet_in_dphi ) {
             h_dphiregion_hth_all_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 1 ) h_dphiregion_hth_mhtc_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 2 ) h_dphiregion_hth_mht1_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 3 ) h_dphiregion_hth_mht2_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 4 ) h_dphiregion_hth_mht3_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 5 ) h_dphiregion_hth_mht4_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         } else {
             h_dphiregion_hth_all_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 1 ) h_dphiregion_hth_mhtc_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 2 ) h_dphiregion_hth_mht1_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 3 ) h_dphiregion_hth_mht2_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 4 ) h_dphiregion_hth_mht3_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 5 ) h_dphiregion_hth_mht4_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         }

      } // hth?


    //------

      if ( is_htm ) {

         h_mdp_htm_all -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 1 ) h_mdp_htm_mhtc -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 2 ) h_mdp_htm_mht1 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 3 ) h_mdp_htm_mht2 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 4 ) h_mdp_htm_mht3 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 5 ) h_mdp_htm_mht4 -> Fill( min_delta_phi, hw ) ;

         if ( badjet_in_dphi ) {
            h_mdp_htm_all_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 1 ) h_mdp_htm_mhtc_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 2 ) h_mdp_htm_mht1_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 3 ) h_mdp_htm_mht2_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 4 ) h_mdp_htm_mht3_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 5 ) h_mdp_htm_mht4_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         } else {
            h_mdp_htm_all_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 1 ) h_mdp_htm_mhtc_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 2 ) h_mdp_htm_mht1_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 3 ) h_mdp_htm_mht2_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 4 ) h_mdp_htm_mht3_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 5 ) h_mdp_htm_mht4_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         }



         h_dphiregion_htm_all -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 1 ) h_dphiregion_htm_mhtc -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 2 ) h_dphiregion_htm_mht1 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 3 ) h_dphiregion_htm_mht2 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 4 ) h_dphiregion_htm_mht3 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 5 ) h_dphiregion_htm_mht4 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;

         if ( badjet_in_dphi ) {
             h_dphiregion_htm_all_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 1 ) h_dphiregion_htm_mhtc_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 2 ) h_dphiregion_htm_mht1_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 3 ) h_dphiregion_htm_mht2_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 4 ) h_dphiregion_htm_mht3_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 5 ) h_dphiregion_htm_mht4_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         } else {
             h_dphiregion_htm_all_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 1 ) h_dphiregion_htm_mhtc_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 2 ) h_dphiregion_htm_mht1_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 3 ) h_dphiregion_htm_mht2_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 4 ) h_dphiregion_htm_mht3_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 5 ) h_dphiregion_htm_mht4_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         }

      } // htm?

    //------

      if ( is_htl ) {

         h_mdp_htl_all -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 1 ) h_mdp_htl_mhtc -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 2 ) h_mdp_htl_mht1 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 3 ) h_mdp_htl_mht2 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 4 ) h_mdp_htl_mht3 -> Fill( min_delta_phi, hw ) ;
         if ( bi_mht == 5 ) h_mdp_htl_mht4 -> Fill( min_delta_phi, hw ) ;

         if ( badjet_in_dphi ) {
            h_mdp_htl_all_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 1 ) h_mdp_htl_mhtc_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 2 ) h_mdp_htl_mht1_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 3 ) h_mdp_htl_mht2_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 4 ) h_mdp_htl_mht3_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 5 ) h_mdp_htl_mht4_badj_in_dphi -> Fill( min_delta_phi, hw ) ;
         } else {
            h_mdp_htl_all_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 1 ) h_mdp_htl_mhtc_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 2 ) h_mdp_htl_mht1_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 3 ) h_mdp_htl_mht2_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 4 ) h_mdp_htl_mht3_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
            if ( bi_mht == 5 ) h_mdp_htl_mht4_badj_not_in_dphi -> Fill( min_delta_phi, hw ) ;
         }



         h_dphiregion_htl_all -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 1 ) h_dphiregion_htl_mhtc -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 2 ) h_dphiregion_htl_mht1 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 3 ) h_dphiregion_htl_mht2 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 4 ) h_dphiregion_htl_mht3 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         if ( bi_mht == 5 ) h_dphiregion_htl_mht4 -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;

         if ( badjet_in_dphi ) {
             h_dphiregion_htl_all_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 1 ) h_dphiregion_htl_mhtc_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 2 ) h_dphiregion_htl_mht1_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 3 ) h_dphiregion_htl_mht2_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 4 ) h_dphiregion_htl_mht3_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 5 ) h_dphiregion_htl_mht4_badj_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         } else {
             h_dphiregion_htl_all_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 1 ) h_dphiregion_htl_mhtc_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 2 ) h_dphiregion_htl_mht1_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 3 ) h_dphiregion_htl_mht2_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 4 ) h_dphiregion_htl_mht3_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
             if ( bi_mht == 5 ) h_dphiregion_htl_mht4_badj_not_in_dphi -> Fill( recalc_in_ldp ? 0 : 1 , hw ) ;
         }

      } // htl?

     //-----



      h_mht_vs_ht_nominal -> Fill( recalc_ht, recalc_mht, hw ) ;

      if ( !recalc_in_ldp ) {
         h_hdp -> Fill( bi_global, hw ) ;
         h_nbsum_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 0 ) h_nb0_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 1 ) h_nb1_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 2 ) h_nb2_hdp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags >= 3 ) h_nb3_hdp -> Fill( bi_nbsum_global, hw ) ;
      } else {
         h_ldp -> Fill( bi_global, hw ) ;
         h_nbsum_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 0 ) h_nb0_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 1 ) h_nb1_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags == 2 ) h_nb2_ldp -> Fill( bi_nbsum_global, hw ) ;
         if ( BTags >= 3 ) h_nb3_ldp -> Fill( bi_nbsum_global, hw ) ;
      }








      n_dump++ ;
      if ( max_dump > 0 && n_dump > max_dump ) break ;

   } // jentry

   printf("\n\n\n Done.\n\n\n") ;
   saveHist( "outputfiles/syst-2015-v2.root", "h*" ) ;

} // Loop

//==================================================================================================

double calc_dphi( double phi1, double phi2 ) {

   double rv = phi1 - phi2 ;
   if ( rv >  3.1415926 ) rv -= 2*3.1415926 ;
   if ( rv < -3.1415926 ) rv += 2*3.1415926 ;
   return rv ;

} // calc_dphi

//==================================================================================================

double syst_2015_v2::dr_match_rec_to_genjet( int rji, int& match_gji ) {

   match_gji = -1 ;

   if ( rji < 0 ) return 99. ;

   int best_gji(-1) ;
   float min_dr(99.) ;
   for ( unsigned int gji=0; gji<GenJets->size(); gji++ ) {
      if ( GenJets -> at(gji).Pt() < 15 ) continue ;
      float dphi = calc_dphi( GenJets -> at(gji).Phi(),  Jets -> at(rji).Phi() ) ;
      float deta = GenJets -> at(gji).Eta() - Jets -> at(rji).Eta() ;
      float dr = sqrt( dphi*dphi + deta*deta ) ;
      if ( dr < min_dr ) {
         best_gji = gji ;
         min_dr = dr ;
      }
   } // gji

   match_gji = best_gji ;
   return min_dr ;

} // dr_match_rec_to_genjet

//==================================================================================================

   void syst_2015_v2::setup_bins() {

      int bi ;

      bi = 0 ;
      bin_edges_nj[bi] = 2.5 ; bi++ ;
      bin_edges_nj[bi] = 4.5 ; bi++ ;
      bin_edges_nj[bi] = 6.5 ; bi++ ;
      bin_edges_nj[bi] = 8.5 ; bi++ ;
      bin_edges_nj[bi] = 99.5 ;
      nb_nj = bi ;

      bi = 0 ;
      bin_edges_nb[bi] = -0.5 ; bi++ ;
      bin_edges_nb[bi] =  0.5 ; bi++ ;
      bin_edges_nb[bi] =  1.5 ; bi++ ;
      bin_edges_nb[bi] =  2.5 ; bi++ ;
      bin_edges_nb[bi] =  99.5 ;
      nb_nb = bi ;

      bi = 0 ;
      bin_edges_mht[bi] =   250. ; bi++ ;
      bin_edges_mht[bi] =   300. ; bi++ ;
      bin_edges_mht[bi] =   350. ; bi++ ;
      bin_edges_mht[bi] =   500. ; bi++ ;
      bin_edges_mht[bi] =   750. ; bi++ ;
      bin_edges_mht[bi] = 20000. ;
      nb_mht = bi ;

      int mbi(1) ;
      nb_htmht = 0 ;
      //--- MHT bin 0
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   300. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 1
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   300. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 2
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   350. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 3
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 4
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   750. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1500. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;



      nb_global = nb_htmht * nb_nb * nb_nj ;

      printf("\n\n") ;
      printf("  nj bins : %d\n", nb_nj ) ; for ( int i=1; i<=nb_nj; i++ ) { printf("    nj bin %2d : [%.1f,%.1f]\n", i, bin_edges_nj[i-1], bin_edges_nj[i] ) ; }
      printf("  nb bins : %d\n", nb_nb ) ; for ( int i=1; i<=nb_nb; i++ ) { printf("    nb bin %2d : [%.1f,%.1f]\n", i, bin_edges_nb[i-1], bin_edges_nb[i] ) ; }
      printf("  mht bins : %d\n", nb_mht ) ; for ( int i=1; i<=nb_mht; i++ ) { printf("    mht bin %2d : [%.1f,%.1f]\n", i, bin_edges_mht[i-1], bin_edges_mht[i] ) ; }
      for ( int mbi=1; mbi<=nb_mht; mbi++ ) {
         printf("  ht bins, mht%d : %d\n", mbi, nb_ht[mbi] ) ; for ( int i=1; i<=nb_ht[mbi]; i++ ) { printf("    ht bin %2d : [%.1f,%.1f]\n", i, bin_edges_ht[mbi][i-1], bin_edges_ht[mbi][i] ) ; }
      }
      printf("\n\n") ;

      TH1F* h_nj_bin_index = new TH1F( "h_nj_bin_index", "Njet bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_nb_bin_index = new TH1F( "h_nb_bin_index", "Nb bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_htmht_bin_index = new TH1F( "h_htmht_bin_index", "HTMHT bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_ht_bin_index = new TH1F( "h_ht_bin_index", "HT bin index", nb_global, 0.5, nb_global+0.5 ) ;
      TH1F* h_mht_bin_index = new TH1F( "h_mht_bin_index", "MHT bin index", nb_global, 0.5, nb_global+0.5 ) ;

      printf("\n\n") ;
      for ( int bi=1; bi<=nb_global; bi++ ) {
         int tbi_nj; int tbi_nb; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin( bi, tbi_nj, tbi_nb, tbi_htmht, tbi_ht, tbi_mht ) ;
         h_nj_bin_index -> SetBinContent( bi, tbi_nj ) ;
         h_nb_bin_index -> SetBinContent( bi, tbi_nb ) ;
         h_htmht_bin_index -> SetBinContent( bi, tbi_htmht ) ;
         h_ht_bin_index -> SetBinContent( bi, tbi_ht ) ;
         h_mht_bin_index -> SetBinContent( bi, tbi_mht ) ;
      } // bi


      TH1F* h_nbsum_nj_bin_index = new TH1F( "h_nbsum_nj_bin_index", "Njet bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;
      TH1F* h_nbsum_htmht_bin_index = new TH1F( "h_nbsum_htmht_bin_index", "HTMHT bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;
      TH1F* h_nbsum_ht_bin_index = new TH1F( "h_nbsum_ht_bin_index", "HT bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;
      TH1F* h_nbsum_mht_bin_index = new TH1F( "h_nbsum_mht_bin_index", "MHT bin index", nb_global/nb_nb, 0.5, nb_global/nb_nb+0.5 ) ;

      printf("\n\n") ;
      for ( int bi=1; bi<=nb_global/nb_nb; bi++ ) {
         int tbi_nj; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin_nbsum( bi, tbi_nj, tbi_htmht, tbi_ht, tbi_mht ) ;
         h_nbsum_nj_bin_index -> SetBinContent( bi, tbi_nj ) ;
         h_nbsum_htmht_bin_index -> SetBinContent( bi, tbi_htmht ) ;
         h_nbsum_ht_bin_index -> SetBinContent( bi, tbi_ht ) ;
         h_nbsum_mht_bin_index -> SetBinContent( bi, tbi_mht ) ;
      } // bi

      printf("\n\n") ;

   } // setup_bins

//=====================================================================================================

   void syst_2015_v2::set_bi( int nj, int nb, double ht, double mht,
                             int& bi_nj, int& bi_nb, int& bi_ht, int& bi_mht,
                             int& bi_htmht,
                             int& bi_global, int& bi_nbsum_global
                           ) {

      bi_ht = 0 ; bi_mht = 0 ; bi_nj = 0 ; bi_nb = 0 ;
      for ( int bi = 1; bi <= nb_nj ; bi++ ) { if ( nj  >= bin_edges_nj[bi-1]  && nj < bin_edges_nj[bi]  ) { bi_nj  = bi ; break ; } }
      for ( int bi = 1; bi <= nb_nb ; bi++ ) { if ( nb  >= bin_edges_nb[bi-1]  && nb < bin_edges_nb[bi]  ) { bi_nb  = bi ; break ; } }
      for ( int bi = 1; bi <= nb_mht; bi++ ) { if ( mht >= bin_edges_mht[bi-1] && mht < bin_edges_mht[bi] ) { bi_mht = bi ; break ; } }
      if ( bi_mht>=1 && bi_mht<=nb_mht ) {
         for ( int bi = 1; bi <= nb_ht[bi_mht] ; bi++ ) { if ( ht  >= bin_edges_ht[bi_mht][bi-1]  && ht < bin_edges_ht[bi_mht][bi]  ) { bi_ht  = bi ; break ; } }
      }

      bi_htmht = 0 ;
      if ( bi_ht >=1 && bi_mht >=1 ) {
         bi_htmht = bi_ht ;
         if ( bi_mht >= 2 ) bi_htmht += 3 ;
         if ( bi_mht >= 3 ) bi_htmht += 3 ;
         if ( bi_mht >= 4 ) bi_htmht += 3 ;
         if ( bi_mht >= 5 ) bi_htmht += 2 ;
      }

      bi_global = 0 ;
      bi_nbsum_global = 0 ;
      if ( bi_nj >= 1 && bi_nb >= 1 && bi_htmht >= 1 ) {
         bi_global = (bi_nj-1) * nb_nb * nb_htmht  +  (bi_nb-1) * nb_htmht   + bi_htmht ;
         bi_nbsum_global = (bi_nj-1) * nb_htmht  +  bi_htmht ;
      }

   } // set_bi

//=====================================================================================================

   void syst_2015_v2::translate_global_bin( int gbi, int& tbi_nj, int& tbi_nb, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) {

      if ( gbi <= 0 ) { printf("\n\n *** translate_global_bin : illegal global bin index %d\n\n", gbi ) ; gSystem -> Exit(-1) ; }

      tbi_nj = (gbi-1) / ( nb_nb*nb_htmht) + 1 ;

      tbi_nb =  ( ((gbi-1)/nb_htmht) + 1 )  % nb_nb  ;
      if ( tbi_nb == 0 ) tbi_nb = nb_nb ;

      tbi_htmht = gbi % nb_htmht ;
      if ( tbi_htmht == 0 ) tbi_htmht = nb_htmht ;

      tbi_ht = ht_bi_from_htmht( tbi_htmht ) ;
      tbi_mht = mht_bi_from_htmht( tbi_htmht ) ;

      printf("  global bin %3d :  bi_nj = %2d ,  bi_nb = %2d ,  bi_htmht = %2d ,  bi_ht = %2d ,  bi_mht = %d\n",
          gbi, tbi_nj, tbi_nb, tbi_htmht, tbi_ht, tbi_mht ) ;

   } // translate_global_bin

//=====================================================================================================

   void syst_2015_v2::translate_global_bin_nbsum( int gbi_nbsum, int& tbi_nj, int& tbi_htmht, int& tbi_ht, int& tbi_mht ) {

      if ( gbi_nbsum <= 0 ) { printf("\n\n *** translate_global_bin_nbsum : illegal global bin index %d\n\n", gbi_nbsum ) ; gSystem -> Exit(-1) ; }

      tbi_nj = (gbi_nbsum-1) / (nb_htmht) + 1 ;

      tbi_htmht = gbi_nbsum % nb_htmht ;
      if ( tbi_htmht == 0 ) tbi_htmht = nb_htmht ;

      tbi_ht = ht_bi_from_htmht( tbi_htmht ) ;
      tbi_mht = mht_bi_from_htmht( tbi_htmht ) ;

      printf("  qcd nbsum global bin %3d :  bi_nj = %2d ,  bi_htmht = %2d ,  bi_ht = %2d ,  bi_mht = %d\n",
          gbi_nbsum, tbi_nj, tbi_htmht, tbi_ht, tbi_mht ) ;

   } // translate_global_bin_nbsum

//=====================================================================================================

   int  syst_2015_v2::ht_bi_from_htmht( int abi_htmht ) {

      if ( abi_htmht <= 0 ) { printf("\n\n *** ht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      if ( abi_htmht > nb_htmht ) { printf("\n\n *** ht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      int rv(0) ;
      int sum1(0) ;
      int sum2(0) ;
      for ( int mbi=1; mbi<=nb_mht; mbi++ ) {
         sum2 += nb_ht[mbi] ;
         if ( abi_htmht > sum1 && abi_htmht <= sum2 ) {
            rv = abi_htmht - sum1 ;
            return rv ;
         }
         sum1 = sum2 ;
      }
      return rv ;

   } // ht_bi_from_htmht

   //===================================

   int  syst_2015_v2::mht_bi_from_htmht( int abi_htmht ) {

      if ( abi_htmht <= 0 ) { printf("\n\n *** mht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      if ( abi_htmht > nb_htmht ) { printf("\n\n *** mht_bi_from_htmht: Illegal htmht bin index: %d\n\n", abi_htmht ) ; gSystem -> Exit(-1) ; }
      int rv(0) ;
      int sum1(0) ;
      int sum2(0) ;
      for ( int mbi=1; mbi<=nb_mht; mbi++ ) {
         sum2 += nb_ht[mbi] ;
         if ( abi_htmht > sum1 && abi_htmht <= sum2 ) {
            rv = mbi ;
            return rv ;
         }
         sum1 = sum2 ;
      }
      return rv ;

   } // ht_bi_from_htmht

   //===================================

   void syst_2015_v2::set_bin_labels( TH1F* hp ) {

      if ( hp == 0x0 ) return ;

      for ( int bi=1; bi<=nb_global/nb_nb; bi++ ) {
         int tbi_nj; int tbi_htmht; int tbi_ht; int tbi_mht ;
         translate_global_bin_nbsum( bi, tbi_nj, tbi_htmht, tbi_ht, tbi_mht ) ;
         char binlabel[100] ;
         //---------
         // sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, tbi_mht, tbi_ht, tbi_htmht, bi ) ;
         //---------
         int print_bi_htmht = tbi_htmht - nb_ht[1] ;
         int print_bi_mht = tbi_mht - 1 ;
         int print_bi = bi - nb_ht[1]*tbi_nj ;
         int control_bi = (tbi_nj-1)*nb_ht[1] + tbi_ht ;
         if ( tbi_mht == 1 ) {
            sprintf( binlabel, "Nj%d-MHTC-HT%d (C%1d)   C%2d", tbi_nj, tbi_ht, tbi_ht, control_bi ) ;
         } else {
            sprintf( binlabel, "Nj%d-MHT%d-HT%d (%2d)   %2d", tbi_nj, print_bi_mht, tbi_ht, print_bi_htmht, print_bi ) ;
         }
         //---------
         hp -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;
      } // bi

      hp -> GetXaxis() -> LabelsOption( "v" ) ;

   } // set_bin_labels

   //===================================






#include "histio.c"


   void data_turnon1() {

      TChain ch_ht("tree") ;
      TChain ch_met("tree") ;

     //-----------
      ///////ch_ht.Add("fnal-prod-v7-skims-slimmed/tree_LDP/tree_JetHT_2016B-2.6ifb-slimskim.root") ;
      ///////ch_met.Add("fnal-prod-v7-skims-slimmed/tree_LDP/tree_MET_2016B-2.6ifb-slimskim.root") ;
     //-----------
      ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016B-slimskim.root") ;
      ch_ht.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_JetHT_2016C-slimskim.root") ;
      ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016B-slimskim.root") ;
      ch_met.Add("fnal-prod-v9-skims-slimmed/tree_LDP/tree_MET_2016C-slimskim.root") ;
     //-----------

      //int nbin(35) ;
      //float xl(0.) ;
      //float xh(700.) ;

      int nbin(10) ;
      double xbins[11] = {0., 200., 220, 240, 260, 280, 300, 340, 400, 500, 700 } ;

      TH1F* h_mht_nj1_denom = new TH1F( "h_mht_nj1_denom", "MHT, Nj1, HT800" ,  nbin, xbins ) ;
      TH1F* h_mht_nj1_numer = new TH1F( "h_mht_nj1_numer", "MHT, Nj1, MET100",  nbin, xbins ) ;
      TH1F* h_eff_nj1 = new TH1F( "h_eff_nj1", "Eff, Nj1", nbin, xbins ) ;

      TH1F* h_mht_nj2_denom = new TH1F( "h_mht_nj2_denom", "MHT, nj2, HT800" ,  nbin, xbins ) ;
      TH1F* h_mht_nj2_numer = new TH1F( "h_mht_nj2_numer", "MHT, nj2, MET100",  nbin, xbins ) ;
      TH1F* h_eff_nj2 = new TH1F( "h_eff_nj2", "Eff, nj2", nbin, xbins ) ;

      TH1F* h_mht_nj3_denom = new TH1F( "h_mht_nj3_denom", "MHT, nj3, HT800" ,  nbin, xbins ) ;
      TH1F* h_mht_nj3_numer = new TH1F( "h_mht_nj3_numer", "MHT, nj3, MET100",  nbin, xbins ) ;
      TH1F* h_eff_nj3 = new TH1F( "h_eff_nj3", "Eff, nj3", nbin, xbins ) ;

      TH1F* h_mht_nj4_denom = new TH1F( "h_mht_nj4_denom", "MHT, nj4, HT800" ,  nbin, xbins ) ;
      TH1F* h_mht_nj4_numer = new TH1F( "h_mht_nj4_numer", "MHT, nj4, MET100",  nbin, xbins ) ;
      TH1F* h_eff_nj4 = new TH1F( "h_eff_nj4", "Eff, nj4", nbin, xbins ) ;


      ch_ht.Draw("MHT>>h_mht_nj1_denom","HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=3 && NJets<=4 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig" ) ;
      ch_met.Draw("MHT>>h_mht_nj1_numer", "HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=3 && NJets<=4 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig && (pass_pfmet100_trig||pass_pfmetnomu100_trig)" ) ;

      ch_ht.Draw("MHT>>h_mht_nj2_denom","HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=5 && NJets<=6 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig" ) ;
      ch_met.Draw("MHT>>h_mht_nj2_numer", "HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=5 && NJets<=6 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig && (pass_pfmet100_trig||pass_pfmetnomu100_trig)" ) ;

      ch_ht.Draw("MHT>>h_mht_nj3_denom","HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=7 && NJets<=8 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig" ) ;
      ch_met.Draw("MHT>>h_mht_nj3_numer", "HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=7 && NJets<=8 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig && (pass_pfmet100_trig||pass_pfmetnomu100_trig)" ) ;

      ch_ht.Draw("MHT>>h_mht_nj4_denom","HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=9 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig" ) ;
      ch_met.Draw("MHT>>h_mht_nj4_numer", "HBHEIsoNoiseFilter == 1 && HBHENoiseFilter == 1 && eeBadScFilter == 1 && EcalDeadCellTriggerPrimitiveFilter == 1 && NVtx > 0 && JetID && NJets>=9 && HT>900 && (MET>100 && CaloMET>80 && MET/CaloMET<5) && pass_ht800_trig && (pass_pfmet100_trig||pass_pfmetnomu100_trig)" ) ;


      for ( int bi=1; bi<=nbin; bi++ ) {

         float eff_val(0.) ;
         float eff_err(0.) ;
         float denom, numer ;

         denom = h_mht_nj1_denom -> GetBinContent( bi ) ;
         numer = h_mht_nj1_numer -> GetBinContent( bi ) ;
         if (denom>0) {
            eff_val = numer/denom ;
            eff_err = sqrt( eff_val*(1.-eff_val)/denom ) ;
         }
         h_eff_nj1 -> SetBinContent( bi, eff_val ) ;
         h_eff_nj1 -> SetBinError( bi, eff_err ) ;


         denom = h_mht_nj2_denom -> GetBinContent( bi ) ;
         numer = h_mht_nj2_numer -> GetBinContent( bi ) ;
         if (denom>0) {
            eff_val = numer/denom ;
            eff_err = sqrt( eff_val*(1.-eff_val)/denom ) ;
         }
         h_eff_nj2 -> SetBinContent( bi, eff_val ) ;
         h_eff_nj2 -> SetBinError( bi, eff_err ) ;


         denom = h_mht_nj3_denom -> GetBinContent( bi ) ;
         numer = h_mht_nj3_numer -> GetBinContent( bi ) ;
         if (denom>0) {
            eff_val = numer/denom ;
            eff_err = sqrt( eff_val*(1.-eff_val)/denom ) ;
         }
         h_eff_nj3 -> SetBinContent( bi, eff_val ) ;
         h_eff_nj3 -> SetBinError( bi, eff_err ) ;


         denom = h_mht_nj4_denom -> GetBinContent( bi ) ;
         numer = h_mht_nj4_numer -> GetBinContent( bi ) ;
         if (denom>0) {
            eff_val = numer/denom ;
            eff_err = sqrt( eff_val*(1.-eff_val)/denom ) ;
         }
         h_eff_nj4 -> SetBinContent( bi, eff_val ) ;
         h_eff_nj4 -> SetBinError( bi, eff_err ) ;

      } // bi

      gStyle -> SetOptStat(0) ;

      h_eff_nj1 -> SetMaximum( 1.05 ) ;
      h_eff_nj1 -> SetMinimum( 0.5 ) ;

      h_eff_nj2 -> SetMaximum( 1.05 ) ;
      h_eff_nj2 -> SetMinimum( 0.5 ) ;

      h_eff_nj3 -> SetMaximum( 1.05 ) ;
      h_eff_nj3 -> SetMinimum( 0.5 ) ;

      h_eff_nj4 -> SetMaximum( 1.05 ) ;
      h_eff_nj4 -> SetMinimum( 0.5 ) ;


      h_eff_nj1 -> SetMarkerStyle(20) ;
      h_eff_nj2 -> SetMarkerStyle(21) ;
      h_eff_nj3 -> SetMarkerStyle(22) ;
      h_eff_nj4 -> SetMarkerStyle(23) ;

      h_eff_nj2 -> SetMarkerColor(2) ;
      h_eff_nj3 -> SetMarkerColor(3) ;
      h_eff_nj4 -> SetMarkerColor(4) ;

      h_eff_nj1 -> Draw() ;
      gPad -> SetGridx(1) ;
      gPad -> SetGridy(1) ;
      h_eff_nj2 -> Draw("same") ;
      h_eff_nj3 -> Draw("same") ;
      h_eff_nj4 -> Draw("same") ;


      saveHist("outputfiles/data-turnon.root", "h*") ;



   } // data_turnon1

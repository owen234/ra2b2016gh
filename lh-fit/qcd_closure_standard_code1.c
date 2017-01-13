

//
//  Copied from Simon's version on Jan 13th, 2017
//    https://github.com/skurz/Lost_Lepton_delphiClass/blob/master/Plot_searchBin_full_extended2016.C
//
//

#include "TPad.h"
#include "TExec.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "../binning.h"
#include "../histio.c"
#include "../get_hist.h"

#include <vector>
#include <cstdio>
#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string
#include "tdrstyle.C"
#include "CMS_lumi.C"


using namespace std;

void shift_bin(TH1* input, TH1* output){

  char tempname[200];  
  char temptitle[200];  
  output->SetName(tempname);
  output->SetTitle(temptitle);
  output->SetBins(input->GetNbinsX(),input->GetBinLowEdge(1)-0.5,input->GetBinLowEdge(input->GetNbinsX()+1)-0.5);
  //input->Print("all");
  //output = new TH1F(tempname,temptitle,input->GetNbinsX(),input->GetBinLowEdge(1)-0.5,input->GetBinLowEdge(input->GetNbinsX()+1)-0.5); 
  // 0: underflow
  // 1: first bin [Use the lowedge of this bin]
  // input->GetNbinsX(): highest bin 
  // input->GetNbinsX()+1: overflow bin [use the lowedge of this bin]
  //

  for (int ibin=1;ibin<=input->GetNbinsX();ibin++){
    output->SetBinContent(ibin,input->GetBinContent(ibin));    
    output->SetBinError(ibin,input->GetBinError(ibin));    
    //std::cout << input->GetBinContent(ibin) << std::endl;
  }

}

void qcd_closure_standard_code1(string option="", int pull=0){ // string option="QCD"

  // Use option="QCD" to produce plots in QCD binning

  char tempname[200];
  // Open root file
  ///////// sprintf(tempname,"outputfiles/model-ratio-hist4.root");

  // true: do closure test (MC prediction vs MC truth)
  // false: do data driven prediction and compare to MC truth
  bool doDataVsMC = false;

  // Add systematics in quadrature to stat. uncertainty on prediction
  // Non-closure systematic not included yet!
  bool showSystematics = false;

  bool doClosurewoIsoTrackVeto = false;

  ///////////////////////////////////////////////////////////////////////////////////////////
  ////Some cosmetic work for official documents.
  //
  // Set basic style
  //
  setTDRStyle();
  gStyle->SetPalette(1) ; // for better color output

  //
  // Canvas size
  int W = 1200;
  int H = 740;
  int H_ref = 740;
  int W_ref = 800;
  float T = 0.10*H_ref;
  float B = 0.06*H_ref;
  float L = 0.16*W_ref;
  float R = 0.04*W_ref;

  //
  // Various vertical line coordinates
  float ymax_top = 300000.;
  float ymin_top = 0.015;

  float ymax2_top = 1000.;
  float ymax3_top = 200.;
  float ymax4_top = 30.;
  float ymax5_top = 5.;





  ////////////////// float ymax_bottom = 1.99;
  ////////////////// float ymin_bottom = 0.01;

  ////////////////// float ymax2_bottom = 2.15;
  ////////////////// float ymax3_bottom = 2.15;
  ////////////////// float ymax4_bottom = 2.15;


  float ymax_bottom = 4.99;
  float ymin_bottom = 0.01;

  float ymax2_bottom = 5.15;
  float ymax3_bottom = 5.15;
  float ymax4_bottom = 5.15;






  //
  // Luminosity information for scaling
  double lumi     = 36.3; // normaliza to this lumi (fb-1)
  double lumi_ref = 36.3; // normaliza to 3 (fb-1)
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  //
  // More specific style set, opening input files etc

  gStyle->SetOptStat(0);  ///to avoid the stat. on the plots
  //gStyle->SetErrorX(0);
  char xtitlename[200];
  char ytitlename[200];

  //TFile * LLFile = new TFile(tempname,"R");
  //printf("Opened %s\n",tempname);

   //+++++++++++ Owen code (

      bool do_text_binlabels(false) ;
      setup_bins();

      gDirectory -> Delete( "h*" ) ;
      loadHist( "outputfiles/model-ratio-hist4.root" ) ;

      TH1F* h_ldp_search_bins_qcdmc = get_hist( "h_ldp_search_bins_qcdmc" ) ;
      TH1F* h_hdp_search_bins_qcdmc = get_hist( "h_hdp_search_bins_qcdmc" ) ;
      TH1F* h_ratio_all = get_hist( "h_ratio_all" ) ;
      TH1F* h_ratio_qcdmc_minus_model = get_hist( "h_ratio_qcdmc_minus_model" ) ;

      int nbins = h_ldp_search_bins_qcdmc->GetNbinsX() ;

      //////TH1F* h_hdp_evts = h_hdp_search_bins_qcdmc ;
      TH1F* h_hdp_evts_orig = h_hdp_search_bins_qcdmc ;
      TH1F* h_hdp_model = new TH1F( "h_hdp_model", "HDP, QCD model", nbins, 0.5, nbins+0.5 ) ;
      TH1F* h_hdp_evts_over_model = new TH1F( "h_hdp_evts_over_model", "HDP, QCD MC events / model", nbins, 0.5, nbins+0.5 ) ;
      TH1F* h_hdp_model_over_model = new TH1F( "h_hdp_model_over_model", "HDP, QCD model / model", nbins, 0.5, nbins+0.5 ) ;

      TH1F* h_hdp_evts = new TH1F( "h_hdp_evts","", nb_global_after_exclusion, 0.5, nb_global_after_exclusion + 0.5 ) ;

      for ( int bi=1; bi<=nbins; bi++ ) {

         h_hdp_evts -> SetBinContent( bi, h_hdp_evts_orig -> GetBinContent( bi ) ) ;
         h_hdp_evts -> SetBinError( bi, h_hdp_evts_orig -> GetBinError( bi ) ) ;

         char label[100] ;
         sprintf( label, "%s", h_ldp_search_bins_qcdmc -> GetXaxis() -> GetBinLabel( bi ) ) ;

         float ldp_val = h_ldp_search_bins_qcdmc -> GetBinContent( bi ) ;
         float ldp_err = h_ldp_search_bins_qcdmc -> GetBinError( bi ) ;

         float hdp_val = h_hdp_search_bins_qcdmc -> GetBinContent( bi ) ;
         float hdp_err = h_hdp_search_bins_qcdmc -> GetBinError( bi ) ;

         float ratio_model_val = h_ratio_all -> GetBinContent( bi ) ;
         float ratio_model_err = h_ratio_all -> GetBinError( bi ) ;

         float ratio_diff_val = h_ratio_qcdmc_minus_model -> GetBinContent( bi ) ;
         float ratio_diff_err = h_ratio_qcdmc_minus_model -> GetBinError( bi ) ;




         ///float ratio_val = ratio_model_val + ratio_diff_val ;
         ///float ratio_err = sqrt( pow( ratio_model_err, 2. ) + pow( ratio_diff_err, 2. ) ) ;

         /////float ratio_val = ratio_model_val  ;
         /////float ratio_err = sqrt( pow( ratio_model_err, 2. ) + pow( ratio_diff_err, 2. ) ) ;

         float ratio_val = ratio_model_val ;
         float ratio_err = ratio_model_err ;





         float hdp_pred_val = ldp_val * ratio_val ;
         float hdp_pred_err(0.) ;
         if ( ldp_val > 0 && ratio_val > 0 ) {
            hdp_pred_err = hdp_pred_val * sqrt( pow( ldp_err/ldp_val, 2. ) + pow( ratio_err/ratio_val, 2. ) ) ;
         }
         h_hdp_model -> SetBinContent( bi, hdp_pred_val ) ;
         h_hdp_model -> SetBinError( bi, hdp_pred_err ) ;
         h_hdp_model -> GetXaxis() -> SetBinLabel( bi, label ) ;

         if ( hdp_pred_val > 0 ) {
            h_hdp_evts_over_model -> SetBinContent( bi, (hdp_val/hdp_pred_val) ) ;
            h_hdp_evts_over_model -> SetBinError( bi, (hdp_err/hdp_pred_val) ) ;
            h_hdp_model_over_model -> SetBinContent( bi, 1. ) ;
            h_hdp_model_over_model -> SetBinError( bi, (hdp_pred_err/hdp_pred_val) ) ;
         }
         if ( do_text_binlabels ) {
            h_hdp_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
            h_hdp_model_over_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
            h_hdp_evts_over_model -> GetXaxis() -> SetBinLabel( bi, label ) ;
         }

      } // bi

      h_hdp_model -> GetXaxis() -> LabelsOption( "v" ) ;

   //+++++++++++ Owen code )






  //
  // Define legend
  //
  Float_t legendX1 = .655; //.50;
  Float_t legendX2 = .955; //.70;
  Float_t legendY1 = .58; //.65;
  Float_t legendY2 = .78;

  TLegend* catLeg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);
  //catLeg1->SetTextSize(0.060);
  catLeg1->SetTextSize(0.044);
  catLeg1->SetTextFont(42);
  catLeg1->SetFillColor(0);
  catLeg1->SetLineColor(1);
  catLeg1->SetBorderSize(1);

  //
  // Define canvas
  //
  TCanvas *canvas = new TCanvas("canvas","canvas",10,10,W,H);

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  canvas->Divide(1, 2);
  
  //
  // Define pads
  //
  TPad* canvas_up = (TPad*) canvas->GetListOfPrimitives()->FindObject("canvas_1");
  TPad* canvas_dw = (TPad*) canvas->GetListOfPrimitives()->FindObject("canvas_2");

  //
  // define the size
  double up_height     = 0.8;  // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.30; // please tune so that the smaller canvas size will work in your environment
  double font_size_dw  = 0.1;  // please tune the font size parameter for bottom figure
  double dw_height     = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.04; // KH, added to put the bottom one closer to the top panel

  //
  // set pad size
  canvas_up->SetPad(0., 1 - up_height,    1., 1.00);
  canvas_dw->SetPad(0., 0.,               1., dw_height+dw_height_offset);
  //
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_up->SetTopMargin(0.12);
  canvas_up->SetLeftMargin(0.1);
  //
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_dw->SetBottomMargin(0.35);
  canvas_dw->SetTopMargin(0);
  canvas_dw->SetLeftMargin(0.1);
  
  //
  // draw top figure
  canvas_up->cd();

  TH1F * GenHist, * EstHist,* thist;
  TH1F * GenHistTemp, * EstHistTemp;
  TH1F * GenHistD, * EstHistD;
  TH1F * GenHistDTemp, * EstHistDTemp;
  TH1F * EstSystematics;
  TH1F * histTemplate;

  double HT_x_max=2500.;
  double HT_x_min=400.;
  double MHT_x_max=1000.;
  double NJet_x_max=15.;
  double NBtag_x_max=4.;
  double search_x_max=174.+0.5;
  if(option.find("QCD")!=string::npos)search_x_max=223.+0.5;
  double search_x_min=1.-0.5;

  //////////  TDirectory *dPre = 0;
  //////////  TDirectory *dExp = (TDirectory*)LLFile->Get("Expectation");

  //////////  if(doDataVsMC){
  //////////    dPre = (TDirectory*)LLFile->Get("Prediction_data");
  //////////  }else{
  //////////    dPre = (TDirectory*)LLFile->Get("Prediction_MC");
  //////////  }  

  //////////  if(doDataVsMC){    
  //////////      EstHistTemp=(TH1F*) dPre->Get("totalPred_LL")->Clone();
  //////////      EstHistDTemp=(TH1F*) dPre->Get("totalPred_LL")->Clone();
  //////////  }else{
  //////////    if(doClosurewoIsoTrackVeto){
  //////////      EstHistTemp=(TH1F*) dPre->Get("totalPred_woIsoTrack_LL_MC")->Clone();
  //////////      EstHistDTemp=(TH1F*) dPre->Get("totalPred_woIsoTrack_LL_MC")->Clone();
  //////////    }else{
  //////////      EstHistTemp=(TH1F*) dPre->Get("totalPred_LL_MC")->Clone();
  //////////      EstHistDTemp=(TH1F*) dPre->Get("totalPred_LL_MC")->Clone();
  //////////    }
  //////////  }

  EstHistTemp = h_hdp_model ;
  EstHistDTemp = h_hdp_model ;


  //////////   if(doClosurewoIsoTrackVeto){
  //////////       GenHistTemp=(TH1F*) dExp->Get("totalExp_woIsoTrack_LL")->Clone();
  //////////       GenHistDTemp=(TH1F*) dExp->Get("totalExp_woIsoTrack_LL")->Clone();;
  //////////   }else{
  //////////       GenHistTemp=(TH1F*) dExp->Get("totalExp_LL")->Clone();
  //////////       GenHistDTemp=(TH1F*) dExp->Get("totalExp_LL")->Clone();
  //////////   }

  GenHistTemp = h_hdp_evts ;
  GenHistDTemp = h_hdp_evts ;




  //////////    if(showSystematics){
  //////////      TDirectory *dSys = (TDirectory*)LLFile->Get("AdditionalContent");
  //////////      if(doDataVsMC){
  //////////        EstSystematics=(TH1F*) dSys->Get("totalPropSysUp_LL")->Clone();
  //////////      }else{
  //////////        EstSystematics=(TH1F*) dSys->Get("totalPropSysUp_LL_MC")->Clone();
  //////////      }
  //////////    }

  if(EstHistTemp->GetNbinsX() != GenHistTemp->GetNbinsX()) std::cout<<"NbinsX of Expectation and Prediction don't agree!"<<std::endl;

  EstHist = new TH1F("Exp", "Exp", EstHistTemp->GetNbinsX(), 0.5, EstHistTemp->GetNbinsX()+0.5);
  EstHistD = new TH1F("ExpD", "ExpD", EstHistDTemp->GetNbinsX(), 0.5, EstHistDTemp->GetNbinsX()+0.5);
  GenHist = new TH1F("Pred", "Pred", GenHistTemp->GetNbinsX(), 0.5, GenHistTemp->GetNbinsX()+0.5);
  GenHistD = new TH1F("PredD", "PredD", GenHistDTemp->GetNbinsX(), 0.5, GenHistDTemp->GetNbinsX()+0.5);

  for(int i = 0; i <= EstHistTemp->GetNbinsX()+1; i++){
    EstHist->SetBinContent(i, EstHistTemp->GetBinContent(i));    
    EstHistD->SetBinContent(i, EstHistDTemp->GetBinContent(i));

    if(showSystematics){
      EstHist->SetBinError(i, std::sqrt(EstHistTemp->GetBinError(i)*EstHistTemp->GetBinError(i)+EstSystematics->GetBinContent(i)*EstSystematics->GetBinContent(i)*EstHistTemp->GetBinError(i)*EstHistTemp->GetBinError(i)));
      EstHistD->SetBinError(i, std::sqrt(EstHistDTemp->GetBinError(i)*EstHistDTemp->GetBinError(i)+EstSystematics->GetBinContent(i)*EstSystematics->GetBinContent(i)*EstHistDTemp->GetBinError(i)*EstHistDTemp->GetBinError(i)));
    }else{
      EstHist->SetBinError(i, EstHistTemp->GetBinError(i));
      EstHistD->SetBinError(i, EstHistDTemp->GetBinError(i));
    }    

    GenHist->SetBinContent(i, GenHistTemp->GetBinContent(i));
    GenHist->SetBinError(i, GenHistTemp->GetBinError(i));
    GenHistD->SetBinContent(i, GenHistDTemp->GetBinContent(i));
    GenHistD->SetBinError(i, GenHistDTemp->GetBinError(i));
  }

  GenHist->SetLineColor(4);
  EstHist->SetLineColor(4);
  //GenHist->GetXaxis()->SetLabelFont(42);
  //GenHist->GetXaxis()->SetLabelOffset(0.007);
  //GenHist->GetXaxis()->SetLabelSize(0.04);
  //GenHist->GetXaxis()->SetTitleSize(0.05);
  //GenHist->GetXaxis()->SetTitleOffset(0.9);
  //GenHist->GetXaxis()->SetTitleOffset(0.5);
  //GenHist->GetXaxis()->SetTitleFont(42);
  //GenHist->GetYaxis()->SetLabelFont(42);
  //GenHist->GetYaxis()->SetLabelOffset(0.007);
  //GenHist->GetYaxis()->SetLabelSize(0.04);
  GenHist->GetYaxis()->SetLabelSize(0.045*1.15);
  GenHist->GetYaxis()->SetTitleSize(0.06*1.15);
  GenHist->GetYaxis()->SetTitleOffset(0.6);
  GenHist->GetYaxis()->SetTitleFont(42);


  //EstHist->GetXaxis()->SetLabelFont(42);
  //EstHist->GetXaxis()->SetLabelOffset(0.007);
  //EstHist->GetXaxis()->SetLabelSize(0.04);
  //EstHist->GetXaxis()->SetTitleSize(0.05);
  //EstHist->GetXaxis()->SetTitleOffset(0.9);
  //EstHist->GetXaxis()->SetTitleFont(42);
  //EstHist->GetYaxis()->SetLabelFont(42);
  //EstHist->GetYaxis()->SetLabelOffset(0.007);
  //EstHist->GetYaxis()->SetLabelSize(0.04);
  //EstHist->GetYaxis()->SetTitleSize(0.08);
  //EstHist->GetYaxis()->SetTitleOffset(2.0);
  //EstHist->GetYaxis()->SetTitleFont(42);
  sprintf(xtitlename,"Search region bin number");
  sprintf(ytitlename,"Events");
  gPad->SetLogy();
  GenHist->SetMaximum(ymax_top);
  GenHist->SetMinimum(ymin_top);
  GenHist->GetXaxis()->SetRangeUser(search_x_min,search_x_max);

  //GenHist->GetYaxis()->SetTickLength(0.015);
  //GenHist->GetXaxis()->SetTickLength(0.02);

  //gPad->SetGridx(1);
  TExec *ex1 = new TExec("ex1","gStyle->SetErrorX(0);");
  TExec *ex2 = new TExec("ex2","gStyle->SetErrorX(0.5);");

  GenHist->SetTitle("");
  GenHist->SetMarkerStyle(20);
  GenHist->SetMarkerSize(1.2);
  GenHist->SetLineColor(1);
  GenHist->GetXaxis()->SetTitle(xtitlename);
  GenHist->GetYaxis()->SetTitle(ytitlename);
  GenHist->Scale(lumi/lumi_ref);
  EstHist->Scale(lumi/lumi_ref);
  TH1F * GenHist_Normalize = static_cast<TH1F*>(GenHist->Clone("GenHist_Normalize"));
  GenHist_Normalize->SetMaximum(ymax_top);
  GenHist_Normalize->SetMinimum(ymin_top);
  ex1->Draw();
  //GenHist_Normalize->GetListOfFunctions()->Add(ex1);
  GenHist_Normalize->DrawCopy("e");

  EstHist->SetFillStyle(3144);
  EstHist->SetFillColor(kRed-10);
  EstHist->SetMarkerStyle(20);
  EstHist->SetMarkerSize(0.0001);
  TH1F * EstHist_Normalize = static_cast<TH1F*>(EstHist->Clone("EstHist_Normalize"));
  ex2->Draw();
  //EstHist_Normalize->GetListOfFunctions()->Add(ex2);
  EstHist_Normalize->DrawCopy("e2same");
  //EstHist_Normalize->DrawCopy("esame");

  TH1F *EstHist_Normalize_Clone = (TH1F*)EstHist_Normalize->Clone();
  for(int i=1; i<EstHist_Normalize_Clone->GetNbinsX(); i++) {
    /////////////EstHist_Normalize_Clone->SetBinError(i,0);
    EstHist_Normalize_Clone->SetBinError(i,0.00000001);
  }
  EstHist_Normalize_Clone->SetFillColor(kWhite);
  EstHist_Normalize_Clone->Draw("esame");

  GenHist->Print("all");
  EstHist->Print("all");
  
  //
  // Re-draw to have "expectation" on top of "prediction"
  ex1->Draw();
  GenHist_Normalize->DrawCopy("esame");
  //

  TString line = "";
  sprintf(tempname,"%8.1f",lumi);
  line+=tempname;
  line+=" fb^{-1} (13 TeV)";
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=0;
    
  writeExtraText = true;  
  if(doDataVsMC) extraText   = "        Preliminary";
  else extraText   = "        Simulation";
  //float extraTextFont = 52;  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  //lumiTextSize     = 0.5;
  //float lumiTextOffset   = 0.2;
  //cmsTextSize      = 0.65;
  //float cmsTextOffset    = 0.1;  // only used in outOfFrame version
  
  //relPosX    = 0.045;
  //relPosY    = 0.035;
  //relExtraDY = 1.2;
  
  // ratio of "CMS" and extra text size
  //float extraOverCmsTextSize  = 0.76;
    
  //TString lumi_13TeV = "20.1 fb^{-1}";
  //TString lumi_8TeV  = "19.7 fb^{-1}";
  //TString lumi_7TeV  = "5.1 fb^{-1}";
  TString lumi_sqrtS = line;

  //
  if(option.find("QCD")==string::npos ){
    
    //-----------------------------------------------------------
    // Putting lines and labels explaining search region definitions
    //-----------------------------------------------------------

    //TString CMSlabel = "";
    //cmsText = "#bf{CMS} #it{Simulation}";
    //CMSlabel += "#splitline{#bf{CMS}}{#scale[0.6]{#it{Simulation}}}";

    /*
    double x0 = gStyle->GetPadLeftMargin();
    double x1 = 1.-gStyle->GetPadRightMargin();
    double y0 = 1.005-gStyle->GetPadTopMargin();
    double y1 = 0.96;
    TPaveText *Lumitxt = new TPaveText(x0,y0,x1,y1,"NDC");
    Lumitxt->SetBorderSize(0);
    Lumitxt->SetFillColor(0);
    Lumitxt->SetTextFont(42);
    Lumitxt->SetTextAlign(31);
    Lumitxt->SetTextSize(1.2*gStyle->GetPadTopMargin());
    Lumitxt->SetMargin(0.);
    Lumitxt->AddText(line);
    //Lumitxt->Draw("same");

    x0 = gStyle->GetPadLeftMargin()+0.03;
    x1 = gStyle->GetPadLeftMargin()+0.13;
    y0 = 0.905-gStyle->GetPadTopMargin();
    y1 = 0.88;
    TPaveText *CMStxt = new TPaveText(x0,y0,x1,y1,"NDC");
    CMStxt->SetBorderSize(0);
    CMStxt->SetFillColor(0);
    CMStxt->SetTextFont(42);
    CMStxt->SetTextAlign(11);
    CMStxt->SetTextSize(1.2*gStyle->GetPadTopMargin());
    CMStxt->SetMargin(0.);
    CMStxt->AddText(CMSlabel);
    //CMStxt->Draw("same");
    */

    // Njet separation lines
    TLine *tl_njet = new TLine();
    tl_njet->SetLineStyle(2);
    tl_njet->DrawLine(31.-0.5,ymin_top,31.-0.5,ymax_top); 
    tl_njet->DrawLine(71.-0.5,ymin_top,71.-0.5,ymax_top); 
    tl_njet->DrawLine(111.-0.5,ymin_top,111.-0.5,ymax_top); 
    tl_njet->DrawLine(143.-0.5,ymin_top,143.-0.5,ymax_top); 

    // Njet labels
    TLatex * ttext_njet = new TLatex();
    ttext_njet->SetTextFont(42);
    ttext_njet->SetTextSize(0.04);
    ttext_njet->SetTextAlign(22);
    ttext_njet->DrawLatex(15.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} = 2");
    ttext_njet->DrawLatex(50.-0.5 , ymax_top/4. , "3 #leq N_{#scale[0.2]{ }jet} #leq 4");
    ttext_njet->DrawLatex(90.-0.5 , ymax_top/4. , "5 #leq N_{#scale[0.2]{ }jet} #leq 6");
    ttext_njet->DrawLatex(127.-0.5 , ymax_top/4. , "7 #leq N_{#scale[0.2]{ }jet} #leq 8");
    ttext_njet->DrawLatex(158.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 9");

    // Nb separation lines
    TLine *tl_nb = new TLine();
    tl_nb->SetLineStyle(3);
    tl_nb->DrawLine(11.-0.5,ymin_top,11.-0.5,ymax2_top); 
    tl_nb->DrawLine(21.-0.5,ymin_top,21.-0.5,ymax2_top); 
    tl_nb->DrawLine(31.-0.5,ymin_top,31.-0.5,ymax2_top);
    tl_nb->DrawLine(41.-0.5,ymin_top,41.-0.5,ymax2_top);
    tl_nb->DrawLine(51.-0.5,ymin_top,51.-0.5,ymax2_top); 
    tl_nb->DrawLine(61.-0.5,ymin_top,61.-0.5,ymax2_top); 
    tl_nb->DrawLine(71.-0.5,ymin_top,71.-0.5,ymax2_top);
    tl_nb->DrawLine(81.-0.5,ymin_top,81.-0.5,ymax3_top);
    tl_nb->DrawLine(91.-0.5,ymin_top,91.-0.5,ymax3_top); 
    tl_nb->DrawLine(101.-0.5,ymin_top,101.-0.5,ymax3_top); 
    tl_nb->DrawLine(111.-0.5,ymin_top,111.-0.5,ymax3_top); 
    tl_nb->DrawLine(119.-0.5,ymin_top,119.-0.5,ymax4_top); 
    tl_nb->DrawLine(127.-0.5,ymin_top,127.-0.5,ymax4_top); 
    tl_nb->DrawLine(135.-0.5,ymin_top,135.-0.5,ymax4_top); 
    tl_nb->DrawLine(143.-0.5,ymin_top,143.-0.5,ymax4_top);
    tl_nb->DrawLine(151.-0.5,ymin_top,151.-0.5,ymax5_top);
    tl_nb->DrawLine(159.-0.5,ymin_top,159.-0.5,ymax5_top);
    tl_nb->DrawLine(167.-0.5,ymin_top,167.-0.5,ymax5_top);

    
    // Nb labels
    TLatex * ttext_nb = new TLatex();
    ttext_nb->SetTextFont(42);
    ttext_nb->SetTextSize(0.04);
    ttext_nb->SetTextAlign(22);
    
    ttext_nb->DrawLatex(9.-0.5 , ymax_top/12. , "N_{#scale[0.2]{ }b-jet}");
    ttext_nb->DrawLatex(6.-0.5 , ymax_top/40. , "0");
    ttext_nb->DrawLatex(16.-0.5 , ymax_top/40. , "1");
    ttext_nb->DrawLatex(26.-0.5 , ymax_top/40. , "2");
    ttext_nb->DrawLatex(36.-0.5 , ymax_top/40. , "0");
    ttext_nb->DrawLatex(46.-0.5 , ymax_top/40. , "1");
    ttext_nb->DrawLatex(56.-0.5 , ymax_top/40. , "2");
    ttext_nb->DrawLatex(66.-0.5 , ymax_top/40. , "#geq 3");

    //
  } else {
    
    //-----------------------------------------------------------
    // Putting lines and labels explaining search region definitions
    //-----------------------------------------------------------

    // Njet separation lines
    TLine *tl_njet = new TLine();
    tl_njet->SetLineStyle(2);
    tl_njet->DrawLine(40.-0.5,ymin_top,40.-0.5,ymax_top); 
    tl_njet->DrawLine(92.-0.5,ymin_top,92.-0.5,ymax_top); 
    tl_njet->DrawLine(144.-0.5,ymin_top,144.-0.5,ymax_top); 
    tl_njet->DrawLine(184.-0.5,ymin_top,184.-0.5,ymax_top); 

    // Njet labels
    TLatex * ttext_njet = new TLatex();
    ttext_njet->SetTextFont(42);
    ttext_njet->SetTextSize(0.04);
    ttext_njet->SetTextAlign(22);
    ttext_njet->DrawLatex(20.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} = 2");
    ttext_njet->DrawLatex(66.-0.5 , ymax_top/4. , "3 #leq N_{#scale[0.2]{ }jet} #leq 4");
    ttext_njet->DrawLatex(118.-0.5 , ymax_top/4. , "5 #leq N_{#scale[0.2]{ }jet} #leq 6");
    ttext_njet->DrawLatex(165.-0.5 , ymax_top/4. , "7 #leq N_{#scale[0.2]{ }jet} #leq 8");
    ttext_njet->DrawLatex(202.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 9");

    // Nb separation lines
    TLine *tl_nb = new TLine();
    tl_nb->SetLineStyle(3);
    tl_nb->DrawLine(14.-0.5,ymin_top,14.-0.5,ymax2_top); 
    tl_nb->DrawLine(27.-0.5,ymin_top,27.-0.5,ymax2_top); 
    tl_nb->DrawLine(40.-0.5,ymin_top,40.-0.5,ymax2_top);
    tl_nb->DrawLine(53.-0.5,ymin_top,53.-0.5,ymax2_top);
    tl_nb->DrawLine(66.-0.5,ymin_top,66.-0.5,ymax2_top); 
    tl_nb->DrawLine(79.-0.5,ymin_top,79.-0.5,ymax2_top); 
    tl_nb->DrawLine(92.-0.5,ymin_top,92.-0.5,ymax2_top); 
    tl_nb->DrawLine(105.-0.5,ymin_top,105.-0.5,ymax3_top); 
    tl_nb->DrawLine(118.-0.5,ymin_top,118.-0.5,ymax3_top); 
    tl_nb->DrawLine(131.-0.5,ymin_top,131.-0.5,ymax3_top); 
    tl_nb->DrawLine(144.-0.5,ymin_top,144.-0.5,ymax3_top);
    tl_nb->DrawLine(154.-0.5,ymin_top,154.-0.5,ymax4_top);
    tl_nb->DrawLine(164.-0.5,ymin_top,164.-0.5,ymax4_top);
    tl_nb->DrawLine(174.-0.5,ymin_top,174.-0.5,ymax4_top);
    tl_nb->DrawLine(184.-0.5,ymin_top,184.-0.5,ymax4_top);
    tl_nb->DrawLine(194.-0.5,ymin_top,194.-0.5,ymax5_top);
    tl_nb->DrawLine(204.-0.5,ymin_top,204.-0.5,ymax5_top);
    tl_nb->DrawLine(214.-0.5,ymin_top,214.-0.5,ymax5_top);
    
    // Nb labels
    TLatex * ttext_nb = new TLatex();
    ttext_nb->SetTextFont(42);
    ttext_nb->SetTextSize(0.04);
    ttext_nb->SetTextAlign(22);
    
    ttext_nb->DrawLatex(9.-0.5 , ymax_top/12. , "N_{#scale[0.2]{ }b-jet}");
    ttext_nb->DrawLatex(8.-0.5 , ymax_top/40. , "0");
    ttext_nb->DrawLatex(20.-0.5 , ymax_top/40. , "1");
    ttext_nb->DrawLatex(33.-0.5 , ymax_top/40. , "2");
    ttext_nb->DrawLatex(46.-0.5 , ymax_top/40. , "0");
    ttext_nb->DrawLatex(59.-0.5 , ymax_top/40. , "1");
    ttext_nb->DrawLatex(72.-0.5 , ymax_top/40. , "2");
    ttext_nb->DrawLatex(85.-0.5 , ymax_top/40. , "#geq 3");

  }

  // Legend & texts
  sprintf(tempname,"QCD background");
  catLeg1->SetHeader(tempname);
  //sprintf(tempname,"#tau_{hadronic} BG expectation (MC truth)");
  sprintf(tempname,"Direct from simulation");
  catLeg1->AddEntry(GenHist,tempname,"pe");
  //sprintf(tempname,"Prediction from MC");
  if(doDataVsMC) sprintf(tempname,"Prediction from data");
  else sprintf(tempname,"Treat simulation like data");
  catLeg1->AddEntry(EstHist,tempname);
  catLeg1->Draw();

  gPad->RedrawAxis();

  //
  // Bottom ratio plot
  //
  // ----------

    //
    // Preparing ratio histograms
      TH1F * numerator   = static_cast<TH1F*>(GenHist->Clone("numerator"));
      TH1F * numerator_fullstaterr   = static_cast<TH1F*>(GenHist->Clone("numerator_fullstaterr"));
      TH1F * denominator = static_cast<TH1F*>(EstHist->Clone("denominator"));

      TH1F * GenHist_Clone = static_cast<TH1F*>(GenHist->Clone("GenHist_Clone"));
      TH1F * EstHist_Clone = static_cast<TH1F*>(EstHist->Clone("EstHist_Clone"));
      TH1F * EstHist_NoError = static_cast<TH1F*>(EstHist->Clone("EstHist_NoError"));
      TH1F * One_NoError = static_cast<TH1F*>(EstHist->Clone("EstHist_NoError"));
      for (int ibin=0; ibin<EstHist_NoError->GetNbinsX()+2; ibin++){ // scan including underflow and overflow bins
        /////////////////EstHist_NoError->SetBinError(ibin,0.);
        EstHist_NoError->SetBinError(ibin,0.0000001);
        One_NoError->SetBinContent(ibin,1.);
        /////////////////One_NoError->SetBinError(ibin,0.);
        One_NoError->SetBinError(ibin,0.0000001);
      }

      //EstHistD->Add(GenHistD,-1);
      numerator->Divide(GenHist_Clone,EstHist_NoError,1,1,"");
      denominator->Divide(EstHist_Clone,EstHist_NoError,1,1,"");

//      if(pull==1){
//      		numerator_fullstaterr->Add(GenHist_Clone,EstHist_Clone,1,-1); // Expectation-Prediction
//      }else{
	      numerator_fullstaterr->Divide(GenHist_Clone,EstHist_Clone,1,1,"");  // Expectation/Prediction
      	  numerator_fullstaterr->Add(One_NoError,-1.);                        // Expectation/Prediction-1
//      }



      // draw bottom figure
      canvas_dw->cd();
      // font size
      numerator->GetXaxis()->SetLabelSize(font_size_dw);
      numerator->GetXaxis()->SetTitleSize(font_size_dw);
      numerator->GetYaxis()->SetLabelSize(font_size_dw);
      numerator->GetYaxis()->SetTitleSize(font_size_dw);

      //
      // Horizontal Lines
      TLine *tline  = new TLine(search_x_min,1.,search_x_max,1.);
      TLine *tline0 = new TLine(search_x_min,0.,search_x_max,0.);

      //
      // Common to all bottom plots
      //
      //sprintf(ytitlename,"#frac{Estimate - #tau_{had} BG}{#tau_{had} BG} ");
      sprintf(ytitlename,"#frac{Direct}{Prediction} ");
      numerator->SetMaximum(ymax_bottom);
      numerator->SetMinimum(ymin_bottom);

      //
      // Specific to each bottom plot
      //
      // Setting style
      //numerator->SetMaximum(1.4);
      //numerator->GetXaxis()->SetLabelFont(42);
      //numerator->GetXaxis()->SetLabelOffset(0.007);
      numerator->GetXaxis()->SetLabelSize(0.18*0.045/0.06);
      numerator->GetXaxis()->SetTitleSize(0.18);
      numerator->GetXaxis()->SetTitleOffset(0.9);
      numerator->GetXaxis()->SetTitleFont(42);
      //numerator->GetYaxis()->SetLabelFont(42);
      //numerator->GetYaxis()->SetLabelOffset(0.007);
      numerator->GetYaxis()->SetLabelSize(0.18*0.045/0.06);
      numerator->GetYaxis()->SetTitleSize(0.18);
      //numerator->GetYaxis()->SetTitleOffset(0.5);
      numerator->GetYaxis()->SetTitleOffset(0.25);
      numerator->GetYaxis()->SetTitleFont(42);

      numerator->GetXaxis()->SetTitle(xtitlename);
      numerator->GetYaxis()->SetTitle(ytitlename);

      //gPad->SetGridx(1);


      if (pull==1){

        //sprintf(ytitlename,"#frac{Exp - Pre}{Stat Error} ");
        sprintf(ytitlename,"pull(Dir./Pred.)");
        numerator_fullstaterr->SetMaximum(3.9);
        numerator_fullstaterr->SetMinimum(-3.9);
        
        //
        // Specific to each bottom plot
        //
        // Setting style

        for (int ibin=0; ibin<numerator_fullstaterr->GetNbinsX()+2; ibin++){ // scan including underflow and overflow bins
          if(numerator_fullstaterr->GetBinError(ibin) < 0.00000001){
            numerator_fullstaterr->SetBinContent(ibin, 0);
          }else{
            numerator_fullstaterr->SetBinContent(ibin,numerator_fullstaterr->GetBinContent(ibin)/numerator_fullstaterr->GetBinError(ibin));
          }         
          /////////////numerator_fullstaterr->SetBinError(ibin,0.);
          numerator_fullstaterr->SetBinError(ibin,0.00000001);
        }

        numerator_fullstaterr->Print("all");
        
        numerator_fullstaterr->GetXaxis()->SetLabelSize(font_size_dw);
        numerator_fullstaterr->GetXaxis()->SetTitleSize(font_size_dw);
        numerator_fullstaterr->GetYaxis()->SetLabelSize(font_size_dw);
        numerator_fullstaterr->GetYaxis()->SetTitleSize(font_size_dw);

        numerator_fullstaterr->GetXaxis()->SetLabelSize(0.18*0.045/0.06);
        numerator_fullstaterr->GetXaxis()->SetTitleSize(0.18);
        numerator_fullstaterr->GetXaxis()->SetTitleOffset(0.9);
        numerator_fullstaterr->GetXaxis()->SetTitleFont(42);
        numerator_fullstaterr->GetYaxis()->SetLabelSize(0.18*0.045/0.06);
        numerator_fullstaterr->GetYaxis()->SetTitleSize(0.15);
        numerator_fullstaterr->GetYaxis()->SetTitleOffset(0.25);
        numerator_fullstaterr->GetYaxis()->SetTitleFont(42);
        
        numerator_fullstaterr->GetXaxis()->SetTitle(xtitlename);
        numerator_fullstaterr->GetYaxis()->SetTitle(ytitlename);
        //numerator_fullstaterr->SetFillColor(kGreen-3);
        numerator_fullstaterr->SetFillColor(kRed-10);
        numerator_fullstaterr->DrawCopy("");

        //
        // Drawing lines
        tline0->SetLineStyle(2);
        //tline0->Draw();

      }
      else {

      //
      // Plotting
      numerator->GetYaxis()->SetNdivisions(505);
      numerator->GetYaxis()->SetTickLength(0.015);
      numerator->GetXaxis()->SetTickLength(0.08);
      numerator->SetTitle("");
      ex1->Draw();
      numerator->DrawCopy();

      ex2->Draw();
      denominator->DrawCopy("e2same");
      //denominator->DrawCopy("same");

      TH1F *denominator_Clone = (TH1F*)denominator->Clone();
      //////////denominator_Clone->SetFillColor(kWhite);
      //////////denominator_Clone->Draw("hist same");

      ex1->Draw();
      numerator->DrawCopy("same");

      numerator->Print("all");
      denominator->Print("all");
      numerator_fullstaterr->Print("all");

      //
      // Drawing lines
      tline->SetLineStyle(2);
      //tline->Draw();

      }
      

      //
      if(option.find("QCD")==string::npos ){
  
        // Njet separation lines
        TLine *tl_njet = new TLine();
        tl_njet->SetLineStyle(2);
        tl_njet->DrawLine(31.-0.5,ymin_bottom,31.-0.5,ymax_bottom); 
        tl_njet->DrawLine(71.-0.5,ymin_bottom,71.-0.5,ymax_bottom); 
        tl_njet->DrawLine(111.-0.5,ymin_bottom,111.-0.5,ymax_bottom); 
        tl_njet->DrawLine(143.-0.5,ymin_bottom,143.-0.5,ymax_bottom); 

        // Nb separation lines
        TLine *tl_nb = new TLine();
        tl_nb->SetLineStyle(3);
        tl_nb->DrawLine(11.-0.5,ymin_bottom,11.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(21.-0.5,ymin_bottom,21.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(31.-0.5,ymin_bottom,31.-0.5,ymax2_bottom);
        tl_nb->DrawLine(41.-0.5,ymin_bottom,41.-0.5,ymax2_bottom);
        tl_nb->DrawLine(51.-0.5,ymin_bottom,51.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(61.-0.5,ymin_bottom,61.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(71.-0.5,ymin_bottom,71.-0.5,ymax2_bottom);
        tl_nb->DrawLine(81.-0.5,ymin_bottom,81.-0.5,ymax2_bottom);
        tl_nb->DrawLine(91.-0.5,ymin_bottom,91.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(101.-0.5,ymin_bottom,101.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(111.-0.5,ymin_bottom,111.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(119.-0.5,ymin_bottom,119.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(127.-0.5,ymin_bottom,127.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(135.-0.5,ymin_bottom,135.-0.5,ymax2_bottom); 
        tl_nb->DrawLine(143.-0.5,ymin_bottom,143.-0.5,ymax2_bottom);
        tl_nb->DrawLine(151.-0.5,ymin_bottom,151.-0.5,ymax2_bottom);
        tl_nb->DrawLine(159.-0.5,ymin_bottom,159.-0.5,ymax2_bottom);
        tl_nb->DrawLine(167.-0.5,ymin_bottom,167.-0.5,ymax2_bottom);

      } else {
  
      // Njet separation lines
      TLine *tl_njet = new TLine();
      tl_njet->SetLineStyle(2);
      tl_njet->DrawLine(40.-0.5,ymin_bottom,40.-0.5,ymax_bottom); 
      tl_njet->DrawLine(92.-0.5,ymin_bottom,92.-0.5,ymax_bottom); 
      tl_njet->DrawLine(144.-0.5,ymin_bottom,144.-0.5,ymax_bottom); 
      tl_njet->DrawLine(184.-0.5,ymin_bottom,184.-0.5,ymax_bottom);

      // Nb separation lines
      TLine *tl_nb = new TLine();
      tl_nb->SetLineStyle(3);
      tl_nb->DrawLine(14.-0.5,ymin_bottom,14.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(27.-0.5,ymin_bottom,27.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(40.-0.5,ymin_bottom,40.-0.5,ymax2_bottom);
      tl_nb->DrawLine(53.-0.5,ymin_bottom,53.-0.5,ymax2_bottom);
      tl_nb->DrawLine(66.-0.5,ymin_bottom,66.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(79.-0.5,ymin_bottom,79.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(92.-0.5,ymin_bottom,92.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(105.-0.5,ymin_bottom,105.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(118.-0.5,ymin_bottom,118.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(131.-0.5,ymin_bottom,131.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(144.-0.5,ymin_bottom,144.-0.5,ymax2_bottom);
      tl_nb->DrawLine(154.-0.5,ymin_bottom,154.-0.5,ymax2_bottom);
      tl_nb->DrawLine(164.-0.5,ymin_bottom,164.-0.5,ymax2_bottom);
      tl_nb->DrawLine(174.-0.5,ymin_bottom,174.-0.5,ymax2_bottom);
      tl_nb->DrawLine(184.-0.5,ymin_bottom,184.-0.5,ymax2_bottom);
      tl_nb->DrawLine(194.-0.5,ymin_bottom,194.-0.5,ymax2_bottom);
      tl_nb->DrawLine(204.-0.5,ymin_bottom,204.-0.5,ymax2_bottom);
      tl_nb->DrawLine(214.-0.5,ymin_bottom,214.-0.5,ymax2_bottom);
    
      }

      gPad->RedrawAxis();

      //
      //

          TLine* blue_line = new TLine() ;
          blue_line -> SetLineColor(4) ;
          blue_line -> DrawLine( 0.5, 1., nbins+0.5, 1. ) ;

  CMS_lumi(canvas, iPeriod, iPos, lumi_sqrtS);

  ////////// if(doDataVsMC){
  //////////   sprintf(tempname,"DataMC_Full_Plot.pdf");
  //////////   if (pull==1)    sprintf(tempname,"DataMCPull_Full_Plot.pdf");
  ////////// }else{
  //////////   if(doClosurewoIsoTrackVeto){
  //////////     if(option.find("QCD")!=string::npos) sprintf(tempname,"Closure_QCD_HDP_woIsoTrack_Full_Plot.pdf");
  //////////       else sprintf(tempname,"Closure_woIsoTrack_Full_Plot.pdf");
  //////////     if (pull==1)    sprintf(tempname,"ClosurePull_woIsoTrack_Full_Plot.pdf");
  //////////   }else{
  //////////     if(option.find("QCD")!=string::npos) sprintf(tempname,"Closure_QCD_HDP_Full_Plot.pdf");
  //////////       else sprintf(tempname,"Closure_Full_Plot.pdf");
  //////////     if (pull==1)    sprintf(tempname,"ClosurePull_Full_Plot.pdf");
  //////////   }
  ////////// }

  canvas->Print( "outputfiles/qcd-closure.pdf" );

  TH1F *h_ratioDist = 0;
  if(pull==0){
    h_ratioDist = new TH1F("h_ratioDist", "h_ratioDist", 39, 0.025, 1.975);

    for(int ibin=1; ibin<numerator->GetNbinsX()+1; ibin++){
      h_ratioDist->Fill(numerator->GetBinContent(ibin));
      h_ratioDist->GetYaxis()->SetRangeUser(0., h_ratioDist->GetBinContent(20)*1.5);
    }
  }
  else{
    h_ratioDist = new TH1F("h_pullDist", "h_pullDist", 19, -2.85, 2.85);

    for(int ibin=1; ibin<numerator_fullstaterr->GetNbinsX()+1; ibin++){
      h_ratioDist->Fill(numerator_fullstaterr->GetBinContent(ibin));
      h_ratioDist->GetYaxis()->SetRangeUser(0., h_ratioDist->GetBinContent(10)*1.5);
    }
  }

  gROOT->SetBatch(true);
  gStyle->SetOptStat(2211);
  gStyle->SetStatW(0.5);
  gStyle->SetStatH(0.5);

  TString name_(h_ratioDist->GetName());
  TCanvas *c1 = new TCanvas(name_,h_ratioDist->GetTitle(),1);
  //c1->SetLogy();
  //////////////////c1->cd();
  //////////////////h_ratioDist->Draw();
  //////////////////c1->SaveAs(name_+".pdf");
  delete c1;
  gROOT->SetBatch(false);

}


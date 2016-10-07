#include "slim-code/run_slimskim.c"

#include "data_turnon1.c"
#include "fill_hists_loop_v2d.c"
#include "fill_data_hists_loop_v2d.c"
#include "make_qcdmc_input_files1.c"

#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include "modelfit3.c"
#include "make_data_input_files1.c"
#include "make_lostlep_input_files1.c"
#include "make_hadtau_input_files1.c"
#include "make_znunu_input_files1.c"
#include "syst_2015_v2.c"
#include "draw_qcd_ratio_v3.c"
#include "draw_badjet_cat_v3.c"
#include "gen_modelfit_input1.c"
#include "run_modelfit3_on_data.c"
#include "create_model_ratio_hist1.c"

void run_all ( TString skim_slim_input_dir = "" )

{

   //Adding necessary dictionaries

   gROOT->ProcessLine(".L loader.C+"); 

   //Run the slim code
   if ( skim_slim_input_dir != "" )
   {

      run_slimskim(skim_slim_input_dir + "/tree_signal");
      run_slimskim(skim_slim_input_dir+ "/tree_LDP");

      gSystem -> Exec("mkdir -p ./fnal-prod-v9-skims-slimmed/tree_signal");
      gSystem -> Exec("mkdir -p ./fnal-prod-v9-skims-slimmed/tree_LDP");

      gSystem -> Exec("mv " + skim_slim_input_dir + "/tree_signal/slim/* ./fnal-prod-v9-skims-slimmed/tree_signal");
      gSystem -> Exec("mv " + skim_slim_input_dir + "/tree_LDP/slim/* ./fnal-prod-v9-skims-slimmed/tree_LDP");

   }
   data_turnon1();

   fill_hists_loop_v2d f1;
   f1.Loop();

   fill_data_hists_loop_v2d f2;
   f2.Loop();

   make_qcdmc_input_files1();
   modelfit3();

   make_data_input_files1();
   make_lostlep_input_files1();
   make_hadtau_input_files1();
   make_znunu_input_files1();

   gen_modelfit_input1();
   run_modelfit3_on_data();
   syst_2015_v2 f3;
   f3.Loop();
   draw_qcd_ratio_v3();
   draw_badjet_cat_v3();
   create_model_ratio_hist1();
}

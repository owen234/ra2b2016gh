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

}

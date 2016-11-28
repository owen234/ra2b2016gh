#include <stdio.h>
#include <sys/stat.h>

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
#include "make_lostlep_input_files2.c"
#include "make_hadtau_input_files2.c"
#include "make_znunu_input_files2.c"
#include "syst_2015_v2.c"
#include "draw_badjet_cat_v3.c"
#include "gen_modelfit_input1.c"
#include "run_modelfit3_on_data.c"
#include "gen_combine_input2.c"
#include "dump_qcdmc_vals.c"
#include "closure_sums3.c"
#include "draw_closure_sums1.c"
#include "create_model_ratio_hist1.c"
#include "closure_v4.c"
#include "draw_qcd_ratio_v3.c"
#include "draw_model_vs_mc.c"

bool does_outputfiles_exist();


void run_all ( TString skim_slim_input_dir = "" )

{

   //Adding necessary dictionaries

   gROOT->ProcessLine(".L loader.C+"); 

   //Run the slim code
   if ( skim_slim_input_dir != "" )
   {

      run_slimskim(skim_slim_input_dir + "/tree_signal");
      run_slimskim(skim_slim_input_dir+ "/tree_LDP");

      gSystem -> Exec("mkdir -p ./fnal-prod-v10-skims-slimmed/tree_signal");
      gSystem -> Exec("mkdir -p ./fnal-prod-v10-skims-slimmed/tree_LDP");

      gSystem -> Exec("mv " + skim_slim_input_dir + "/tree_signal/slim/* ./fnal-prod-v10-skims-slimmed/tree_signal");
      gSystem -> Exec("mv " + skim_slim_input_dir + "/tree_LDP/slim/* ./fnal-prod-v10-skims-slimmed/tree_LDP");

   }

   if ( !does_outputfiles_exist() ) gSystem -> Exec("mkdir ./outputfiles");

   std::cout << "\n-----------   running code data_turnon1   -----------\n" << std::endl;
   data_turnon1();

   std::cout << "\n-----------   running code fill_hists_loop_v2d   -----------\n" << std::endl;
   fill_hists_loop_v2d f1;
   f1.Loop();

   std::cout << "\n-----------   running code fill_data_hists_loop_v2d   -----------\n" << std::endl;
   fill_data_hists_loop_v2d f2;
   f2.Loop();


   std::cout << "\n-----------   running code make_qcdmc_input_files1   -----------\n" << std::endl;
   make_qcdmc_input_files1();
   modelfit3();

   std::cout << "\n-----------   running code make_data_input_files1   -----------\n" << std::endl;
   make_data_input_files1();

   std::cout << "\n-----------   running code make_lostlep_input_files2   -----------\n" << std::endl;
   make_lostlep_input_files2();

   std::cout << "\n-----------   running code make_hadtau_input_files2   -----------\n" << std::endl;
   make_hadtau_input_files2();

   std::cout << "\n-----------   running code make_znunu_input_files2   -----------\n" << std::endl;   
   make_znunu_input_files2();

   std::cout << "\n-----------   running code gen_modelfit_input1   -----------\n" << std::endl;    
   gen_modelfit_input1();

   std::cout << "\n-----------   running code run_modelfit3_on_data   -----------\n" << std::endl;
   run_modelfit3_on_data();

   std::cout << "\n-----------   running code syst_2015_v2   -----------\n" << std::endl;   
   syst_2015_v2 f3;
   f3.Loop();

   std::cout << "\n-----------   running code draw_qcd_ratio_v3   -----------\n" << std::endl;   
   draw_qcd_ratio_v3();

   std::cout << "\n-----------   running code draw_badjet_cat_v3   -----------\n" << std::endl;   
   draw_badjet_cat_v3();

   std::cout << "\n-----------   running code create_model_ratio_hist1   -----------\n" << std::endl;   
   create_model_ratio_hist1();

   std::cout << "\n-----------   running code gen_combine_input2   -----------\n" << std::endl;   
   gen_combine_input2();

   std::cout << "\n-----------   running code dump_qcdmc_vals   -----------\n" << std::endl;
   dump_qcdmc_vals();

   std::cout << "\n-----------   running code closure_sums3   -----------\n" << std::endl;   
   closure_sums3();

   std::cout << "\n-----------   running code draw_closure_sums1   -----------\n" << std::endl;
   draw_closure_sums1("njet");
   draw_closure_sums1("mht");
   draw_closure_sums1("ht");
   draw_closure_sums1("nb");
   draw_closure_sums1("10boxes");

   std::cout << "\n-----------   running code closure_v4   -----------\n" << std::endl;   

   closure_v4();

   std::cout << "\n-----------   running code draw_model_vs_mc   -----------\n" << std::endl;

   draw_model_vs_mc();

}




bool does_outputfiles_exist()
{
    const char* folderr;
    folderr = "./outputfiles";
    struct stat sb;

    if (stat(folderr, &sb) == 0 && S_ISDIR(sb.st_mode))
    {
       return true;
    }
    else
    {
       return false;
    }
}// does_outputfiles_exist

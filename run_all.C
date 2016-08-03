#include "slim-code/run_slimskim.c"
#include "data_turnon1.c"

#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"



void run_all (bool do_skim_slim = false )

{

   //Adding necessary dictionaries
   gROOT->ProcessLine(".L loader.C+");

   //Run the slim code
   if ( do_skim_slim == true )
   {

      run_slimskim("./fnal-prod-v9-skims/tree_signal");
      run_slimskim("./fnal-prod-v9-skims/tree_LDP");

      gSystem -> Exec("mkdir -p ./fnal-prod-v9-skims-skimmed/tree_signal");
      gSystem -> Exec("mkdir -p ./fnal-prod-v9-skims-skimmed/tree_LDP");

      gSystem -> Exec("mv ./fnal-prod-v9-skims/tree_signal/slim/* ./fnal-prod-v9-skims-slimmed/tree_signal");
      gSystem -> Exec("mv ./fnal-prod-v9-skims/tree_LDP/slim/* ./fnal-prod-v9-skims-slimmed/tree_LDP");

   }
   data_turnon1();

}

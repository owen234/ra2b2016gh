#include "slim-code/run_slimskim.c"

#include "data_turnon1.c"
#include "fill_hists_loop_v2d.c"
#include "fill_data_hists_loop_v2d.c"

#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"

void run_all (bool do_skim_slim = false )

{

   //Adding necessary dictionaries
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class std::vector<bool>+;
#endif

   //Run the slim code
   if ( do_skim_slim == true )
   {

      run_slimskim("./fnal-prod-v9-skims/tree_signal");
      run_slimskim("./fnal-prod-v9-skims/tree_LDP");

      gSystem -> Exec("mkdir -p ./fnal-prod-v9-skims-slimmed/tree_signal");
      gSystem -> Exec("mkdir -p ./fnal-prod-v9-skims-slimmed/tree_LDP");

      gSystem -> Exec("mv ./fnal-prod-v9-skims/tree_signal/slim/* ./fnal-prod-v9-skims-slimmed/tree_signal");
      gSystem -> Exec("mv ./fnal-prod-v9-skims/tree_LDP/slim/* ./fnal-prod-v9-skims-slimmed/tree_LDP");

   }
   data_turnon1();

   fill_hists_loop_v2d f1;
   f1.Loop();

   fill_data_hists_loop_v2d f2;
   f2.Loop();


}

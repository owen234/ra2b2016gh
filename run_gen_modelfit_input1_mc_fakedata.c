
#include "TSystem.h"

#include "gen_modelfit_input1.c"

   void run_gen_modelfit_input1_mc_fakedata() {

      gSystem -> Exec( "mkdir -p outputfiles/files-modelfit-input-mc-fakedata" ) ;

      gen_modelfit_input1(  "outputfiles/mc-nbsum-input-fakedata.txt",
                            "outputfiles/mc-nbsum-input-lostlep.txt",
                            "outputfiles/mc-nbsum-input-hadtau.txt",
                            "outputfiles/mc-nbsum-input-znunu.txt",
                            "outputfiles/modelfit-input-mc-fakedata.root" ) ;

      gSystem -> Exec("mv outputfiles/modelfit-input*.* outputfiles/files-modelfit-input-mc-fakedata" ) ;

   }


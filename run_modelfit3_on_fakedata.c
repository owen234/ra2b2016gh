

#include "modelfit3.c"

   void run_modelfit3_on_fakedata() {

      modelfit3( "outputfiles/files-modelfit-input-mc-fakedata/modelfit-input-mc-fakedata.root",
                 "outputfiles/fakedata-chi2-fit", "h_ratio", -0.1, 0.8 ) ;

   }





#include "TSystem.h"
#include "modelfit3.c"

   void run_modelfit3_on_data() {

      modelfit3( "outputfiles/modelfit-input-data.root", "outputfiles/data-chi2-fit", "h_ratio", -0.1, 0.4 ) ;
      gSystem -> Exec("mv outputfiles/data-chi2-fit-modelfit.pdf outputfiles/data-chi2-fit-modelfit-range1.pdf" ) ;

      modelfit3( "outputfiles/modelfit-input-data.root", "outputfiles/data-chi2-fit", "h_ratio", -0.1, 0.8 ) ;
      gSystem -> Exec("mv outputfiles/data-chi2-fit-modelfit.pdf outputfiles/data-chi2-fit-modelfit-range2.pdf" ) ;

      modelfit3( "outputfiles/modelfit-input-data.root", "outputfiles/data-chi2-fit", "h_ratio", -0.1, 1.5 ) ;
      gSystem -> Exec("mv outputfiles/data-chi2-fit-modelfit.pdf outputfiles/data-chi2-fit-modelfit-range3.pdf" ) ;

   }



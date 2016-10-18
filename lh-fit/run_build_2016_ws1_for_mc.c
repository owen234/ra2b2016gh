
#include "build_2016_ws1.c"


void run_build_2016_ws1_for_mc() {

   build_2016_ws1( "outputfiles/ws-lhfit-test-mc.root",
                   "../outputfiles/mc-nbsum-input-lostlep.txt",
                   "../outputfiles/mc-nbsum-input-hadtau.txt",
                   "../outputfiles/mc-nbsum-input-znunu.txt",
                   "../outputfiles/mc-nbsum-input-fakedata.txt",
                   "../outputfiles/nbsum-input-T1bbbbH.txt" ) ;

}


#include "TH1.h"
#include "TFile.h"

#ifndef input_files_h
#define input_files_h


TH1F* get_hist( TFile* tf, const char* hname ) {

   TH1F* hp   = (TH1F*) tf -> Get( hname ) ;
   if ( hp == 0x0 ) {

      printf("\n\n *** Can't find hist %s.\n\n", hname ) ;
      gSystem -> Exit(-1) ;

   }

   return hp ;

} // get_hist

#endif

#ifndef binning_h
#define binning_h

#include <iostream>
#include "TSystem.h"
//
//  Organizing principles of binning.
//
//    For each dimension
//      bin index 0 means below bins used in analysis (e.g. HT<500 is HT bin 0).
//      bin index 1 is the first bin used
//      This also goes for btagging, so btag bin 1 is nb=0.
//      In the bin edges arrays, it goes like this
//        bin_edge_x[0] = 500. ; // low edge of bin 1.
//        bin_edge_x[1] = 750. ; // high edge of bin 1, low edge of bin 2.
//        ...
//


void htmht_bin_to_ht_and_mht_bins(int bin_htmht, int& bin_ht, int& bin_mht);
bool is_this_bin_excluded(int bin_nj, int bin_nb, int bin_ht, int bin_mht);
bool is_this_bin_excluded(int bin_nj, int bin_nb, int bin_htmht);

int global_bin_with_mhtc( int nji, int nbi, int htmhti ) ;
int global_search_bin ( int arg_nji, int arg_nbi, int arg_htmhti ) ;
bool is_this_bin_excluded(int bin_no);
void fill_gbi ();
bool translate_qcd_bin_to_nj_nb_ht_mht( int qcd_bin_no, int& arg_nj, int& arg_nb, int& arg_ht, int& arg_mht );
bool translate_search_bin_to_nj_nb_ht_mht( int search_bin_no, int& arg_nj, int& arg_nb, int& arg_ht, int& arg_mht );

using namespace std ;


   float bin_edges_nj[100]  ;
   float bin_edges_nb[100]  ;
   float bin_edges_mht[100]  ;
   float bin_edges_ht[100][100]  ;
   int   nb_nj ;
   int   nb_nb ;
   int   nb_ht[100] ;
   int   nb_mht ;

   int   nb_htmht ;
   int   nBinsHT = 3; 

   int   nb_global ;
   int   nb_global_after_exclusion ;   
   int   nb_global_w_exclusion_w_mhtc ;
   int   bi_global ;
   int   bi_nbsum_global ;

   int   njet_bin_to_be_fixed_in_qcd_model_fit;

   int   no_bin_ht  [10] = {};
   int   no_bin_mht [10] = {};
   int   no_bin_bjet[10] = {};
   int   no_bin_njet[10] = {};

   int   no_bin_ht_w_exclusion_w_mhtc  [10] = {};
   int   no_bin_mht_w_exclusion_w_mhtc [10] = {};
   int   no_bin_bjet_w_exclusion_w_mhtc[10] = {};
   int   no_bin_njet_w_exclusion_w_mhtc[10] = {};

   bool  gbi_array_ready(false) ;
   int   gbi_with_mhtc[6][5][14] ;
   int   gbi_search_bins  [6][5][14] ;

//=====================================================================================================

   void setup_bins() {

      int bi ;

      bi = 0 ;

//      bin_edges_nj[bi] = 1.5 ; bi++ ; // *** adding NJets=2 binning
      njet_bin_to_be_fixed_in_qcd_model_fit  = bi; // Setting bin NJets = {3,4} as the bin to be fixed to one while fitting QCD toymodel parameters
      bin_edges_nj[bi] = 2.5 ; bi++ ;
      bin_edges_nj[bi] = 4.5 ; bi++ ;
      bin_edges_nj[bi] = 6.5 ; bi++ ;
      bin_edges_nj[bi] = 8.5 ; bi++ ;
      bin_edges_nj[bi] = 99.5 ;
      nb_nj = bi ;

      bi = 0 ;
      bin_edges_nb[bi] = -0.5 ; bi++ ;
      bin_edges_nb[bi] =  0.5 ; bi++ ;
      bin_edges_nb[bi] =  1.5 ; bi++ ;
      bin_edges_nb[bi] =  2.5 ; bi++ ;
      bin_edges_nb[bi] =  99.5 ;
      nb_nb = bi ;

      bi = 0 ;
      bin_edges_mht[bi] =   250. ; bi++ ;
      bin_edges_mht[bi] =   300. ; bi++ ;
      bin_edges_mht[bi] =   350. ; bi++ ;
      bin_edges_mht[bi] =   500. ; bi++ ;
      bin_edges_mht[bi] =   750. ; bi++ ;
      bin_edges_mht[bi] = 20000. ;
      nb_mht = bi ;

      int mbi(1) ;
      nb_htmht = 0 ;
      //--- MHT bin 0
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   300. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 1
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   300. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 2
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   350. ; bi++ ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 3
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   500. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1000. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;
      mbi++ ;
      //--- MHT bin 4
      bi = 0 ;
      bin_edges_ht[mbi][bi] =   750. ; bi++ ;
      bin_edges_ht[mbi][bi] =  1500. ; bi++ ;
      bin_edges_ht[mbi][bi] = 20000. ;
      nb_ht[mbi] = bi ;
      nb_htmht += bi ;


      nb_global = nb_htmht * nb_nb * nb_nj ;
      nb_global_after_exclusion = 0;
      nb_global_w_exclusion_w_mhtc = 0;

      for ( int bin_nj = 0; bin_nj < nb_nj;    bin_nj++)
         for ( int bin_nb = 0;    bin_nb < nb_nb;    bin_nb++)
            for ( int bin_htmht = 1; bin_htmht <= nb_htmht; bin_htmht++)
	    {

               int bin_ht, bin_mht;
	       htmht_bin_to_ht_and_mht_bins (bin_htmht, bin_ht, bin_mht);
               if ( !is_this_bin_excluded(bin_nj, bin_nb, bin_ht-1, bin_mht-1 ) ) 
	       {
		  if ( bin_htmht > 3) 
		  {
                     nb_global_after_exclusion++;
                     no_bin_ht  [bin_ht-1  ]++;
                     no_bin_mht [bin_mht-1 ]++;		  
                     no_bin_bjet[bin_nb    ]++;
                     no_bin_njet[bin_nj    ]++;

		  }

                  nb_global_w_exclusion_w_mhtc++;
                  no_bin_ht_w_exclusion_w_mhtc  [bin_ht-1  ]++;
                  no_bin_mht_w_exclusion_w_mhtc [bin_mht-1 ]++;
                  no_bin_bjet_w_exclusion_w_mhtc[bin_nb    ]++;
                  no_bin_njet_w_exclusion_w_mhtc[bin_nj    ]++;

	       }//if is_this_bin_excluded
            }//bin_nj

      fill_gbi();

   } // setup_bins


bool is_this_bin_excluded(int bin_nj, int bin_nb, int bin_ht, int bin_mht)
{
  // all variables start from zero
  // bin_mht = 0 corresponds to MHTC bins
/*
  if ( bin_edges_nj[0] == 1.5 ) // If NJets = 2 binning is considered
  {	
     if ( bin_nj == 0 && bin_nb == 3 ) return true; // Removing bins NJets0_BTags3_...
  } 
  if ( bin_nj == nb_nj-2 && bin_ht == 0 && bin_mht == 0) return true; // Removing bins NJets3_BTagsX_MHTC_HT0
  if ( bin_nj == nb_nj-2 && bin_ht == 0 && bin_mht == 1) return true; // Removing bins NJets3_BTagsX_MHT0_HT0
  if ( bin_nj == nb_nj-2 && bin_ht == 0 && bin_mht == 2) return true; // Removing bins NJets3_BTagsX_MHT1_HT0

  if ( bin_nj == nb_nj-1 && bin_ht == 0 && bin_mht == 0) return true; // Removing bins NJets4_BTagsX_MHTC_HT0
  if ( bin_nj == nb_nj-1 && bin_ht == 0 && bin_mht == 1) return true; // Removing bins NJets4_BTagsX_MHT0_HT0
  if ( bin_nj == nb_nj-1 && bin_ht == 0 && bin_mht == 2) return true; // Removing bins NJets4_BTagsX_MHT1_HT0
*/
  return false;  //don't exclude this bin
}


bool is_this_bin_excluded(int bin_nj, int bin_nb, int bin_htmht)
{
   int bin_ht, bin_mht;
   htmht_bin_to_ht_and_mht_bins (bin_htmht+1, bin_ht, bin_mht);
   return is_this_bin_excluded(bin_nj, bin_nb, bin_ht-1, bin_mht-1);

}

bool is_this_bin_excluded(int bin_no){
      int gbi(0) ;
      for ( int nji=1; nji<=nb_nj; nji++ ) { 
         for ( int nbi=1; nbi<=nb_nb; nbi ++ ) {
            for ( int htmhti=1; htmhti<=nb_htmht; htmhti++ ) {
               gbi++;
               if ( gbi == bin_no ){
                  return is_this_bin_excluded( nji-1, nbi-1, htmhti-1 ) ;
               }//if ( gbi == bin_no )
	    } //htmhti
         } // nbi
      } // mhi
      return false;
}  //is_this_bin_excluded


void htmht_bin_to_ht_and_mht_bins(int bin_htmht, int& bin_ht, int& bin_mht)
{

   if ( bin_htmht < 1 || bin_htmht > nb_htmht )
   {
      std::cout << "Fatal error: ht_mht bin number out of range: " << bin_htmht << std::endl;
      gSystem -> Exit(-1) ;
   }
  
   if ( bin_htmht == 1 ) { bin_ht = 1; bin_mht = 1; }
   if ( bin_htmht == 2 ) { bin_ht = 2; bin_mht = 1; }
   if ( bin_htmht == 3 ) { bin_ht = 3; bin_mht = 1; }

   if ( bin_htmht == 4 ) { bin_ht = 1; bin_mht = 2; }
   if ( bin_htmht == 5 ) { bin_ht = 2; bin_mht = 2; }
   if ( bin_htmht == 6 ) { bin_ht = 3; bin_mht = 2; }

   if ( bin_htmht == 7 ) { bin_ht = 1; bin_mht = 3; }
   if ( bin_htmht == 8 ) { bin_ht = 2; bin_mht = 3; }
   if ( bin_htmht == 9 ) { bin_ht = 3; bin_mht = 3; }

   if ( bin_htmht ==10 ) { bin_ht = 2; bin_mht = 4; }
   if ( bin_htmht ==11 ) { bin_ht = 3; bin_mht = 4; }

   if ( bin_htmht ==12 ) { bin_ht = 2; bin_mht = 5; }
   if ( bin_htmht ==13 ) { bin_ht = 3; bin_mht = 5; }

}
//=========================================

void fill_gbi ()
{

   if ( !gbi_array_ready ) {
      int gbi(0) ;
      int gbi_no_mhtc(0) ;      
      for ( int nji=1; nji<=nb_nj; nji++ ) {
         for ( int nbi=1; nbi<=nb_nb; nbi ++ ) {
            for ( int htmhti=1; htmhti<=nb_htmht; htmhti++ ) {
               int hti, mhti ;
               htmht_bin_to_ht_and_mht_bins( htmhti, hti, mhti ) ;
               bool excluded = is_this_bin_excluded( nji-1, nbi-1, hti-1, mhti-1 ) ;
               if ( !excluded ) {
                  gbi++ ;
                  gbi_with_mhtc[nji][nbi][htmhti] = gbi ;

                  if ( htmhti > 3 )
		  {
                     gbi_no_mhtc++ ;
                     gbi_search_bins[nji][nbi][htmhti] = gbi_no_mhtc ;
		  }//if htmhti
		  else                                   //else for if ( htmhti > 3 )
                     gbi_search_bins[nji][nbi][htmhti] = -1 ;

	       } else {                                 // else for if ( !excluded ) {
                  gbi_with_mhtc[nji][nbi][htmhti] = -1 ;
                  gbi_search_bins[nji][nbi][htmhti] = -1;
               } // else if !excluded
            } //htmhti
         } // nbi
      } // mhi
      gbi_array_ready = true ;
   }//if (!gbi_array_ready)

}//fill_gbi

int global_bin_with_mhtc ( int arg_nji, int arg_nbi, int arg_htmhti ) {

//all variables start from one
 
   if ( arg_nji    < 1 || arg_nji    > nb_nj    ) {std::cout << "Error: Njet argument out of range"   << std::endl; return -1; }
   if ( arg_nbi    < 1 || arg_nbi    > nb_nb    ) {std::cout << "Error: Nbjet argument out of range"  << std::endl; return -1; }
   if ( arg_htmhti < 1 || arg_htmhti > nb_htmht ) {std::cout << "Error: HT_MHT argument out of range" << std::endl; return -1; }


   if ( !gbi_array_ready ) {
      fill_gbi ();
      gbi_array_ready = true ;
   }

   return gbi_with_mhtc[arg_nji][arg_nbi][arg_htmhti] ;

} // global_bin_with_mhtc


int global_search_bin ( int arg_nji, int arg_nbi, int arg_htmhti ) {

//all variables start from one

   if ( arg_nji    < 1 || arg_nji    > nb_nj    ) {std::cout << "Error: Njet argument out of range"   << std::endl; return -1; }
   if ( arg_nbi    < 1 || arg_nbi    > nb_nb    ) {std::cout << "Error: Nbjet argument out of range"  << std::endl; return -1; }
   if ( arg_htmhti < 1 || arg_htmhti > nb_htmht ) {std::cout << "Error: HT_MHT argument out of range" << std::endl; return -1; }


   if ( !gbi_array_ready ) {
         fill_gbi ();
       gbi_array_ready = true ;
   }
   return gbi_search_bins[arg_nji][arg_nbi][arg_htmhti] ;

} // global_bin_with_mhtc



bool translate_search_bin_to_nj_nb_ht_mht( int search_bin_no, int& arg_nj, int& arg_nb, int& arg_ht, int& arg_mht )
{

   if ( search_bin_no < 1 || search_bin_no > nb_global_after_exclusion ) {std::cout << "Error: Search bin index out of range" << std::endl; return -1; }

   for ( int nji=1; nji<=nb_nj; nji++ )
   {
      for ( int nbi=1; nbi<=nb_nb; nbi ++ )
      {
         for ( int htmhti=1; htmhti<=nb_htmht; htmhti++ )
         {
            if ( gbi_search_bins[nji][nbi][htmhti] == search_bin_no )
            {
               arg_nj = nji;
               arg_nb = nbi;
               htmht_bin_to_ht_and_mht_bins(htmhti,arg_ht,arg_mht);
               return 1; //successfully found the bin
            }
         }//htmhti
      }//nbi
   }//nji

   arg_nj  = -1;
   arg_nb  = -1;
   arg_ht  = -1;
   arg_mht = -1;
   return 0; // Didn't find such a bin
}// translate_search_bin_to_nj_nb_ht_mht



bool translate_qcd_bin_to_nj_nb_ht_mht( int qcd_bin_no, int& arg_nj, int& arg_nb, int& arg_ht, int& arg_mht )
{

   if ( qcd_bin_no < 1 || qcd_bin_no > nb_global_w_exclusion_w_mhtc ) {std::cout << "Error: QCD bin index out of range" << std::endl; return -1; }

   for ( int nji=1; nji<=nb_nj; nji++ )
   {
      for ( int nbi=1; nbi<=nb_nb; nbi ++ )
      {
         for ( int htmhti=1; htmhti<=nb_htmht; htmhti++ )
         {
            if ( gbi_with_mhtc[nji][nbi][htmhti] == qcd_bin_no )
            {
               arg_nj = nji;
               arg_nb = nbi;
               htmht_bin_to_ht_and_mht_bins(htmhti,arg_ht,arg_mht);
               return 1; //successfully found the bin
            }
         } //htmhti
      } //nbi
   } //nji
   arg_nj  = -1;
   arg_nb  = -1;
   arg_ht  = -1;
   arg_mht = -1;
   return 0; // Didn't find such a bin
}//translate_qcd_bin_to_nj_nb_ht_mht


TString num_to_str(int value)
{
   
   TString str;
   str.Form("%d",value);
   return str;

}

TString num_to_str(double value)
{

   TString str; 
   str.Form("%f",value);
   return str;

}

TString num_to_str(float value)
{

   TString str;
   str.Form("%f",value);
   return str;

}


//=========================================




#endif

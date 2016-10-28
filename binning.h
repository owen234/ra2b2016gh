#ifndef binning_h
#define binning_h

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


void htmht_bin_to_ht_mht_bins(int bin_htmht, int& bin_ht, int& bin_mht);
bool exclude_this_bin(int bin_nj, int bin_nb, int bin_ht, int bin_mht);

using namespace std ;


   float bin_edges_nj[100]  ;
   float bin_edges_nb[100]  ;
   float bin_edges_mht[100]  ;
   float bin_edges_ht[100][100]  ;
   int   nb_nj ;
   int   nb_nb ;
   int   nb_ht[100] ;
   int   nb_mht ;
   int   bi_nj ;
   int   bi_nb ;
   int   bi_ht ;
   int   bi_mht ;

   int   nb_htmht ;
   int   bi_htmht ;
   int   nBinsHT = 3; 

   int   nb_global ;
   int   nb_global_after_exclusion ;   
   int   nb_global_w_mhtc ;
   int   bi_global ;
   int   bi_nbsum_global ;

   int   njet_bin_to_fix_in_qcd_model_fit;


   int   no_bin_ht  [10]{};
   int   no_bin_mht [10]{};
   int   no_bin_bjet[10]{};
   int   no_bin_njet[10]{};

   int   no_bin_ht_w_mhtc  [10]{};
   int   no_bin_mht_w_mhtc [10]{};
   int   no_bin_bjet_w_mhtc[10]{};
   int   no_bin_njet_w_mhtc[10]{};


//=====================================================================================================

   void setup_bins() {

      int bi ;

      bi = 0 ;

      bin_edges_nj[bi] = 1.5 ; bi++ ; // *** adding NJets=2 binning
      njet_bin_to_fix_in_qcd_model_fit  = bi; // Setting bin NJets = {3,4} as the bin to be fixed to one while fitting QCD toymodel parameters
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
      nb_global_w_mhtc = 0;

      for ( int bin_htmht = 1; bin_htmht <= nb_htmht; bin_htmht++)
         for ( int bin_nb = 0;    bin_nb < nb_nb;    bin_nb++)
            for ( int bin_nj = 0; bin_nj < nb_nj;    bin_nj++)
	    {

               int bin_ht, bin_mht;
	       htmht_bin_to_ht_mht_bins (bin_htmht, bin_ht, bin_mht);
               if ( !exclude_this_bin(bin_nj, bin_nb, bin_ht-1, bin_mht-1 ) ) 
	       {
		  if ( bin_htmht > 3) 
		  {
                     nb_global_after_exclusion++;
                     no_bin_ht  [bin_ht-1  ]++;
                     no_bin_mht [bin_mht-1 ]++;		  
                     no_bin_bjet[bin_nb    ]++;
                     no_bin_njet[bin_nj    ]++;

		  }
		  else
		  {
                     nb_global_w_mhtc++;
                     no_bin_ht_w_mhtc  [bin_ht-1  ]++;
                     no_bin_mht_w_mhtc [bin_mht-1 ]++;
                     no_bin_bjet_w_mhtc[bin_nb    ]++;
                     no_bin_njet_w_mhtc[bin_nj    ]++;
		  }
	       }//if exclude_this_bin
            }//bin_nj

   } // setup_bins


bool exclude_this_bin(int bin_nj, int bin_nb, int bin_ht, int bin_mht)
{
  // all variables start from zero
  // bin_mht 0 corresponds to MHTC bins

if ( bin_edges_nj[0] == 1.5 ) 
  if ( bin_nj == 0 && bin_nb == 3 ) return true; // Removing bins NJets0_BTags3_...
  
  if ( bin_nj == nb_nj-2 && bin_ht == 0 && bin_mht == 0) return true; // Removing bins NJets3_BTags0_MHTC_HT0
  if ( bin_nj == nb_nj-2 && bin_ht == 0 && bin_mht == 1) return true; // Removing bins NJets3_BTags0_MHT0_HT0
  if ( bin_nj == nb_nj-2 && bin_ht == 0 && bin_mht == 2) return true; // Removing bins NJets3_BTags0_MHT1_HT0

  if ( bin_nj == nb_nj-1 && bin_ht == 0 && bin_mht == 0) return true; // Removing bins NJets4_BTags0_MHTC_HT0
  if ( bin_nj == nb_nj-1 && bin_ht == 0 && bin_mht == 1) return true; // Removing bins NJets4_BTags0_MHT0_HT0
  if ( bin_nj == nb_nj-1 && bin_ht == 0 && bin_mht == 2) return true; // Removing bins NJets4_BTags0_MHT1_HT0

  return false;  //don't exclude this bin
}

void htmht_bin_to_ht_mht_bins(int bin_htmht, int& bin_ht, int& bin_mht)
{
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

#endif

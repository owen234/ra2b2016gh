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
   int   bi_global ;
   int   bi_nbsum_global ;


//=====================================================================================================

   void setup_bins() {

      int bi ;

      bi = 0 ;

      bin_edges_nj[bi] = 1.5 ; bi++ ; // *** adding NJets=2 binning
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

   } // setup_bins

#endif

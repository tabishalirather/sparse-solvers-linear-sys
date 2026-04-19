/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
*/

#include <stdlib.h>
#include "itsolver.h"

csr *csr_testmatrices(int number)
{
  int i;
  csr *new;
  
  int rowa[11] = {0,3,5,8,11,13,16,19,22,25,28};
  int cola[28] = {0,      4,      8,
                    1,  3,
                      2,3,    6,
                    1,2,3,
                  0,      4,
                            5,  7,  9,
                      2,      6,7,
                            5,6,7,
                  0,              8,9,
                            5,    8,9};
  
  int rowb[16] = {0,3,8,12,16,21,24,29,36,41,46,53,60,65,70,75};
  int colb[75] = {0,          6,  8,
	  	    1,          7,8,  10,   12,
		      2,      6,7,  9,
       	                3,               11,12,   14,
                          4,        9,10,11,   13,
                            5,                 13,14,
	  	  0,  2,      6,7,8,
                    1,2,      6,7,8,9,10,
                  0,1,        6,7,8,
                      2,  4,    7,  9,10,
		    1,    4,    7,  9,10,11,12,
		        3,4,          10,11,12,13,14,
                    1,  3,            10,11,12,
                          4,5,           11,   13,14,
                        3,  5,           11,   13,14};
		      
  int rowc[10] = {0,9,11,13,15,17,19,21,23,25};
  int colc[25] = {0,1,2,3,4,5,6,7,8,
                  0,1,
                  0,  2,
                  0,    3,
                  0,      4,
                  0,        5,
                  0,          6,
                  0,            7,
                  0,              8};


  /* Bsp. für zerfallenden Matrixgraph */
  int rowd[9] = {0,2,5,8,11,14,17,20,22};
  int cold[22] = {0,7,1,4,6,2,3,5,2,3,5,1,4,6,2,3,5,1,4,6,0,7};

  switch (number)
  {
  case 0: /* Abbildung 2 */
    new = csr_alloc(10, 10, 28);
    for (i=0; i<11; i++) new->row[i] = rowa[i];
    for (i=0; i<28; i++)
    {
      new->col[i] = cola[i];
      new->coeff[i] = 1.0;
    }
    break;  
  
  case 1: /* Saad, Figure 2.10 / Abbildung 3 links */
    new = csr_alloc(15, 15, 75);
    for (i=0; i<16; i++) new->row[i] = rowb[i];
    for (i=0; i<75; i++)
    {
      new->col[i] = colb[i];
      new->coeff[i] = 1.0;
    }
    break;
    
  case 2: /* Saad, Figure 3.4 / Abbildung 3 rechts*/
    new = csr_alloc(9, 9, 25);
    for (i=0; i<10; i++) new->row[i] = rowc[i];
    for (i=0; i<25; i++)
    {
      new->col[i] = colc[i];
      new->coeff[i] = 1.0;
    }
    break;
  case 3:
    new = csr_alloc(8,8,22);
    for (i=0; i<9; i++) new->row[i] = rowd[i];
    for (i=0; i<22; i++)
    {
      new->col[i] = cold[i];
      new->coeff[i] = 1.0;
    }
    break;
    
  default:
    new = csr_alloc(1, 1, 0);
  }
  return new;

} /* csr_testmatrices */

/*! \file
 *  \brief Test program for computing a transposed sparse matrix.
 *
 *   \author Michael Wende
 *   \author Jens-Peter Zemke
 *   \author Sabine Le Borne
 */

#include "../Library/itsolver.h"

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
  csr *A, *AT;

  /*
         / 0 3 1 0 0 \
     A = | 1 0 4 0 5 |
         \ 0 3 2 4 0 /
   */

  int row[4] = {0,2,5,8};
  int col[8] = {1,2,0,2,4,1,2,3};
  double coeff[8] = {3.,1.,1.,4.,5.,3.,2.,4.};

  int i;

  A = csr_alloc(3, 5, 8);

  for (i=0; i<A->nr+1; i++)
    A->row[i] = row[i];
  for (i=0; i<A->nz; i++)
  {
    A->col[i] = col[i];
    A->coeff[i] = coeff[i];
  }

  /*  display in coordinate format and plot (original matrix)   */
csr_ascii(A);
// csr_spy(A, "A.eps");
  /* computing the transposed matrix */
//AT = csr_transpose;

  /*  display in coordinate format and plot (transposed matrix)  */

  csr_free(A);
  csr_free(AT);

  return EXIT_SUCCESS;
} /* main */

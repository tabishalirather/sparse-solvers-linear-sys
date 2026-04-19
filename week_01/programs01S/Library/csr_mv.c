/*! \file
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */


#include "itsolver.h"

void csr_mv(csr *A, double alpha, const double *x, double beta, double *y)
{
  printf("print by Tabish, Library/from csr_mv.c\n");
  int i, k;
  for (i = 0; i < A->nr; i++) {
    // Set y = beta * y allowing uninitialized y in case of beta = 0.
    if (beta == 0.0) {
      *y = 0.0;
    } else {
      *y *= beta;
    }
    // go through every entry of a given row.
    for (k = *A->row++; k < *A->row; k++) {
      *y += alpha * *A->coeff++ * x[*A->col++]; // *A->col gives the correct entry of x for a coefficient of A
    }
    y++;
  }

  /* Move pointers to the first entry. */
  A->row -= A->nr;
  A->col -= A->nz;
  A->coeff -= A->nz;

  /* No need to move the pointer to y back to the first entry,
     as this is only a copy. */
}


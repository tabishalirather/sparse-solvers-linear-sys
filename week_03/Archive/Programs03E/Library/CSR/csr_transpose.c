/*! \file
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <stdio.h>
#include <stdlib.h>
#include "itsolver.h"

csr *csr_transpose(csr *A)
{
  int i, j, q;

  csr *B = malloc(sizeof(csr));
  B->nr = A->nc;
  B->nc = A->nr;
  B->nz = A->nz;
  B->coeff = malloc(A->nz * sizeof(double));
  B->col = malloc(A->nz * sizeof(int));
  B->row = calloc(A->nc + 1,  sizeof(int));

  // Calculate how many entries are in each row and not calculating the pointer for the next row but the one after that
  for (j = 0; j < A->nz; j++) {
    if (A->col[j] + 1 < A->nc) {
      B->row[A->col[j] + 2]++;
    }
  }

  // B->row[i] only gives how many elements are in the row i-2
  for (i = 2; i < A->nc; i++) {
    B->row[i + 1] += B->row[i];
  }

  for (i = 0; i < A->nr; i++) { // go over every row
    for (j = A->row[i]; j < A->row[i + 1]; j++) { // go over every entry of one row, p is overall index of entry
      q = B->row[A->col[j] + 1]++; // A->col[j] column index of entry j --> pointer on the beginn of the following row is moved one up, q is the overall index for B corresponding to j
      B->col[q] = i;
      B->coeff[q] = A->coeff[j];
    }
  }

  return B;
}

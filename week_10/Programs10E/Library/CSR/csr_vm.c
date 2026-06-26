/*! \file
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <stdio.h>
#include <stdlib.h>
#include "itsolver.h"

void csr_vm(csr *A, double alpha, double *x, double beta, double *y)
{
  int i, k;

  for(i=0; i<A->nc; i++)
    *y++ *= beta;

  y -= A->nc;

  for (i=0; i<A->nr; i++)
  {
    for (k=*A->row++; k<*A->row; k++)
    {
      y[*A->col++] += alpha * *A->coeff++ * x[i];
    }
  }

  /* move pointers to the first entry */

  A->row -= A->nr;
  A->col -= A->nz;
  A->coeff -= A->nz;

} /* csr_vm */

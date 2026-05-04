/*! \file
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "itsolver.h"

#define max(A,B) ((A) > (B) ? (A) : (B))

void csr_ascii(csr * A)
{
  const int *Arow;
  const int *Acol;
  const double *coeff;
  int i, r;

  assert(A != NULL);

  printf("Matrix output: %d rows and %d columns, nnz = %d.\n",
     A->nr,A->nc,A->nz);
  printf("Average of %2.2f non-zero entries per row.\n",
    (1.0*A->nz)/(1.0*A->nr));
  printf("Portion of non-zero entries: %2.2f%%.\n",
    (100.0*A->nz)/(1.0*A->nr*A->nc));
  printf("\n");


  printf("Commented coordinate format:\n");

  Arow = A->row;
  Acol = A->col;
  coeff = A->coeff;

  for(i=0; i<A->nr; i++)
  {
    for(r=Arow[i]; r<Arow[i+1]; r++)
    {
      printf("Row: %3d, Column: %3d, Value: %g\n",
                       i,       Acol[r],    coeff[r]);
    }
  }
  printf("\n");

} /* csr_ascii */

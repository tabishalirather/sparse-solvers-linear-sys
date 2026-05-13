/*! \file
 *  \brief Spyplots of a sparse matrix discretizing a Laplacian on unit cubes
 *         in \f$ d = 1, 2, 3 \f$ dimensions.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "itsolver.h"

int main() {
  int n, d;
  char *filename;
  csr *A;

  n = 5;
  A = csr_laplace1d(n);
  filename = (char *) malloc(100);
  sprintf(filename, "testprog_03_myLaplace_1D.eps");
  csr_ascii(A);
  csr_spy(A, filename, 's');
  free(filename);
  csr_free(A);

  for (d = 1; d <= 3; ++d) {
    A = csr_laplace(n, d);
    filename = (char *) malloc(100);
    sprintf(filename, "testprog_03_laplace_n%d_d%d.eps", n, d);
    csr_spy(A, filename, 's');
    free(filename);
    csr_free(A);
  }
  // free(A);

  return 0;
}

/*! \file
 *  \brief How to initialize a CSR matrix by hand.
 *
 *   \author Jens-Peter Zemke, 
 *   \author Sabine Le Borne
 */

#include "itsolver.h"

int main() {

  /* Number of rows, columns and nonzeros. */
  int nr = 3;
  int nc = 5;
  int nnz = 8;
  csr *A = csr_alloc(nr, nc, nnz);

  /* Fill in coefficients and column indices. */
  A->row[0] = 0;

  /* First row. */
  A->coeff[0] = 3.0;
  A->col[0] = 1;
  A->coeff[1] = 1.0;
  A->col[1] = 2;
  A->row[1] = 2;

  /* Second row. */
  A->coeff[2] = 1.0;
  A->col[2] = 0;
  A->coeff[3] = 4.0;
  A->col[3] = 2;
  A->coeff[4] = 5.0;
  A->col[4] = 4;
  A->row[2] = 5;

  /* Third row. */
  A->coeff[5] = 3.0;
  A->col[5] = 1;
  A->coeff[6] = 2.0;
  A->col[6] = 2;
  A->coeff[7] = 4.0;
  A->col[7] = 3;
  A->row[3] = 8;

  /* Print to stdout. */
  csr_ascii(A);

  /* Visualize sparsity pattern with squares */
  /* csr_spy(A, "MatrixWithCircles.eps", 'c'); */
  csr_spy(A, "MatrixWithSquares.eps", 's');

  /* Free memory. */
  csr_free(A);

  return 0;
}

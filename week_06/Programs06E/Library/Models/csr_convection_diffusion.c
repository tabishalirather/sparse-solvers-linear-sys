 /*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/* Some example coefficients. */
void b_func1(int n, int d, int i, double *b) {
  int n1 = n - 1;
  double h = 1.0 / n;
  int iTmp = i;
  int k;
  for (k = 0; k < d; ++k) {
    double xk = ((iTmp % n1) + 1) * h;
    b[k] = 0.5 * (xk + 1.0);
    i /= n1;
  }
}

/* Any-dimensional case where the convection coefficients are constant and
 * the same for each dimension. */
static double b_const;
static void b_const_func(int n, int d, int i, double *b_val) {
  (void) n;
  (void) d;
  (void) i;
  int k;
  for (k = 0; k < d; ++k) {
    b_val[k] = b_const;
  }
}

csr *csr_const_convection_diffusion(int n, int d, double epsilon, double b_const_local) {
  b_const = b_const_local;
  return csr_convection_diffusion(n, d, epsilon, b_const_func);
}

/* Any-dimensional case. */
csr *csr_convection_diffusion(int n, int d, double epsilon, void (*b_func)(int, int, int, double *)) {
  int n1, nd, nnz;
  int i, j;
  int n1pow;
  int k, p;
  csr *A;
  assert(n > 0);
  assert(d > 0);

  /* Number of interior points in one dimension. */
  n1 = n - 1;

  /* Number of interior points in d dimensions. */
  nd = 1;
  for (k = 0; k < d; ++k) {
    nd *= n1;
  }

  /* Compute the number of non-zeros nnz. */
  nnz = laplace_nnz(n, d);

  /* Storage for convection coefficients. */
  double *bval = malloc(d * sizeof(double));

  /* Compute entries and allocate storage. */
  A = csr_alloc(nd, nd, nnz);
  A->row[0] = 0;
  p = 0;
  for (i = 0; i < nd; ++i) {
    /* Compute convection coefficients b = (b_1, ...., b_d) centered at i = (i_1, ..., i_d). */
    b_func(n, d, i, bval);

    /* Left of the diagonal. */
    n1pow = nd / n1;
    for (k = 0; k < d; ++k) {
      j = i - n1pow;
      if ((i / n1pow) % n1 > 0) {
        A->coeff[p] = -epsilon - 0.5 * bval[k] * n;
        A->col[p] = j;
        ++p;
      }
      n1pow /= n1;
    }
    /* On the diagonal. */
    A->coeff[p] = 2 * d * epsilon;
    A->col[p] = i;
    ++p;
    /* Right of the diagonal. */
    n1pow = 1;
    for (k = 0; k < d; ++k) {
      j = i + n1pow;
      if ((i / n1pow) % n1 < n1 - 1) {
        A->coeff[p] = -epsilon + 0.5 * bval[k] * n;
        A->col[p] = j;
        ++p;
      }
      n1pow *= n1;
    }
    A->row[i + 1] = p;
  }
  assert(p == nnz);
  assert(p == A->row[nd]);

  free(bval);

  return A;
}

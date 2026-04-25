/*! \file
 *  
 *  \author Sabine Le Borne
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "itsolver.h"

csr *csr_laplace1d(int n) {
  printf("csr_laplace1d called\n");
  int n1, nnz;
  int i, p;
  double adiag, aoff;
  csr *A;
  assert(n > 0);
  n1 = n - 1;
  /* Compute the number of non-zeros nnz. */
  switch (n1) {
    case 0:
      nnz = 0;
      break;
    case 1:
      nnz = 1;
      break;
    default:
      nnz = 3 * (n1 - 2) + 4;
  }
  adiag = 2.0;
  aoff = -1.0;
  A = csr_alloc(n1, n1, nnz);
  int temp_row[] = {0, 2, 5, 8};
  A->row[0] = 0;
  for (i = 0; i < n1; ++i) {
    // printf("i: %d\n", i);
    // A->row[i] = temp_row[i];


    /* ******************************************************************** */
    /*                                                                      */
    /*     TODO --- Exercise (Sheet 3, Problem 3)                         */
    /*                                                                      */
    /* ******************************************************************** */
    p = A->row[i];
    int j = i-1;
    if (j>=0) {
      A->coeff[p] = aoff;
      A->col[p] = j;
      p++;
    }

    A->coeff[p] = adiag;
    A->col[p] = i;

    j = i+1;
    if (j<n1) {
      A->coeff[p] = aoff;
      A->col[p] = j;
      p++;
    }
    A->row[i+1]=p;

  }
  assert(p == nnz);
  assert(p == A->row[n1]);
  return A;
}

csr *csr_laplace(int n, int d) {
  int n1, nd, nnz;
  int i, j;
  int n1pow;
  int k, p;
  double adiag, aoff;
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

  /* Compute entries and allocate storage. */
  adiag = 2.0 * d;
  aoff = -1.0;
  A = csr_alloc(nd, nd, nnz);
  A->row[0] = 0;
  for (i = 0; i < nd; ++i) {
    p = A->row[i];
    /* Left of the diagonal. */
    n1pow = nd / n1;
    for (k = 0; k < d; ++k) {
      j = i - n1pow;
      if ((i / n1pow) % n1 > 0) {
        A->coeff[p] = aoff;
        A->col[p] = j;
        ++p;
      }
      n1pow /= n1;
    }
    /* On the diagonal. */
    A->coeff[p] = adiag;
    A->col[p] = i;
    ++p;
    /* Right of the diagonal. */
    n1pow = 1;
    for (k = 0; k < d; ++k) {
      j = i + n1pow;
      if ((i / n1pow) % n1 < n1 - 1) {
        A->coeff[p] = aoff;
        A->col[p] = j;
        ++p;
      }
      n1pow *= n1;
    }
    A->row[i + 1] = p;
  }
  assert(p == nnz);
  assert(p == A->row[nd]);
  return A;
}


void laplace_ev(int n, int d, int *i, double *v) {
  int n1, nd;
  int j, k;

  n1 = n - 1;
  nd = 1;
  for (k = 0; k < d; ++k) {
    nd *= n1;
  }

  for (j = 0; j < nd; ++j) {
    int jTmp = j;
    v[j] = 1.0;
    for (k = 0; k < d; ++k) {
      int jk = jTmp % n1;
      jTmp /= n1;
      v[j] *= sin((i[k] + 1) * (jk + 1) * M_PI / n);
    }
  }
}

void laplace_ev_min(int n, int d, double *v) {
  int n1, nd;
  int j, k;

  n1 = n - 1;
  nd = 1;
  for (k = 0; k < d; ++k) {
    nd *= n1;
  }

  for (j = 0; j < nd; ++j) {
    int jTmp = j;
    v[j] = 1.0;
    for (k = 0; k < d; ++k) {
      int jk = jTmp % n1;
      jTmp /= n1;
      v[j] *= sin((jk + 1) * M_PI / n);
    }
  }
}

void laplace_ev_max(int n, int d, double *v) {
  int n1, nd;
  int j, k;

  n1 = n - 1;
  nd = 1;
  for (k = 0; k < d; ++k) {
    nd *= n1;
  }

  for (j = 0; j < nd; ++j) {
    int jTmp = j;
    v[j] = 1.0;
    for (k = 0; k < d; ++k) {
      int jk = jTmp % n1;
      jTmp /= n1;
      v[j] *= sin(n1 * (jk + 1) * M_PI / n);
    }
  }
}

double laplace_ew(int n, int d, int *i) {
  int k;
  double lambda = 0.0;
  for (k = 0; k < d; ++k) {
    double s = sin((i[k] + 1) * M_PI / (2 * n));
    lambda += 4.0 * s * s;
  }
  return lambda;
}

double laplace_ew_min(int n, int d) {
  double s = sin(M_PI / (2 * n));
  double lambda = 4.0 * d * s * s;
  return lambda;
}

double laplace_ew_max(int n, int d) {
  double s = sin((n - 1) * M_PI / (2 * n));
  double lambda = 4.0 * d * s * s;
  return lambda;
}

double laplace_cond(int n, int d) {
  double lambda_min = laplace_ew_min(n, d);
  double lambda_max = laplace_ew_max(n, d);
  double kappa = lambda_max / lambda_min;
  return kappa;
}

int laplace_nnz(int n, int d) {
  int i;
  int n1, nd, nnz;
  assert(n > 0);
  assert(d > 0);
  n1 = n - 1;
  nd = 1;
  for (i = 0; i < d; ++i) {
    nd *= n1;
  }
  if (n1 == 0) {
    nnz = 0;
  } else {
    nnz = (2 * d + 1) * nd - 2 * d * (nd / n1);
  }
  assert(nnz > 0);
  return nnz;
}

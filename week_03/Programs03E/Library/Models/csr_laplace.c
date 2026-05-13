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
  adiag = 3.0;
  // aoff = -1.0;
  A = csr_alloc(n1, n1, nnz);
  A->row[0] = 0;
  A->col[0] = 0;
  int nr = n;

  // for (i = 0; i < n1; ++i) {

  /* ******************************************************************** */
  /*                                                                      */
  /*     TODO --- Exercise (Sheet 3, Problem 3)                         */
  /*                                                                      */
  /* ******************************************************************** */
  //The goal is to create tri-diagonal matrices with 2 on diagonal and -1 on off diagonal but in CSR format.
  /* For example, for a 4*4 matrix of this type
   * [
   * 2 -1  0  0
   * -1 2 -1  0
   * 0 -1  2 -1
   * 0  0 -1  2
   */
  //Let's fill diagonal elements first. In csr format,
  //Col will be like  0, 1, 2, 3, ....n and row will be 0,1,2,3,... with coeff = 2,2,2,2,...
  // A->row[0] = 0;

  // int col[10] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
  // double coeff[10] = {2, 3, 1, 2, 3, 1, 2, 3, 1, 2};
  double main = adiag;
  puts("populating the tridiagonal matrix");
  for (int i = 1; i <= nr; i++) {
    if (i == 1 ^ i == nr)
      A->row[i] = A->row[i - 1] + 2;
    else {
      A->row[i] = A->row[i - 1] + 3;
    }
  }

  // for (int j = 0; j <= nr; j++) {
  //   printf("%d ", A->row[j]);
  // }
  puts("row matrix setup");
  for (int i = 0; i <= A->nr; i++) {
    // printf("running inside the for loop");
    for (int j = A->row[i]; j < A->row[i + 1]; j++) {
      double lower = -1;
      double upper = -1;
      // printf("running inside the j for loop");
      if (i == 0) {
        //this means first row, i.e. fill main and upper diagonal.
        if (j == A->row[i]) {
          //then we are filling the first element
          A->coeff[j] = main;
          // The general pattern for filling columns in i-1, i, i+1
          A->col[j] = i; // index of first nz col
        } else if (j == A->row[i] + 1) {
          //then we are filling the second element
          A->coeff[j] = upper;
          A->col[j] = i + 1; //index of second nz col.
        }
      } else if (i + 1 == A->nr) {
        //this means last row, i.e. fill lower and main diagonal.
        if (j == A->row[i]) {
          //then we are filling the first element
          A->coeff[j] = lower;
          A->col[j] = i - 1;
        } else if (j == A->row[i] + 1) {
          //then we are filling the second element
          A->coeff[j] = main;
          A->col[j] = i;
        }
      } else {
        //fill lower, then main and then upper diagonal.
        if (j == A->row[i]) {
          //then we are filling the first element
          A->coeff[j] = lower;
          A->col[j] = i - 1;
        } else if (j == A->row[i] + 1) {
          //then we are filling the second element
          A->coeff[j] = main;
          A->col[j] = i;
        } else if (j == A->row[i] + 2) {
          //then we are filling the second element
          A->coeff[j] = upper;
          A->col[j] = i + 1;
        }
      }
    }
  }


  printf("%d \n", n1);
  puts("Matrix A is in week3 csr laplace: ");
  for (int i = 0; i < n1; i++) {
    for (int j = A->row[i]; j < A->row[i + 1]; j++) {
      printf("%f ", A->coeff[j]);
    }
    printf("\n");
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

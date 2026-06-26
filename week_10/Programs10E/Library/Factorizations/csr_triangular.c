/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */


#include "itsolver.h"

csr *csr_extract_unit_lower(csr *LU) {
  csr *L = malloc(sizeof(csr));
  assert(LU->nr == LU->nc);
  int n = LU->nr;
  L->nr = n;
  L->nc = n;
  L->row = malloc((n + 1) * sizeof(int));
  L->row[0] = 0;
  int i, p;
  /* Count coefficients. */
  for (i = 0; i < n; ++i) {
    /* Plus one because of the additional diagonal entry in each row. */
    L->row[i + 1] = L->row[i] + 1;
    for (p = LU->row[i]; p < LU->row[i + 1]; ++p) {
      if (LU->col[p] < i) {
        L->row[i + 1]++;
      }
    }
  }
  L->nz = L->row[n];
  L->col = malloc(L->nz * sizeof(int));
  L->coeff = malloc(L->nz * sizeof(double));
  /* Copy coefficients and add an unit diagonal. */
  int q = 0;
  for (i = 0; i < n; ++i) {
    for (p = LU->row[i]; p < LU->row[i + 1]; ++p) {
      if (LU->col[p] < i) {
        L->col[q] = LU->col[p];
        L->coeff[q] = LU->coeff[p];
        q++;
      }
    }
    L->col[q] = i;
    L->coeff[q] = 1.0;
    q++;
  }
  assert(q == L->nz);
  return L;
}

csr *csr_extract_upper(csr *LU) {
  csr *U = malloc(sizeof(csr));
  assert(LU->nr == LU->nc);
  int n = LU->nr;
  U->nr = n;
  U->nc = n;
  U->row = malloc((n + 1) * sizeof(int));
  U->row[0] = 0;
  int i, p;
  /* Count coefficients. */
  for (i = 0; i < n; ++i) {
    U->row[i + 1] = U->row[i];
    for (p = LU->row[i]; p < LU->row[i + 1]; ++p) {
      if (LU->col[p] >= i) {
        U->row[i + 1]++;
      }
    }
  }
  U->nz = U->row[n];
  U->col = malloc(U->nz * sizeof(int));
  U->coeff = malloc(U->nz * sizeof(double));
  /* Copy coefficients. */
  int q = 0;
  for (i = 0; i < n; ++i) {
    for (p = LU->row[i]; p < LU->row[i + 1]; ++p) {
      if (LU->col[p] >= i) {
        U->col[q] = LU->col[p];
        U->coeff[q] = LU->coeff[p];
        q++;
      }
    }
  }
  assert(q == U->nz);
  return U;
}

void csr_lowersolve(csr *L, double *x) {
  int i, p;
  assert(L->nr == L->nc);

  /*****************************************************************/
  /*                                                               */
  /*   TODO: Exercise sheet 9, problem 2                          */
  /*                                                               */
  /*****************************************************************/
  
}

void csr_uppersolve(csr *U, double *x) {
  int i, p;
  assert(U->nr == U->nc);
  /*****************************************************************/
  /*                                                               */
  /*   TODO: Exercise sheet 9, problem 2                          */
  /*                                                               */
  /*****************************************************************/
}

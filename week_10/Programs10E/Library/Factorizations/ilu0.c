/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void ilu0_inplace_csr_row(csr *A, int i)
{
  assert(i >= 0 && i < A->nr);
  /* For all existing A_{i, j} where j < i. */
  int pij;
  for (pij = A->row[i]; pij < A->row[i + 1] && A->col[pij] < i; ++pij) {
    int j = A->col[pij];
    /* Find A_{j,j}, which we assume/ has to exist. */
    int pjj;
    for (pjj = A->row[j]; pjj < A->row[j + 1]; ++pjj) {
      if (A->col[pjj] == j) {
        assert(A->coeff[pjj] != 0.0);
        A->coeff[pij] /= A->coeff[pjj];
        break;
      }
    }
    /* Assert that we have found A_{j,j}. */
    assert(pjj < A->row[j + 1]);
    /* For all existing A_{i, k} such that k > j. */
    int pik;
    for (pik = pij + 1; pik < A->row[i + 1]; ++pik) {
      int k = A->col[pik];
      /* Try to find A_{j, k}, which might or might not exist. */
      int pjk;
      for (pjk = A->row[j]; pjk < A->row[j + 1]; ++pjk) {
        if (A->col[pjk] == k) {
          A->coeff[pik] -= A->coeff[pij] * A->coeff[pjk];
          break;
        }
      }
    }
  }
  /* We are done with row i, and now A_{i, i} equals U_{i, i}.
     Assert that this is non-zero since we will divide by this
     value in the next step. */
  for (pij = A->row[i]; pij < A->row[i + 1]; ++pij) {
    if (A->col[pij] == i) {
      assert(A->coeff[pij] != 0.0);
      break;
    }
  }
}

void ilu0_inplace_csr(csr *A) {
  /* For all rows i = 0, ..., n - 1. */
  int i;
  for (i = 0; i < A->nr; ++i) {
    ilu0_inplace_csr_row(A, i);
  }
}

csr *ilu0_csr(csr *A) {
  assert(A->nr == A->nc);
  csr *LU = csr_copy(A);
  ilu0_inplace_csr(LU);
  return LU;
}

double ilu0_check_csr(csr *A, csr *LU) {
  int n = A->nr;
  csr *L = csr_extract_unit_lower(LU);
  csr *U = csr_extract_upper(LU);
  csr *LtimesU = csr_multiply(L, U);
  csr_free(L);
  csr_free(U);
  double error = 0.0;
  int i, p, q;
  for (i = 0; i < n; ++i) {
    q = LtimesU->row[i];
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      while (q < LtimesU->row[i + 1] && A->col[p] > LtimesU->col[q]) {
	++q;
      }
      if (q < LtimesU->row[i + 1] && A->col[p] == LtimesU->col[q]) {	
	double tmp = A->coeff[p] - LtimesU->coeff[q];
	error += tmp * tmp;
      }
    }
  }
  csr_free(LtimesU);
  error = sqrt(error);
  return error;
}


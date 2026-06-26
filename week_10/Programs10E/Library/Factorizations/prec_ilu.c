/*! \file
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

static void prec_ilu0_func(void *data, double *r, double *z) {
  csr *LU = data;
  int n = LU->nr;
  vec_copy(n, r, z);
  csr_lowersolve(LU, z);
  csr_uppersolve(LU, z);
}

static void prec_ilu0_delete(prec_t *prec) {
  csr *LU = prec->data;
  csr_free(LU);
  free(prec);
}

prec_t *prec_ilu0_new(csr *A) {
  prec_t *prec = malloc(sizeof(prec_t));
  prec->func = &prec_ilu0_func;
  prec->delete = &prec_ilu0_delete;
  csr *LU = csr_copy(A);
  ilu0_inplace_csr(LU);
  prec->data = LU;
  return prec;
}

static void prec_ilup_func(void *data, double *r, double *z) {
  csr *LU = data;
  int n = LU->nr;
  vec_copy(n, r, z);
  csr_lowersolve(LU, z);
  csr_uppersolve(LU, z);
}

static void prec_ilup_delete(prec_t *prec) {
  csr *LU = prec->data;
  csr_free(LU);
  free(prec);
}

prec_t *prec_ilup_new(csr *A, int p) {
  assert(p >= 0);
  prec_t *prec = malloc(sizeof(prec_t));
  prec->func = &prec_ilup_func;
  prec->delete = &prec_ilup_delete;
  csr *LU = ilup_csr(A, p);
  prec->data = LU;
  return prec;
}

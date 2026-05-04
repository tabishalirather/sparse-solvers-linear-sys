
#include "itsolver.h"

void prec_apply(prec_t *prec, double *r, double *z) {
  prec->func(prec->data, r, z);
}

void prec_delete(prec_t *prec) {
  prec->delete(prec);
}

static void prec_identity_func(void *data, double *r, double *z) {
  csr *A = data;
  assert(A->nc == A->nr);
  int n = A->nc;
  int i;
  for (i = 0; i < n; ++i) {
    z[i] = r[i];
  }
}

static void prec_identity_delete(prec_t *prec) {
  prec->func = NULL;
  prec->delete = NULL;
  prec->data = NULL;
  free(prec);
}

prec_t *prec_identity_new(csr *A) {
  prec_t *prec = malloc(sizeof(prec_t));
  prec->func = &prec_identity_func;
  prec->delete = &prec_identity_delete;
  prec->data = A;
  return prec;
}

static void prec_jacobi_func(void *data, double *r, double *z) {
  int n;
  csr *A = data;
  int i, p;
  assert(A->nc == A->nr);
  n = A->nc;
  for (i = 0; i < n; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] == i) {
	z[i] = r[i] / A->coeff[p];
	continue;
      }
    }
  }
}

static void prec_jacobi_delete(prec_t *prec) {
  prec->func = NULL;
  prec->delete = NULL;
  prec->data = NULL;
  free(prec);
}

prec_t *prec_jacobi_new(csr *A) {
  prec_t *prec = malloc(sizeof(prec_t));
  prec->func = &prec_jacobi_func;
  prec->delete = &prec_jacobi_delete;
  prec->data = A;
  return prec;
}

static void prec_gs_func(void *data, double *r, double *z) {
  int n;
  csr *A = data;
  int i, p;
  assert(A->nc == A->nr);
  n = A->nc;

  /* Copy r to z. */
  for (i = 0; i < n; ++i) {
    z[i] = r[i];
  }

  /* Forward substitution on z. */
  for (i = 0; i < n; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] < i) {
	z[i] -= A->coeff[p] * z[A->col[p]];
      }
    }
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] == i) {
	z[i] /= A->coeff[p];
	continue;
      }
    }
  }
}

static void prec_gs_delete(prec_t *prec) {
  prec->func = NULL;
  prec->delete = NULL;
  prec->data = NULL;
  free(prec);
}

prec_t *prec_gs_new(csr *A) {
  prec_t *prec = malloc(sizeof(prec_t));
  prec->func = &prec_gs_func;
  prec->delete = &prec_gs_delete;
  prec->data = A;
  return prec;
}

static void prec_sgs_func(void *data, double *r, double *z) {
  int n;
  csr *A = data;
  int i, p;
  assert(A->nc == A->nr);
  n = A->nc;

  /* Copy r to z. */
  for (i = 0; i < n; ++i) {
    z[i] = r[i];
  }

  /* Forward substitution on z. */
  for (i = 0; i < n; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] < i) {
	z[i] -= A->coeff[p] * z[A->col[p]];
      }
    }
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] == i) {
	z[i] /= A->coeff[p];
	continue;
      }
    }
  }

  /* Multiply z by diagonal. */
  for (i = 0; i < n; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] == i) {
	z[i] *= A->coeff[p];
	continue;
      }
    }
  }

  /* Backward substitution on z. */
  for (i = n - 1; i >= 0; --i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] > i) {
	z[i] -= A->coeff[p] * z[A->col[p]];
      }
    }
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] == i) {
	z[i] /= A->coeff[p];
	continue;
      }
    }
  }
}

static void prec_sgs_delete(prec_t *prec) {
  prec->func = NULL;
  prec->delete = NULL;
  prec->data = NULL;
  free(prec);
}

prec_t *prec_sgs_new(csr *A) {
  prec_t *prec = malloc(sizeof(prec_t));
  prec->func = &prec_sgs_func;
  prec->delete = &prec_sgs_delete;
  prec->data = A;
  return prec;
}

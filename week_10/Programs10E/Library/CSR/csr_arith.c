
#include "itsolver.h"

double csr_xtAy(csr *A, double *x, double *y) {
  int i, p;
  int nr = A->nr;
  double dyad = 0.0;
  for (i = 0; i < nr; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      dyad += x[i] * A->coeff[p] * y[A->col[p]];
    }
  }
  return dyad;
}

void csr_mv_transposed(csr *A, double alpha, double *x, double beta, double *y) {
  int i, p;
  int nr = A->nr;
  int nc = A->nc;
  for (i = 0; i < nc; ++i) {
    y[i] *= beta;
  }
  for (i = 0; i < nr; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      y[A->col[p]] += alpha * A->coeff[p] * x[i];
    }
  }
}

void csr_residual(csr *A, double *b, double *x, double *r) {
  int i, p;
  for (i = 0; i < A->nr; ++i) {
    r[i] = b[i];
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      r[i] -= A->coeff[p] * x[A->col[p]];
    }
  }
}

double csr_residual_norm(csr *A, double *b, double *x) {
  int i, p;
  double r_at; /* Current entry of residual vector. */
  double r_sq = 0.0;
  for (i = 0; i < A->nr; ++i) {
    r_at = b[i];
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      r_at -= A->coeff[p] * x[A->col[p]];
    }
    r_sq += r_at * r_at;
  }
  return sqrt(r_sq);
}

void csr_residual_transposed(csr *A, double *b, double *x, double *r) {
  int i, n, p;
  assert(A->nr == A->nc);
  n = A->nr;
  for (i = 0; i < n; ++i) {
    r[i] = b[i];
  }
  for (i = 0; i < n; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      r[A->col[p]] -= A->coeff[p] * x[i];
    }
  }
}

csr *csr_multiply(csr *A, csr *B) {
  assert(A->nc == B->nr);
  csr *C = malloc(sizeof(csr));
  C->nr = A->nr;
  C->nc = B->nc;
  int i, p;
  int c = 0;
  for (i = 0; i < A->nr; ++i) {
    int ci = A->row[i + 1] - A->row[i];
    if (c < ci) {
      c = ci;
    }
  }
  int *pick = malloc(c * sizeof(int));
  /* Count the number of non-zeros. */
  C->row = malloc((A->nr + 1) * sizeof(int));
  C->row[0] = 0;
  for (i = 0; i < A->nr; ++i) {
    C->row[i + 1] = C->row[i];
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      pick[p - A->row[i]] = B->row[A->col[p]];
    }
    while (1) {
      int jmin = INT_MAX;
      for (p = A->row[i]; p < A->row[i + 1]; ++p) {
        if (pick[p - A->row[i]] < B->row[A->col[p] + 1]) {
          if (jmin > B->col[pick[p - A->row[i]]]) {
            jmin = B->col[pick[p - A->row[i]]];
          }
        }
      }
      if (jmin == INT_MAX) {
        break;
      }
      C->row[i + 1]++;
      for (p = A->row[i]; p < A->row[i + 1]; ++p) {
        if (pick[p - A->row[i]] < B->row[A->col[p] + 1]) {
          if (jmin == B->col[pick[p - A->row[i]]]) {
            pick[p - A->row[i]]++;
          }
        }
      }
    }
  }
  C->nz = C->row[A->nr];
  /* Compute the coefficients. */
  C->col = malloc(C->nz * sizeof(int));
  C->coeff = malloc(C->nz * sizeof(double));
  int q = 0;
  for (i = 0; i < A->nr; ++i) {
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      pick[p - A->row[i]] = B->row[A->col[p]];
    }
    while (1) {
      int jmin = INT_MAX;
      for (p = A->row[i]; p < A->row[i + 1]; ++p) {
        if (pick[p - A->row[i]] < B->row[A->col[p] + 1]) {
          if (jmin > B->col[pick[p - A->row[i]]]) {
            jmin = B->col[pick[p - A->row[i]]];
          }
        }
      }
      if (jmin == INT_MAX) {
        break;
      }
      C->col[q] = jmin;
      C->coeff[q] = 0.0;
      for (p = A->row[i]; p < A->row[i + 1]; ++p) {
        if (pick[p - A->row[i]] < B->row[A->col[p] + 1]) {
          if (jmin == B->col[pick[p - A->row[i]]]) {
            C->coeff[q] += A->coeff[p] * B->coeff[pick[p - A->row[i]]++];
          }
        }
      }
      q++;
    }
  }
  assert(q == C->nz);
  free(pick);
  return C;
}

void csr_diag_comp(csr *A, double alpha) {
  assert(A->nr == A->nc);
  int n = A->nr;
  int i = 0;
  for (; i < n; ++i) {
    int p = A->row[i];
    for (; p < A->row[i + 1] && A->col[p] != i; ++p) {
    }
    /* Assert that a diagonal element has been found. */
    assert(p < A->row[i + 1]);
    /* Add to diagonal element. */
    A->coeff[p] += alpha;
  }
}

csr *csr_add(csr *a, csr *b) {
  int nr, nc, nz;
  int i, pa, pb, pc;
  csr *c;
  assert(a->nr == b->nr);
  assert(a->nc == b->nc);
  nr = a->nr;
  nc = a->nc;
  nz = 0;
  for (i = 0; i < nr; ++i) {
    pa = a->row[i];
    pb = b->row[i];
    while (1) {
      if (pa < a->row[i + 1] && pb < b->row[i + 1]) {
        if (a->col[pa] == b->col[pb]) {
          ++pa;
          ++pb;
          ++nz;
          continue;
        }
        if (a->col[pa] < b->col[pb]) {
          ++pa;
          ++nz;
          continue;
        }
        if (a->col[pa] > b->col[pb]) {
          ++pb;
          ++nz;
          continue;
        }
        assert(0);
      }
      if (pa == a->row[i + 1] && pb < b->row[i + 1]) {
        ++pb;
        ++nz;
        continue;
      }
      if (pa < a->row[i + 1] && pb == b->row[i + 1]) {
        ++pa;
        ++nz;
        continue;
      }
      assert(pa == a->row[i + 1] && pb == b->row[i + 1]);
      break;
    }
  }
  c = csr_alloc(nr, nc, nz);
  pc = 0;
  for (i = 0; i < nr; ++i) {
    pa = a->row[i];
    pb = b->row[i];
    c->row[i] = pc;
    while (1) {
      if (pa < a->row[i + 1] && pb < b->row[i + 1]) {
        if (a->col[pa] == b->col[pb]) {
          c->col[pc] = a->col[pa];
          c->coeff[pc] = a->coeff[pa] + b->coeff[pb];
          ++pa;
          ++pb;
          ++pc;
          continue;
        }
        if (a->col[pa] < b->col[pb]) {
          c->col[pc] = a->col[pa];
          c->coeff[pc] = a->coeff[pa];
          ++pa;
          ++pc;
          continue;
        }
        if (a->col[pa] > b->col[pb]) {
          c->col[pc] = b->col[pb];
          c->coeff[pc] = b->coeff[pb];
          ++pb;
          ++pc;
          continue;
        }
        assert(0);
      }
      if (pa == a->row[i + 1] && pb < b->row[i + 1]) {
        c->col[pc] = b->col[pb];
        c->coeff[pc] = b->coeff[pb];
        ++pb;
        ++pc;
        continue;
      }
      if (pa < a->row[i + 1] && pb == b->row[i + 1]) {
        c->col[pc] = a->col[pa];
        c->coeff[pc] = a->coeff[pa];
        ++pa;
        ++pc;
        continue;
      }
      assert(pa == a->row[i + 1] && pb == b->row[i + 1]);
      break;
    }
  }
  c->row[nr] = pc;
  return c;
}

csr *csr_substract(csr *a, csr *b) {
  int nr, nc, nz;
  int i, pa, pb, pc;
  csr *c;
  assert(a->nr == b->nr);
  assert(a->nc == b->nc);
  nr = a->nr;
  nc = a->nc;
  nz = 0;
  for (i = 0; i < nr; ++i) {
    pa = a->row[i];
    pb = b->row[i];
    while (1) {
      if (pa < a->row[i + 1] && pb < b->row[i + 1]) {
        if (a->col[pa] == b->col[pb]) {
          ++pa;
          ++pb;
          ++nz;
          continue;
        }
        if (a->col[pa] < b->col[pb]) {
          ++pa;
          ++nz;
          continue;
        }
        if (a->col[pa] > b->col[pb]) {
          ++pb;
          ++nz;
          continue;
        }
        assert(0);
      }
      if (pa == a->row[i + 1] && pb < b->row[i + 1]) {
        ++pb;
        ++nz;
        continue;
      }
      if (pa < a->row[i + 1] && pb == b->row[i + 1]) {
        ++pa;
        ++nz;
        continue;
      }
      assert(pa == a->row[i + 1] && pb == b->row[i + 1]);
      break;
    }
  }
  c = csr_alloc(nr, nc, nz);
  pc = 0;
  for (i = 0; i < nr; ++i) {
    pa = a->row[i];
    pb = b->row[i];
    c->row[i] = pc;
    while (1) {
      if (pa < a->row[i + 1] && pb < b->row[i + 1]) {
        if (a->col[pa] == b->col[pb]) {
          c->col[pc] = a->col[pa];
          c->coeff[pc] = a->coeff[pa] - b->coeff[pb];
          ++pa;
          ++pb;
          ++pc;
          continue;
        }
        if (a->col[pa] < b->col[pb]) {
          c->col[pc] = a->col[pa];
          c->coeff[pc] = a->coeff[pa];
          ++pa;
          ++pc;
          continue;
        }
        if (a->col[pa] > b->col[pb]) {
          c->col[pc] = b->col[pb];
          c->coeff[pc] = -b->coeff[pb];
          ++pb;
          ++pc;
          continue;
        }
        assert(0);
      }
      if (pa == a->row[i + 1] && pb < b->row[i + 1]) {
        c->col[pc] = b->col[pb];
        c->coeff[pc] = b->coeff[pb];
        ++pb;
        ++pc;
        continue;
      }
      if (pa < a->row[i + 1] && pb == b->row[i + 1]) {
        c->col[pc] = a->col[pa];
        c->coeff[pc] = a->coeff[pa];
        ++pa;
        ++pc;
        continue;
      }
      assert(pa == a->row[i + 1] && pb == b->row[i + 1]);
      break;
    }
  }
  c->row[nr] = pc;
  return c;
}

void csr_scale(csr *a, double alpha) {
  vec_scale(a->nz, alpha, a->coeff);
}

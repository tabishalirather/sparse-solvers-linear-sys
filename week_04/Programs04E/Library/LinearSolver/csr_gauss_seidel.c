
#include "itsolver.h"

void csr_gauss_seidel(
  csr *A, double *b, double *x,
  double tol, int stopcrit, int maxit,
  double *esterr, int *iter, double *errvec)
{
  int i, j, p;
  int n;
  double rho = 0.0;
  double *r;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate residual vector. */
  r = (double *) malloc(n * sizeof(double));

  /* Compute start residual. */
  vec_copy(n, b, r);
  csr_mv(A, -1.0, x, 1.0, r);
  if (stopcrit == ITSOLVER_RELRES) {
    rho = vec_norm(n, r);
    errvec[0] = rho;
  }

  /* Start the iteration. */
  j = 0;
  while (j < maxit) {
    if (stopcrit == ITSOLVER_RELRES && errvec[j] < tol * rho) {
      break;
    }
    /* Compute new iterate. */
    for (i = 0; i < n; ++i) {
      for (p = A->row[i]; p < A->row[i + 1]; ++p) {
        if (A->col[p] < i) {
          r[i] -= A->coeff[p] * r[A->col[p]];
        }
      }
      for (p = A->row[i]; p < A->row[i + 1]; ++p) {
        if (A->col[p] == i) {
          r[i] /= A->coeff[p];
          x[i] += r[i];
          /* Skip to the next row. */
          break;
        }
      }
    }

    /* Compute new residual. */
    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    if (stopcrit == ITSOLVER_RELRES) {
      errvec[j + 1] = vec_norm(n, r);
    }

    /* Proceed to the next step. */
    j++;
  }

  /* Output iter and esterr. */
  *iter = j;
  if (stopcrit == ITSOLVER_RELRES) {
    *esterr = errvec[j];
  }

  /* Free memory. */
  free(r);
}

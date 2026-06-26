 /*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

static void csr_richardson_preconditioned(
  csr *A, prec_t *prec, double *b, double *x,
  double tol, int stopcrit, int maxit,
  double *esterr, int *iter, double *errvec)
{
  int i, j;
  int n;
  double rho = 0.0;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate residual and preconditioned residual. */
  double *r = (double *) malloc(n * sizeof(double));
  double *z = (double *) malloc(n * sizeof(double));

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

    /* Apply Preconditioner. */
    prec_apply(prec, r, z);

    /* Add preconditioned residual to iterate. */
    for (i = 0; i < n; ++i) {
      x[i] += z[i];
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
  free(z);
}

void csr_sgs(csr *A, double *b, double *x,
             double tol, int stopcrit, int maxit,
             double *esterr, int *iter, double *errvec)
{
  prec_t *prec = prec_sgs_new(A);
  csr_richardson_preconditioned(A, prec, b, x, tol, stopcrit, maxit, esterr, iter, errvec);
  prec_delete(prec);
}

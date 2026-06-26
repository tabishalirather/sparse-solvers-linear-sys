/*! \file
 *  \brief Biconjugate gradients (BiCG).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void csr_bicg(csr *A, double *b, double *x,
        double tol, int stopcrit, int maxit,
        double *esterr, int *iter, double *errvec) {
  int k;
  int n;
  double done = 1.0;
  double rho = 0.0;
  double alpha, beta;
  double rs0, rs1;
  double *r, *s;
  double *p, *q;
  double *Ap, *qA;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate vectors. */
  r = (double *) malloc(n * sizeof(double));
  s = (double *) malloc(n * sizeof(double));
  p = (double *) malloc(n * sizeof(double));
  q = (double *) malloc(n * sizeof(double));
  Ap = (double *) calloc(n, sizeof(double));
  qA = (double *) calloc(n, sizeof(double));

  /* Right and left residuals r and s. */
  csr_residual(A, b, x, r);
  csr_residual_transposed(A, b, x, s);

  /* Right and left direction vectors p and q. */
  vec_copy(n, r, p);
  vec_copy(n, s, q);

  /* Compute residual norm if requested. */
  if (stopcrit == ITSOLVER_RELRES) {
    rho = vec_norm(n, r);
    errvec[0] = rho;
  }

  /* Scalar products between r and s. */
  rs0 = 0.0;
  rs1 = vec_dot(n, r, s);

  /* Start the iteration. */
  k = 0;
  while (k < maxit) {
    /* Check whether the desired accuracy has been reached. */
    if (stopcrit == ITSOLVER_RELRES && errvec[k] < tol * rho) {
  break;
    }

    /* One step of BiCG. */
    rs0 = rs1;
    csr_mv(A, 1.0, p, 0.0, Ap);
    csr_mv_transposed(A, 1.0, q, 0.0, qA);
    alpha = rs0 / vec_dot(n, q, Ap);
    vec_axpy(n, alpha, p, x);
    alpha *= -1;
    vec_axpy(n, alpha, Ap, r);
    vec_axpy(n, alpha, qA, s);
    rs1 = vec_dot(n, r, s);
    beta = rs1 / rs0;
    vec_scale(n, beta, p);
    vec_axpy(n, done, r, p);
    vec_scale(n, beta, q);
    vec_axpy(n, done, s, q);

    /* Compute residual norm if requested. */
    if (stopcrit == ITSOLVER_RELRES) {
      errvec[k + 1] = vec_norm(n, r);
    }

    /* Proceed to the next step. */
    k++;
  }

  /* Output iter and esterr. */
  *iter = k;
  if (stopcrit == ITSOLVER_RELRES) {
    *esterr = errvec[k];
  }

  /* Free memory. */
  free(r); free(s);
  free(p); free(q);
  free(Ap); free(qA);
}


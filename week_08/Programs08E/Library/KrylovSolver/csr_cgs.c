/*! \file
 *  \brief Conjugate gradients squared (CGS).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void csr_cgs(csr *A, double *b, double *x,
       double tol, int stopcrit, int maxit,
       double *esterr, int *iter, double *errvec) {
  int k;
  int n;
  double done = 1.0;
  double rho = 0.0;
  double alpha, malpha, beta;
  double rs0, rs1;
  double *r, *s, *p, *u, *q, *Ap;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate vectors. */
  r = (double *) malloc(n * sizeof(double));
  s = (double *) malloc(n * sizeof(double));
  p = (double *) malloc(n * sizeof(double));
  u = (double *) malloc(n * sizeof(double));
  q = (double *) malloc(n * sizeof(double));
  Ap = (double *) calloc(n, sizeof(double));

  /* Residual r and vectors s, p, u. */
  csr_residual(A, b, x, r);
  csr_residual_transposed(A, b, x, s);
  vec_copy(n, r, p);
  vec_copy(n, r, u);

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

    /* One step of CGS. */
    rs0 = rs1;
    csr_mv(A, 1.0, p, 0.0, Ap);
    alpha = rs0 / vec_dot(n, s, Ap);
    malpha = -alpha;

    /* q = u - alpha * A * p */
    vec_copy(n, u, q);
    vec_axpy(n, malpha, Ap, q);

    /* u = u + q. */
    vec_axpy(n, done, q, u);

    /* x = x + alpha * (u + q) */
    vec_axpy(n, alpha, u, x);

    /* r = r - alpha * A * (u + q) */
    csr_mv(A, 1.0, u, 0.0, Ap);
    vec_axpy(n, malpha, Ap, r);

    /* beta = <r1,s> / <r0,s> */
    rs1 = vec_dot(n, r, s);
    beta = rs1 / rs0;

    /* u = r + beta * q */
    vec_copy(n, r, u);
    vec_axpy(n, beta, q, u);

    /* p = u + beta * (q + beta * p) */
    vec_axpy(n, beta, p, q);
    vec_copy(n, u, p);
    vec_axpy(n, beta, q, p);

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
  free(r);
  free(s);
  free(p);
  free(u);
  free(q);
  free(Ap);
}


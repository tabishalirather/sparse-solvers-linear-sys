/*! \file
 *  \brief Stabilized biconjugate gradients (BiCGStab).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void csr_bicgstab(csr *A, double *b, double *x,
      double tol, int stopcrit, int maxit,
      double *esterr, int *iter, double *errvec) {
  int k;
  int n;
  double rho = 0.0;
  double alpha, malpha, omega, momega, beta;
  double rtilder0, rtilder1;
  double *r, *tilder, *p, *s, *Ap, *As;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate vectors. */
  r = (double *) malloc(n * sizeof(double));
  tilder = (double *) malloc(n * sizeof(double));
  p = (double *) malloc(n * sizeof(double));
  s = (double *) malloc(n * sizeof(double));
  Ap = (double *) calloc(n, sizeof(double));
  As = (double *) calloc(n, sizeof(double));

  /* Residual r and direction vector p. */
  csr_residual(A, b, x, r);
  csr_residual_transposed(A, b, x, tilder);
  /*
  for (i = 0; i < n; ++i) {
    tilder[i] = drand48();
  }
  */
  vec_copy(n, r, p);

  /* Compute residual norm if requested. */
  if (stopcrit == ITSOLVER_RELRES) {
    rho = vec_norm(n, r);
    errvec[0] = rho;
  }

  /* Scalar products between r and s. */
  rtilder0 = 0.0;
  rtilder1 = vec_dot(n, r, tilder);

  /* Start the iteration. */
  k = 0;
  while (k < maxit) {
    /* Check whether the desired accuracy has been reached. */
    if (stopcrit == ITSOLVER_RELRES && errvec[k] < tol * rho) {
  break;
    }

    /* One step of BiCGStab. */
    rtilder0 = rtilder1;
    csr_mv(A, 1.0, p, 0.0, Ap);
    alpha = rtilder0 / vec_dot(n, Ap, tilder);
    malpha = -alpha;

    /* s = r - alpha * A * p */
    vec_copy(n, r, s);
    vec_axpy(n, malpha, Ap, s);

    /* omega = <As,s> / <As,As> */
    csr_mv(A, 1.0, s, 0.0, As);
    omega = vec_dot(n, As, s) / vec_dot(n, As, As);
    momega = -omega;

    /* x = x + alpha * p + omega * s */
    vec_axpy(n, alpha, p, x);
    vec_axpy(n, omega, s, x);

    /* r = s - omega * A * s */
    vec_copy(n, s, r);
    vec_axpy(n, momega, As, r);

    /* beta = <r1,tilder> / <r0,tilder> * alpha / omega */
    rtilder1 = vec_dot(n, r, tilder);
    beta = rtilder1 / rtilder0 * alpha / omega;

    /* p = r + beta * (p - omega * A * p) */
    vec_copy(n, p, s);
    vec_axpy(n, momega, Ap, s);
    vec_copy(n, r, p);
    vec_axpy(n, beta, s, p);

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
  free(tilder);
  free(p);
  free(s);
  free(Ap);
  free(As);
}


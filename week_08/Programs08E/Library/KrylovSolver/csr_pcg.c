 /*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void csr_pcg(csr *A, prec_t *prec, double *b, double *x,
             double tol, int stopcrit, int maxit,
             double *esterr, int *iter, double *errvec)
{
  int k;
  int n;
  double rho = 0.0;
  double alpha, beta, rz0, rz1;
  double *r, *z, *p, *Ap;
  double *temp;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate vectors. */
  r = (double *) malloc(n * sizeof(double));
  z = (double *) malloc(n * sizeof(double));
  p = (double *) malloc(n * sizeof(double));
  Ap = (double *) calloc(n, sizeof(double));

  /* Residual r^0. */
  vec_copy(n, b, r);
  csr_mv(A, -1.0, x, 1.0, r);

  /* Preconditioned residual z^0 = M \ r^{0}. */
  prec_apply(prec, r, z);

  /* Direction vector p^0 = z^0. */
  vec_copy(n, z, p);

  /* Compute residual norm if requested. */
  if (stopcrit == ITSOLVER_RELRES) {
    rho = vec_norm(n, r);
    errvec[0] = rho;
  }

  /* Scalar product <r^{k}, z^{k}> and <r^{k+1}, z^{k+1}>. */
  rz0 = 0.0; /* Will not be used. */
  rz1 = vec_dot(n, r, z);

  /* Start the iteration. */
  k = 0;
  while (k < maxit) {

    /* Check whether the desired accuracy has been reached. */
    if (stopcrit == ITSOLVER_RELRES && errvec[k] < tol * rho) {
      break;
    }

    /* One step of PCG. */
    rz0 = rz1;
    csr_mv(A, 1.0, p, 0.0, Ap);
    alpha = rz0 / vec_dot(n, p, Ap);
    vec_axpy(n, alpha, p, x);
    alpha *= -1;
    vec_axpy(n, alpha, Ap, r);
    prec_apply(prec, r, z);
    rz1 = vec_dot(n, r, z);
    beta = rz1 / rz0;
    vec_axpy(n, beta, p, z);
    /* Swap pointers p and z. */
    temp = p;
    p = z;
    z = temp;

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
  free(Ap);
  free(p);
  free(z);
  free(r);
}

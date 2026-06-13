/* 
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
*/

#include "itsolver.h"

void csr_cg(csr *A, double *b, double *x,
            double tol, int stopcrit, int maxit,
            double *err, int *iter, double *errvec)
{

  double pAp, alpha, beta;
    
  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  int n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate residual and direction vector. */
  double *r  = malloc(n * sizeof(double));
  double *p  = malloc(n * sizeof(double));
  double *Ap = calloc(n,  sizeof(double));

  /* Compute start residual r^0. */
  vec_copy(n, b, r);
  csr_mv(A, -1.0, x, 1.0, r);

  /* Compute start direction vector p^0 = r^0. */
  vec_copy(n, r, p);

  /* Squares of residual vectors for the current and for the next iteration index. */
  double rsq0 = 0.0; /* Will not be used. */
  double rsq1 = vec_dot(n, r, r);

  /* Compute residual norm. */
  double rho = sqrt(rsq1);
  double rho0 = rho;
  if (errvec != NULL) {
    errvec[0] = rho0;
  }

  /* Start the iteration. */
  int k = 0;
  while (k < maxit) {
    /* Shift. */
    rsq0 = rsq1;

    /* Check whether the desired accuracy has been reached. */
    if (stopcrit == ITSOLVER_RELRES && rho <= tol * rho0) {
      break;
    }

    /* One step of CG. */
    csr_mv(A, 1.0, p, 0.0, Ap);
    pAp = vec_dot(n, p, Ap);
    assert(pAp != 0.0);
    alpha = rsq0 / pAp;
    vec_axpy(n,  alpha,  p, x);
    vec_axpy(n, -alpha, Ap, r);
    rsq1 = vec_dot(n, r, r);
    assert(rsq0 != 0.0);
    beta = rsq1 / rsq0;
    /* Use Ap as a temporary. */
    vec_copy(n, r, Ap);
    vec_axpy(n, beta, p, Ap);
    /* Swap pointers p and Ap. */
    double *temp = p; 
    p = Ap;
    Ap = temp; 

    /* Compute residual norm. */
    rho = vec_norm(n, r); /* == sqrt(rsq1) */
    if (errvec != NULL) {
      errvec[k + 1] = rho;
    }

    /* Proceed to the next step. */
    k++;
  }

  /* Output iter and esterr. */
  if (iter != NULL) {
    *iter = k;
  }
  if (err != NULL) {
    *err = rho;
  }

  /* Free memory. */
  free(Ap);
  free(p);
  free(r);
}

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

  double pAp, alpha, beta, old_r_vec_dot;

  /* Number of rows and columns. */
  assert(A->nr == A->nc);
  int n = A->nr;

  assert(stopcrit == ITSOLVER_NORES || stopcrit == ITSOLVER_RELRES);

  /* Allocate residual and direction vector. */
  double *r = malloc(n * sizeof(double));
  double *p = malloc(n * sizeof(double));
  double *Ap = calloc(n, sizeof(double));

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
  if (errvec != NULL)
  {
    errvec[0] = rho0;
  }

  /* Start the iteration. */
  int k = 0;
  while (k < maxit)
  {
    /* Shift. */
    rsq0 = rsq1;

    /* Check whether the desired accuracy has been reached. */
    if (stopcrit == ITSOLVER_RELRES && rho <= tol * rho0)
    {
      break;
    }

    /* One step of CG. */
    
    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 2          */
    /* Implement one step of the Conjugate  */
    /* Gradient method.                     */
    /*                                      */
    /* ************************************ */
    //setup vector Ap to be used in alpha calculation
    csr_mv(A,1,p,0,Ap);
    alpha = vec_dot(n, r, r) / vec_dot(n, Ap, p);

    old_r_vec_dot = vec_dot(n, r,r); //we save the dot product of r,r before updating it
    vec_axpy(n, alpha, p, x); //here we do x <-- x + alpha*p
    csr_mv(A, -alpha, p,1, r); // r <-- r - alpha*A*p

    beta = vec_dot(n, r,r) /old_r_vec_dot;   // beta = vec_dot(n, new_r, new_r) / old_r_vec_dot; denominator is the product given by old r

    vec_scale(n, beta, p);
    vec_axpy(n, 1, r, p);

    /* Compute residual norm. */
    rho = vec_norm(n, r); /* == sqrt(rsq1) */
    if (errvec != NULL)
    {
      errvec[k + 1] = rho;
    }

    /* Proceed to the next step. */
    k++;
  }

  /* Output iter and esterr. */
  if (iter != NULL)
  {
    *iter = k;
  }
  if (err != NULL)
  {
    *err = rho;
  }

  /* Free memory. */
  free(Ap);
  free(p);
  free(r);
}

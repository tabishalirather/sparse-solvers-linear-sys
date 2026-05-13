/*
 * Jens-Peter Zemke
 * Sabine Le Borne
*/
/*
 * Richardson Method for solving Ax = b
 *
 * The iteration formula is:
 *   x_new = x_old + omega * (b - A*x_old)
 * where (b - A*x) is the residual vector
 *
 * Two stopping criteria:
 *
 * NORES (No Residual check):
 *   - Just run maxit iterations, don't check anything
 *   - Fast but you don't know how good the solution is
 *
 * RELRES (Relative Residual):
 *   - Compute residual at each iteration
 *   - Stop early when: (current residual) / (initial residual) < tolerance
 *   - Smarter because it accounts for how bad the initial guess was
 *
 * Outputs:
 *   - x: the approximate solution (updated)
 *   - iter: which iteration we stopped at
 *   - esterr: final relative error (final residual / initial residual)
 *   - errvec: array storing all residual norms (RELRES only)
 */

#include "itsolver.h"

/*
  A, b, x    : data for every iterative scheme:
               A - system matrix in CSR format;
               b - right-hand side;
               x - starting guess / approximate solution.
  stopcrit   : int parameter to choose the stopping criterion:
               ITSOLVER_PREC - for use as a preconditioner;
               ITSOLVER_NORES - no computation of a residual;
               ITSOLVER_RELRES - relative residual;
               ITSOLVER_BACKERR - backward error.
  tol, maxit : stopping criteria for every scheme;
  omega      : parameter for relaxed Richardson's iteration;
  esterr     : used by tol for stopping criterion;
  iter       : iteration number reached;
  errvec     : vector of all esterr of size maxit+1
               or untouched (ITSOLVER_NORES);
*/

void csr_richardson(csr *A, double *b, double *x,
  double tol, int stopcrit, int maxit,
  double omega,
  double *esterr, int *iter, double *errvec)
{
  int i, n;
  double *r;
  double rho;

  n = A->nc;

  switch (stopcrit)
  {
  case ITSOLVER_NORES:
    r = (double *) malloc(n*sizeof(double));
    /* for i=0,...,maxit
         r <- b;
         r <- r-A*x;
         x <- x+omega*r;
       end for
     */
    for (i=0; i<maxit; i++)
    {
      vec_copy(n, b, r);
      csr_mv(A,-1.0,x,1.0,r);
      vec_axpy(n, omega, r, x);
    }
    free(r);
    break;
  case ITSOLVER_RELRES:
    r = (double *) malloc(n*sizeof(double));
    /* r <- b;
       r <- r-A*x;
       compute norm(r,2);
       x <- x+omega*r;
       for i=1,...,maxit
         r <- b;
         r <- r-A*x;
         compute norm(r,2);
         x <- x+omega*r;
       end for
     */
    vec_copy(n, b, r);
    csr_mv(A,-1.0,x,1.0,r);
    rho = vec_norm(n, r);
    *errvec++ = rho;
    tol = tol*rho;
    *iter = maxit;
    for (i=0; i<maxit; i++)
    {
      vec_axpy(n, omega, r, x);
      vec_copy(n, b, r);
      csr_mv(A,-1.0,x,1.0,r);
      *errvec = vec_norm(n, r);
      if ( *errvec++ <= tol )
      {
        *iter = i+1;
        break;
      }
    }
    *esterr = *(--errvec)/rho;

    /* no need to move back errvec, is a copy */

    free(r);
    break;
  }

} /* csr_richardson */

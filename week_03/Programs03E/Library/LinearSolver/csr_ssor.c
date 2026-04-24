/*
 * Jens-Peter Zemke
 * Sabine Le Borne
*/

#include "itsolver.h"

/*
  A, b, x    : data for every iterative scheme:
               A - system matrix in CSR format;
               b - right-hand side;
               x - starting guess / approximate solution.
  stopcrit   : int parameter to choose the stopping criterion:
               ITSOLVER_NORES - no computation of a residual;
               ITSOLVER_RELRES - relative residual;
               ITSOLVER_BACKERR - backward error.
  tol, maxit : stopping criteria for every scheme;
  esterr     : used by tol for stopping criterion;
  iter       : iteration number reached;
  errvec     : vector of all esterr of size maxit+1
               or untouched (ITSOLVER_NORES);
*/

void csr_ssor(csr *A, double *b, double *x,
  double tol, int stopcrit, int maxit, double omega,
  double *esterr, int *iter, double *errvec)
{
  int i, j, k, n;
  double *r;
  double oneminusomega;
  double rho, Aii, xtmp;
  double normA, normb;

  n = A->nc;
  oneminusomega = 1.0-omega;

  /* our implementation relies on the assumption that Aii is always
     the first element in a row */
  csr_sort_diagFirst(A); 
  
  switch (stopcrit)
  {
  case ITSOLVER_NORES:
    for (j=0; j<maxit; j++)
    {

      for (i=0; i<A->nr; i++)
      {

        xtmp = x[i];
        x[i] = *b++;

        /* we rely on the assumption that Aii is always
           the first element in a row */

        Aii = *A->coeff++; A->col++;
        for (k=1+*A->row++; k<*A->row; k++)
        {
          x[i] -= *A->coeff++ * x[*A->col++];
        }

        /* compute relaxed update
           x = r + omega*(x-r)
             = (1-omega)*r + omega*x */

        x[i] *= (omega/Aii);
        x[i] += (oneminusomega * xtmp);

      }

      /* move pointers to the first entry */

      A->row -= A->nr;
      A->col -= A->nz;
      A->coeff -= A->nz;
      b -= A->nr;

      for (i=A->nr-1; i>-1; i--)
      {

        xtmp = x[i];
        x[i] = b[i];

        /* we rely on the assumption that Aii is always
           the first element in a row */

        Aii = A->coeff[A->row[i]];
        for (k=A->row[i]+1; k<A->row[i+1]; k++)
        {
          x[i] -= A->coeff[k] * x[A->col[k]];
        }

        /* compute relaxed update
           x = r + omega*(x-r)
             = (1-omega)*r + omega*x */

        x[i] *= (omega/Aii);
        x[i] += (oneminusomega * xtmp);

      }

    } /* for j */

    break;
  case ITSOLVER_RELRES:
    r = (double *) malloc(n*sizeof(double));

    /* compute zeroth residual */

    vec_copy(n, b, r);
    csr_mv(A,-1.0,x,1.0,r);
    rho = vec_norm(n, r);
    *errvec++ = rho;
    tol = tol*rho;
    *iter = maxit;

    for (j=0; j<maxit; j++)
    {

      for (i=0; i<A->nr; i++)
      {

        xtmp = x[i];
        x[i] = *b++;

        /* we rely on the assumption that Aii is always
           the first element in a row */

        Aii = *A->coeff++; A->col++;
        for (k=1+*A->row++; k<*A->row; k++)
        {
          x[i] -= *A->coeff++ * x[*A->col++];
        }

        /* compute relaxed update
           x = r + omega*(x-r)
             = (1-omega)*r + omega*x */

        x[i] *= (omega/Aii);
        x[i] += (oneminusomega * xtmp);

      }

      /* move pointers to the first entry */

      A->row -= A->nr;
      A->col -= A->nz;
      A->coeff -= A->nz;
      b -= A->nr;

      for (i=A->nr-1; i>-1; i--)
      {

        xtmp = x[i];
        x[i] = b[i];

        /* we rely on the assumption that Aii is always
           the first element in a row */

        Aii = A->coeff[A->row[i]];
        for (k=A->row[i]+1; k<A->row[i+1]; k++)
        {
          x[i] -= A->coeff[k] * x[A->col[k]];
        }

        /* compute relaxed update
           x = r + omega*(x-r)
             = (1-omega)*r + omega*x */

        x[i] *= (omega/Aii);
        x[i] += (oneminusomega * xtmp);

      }

      vec_copy(n, b, r);
      csr_mv(A,-1.0,x,1.0,r);
      *errvec = vec_norm(n, r);
      if ( *errvec++ <= tol )
      {
        *iter = j+1;
        break;
      }

    } /* for j */

    *esterr = *(--errvec)/rho;

    /* no need to move back errvec, is a copy */

    free(r);
    break;
  case ITSOLVER_BACKERR:
    r = (double *) malloc(n*sizeof(double));
    normb = vec_norm(n, b);
    normA = csr_norm_infty(A);

    /* compute zeroth backward error */

    vec_copy(n, b, r);
    csr_mv(A,-1.0,x,1.0,r);
    *errvec =
      vec_norm(n, r)/(normA*vec_norm(n, x)+normb);
    if ( *errvec++ <= tol )
    {
      maxit = 0;
      *iter = 0;
    }

    for (j=0; j<maxit; j++)
    {

      for (i=0; i<A->nr; i++)
      {

        xtmp = x[i];
        x[i] = *b++;

        /* we rely on the assumption that Aii is always
           the first element in a row */

        Aii = *A->coeff++; A->col++;
        for (k=1+*A->row++; k<*A->row; k++)
        {
          x[i] -= *A->coeff++ * x[*A->col++];
        }

        /* compute relaxed update
           x = r + omega*(x-r)
             = (1-omega)*r + omega*x */

        x[i] *= (omega/Aii);
        x[i] += (oneminusomega * xtmp);

      }

      /* move pointers to the first entry */

      A->row -= A->nr;
      A->col -= A->nz;
      A->coeff -= A->nz;
      b -= A->nr;

      for (i=A->nr-1; i>-1; i--)
      {

        xtmp = x[i];
        x[i] = b[i];

        /* we rely on the assumption that Aii is always
           the first element in a row */

        Aii = A->coeff[A->row[i]];
        for (k=A->row[i]+1; k<A->row[i+1]; k++)
        {
          x[i] -= A->coeff[k] * x[A->col[k]];
        }

        /* compute relaxed update
           x = r + omega*(x-r)
             = (1-omega)*r + omega*x */

        x[i] *= (omega/Aii);
        x[i] += (oneminusomega * xtmp);

      }

      /* compute backward error */

      vec_copy(n, b, r);
      csr_mv(A,-1.0,x,1.0,r);
      *errvec =
        vec_norm(n, r)/(normA*vec_norm(n, x)+normb);
      if ( *errvec++ <= tol )
      {
        *iter = j+1;
        break;
      }

    } /* for j */

    *esterr = *(--errvec);

    /* no need to move back errvec, is a copy */

    free(r);
    break;
  }

} /* csr_ssor */

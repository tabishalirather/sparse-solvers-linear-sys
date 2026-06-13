/*! \file
 *  \brief Investigation of convergence behaviour of CG when the right hand
 *         side is a combination of eigenvectors.
 *
 * The number of eigenvectors is increased from few to many to all.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NPLOTS (10)

int n, d, nd;
int i;
csr *A;
double *b, *x, *xstart;
int *eigidx, eigidxmin, eigidxmax;
double *eigvec;
double tol;
int stopcrit;
int maxit[NPLOTS];
double esterr;
int iter[NPLOTS];
double *errvec[NPLOTS];
char *legend_entries[NPLOTS];
/* Blas constants. */
int inc_one = 1;
double alpha_one = 1.0;

/* k: Run index = 0, 1, 2. */
void cg_rhs(int k, int low) {

  int eignum;

  assert(0 <= k && k <= n - 2);

  /* Set right hand side to zero. */
  memset(b, 0, nd * sizeof(double));

  /* Number of eigenvectors contained in the right hand side. */
  eignum = 1;
  for (i = 0; i < d; ++i) {
    eignum *= (k + 1);
  }

  if (low != 1) {
    /* use eignum  highest eigenvector indices */
    /* Initialize eigenvector multiindex. */
    eigidxmin = n - 2 - k;
    for (i = 0; i < d; ++i) {
      eigidx[i] = eigidxmin;
    }

    while (1) {
      /* Add eigenvector to right hand side. */
      laplace_ev(n, d, eigidx, eigvec);
      vec_axpy(nd, alpha_one, eigvec, b);
 
      /* Proceed to the next multiindex. */
      for (i = 0; i < d && eigidx[i] == n - 2; ++i) {
        eigidx[i] = eigidxmin;
      }
    
      if (i == d) {
        break;
      }
      ++eigidx[i];
     }
  }
  else {
    /* use eignum lowest eigenvector indices */
    
  /* ******************************************************************** */
  /*                                                                      */
  /*     TODO --- Exercise (Exercise sheet 7, Problem 2)                */
  /*                                                                      */
  /* ******************************************************************** */
   
  }

  legend_entries[k] = malloc(100);
  sprintf(legend_entries[k], "%04d eigen shares", eignum);
  printf("%s\n", legend_entries[k]);
  maxit[k] = 100;
  errvec[k] = (double *) malloc((maxit[k] + 1) * sizeof(double));
  memcpy(x, xstart, nd * sizeof(double));
  csr_cg(A, b, x, tol, stopcrit, maxit[k],
   &esterr, &iter[k], errvec[k]);
}

int main() {
  int k;
  double kappa, kappa_sqrt, rho, rho_log10;

  /* Laplace. */
  n = 11;
  d = 3;
  A = csr_laplace(n, d);
  nd = A->nc;
  printf("n = %d, d = %d, nd = %d\n", n, d, nd);

  /* Convergence rate of CG. */
  kappa = laplace_cond(n, d);
  kappa_sqrt = sqrt(kappa);
  rho = (kappa_sqrt - 1.0) / (kappa_sqrt + 1);
  rho_log10 = log10(rho);
  printf("Predicted convergence rate of CG:\n");
  printf("rho = %g\n", rho);
  printf("log10(rho) = %g\n", rho_log10);

  /* Right hand side. */
  b = (double *) malloc(nd * sizeof(double));
  eigidx = (int *) malloc(d * sizeof(int));
  eigvec = (double *) malloc(nd * sizeof(double));

  /* Zero start vector. */
  x = (double *) malloc(nd * sizeof(double));
  xstart = (double *) malloc(nd * sizeof(double));
  for (i = 0; i < nd; ++i) {
    xstart[i] = 0.0;
  }

  /* Tolerance and stopping criterium. */
  tol = 1e-16;
  stopcrit = ITSOLVER_RELRES;

  /* Run CG multiple times. */
  for (k = 0; k < NPLOTS; ++k) {
    cg_rhs(k, 0);
  }

  /* Combined plot of all the methods. */
  convplot_steps(
    "testprog_06_cg_rhs.eps",
    "CG method with different right hand sides, highest eigenvalues",
    "Iteration steps",
    "Residual norm",
    NPLOTS,
    iter,
    errvec,
    legend_entries);

  /* Free memory. */
  for (k = 0; k < NPLOTS; ++k) {
    free(errvec[k]);
    free(legend_entries[k]);
  }
  free(xstart);
  free(x);
  free(eigvec);
  free(eigidx);
  free(b);
  csr_free(A);

  return 0;
}

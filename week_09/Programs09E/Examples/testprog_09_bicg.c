/*! \file
 *  \brief Comparison of CG, BiCG, CGS and BiCGStab.
 *
 *  The tested matrix comes from an oil reservoir simulation
 *  problem (ORSIRR1 from the Harwell-Boeing collection).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

#define K_MAX 5

int main() {
  char *file_mtx;
  int k;
  int i, n;
  csr *A;
  csr *e;
  double *b;
  double *x, *x0;
  double tol;
  int stopcrit;
  int maxit[K_MAX];
  double esterr;
  int iter[K_MAX];
  double *errvec[K_MAX];
  char *legend_entries[K_MAX];

  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Import the testmatrix orsirr_1, and  */
  /* define the right-hand side.          */
  /*                                      */
  /* ************************************ */
  
  file_mtx = "../../TestMatrixCollection/orsirr_1.mtx";
  A = csr_import_from_mm_mtx(file_mtx);
  e = vec_ones(n);
  if (A == NULL) {
    printf("Matrix import failed\n");
  }

  b = csr
  /* Random start vector. */
  x0 = malloc(n * sizeof(double));
  for (i = 0; i < n; ++i) {
    x0[i] = drand48();
  }

  /* Stopping. */
  tol = 1e-7;
  stopcrit = ITSOLVER_RELRES;

  /* Counter for the current iterative schmeme. */
  k = 0;

  /* CG. */
  printf("CG:\n");
  maxit[k] = 1000;
  errvec[k] = malloc((maxit[k] + 1) * sizeof(double));
  vec_copy(n, x0, x);
  csr_cg(A, b, x, tol, stopcrit, maxit[k], &esterr, &iter[k], errvec[k]);
  legend_entries[k] = "CG";
  printf("  Iterations: %d\n", iter[k]);
  k++;

  /* BiCG. */
  printf("BiCG:\n");
  maxit[k] = 1000;
  errvec[k] = malloc((maxit[k] + 1) * sizeof(double));
  vec_copy(n, x0, x);
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Call the BiCG method in order to     */
  /* solve the linear system.             */
  /*                                      */
  /* ************************************ */
  csr_bicg(A, b)
  legend_entries[k] = "BiCG";
  printf("  Iterations: %d\n", iter[k]);
  k++;

  /* CGS. */
  printf("CGS:\n");
  maxit[k] = 1000;
  errvec[k] = malloc((maxit[k] + 1) * sizeof(double));
  vec_copy(n, x0, x);
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Call the CGS method in order to      */
  /* solve the linear system.             */
  /*                                      */
  /* ************************************ */
  
  legend_entries[k] = "CGS";
  printf("  Iterations: %d\n", iter[k]);
  k++;

  /* BiCGStab. */
  printf("BiCGStab:\n");
  maxit[k] = 1000;
  errvec[k] = malloc((maxit[k] + 1) * sizeof(double));
  vec_copy(n, x0, x);
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Call the BiCGStab method in order    */
  /* to solve the linear system.          */
  /*                                      */
  /* ************************************ */

  legend_entries[k] = "BiCGStab";
  printf("  Iterations: %d\n", iter[k]);
  k++;

  /* GMRes. */
  printf("GMRes:\n");
  maxit[k] = 1000;
  errvec[k] = malloc((maxit[k] + 1) * sizeof(double));
  gmres_t *gmres = gmres_new(A, b, x0);
  gmres_tol(gmres, tol);
  gmres_iter_max(gmres, maxit[k]);
  gmres_solve(gmres);
  iter[k] = gmres_iter(gmres);
  vec_copy(iter[k] + 1, gmres_residual_norms(gmres), errvec[k]);
  gmres_delete(gmres);
  legend_entries[k] = "GMRes";
  printf("  Iterations: %d\n", iter[k]);
  k++;

  /* Plot. */
  assert(k == K_MAX);
  convplot_steps(
     "testprog_09_bicg.eps",
     "Matrix ORSIRR1 (oil reservoir simulation)",
     "Iteration",
     "Residual norm",
     K_MAX,
     iter,
     errvec,
     legend_entries);

  return 0;
}

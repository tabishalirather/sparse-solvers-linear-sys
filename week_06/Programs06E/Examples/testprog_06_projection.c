/*
 * \author Jens-Peter Zemke
 * \author Sabine Le Borne
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "itsolver.h"

int main(void)
{
  csr *A;       // Matrix.
  int i;        // Loop index.
  int j = 0;    // Loop index.
  int k = 6;    // Number of methods.

  // Variables for the plot.
  int *m;
  char *name[] = {"Richardson",
                  "Jacobi",
                  "Gauss-Seidel",
                  "Steepest Descent",
                  "Minimal Residual",
                  "Residual Norm Steepest Descent"};
  double **vecvec;

  // Variables for the methods.
  int iter, maxit = 100;
  double tol = 1e-18;
  double esterr;

  // Allocate space for matrix, solution, and right-hand side.
  int n = 10;
  int d = 3;
  A  = csr_laplace(n + 1, d);
  double *b  = (double *) malloc(A->nc*sizeof(double));
  double *x  = (double *) malloc(A->nc*sizeof(double));

  // allocate space for plot-data.
  m = (int *) malloc(k*sizeof(int));
  vecvec = (double **) malloc(k*sizeof(double *));
  for (i=0; i<k; i++) {
    vecvec[i] = (double *) malloc((maxit+1)*sizeof(double *));
  }

  // set up b = A*e, e vector of ones.
  for (i=0; i<A->nc; i++)
    *x++ = 1.0;
  x -= A->nc;
  csr_mv(A,1.0,x,0.0,b);



  // Richardson for optimal omega.
  
  // Eigenvalues and optimal omega for Richardson's iteration.
  double lambda_max, lambda_min;
  double omega;
  lambda_min = 4 * d * pow(sin(          M_PI / (2 * n)), 2);
  lambda_max = 4 * d * pow(sin((n - 1) * M_PI / (2 * n)), 2);
  omega = 2.0 / (lambda_min + lambda_max);

  printf("\n");
  printf("Richardson: omega = %g.\n", omega);

  m[j] = maxit;
  vec_scale(A->nc, 0.0, x);
  csr_richardson(A, b, x, tol, ITSOLVER_RELRES, m[j], omega,
    &esterr, &iter, vecvec[j]);
  m[j] = iter;
  printf("\n");
  printf("Richardson:\n");
  printf("===========\n");
  printf("%d steps were computed.\n", iter);
  printf("Maximal of steps was %d.\n", maxit);
  printf("Start residual norm = %g.\n", vecvec[j][0]);
  printf("Residual norm = %g.\n", vecvec[j][iter]);
  printf("Relative Residual norm = %g.\n", esterr);
  printf("Tolerance = %g.\n", tol);
  printf("x(1:4) = [%g,%g,%g,%g]\n",x[0],x[1],x[2],x[3]);
  j++;



  // Jacobi.

  m[j] = maxit;
  vec_scale(A->nc, 0.0, x);
  csr_jacobi(A, b, x, tol, ITSOLVER_RELRES, m[j],
    &esterr, &iter, vecvec[j]);
  m[j] = iter;
  printf("\n");
  printf("Jacobi:\n");
  printf("=======\n");
  printf("%d steps were computed.\n", iter);
  printf("Maximal of steps was %d.\n", maxit);
  printf("Start residual norm = %g.\n", vecvec[j][0]);
  printf("Residual norm = %g.\n", vecvec[j][iter]);
  printf("Relative Residual norm = %g.\n", esterr);
  printf("Tolerance = %g.\n", tol);
  printf("x(1:4) = [%g,%g,%g,%g]\n",x[0],x[1],x[2],x[3]);
  j++;



  // Gauß-Seidel.

  m[j] = maxit;
  vec_scale(A->nc, 0.0, x);
  csr_gauss_seidel(A, b, x, tol, ITSOLVER_RELRES, m[j],
    &esterr, &iter, vecvec[j]);
  m[j] = iter;
  printf("\n");
  printf("Gauss-Seidel:\n");
  printf("=============\n");
  printf("%d steps were computed.\n", iter);
  printf("Maximal of steps was %d.\n", maxit);
  printf("Start residual norm = %g.\n", vecvec[j][0]);
  printf("Residual norm = %g.\n", vecvec[j][iter]);
  printf("Relative Residual norm = %g.\n", esterr);
  printf("Tolerance = %g.\n", tol);
  printf("x(1:4) = [%g,%g,%g,%g]\n",x[0],x[1],x[2],x[3]);
  j++;



  // * Steepest descent.

  m[j] = maxit;
  vec_scale(A->nc, 0.0, x);
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Call the function implementing the   */
  /* Steepest Descent method.             */
  /*                                      */
  /* ************************************ */

  csr_steepest_descent(A, b, x, tol, ITSOLVER_RELRES, maxit, &esterr, &iter, vecvec[j]);
  
  m[j] = iter;
  printf("\n");
  printf("Steepest Descent:\n");
  printf("=================\n");
  printf("%d steps were computed.\n", iter);
  printf("Maximal of steps was %d.\n", maxit);
  printf("Start residual norm = %g.\n", vecvec[j][0]);
  printf("Residual norm = %g.\n", vecvec[j][iter]);
  printf("Relative Residual norm = %g.\n", esterr);
  printf("Tolerance = %g.\n", tol);
  printf("x(1:4) = [%g,%g,%g,%g]\n",x[0],x[1],x[2],x[3]);
  j++;


  
  // * Minimal residual.
  
  m[j] = maxit;
  vec_scale(A->nc, 0.0, x);
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Call the function implementing the   */
  /* Minimal Resudial method.             */
  /*                                      */
  /* ************************************ */

  csr_minimal_residual(A, b, x, tol, ITSOLVER_RELRES, maxit, &esterr, &iter, vecvec[j]);


  m[j] = iter;
  printf("\n");
  printf("Minimal Residual:\n");
  printf("=================\n");
  printf("%d steps were computed.\n", iter);
  printf("Maximal of steps was %d.\n", maxit);
  printf("Start residual norm = %g.\n", vecvec[j][0]);
  printf("Residual norm = %g.\n", vecvec[j][iter]);
  printf("Relative Residual norm = %g.\n", esterr);
  printf("Tolerance = %g.\n", tol);
  printf("x(1:4) = [%g,%g,%g,%g]\n",x[0],x[1],x[2],x[3]);
  j++;



  // * Residual norm steepest descent.
  
  m[j] = maxit;
  vec_scale(A->nc, 0.0, x);
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Call the function implementing the   */
  /* Residual Norm Steepest Descent       */
  /* method.                              */
  /*                                      */
  /* ************************************ */
  csr_residual_norm_steepest_descent(A, b, x, tol, ITSOLVER_RELRES, maxit, &esterr, &iter, vecvec[j]);

  m[j] = iter;
  printf("\n");
  printf("Residual Norm Steepest Descent:\n");
  printf("===============================\n");
  printf("%d steps were computed.\n", iter);
  printf("Maximal of steps was %d.\n", maxit);
  printf("Start residual norm = %g.\n", vecvec[j][0]);
  printf("Residual norm = %g.\n", vecvec[j][iter]);
  printf("Relative Residual norm = %g.\n", esterr);
  printf("Tolerance = %g.\n", tol);
  printf("x(1:4) = [%g,%g,%g,%g]\n",x[0],x[1],x[2],x[3]);
  j++;



  // Plot.
  convplot_steps(
    "testprog_projection.eps",
    "Linear iteration methods vs. one-dimensional projection methods",
    "Steps",
    "Residual norms", j, m, vecvec, name);

  // Free space for linear system.
  free(x);
  free(b);
  csr_free(A);

  // Free space for convergence curves.
  free(m);
  for (i=0; i<k; i++) free(vecvec[i]);
  free(vecvec);

  return EXIT_SUCCESS;
} /* main */

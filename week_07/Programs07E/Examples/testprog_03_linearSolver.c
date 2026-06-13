/*! \file
 *
 *   \author Michael Wende 
 *   \author Jens-Peter Zemke 
 *   \author Sabine Le Borne
 */


#include "itsolver.h"

int main() {
  int n, d, nd;
  double omega;
  int i, iter, maxit;
  csr *A;
  double *b, *x, *xstart;
  double tol, esterr, *errvec;
  int stopcrit;

  /* Symmetric positive definite matrix. */
  n = 10;
  d = 3;
  A = csr_laplace(n, d);
  nd = A->nc;

  /* Zero right hand side. */
  b = (double *) malloc(nd * sizeof(double));
  for (i = 0; i < nd; ++i) {
    b[i] = 0.0;
  }

  /* All ones start vector. */
  x = (double *) malloc(nd * sizeof(double));
  xstart = (double *) malloc(nd * sizeof(double));
  for (i = 0; i < nd; ++i) {
    xstart[i] = 1.0;
  }

  /* Stopping criterion. */
  tol = 1e-4;
  stopcrit = ITSOLVER_RELRES;
  maxit = 100;
  errvec = malloc((maxit + 1) * sizeof(double));

  /* Create plot. */
  conv_plot_t *conv_plot = conv_plot_new();
  conv_plot_titlef(conv_plot, "Linear iteration methods (3d Laplace operator, n=%d)", nd);
  conv_plot_xlabel(conv_plot, "Iteration steps");
  conv_plot_ylabel(conv_plot, "Absolute residual norm");

  { /* Jacobi. */
    vec_copy(nd, xstart, x);
    vec_scale(maxit, 0.0, errvec);
    csr_jacobi(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "Jacobi");
  }

  { /* Gauss-Seidel (vorwaerts). */
    vec_copy(nd, xstart, x);
    vec_scale(maxit, 0.0, errvec);
    csr_gauss_seidel(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "GS");
  }

  { /* Gauss-Seidel (rueckwaerts). */
    /*
    vec_copy(nd, xstart, x);
    vec_scale(maxit, 0.0, errvec);
    csr_gauss_seidel_reverse(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "rGS");
    */
  }

  { /* Gauss-Seidel (symmetrisch). */
    vec_copy(nd, xstart, x);
    vec_scale(maxit, 0.0, errvec);
    csr_sgs(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "SGS");
  }

  { /* SOR (omega). */
    omega = 1.5;
    vec_copy(nd, xstart, x);
    vec_scale(maxit, 0.0, errvec);
    csr_sor(A, b, x, tol, stopcrit, maxit, omega, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "SOR (omega=%.4f)", omega);
  }

  { /* SSOR (omega). */
    omega = 1.5;
    vec_copy(nd, xstart, x);
    vec_scale(maxit, 0.0, errvec);
    csr_ssor(A, b, x, tol, stopcrit, maxit, omega, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "SSOR (omega=%.4f)", omega);
  }

  /* Save plot to file. */
  conv_plot_save(conv_plot, "testprog_linearSolver");
  conv_plot_delete(conv_plot);

  /* Free memory. */
  free(errvec);
  free(x);
  free(b);
  free(xstart);
  csr_free(A);

  return 0;
}

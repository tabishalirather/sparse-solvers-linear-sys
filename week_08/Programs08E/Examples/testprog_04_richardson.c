/*! \file
 *
 *   \author Michael Wende 
 *   \author Jens-Peter Zemke 
 *   \author Sabine Le Borne
 */

#include "itsolver.h"

int main() {
  int n, d, nd;
  double lambda_min, lambda_max;
  double omega_opt, omega_crit;
  int i;
  csr *a;
  double *b, *x, *xstart;
  double tol;
  int stopcrit;

  /* Create plot. */
  conv_plot_t *conv_plot = conv_plot_new();
  conv_plot_title(conv_plot, "Richardson method applied to the discretization of the Laplace operator");
  conv_plot_xlabel(conv_plot, "Iteration steps");
  conv_plot_ylabel(conv_plot, "Absolute residual norm");

  /* Symmetric positive definite matrix. */
  n = 10;
  d = 3;
  a = csr_laplace(n, d);
  nd = a->nc;

  /* Compute optimal relaxation parameter. */
  lambda_min = 4 * d * pow(sin(          M_PI / (2 * n)), 2);
  lambda_max = 4 * d * pow(sin((n - 1) * M_PI / (2 * n)), 2);
  omega_opt = 2.0 / (lambda_min + lambda_max);
  omega_crit = 2.0 / lambda_max;
  printf("n = %d, d = %d, nd = %d\n", n, d, nd);
  printf("Optimal convergence for omega = %g\n", omega_opt);
  printf("No convergence for omega >= %g\n", omega_crit);

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

  /* Perform Richardson's iteration. */
  tol = 1e-6;
  stopcrit = ITSOLVER_RELRES;

  { /* Omega too small. */
    double omega = 0.8 * omega_opt;
    int maxit = 1000;
    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(nd, xstart, x);
    double esterr;
    int iter;
    csr_richardson(a, b, x, tol, stopcrit, maxit, omega,
                   &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "omega = %.4f", omega);
    free(errvec);
  }

  { /* Omega just right. */
    double omega = 1.0 * omega_opt;
    int maxit = 1000;
    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(nd, xstart, x);
    double esterr;
    int iter;
    csr_richardson(a, b, x, tol, stopcrit, maxit, omega,
                   &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "omega = %.4f", omega);
    free(errvec);
  }

  { /* Omega too large. */
    double omega = 1.05 * omega_opt;
    int maxit = 150; /* Prevent plot window from getting too large. */
    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(nd, xstart, x);
    double esterr;
    int iter;
    csr_richardson(a, b, x, tol, stopcrit, maxit, omega,
                   &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "omega = %.4f", omega);
    free(errvec);
  }

  /* Save plot to file. */
  conv_plot_save(conv_plot, "testprog_richardson");
  conv_plot_delete(conv_plot);

  /* Free memory. */
  free(x);
  free(b);
  free(xstart);
  csr_free(a);

  return 0;
}

#include "itsolver.h"

int main() {
  int n, d, nd;
  double lambda_min, lambda_max;
  double omega_opt;
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

  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /* Use the Richardson method on the     */
  /* Laplace matrix for different damping */
  /* parameters.                          */
  /*                                      */
  /* ************************************ */

  // TODO: Generate discrete Laplace matrix.

  // TODO: Compute optimal damping parameter.
  
  // TODO: Zero right hand side.

  // TODO: All ones start vector.



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

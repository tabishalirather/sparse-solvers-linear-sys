/*! \file
 *  \brief Iterative solution of convection diffusion test problem.
 *
 *  Comparision of linear iterations (Jacobi, Gauss-Seidel, ...).
 *
 *  \author Michael Wende 
 *  \author Jens Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

static void test_convection_diffusion(int n, int epsilon_mlog, double convection_b) {

  /* Compute epsilon from epsilon_mlog where mlog is intended to mean
   * the negative logarithm base two.
   * We multiply by n * n here because the entries of our convection diffusion
   * test matrix are not divided by h * h where h = 1 / n is the distance
   * between nodes.
   */
  double epsilon = n * n * pow(0.5, epsilon_mlog);

  /* Create a new plot. */
  conv_plot_t *conv_plot = conv_plot_new();
  conv_plot_xlabel(conv_plot, "Iteration step");
  conv_plot_ylabel(conv_plot, "Residual norm");
  conv_plot_titlef(conv_plot,
     "Convection-diffusion matrix with n = %d, b = (%.1f, %.1f), epsilon = 0.5 ** %d",
     n, convection_b, convection_b, epsilon_mlog);

  /* Some constants. */
  int nsq = n * n;
  int maxit = 100;
  double tol = 0.0; /* Zero tolerance since we want maxit iterations to be performed. */
  int stopcrit = ITSOLVER_RELRES;

  /* Test matrix in two dimensions. */
  csr *A = csr_const_convection_diffusion(n, 2, epsilon, convection_b);

  /* Prescribed solution consisting of all ones. */
  double *xstar = malloc(nsq * sizeof(double));
  int i;
  for (i = 0; i < nsq; ++i) {
    xstar[i] = 1.0;
  }

  /* Corresponding right hand side. */
  double *b = malloc(nsq * sizeof(double));
  for (i = 0; i < nsq; ++i) {
    b[i] = 0.0;
  }
  csr_mv(A, 1.0, xstar, 0.0, b);

  /* All zero start solution. */
  double *x0 = malloc(nsq * sizeof(double));
  for (i = 0; i < nsq; ++i) {
    x0[i] = 0.0;
  }

  /* Iteration vector for every method. */
  double *x = malloc(nsq * sizeof(double));

  { /* Jacobi. */
    vec_copy(nsq, x0, x);
    int iter;
    double esterr;
    double *errvec = malloc((maxit + 1) * sizeof(double));
    csr_jacobi(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "Jacobi");
    free(errvec);
  }

  { /* Gauss-Seidel (vorwaerts). */
    vec_copy(nsq, x0, x);
    int iter;
    double esterr;
    double *errvec = malloc((maxit + 1) * sizeof(double));
    csr_gauss_seidel(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "GS");
    free(errvec);
  }

  { /* Gauss-Seidel (rueckwaerts). */
    vec_copy(nsq, x0, x);
    int iter;
    double esterr;
    double *errvec = malloc((maxit + 1) * sizeof(double));
    csr_gauss_seidel_reverse(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "bGS");
    free(errvec);
  }

  { /* SGS. */
    vec_copy(nsq, x0, x);
    int iter;
    double esterr;
    double *errvec = malloc((maxit + 1) * sizeof(double));
    csr_sgs(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "SGS");
    free(errvec);
  }

  /* Save plot to file and delete plot. */
  char *file_name = malloc(100);
  snprintf(file_name, 100, "testprog_convection_diffusion_n%d_b%.1f_%.1f_epsilon_mlog%d",
    n, convection_b, convection_b, epsilon_mlog);
  conv_plot_save(conv_plot, file_name);
  free(file_name);
  conv_plot_delete(conv_plot);

  /* Free memory. */
  free(x);
  free(x0);
  free(b);
  free(xstar);
  csr_free(A);
}

/* Testprogram investigating the speed of convergence for the
 * convection diffusion matrix for different values of n and epsilon.
 */
int main() {

   { /* Some convection and different values of epsilon. */
    int n = 10;
    double b_const = 0.0;
    int epsilon_mlog;
    for (epsilon_mlog = 0; epsilon_mlog <= 8; ++epsilon_mlog) {
      test_convection_diffusion(n, epsilon_mlog, b_const);
    }
  }
  { /* Some convection and different values of epsilon. */
    int n = 20;
    double b_const = 1.0;
    int epsilon_mlog;
    for (epsilon_mlog = 0; epsilon_mlog <= 8; ++epsilon_mlog) {
      test_convection_diffusion(n, epsilon_mlog, b_const);
    }
  }

  return 0;
}

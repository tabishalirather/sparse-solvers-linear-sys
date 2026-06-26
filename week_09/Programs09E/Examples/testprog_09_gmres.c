/*! \file
 *  \brief Comparison of restarted GMRes for a symmetric postive definite
 *         matrix.
 *
 *  Comparison of full GMRes without restart and restarted GMRes(m) with restart
 *  parameters m = 1, 2, ..., 8. The Minimal Residual method is included in the
 *  comparision since the results are expected to be the same as for GMRes(1).
 *  Note that the Minimal Residual Method requires a symmetric positive
 *  definite matrix, which is why we use the convection diffusion matrix
 *  with zero diffusion.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

int main(int argc, char **argv)
{

  /* Number of intervals. */
  int n = 20;
  
  /* Spatial dimension. */
  int d = 3;
  printf("n = %d, d = %d\n", n, d);

  /* convection diffusion model problem. */
  double epsilon = 1.0;
  double convection = 0.0 / n;
  csr *A = csr_const_convection_diffusion(n, d, epsilon, convection);

  int nd = A->nr;

  /* Constant right hand side. */
  double *b = vec_ones(nd);

  /* Maximum number of iterations and tolerance. */
  int iter_max = 80 * n;
  double tol = 1e-8;

  /* Create plot. */
  conv_plot_t *conv_plot = conv_plot_new();

  { /* Minimal Residual. */
    printf("Minimal Residual\n");
    double *x = vec_zeros(nd);
    double esterr;
    int iter;
    double *errvec = vec_zeros(iter_max + 1);
    csr_minimal_residual(A, b, x,
                         tol, ITSOLVER_RELRES, iter_max,
                         &esterr, &iter, errvec);
    printf("  iter = %d\n", iter);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_name(conv_plot, "MR");
    free(x);
    free(errvec);
  }

  /* Halve iter_max to make sure that the curve corresponding to GMRes(0) is
   * not completely hidden by the MR curve in the convergence plot. */
  iter_max /= 2;

  /* Compare different values of the restart parameter where restart = 0 means
   * not to restart ("full" GMRes). */
  for (int restart = 8; restart >= 0; restart -= 1) {

    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 3          */
    /* Configure and run GMRES; and add the */
    /* convergence line to the plot.        */
    /*                                      */
    /* ************************************ */
    
  }

  /* Save and delete plot. */
  conv_plot_save(conv_plot, "testprog_09_gmres");
  conv_plot_delete(conv_plot);

  /* Free A, b. */
  free(b);
  csr_free(A);

  return 0;
}

/**
 * @file testprog_04_tridiagonal.c
 *
 * @author Jonas Grams
 *
 */

#include "itsolver.h"

int main()
{
    int n;
    double alpha1, alpha2, alpha3;

    n = 1000;
    alpha1 = 2.0;
    alpha2 = 3.0;
    alpha3 = 1.9;

    /* Create plot */
    conv_plot_t *conv_plot = conv_plot_new();
    conv_plot_title(conv_plot, "Richardson method applied to the tridiagonal matrix for different diagonal entries");
    conv_plot_xlabel(conv_plot, "Iteration steps");
    conv_plot_ylabel(conv_plot, "Absolute residual norm");

    /* Zero right hand side. */
    double *b = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
    {
        b[i] = 0.0;
    }

    /* All ones start vector. */
    double *x = (double *)malloc(n * sizeof(double));
    double *xstart = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
    {
        xstart[i] = 1.0;
    }

    /* Perform Richardson's iteration. */
    double tol = 1e-6;
    int stopcrit = ITSOLVER_RELRES;
    int maxit = 1000;

    { /* Laplace discretization (alpha = 2)*/
        csr *A = csr_tridiagonal(n, alpha1);

        // Relaxation parameter 1/lambda_max
        double omega = 1.0 / tridiagonal_eigval_max(n, alpha1);

        // Perform richardson iteration
        double *errvec = malloc((maxit + 1) * sizeof(double));
        vec_copy(n, xstart, x);
        double esterr;
        int iter;
        csr_richardson(A, b, x, tol, stopcrit, maxit, omega,
                       &esterr, &iter, errvec);

        // Plot residual norm
        conv_plot_line(conv_plot, iter, errvec);
        conv_plot_namef(conv_plot, "alpha = %2.1f", alpha1);
        free(errvec);
        csr_free(A);
    }

    { /* Positive definite case (alpha = 3)*/
        csr *A = csr_tridiagonal(n, alpha2);

        // Relaxation parameter 1/lambda_max
        double omega = 1.0 / tridiagonal_eigval_max(n, alpha2);

        // Perform richardson iteration
        double *errvec = malloc((maxit + 1) * sizeof(double));
        vec_copy(n, xstart, x);
        double esterr;
        int iter;
        csr_richardson(A, b, x, tol, stopcrit, maxit, omega,
                       &esterr, &iter, errvec);

        // Plot residual norm
        conv_plot_line(conv_plot, iter, errvec);
        conv_plot_namef(conv_plot, "alpha = %2.1f", alpha2);
        free(errvec);
        csr_free(A);
    }

    { /* Indefinite case (alpha = 2)*/
        csr *A = csr_tridiagonal(n, alpha3);

        // Relaxation parameter 1/lambda_max
        double omega = 1.0 / tridiagonal_eigval_max(n, alpha3);

        // Perform richardson iteration
        double *errvec = malloc((maxit + 1) * sizeof(double));
        vec_copy(n, xstart, x);
        double esterr;
        int iter;
        csr_richardson(A, b, x, tol, stopcrit, maxit, omega,
                       &esterr, &iter, errvec);

        // Plot residual norm
        conv_plot_line(conv_plot, iter, errvec);
        conv_plot_namef(conv_plot, "alpha = %2.1f", alpha3);
        free(errvec);
        csr_free(A);
    }

    /* Save plot to file. */
    conv_plot_save(conv_plot, "testprog_04_tridiagonal");
    conv_plot_delete(conv_plot);

    /* Free memory. */
    free(x);
    free(b);
    free(xstart);

    return 0;
}
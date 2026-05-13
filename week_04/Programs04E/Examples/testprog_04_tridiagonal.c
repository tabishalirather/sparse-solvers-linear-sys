/**
 * @file testprog_04_tridiagonal.c
 *
 * @author Jonas Grams
 *
 */

#include <iso646.h>

#include "itsolver.h"

void print_csr(csr* A) {
  for (int i = 0; i < A->nr; i++) {
    // puts("running inside the print for loop");
    for (int j = A->row[i]; j < A->row[i + 1]; j++) {
      // puts("running inside the j print for loop");
      printf("%f ", A->coeff[j]);
    }
    printf("\n");
  }
}
int main() {
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
  double *b = (double *) malloc(n * sizeof(double));
  for (int i = 0; i < n; ++i) {
    b[i] = 0.0;
  }

  /* All ones start vector. */
  double *x = (double *) malloc(n * sizeof(double));
  double *xstart = (double *) malloc(n * sizeof(double));
  for (int i = 0; i < n; ++i) {
    xstart[i] = 1.0;
  }

  /* Perform Richardson's iteration. */
  double tol = 1e-6;
  int stopcrit = ITSOLVER_RELRES;
  int maxit = 1000;

  {
    /* Laplace discretization (alpha = 2) */

    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(n, xstart, x);
    double esterr;
    int iter;

    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 1          */
    /*                                      */
    /* Set up the tridiagonal matrix with   */
    /* damping parameter alpha = 2.         */
    /*                                      */
    /* ************************************ */

    // TODO: Set up the Tridiagonal matrix.
    csr *A = csr_tridiagonal(n, alpha1);
    // TODO: Get first and last row out of the for loop and then we only  need one for loop.

    /*//in general we will have 2elements in row 1 and 2elements  in last row, then then for (n-2) rows we have 3 elements, or 3(n-2) + 4 = 3n + 2
    // csr* A = csr_alloc(nr, nc, nz);
    // printf("%f\n", A->coeff[4]);

    // int row[5] = {0, 2, 5, 8, 10};
    // int row[5] = {0};
    // A->row[0] = 0;
    // int col[10] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    // double coeff[10] = {2, 3, 1, 2, 3, 1, 2, 3, 1, 2};
    // double main = alpha2;
    // puts("populating the tridiagonal matrix");
    // for (int i=1; i<=nr; i++) {
    //     if (i==1 ^ i==nr)
    //     A->row[i] = A->row[i-1] + 2;
    //     else {
    //         A->row[i] = A->row[i-1] + 3;
    //     }
    // }
    //
    // for (int j=0; j<=nr; j++) {
    //     printf("%d ", A->row[j]);
    // }
    // puts("row matrix setup");
    // for (int i=0; i<=A->nr; i++) {
    //     // printf("running inside the for loop");
    //     for (int j=A->row[i]; j<A->row[i+1]; j++) {
    //         double lower = -1;
    //         double upper = -1;
    //         // printf("running inside the j for loop");
    //
    //         if (i==0) {
    //             //this means first row, i.e. fill main and upper diagonal.
    //             if (j==A->row[i]) { //then we are filling the first element
    //                 A->coeff[j] = main;
    //                 // The general pattern for filling columns in i-1, i, i+1
    //                 A->col[j] = i; // index of first nz col
    //             }else if (j==A->row[i]+1) { //then we are filling the second element
    //                 A->coeff[j] = upper;
    //                 A->col[j] = i+1; //index of second nz col.
    //             }
    //         }
    //         else if (i+1 == A->nr){
    //             //this means last row, i.e. fill lower and main diagonal.
    //             if (j==A->row[i]) { //then we are filling the first element
    //                 A->coeff[j] = lower;
    //                 A->col[j] = i-1;
    //             }else if (j==A->row[i]+1) { //then we are filling the second element
    //                 A->coeff[j] = main;
    //                 A->col[j] = i;
    //             }
    //         } else{
    //             //fill lower, then main and then upper diagonal.
    //             if (j==A->row[i]) { //then we are filling the first element
    //                 A->coeff[j] = lower;
    //                 A->col[j] = i-1;
    //             }else if (j==A->row[i]+1) { //then we are filling the second element
    //                 A->coeff[j] = main;
    //                 A->col[j] = i;
    //             }else if (j==A->row[i]+2) { //then we are filling the second element
    //                 A->coeff[j] = upper;
    //                 A->col[j] = i+1;
    //             }
    //         }
    //     }
    // }*/

    // puts("\n Printing tridiag matrix");

    // TODO: Set the damping parameter.
    double damping_param = 1 / tridiagonal_eigval_max(n, alpha1);


    // TODO: Start the Richardson method.
    csr_richardson(A,
                   b,
                   x,
                   tol,
                   stopcrit,
                   maxit,
                   damping_param,
                   &esterr,
                   &iter,
                   errvec);
    puts("For alpha1 = 2");
    print_csr(A);
    // Plot residual norm
    conv_plot_line(conv_plot, iter, errvec);
    puts("after plot line");
    conv_plot_namef(conv_plot, "alpha = %2.1f", alpha1);
    puts("after plot namef line");

    free(errvec);
    csr_free(A);
  }

  {
    /* Positive definite case (alpha = 3)*/

    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(n, xstart, x);
    double esterr;
    int iter;
    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 1          */
    /*                                      */
    /* Set up the tridiagonal matrix with   */
    /* damping parameter alpha = 3.         */
    /*                                      */
    /* ************************************ */

    // TODO: Set up the Tridiagonal matrix.
    // csr* A = 0;
    csr *A = csr_tridiagonal(n, alpha2);
    // TODO: Set the damping parameter.
    double damping_param = 1/tridiagonal_eigval_max(n, alpha2);


    // TODO: Start the Richardson method.
    csr_richardson(A,
                   b,
                   x,
                   tol,
                   stopcrit,
                   maxit,
                   damping_param,
                   &esterr,
                   &iter,
                   errvec);


    // Plot residual norm
    puts("For alpha2 = 3");
    print_csr(A);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "alpha = %2.1f", alpha2);
    free(errvec);
    csr_free(A);
  }

  {
    /* Indefinite case (alpha = 1.9)*/

    // Perform richardson iteration
    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(n, xstart, x);
    double esterr;
    int iter;

    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 1          */
    /*                                      */
    /* Set up the tridiagonal matrix with   */
    /* damping parameter alpha = 1.9.       */
    /*                                      */
    /* ************************************ */

    // TODO: Set up the Tridiagonal matrix.
    csr *A = csr_tridiagonal(n, alpha3);

    // TODO: Set the damping parameter.
    double damping_param = 1/tridiagonal_eigval_max(n, alpha3);

    // TODO: Start the Richardson method.
    csr_richardson(A,
                   b,
                   x,
                   tol,
                   stopcrit,
                   maxit,
                   damping_param,
                   &esterr,
                   &iter,
                   errvec);
    // Plot residual norm
    puts("For alpha3 = 1.9");
    print_csr(A);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "alpha = %2.1f", alpha3);
    free(errvec);
    csr_free(A);
  }

  {
    /* Positive definite case (alpha = 3)*/

    double *errvec = malloc((maxit + 1) * sizeof(double));
    vec_copy(n, xstart, x);
    double esterr;
    int iter;
    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 1          */
    /*                                      */
    /* Set up the tridiagonal matrix with   */
    /* damping parameter alpha = optiomal.: 1/(max_eigen + min_eigen)         */
    /*                                      */
    /* ************************************ */

    // TODO: Set up the Tridiagonal matrix.
    // csr* A = 0;
    csr *A = csr_tridiagonal(n, alpha1);
    // TODO: Set the damping parameter.
    double damping_param = 2/(tridiagonal_eigval_max(n, alpha1)+tridiagonal_eigval_min(n, alpha1));


    // TODO: Start the Richardson method.
    csr_richardson(A,
                   b,
                   x,
                   tol,
                   stopcrit,
                   maxit,
                   damping_param,
                   &esterr,
                   &iter,
                   errvec);


    // Plot residual norm
    puts("For alpha2 = 3");
    print_csr(A);
    conv_plot_line(conv_plot, iter, errvec);
    conv_plot_namef(conv_plot, "damping_param_optimal = %2.1f", damping_param);
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

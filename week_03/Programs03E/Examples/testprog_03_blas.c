/*! \file
 *
 *   \author Michael Wende 
 *   \author Jens-Peter Zemke 
 *   \author Sabine Le Borne
 *   \author Lina Fesefeldt 
 */

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <itsolver.h>
#include <cblas.h>

double cblas_vec_dot(int n, const double *x, const double *y) {
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /*                                      */
  /*  Implement a function that computes  */
  /*  the inner product of two vectors.   */
  /* If it doens't run, use cblas blas instead of openblas */
  /* ************************************ */
  // printf("my_vec_dot called\n");
  double inner_product_result = cblas_ddot(n, x, 1, y, 1);
  return inner_product_result;
}


void rand_vector(int size, double *x) {
  // printf("100*((double) (rand() / RAND_MAX): %f\n", 100*((double) rand() / RAND_MAX));
  for (int i = 0; i < size; i++) {
    // printf("rand() / RAND_MAX;%f ", 100*((double)rand() / RAND_MAX));
    x[i] = 10.0 * ((double) rand() / RAND_MAX);
    // printf("x[i]: %f", x[i]);
  }
}

void run_cblas_vector_mult(double *x, double *y, double *inner_product_result, int n, int product_count,
                           struct timeval tv1, struct timeval tv2, struct timezone tz) {
  gettimeofday(&tv1, &tz);
  // gettimeofday(&tv2, &tz);

  //This code generates inner product of n random vectors and stores that in inner_product_result
  for (int i = 0; i < product_count; i++) {
    rand_vector(n, x);
    rand_vector(n, y);
    inner_product_result[i] = cblas_vec_dot(n, x, y);
    // printf(("inner_product_result[%d]: %f\n"),i, inner_product_result);
  }

  for (int i = 0; i < product_count; i++) {
    // inner_product_result[i] = my_vec_dot(n, x, y);
    printf(("inner_product_result[%d]: %f\n"), i, inner_product_result[i]);
  }
}


void vector_mult_manual(int n, int product_count, double *x, double *y, double *inner_product_result) {
  for (int j = 0; j < product_count; j++) {
    rand_vector(n, x);
    rand_vector(n, y);
    inner_product_result[j] = 0.0;
    for (int i = 0; i < n; i++) {
      inner_product_result[j] += x[i] * y[i];
    }
  }
}

// void run_manual_mult(struct timeval tv3,struct timezone tz, int n, double *vector_A, double *vector_B, double *result) {
//   gettimeofday(&tv3, &tz);
//   vector_mult_manual(n, vector_A, vector_B, result);
//
// }


int main() {
  srand(time(NULL));

  struct timeval tv1, tv2;
  struct timezone tz;

  struct timeval tv3, tv4;
  struct timeval tz_1;

  int n = 10000000;
  int product_count = 100;

  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
  double *inner_product_result = malloc(product_count * sizeof(double));


  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /*                                      */
  /*  Implement a loop that computes 100  */
  /*  scalar products of random vectors   */
  /*  using my_vec_dot.                   */
  /*                                      */
  /* ************************************ */


  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /*                                      */
  /*  Implement a loop that uses BLAS to  */
  /*  compute 100 scalar products of      */
  /*  random vectors.                     */
  /*                                      */
  /* ************************************ */
  gettimeofday(&tv1, &tz);
  run_cblas_vector_mult(x, y, inner_product_result, n, product_count, tv1, tv2, tz);
  gettimeofday(&tv2, &tz);
  double time_scalar_prod_cblas = (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec - tv1.tv_usec) / 1000.0;
  printf("time_scalar_prod_cblas: %f\n", time_scalar_prod_cblas);
  gettimeofday(&tv3, &tv4);
  free(x);
  free(y);
  free(inner_product_result);
  x = NULL;
  y = NULL;
  inner_product_result = NULL;


  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  inner_product_result = malloc(product_count * sizeof(double));

  gettimeofday(&tv3, &tz_1);
  vector_mult_manual(n, product_count, x, y, inner_product_result);
  gettimeofday(&tv3, &tz_1);
  double time_scalar_prod_own = (tv4.tv_sec - tv3.tv_sec) * 1000.0 + (tv4.tv_usec - tv3.tv_usec) / 1000.0;
  printf("time_scalar_prod_own: %f\n", time_scalar_prod_own);

  free(x);
  free(y);
  free(inner_product_result);
  return 0;
}

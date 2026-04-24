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

void rand_vector(int size, double *x) {
  // printf("100*((double) (rand() / RAND_MAX): %f\n", 100*((double) rand() / RAND_MAX));
  for (int i = 0; i < size; i++) {
    // printf("rand() / RAND_MAX;%f ", 100*((double)rand() / RAND_MAX));
    x[i] = 10.0 * ((double) rand() / RAND_MAX);
    // printf("x[i]: %f", x[i]);
  }
}

double my_vec_dot(int n, const double *x, const double *y) {
  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /*                                      */
  /*  Implement a function that computes  */
  /*  the inner product of two vectors.   */
  /* If it doens't run, use cblas blas instead of openblas */
  /* ************************************ */
  printf("my_vec_dot called\n");
  double inner_product_result = cblas_ddot(n, x, 1, y, 1);
  // printf("inner_product_result: %f\n", inner_product_result);

  return inner_product_result;
}


int main() {
  srand(time(NULL));

  struct timeval tv1, tv2;
  struct timezone tz;

  int n = 10000000;
  int product_count = 100;

  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
  double *inner_product_result = malloc(product_count * sizeof(double));

  gettimeofday(&tv1, &tz);

  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /*                                      */
  /*  Implement a loop that computes 100  */
  /*  scalar products of random vectors   */
  /*  using my_vec_dot.                   */
  /*                                      */
  /* ************************************ */


  gettimeofday(&tv2, &tz);
  double time_scalar_prod_own = (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec - tv1.tv_usec) / 1000.0;

  gettimeofday(&tv1, &tz);

  /* ************************************ */
  /*                                      */
  /*          TODO --- Problem 2          */
  /*                                      */
  /*  Implement a loop that uses BLAS to  */
  /*  compute 100 scalar products of      */
  /*  random vectors.                     */
  /*                                      */
  /* ************************************ */
  for (int i=0; i<product_count; i++) {
    rand_vector(n, x);
    rand_vector(n, y);
    inner_product_result[i] = my_vec_dot(n, x, y);
    // printf(("inner_product_result[%d]: %f\n"),i, inner_product_result);
  }

  for (int i=0; i<product_count; i++) {
    // inner_product_result[i] = my_vec_dot(n, x, y);
    printf(("inner_product_result[%d]: %f\n"),i, inner_product_result[i]);
  }
  gettimeofday(&tv2, &tz);
  printf("time_scalar_prod_own: %f\n", time_scalar_prod_own);
  free(x);
  free(y);
  free(inner_product_result);

  return 0;
}

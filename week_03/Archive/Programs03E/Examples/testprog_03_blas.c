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



double my_vec_dot(int n, const double *x, const double *y) {

    /* ************************************ */
    /*                                      */
    /*          TODO --- Problem 2          */
    /*                                      */
    /*  Implement a function that computes  */
    /*  the inner product of two vectors.   */
    /*                                      */
    /* ************************************ */
    
}



int main() {
    struct timeval tv1, tv2;
    struct timezone tz;

    int n = 10000000;
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

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
    
    gettimeofday(&tv2, &tz);

    double time_scalar_prod_blas = (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec - tv1.tv_usec) / 1000.0;

    printf("Speedup: %lf\n", time_scalar_prod_blas/time_scalar_prod_own);

    free(x);
    free(y);

    return 0;
}


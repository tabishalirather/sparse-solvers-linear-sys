/*! \file
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */


#include "itsolver.h"

void csr_mv(csr *A, double alpha, const double *x, double beta, double *y)
{
  /* ******************************************************************** */
  /*                                                                      */
  /*     TODO --- Exercise 1, Problem 3                                   */
  /*                                                                      */
  /* ******************************************************************** */
    printf("\nPrinting in Examples/csr_mv exercise code");
    for (int i=0; i<A->nr; i++) {
        double sum = 0;
        for (int j=A->row[i]; j<A->row[i+1]; j++) {
            sum += A->coeff[j] * x[A->col[j]];
        }
        y[i] = alpha*sum + beta*y[i];
    }

    // printf("Results of matrix mult are: \n");
    // for (int i=0; i<A->nr; i++) {
    //     printf("%f\n", y[i]);
    // }

}

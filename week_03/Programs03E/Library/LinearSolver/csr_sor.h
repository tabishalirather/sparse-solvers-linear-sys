/*! \file
 *  \brief Successive overrelaxation for solving a linear system with a CSR
 *         matrix.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_SOR_H
#define CSR_SOR_H

void csr_sor(csr *A, double *b, double *x,
            double tol, int stopcrit, int maxit, double omega,
            double *esterr, int *iter, double *errvec);

#endif

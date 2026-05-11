/*! \file
 *  \brief Jacobi iteration for solving a linear system with
 *         a CSR matrix.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_JACOBI_H
#define CSR_JACOBI_H

void csr_jacobi(csr *A, double *b, double *x, 
		double tol, int stopcrit, int maxit,
		double *esterr, int *iter, double *errvec);

#endif

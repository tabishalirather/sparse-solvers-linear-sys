/*! \file
 *  \brief Reverse Gauss-Seidel iteration for solving a linear system with
 *         a CSR matrix.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_GAUSS_SEIDEL_REVERSE_H
#define CSR_GAUSS_SEIDEL_REVERSE_H

void csr_gauss_seidel_reverse(csr *A, double *b, double *x, 
			      double tol, int stopcrit, int maxit,
			      double *esterr, int *iter, double *errvec);
#endif

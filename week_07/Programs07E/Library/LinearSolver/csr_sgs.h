/*! \file
 *  \brief Symmetric Gauss-Seidel iteration for solving a linear system with a CSR  matrix.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_SGS_H
#define CSR_SGS_H

void csr_sgs(csr *A, double *b, double *x, 
	     double tol, int stopcrit, int maxit,
	     double *esterr, int *iter, double *errvec);

#endif

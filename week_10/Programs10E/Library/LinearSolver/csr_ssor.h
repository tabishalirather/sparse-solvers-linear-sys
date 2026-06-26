/*! \file
 *  \brief Symmetric successive overrelaxation (SSOR) for solving a linear
 *         system with a CSR matrix.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_SSOR_H
#define CSR_SSOR_H

void csr_ssor(csr *A, double *b, double *x, 
	      double tol, int stopcrit, int maxit, double omega,
	      double *esterr, int *iter, double *errvec);

#endif

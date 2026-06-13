/*! \file
 *  \brief Conjugate gradients squared (CGS).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_CGS_H
#define CSR_CGS_H

/*! \brief Implementation of the CGS algorithm.
 */
void csr_cgs(csr *A, double *b, double *x,
	     double tol, int stopcrit, int maxit,
	     double *esterr, int *iter, double *errvec);

#endif

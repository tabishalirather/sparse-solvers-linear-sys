/*! \file
 *  \brief Biconjugate gradients (BiCG).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_BICG_H
#define CSR_BICG_H

/*! \brief Implementation of the BiCG algorithm.
 */
void csr_bicg(csr *A, double *b, double *x,
	      double tol, int stopcrit, int maxit,
	      double *esterr, int *iter, double *errvec);

#endif

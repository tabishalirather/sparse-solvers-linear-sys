/*! \file
 *  \brief Stabilized biconjugate gradients (BiCGStab).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_BICGSTAB_H
#define CSR_BICGSTAB_H

/*! \brief Implementation of the BiCGStab algorithm.
 */
void csr_bicgstab(csr *A, double *b, double *x,
	          double tol, int stopcrit, int maxit,
	          double *esterr, int *iter, double *errvec);

#endif

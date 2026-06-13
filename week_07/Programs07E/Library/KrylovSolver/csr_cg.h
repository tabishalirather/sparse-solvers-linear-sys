/*! \file
 *  \brief Conjugate gradients (CG).
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_CG_H
#define CSR_CG_H

/*! \brief Implementation of CG algorithm.
 */
void csr_cg(csr *A, double *b, double *x,
            double tol, int stopcrit, int maxit,
            double *esterr, int *iter, double *errvec);

#endif

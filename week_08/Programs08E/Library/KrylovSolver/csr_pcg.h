/*! \file
 *  \brief Preconditioned CG iteration.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_PCG_H
#define CSR_PCG_H

/*! \brief Preconditioned conjugate grandient method.
 *
 *  \param A Symmetric positive definite CSR matrix.
 *  \param prec Symmetric postive definite preconditioner.
 *  \param b Right hand side.
 *  \param x Initial guess (input) and computed solution (output).
 *  \param tol Tolerance.
 */
void csr_pcg(csr *A, prec_t *prec, double *b, double *x,
             double tol, int stopcrit, int maxit,
             double *esterr, int *iter, double *errvec);

#endif

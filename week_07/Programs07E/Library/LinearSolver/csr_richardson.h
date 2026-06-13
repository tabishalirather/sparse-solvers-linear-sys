/*! \file
 *  \brief Richardson's iteration for solving a linear system with a CSR
 *         matrix.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_RICHARDSON_H
#define CSR_RICHARDSON_H

/*! Richardson iteration for solving a linear system A x = b with a positive definite matrix A.
 *  In any case not more than maxit iterations are performed.
 *  A lesser number of iterations may be performed if one of the stopping criteria is used.
 *
 *    \param A Symmetric positive definite matrix.
 *    \param b Right hand side.
 *    \param x Start solution.
 *    \param stopcrit Stopping criterion to be used.
 *                    For \c ITSOLVER_NORES no residual is computed and maxit iterations are performed in any case.
 *                    With \c ITSOLVER_RELRES the current residual is compared to the initial residual times the tolerance.
 *    \param tol Tolerance to be used with the selected stopping criterion.
 *    \param maxit Upper bound on the number of iterations to be performed.
 *    \param omega Parameter of the Richardson iteration.
 *    \param esterr Estimated error (residual norm) after the last step.
 *    \param iter Number of steps that where actually performed.
 *    \param errvec Estimated error (residual norm) initially and after each step.
 */

void csr_richardson(csr *A, double *b, double *x, 
		    double tol, int stopcrit, int maxit,
		    double omega,
		    double *esterr, int *iter, double *errvec);

#endif

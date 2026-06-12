/*! \file
 *  \brief One dimensional projection methods.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_PROJECTION_H
#define CSR_PROJECTION_H

void csr_steepest_descent(csr *A, double *b, double *x,
                          double tol, int stopcrit, int maxit,
                          double *esterr, int *iter, double *errvec);

void csr_minimal_residual(csr *A, double *b, double *x,
                          double tol, int stopcrit, int maxit,
                          double *esterr, int *iter, double *errvec);

void csr_residual_norm_steepest_descent(csr *A, double *b, double *x,
                                        double tol, int stopcrit, int maxit,
                                        double *esterr, int *iter, double *errvec);

#endif

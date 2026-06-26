/*! \file
 *  \brief Some arithmetic operations involving CSR matrices.
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_ARITH_H
#define CSR_ARITH_H

#include "csr.h"

/*! \brief Perform two-sided matrix-vector multiplication.
 *
 *    \param A Matrix in CSR format.
 *    \param x Vector that will be multiplied from the left.
 *    \param y Vector that will be multiplied from the right.
 *    \return Two-sided product \f$ x^T A y \f$.
 */
double csr_xtAy(csr *A, double *x, double *y);

/*! \brief Multiplication of vector by a transposed matrix.
 *
 *  For a given matrix \f$ A \f$ a multiplication with the transposed matrix will
 *  be performed without actually forming and storing \f$ A^T \f$.
 *  The vector \f$ y \f$ will be updated according to \f$ y = \alpha A^T x + \beta y \f$. 
 *
 *    \param A Matrix with \f$ m \f$ rows and \f$ n \f$ columns in CSR-Format.
 *    \param alpha Scalar multiplier.
 *    \param x Vector with \f$ m \f$ entries.
 *    \param beta Scalar multiplier.
 *    \param y Vector with \f$ n \f$ entries.
 */
void csr_mv_transposed(csr *A, double alpha, double *x, double beta, double *y);

/*! \brief Compute the residual \f$ r = b - A x \f$.
 *
 *    \param A Matrix.
 *    \param b Right hand side.
 *    \param x Approximate solution.
 *    \param r Residual.
 */
void csr_residual(csr *A, double *b, double *x, double *r);

/*! \brief Compute the Euclidean norm of the residual \f$ r = b - A x \f$.
 *
 *    \param A Matrix.
 *    \param b Right hand side.
 *    \param x Approximate solution.
 *    \return Euclidean norm of the residual.
 */
double csr_residual_norm(csr *A, double *b, double *x);

/*! \brief Compute the residual \f$ r = b - A^T x \f$ using the transposed matrix.
 *
 *    \param A Matrix to be transposed.
 *    \param b Right hand side.
 *    \param x Approximate solution.
 *    \param r Residual.
 */
void csr_residual_transposed(csr *A, double *b, double *x, double *r);

/*! \brief CSR multiplication.
 *
 *    \param A
 *    \param B
 *    \return Product \f$ C = A * B \f$.
 */
csr *csr_multiply(csr *A, csr *B);

/*! \brief Add a multiple of an identity matrix to a CSR matrix.
 *    \param A Square matrix in CSR format.
 *    \param alpha Factor by which the identity matrix that is added to given
 *                 matrix is scaled.
 */
void csr_diag_comp(csr *A, double alpha);

/*! \brief Produce a (newly allocated) matrix \f$ C = A + B \f$.
 *
 *    \param a Matrix \f$ A \f$.
 *    \param b Matrix \f$ B \f$, possibly with non-zero pattern different from \f$ A \f$.
 *    \return Sum \f$ C = A + B \f$ with a non-zero pattern corresponding to the union
 *            of the non-zero patterns of \f$ A \f$ and \f$ B \f$.
 */
csr *csr_add(csr *a, csr *b);

/*! \brief Produce a (newly allocated) matrix \f$ C = A - B \f$.
 *
 *    \param a Matrix \f$ A \f$.
 *    \param b Matrix \f$ B \f$, possibly with non-zero pattern different from \f$ A \f$.
 *    \return Difference \f$ C = A - B \f$ with a non-zero pattern corresponding to the union
 *            of the non-zero patterns of \f$ A \f$ and \f$ B \f$.
 */
csr *csr_substract(csr *a, csr *b);

/*! \brief Scale CSR matrix by a given factor.
 *
 *    \param A Matrix to be scaled.
 *    \param alpha Factor by which to scale.
 */
void csr_scale(csr *a, double alpha);

#endif

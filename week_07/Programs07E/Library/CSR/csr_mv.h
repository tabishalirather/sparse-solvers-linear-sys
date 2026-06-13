/*! \file
 *  \brief Multiplication of a matrix in CSR format by a vector from the right.
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_MV_H
#define CSR_MV_H

/*! \brief Multiplication of a matrix in CSR format by a vector from the right.
 *
 *  The vector y is updated according to \f$ y \gets \alpha A x + \beta y \f$ .
 *
 *    \param A Matrix \f$ A \in \mathbb{R}^{m \times n} \f$.
 *    \param alpha Scalar.
 *    \param x Vector \f$ x \in \mathbb{R}^{n} \f$.
 *    \param beta Scalar.
 *    \param y Vector \f$ y \in \mathbb{R}^{m} \f$
 */
void csr_mv(csr *A, double alpha, const double *x, double beta, double *y);

#endif



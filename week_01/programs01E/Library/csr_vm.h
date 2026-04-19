/*! \file
 *  \brief Multiplication of a matrix in CSR format by a vector from the left.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_VM_H
#define CSR_VM_H

/*! \brief Multiplication of a matrix in CSR format by a vector from the left.
 *
 *  The vector y is updated according to \f$ y^T \gets \alpha x^T A + \beta y^T \f$ .
 *
 *    \param A Matrix \f$ A \in \mathbb{R}^{m \times n} \f$.
 *    \param alpha Scalar.
 *    \param x Vector \f$ x \in \mathbb{R}^{m} \f$.
 *    \param beta Scalar.
 *    \param y Vector \f$ y \in \mathbb{R}^{n} \f$
 */
void csr_vm(csr *A, double alpha, double *x, double beta, double *y);

#endif

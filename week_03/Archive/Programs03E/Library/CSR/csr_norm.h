/*! \file
 *  \brief Infinity and Frobenius norms of CSR matrices.
 *  \author Sabine Le Borne
 */

#ifndef CSR_NORM_H
#define CSR_NORM_H

#include "csr.h"

/*! \brief Compute the infinity norm of a CSR matrix.
 *    \param A CSR matrix \f$ A \f$.
 *    \return Frobenius norm \f$ \| A \|_{\infty} \f$.
 */
double csr_norm_infty(const csr *A);

/*! \brief Compute the Frobenius norm of a CSR matrix.
 *    \param A CSR matrix \f$ A \f$.
 *    \return Frobenius norm \f$ \| A \|_F \f$.
 */
double csr_norm_frob(const csr *A);

#endif

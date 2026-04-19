/*! \file
 *  \brief Compute the transpose of a sparse matrix.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_TRANSPOSE_H
#define CSR_TRANSPOSE_H

#include "csr.h"

/*! \brief Compute the transpose of a sparse matrix with a
 *         complexity proportional to the number of non-zeros.
 *
 *    \param A CSR matrix.
 *    \return Transpose of A.
 */
csr *csr_transpose(csr *A);

#endif



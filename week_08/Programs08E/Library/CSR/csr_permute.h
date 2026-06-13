/*! \file
 *  \brief Permute rows and columns of a matrix in CSR format.
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_PERMUTE_H
#define CSR_PERMUTE_H

#include "csr.h"
#include "Permutations/perm.h"

/*! \brief Permute rows and columns of a matrix in CSR format.
 *
 *  The input matrix A will be permuted to produce B = P A Q where
 *  P and Q are permutation matrices.
 *  Let A have m rows and n columns.
 *  The permutation matrix P is represented by a vector
 *    p = (p_1, ..., p_m)
 *  of length m from which P is defined as
 *    P = (e_{p_1}, ..., e_{p_m})
 *  where e_{i} indicates the i-th unit vector of length m.
 *  The permutation matrix Q is represented by a vector q of length n
 *  in the same way. While the vector q is given directly to this function,
 *  permuting the rows of A is easier to implement when the inverse
 *  permutation of P is known. Therefore, this function expects pinv to be
 *  given instead of p where pinv is the vector representation of the inverse
 *  permutation of p.
 *
 *    \param A Matrix in CSR format, possibly rectangular.
 *    \param pinv Vector representation of the inverse of the permutation applied to rows.
 *    \param q Vector representation of permutation applied to columns.
 *
 *    \return Permuted matrix.
 */
csr *csr_permute(csr *A, const int *pinv, const int *q);

csr *csr_perm(csr *A, perm_t *row_perm, perm_t *col_perm);

csr *csr_perm_symm(csr *A, perm_t *perm);

#endif

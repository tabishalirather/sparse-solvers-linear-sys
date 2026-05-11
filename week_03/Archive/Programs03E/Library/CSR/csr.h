/*! \file
 *  \brief Matrix in compressed sparse row (CSR) storage.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_H
#define CSR_H

#include <stdlib.h> /* size_t */

/*! \brief Structure representing a CSR matrix.
 */
typedef struct csr_matrix
{
  int nr;          /*!< Number of rows. */
  int nc;          /*!< Number of columns. */
  int nz;          /*!< Number of non-zero entries. */
  int *row;        /*!< Pointers to the beginning of each row. */
  int *col;        /*!< Column indices of non-zero coefficients. */
  double *coeff;   /*!< Non-zero coefficient values. */
} csr;

/*! \brief Allocate storage for a CSR matrix.
 *
 *    \param rows Number of rows.
 *    \param cols Number of columns.
 *    \param elements Number of non-zeros.
 */
csr *csr_alloc(int rows, int cols, int elements);

/*! \brief Free storage associated with a CSR matrix.
 *
 *    \param A Matrix in CSR format.
 */
void csr_free(csr *A);

/*! \brief Produce a new matrix that is a copy of the given matrix.
 *
 *    \param A Matrix to be copied.
 *    \return Newly allocated copy.
 */
csr *csr_copy(csr *A);

/*! \brief Produce a random CSR matrix.
 *
 *  The number of non-zero elements per row is random but bounded
 *  by a constant.
 *  The coefficient values are small integers.
 *
 *   \param nrows Number of rows.
 *   \param ncols Number of cols.
 *   \return Random matrix.
 */
csr *csr_rand(int nrows, int ncols);

/*! \brief Compute the bandwith of a sparse matrix.
 *
 *  In this case the bandwidth is defined to be the largest value of \f$ | i - j | \f$
 *  for any non-zero entry \f$ a_{i,j} \ne 0 \f$.
 *
 *    \param A CSR Matrix.
 *    \return Bandwith of A.
 */
int csr_bandwidth(csr *A);

/*! \brief Compute the storage required for a matrix in CSR format.
 *
 *  The storage comprises the structure for A itself as well as the
 *  referenced arrays for storing coefficients, column indices and pointers.
 *
 *    \param A CSR Matrix.
 *    \return Total memory consumption in bytes.
 */
size_t csr_bytes(csr *A);

/*! \brief Sort entries in each row by column indices.
 *
 *    \param A CSR Matrix.
 */
void csr_sort_col(csr *A);

/*! \brief Order diagonal entry in each row first.
 *
 *    \param A CSR Matrix.
 */
void csr_sort_diagFirst(csr *A);


#endif

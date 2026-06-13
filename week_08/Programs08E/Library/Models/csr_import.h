/*! \file
 *  \brief Functions for importing sparse matrices and right hand sides
 *         from the Matrix Market repository.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_IMPORT_H
#define CSR_IMPORT_H

/*! \brief Import sparse matrix from a file.
 *
 *  \param file_name File name with suffix .mtx holding a sparse matrix
 *                   from the Matrix Market repository.
 *  \return Imported matrix in CSR format.
 */
csr *csr_import_from_mm_mtx(const char *file_name);

/*! \brief Import right hand side from a file.
 *
 *  \param file_name File name with suffix .mtx holding a (dense) right hand
 *                   side from the Matrix Market repository.
 *  \return Imported right hand side vector.
 */
double *csr_import_from_mm_rhs(const char *file_name);

#endif

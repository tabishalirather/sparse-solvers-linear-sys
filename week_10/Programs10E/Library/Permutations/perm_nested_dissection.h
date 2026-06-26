/*! \file
 *  \brief Nested dissection permutation.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef PERM_NESTED_DISSECTION_H
#define PERM_NESTED_DISSECTION_H

/*! \brief Nested dissection node reordering.
 *
 *  \param A Square CSR matrix.
 *  \return Vector respresentation of nested dissection permutation.
 */
perm_t *perm_nested_dissection(csr *A);

#endif

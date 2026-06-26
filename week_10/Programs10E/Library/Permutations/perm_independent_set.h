/*! \file
 *  \brief Independent set permutation.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef PERM_INDEPENDENT_SET_H
#define PERM_INDEPENDENT_SET_H

/*! \brief Independent set node reordering.
 *
 *  \param A Square CSR matrix.
 *  \return Vector respresentation of independent set permutation.
 */
perm_t *perm_independent_set(csr *A);

#endif

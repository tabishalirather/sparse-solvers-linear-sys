/*! \file
 *  \brief Level set permutations BFS, CMK and RCM.
 *
 *  \author Sabine Le Borne
 */

#ifndef PERM_LEVEL_SET_H
#define PERM_LEVEL_SET_H

enum LEVEL_SET_VARIANT {
  LEVEL_SET_VARIANT_BFS, /* Breadth first search.  */
  LEVEL_SET_VARIANT_CMK, /* Cuthill-McKee.         */
  LEVEL_SET_VARIANT_RCM  /* Reverse Cuthill-McKee. */
};

/*! \brief Function that can be used to for different variants of
 *         level set reorderings, especially BFS and CMK.
 *
 *  \param A CSR matrix.
 *  \param start Index of a start node.
 *  \param level_set_variant Variant of the level set algorithm.
 *  \return Vector representation of the computed level set permutation.
 */
perm_t *perm_level_set(csr *A, int start, enum LEVEL_SET_VARIANT level_set_variant);

/*! \brief Breadth first search (BFS) permutation of a given CSR matrix.
 *
 *  \param A Square CSR matrix.
 *  \param start Start node or -1 to let the implementation choose.
 *  \return Vector representation of BFS permutation.
 */
perm_t *perm_bfs(csr * A, int start);

/*! \brief Cuthill McKee (CMK) permutation of a given CSR matrix.
 *
 *  \param A Square CSR matrix.
 *  \param start Start node or -1 to let the implementation choose.
 *  \return Vector representation of the CMK permutation, which is the
 *          reverse (not the inverse) of the RCM permutation.
 */
perm_t *perm_cmk(csr * A, int start);

/*! \brief Reverse Cuthill-McKee (RCM) permutation of a given CSR matrix.
 *
 *  \param A Square CSR matrix.
 *  \param start Start node or -1 to let the implementation choose.
 *  \return Vector representation of the RCM permutation, which is the
 *          reverse (not the inverse) of the CMK permutation.
 */
perm_t *perm_rcm(csr * A, int start);

#endif

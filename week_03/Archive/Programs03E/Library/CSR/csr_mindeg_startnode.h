/*! \file
 *  \brief Compute node to be used as a starting point for reordering
 *         algorithms.
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_MINDEG_STARTNODE_H
#define CSR_MINDEG_STARTNODE_H

/*! \brief Compute start node to be used for reordering matrix indices.
 *
 *  The computed start node will be of minimum degree among all nodes in
 *  the connectivity graph of the given matrix.
 *
 *    \param A CSR matrix, possibly non-square.
 *
 *    \return Index of start node.
 */
int csr_mindeg_startnode(csr *A);

#endif

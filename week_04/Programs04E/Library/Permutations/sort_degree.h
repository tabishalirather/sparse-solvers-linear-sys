/*! \file
 *  \brief Sort nodes according to their degree.
 *
 *  \author Jens-Peter Zemke
 */

#ifndef SORT_DEGREE_H
#define SORT_DEGREE_H

/*! \brief Sort a vector of indices where each index represents a node in a
 *         matrix graph. A second vector specifies the degree of each node
 *         in the graph. Both vectors will be sorted in ascending order with
 *         respect to degree.
 *         .
 *    \param v Vector of node indices to be sorted.
 *    \param deg Vector of node degrees which will be sorted accordingly.
 *    \param n Length of both vectors.
 */
void sort_degree(int *v, int *deg, int n);

#endif

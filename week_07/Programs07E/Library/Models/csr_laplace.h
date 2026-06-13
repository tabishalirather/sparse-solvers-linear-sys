/*! \file
 *  \brief Discretisation of the (negative) Laplacian operator on the unit cube in arbitrary dimensions.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_LAPLACE_H
#define CSR_LAPLACE_H

#include "CSR/csr.h"

/*! \brief Discretisation of the (negative) Laplacian operator on the unit cube in arbitrary dimensions.
 *
 *  The returned matrix \f$ A_{n, d} \in \mathbb{R}^{m \times m} \f$ consists of \f$ m = (n - 1)^d \f$
 *  rows and columns.
 *
 *    \param n Number of intervals in each dimension.
 *    \param d Space dimension.
 *    \return The computed matrix in CSR format.
 */
csr *csr_laplace(int n, int d);

/*! \brief Discretisation of the (negative) Laplacian operator on the unit interval.
 *
 *  The matrix returned by csr_laplace1d(n) should be the same as for csr_laplace(n, 1).
 *
 *    \param n Number of intervals in each dimension.
 *    \return Matrix with \f$ n - 1 \f$ rows and columns.
 */
csr *csr_laplace1d(int n);

/*! \brief Compute the eigenvector of \f$ A_{n,d} \f$
 *         corresponding to indices \f$ (i_1, \dotsc, i_d) \f$.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \param i Multiindex to select the eigenvector.
 *    \param v Preallocated storage for the computed eigenvector.
 */
void laplace_ev(int n, int d, int *i, double *v);

/*! \brief Compute the eigenvector of \f$ A_{n,d} \f$
 *         corresponding to the smallest eigenvalue.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \param v Preallocated storage for the computed eigenvector.
 */
void laplace_ev_min(int n, int d, double *v);

/*! \brief Compute the eigenvector of \f$ A_{n,d} \f$
 *         corresponding to the largest eigenvalue.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \param v Preallocated storage for the computed eigenvector.
 */
void laplace_ev_max(int n, int d, double *v);

/*! \brief Compute the eigenvalue of the matrix \f$ A_{n,d} \f$ that
 *  corresponds to the indices \f$ (i_1, \dotsc, i_d) \f$.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \param i Multiindex to select the eigenvalue.
 *    \return Computed eigenvalue.
 */
double laplace_ew(int n, int d, int *i);

/*! \brief Compute the minimal eigenvalue of \f$ A_{n, d} \f$.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \return Smallest eigenvalue.
 */
double laplace_ew_min(int n, int d);

/*! \brief Compute the maximal eigenvalue of \f$ A_{n, d} \f$.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \return Largest eigenvalue.
 */
double laplace_ew_max(int n, int d);


/*! \brief Compute the condition number of \f$ A_{n, d} \f$.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 *    \return Condition number \f$ \kappa = \lambda_{max} / \lambda_{min} \f$.
 */
double laplace_cond(int n, int d);

/*! \brief Compute th number of non-zero entries in \f$ A_{n,d} \f$.
 *
 *    \param n Number of intervals in one dimension.
 *    \param d Dimension.
 */
int laplace_nnz(int n, int d);

#endif


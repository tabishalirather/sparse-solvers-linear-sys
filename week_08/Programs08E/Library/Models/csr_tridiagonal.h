/**
 * \file csr_tridiagonal.c
 * \brief Tridiagonal test matrix with -1 on the off-diagonals, and some real value \f$\alpha\f$
 *        on the main diagonal.
 *
 * \author Jonas Grams
 *
 */

/**
 * \brief Tridiagonal test matrix \f$A_{\alpha}\in\mathbb{R}^{n\times n\f$
 *        with -1 on the off-diagonals, and some real value \f$\alpha\f$
 *        on the main diagonal with.
 * 
 * @param n Matrix dimension.
 * @param alpha Diagonal value. For \f$\alpha=1\f$, the discretization
 *  of the (negative) 1D Laplace operator is returned.
 * @return csr* Matrix in the CSR format.
 */
csr *csr_tridiagonal(int n, double alpha);

/**
 * \brief Compute the minimal eigenvalue of the tridiagonal matrix 
 *  created by csr_tridiagonal.
 * 
 * @param n Matrix dimension.
 * @param alpha Diagonal entry.
 * @return double Minimal eigenvalue of the \f$n\times n\f$ tridiagonal matrix 
 *  with diagonal entry alpha, created by csr_tridiagonal
 *  
 */
double tridiagonal_eigval_min(int n, double alpha);

/**
 * \brief Compute the maximal eigenvalue of the tridiagonal matrix 
 *  created by csr_tridiagonal.
 * 
 * @param n Matrix dimension.
 * @param alpha Diagonal entry.
 * @return double Maximal eigenvalue of the \f$n\times n\f$ tridiagonal matrix 
 *  with diagonal entry alpha, created by csr_tridiagonal
 *  
 */
double tridiagonal_eigval_max(int n, double alpha);
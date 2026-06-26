/*! \file
 *  \brief Incomplete LU factorisation.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef ILU0_H
#define ILU0_H

/*! \brief Compute a single row in the ILU(0) factorisation of a CSR matrix
 *         assuming all previous rows have already been computed.
 *
 *  \param A Matrix to be overwritten with its ILU factorisation.
 *  \param i Row to be computed.
 */
void ilu0_inplace_csr_row(csr *A, int i);

/*! \brief Incomplete LU factorisations without fill-in (ILU(0)).
 *
 *  Can we have an rectangular input matrix?
 *  Maybe check whether the element we devide by is zero or small.
 *
 *    \param A Matrix $\f A \f$ on input, ILU(0) factors \f$ L, U \f$ on output.
 *             The strictly lower part will contain \f$ L \f$ with its unit diagonal stripped.
 *             The upper part will contain \f$ U \f$.
 */
void ilu0_inplace_csr(csr *A);

/*! \brief Compute and return ILU(0) factorisation of a CSR matrix.
 *
 *  \param A Matrix to be factorised which will not be altered.
 *  \return Newly allocated matrix containing the computed factors L and U.
 */
csr *ilu0_csr(csr *A);

/*! \brief Compute the Frobenius norm in that part of the residual R = A - L * U
 *         where the coefficients of A are non-zero.
 *
 *  Since only the pattern of A is considered, this is not the norm of R.
 *  In absence of rounding errors, the return value should be zero if L and U come
 *  from the ILU(0) decomposition of A.
 *  The function is used to check that the implementation of ILU(0) works.
 *
 *    \param A
 *    \param LU Combination of L in the strictly lower and U in the upper part.
 *    \return Mean square error.
 */
double ilu0_check_csr(csr *A, csr *LU);

#endif

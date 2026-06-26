/*! \file
 *  \brief Operations actions the upper/ lower part of a CSR matrix.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_TRIANGULAR_H
#define CSR_TRIANGULAR_H

/*! \brief Extract the strict lower part of a matrix and add a unit diagonal.
 *    \param LU
 *    \return Newly allocated matrix.
 */
csr *csr_extract_unit_lower(csr *LU);

/*! \brief Extract the upper part of a matrix.
 *    \param LU
 *    \return Newly allocated matrix.
 */
csr *csr_extract_upper(csr *LU);

/*! \brief Perform forward substitution using the lower part of a matrix
 *         and a unit diagonal.
 *
 *    \param L Matrix of which only the strict lower part is used.
 *    \param x Right hand side on input, solution on output.
 */
void csr_lowersolve(csr *L, double *x);

/*! \brief Perform backward substitution using the upper part of a matrix
 *         (with non-zero diagonal elements, of course).
 *
 *    \param U Matrix of which only the upper part is used.
 *    \param x Right hand side on input, solution on output.
 */
void csr_uppersolve(csr *U, double *x);

#endif

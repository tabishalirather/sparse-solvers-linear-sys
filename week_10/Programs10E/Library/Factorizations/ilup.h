/*! \file
 *  \brief Incomplete LU factorisation with level of fill.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef ILUP_H
#define ILUP_H

#include "levmat.h"

/*! \brief Compute a single row of an ILU(p) factorisation.
 *
 *  \param lm Linked list matrix representation to be modidified.
 *  \param p Fill parameter.
 *  \param i Row to be computed.
 */
void ilup_levmat_row(levmat_t *lm, int p, int i);

/*! \brief Compute the ILU(p) factorisation.
 *
 *    \param A Matrix \f$ A \f$ to be factorised.
 *    \param p Allowable level of fill defining the non-zero pattern of
 *             the factorisation.
 *    \return Lower and upper triangular factor combined.
 *            The non-zero pattern of this matrix is a superset
 *            of the non-zero pattern of A.
 */
csr *ilup_csr(csr *A, int p);

#endif

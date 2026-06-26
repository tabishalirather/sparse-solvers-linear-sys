/*! \file
 *  \brief Preconditioners for incomplete factorisations (ILU).
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef PREC_ILU_H
#define PREC_ILU_H

prec_t *prec_ilu0_new(csr *a);
prec_t *prec_ilup_new(csr *a, int p);

#endif

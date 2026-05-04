/*! \file
 *  \brief Preconditioners corresponding to linear iterations.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef PREC_H
#define PREC_H

struct prec_s;
typedef struct prec_s prec_t;

struct prec_s {
  void (*func)(void *data, double *r, double *z);
  void (*delete) (prec_t *prec);
  void *data;
};

void prec_apply(prec_t *prec, double *r, double *z);
void prec_delete(prec_t *prec);

prec_t *prec_identity_new(csr *A);
prec_t *prec_jacobi_new(csr *A);
prec_t *prec_gs_new(csr *A);
prec_t *prec_sgs_new(csr *A);

#endif

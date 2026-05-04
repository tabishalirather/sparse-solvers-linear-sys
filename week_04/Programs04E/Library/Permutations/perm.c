/*! \file
 *  
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

perm_t *perm_new(int n, int *idx, int *idx_inv)
{
  perm_t *perm = malloc(sizeof(perm_t));
  perm->size = n;
  perm->idx = idx;
  if (idx_inv != NULL)
    perm->idx_inv = idx_inv;
  else
    perm->idx_inv = perminv(n, idx);
  return perm;
}

void perm_delete(perm_t *perm)
{
  free(perm->idx);
  free(perm->idx_inv);
  free(perm);
}

perm_t *perm_identity(int n) {
  int *idx     = malloc(n * sizeof(int));
  int *idx_inv = malloc(n * sizeof(int));
  int i;
  for (i = 0; i < n; ++i) {
    idx[i] = idx_inv[i] = i;
  }
  return perm_new(n, idx, idx_inv);
}

void perm_print(perm_t *perm)
{
  int i;
  printf("size = %d\n", perm->size);

  /* Permutation vector. */
  printf("idx     :");
  for (i = 0 ; i < perm->size; i++) {
    printf(" %d", perm->idx[i]);
  }
  printf("\n");

  /* Inverse permutation vector. */
  printf("idx_inv :");
  for (i = 0 ; i < perm->size; i++) {
    printf(" %d", perm->idx_inv[i]);
  }
  printf("\n");
}

bool perm_check(perm_t *perm)
{
  /* Compute inverse permutation. */
  int n = perm->size;
  int i;
  int *idx_inv_ref;

  idx_inv_ref = perminv(n, perm->idx);
  
  /* Check that idx_inv_ref equals to perm->idx_inv. */
  for (i = 0; i < n; ++i) {
    if (perm->idx_inv[i] != idx_inv_ref[i]) {
      free(idx_inv_ref);
      return false;
    }
  }
  free(idx_inv_ref);
  return true;
}

int *perminv(int n, const int *p)
{
  int i, *pinv;

  pinv = (int *) malloc(n * sizeof(int));
  
  for (i=0; i<n; i++)
    pinv[p[i]] = i;

  return pinv;

} /* perminv */

double *vec_permute(const int n, const double *b, const int *pinv)
{
  int i;
  double *Pb;

  Pb = (double *) malloc(n*sizeof(double));

  for (i=0; i<n; i++)
    Pb[i] = b[pinv ? pinv[i] : i];

  /* Here you can also use p instead of pinv:   Pb[p ? p[i] : i] = b[i] */

  return Pb;

} /* vec_permute */

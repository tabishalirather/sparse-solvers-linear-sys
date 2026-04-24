/*
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
*/

#include "itsolver.h"

csr *csr_permute(csr *A, const int *pinv, const int *q)
{
  csr *B = csr_alloc(A->nr, A->nc, A->nz);
  int ne = 0;
  int i;
  for (i = 0; i < B->nr; i++)
  {
    B->row[i] = ne;
    /* Use identity for pinv in case of NULL pointer. */
    /* The pinv[i]-th row of A becomes the i-th row of B. */
    int k = pinv ? (pinv[i]) : i;
    int ell;
    for (ell = A->row[k]; ell < A->row[k+1]; ell++)
    {
      B->coeff[ne] = A->coeff[ell];
      /* Use identity for q in case of NULL pointer. */
      /* The j-th column of A becomes the q[j]-th column of B. */
      B->col[ne++] = q ? (q[A->col[ell]]) : A->col[ell];
      /* Sort column indices. */
      int h;
      for (h = ne - 2; h >= B->row[i]; h--) {
        if (B->col[h] > B->col[h + 1]) {
          int j_tmp = B->col[h + 1];
          B->col[h + 1] = B->col[h];
          B->col[h] = j_tmp;
          double coeff_tmp = B->coeff[h + 1];
          B->coeff[h + 1] = B->coeff[h];
          B->coeff[h] = coeff_tmp;
        }
      }
    }
  }

  B->row[B->nr] = ne;

  return B;
}

csr *csr_perm(csr *A, perm_t *row_perm, perm_t *col_perm)
{
  return csr_permute(A, row_perm->idx_inv, col_perm->idx);
}

csr *csr_perm_symm(csr *A, perm_t *perm)
{
  return csr_perm(A, perm, perm);
}


/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void ilup_levmat_row(levmat_t *lm, int p, int i)
{
  assert(lm->nr == lm->nc);
  assert(p >= 0);
  assert(i >= 0 && i < lm->nr);
  /* For all existing A_{i, j} where j < i. */
  levent_t *aij = lm->entries[i];
  while ((aij = aij->next) != lm->entries[i] && aij->col < i)
  {
    int j = aij->col;
    /* For all existing A_{j, k} where k > j. */
    levent_t *ajk = lm->entries[j];
    /* Skip all k < j. */
    while ((ajk = ajk->next) != lm->entries[j] && ajk->col < j)
      ;
    /* Assume there is k = j. */
    assert(ajk != lm->entries[j]);
    assert(ajk->col == j);
    assert(ajk->val != 0.0);
    assert(ajk->val > 0.0); /* If A is SPD. */
    aij->val /= ajk->val;

    /* Proceed to k > j. */
    ajk = ajk->next;

    levent_t *aik = aij->next;
    while (ajk != lm->entries[j])
    {
      if (aik == lm->entries[i] || ajk->col < aik->col)
      {
        /* Insert a new entry. */
        int lev = aij->lev + ajk->lev + 1;
        if (lev <= p)
        {
          levent_t *aik_new = malloc(sizeof(levent_t));
          aik_new->col = ajk->col;
          aik_new->lev = lev;
          aik_new->val = -aij->val * ajk->val;
          levent_insert_before(aik, aik_new);
        }
        ajk = ajk->next;
      }
      else
      {
        if (ajk->col > aik->col)
        {
          /* Skip A_{i, k}. */
          aik = aik->next;
        }
        else
        {
          /* Update level in A_{i, k}.*/
          assert(ajk->col == aik->col);
          int lev = aij->lev + ajk->lev + 1;
          if (lev < aik->lev)
          {
            aik->lev = lev;
          }
          aik->val -= aij->val * ajk->val;
          aik = aik->next;
          ajk = ajk->next;
        }
      }
    }
  }
}

csr *ilup_csr(csr *A, int p)
{
  assert(A->nr == A->nc);
  assert(p >= 0);
  int n = A->nr;
  levmat_t *lm = levmat_from_csr(A);
  /* For all rows except the first. */
  int i;
  for (i = 1; i < n; ++i)
  {
    ilup_levmat_row(lm, p, i);
  }
  csr *LU = csr_from_levmat(lm);
  levmat_delete(lm);
  return LU;
}

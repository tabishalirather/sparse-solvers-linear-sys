/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void levent_insert_before(levent_t *list, levent_t *node) {
  node->prev = list->prev;
  node->next = list;
  list->prev->next = node;
  list->prev = node;
}

void levmat_print(levmat_t *lm) {
  printf("nr = %d, nc = %d\n", lm->nr, lm->nc);
  int i;
  for (i = 0; i < lm->nr; ++i) {
    printf("row=%d:", i);
    levent_t *curr = lm->entries[i];
    while ((curr = curr->next) != lm->entries[i]) {
      printf(" (col=%d, lev=%d, val=%f)", curr->col, curr->lev, curr->val);
    }
    printf("\n");
  }
}

levmat_t *levmat_from_csr(csr *A) {
  (void) levmat_print;
  int nr = A->nr;
  int nc = A->nc;
  levent_t **entries = malloc(nr * sizeof(levent_t *));
  int i, p;
  for (i = 0; i < A->nr; ++i) {
    entries[i] = malloc(sizeof(levent_t));
    entries[i]->next = entries[i];
    entries[i]->prev = entries[i];
    entries[i]->lev = -1;
    entries[i]->col = -1;
    entries[i]->val = NAN;
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      levent_t *le = malloc(sizeof(levent_t));
      le->lev = 0;
      le->col = A->col[p];
      le->val = A->coeff[p];
      le->next = NULL;
      le->prev = NULL;
      levent_insert_before(entries[i], le);
    }
  }
  levmat_t *lm = malloc(sizeof(levmat_t));
  lm->nr = nr;
  lm->nc = nc;
  lm->entries = entries;
  return lm;
}

csr *csr_from_levmat(levmat_t *lm) {

  int i;
  int nr = lm->nr;
  int nc = lm->nc;
  csr *A = malloc(sizeof(csr));
  A->nr = nr;
  A->nc = nc;
  A->row = malloc((nr + 1) * sizeof(int));
  A->row[0] = 0;
  for (i = 0; i < nr; ++i) {
    A->row[i + 1] = A->row[i];
    levent_t *curr = lm->entries[i];
    while ((curr = curr->next) != lm->entries[i]) {
      A->row[i + 1]++;
    }
  }
  int nz = A->row[nr];
  A->nz = nz;
  A->col = malloc(nz * sizeof(int));
  A->coeff = malloc(nz * sizeof(double));
  /* Copy column indices and coefficient values. */
  int p = 0;
  for (i = 0; i < nr; ++i) {
    levent_t *curr = lm->entries[i];
    while ((curr = curr->next) != lm->entries[i]) {
      A->col[p] = curr->col;
      A->coeff[p] = curr->val;
      ++p;
    }
  }
  assert(p == nz);
  return A;
}

void levmat_delete(levmat_t *lm) {
  int i;
  for (i = 0; i < lm->nr; ++i) {
    levent_t *curr = lm->entries[i];
    levent_t *curr_next;
    do {
      curr_next = curr->next;
      free(curr);
      curr = curr_next;
    } while (curr != lm->entries[i]);
  }
  free(lm->entries);
  free(lm);
}


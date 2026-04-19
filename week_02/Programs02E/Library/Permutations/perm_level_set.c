/*! \file
 *  
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

perm_t *perm_level_set(csr *A, int start, enum LEVEL_SET_VARIANT level_set_variant)
{
  assert(A->nr == A->nc);
  int n = A->nr;
  assert(n > 0);

  /* Level set permutation to be returned. */
  perm_t *perm = perm_identity(n);

  /* We want to assume n to be positive later on. */
  if (n == 0) {
    return perm;
  }

  /* Reference vector representation of the permutation and its inverse. */
  int *idx     = perm->idx_inv;
  int *idx_inv = perm->idx;

  /* Choose either the given node as a starting point or start from node zero
   * if the value of the parameter is out of range. */
  if (start >= 0 && start < n) {
    idx[0] = start;
    idx[start] = 0;
    idx_inv[0] = start;
    idx_inv[start] = 0;
  }

  /* Actual level set algorithm. */
  int i = 0;
  int j = 1;
  while (j < n) {
    /* Within the loop successors of node k are added to the permutation. */
    int k = idx[i];
    int p;
    /* Build a list of all umarked nodes reachable from node k. */
    /* First count the number of unmarked nodes reachable from k. */
    int list_len = 0;
    for (p = A->row[k]; p < A->row[k + 1]; ++p) {
      int l = A->col[p];
      /* Check whether node l is unmarked. */
      if (idx_inv[l] >= j) {
        ++list_len;
      }
    }
    /* Allocate and build list of unmarked nodes reachable from k. */
    int *list = malloc(list_len * sizeof(int));
    int *list_end = list + list_len;
    int *list_cur = list;
    for (p = A->row[k]; p < A->row[k + 1]; ++p) {
      int l = A->col[p];
      /* Check whether node l is unmarked. */
      if (idx_inv[l] >= j) {
        *list_cur++ = l;
      }
    }
    assert(list_cur == list_end);

    if (level_set_variant == LEVEL_SET_VARIANT_CMK ||
        level_set_variant == LEVEL_SET_VARIANT_RCM) {
      /* Sort list such that nodes nodes of lower degree appear before nodes of
       * higher degree. */
      /* Insertion sort. */
      for (list_cur = list + 1; list_cur < list_end; ++list_cur) {
        int *list1 = list_cur;
        do {
          int *list0 = list1 - 1;
          if (list0 >= list) {
            int idx0 = *list0;
            int idx1 = *list1;
            int deg0 = A->row[idx0 + 1] - A->row[idx0];
            int deg1 = A->row[idx1 + 1] - A->row[idx1];
            if (deg0 > deg1) {
              *list0 = idx1;
              *list1 = idx0;
            } else {
              break;
            }
          } else {
            break;
          }
          --list1;
        } while (true);
      }
    }

    /* Add each node from the list to the level set permutation. */
    for (list_cur = list; list_cur < list_end; ++list_cur) {
      int l = *list_cur;
      idx[idx_inv[l]] = idx[j];
      idx_inv[idx[j]] = idx_inv[l];
      idx[j] = l;
      idx_inv[l] = j;
      ++j;
    }

    /* Deallocate list. */
    free(list);

    /* Go on to check succesors of next node. */
    ++i;
  }

  /* Reverse permutation in case of RCM. */
  if (level_set_variant == LEVEL_SET_VARIANT_RCM) {

  /* ******************************************************************** */
  /*                                                                      */
  /*     TODO --- Problem 3                                               */
  /*                                                                      */
  /* ******************************************************************** */  

    /* Recompute the inverse permutation.
     * This cannot be computed by reversing the vector representation of the
     * inverse permutation. */
    for (i = 0; i < n; ++i) {
      idx_inv[idx[i]] = i;
    }
  }

  return perm;
}

perm_t *perm_bfs(csr *A, int start)
{
  return perm_level_set(A, start, LEVEL_SET_VARIANT_BFS);
}

perm_t *perm_cmk(csr *A, int start)
{
  return perm_level_set(A, start, LEVEL_SET_VARIANT_CMK);
}

perm_t *perm_rcm(csr *A, int start)
{
  return perm_level_set(A, start, LEVEL_SET_VARIANT_RCM);
}


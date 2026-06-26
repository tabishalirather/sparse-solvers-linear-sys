/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

static void perm_nested_dissection_rec(csr *A, int *queue, int *queue_inv, int *level, int rec_begin, int rec_end) {

  int components = 0;
  int head = rec_begin;
  int foot = rec_begin;
  int head_level = -1;
  while (foot < rec_end) {
    assert(head <= foot);
    if (head == foot) {
      level[foot] = ++head_level;
      foot++;
      components++;
      continue;
    }
    head_level = level[head];
    int p;
    for (p = A->row[queue[head]]; p < A->row[queue[head] + 1]; ++p) {
      int new = A->col[p];
      if (queue_inv[new] >= foot && queue_inv[new] < rec_end) {
        /* Swap in queue. */
        int idx = queue_inv[new];
        int tmp = queue[foot];

        assert(queue[idx] == new);
        assert(queue_inv[tmp] == foot);

        queue[idx] = tmp;
        queue[foot] = new;

        queue_inv[new] = foot;
        queue_inv[tmp] = idx;

        level[foot] = head_level + 1;
        foot++;
      }
    }
    head++;
  }
  /* Assert the graph is connected. */
  //assert(components <= 1);

  int level_max = -1;
  if (rec_end - rec_begin > 0) {
    level_max = level[rec_end - 1];
  }

  /* Choose separator D by minimizing abs(#S1 - #S2). */
  /* Assert we have at least three distinct levels 0, 1, 2. */
  if (level_max >= 2) {

    /* Find the largest index i such that
       the number of elements on all levels <= i is
       smaller than the number of all elements of levels > i. */
    int ptr_a = rec_begin;
    int ptr_b = rec_begin;
    int level_tmp = level[ptr_b];
    while (level_tmp == level[++ptr_b]);
    int ptr_c = ptr_b;
    level_tmp = level[ptr_c];
    while (level_tmp == level[++ptr_c]);
    do {
      /* Advance all three pointers. */

      level_tmp = level[ptr_a];
      while (level_tmp == level[++ptr_a]);
      level_tmp = level[ptr_b];
      while (level_tmp == level[++ptr_b]);
      level_tmp = level[ptr_c];
      while (level_tmp == level[++ptr_c]);
    } while (ptr_b < rec_end - ptr_c);

    int s1_begin = rec_begin;
    int d_begin = ptr_a;
    int s2_begin = ptr_b;
    if (ptr_b - (rec_end - ptr_c) < (rec_end - ptr_b) - ptr_a) {
      d_begin = ptr_b;
      s2_begin = ptr_c;
    }
    int cardS1 = d_begin - s1_begin;
    int cardS2 = rec_end - s2_begin;

    /* Move D to the end of the permuation array. */
    int idx = d_begin;
    while (idx < s2_begin) {

      int idx_swp = rec_end - 1 - (idx - d_begin);

      int tmp = queue[idx];
      int tmp_swp = queue[idx_swp];

      queue[idx] = tmp_swp;
      queue[idx_swp] = tmp;

      queue_inv[tmp] = idx_swp;
      queue_inv[tmp_swp] = idx;

      idx++;
    }

    /* Recursion. */
    /* S1. */
    perm_nested_dissection_rec(A, queue, queue_inv, level, rec_begin, rec_begin + cardS1);
    /* S2. */
    perm_nested_dissection_rec(A, queue, queue_inv, level, rec_begin + cardS1, rec_begin + cardS1 + cardS2);
  }
}

perm_t *perm_nested_dissection(csr *A) {
  assert(A->nr == A->nc);
  int n = A->nr;
  perm_t *perm = perm_identity(n);

  /* Call recursive version of this function with level as an work array. */
  int *level = malloc(n * sizeof(int));
  perm_nested_dissection_rec(A, perm->idx_inv, perm->idx, level, 0, n);
  free(level);

  return perm;
}

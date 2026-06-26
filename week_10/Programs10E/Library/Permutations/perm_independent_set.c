/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

perm_t *perm_independent_set(csr *A) {
  assert(A->nr == A->nc);
  int n = A->nr;
  int *color_used = malloc(n * sizeof(int));
  int *color = malloc(n * sizeof(int));
  int i;
  /* An entry color_used[c] = i means that a neighbor of node i
   * is of color c. Initialize to -1 since i = 0, ..., n-1. */
  for (i = 0; i < n; ++i) {
    color_used[i] = -1;
  }
  /* Assign colors to all nodes. */
  for (i = 0; i < n; ++i) {
    /* Find out which color cannot be used for i. */
    int p;
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      int j = A->col[p];
      if (j < i) {
        color_used[color[j]] = i;
      }
    }
    /* Pick first available color for i. */
    int c = 0;
    while (c < n && color_used[c] == i) {
      c++;
    }
    assert(c < n);
    color[i] = c;
  }
  /* Set color_used to zero. */
  for (i = 0; i < n; ++i) {
    color_used[i] = 0;
  }
  /* For each color c store in color_used[c] the number of
   * nodes which have color c-1. */
  for (i = 0; i < n; ++i) {
    if (color[i] < n - 1) {
      ++color_used[color[i] + 1];
    }
  }
  /* For each color c store in color_used[c] the number of
   * nodes which have color < c. */
  for (i = 1; i < n; ++i) {
    color_used[i] += color_used[i - 1];
  }

  /* Collect nodes of each color. */
  int *queue = malloc(n * sizeof(int));
  int *queue_inv = malloc(n * sizeof(int));
  for (i = 0; i < n; ++i) {
    int j = color_used[color[i]]++;
    queue[j] = i;
    queue_inv[i] = j;
  }

  free(color);
  free(color_used);

  return perm_new(n, queue_inv, queue);
}

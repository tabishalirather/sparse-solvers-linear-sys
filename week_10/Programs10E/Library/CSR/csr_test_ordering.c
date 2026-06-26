/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "itsolver.h"

void csr_test_ordering(csr *A, char *name, int start)
{
  /* create buffer for the filename */

  int slen = strlen(name) + 20;
  char *filename = malloc(slen);

  /* show A before symmetric permutations */
  snprintf(filename, slen, "%s_vorher.eps", name);
  csr_spy(A, filename, 's');

  /* compute startnode if necessary */
  if (start == -1) {
    start = csr_mindeg_startnode(A);
  }

  /* Level set algorithms BFS, CMK and RCM. */
  char *level_set_string[] = {"BFS", "CMK", "RCM"};
  int level_set_variant;
  for (level_set_variant = 0; level_set_variant < 3; ++level_set_variant) {

    /* Compute permutation and its inverse. */
    perm_t *perm = perm_level_set(A, start, level_set_variant);

    /* Display permutation. */
    printf("%s-Permutation:\n", level_set_string[level_set_variant]);
    perm_print(perm);

    /* Check permutation for validity. */
    if (!perm_check(perm)) {
      printf("Fehler: Indexvektoren stellen keine Permutation und ihre Inverse dar.\n");
    }

    /* Permute A symmetrically from both sides. */
    csr *B = csr_perm_symm(A, perm);

    /* Plot B's sparsity structure. */
    snprintf(filename, 100, "%s_nach%s.eps", name,
             level_set_string[level_set_variant]);
    csr_spy(B, filename, 's');

    /* Free memory. */
    csr_free(B);
    perm_delete(perm);
  }

  free(filename);
}

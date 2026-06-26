/*! \file
 *  \brief ILU(p) factorisation of a matrix whose unknowns are permuted in different ways.
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

int has_diagonal(csr *A) {
  int has = 1;
  int i;
  for (i = 0; i < A->nr; ++i) {
    int p;
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      if (A->col[p] == i) {
        break;
      }
    }
    if (p == A->row[i + 1]) {
      has = 0;
      break;
    }
  }
  return has;
}

int main(int argc, char **argv) {
  char filename[100];

  int d = 3;
  if (argc > 1) {
    d = atoi(argv[1]);
  }
  int n = 10;
  if (argc > 2) {
    n = atoi(argv[2]);
  }
  int p = 5;
  if (argc > 3) {
    p = atoi(argv[3]);
  }
  printf("n = %d, d = %d, p = %d\n", n, d, p);

  /* Testmatrices. You can observe that ILU(p) works well for the Laplace 
  problem but fails for the orsirr matrix */
  csr *A = csr_laplace(n, d);
  // csr *A = csr_import_from_mm_mtx("../../TestMatrixCollection/orsirr_1.mtx");
  csr_spy(A, "testprog_ilup_perm_initial_matrix.eps", 'c');

  /* Otherwise ILU(p) won't succeed. */
  printf("A has a complete diagonal: %s\n", has_diagonal(A) ? "yes" : "no");

  /* Declare variables. */
  perm_t *perm;
  csr *B, *C;
  

  // * Identity.

  // Compute the permutation.
  perm = perm_identity(A->nr);

  // Permute A.
  B = csr_perm_symm(A, perm);
  csr_spy(B, "A_perm_ID.eps", 'c');

  // Apply ILU(p) on the permuted matrix.
  C = ilup_csr(B, p);
  sprintf(filename, "ILU(%u)_ID.eps", p);
  csr_spy(C, filename, 'c');

  
  /* ************************************ */
  /*                                      */
  /*         TODO --- Problem 1b          */
  /* Apply the permutations for Cuthill   */
  /* McKee, Nested Dissection, and        */
  /* Independent Set, respectively.       */
  /*                                      */
  /* ************************************ */


  // * Cuthill McKee.

  // Compute the permutation.
  perm = perm_cmk(A, -1);


  // Permute A.
  B = csr_perm_symm(A, perm);
  csr_spy(B, "A_perm_CMK.eps", 'c');

  // Apply ILU(p) on the permuted matrix.
  C = ilup_csr(B, p);
  sprintf(filename, "ILU(%u)_CMK.eps", p);
  csr_spy(C, filename, 'c');

  
  // * Nested Dissection.

  // Compute the permutation.
  perm = perm_nested_dissection(A);

  // Permute A.
  B = csr_perm_symm(A, perm);
  csr_spy(B, "A_perm_ND.eps", 'c');

  // Apply ILU(p) on the permuted matrix.
  C = ilup_csr(B, p);
  sprintf(filename, "ILU(%u)_ND.eps", p);
  csr_spy(C, filename, 'c');

  
  // * Independent Set.

  // Compute the permutation.
  perm = perm_independent_set(A);

  // Permute A.
  B = csr_perm_symm(A, perm);
  csr_spy(B, "A_perm_IS.eps", 'c');

  // Apply ILU(p) on the permuted matrix.
  C = ilup_csr(B, p);
  sprintf(filename, "ILU(%u)_IS.eps", p);
  csr_spy(C, filename, 'c');



  /* Free memory. */
  csr_free(C);
  csr_free(B);
  csr_free(A);
  perm_delete(perm);

  return 0;
}

/*! \file
 *  \brief Amount of fillin generated with ILU(p) for different values of p.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

void print_matrix_csr(csr* matrix_p) {
  printf("\nPrinting in print_matrix_csr fxn\n");

  for (int i=0; i<matrix_p->nr; i++) {
    // printf("%d ", matrix_p->row[i]);
    for (int j=matrix_p->row[i]; j<matrix_p->row[i+1]; j++) {
      printf("%0.2f ", matrix_p->coeff[j]);
    }
    printf("\n");
  }
}



int main(int argc, char **argv)
{
  int d = 3;
  if (argc > 1)
  {
    d = atoi(argv[1]);
  }
  int n = 10;
  if (argc > 2)
  {
    n = atoi(argv[2]);
  }
  int pmax = 20;
  if (argc > 3)
  {
    pmax = atoi(argv[3]);
  }
  printf("n = %d, d = %d, pmax = %d\n", n, d, pmax);
  int nd = pow(n - 1, d);
  // printf("just after nd");

  /* Test matrix. */
  csr *A = csr_laplace(n, d);

  /* Permutation. */
  perm_t *perms[4];
  double *nnz_rel[4];

  perms[0] = perm_identity(nd);
  perms[1] = perm_cmk(A, -1);
  perms[2] = perm_nested_dissection(A);
  perms[3] =  perm_independent_set(A);
  // printf("%d, in for loop:", 2);

  /* ************************************ */
  /*                                      */
  /*         TODO --- Problem 1a          */
  /* Apply permutation and store the      */
  /* relative number of nonzero entries.  */
  /*                                      */
  /* ************************************ */

  for (int i = 0; i < 4; i++)
  {
    // printf("%d, in for loop: \n", i);

    /* Reserve storage for the relative number of nonzero entries */
    nnz_rel[i] = (double *) malloc((pmax+1) * sizeof(double));
    
    /* Apply a permutation */
    csr *B = csr_perm_symm(A, perms[i]);
    // print_matrix_csr(B);
    for (int p = 0; p <= pmax; p++)
    {
      printf("%d, in p loop:\n", p);
      /* ILU(p). */
      csr *C = ilup_csr(B,p);

      /* Compute and store relative number of non-zeros. */
      nnz_rel[i][p] = C->nz - A->nz ;

      /* Free memory.*/
      csr_free(C);
    }
  }

  /* Free memory.*/
  csr_free(A);
  for (int i = 0; i < 4; i++)
    perm_delete(perms[i]);

  /* Output File */
  FILE *OUT;
  // File is overwritten (and created if not existing)
  OUT = fopen("testprog_10_ilup_fillin.txt", "w");

  /* Write data to the output file */ 
  for (int p = 0; p <= pmax; p++)
  {

    /* Print p */
    fprintf(OUT, "%2d", p);

    /* Print relative number of non-zero entries for each permutation */
    for (int i = 0; i < 4; i++)
      fprintf(OUT, " %1.6e", nnz_rel[i][p]);
    fprintf(OUT, "\n");
  }

  /* Close output file */
  fclose(OUT);

  /* Free memory*/
  for (int i = 0; i < 4; i++)
    free(nnz_rel[i]);

  return 0;
}

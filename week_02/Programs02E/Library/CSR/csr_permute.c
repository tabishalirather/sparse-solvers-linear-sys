/*
 *
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

void print_vec_2(int n, const int* vec) {
  if (vec == NULL) {
    printf("vec is null");
  }

  printf("[ ");
  for (int i=0; i<n; i++) {

    printf("%d ", vec[i]);
  }
  printf("]");
  printf("\n");
}
//
// void print_vec_db(int n, double* vec) {
//   for (int i=0; i<n; i++) {
//     printf("%f ", vec[i]);
//   }
//   // printf("\n");
// }

csr *csr_permute(csr *A, const int *pinv, const int *q)
{
  csr *B = csr_alloc(A->nr, A->nc, A->nz);
  int ne = 0;
  int i;
  /* ******************************************************************** */
  /*                                                                      */
  /*     TODO --- Problem 1                                               */
  /*                                                                      */
  /* ******************************************************************** */
    printf("Inside pinv not null csr-permute");
    puts("------------------------Start csr_permnute-------------------------");
    printf("Print matrix_A below");
    print_matrix_csr(A);
    // if (q==NULL) {
    //   printf("q is null\n");
    // }
  //Can I print the matrix by columns. Select the col index in a row, and then skip by number of elements in the row and the print again.
  // int s=0;
  puts("Trying to print by column");
  for (int s=0; s<A->nr; s++) {
  for (int i=0; i<A->nr;i++) {

    for (int j=A->row[i]+s; j<A->row[i]+s+1; j++) {
      printf("%.1f ", A->coeff[j]);
    }
  }
    printf("\n");

  }

  /*for (int i=0; i<A->nc;i++) {
    // int k=0;
    for (int j=A->row[i]+1; j<A->row[i]+2; j++) {
      printf("%.1f ", A->coeff[j]);
    }
  }
  printf("\n");
  for (int i=0; i<A->nc;i++) {
    // int k=0;
    for (int j=A->row[i]+2; j<A->row[i]+3; j++) {
      printf("%.1f ", A->coeff[j]);
    }
  }*/
  printf("\n");
  printf("\n");

  int k=0;
  if (pinv!=NULL && q==NULL) {
    puts("vector pinv is: ");
    print_vec_2(3, pinv);
    puts("This returns PA");
    k=0;
    int row_test_idx = 0;
    for (int i=0; i<A->nr; i++) {
      for (int j=A->row[pinv[i]]; j<A->row[pinv[i]+1]; j++) {
        printf("\n");
        B->row[0] = 0;
        B->row[i+1] = row_test_idx++;
        B->coeff[j] = A->coeff[k++];
      }
      B->row[i+1]++;
    }
    for (int i=0; i<A->nz; i++) {
      B->col[i] = A->col[i];
    }
  }
  else if (pinv == NULL && q!=NULL) {
    //What we are saying now is that this should a column transpose, not a row transpose.
    puts("vector q is: ");
    print_vec_2(3, q);
    puts("This returns AP");
    for (int i=0; i<A->nr;i++) {
      for (int j=A->row[q[i]]; j<A->row[q[i]+1]; j++){
        B->coeff[k++] = A->coeff[j];
        printf("%0.1f ", A->coeff[j]);       // printf(" (k,i):(%d,%d) ", k,i);
      }
      // printf(" (k: %d) ", k);
      B->row[0] = 0;
      B->row[i+1] = k;
      // printf("\n");
    }
    for (int i=0; i<A->nz; i++) {
      B->col[i] = A->col[i];
    }
  }
     //apply above algorithm to A matrix now.
if (pinv != NULL ^ q!=NULL) {

  puts("B.coeff is:");
  for (int i=0; i<B->nz;i++) {
    printf("%.1f ", B->coeff[i]);
  }

  printf("\n");
  puts("permuted_row is: ");
  for (int i=0; i<=B->nr; i++) {
    printf("%d ", B->row[i]);
  }

  printf("\n");
  printf("Matrix B printed below");
  print_matrix_csr(B);
  puts("\n------------------------End   csr_permnute-------------------------");

  printf("**********************Outside pinv not null csr-permute****************************");
}
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


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
    if (q==NULL) {
      printf("q is null\n");
    }
    // // printf("pinv is: \n");
    if (pinv != NULL) {
      print_vec_2(3, pinv);
    }else {
      pinv = q;
    }
  if (pinv!=NULL && q==NULL) {
    puts("This returns PA");
  } else if (q!=NULL) {
    puts("This returns AP");
  }
    //apply above algorithm to A matrix now.
    int k=0;
    for (int i=0; i<A->nr;i++) {
      for (int j=A->row[pinv[i]]; j<A->row[pinv[i]+1]; j++){
        B->coeff[k++] = A->coeff[j];
        printf("%0.1f ", A->coeff[j]);       // printf(" (k,i):(%d,%d) ", k,i);
      }
      // printf(" (k: %d) ", k);
      B->row[0] = 0;
      B->row[i+1] = k;
      printf("\n");
    }
    for (int i=0; i<A->nz; i++) {
      B->col[i] = A->col[i];
    }

    /*// for (int i=0; i<3;i++) {
    //   for (int j=A->row[pinv[i]]; j<A->row[pinv[i]+1]; j++){
    //     B->coeff[k++] = A->coeff[j];
    //     printf("%d ", B->coeff[k++]);
    //   }
    //   printf("\n");
    // }*/
    printf("permuted vector is: \n");
    for (int i=0; i<A->nz; i++) {
      printf("%.f ", B->coeff[i]);
    }
    printf("\n");
    printf("permuted row is: \n");
    for (int i=0; i<A->nr+1; i++) {
      printf("%d ", B->row[i]);
    }
    printf("\n");
    printf("permuted col is: \n");
    for (int i=0; i<B->nz; i++) {

      printf("%d ", B->col[i]);
    }
    printf("\n");

    //matrix permuted, now we need to figure out how to sort out the row variable to suit the matrix.
    printf("Matrix B printed below");
    print_matrix_csr(B);
    puts("\n------------------------End   csr_permnute-------------------------");
    return B;
  printf("**********************Outside pinv not null csr-permute****************************");

}
csr *csr_perm(csr *A, perm_t *row_perm, perm_t *col_perm)
{
  return csr_permute(A, row_perm->idx_inv, col_perm->idx);
}

csr *csr_perm_symm(csr *A, perm_t *perm)
{
  return csr_perm(A, perm, perm);
}


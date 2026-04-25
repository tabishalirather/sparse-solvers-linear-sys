/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
*/

#include "itsolver.h"

void print_matrix_csr(csr *matrix_p) {
  printf("\nPrinting in print_matrix_csr fxn\n");

  for (int i = 0; i < matrix_p->nr; i++) {
    // printf("%d \n", matrix_p->row[i]);
    for (int j = matrix_p->row[i]; j < matrix_p->row[i + 1]; j++) {
      printf("%0.2f ", matrix_p->coeff[j]);
      printf("%d ", matrix_p->col[j]);
      // printf("\n(i,j):(%d,%d\n ", i, j);
    }
    printf("\n");
  }
}

void print_vec_2(int n, const int *vec) {
  if (vec == NULL) {
    printf("vec is null");
  }

  printf("[ ");
  for (int i = 0; i < n; i++) {
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

csr *csr_permute(csr *A, const int *pinv, const int *q) {
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
  int nr = 3;
  int nc = 3;
  int nz = 5;

  int row[4] = {0, 1, 3, 5};
  int col[5] = {0, 1, 2, 0, 2};
  // From this we have 4 belongs
  // 4 and 1 belong to col 0
  // 7 and 9 belong to col 2
  // 5 belongs to col 1
  // what we can do now is store coeff in the
  double coeff[5] = {4, 7, 5, 1, 9};

  int perm[3] = {2, 0, 1};
  int sorted_col[5] = {0};
  int to_sort[4] = {3, 4, 2, 1};
  int to_sort_sorted[4] = {0};

  // puts("Before sorting, the var sorted_col is: ");
  //
  // printf("\n");
  //
  // printf("\n");
  //
  // printf("\n");
  // printf("\n");

  int k = 0;
  if (pinv != NULL && q == NULL) {
    puts("vector pinv is: ");
    print_vec_2(3, pinv);
    puts("This returns PA");
    k = 0;
    int row_test_idx = 0;
    for (int i = 0; i < A->nr; i++) {
      for (int j = A->row[pinv[i]]; j < A->row[pinv[i] + 1]; j++) {
        printf("\n");
        B->row[0] = 0;
        B->row[i + 1] = row_test_idx++;
        B->coeff[j] = A->coeff[k++];
      }
      B->row[i + 1]++;
    }
    for (int i = 0; i < A->nz; i++) {
      B->col[i] = A->col[i];
    }
  }

  else if (pinv == NULL && q != NULL) {
    puts("This returns AP\n");
    //What we are saying now is that this should a column transpose, not a row transpose.
    puts("vector q is:");
    print_vec_2(3, q);

    puts("A.col is:");
    print_vec_2(A->nz, A->col);
    // puts("B.coeff is:");
    for (int i = 0; i < A->nz; i++) {
      // printf("i is: %d ", i);
      // printf("B.col[i] is: %d  ", B->col[i]);
      B->col[i] = q[A->col[i]];
      B->coeff[i] = A->coeff[i];
    }
    for (int i=0; i<=B->nr; i++) {
      B->row[i] = A->row[i];
    }
    printf("\n");
    puts("B.row is: ");
    for (int i=0; i<=B->nr; i++) printf("%d ", B->row[i]);

    puts("\nA.row is: ");
    for (int i=0; i<=A->nr; i++) printf("%d ", A->row[i]);

    printf("\n");
    puts("Permuted col: B.col[i]");
    for (int i = 0; i < B->nz; i++) {
      printf("%d ", B->col[i]);
    }
    printf("\n");
    puts("Original col: A.col[i]");
    for (int i = 0; i < A->nz; i++) {
      printf("%d ", A->col[i]);
    }
    printf("\n");
    puts("B.coeff[i] is:");
    for (int i=0; i<B->nz; i++) {
      printf("%f ", B->coeff[i]);
    }



    printf("\n");
    printf("Matrix B printed below");
    print_matrix_csr(B);
    puts("\n------------------------End   csr_permnute-------------------------");

    printf("**********************Outside pinv not null csr-permute****************************");
  }
  return B;
}

csr *csr_perm(csr *A, perm_t *row_perm, perm_t *col_perm) {
  return csr_permute(A, row_perm->idx_inv, col_perm->idx);
}

csr *csr_perm_symm(csr *A, perm_t *perm) {
  return csr_perm(A, perm, perm);
}

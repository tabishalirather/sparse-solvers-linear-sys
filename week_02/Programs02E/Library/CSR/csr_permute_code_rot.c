/*
/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
#1#

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
  // puts("calling print_vec_2:\n");
  // printf("vec[0] is: %d ", vec[0]);
  // printf("%d \n", n);
  printf("[ ");
  for (int i=0; i<n; i++) {
    // puts("inside loop");
    // printf("i: %d \n", i);
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
  /* ******************************************************************** #1#
  /*                                                                      #1#
  /*     TODO --- Problem 1                                               #1#
  /*                                                                      #1#
  /* ******************************************************************** #1#
    printf("Inside pinv not null csr-permute");
    puts("------------------------Start csr_permnute-------------------------");
    printf("Print matrix_A below");
    print_matrix_csr(A);
    if (q==NULL) {
      printf("q is null\n");
    }

    printf("pinv is: \n");
    if (pinv != NULL) {
      print_vec_2(3, pinv);
    }else {
      pinv = q;
    }

    for (int i=0; i<A->nz; i++) {
      // B->coeff[i] = A->coeff[i];
      for (int j=0; j<A->nr; j++) {
        // B->row[i] = A->row[i];
        // B->col[i] = A->col[i];
      }
    }

    // printf("A.nr: %d\n", A->nr);
    // printf("A.nc: %d\n", A->nc);


    //working with csr format
    puts("testing stuff with csr format");
    int nr = 3;
    int nc = 3;
    int nz = 5;
    int row[4] = {0, 2, 3, 5};
    int col[5] = {0, 2, 1, 0, 2};
    int coeff[5] = {4, 7, 5, 1, 9};
    int permuted_coeff[5] = {0};
    int permuted_row[4] = {0};
    // int p[3] = {1, 2, 0};
    int k=0;
    /#1#/ for (int i=0; i<nr; i++) {
      // printf(" ", i);
      // for (int j=row[i]; j<row[i+1]; j++) {
        // this printf, prints a row, in other words selects a row. Now we just need to permute the selected row.
        // printf("(i,j:%d, %d): ",i, j);
        // printf("%d \n", coeff[j]);
        // permuted_coeff[p[i]] = coeff[j];

        // permuted_index[i] += j;
        //What we want is this: for p = (1 2 0). Select row with index i and store in index(p(i))
        //manually, this is what I'd do:
        //first select row with index p(0) = 1 which means the second row. and then store this in the first position in a vector
        //next select row with index p(1) = 2 which means the third row and store this in  the second position in a vector
        //next select row with index p(2) = 0 which means the first row and store this in the third position.
        // Right now each run of the for-j loop, selects a row. The challenge is to figure out how to store these elements in the place determined by the by index of the permutation element.
        // For example here, in the first iteration, i=0, which means the first row. p(0) = 1

  //(i,j:0, 0): 4
  //(i,j:0, 1): 7
  // new--row
  //(i,j:1, 2): 5
  //new--row
  //(i,j:2, 3): 1
  //(i,j:2, 4): 9
  // new--row
      // }
      // permuted_index[i] += i;
      // printf("new--row\n");
    // }#1#

    // for (int i=0; i<3;i++) {
    //   for (int j=row[p[i]]; j<row[p[i]+1]; j++){
    //     permuted_coeff[k++] = coeff[j];
    //     printf("%d ", coeff[j]);
    //   }
      // printf(" (k: %d) ", k);
    //   permuted_row[i+1] = k;
    //   printf("\n");
    // }
    //apply above algorithm to A matrix now.
    for (int i=0; i<nr;i++) {
      for (int j=A->row[pinv[i]]; j<A->row[pinv[i]+1]; j++){
        B->coeff[k++] = A->coeff[j];
        printf("%0.1f ", A->coeff[j]);
        // printf(" (k,i):(%d,%d) ", k,i);
      }
      // printf(" (k: %d) ", k);
      B->row[0] = 0;
      B->row[i+1] = k;
      printf("\n");
    }
    for (int i=0; i<A->nz; i++) {
      B->col[i] = A->col[i];
    }

    /#1#/ for (int i=0; i<3;i++) {
    //   for (int j=A->row[pinv[i]]; j<A->row[pinv[i]+1]; j++){
    //     B->coeff[k++] = A->coeff[j];
    //     printf("%d ", B->coeff[k++]);
    //   }
    //   printf("\n");
    // }#1#
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
*/


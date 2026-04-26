/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
*/

#include "itsolver.h"



void sort_csr_matrix_col(csr *matrix) {
  puts("Calling sort_csr_matrix_col fxn");
  for (int i = 0; i < matrix->nr; i++) {
    int start = matrix->row[i];
    int end = matrix->row[i + 1];

    for (int j = start; j < end - 1; j++) {
      for (int k = j + 1; k < end; k++) {
        if (matrix->col[k] < matrix->col[j]) {
          int temp_col = matrix->col[j];
          matrix->col[j] = matrix->col[k];
          matrix->col[k] = temp_col;

          double temp_coeff = matrix->coeff[j];
          matrix->coeff[j] = matrix->coeff[k];
          matrix->coeff[k] = temp_coeff;
        }
      }
    }
  }
  // printf("\n");
}

void return_AP_APT(csr *A, csr *B, const int *q) {
  for (int i = 0; i < A->nz; i++) {
    B->col[i] = q[A->col[i]];
    B->coeff[i] = A->coeff[i];
  }
  for (int i = 0; i <= B->nr; i++) {
    B->row[i] = A->row[i];
  }
  sort_csr_matrix_col(B);

}
void return_PA_PTA(csr *A, csr *B, const int *pinv) {
  //this is row permutation/
  puts("Printing from fxn return_PA");
  int k = 0;
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


void print_matrix_csr(csr *matrix_p) {
  printf("\nPrinting in print_matrix_csr fxn\n");

  for (int i = 0; i < matrix_p->nr; i++) {
    for (int j = matrix_p->row[i]; j < matrix_p->row[i + 1]; j++) {
      printf("%0.2f ", matrix_p->coeff[j]);
      // printf("%d ", matrix_p->col[j]);
    }
    printf("\n");
  }
}

void print_vec_int(int n, const int *vec) {
  // printf(title)
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

void print_vec_double(int n, double* vec) {
  // printf(title)
  if (vec == NULL) {
    printf("double* vec is null");
  }

  printf("[ ");
  for (int i = 0; i < n; i++) {
    printf("%f ", vec[i]);
  }
  printf("]");
  printf("\n");
}

void print_attrs_B(csr *B) {
  puts("B.row is: ");
  print_vec_int(B->nr + 1, B->row);

  puts("Permuted B.col[i] is: ");
  print_vec_int(B->nz, B->col);

  puts("B.coeff[i] is:");
  print_vec_double(B->nz, B->coeff);

  puts("Permuted B.coeff[i] is:");
  print_vec_double(B->nz, B->coeff);
}
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

  if (pinv != NULL && q == NULL) {
    puts("vector pinv is: ");
    print_vec_int(3, pinv);
    puts("This returns PA and PTA");
    //  PA = csr_permute(A, pinv=iperm, q=NULL);
    // PTA   = csr_permute(A, pinv=perm, q=NULL);
    // The idea we use inverse of pinv argument, so if arg has perm i.e. P, then code uses PT (iperm) for permutation. Therefore if you want PA, then use p=iperm as argument

    // q = NULL and pinv = iperm for PA and pinv = perm for PTA
    return_PA_PTA(A, B, pinv);
  } if (pinv == NULL && q != NULL) {
    puts("This returns AP\n");
    //What we are saying now is that this should a column permutation not a row permutation.
    puts("vector q is:");
    print_vec_int(3, q);
    //  AP = csr_permute(A, pinv=NULL, q=perm); When using q=perm, this is not inverted and is used as if. I you want APT, then pass q=iperm like below:
    //   APT   = csr_permute(A, NULL, iperm);

    //pinv = NULL and q = iperm for APT and q = perm for AP
    return_AP_APT(A, B, q);
    // printf("Matrix B printed below");
    // print_matrix_csr(B);
  } if (pinv!=NULL && q != NULL) {
    puts("This returns dual-permutations");
    //To get PAP:   PAP   = csr_permute(A, pinv=iperm, q=perm);
    // First we do (PA) and then (PA)*P. We use pinv to get PA as the code is setup such that for PA we must pass pinv and for PTA,we must pass p(or here p=q i.e q = (pinv)^-1 = p)

    // pinv = ipermn and q = perm
    return_PA_PTA(A, B, pinv);
    return_AP_APT(B, B, q);

    //To get PAPT:
    //Here we need PA first, so we do with pinv as that is converted to P, in the second fxn call we use pinv as well to do PT as the second fxn uses argumetn directly for column permutations without inversion nagain.

    //pinv = iperm and q = iperm
    return_PA_PTA(A, B, pinv);
    return_AP_APT(B, B, q);

    // To get PTAP:
    // Since here we need PTA first, we must use q, as q, if we use pinv, then code will convert that to p and we will get PA, not PTA, to get PTA, we have to pass p, then code will use that to build PT. In simple, when pinv = perm, then we get PTA and when pinv = iperm, we get PA, and when q = perm, we get AP and when q = iperm, we get APT
    //since here want PT first, then pinv=perm and q=perm

    //pinv = perm and q = perm
    return_PA_PTA(A, B, pinv);
    return_AP_APT(B, B, q);

    //Now to get PTAPT
    // now we want PTA first, so

    // pinv = perm and q=ipem
    return_PA_PTA(A, B, pinv);
    return_AP_APT(B, B, q);

  }
    puts("\n------------------------End   csr_permnute-------------------------");
  return B;
}

csr *csr_perm(csr *A, perm_t *row_perm, perm_t *col_perm) {
  return csr_permute(A, row_perm->idx_inv, col_perm->idx);
}

csr *csr_perm_symm(csr *A, perm_t *perm) {
  return csr_perm(A, perm, perm);
}

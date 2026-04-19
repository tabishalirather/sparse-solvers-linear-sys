/*/*! \file
 *
 *  \author Sabine Le Borne
 #1#
/*
#include <iso646.h>

#include "itsolver.h"

perm_t *perm_new(int n, int *idx, int *idx_inv)
{
  perm_t *perm = malloc(sizeof(perm_t));
  perm->size = n;
  perm->idx = idx;
  if (idx_inv != NULL)
    perm->idx_inv = idx_inv;
  else
    perm->idx_inv = perminv(n, idx);
  return perm;
}


void perm_delete(perm_t *perm)
{
  free(perm->idx);
  free(perm->idx_inv);
  free(perm);
}

perm_t *perm_identity(int n) {
  int *idx     = malloc(n * sizeof(int));
  int *idx_inv = malloc(n * sizeof(int));
  int i;
  for (i = 0; i < n; ++i) {
    idx[i] = idx_inv[i] = i;
  }
  return perm_new(n, idx, idx_inv);
}

void print_vec(int n, int* vec) {
  for (int i=0; i<n; i++) {
    printf("%d ", vec[i]);
  }
  printf("\n");
}

void verify_pinv(int n, int* p, int* pinv) {
  int pass_count = 0;
  for (int i=0; i<n; i++) {
    // printf("p[pinv[%d]: %d \n", i, p[pinv[i]]);
    if (p[pinv[i]] == i and pinv[p[i]] == i) {
      pass_count++;
    }
  }
  // printf("pass count: %d", pass_count);
  if (pass_count == n) {
    printf("Good job, all test passed.");
  }
}

void print_matrix(int n, int* matrix_p) {
  for (int i =0; i<n; i++) {
    for (int j=0; j<n; j++){
      printf("%d ", *(matrix_p + i*n + j));

    }
    printf("\n");
  }
}

void make_pinv(int n, int* pinv, int* p) {
  for (int i=0; i<n; i++) {
    pinv[p[i]] = i;
  }
  printf("\npinv is:");
  print_vec(n, pinv);
}

void make_p_matrix_t(int n, int* matrix_p_t, int* matrix_p) {
  for (int i=0; i<n; i++) {
     for (int j=0; j<n; j++){
      *(matrix_p_t+ i*n +j) = *(matrix_p+ j*n +i);
      // printf("%d ", *(matrix_p_t+ i*n +j));
    }
    // printf("\n");
  }
  printf("matrix p_t is: \n");
  print_matrix(n, matrix_p_t);
}

void make_p_matrix(int n, int* matrix_p, int* p) {
  for (int i =0; i<n; i++) {
    for (int j=0; j<n; j++){
      // printf("j: [%d]; p[%d]: %d \n", j, i, p[i]);
      if (j == p[i])
        *(matrix_p + i*n + j) = 1;
    }
  }
  printf("matrix p is: \n");
  print_matrix(n, matrix_p);
}

void make_pinv_from_matrix(int n, int* matrix_p_t, int* pinv) {
  for (int i =0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (*(matrix_p_t +i*n +j) != 0) {
        pinv[i] = j;
      }
    }
  }
}


void perm_print(perm_t *perm)
{
  int i;
  printf("size = %d\n", perm->size);

  /* Permutation vector. #1#
  printf("idx     :");
  for (i = 0 ; i < perm->size; i++) {
    printf(" %d", perm->idx[i]);
  }
  printf("\n");

  /* Inverse permutation vector. #1#
  printf("idx_inv :");
  for (i = 0 ; i < perm->size; i++) {
    printf(" %d", perm->idx_inv[i]);
  }
  printf("\n");
}

bool perm_check(perm_t *perm)
{
  /* Compute inverse permutation. #1#
  int n = perm->size;
  int i;
  int *idx_inv_ref;

  idx_inv_ref = perminv(n, perm->idx);

  /* Check that idx_inv_ref equals to perm->idx_inv. #1#
  for (i = 0; i < n; ++i) {
    if (perm->idx_inv[i] != idx_inv_ref[i]) {
      free(idx_inv_ref);
      return false;
    }
  }
  free(idx_inv_ref);
  return true;
}

int *perminv(int n, const int *p) {
  int i, *pinv;
  printf("calling pinv from Library/perm.c \n");
  pinv = (int *) malloc(n * sizeof(int));
  printf("Intial pinv is:");
  print_vec(n, pinv);
  int *matrix_p = calloc(n*n, sizeof(int));
  int *matrix_p_t = calloc(n*n, sizeof(int));

  /* ******************************************************************** #1#
  /*                                                                      #1#
  /*     TODO ---  Problem 1                                              #1#
  /*                                                                      #1#
  /* ******************************************************************** #1#
  //

  //Ex: p = 1 2 0 means 1's are at: (0,1), (1,2) and (2,0)
  make_p_matrix(n, matrix_p, p);


  //write transpose of the matrix;
  // make_p_matrix_t(n, matrix_p_t, matrix_p);
  // make_pinv_from_matrix(n, matrix_p_t, pinv);
  //very this is indeed inverse.
  verify_pinv(n, p, pinv);

  printf("\n");
  printf("p: "); print_vec(n, p);
  printf("pinv: "); print_vec(n, pinv);

  printf("Efficient method:");
  //let's try to make this efficient now.
  make_pinv(n, pinv, p);


  free(matrix_p);
  free(matrix_p_t);
  return pinv;

} /* perminv #1#

double *vec_permute(const int n, const double *b, const int *pinv)
{
  int i;
  double *Pb;

  Pb = (double *) malloc(n*sizeof(double));


  /* ******************************************************************** #1#
  /*                                                                      #1#
  /*     TODO --- Problem 1                                               #1#
  /*                                                                      #1#
  /* ******************************************************************** #1#

  return Pb;

} /* vec_permute #1#
#1#*/
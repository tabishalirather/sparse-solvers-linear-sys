/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */


#include "csr.h"

#include <stdlib.h>

csr *csr_alloc(int rows, int cols, int elements)
{
  csr *new;
  new = (csr *) malloc(sizeof(csr));

  new->nr = rows;
  new->nc = cols;
  new->nz = elements;

  new->row     = (int *)  malloc((rows+1)*sizeof(int));
  new->col     = (int *)  malloc(elements*sizeof(int));
  new->coeff = (double *) malloc(elements*sizeof(double));

  return new;
}

void csr_free(csr *A)
{
  free(A->row);
  free(A->col);
  free(A->coeff);
  free(A);
}

csr *csr_copy(csr *A) {
  int i, j;
  int nr = A->nr;
  csr *B = csr_alloc(nr, A->nc, A->nz);
  for (i = 0; i <= nr; ++i) {
    B->row[i] = A->row[i];
  }
  for (i = 0; i < nr; ++i) {
    for (j = A->row[i]; j < A->row[i + 1]; ++j) {
      B->coeff[j] = A->coeff[j];
      B->col[j] = A->col[j];
    }
  }
  return B;
}

csr *csr_rand(int nr, int nc) {
  int i, j, k, q;
  int nnz_row_max = 7;
  int val_max = 3;
  csr *A = malloc(sizeof(csr));
  A->nr = nr;
  A->nc = nc;
  A->row = malloc((nr + 1) * sizeof(int));
  A->row[0] = 0;
  for (i = 0; i < nr; ++i) {
    int nnz_row = rand() % nnz_row_max + 1;
    if (nnz_row > nc) {
      nnz_row = nc;
    }
    A->row[i + 1] = A->row[i] + nnz_row;
  }
  A->nz = A->row[nr];
  A->col = malloc(A->nz * sizeof(int));
  A->coeff = malloc(A->nz * sizeof(double));
  for (i = 0; i < nr; ++i) {
    for (j = A->row[i]; j < A->row[i + 1]; ++j) {
      do {
	k = rand() % nc;
	for (q = A->row[i]; q < j; ++q) {
	  if (A->col[q] == k) {
	    k = -1;
	    break;
	  }
	}
      } while (k == -1);
      for (q = j; q > A->row[i] && A->col[q - 1] > k; --q) {
	A->col[q] = A->col[q - 1];
	A->coeff[q] = A->coeff[q - 1];
      }
      A->col[q] = k;
      A->coeff[q] = rand() % val_max + 1;
    }
  }
  return A;
}

int csr_bandwidth(csr *A) {
  int bw = 0;
  int i, j, bw_update;
  for (i = 0; i < A->nr; ++i) {
    for (j = A->row[i]; j < A->row[i + 1]; ++j) {
      bw_update = i < A->col[j] ? A->col[j] - i : i - A->col[j];
      if (bw < bw_update) {
	bw = bw_update;
      }
    }
  }
  return bw;
}

size_t csr_bytes(csr *A) {
  size_t bytes = sizeof(csr); /* A */
  bytes += (A->nr + 1) * sizeof(int); /* A->row */
  bytes += A->nz * sizeof(int); /* A->col */
  bytes += A->nz * sizeof(double); /* A->coeff */
  return bytes;
}

void csr_sort_col(csr *A) {
  int c, r, j, q;
  double val;
  for (r = 0; r < A->nr; ++r) {
    for (j = A->row[r] + 1; j < A->row[r + 1]; ++j) {
      c = A->col[j];
      val = A->coeff[j];
      q = j;
      for (; q >= A->row[r] + 1 && A->col[q - 1] > A->col[q]; --q) {
	A->col[q] = A->col[q - 1];
	A->coeff[q] = A->coeff[q - 1];
      }
      A->col[q] = c;
      A->coeff[q] = val;
    }
  }
}

void csr_sort_diagFirst(csr *A) {
  int c, r, j;
  double val;
  for (r = 0; r < A->nr; ++r) {
    if (A->col[A->row[r]] != r) {
      for (j = A->row[r]+1; j < A->row[r + 1]; ++j) {
        c = A->col[j];
        if (c == r ) {
          val = A->coeff[j];
	  A->col[j] = A->col[A->row[r]];
	  A->coeff[j] = A->coeff[A->row[r]];
	  A->col[A->row[r]] = r;
	  A->coeff[A->row[r]] = val;
	}
      }
    }
  }
}

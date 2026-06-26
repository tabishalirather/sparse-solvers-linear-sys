/* 
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

/* Matrix market IO header. */
#include "mmio.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "itsolver.h"

csr *csr_import_from_mm_mtx(const char *file_name) {
  MM_typecode matcode;
  FILE *f;
  long fpos;
  csr *A;
  int nr, nc, nz;
  int i, j, p;
  int i_read, j_read;
  double val;

  /* Try to open the file with the given name. */
  if ((f = fopen(file_name, "rb")) == NULL) {
    printf("File %s could not be opened for reading.\n", file_name);
    return NULL;
  }

  /* Read the first line of the file describing its contents. */
  if (mm_read_banner(f, &matcode)) {
    fclose(f);
    return NULL;
  }

  /* Make sure this is a sparse matrix with real entries in coordinate format. */
  if (!mm_is_matrix(matcode) ||
      !mm_is_sparse(matcode) ||
      !mm_is_coordinate(matcode) ||
      !mm_is_real(matcode)) {
    fclose(f);
    return NULL;
  }

  /* Get the number of rows, columns and non-zeros. */
  if (mm_read_mtx_crd_size(f, &nr, &nc, &nz)) {
    fclose(f);
    return NULL;
  }

  /* Allocate matrix in CSR-format. */
  A = csr_alloc(nr, nc, nz);

  /* Initialize row pointers to zero and use A->row[i + 1] to count
     the number of entries in row 0 <= i < nr. */
  for (i = 0; i <= nr; ++i) {
    A->row[i] = 0;
  }

  /* Remember position in stream. */
  fpos = ftell(f);

  /* Store the number of all elements in row i in A->row[i + 2] for
     i = 0, ..., nr - 2. */
  for (p = 0; p < nz; ++p) {
    fscanf(f, "%d %d %lg\n", &i_read, &j_read, &val);
    i = i_read - 1;
    if (i <= nr - 2) {
      A->row[i + 2] += 1;
    }
  }

  /* Sum up entries in A->row. */
  for (i = 3; i <= nr; ++i) {
    A->row[i] += A->row[i - 1];
  }

  /* Go back to prior position. */
  if (fseek(f, fpos, SEEK_SET)) {
    csr_free(A);
    fclose(f);
    return NULL;
  }

  /* Read matrix entries. */
  for (p = 0; p < nz; ++p) {

    fscanf(f, "%d %d %lg\n", &i_read, &j_read, &val);

    /* Adjust from one-based to zero-based. */
    i = i_read - 1;
    j = j_read - 1;

    assert(0 <= i && i < nr);
    assert(0 <= j && j < nc);

    A->coeff[A->row[i + 1]] = val;
    A->col[A->row[i + 1]] = j;
    A->row[i + 1] += 1;
  }

  /* Close the opened file. */
  fclose(f);

  /* Return matrix. */
  return A;
}

double *csr_import_from_mm_rhs(const char *file_name) {
  MM_typecode matcode;
  FILE *f;
  double *b;
  int nr, nc;
  int i;
  double val;

  /* Try to open the file with the given name. */
  if ((f = fopen(file_name, "rb")) == NULL) {
    printf("File %s could not be opened for reading.\n", file_name);
    return NULL;
  }

  /* Read the first line of the file describing its contents. */
  if (mm_read_banner(f, &matcode)) {
    fclose(f);
    return NULL;
  }

  /* Make sure this is a matrix with real entries in array format. */
  if (!mm_is_matrix(matcode) ||
      !mm_is_array(matcode) ||
      !mm_is_real(matcode)) {
    fclose(f);
    return NULL;
  }

  /* Get the number of rows, columns. */
  if (mm_read_mtx_array_size(f, &nr, &nc)) {
    fclose(f);
    return NULL;
  }

  /* Check that this is really a vector, i.e. a single column matrix. */
  if (nc != 1) {
    fclose(f);
    return NULL;
  }

  /* Allocate vector. */
  b = (double *) malloc(nr * sizeof(double));

  /* Read vector entries. */
  for (i = 0; i < nr; ++i) {
    fscanf(f, "%lg\n", &val);
    b[i] = val;
  }

  /* Close the opened file. */
  fclose(f);

  /* Return vector. */
  return b;
}

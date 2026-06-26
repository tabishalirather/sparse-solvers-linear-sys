/*! \file
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"


int main() {

  csr *A = NULL;
  char *file_mtx;
  char *file_mtx_spy;

  
  file_mtx = "../../TestMatrixCollection/orsirr_1.mtx";
  file_mtx_spy = "orsirr_1.eps";
  
  printf("Matrix from Matrix Market\n");
  printf("---------------------------------------------\n");
  A = csr_import_from_mm_mtx(file_mtx);
  if (A == NULL) {
    printf("Matrix import failed\n");
  }
  csr_spy(A, file_mtx_spy, 's');

  csr_free(A);
  
  return 0;
}

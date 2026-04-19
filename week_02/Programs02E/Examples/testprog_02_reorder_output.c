/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 *
*/

#include <stdlib.h>
#include "itsolver.h"

int main(void)
{
  csr *A;

  A = csr_testmatrices(0);
  csr_test_ordering(A, "PermTridiag", -1);
  csr_free(A);

  A = csr_testmatrices(1);
  csr_test_ordering(A, "Figure2_10a", 2);
  csr_free(A);

  A = csr_testmatrices(1);
  csr_test_ordering(A, "Figure2_10b", 0);
  csr_free(A);

  A = csr_testmatrices(2);
  csr_test_ordering(A, "Figure3_4a", 0);
  csr_free(A);

  A = csr_testmatrices(2);
  csr_test_ordering(A, "Figure3_4b", -1);
  csr_free(A);

  A = csr_testmatrices(3);
  csr_test_ordering(A, "disjunct", -1);
  csr_free(A);

  return EXIT_SUCCESS;
} /* main */

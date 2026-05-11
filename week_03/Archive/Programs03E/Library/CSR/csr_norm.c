
#include "itsolver.h"

#include <assert.h>
#include <math.h>

double csr_norm_infty(const csr *A)
{
  int i, p;
  double absrowsum, normA = 0.0;
  for (i = 0; i < A->nr; i++)
  {
    absrowsum = 0.0;
    for (p = A->row[i]; p < A->row[i + 1]; ++p)
    {
      absrowsum += fabs(A->coeff[p]);
    }
    if (normA < absrowsum) {
      normA = absrowsum;
    }
  }
  return normA;
}

double csr_norm_frob(const csr *A) {
  double val = 0.0;
  int i;
  for (i = 0; i < A->nr; ++i) {
    int p;
    for (p = A->row[i]; p < A->row[i + 1]; ++p) {
      double entry = A->coeff[p];
      val += entry * entry;
    }
  }
  return sqrt(val);
}

/**
 * \file csr_tridiagonal.c
 *
 * \author Jonas Grams
 *
 */

#include "itsolver.h"
#include <assert.h>
#include <math.h>


csr *csr_tridiagonal(int n, double alpha)
{
    csr *A;

    int nnz;

    // First and last row with two entries,
    // all remaining rows three entries.
    nnz = 3 * (n - 2) + 4;

    A = csr_alloc(n, n, nnz);

    double adiag, aoff;
    adiag = alpha;
    aoff = -1.0;

    A->row[0] = 0;
    for (int i = 0; i < n; i++)
    {
        int p = A->row[i];

        // Left off-diagonal, except for first row
        if (i - 1 >= 0)
        {
            A->col[p] = i - 1;
            A->coeff[p] = aoff;
            p++;
        }

        // Diagonal entry
        A->col[p] = i;
        A->coeff[p] = adiag;
        p++;

        // Right off-diagonal entry, except for last row
        if (i + 1 < n)
        {
            A->col[p] = i + 1;
            A->coeff[p] = aoff;
            p++;
        }

        A->row[i + 1] = p;
    }

    assert(nnz == A->row[n]);
    return A;
}

double tridiagonal_eigval_min(int n, double alpha)
{
    return alpha - 2 * cos(M_PI / (n + 1));
}

double tridiagonal_eigval_max(int n, double alpha)
{
    int n1 = n - 1;

    return alpha - 2 * cos(n1 * M_PI / (n + 1));
}
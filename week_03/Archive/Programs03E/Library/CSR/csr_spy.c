/*! \file
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "itsolver.h"

#define max(A,B) ((A) > (B) ? (A) : (B))

void csr_spy(csr * A, const char *filename, const char style)
{
  const int *Arow;
  const int *Acol;
  const double *coeff;
  int nr, nc, nz;
  FILE *out;
  double minval, maxval, scale, translate;
  int i, r;

  assert(A != NULL);

  printf("\n");
  printf("Matrix output: %d rows and %d columns, nnz = %d.\n",
     A->nr,A->nc,A->nz);
  printf("Average of %2.2f non-zero entries per row.\n",
    (1.0*A->nz)/(1.0*A->nr));
  printf("Portion of non-zero entries: %2.2f%%.\n",
    (100.0*A->nz)/(1.0*A->nr*A->nc));
  printf("\n");

  nr = A->nr;
  nc = A->nc;
  nz = A->nz;
  Arow = A->row;
  Acol = A->col;
  coeff = A->coeff;

  minval = maxval = 0.0;
  if (nz > 0)
  {
    minval = maxval = coeff[0];
    for (i=1; i<nz; i++)
    {
      if (minval > coeff[i]) minval = coeff[i];
      if (maxval < coeff[i]) maxval = coeff[i];
    }
  }

  out = fopen(filename, "w");
  assert(out != NULL);

  scale = 500.0 / max(nr, nc);
  switch (style)
  {
    case 's': translate = 0.1;
      (void) fprintf(out,
        "%%!PS-Adobe-2.0 EPSF-2.0\n"
        "%%%%BoundingBox: 0 0 %d %d\n"
        "%f dup scale\n"
        "%f dup translate\n"
        "0 setgray\n"
        "/bx {dup 0 exch 1 exch sub setrgbcolor "
        "moveto 0.8 0.0 "
        "rlineto 0.0 0.8 "
        "rlineto -0.8 0.0 "
        "rlineto closepath fill} def\n",
        (int) (scale * nc), (int) (scale * nr), scale, translate);
    break; /* squares */
    case 'c': translate = 0.5;
      (void) fprintf(out,
        "%%!PS-Adobe-2.0 EPSF-2.0\n"
        "%%%%BoundingBox: 0 0 %d %d\n"
        "%f dup scale\n"
        "%f dup translate\n"
        "0 setgray\n"
        "/bx {dup 0 exch 1 exch sub setrgbcolor "
        "0.45 0 360 arc closepath fill} def\n",
        (int) (scale * nc), (int) (scale * nr), scale, translate);
    break; /* circles */
  }

  scale = maxval-minval;
  if (scale != 0.0) scale = 1.0/scale;

  for (i=0; i<nr; i++)
  {
    for (r=Arow[i]; r<Arow[i+1]; r++)
    {
      (void) fprintf(out,"%d %d %f bx\n", Acol[r], nr-i-1,
         (coeff[r]-minval)*scale);
    }
  }

  (void) fprintf(out,"showpage\n");
  (void) fclose(out);
}

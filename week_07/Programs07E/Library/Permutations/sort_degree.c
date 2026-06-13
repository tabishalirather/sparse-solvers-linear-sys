/*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <stdlib.h>
#include "itsolver.h"

/* sort n entries of vector v by increasing degree (stored in deg) */

void sort_degree(int *v, int *deg, int n)
{
  int *b;
  int k,l,tmp;

  b = (int *) malloc(n*sizeof(int));

  b[0] = 0;
  for (k = 1; k<n; k++)
  {
    if (deg[k]==deg[k-1]) b[k] = b[k-1]+1;
    else if (deg[k]>deg[k-1]) b[k] = 0;
    else
    {
      l = k-1;
      while ((deg[l+1]<deg[l]) && (l>=0))
      {
        tmp = deg[l+1];
        deg[l+1] = deg[l];
        b[l+1] = b[l];
        deg[l-b[l]] = tmp;
        tmp = v[l+1];
        v[l+1] = v[l];
        v[l-b[l]] = tmp;
        l = l-b[l]-1;
      }
    }
  } /* for k */

  free(b);


} /* sort_degree */


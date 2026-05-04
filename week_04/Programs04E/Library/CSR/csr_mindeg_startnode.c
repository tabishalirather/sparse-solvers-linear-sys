/*
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 * 
*/
#include "itsolver.h"

int csr_mindeg_startnode(csr *A)
{
  int i, start = 0;
  int deg, mindeg, Arow, Arownew;

  Arow = *A->row++;
  Arownew = *A->row++;
  deg = Arownew-Arow;
  mindeg = deg;

  for (i=1; i<A->nr; i++)
  {
    Arow = Arownew;
    Arownew = *A->row++;
    deg = Arownew-Arow;
    if (mindeg > deg)
    {
      mindeg = deg;
      start = i;
    }
  }
  A->row -= A->nr+1;

  return start;

} /* csr_mindeg_startnode */

/*! \file
 *  \brief Test program showing how to multiply a sparse matrix and a vector.
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
  csr A;

  /*
         / 0 3 1 0 0 \
     A = | 1 0 4 0 5 |
         \ 0 3 2 4 0 /
   */

  int row[4] = {0,2,5,8};
  int col[8] = {1,2,0,2,4,1,2,3};
  double coeff[8] = {3.,1.,1.,4.,5.,3.,2.,4.};

  double x[5] = {1.,2.,3.,4.,5.};
  double y[3] = {1.,2.,3.};

  A.nr = 3;
  A.nc = 5;
  A.nz = 8;
  A.row = row;
  A.col = col;
  A.coeff = coeff;

  /* testing the matrix-vector product routine */
  printf("\n");
  printf("  y = [%3g;%3g;%3g]\n",y[0],y[1],y[2]);
  printf("  x = [%3g;%3g;%3g;%3g;%3g]\n\n",x[0],x[1],x[2],x[3],x[4]);
  printf("Test of the matrix-vector operation y := y + A*x: ");
  /* y := y + A*x */
  csr_mv(&A, 1.0, x, 1.0, y);
  printf("  y = [%3g;%3g;%3g]\n",y[0],y[1],y[2]);
  printf("Test of the matrix-vector operation y := -y + A*x: ");
  /* y := -y + A*x  */
  csr_mv(&A, 1.0, x, -1.0, y);
  printf("  y = [%3g;%3g;%3g]\n",y[0],y[1],y[2]);
  printf("\n");

  /* setting back -y to y */
  csr_mv(&A, 0.0, x, -1.0, y);

  /* testing the vector-matrix product routine */
  printf("\n");
  printf("  y = [%3g;%3g;%3g]\n",y[0],y[1],y[2]);
  printf("  x = [%3g;%3g;%3g;%3g;%3g]\n",x[0],x[1],x[2],x[3],x[4]);
  printf("Test of the matrix-vector operation x' = x' + y'*A: ");
  /* x' = x' + y'*A */
  csr_vm(&A, 1.0, y, 1.0, x);
  printf("  x = [%3g;%3g;%3g;%3g;%3g]\n",x[0],x[1],x[2],x[3],x[4]);
  printf("Test of the matrix-vector operation x' = -x' + y'*A: ");
  /* x' = -x' + y'*A */
  csr_vm(&A, 1.0, y, -1.0, x);
  printf("  x = [%3g;%3g;%3g;%3g;%3g]\n",x[0],x[1],x[2],x[3],x[4]);
  printf("\n");
  printf("This is test print from mv.c");

  return EXIT_SUCCESS;
} /* main */

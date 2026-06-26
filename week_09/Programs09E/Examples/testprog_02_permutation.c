/*
 *
 *  \author Sabine Le Borne
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include "itsolver.h"

int main(void)
{
  /*

    Test the permutation routines with

           / 0 0 1 \
       P = | 1 0 0 | = [e_2,e_3,e_1] <=> p = [2 3 1].
           \ 0 1 0 /

                       / e_3^T \
    pinv = [3 1 2]^T = | e_1^T |.
                       \ e_2^T /

    In C, p -> [1 2 0] = perm.

   */

  csr *A, *PA, *AP, *PTA, *APT, *PAP, *PTAP, *PAPT, *PTAPT;
  int Arow[4] = {0,3,6,9};
  int Acol[9] = {0,1,2,0,1,2,0,1,2};
  double Acoeff[9] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};

  double b[3] = {1.0,-2.0,1.0};
  double *Pb, *PTb;

  int perm[3] = {1,2,0};
  int *iperm;

  int i,j;

  /* set up A */
  A = csr_alloc(3,3,9);

  for (i=0; i<4; i++) A->row[i] = Arow[i];
  for (i=0; i<9; i++) A->col[i] = Acol[i];
  for (i=0; i<9; i++) A->coeff[i] = Acoeff[i];

  /* compute inverse permutation */
  iperm = perminv(3, perm);

  /* compute permuted vectors */
  Pb  = vec_permute(3, b, iperm);
  PTb = vec_permute(3, b, perm);

  /* produce color plots of all permuted matrices */
  csr_spy(A,"A.eps",'s');

  PA    = csr_permute(A, iperm, NULL);
  csr_spy(PA, "PA.eps", 's');

  AP    = csr_permute(A, NULL, perm);
  csr_spy(AP, "AP.eps", 's');

  PTA   = csr_permute(A, perm, NULL);
  csr_spy(PTA, "PTA.eps", 's');

  APT   = csr_permute(A, NULL, iperm);
  csr_spy(APT, "APT.eps", 's');

  PAP   = csr_permute(A, iperm, perm);
  csr_spy(PAP,"PAP.eps",'s');

  PAPT  = csr_permute(A, iperm, iperm);
  csr_spy(PAPT, "PAPT.eps", 's');

  PTAP  = csr_permute(A, perm, perm);
  csr_spy(PTAP, "PTAP.eps", 's');

  PTAPT = csr_permute(A, perm, iperm);
  csr_spy(PTAPT, "PTAPT.eps", 's');

  /* output of permutation vectors */
  printf("Output of the permutation vectors:\n");
  printf("==================================\n");

  printf(" perm =");
  for (i=0;i<3; i++)
    printf(" %d", perm[i]+1);
  printf("\n");

  printf("iperm =");
  for (i=0;i<3; i++)
    printf(" %d", iperm[i]+1);
  printf("\n");

  /* output of permuted matrices */
  printf("\n");
  printf("Output of the permuted matrices:\n");
  printf("================================\n");

  printf("A =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", A->coeff[i*3+j]);

    printf("\n");
  }
  printf("\n");

  printf("PA =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", PA->coeff[i*3+j]);
   
    printf("\n");
  }
  printf("\n");

  printf("AP =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", AP->coeff[i*3+AP->col[j]]);
   
    printf("\n");
  }
  printf("\n");

  printf("PTA =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", PTA->coeff[i*3+j]);
   
    printf("\n");
  }
  printf("\n");

  printf("APT =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", APT->coeff[i*3+APT->col[j]]);
   
    printf("\n");
  }
  printf("\n");

  printf("PAP =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", PAP->coeff[i*3+PAP->col[j]]);
   
    printf("\n");
  }
  printf("\n");

  printf("PAPT =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", PAPT->coeff[i*3+PAPT->col[j]]);
   
    printf("\n");
  }
  printf("\n");

  printf("PTAP =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", PTAP->coeff[i*3+PTAP->col[j]]);
   
    printf("\n");
  }
  printf("\n");

  printf("PTAPT =\n");
  for (i=0;i<3; i++)
  {
    for (j=0;j<3; j++)
      printf(" %g", PTAPT->coeff[i*3+PTAPT->col[j]]);
   
    printf("\n");
  }
  printf("\n");

  printf("Output of the permuted vectors:\n");
  printf("===============================\n");

  printf("b   =");
  for (i=0;i<3; i++)
    printf(" %+g", b[i]);
  printf("\n");

  printf("Pb  =");
  for (i=0;i<3; i++)
    printf(" %+g", Pb[i]);
  printf("\n");

  printf("PTb =");
  for (i=0;i<3; i++)
    printf(" %+g", PTb[i]);
  printf("\n\n");

  free(iperm);

  csr_free(A);
  csr_free(PA);
  csr_free(AP);
  csr_free(PTA);
  csr_free(APT);
  csr_free(PAP);
  csr_free(PAPT);
  csr_free(PTAP);
  csr_free(PTAPT);

  free(Pb);
  free(PTb);

  return EXIT_SUCCESS;
} /* main */

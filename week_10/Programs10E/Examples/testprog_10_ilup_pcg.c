/*! \file
 *  \brief Number of PCG iterations with an ILU(p) preconditioner for different values of p.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

int main(int argc, char **argv) {

  /* Discretisation parameter. */
  int n = 200;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  printf("n = %d\n", n);

  /* Spatial dimension. */
  int d = 2;
  if (argc > 2) {
    d = atoi(argv[2]);
  }
  printf("d = %d\n", d);

  /* ILU(p) will be used with p = 0, ..., p_max. */
  int p_max = 7;
  if (argc > 3) {
    p_max = atoi(argv[3]);
  }
  printf("p_max = %d\n", p_max);

  /* Testmatrix. */
  csr *A = csr_laplace(n, d);
  int m = A->nr;

  /* Exact solution of all ones and corresponding right hand side. */
  double *x0 = vec_zeros(m);
  double *xstar = vec_ones(m);
  double *b = vec_zeros(m);
  csr_mv(A, 1.0, xstar, 0.0, b);

  /* Iteration paramters. */
  double *x = vec_raw(m);
  int stopcrit = ITSOLVER_RELRES;
  double tol = 1e-6;
  int maxit = 1000;
  double esterr;
  int iter;
  double *errvec = malloc((maxit + 1) * sizeof(double));

  /* CG */
  printf("CG:\n");
  vec_copy(m, x0, x);

  clock_t cg_time = clock();
  csr_cg(A, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
  cg_time = clock() - cg_time;
  printf("   CG time      : %7.3f ms\n", (double) cg_time / CLOCKS_PER_SEC * 1000);

  printf("   CG iterations: %d\n", iter);

  /* PCG with ILU(p) and p = 0, 1, ... */
  int p;
  for (p = 0; p <= p_max; ++p) {
    /* PCG. */
    printf("PCG with ILU(p), p = %d:\n", p);
    vec_copy(m, x0, x);

    clock_t ilup_time = clock();
    prec_t *prec = prec_ilup_new(A, p);
    ilup_time = clock() - ilup_time;
    printf("     ILU(p) time: %7.3f ms\n", (double) ilup_time / CLOCKS_PER_SEC * 1000);

    clock_t pcg_time = clock();
    csr_pcg(A, prec, b, x, tol, stopcrit, maxit, &esterr, &iter, errvec);
    pcg_time = clock() - pcg_time;
    printf("        PCG time: %7.3f ms\n", (double) pcg_time / CLOCKS_PER_SEC * 1000);

    printf("      Total time: %7.3f ms\n", (double) (ilup_time + pcg_time) / CLOCKS_PER_SEC * 1000);

    prec_delete(prec);
    printf("  PCG iterations: %d\n", iter);
  }

  /* Free memory. */
  csr_free(A);
  free(x0);
  free(xstar);
  free(b);
  free(x);
  free(errvec);

  return 0;
}

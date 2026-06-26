/* Autor: Sabine Le Borne */

#include "itsolver.h"

int main(int argc, char **argv) {

  double lambda_min, lambda_max, kappa;
  int i;
  double omega_opt;
  double *b, *x, *xstart;
  double erra_start;
  double erra_end_ri, erra_relred_ri;
  double erra_end_sd, erra_relred_sd;
  double erra_end_cg, erra_relred_cg;
  double rho_ri, rho_sd, rho_cg;
  double tol;
  int stopcrit;
  int maxit[3];
  double esterr;
  int iter[3];
  double *errvec[3];
  char *legend_entries[3];

  /* Input parameters. */
  int n = 20;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  printf("n = %d\n", n);
  int d = 2;
  if (argc > 2) {
    d = atoi(argv[2]);
  }
  printf("d = %d\n", d);

  /* Laplace model problem. */
  csr *A = csr_laplace(n, d);
  int nd = A->nc;
  printf("Memory CSR (coefficients only): %.3e\n", (double) laplace_nnz(n, d) * 8);

  /* Compute optimal relaxation parameter for Richardson. */
  lambda_min = 4 * d * pow(sin(          M_PI / (2 * n)), 2);
  lambda_max = 4 * d * pow(sin((n - 1) * M_PI / (2 * n)), 2);
  omega_opt = 2.0 / (lambda_min + lambda_max);
  kappa = lambda_max / lambda_min;
  rho_ri = (kappa - 1.0) / (kappa + 1.0);
  rho_sd = rho_ri;
  rho_cg = (sqrt(kappa) - 1.0) / (sqrt(kappa) + 1.0);
  printf("n = %d, d = %d, nd = %d\n", n, d, nd);

  /* Zero right hand side. */
  b = (double *) malloc(nd * sizeof(double));
  vec_scale(nd, 0.0, b);

  /* All ones start vector. */
  x = (double *) malloc(nd * sizeof(double));
  xstart = (double *) malloc(nd * sizeof(double));
  for (i = 0; i<nd; i++)
    xstart[i] = 1.0;
  
  /* Start error in A-norm. */
  erra_start = sqrt(csr_xtAy(A, xstart, xstart));
  printf("Start error in A-Norm: %g\n", erra_start);

  /* Tolerance and stopping criterium. */
  tol = 1e-16;
  stopcrit = ITSOLVER_RELRES;

  /* Richardson's iteration with optimal omega. */
  printf("\n");
  printf("Richardson\n");
  printf("----------\n");
  maxit[0] = 1000;
  errvec[0] = (double *) malloc((maxit[0] + 1) * sizeof(double));
  memcpy(x, xstart, nd * sizeof(double));
 
  legend_entries[0] = "Richardson";
  csr_richardson(A, b, x, tol, stopcrit, maxit[0], omega_opt, &esterr,
		 &iter[0], errvec[0]);
		 
  /* Error reduction in A-norm. */
  erra_end_ri = sqrt(csr_xtAy(A, x, x));

  erra_relred_ri = erra_end_ri / erra_start;
  printf("Theoretical convergence rate (log10): %g\n", log10(rho_ri));
  printf("End error in A-Norm: %g\n", erra_end_ri);
  printf("Number of iteration steps: %d\n", iter[0]);
  printf("Relative error reduction in A-Norm: %g\n", erra_relred_ri);

  /* Steepest descent. */
  printf("\n");
  printf("Steepest Descent\n");
  printf("----------------\n");
  maxit[1] = 1000;
  errvec[1] = (double *) malloc((maxit[1] + 1) * sizeof(double));
  memcpy(x, xstart, nd * sizeof(double));
  
  legend_entries[1] = "Steepest descent";
  csr_steepest_descent(A, b, x, tol, stopcrit, maxit[1], &esterr,
		 &iter[1], errvec[1]);
  
  /* Error reduction in A-norm. */
  erra_end_sd = sqrt(csr_xtAy(A, x, x));

  erra_relred_sd = erra_end_sd / erra_start;
  printf("Theoretical convergence rate (log10): %g\n", log10(rho_sd));
  printf("End error in A-Norm: %g\n", erra_end_sd);
  printf("Number of iteration steps: %d\n", iter[1]);
  printf("Relative error reduction in A-Norm: %g\n", erra_relred_sd);

  /* Conjugate gradients */
  printf("\n");
  printf("CG\n");
  printf("--\n");
  maxit[2] = 1000;
  errvec[2] = (double *) malloc((maxit[2] + 1) * sizeof(double));
  memcpy(x, xstart, nd * sizeof(double));
  
  csr_cg(A, b, x, tol, stopcrit, maxit[2], &esterr,
		 &iter[2], errvec[2]);
  
  legend_entries[2] = "CG";

  /* Error reduction in A-norm. */
  erra_end_cg = sqrt(csr_xtAy(A, x, x));

  erra_relred_cg = erra_end_cg / erra_start;
  printf("Theoretical convergence rate (log10): %g\n", log10(rho_cg));
  printf("End error in A-Norm: %g\n", erra_end_cg);
  printf("Number of iteration steps: %d\n", iter[2]);
  printf("Relative error reduction in A-Norm: %g\n", erra_relred_cg);

  /* Combined plot of all the methods. */
  convplot_steps(
    "testprog_06_cg.eps",
    "CG method applied to the discretization of the Laplace operator",
    "Iteration step",
    "Residual norm",
    3,
    iter,
    errvec,
    legend_entries);

  /* Individual plot for CG. */
  char *file_name = malloc(100);
  char *title = malloc(100);
  snprintf(file_name, 100, "testprog_06_cg_only_n%d_d%d.eps", n, d);
  snprintf(title, 100, "Laplace model problem with n = %d, d = %d", n, d);
  convplot_steps(
    file_name,
    title,
    "Iteration step",
    "Residual norm",
    1,
    iter + 2,
    errvec + 2,
    legend_entries + 2);
  free(file_name);
  free(title);

  /* Free memory. */
  free(errvec[2]);
  free(errvec[1]);
  free(errvec[0]);
  free(xstart);
  free(x);
  free(b);
  csr_free(A);

  return 0;
}

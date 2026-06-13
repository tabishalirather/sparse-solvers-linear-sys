
#include "itsolver.h"

gmres_t *gmres_new(csr *A, double *b, double *x0)
{
  // Assert that A is square and b is given.
  // If x0 is not given the zero vector will be used.
  assert(A);
  assert(A->nr == A->nc);
  assert(b);

  // Initialize GMRes stucture.
  gmres_t *gmres = malloc(sizeof(gmres_t));
  gmres->n  = A->nr;
  gmres->A  = A;
  gmres->b  = b; 

  // Set x0 to zero if not specified.
  gmres->x0 = calloc(gmres->n, sizeof(double));
  if (x0) {
    vec_copy(gmres->n, x0, gmres->x0);
  }

  // Copy x0 to x.
  gmres->x = malloc(gmres->n * sizeof(double));
  vec_copy(gmres->n, gmres->x0, gmres->x);

  gmres->iter = 0;  // Current number of total iterations performed.
  gmres->iter_max = 1; // One iteration only.
  gmres->restart = -1; // Meaning no restart.
  gmres->tol = 1e-8;
  gmres->residual_norms = NULL; // Will be allocated in gmres_solve.
  gmres->n_ortho = 2; // Orthogonalize twice by default.

  return gmres;
}

void gmres_iter_max(gmres_t *gmres, int iter_max)
{
  assert(gmres);
  assert(iter_max >= 0);
  gmres->iter_max = iter_max;
}

void gmres_restart(gmres_t *gmres, int restart)
{
  assert(gmres);
  if (restart < 0) {
    gmres->restart = 0;
  } else {
    gmres->restart = restart;
  }
}

void gmres_tol(gmres_t *gmres, double tol)
{
  assert(gmres);
  assert(tol >= 0.0 || tol <= 1.0);
  gmres->tol = tol;
}

static void gmres_compute_x(gmres_t *gmres, int m, int i, double *h, double *v, double *y)
{
  assert(m >= 0);
  assert(i >= 0 && i <= m);

  // Set x = x0.
  vec_copy(gmres->n, gmres->x0, gmres->x);

  // Update x += v * z where h * z = y.
  if (i > 0) {

    // Solve h * z = y for z.
    double *z = malloc(i * sizeof(double));
    int j;
    for (j = i - 1; j >= 0; --j) {
      z[j] = y[j];
      int k;
      for (k = i - 1; k > j; --k) {
        z[j] -= h[j + k * (m + 1)] * z[k];
      }
      z[j] /= h[j + j * (m + 1)];
    }

    // Update x += v * z.
    for (j = 0; j < i; ++j) {
      vec_axpy(gmres->n, z[j], v + j * gmres->n, gmres->x);
    }
    free(z);
  }
}

void gmres_solve(gmres_t *gmres)
{
  // Do not expect that this function is called more than once.
  assert(gmres->iter == 0);

  // Flag used in the outer loop to determine whether the algorithm should be
  // stopped indefinitely or to continue with a restart.
  bool full_stop;

  // Determine the maximum number of inner iterations within each outer
  // iteration, i.e. before a restart is performed.
  // At most m + 1 basis vectors are stored at any time.
  int m = gmres->iter_max;
  if (gmres->restart > 0) { // Non-postive value means no restart.
    if (gmres->restart < m) { // Assure that m will be no larger than iter_max.
      m = gmres->restart;
    }
  }
  printf("m = %d\n", m);

  // Vector of residual norms.
  gmres->residual_norms = calloc(gmres->iter_max + 1, sizeof(double));

  // Outer loop.
  while (true) {

    // Hessenberg matrix which is transformed into an upper triangular matrix by
    // applying the necessary Givens rotations immediately after each step.
    double *h = calloc((m + 1) * m, sizeof(double));

    // Matrix of up to m + 1 basis vectors.
    double *v = malloc(gmres->n * (m + 1) * sizeof(double));

    // Copy current solution x to initial solution x0 after a restart.
    if (gmres->iter > 0) {
      vec_copy(gmres->n, gmres->x, gmres->x0);
    }

    // Compute the initial residual at the beginning of the inner loop.
    vec_copy(gmres->n, gmres->b, v);
    csr_mv(gmres->A, -1.0, gmres->x, 1.0, v);

    // Norm of the initial residual.
    double beta = vec_norm(gmres->n, v);
    vec_scale(gmres->n, 1.0 / beta, v);

    // Right hand side of the least squares problem.
    double *y = calloc(m + 1, sizeof(double));
    y[0] = beta;

    // Givens rotations.
    double *c = malloc(m * sizeof(double)); // cosine
    double *s = malloc(m * sizeof(double)); // sine

    // Inner iteration index.
    int i = 0;

    // Inner loop.
    while (true) {

      // Make sure that the last residual norm before a restart is not stored
      // a second time as the first residual norm of the next outer iteration.
      if (gmres->iter == 0 || i > 0) {

        // Estimate of the residual norm obtained from the last entry of
        // the modified right hand side of the least squares problem.
        gmres->residual_norms[gmres->iter] = fabs(y[i]);
      }

      // Check whether the desired accuracy or the maximum number of iterations
      // have been reached. In this case there will be no restart and setting
      // this flag to true will cause the outer loop to terminate as well.
      full_stop = true;
      if (gmres->iter >= gmres->iter_max ||
          gmres->residual_norms[gmres->iter] <= gmres->tol * gmres->residual_norms[0]) {
        break;
      }

      // Check for a restart and set the flag correspondingly.
      if (i >= m) {
        full_stop = false;
        break;
      }

      // Matrix-vector-multiplication.
      csr_mv(gmres->A, 1.0, v + i * gmres->n, 0.0, v + (i + 1) * gmres->n);

      { // Orthogonalize a number of times using Modified Gram-Schmidt (MGS).
        int ortho;
        for (ortho = 0; ortho < gmres->n_ortho; ++ortho) {
          int j;
          for (j = 0; j <= i; ++j) {
            // Compute dot product.
            double h_update = vec_dot(gmres->n, v + j * gmres->n, v + (i + 1) * gmres->n);
            // Update coefficient of Hessenberg matrix.
            h[j + i * (m + 1)] += h_update;
            // Substract from new basis vector.
            vec_axpy(gmres->n, -h_update, v + j * gmres->n, v + (i + 1) * gmres->n);
          }
        }
      }

      // Normalize the new basis vector.
      beta = vec_norm(gmres->n, v + (i + 1) * gmres->n);
      h[i + 1 + i * (m + 1)] = beta;
      vec_scale(gmres->n, 1.0 / beta, v + (i + 1) * gmres->n);

      { // Apply all previous Givens rotations to the new column of the Hessenberg matrix.
        int j;
        for (j = 0; j < i; ++j) {
          // The leading dimensions/ increments are irrelevant since
          // the rotation is applied to a single column only.
          givens_rot(1, h + j + i * (m + 1), m + 1, h + j + 1 + i * (m + 1), m + 1, c[j], s[j]);
        }
      }

      // Compute new Givens rotation.
      double *h0 = h + i + i * (m + 1);
      double *h1 = h0 + 1;

      // The norm of [h0 h1], i.e. r = sqrt(h0 ** 2 + h1 ** 2) is stored in h0.
      givens_rotg(h0, h1, c + i, s + i);

      // Apply the newly computed Givens rotation to the right hand side.
      givens_rot(1, y + i, m + 1, y + i + 1, m + 1, c[i], s[i]);

      // Increment counts for the total number of inner iterations and for the
      // number of inner iterations within the current outer iteration.
      ++gmres->iter;
      ++i;
    }

    // Compute solution.
    gmres_compute_x(gmres, m, i, h, v, y);

    // Free memory.
    free(h);
    free(v);
    free(y);
    free(c);
    free(s);

    // Quit or restart.
    if (full_stop) {
      break;
    }
  }
}

double *gmres_solution(gmres_t *gmres)
{
  assert(gmres);
  return gmres->x;
}

int gmres_iter(gmres_t *gmres)
{
  assert(gmres);
  return gmres->iter;
}

double *gmres_residual_norms(gmres_t *gmres)
{
  assert(gmres);
  return gmres->residual_norms;
}

void gmres_delete(gmres_t *gmres)
{
  assert(gmres);
  // Vector of residual norms may be null.
  free(gmres->residual_norms);
  assert(gmres->x0);
  free(gmres->x0);
  assert(gmres->x);
  free(gmres->x);
  free(gmres);
}

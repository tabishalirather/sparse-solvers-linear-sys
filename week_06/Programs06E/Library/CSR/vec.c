
#include "itsolver.h"

double *vec_rand(int n)
{
  double *x = malloc(n * sizeof(double));
  int i;
  for (i = 0; i < n; ++i) {
    x[i] = rand() / (RAND_MAX + 1.0);
  }
  return x;
}

double *vec_raw(int n)
{
  double *x = malloc(n * sizeof(double));
  return x;
}

double *vec_ones(int n)
{
  double *e = malloc(n * sizeof(double));
  int i;
  for (i = 0; i < n; ++i) {
    e[i] = 1.0;
  }
  return e;
}

double *vec_zeros(int n)
{
  double *o = calloc(n, sizeof(double));
  return o;
}

void vec_copy(int n, const double *x, double *y)
{
#if USE_BLAS
  int one = 1;
  dcopy_(&n, (double *) x, &one, y, &one);
#else
  int i;
  for (i = 0; i < n; ++i) {
    y[i] = x[i];
  }
#endif
}

void vec_print(int n, const double *x)
{
  int i;
  for (i = 0; i < n; ++i) {
    printf("%e\n", x[i]);
  }
}

double vec_norm(int n, const double *x)
{
#if USE_BLAS
  int one = 1;
  return dnrm2_(&n, (double *) x, &one);
#else
  double val = 0;
  int i;
  for (i = 0; i < n; ++i) {
    val += x[i] * x[i];
  }
  return sqrt(val);
#endif
}

double vec_norm_inf(int n, const double *x)
{
  double max = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    double abs = fabs(x[i]);
    if (abs > max) {
      max = abs;
    }
  }
  return max;
}

double vec_dot(int n, const double *x, const double *y)
{
#if USE_BLAS
  int one = 1;
  return ddot_(&n, (double *) x, &one, (double *) y, &one);
#else
  double val = 0;
  int i;
  for (i = 0; i < n; ++i) {
    val += x[i] * y[i];
  }
  return val;
#endif
}

void vec_scale(int n, double alpha, double *x)
{
#if USE_BLAS
  int one = 1;
  dscal_(&n, &alpha, x, &one);
#else
  int i;
  for (i = 0; i < n; ++i) {
    x[i] *= alpha;
  }
#endif
}

void vec_axpy(int n, double alpha, const double *x, double *y)
{
#if USE_BLAS
  int one = 1;
  daxpy_(&n, &alpha, (double *) x, &one, y, &one);
#else
  int i;
  for (i = 0; i < n; ++i) {
    y[i] += alpha * x[i];
  }
#endif
}



#include "itsolver.h"

#if USE_BLAS
void drotg_(double *a, double *b, double *c, double *s);
#endif

void givens_rotg(double *x, double *y, double *c, double *s)
{
#if USE_BLAS
  drotg_(x, y, c, s);
#else
  double r = sqrt(*x * *x + *y * *y);
  *c = *x / r;
  *s = *y / r;
  *x = r;
  if (fabs(*c) > fabs(*s)) {
    *y = *s;
  } else {
    *y = (*c != 0.0) ? 1.0 / *c : 1.0;
  }
#endif
}

#if USE_BLAS
void drot_(const int *N, double *X, const int *incX,
           double *Y, const int *incY, const double *c, const double *s);
#endif

void givens_rot(int n, double *x, int incx, double *y, int incy, double c, double s)
{
#if USE_BLAS
  drot_(&n, x, &incx, y, &incy, &c, &s);
#else
  int i;
  for (i = 0; i < n; ++i) {
    double x_tmp = x[i * incx];
    double y_tmp = y[i * incy];
    x[i * incx] =  c * x_tmp + s * y_tmp;
    y[i * incx] = -s * x_tmp + c * y_tmp;
  }
#endif
}

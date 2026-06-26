/*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <stdio.h>
#include <stdlib.h>
#include "itsolver.h"

inline double
csr_steepest_descent_step(csr *A, double *r, double *p, double *x)
{
  int n;
  double alpha;

  n = A->nc;

  alpha = vec_dot(n, r, r) / vec_dot(n, p, r);
  vec_axpy(n, alpha, r, x);
  alpha = -alpha;
  vec_axpy(n, alpha, p, r);
  csr_mv(A, 1.0, r, 0.0, p);

  return vec_norm(n, r);
}

void csr_steepest_descent(csr *A, double *b, double *x,
                          double tol, int stopcrit, int maxit,
                          double *esterr, int *iter, double *errvec)
{
  int i, n;
  double *r, *p;
  double rho;
  double normA, normb;

  n = A->nc;

  switch (stopcrit)
  {
  case ITSOLVER_NORES:
    r = (double *)malloc(n * sizeof(double));
    p = (double *)calloc(n, sizeof(double));
    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    csr_mv(A, 1.0, r, 0.0, p);
    for (i = 0; i < maxit; i++)
    {
      csr_steepest_descent_step(A, r, p, x);
    }
    free(r);
    free(p);
    break;
  case ITSOLVER_RELRES:
    r = (double *)malloc(n * sizeof(double));
    p = (double *)calloc(n, sizeof(double));

    /* compute zeroth residual */

    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    rho = vec_norm(n, r);
    *errvec++ = rho;
    tol = tol * rho;
    *iter = maxit;
    csr_mv(A, 1.0, r, 0.0, p);

    for (i = 0; i < maxit; i++)
    {
      double resnorm;
      resnorm = csr_steepest_descent_step(A, r, p, x);

      /* estimate residual */

      *errvec = resnorm;
      if (*errvec++ <= tol)
      {
        *iter = i + 1;
        break;
      }
    }

    *esterr = *(--errvec) / rho;

    /* no need to move back errvec, is a copy */

    free(r);
    free(p);
    break;
  case ITSOLVER_BACKERR:
    r = (double *)malloc(n * sizeof(double));
    normb = vec_norm(n, b);
    normA = csr_norm_infty(A);
    p = (double *)calloc(n, sizeof(double));

    /* compute zeroth backward error */

    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    *errvec =
        vec_norm(n, r) / (normA * vec_norm(n, x) + normb);
    if (*errvec++ <= tol)
    {
      maxit = 0;
      *iter = 0;
    }
    csr_mv(A, 1.0, r, 0.0, p);

    for (i = 0; i < maxit; i++)
    {
      double resnorm;
      resnorm = csr_steepest_descent_step(A, r, p, x);

      /* compute backward error */

      *errvec = resnorm / (normA * vec_norm(n, x) + normb);
      if (*errvec++ <= tol)
      {
        *iter = i + 1;
        break;
      }
    }

    *esterr = *(--errvec);

    /* no need to move back errvec, is a copy */

    free(r);
    break;
  }
}

inline double
csr_minimal_residual_step(csr *A, double *r, double *p, double *x)
{
  int n;
  double alpha;

  n = A->nc;

  alpha = vec_dot(n, r, p) / vec_dot(n, p, p);
  vec_axpy(n, alpha, r, x);
  alpha = -alpha;
  vec_axpy(n, alpha, p, r);
  csr_mv(A, 1.0, r, 0.0, p);

  return vec_norm(n, r);
}

void csr_minimal_residual(csr *A, double *b, double *x,
                          double tol, int stopcrit, int maxit,
                          double *esterr, int *iter, double *errvec)
{
  int i, n;
  double *r, *p;
  double rho;
  double normA, normb;

  n = A->nc;

  switch (stopcrit)
  {
  case ITSOLVER_NORES:
    r = (double *)malloc(n * sizeof(double));
    p = (double *)calloc(n, sizeof(double));
    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    csr_mv(A, 1.0, r, 0.0, p);
    for (i = 0; i < maxit; i++)
    {
      csr_minimal_residual_step(A, r, p, x);
    }
    free(r);
    free(p);
    break;
  case ITSOLVER_RELRES:
    r = (double *)malloc(n * sizeof(double));
    p = (double *)calloc(n, sizeof(double));

    /* compute zeroth residual */

    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    rho = vec_norm(n, r);
    *errvec++ = rho;
    tol = tol * rho;
    *iter = maxit;
    csr_mv(A, 1.0, r, 0.0, p);

    for (i = 0; i < maxit; i++)
    {
      double resnorm;
      resnorm = csr_minimal_residual_step(A, r, p, x);

      /* estimate residual */

      *errvec = resnorm;
      if (*errvec++ <= tol)
      {
        *iter = i + 1;
        break;
      }
    }

    *esterr = *(--errvec) / rho;

    /* no need to move back errvec, is a copy */

    free(r);
    free(p);
    break;
  case ITSOLVER_BACKERR:
    r = (double *)malloc(n * sizeof(double));
    normb = vec_norm(n, b);
    normA = csr_norm_infty(A);
    p = (double *)calloc(n, sizeof(double));

    /* compute zeroth backward error */

    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    *errvec =
        vec_norm(n, r) / (normA * vec_norm(n, x) + normb);
    *iter = maxit;
    if (*errvec++ <= tol)
    {
      maxit = 0;
      *iter = 0;
    }
    csr_mv(A, 1.0, r, 0.0, p);

    for (i = 0; i < maxit; i++)
    {
      double resnorm;
      resnorm = csr_minimal_residual_step(A, r, p, x);

      /* compute backward error */

      *errvec = resnorm / (normA * vec_norm(n, x) + normb);
      if (*errvec++ <= tol)
      {
        *iter = i + 1;
        break;
      }
    }

    *esterr = *(--errvec);

    /* no need to move back errvec, is a copy */

    free(r);
    break;
  }
}

inline double
csr_residual_norm_steepest_descent_step(csr *A, double *r, double *v, double *w, double *x)
{
  int n;
  double alpha;

  n = A->nc;

  csr_vm(A, 1.0, r, 0.0, v);
  csr_mv(A, 1.0, v, 0.0, w);
  alpha = vec_dot(n, v, v) / vec_dot(n, w, w);
  vec_axpy(n, alpha, v, x);
  alpha = -alpha;
  vec_axpy(n, alpha, w, r);

  return vec_norm(n, r);
}

void csr_residual_norm_steepest_descent(
    csr *A, double *b, double *x,
    double tol, int stopcrit, int maxit,
    double *esterr, int *iter, double *errvec)
{
  int i, n;
  double *r, *v, *w;
  double rho;
  double normA, normb;

  n = A->nc;

  switch (stopcrit)
  {
  case ITSOLVER_NORES:
    r = (double *)malloc(n * sizeof(double));
    v = (double *)calloc(n, sizeof(double));
    w = (double *)calloc(n, sizeof(double));
    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    for (i = 0; i < maxit; i++)
    {
      csr_residual_norm_steepest_descent_step(A, r, v, w, x);
    }
    free(r);
    free(v);
    free(w);
    break;
  case ITSOLVER_RELRES:
    r = (double *)malloc(n * sizeof(double));
    v = (double *)calloc(n, sizeof(double));
    w = (double *)calloc(n, sizeof(double));

    /* compute zeroth residual */

    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    rho = vec_norm(n, r);
    *errvec++ = rho;
    tol = tol * rho;
    *iter = maxit;

    for (i = 0; i < maxit; i++)
    {
      double resnorm;
      resnorm = csr_residual_norm_steepest_descent_step(A, r, v, w, x);

      /* estimate residual */

      *errvec = resnorm;
      if (*errvec++ <= tol)
      {
        *iter = i + 1;
        break;
      }
    }

    *esterr = *(--errvec) / rho;

    /* no need to move back errvec, is a copy */

    free(r);
    free(v);
    free(w);
    break;
  case ITSOLVER_BACKERR:
    r = (double *)malloc(n * sizeof(double));
    normb = vec_norm(n, b);
    normA = csr_norm_infty(A);
    v = (double *)calloc(n, sizeof(double));
    w = (double *)calloc(n, sizeof(double));

    /* compute zeroth backward error */

    vec_copy(n, b, r);
    csr_mv(A, -1.0, x, 1.0, r);
    *errvec =
        vec_norm(n, r) / (normA * vec_norm(n, x) + normb);
    *iter = maxit;
    if (*errvec++ <= tol)
    {
      maxit = 0;
      *iter = 0;
    }

    for (i = 0; i < maxit; i++)
    {
      double resnorm;
      resnorm = csr_residual_norm_steepest_descent_step(A, r, v, w, x);


      /* compute backward error */

      *errvec = resnorm / (normA * vec_norm(n, x) + normb);
      if (*errvec++ <= tol)
      {
        *iter = i + 1;
        break;
      }
    }

    *esterr = *(--errvec);

    /* no need to move back errvec, is a copy */

    free(r);
    free(v);
    free(w);
    break;
  }
}

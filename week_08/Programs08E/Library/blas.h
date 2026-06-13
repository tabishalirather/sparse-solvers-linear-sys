/*! \file
 *  \brief Declarations for a selection of BLAS routines.
 *
 *  See <a href=http://netlib.org/blas>http://netlib.org/blas</a>.
 */

#ifndef BLAS_H
#define BLAS_H

/*! \brief Swap vectors \f$ x \f$ and \f$ y \f$.
 *
 *  Update \f$ x \leftrightarrow y \f$ .
 *
 *    \param N Number of entries in each vector.
 *    \param X First vector.
 *    \param INCX Increment between subsequent entries in \p X.
 *    \param Y Second vector.
 *    \param INCY Increment between subsequent entries in \p Y.
 */
extern void dswap_(int *N, double *X, int *INCX, double *Y, int *INCY);

/*! \brief Copy vector \f$ x \f$ into \f$ y \f$.
 *
 *  Update \f$ y \gets x \f$ .
 *
 *    \param N Number of entries in each vector.
 *    \param X Vector to be duplicated.
 *    \param INCX Increment between subsequent entries in \p X.
 *    \param Y Vector to be overwritten.
 *    \param INCY Increment between subsequent entries in \p Y.
 */
extern void dcopy_(int *N, double *X, int *INCX, double *Y, int *INCY);

/*! \brief Scale vector \f$ x \f$ by a factor \f$ \alpha \f$.
 *
 *  Update \f$ x \gets \alpha x \f$ .
 *
 *    \param N Number of entries in \p X .
 *    \param ALPHA Scalar multiplier \f$ \alpha \f$ .
 *    \param X Vector being scaled. 
 *    \param INCX Increment between subsequent entries in \p X .
 */
extern void dscal_(int *N, double *ALPHA, double *X, int* INCX) ;

/*! \brief Scale a vector \f$ x \f$ by \f$ \alpha \f$ and add to \f$ y \f$. 
 *
 *  Update \f$ y \gets \alpha x + y \f$ .
 *
 *    \param N Number of entries in each vector.
 *    \param ALPHA Scalar multiplier \f$ \alpha \f$ .
 *    \param X Vector to be sclaled.
 *    \param INCX Increment between subsequent entries in \p X.
 *    \param Y Vector to be updated. 
 *    \param INCY Increment between subsequent entries in \p Y.
 */
extern void daxpy_(int *N, double *ALPHA,
                   double *X, int *INCX,
                   double *Y, int *INCY);

/*! \brief Compute the scalar product \f$ x^T y \f$.
 *
 *    \param N Number of entries in each vector.
 *    \param X First vector.
 *    \param INCX Increment between subsequent entries in \p X.
 *    \param Y Second vector.
 *    \param INCY Increment between subsequent entries in \p Y.
 *    \return Scalar product of \p X and \p Y .
 */
extern double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);

/*! \brief Compute the Euclidean norm \f$ \| x \|_2 f$.
 *
 *    \param N Number of entries in \p X .
 *    \param X Vector. 
 *    \return The norm of \p X .
 */
extern double dnrm2_(int *N, double *X, int *INCX);

/*! \brief Scale vector by some factor.
 */
extern void dscal_(int *N, double *ALPHA, double *X, int *INCX);

#endif

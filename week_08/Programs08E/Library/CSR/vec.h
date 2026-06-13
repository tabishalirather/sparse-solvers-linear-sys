/*! \file
 *  \brief Vector allocation and initialisation.
 *
 *  \author Michael Wende 
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef VEC_H
#define VEC_H

/*! \brief Produce a vector with entries uniformly distributed in [0, 1).
 *
 *    \param n Length.
 *    \return Newly allocated.
 */
double *vec_rand(int n);

/*! \brief Allocate an uninitialized vector.
 *
 *    \param n Length.
 *    \return Allocated memory.
 */
double *vec_raw(int n);

/*! \brief Produce a vector with unit entries.
 *
 *    \param n Number of entries.
 *    \return Newly allocated.
 */
double *vec_ones(int n);

/*! \brief Produce a vector with zero entries.
 *
 *    \param n Number of entries.
 *    \return Newly allocated.
 */
double *vec_zeros(int n);

/*! \brief Copy one vector into another.
 *
 *    \param n Common length.
 *    \param x Source vector.
 *    \param y Destination vector.
 */
void vec_copy(int n, const double *x, double *y);

/*! \brief Print a vector to standard output.
 *
 *    \param n Length.
 *    \param x Vector.
 */
void vec_print(int n, const double *x);

/*! \brief Euclidean norm of a vector.
 *
 *    \param n Length.
 *    \param x Vector.
 */
double vec_norm(int n, const double *x);

/*! \brief Infinity norm of a vector.
 *
 *    \param n Length.
 *    \param x Vector.
 */
double vec_norm_inf(int n, const double *x);

/*! \brief Dot product between vectors.
 *
 *    \param n Length.
 *    \param x Vector.
 *    \param y Vector.
 */
double vec_dot(int n, const double *x, const double *y);

/*! \brief Scale a vector by some factor.
 *
 *    \param n Length.
 *    \param alpha Factor by which to scale.
 *    \param x Vector to be scaled.
 */
void vec_scale(int n, double alpha, double *x);

/*! \brief Scaled vector update \f$ y \gets \alpha x + y \f$.
 *
 *    \param n Length of vectors x and y.
 *    \param alpha Factor by which x is scaled.
 *    \param x Input vector to be scaled.
 *    \param y Vector to be updated.
 */
void vec_axpy(int n, double alpha, const double *x, double *y);

#endif

/*! \file
 *  \brief BLAS-like Givens rotations.
 *  \author Jens-Peter Zemke.
 *  \author Sabine Le Borne.
 */

/*! \brief Generate Givens rotation.
 *
 *  Generate a Givens rotation which will can be used to eliminate the
 *  second entry in a two-by-one vector z = (x, y)^T.
 *  The corresponding sine and cosine are returned and the vector z is
 *  modified in place.
 *
 *  \param x Will be overwriten with the norm of z.
 *  \param y Will also be overwritten, but not necessarily with zero.
 *  \param c Will be set to the cosine value of the Givens rotation.
 *  \param s Will be set to the sine value of the Givens rotation.
 */
void givens_rotg(double *a, double *b, double *c, double *s);

/*! \brief Apply a single Givens rotation multiple times.
 *
 *  The Givens rotation defined by the given cosine and sine values is applied
 *  to a sequence of two-by-one column vector (x_i, y_i)^T for i = 1, ..., n.
 *  The values x_i and y_i are elements of given vectors x and y, not necessarily
 *  adjacent but with a given increment between them.
 *
 *  \param n Number of applications.
 *  \param x First element in vector x.
 *  \param incx Increment between positions in x.
 *  \param y Left-most element in the lower row.
 *  \param incy First element in vector y.
 *  \param c Cosine value of the Givens rotation.
 *  \param s Sine value of the Givens rotation.
 */
void givens_rot(int n, double *x, int incx, double *y, int incy, double c, double s);

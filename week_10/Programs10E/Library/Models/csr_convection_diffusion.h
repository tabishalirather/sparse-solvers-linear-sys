/*! \file
 *  \brief Convection diffusion test problem in the unit cube in arbitrary dimension.
 *
 *  The functions in this file set up the matrix \f$ A \f$ from the
 *  discretisation of the convection diffusion equation.
 *
 *  The differential equation of the convection diffusion problem is
 *    \f[ - \epsilon \Delta u + b^T \nabla u = f \f]
 *  where \f$ u \f$ is the function to be determined and \f$ f \f$ is
 *  a given right hand side.
 *  The convection part is described by the vector \f$ b \in \mathbb{R}^{d} \f$
 *  which can be constant or dependent on the location in the unit square.
 *  The diffusion part is described by a constant \f$ \epsilon > 0 \f$.
 *
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_CONVECTION_DIFFUSION_H
#define CSR_CONVECTION_DIFFUSION_H

/*! \brief Convection diffusion equation where all entries of the vector \f$ b \f$
 *         are equal to a constant (no spatial dependency).
 *
 *    \param n Number of intervals into which the the unit cube is subdivided
 *             along each dimension such that the distance between neighboring
 *             nodes will be h = 1 / n.
 *    \param d Spatial dimensions.
 *    \param epsilon
 *    \param b_const A single value defining all of the \f$ d \f$ entries in \f$ b \f$.
 */
csr *csr_const_convection_diffusion(int n, int d, double epsilon, double b_const);

/*! \brief Convection diffusion equation where the vector \f$ b \f$
 *         may be a function depending on the location in the unit square.
 *    \param n
 *    \param d
 *    \param epsilon
 *    \param b_func Pointer to a function that computes the vector \f$ b \in \mathbb{R}^{d} \f$
 *                  given the parameters n, d and the index of the grid point where \f$ b \f$
 *                  has to be evaluated.
 *                  The last argument to the function is the output.
 */
csr *csr_convection_diffusion(int n, int d, double epsilon, void (*b_func)(int, int, int, double *));

#endif

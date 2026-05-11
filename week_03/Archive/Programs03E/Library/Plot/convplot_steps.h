/*! \file
 *  \brief Plotting routine for convergence diagrams.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CONVPLOT_STEPS_H
#define CONVPLOT_STEPS_H

/*! \brief Plot convergence curves of multiple methods into one diagram.
 *
 *    \param filename Name of the file for an eps picture.
 *    \title title Title string to show up on top of the plot.
 *    \param k Number of methods.
 *    \param m Number of steps for each method.
 *    \param y For each method one initial value and one value after each step.
 *    \param name Legend entry for each method.
 */
void convplot_steps(const char *filename,
		    const char *title, const char *xlabel, const char *ylabel,
		    const int k, const int *m, double **y, char **name);

#endif

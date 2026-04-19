/*! \file
 *  \brief Draw a picture of the non-zero pattern of a sparse matrix.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_SPY_H
#define CSR_SPY_H

/*! \brief Produce a postscript file with a spyplot of a
 *         sparse matrix showing the locations and magnitudes
 *         of its non-zero entries.
 *
 *    \param A CSR matrix.
 *    \param filename A file with this name will be generated in the
 *                    current directory.
 *                    An existing file will be overwritten and a filename
 *                    suffix is not appended.
 *    \param style A single character 's' or 'c' determining whether
 *                 non-zero entries are to be shown as squares or circles.
 */
void csr_spy(csr *A, const char *filename, char style);

#endif



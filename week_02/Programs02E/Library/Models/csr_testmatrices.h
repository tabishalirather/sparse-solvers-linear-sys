/*! \file
 *  \brief Routine to construct test matrices which can be used to
 *         compare the effect of different reordering strategies.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
*/

#ifndef CSR_TESTMATRICES_H
#define CSR_TESTMATRICES_H

/*! \brief Beispielmatrizen zum Testen der Umordnungsalgorithmen.
 *
 *    \param number Auswahl einer Beispielmatrix 0, 1, 2 oder 3.
 *
 *    \return Beispielmatrix.
 */
csr *csr_testmatrices(int number);

#endif

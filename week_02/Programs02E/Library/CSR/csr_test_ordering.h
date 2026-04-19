/*! \file
 *  \brief Apply multiple reorderings to a given matrix and produce
 *         plots of the resulting non-zero patterns.
 *
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef CSR_TEST_ORDERING_H
#define CSR_TEST_ORDERING_H

/*! \brief Test der Umordungsalgorithmen BFS, CMK und RCM fuer eine
 *         vorgegebene Matrix und Startknoten.
 *
 *  Die Besetztheitsstruktur der umgeordneten Matrizen wird in Bilddateien geschrieben.
 *
 *    \param A Quadratische Matrix.
 *    \param name Grundname fuer die erzeugten Bilddateien.
 *    \param start Index des Startknotens. Bei Angabe von -1 wird
 *                 ein Knoten minmalen Grades verwendet.
 */
void csr_test_ordering(csr *A, char *name, int start);

#endif

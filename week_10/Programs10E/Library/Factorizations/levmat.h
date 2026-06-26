/*! \file
 *  \brief Linked list representation of a sparse matrix intended for ILU(p).
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#ifndef LEVMAT_H
#define LEVMAT_H

struct levmat_s;
typedef struct levmat_s levmat_t;

struct levent_s;
typedef struct levent_s levent_t;

/*! \brief Sparse matrix where each row is represented as a linked list.
 *         For each entry the level of fill is stored.
 *
 *  The matrix will be square when this structure is used to implement the
 *  ILU(p) factorisation as intendend.
 *  Here we allow separate number of rows and columns anyway.
 */
struct levmat_s{
  int nr;             /*!< Number of rows.                     */
  int nc;             /*!< Number of columns.                  */
  levent_t **entries; /*!< Array of linked lists for each row. */
};

/*! \brief Entry of a level matrix.
 */
struct levent_s {
  int lev;        /*!< Level of fill.                  */
  int col;        /*!< Column index.                   */
  double val;     /*!< Coefficient value.              */
  levent_t *next; /*!< Next entry within this row.     */
  levent_t *prev; /*!< Previous entry within this row. */
};

/*! \brief Insert level entry into list.
 *
 *    \param list Level entry to which the new node will be prepended.
 *    \param node New level entry to insert.
 */
void levent_insert_before(levent_t *list, levent_t *node);

/*! \brief Print level matrix to standard output.
 *  \param lm Matrix to be printed.
 */
void levmat_print(levmat_t *lm);

/*! \brief Convert a CSR matrix into a level matrix.
 *
 *  \param A CSR matrix.
 *  \return Newly allocated level matrix.
 */
levmat_t *levmat_from_csr(csr *A);

/*! \brief Convert a level matrix to CSR format.
 *
 *    \param lm Level matrix to be converted.
 *    \return Newly allocated CSR matrix.
 */
csr *csr_from_levmat(levmat_t *lm);

/*! \brief Delete a level matrix.
 *
 *  \param lm Level matrix to delete.
 */
void levmat_delete(levmat_t *lm);

#endif

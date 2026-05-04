/*! \file
 *  \author Sabine Le Borne
 */

#ifndef PERM_H
#define PERM_H

struct perm_s;
typedef struct perm_s perm_t;
struct perm_s {
  int size;
  int *idx;
  int *idx_inv;
};

perm_t *perm_new(int n, int *idx, int *idx_inv);

/*! \brief Delete permutation and free memory.
 *
 * \param perm Permutation to be deleted.
 */
void perm_delete(perm_t *perm);

/*! \brief Print permutation to standard output.
 *
 *  \param perm Permutation.
 */
/*! \brief Identity permutation.
 *
 *  \param n Size of the permutation.
 *  \return Identity permutation.
 */
perm_t *perm_identity(int n);

/*! \brief Create a new permutation intialized with given index vectors.
 *
 *  \param n Size of the permutation.
 *  \param idx Index vector describing the permutation.
 *  \param idx_inv Index vector describing the inverse permutation.
 *  \return New allocation object.
 */

void perm_print(perm_t *perm);

/*! \brief Check that a permutations internal index vectors represent a
 *         permutation and its inverse as they should.
 *
 *  \param perm Permutation to check.
 *  \return True if index vectors are set up correctly, otherwise false.
 */
bool perm_check(perm_t *perm);

/*! \brief Compute the inverse of a given permutation vector.
 *
 *  \param n Number of indices.
 *  \param p Vector of n indices representing a permutation.
 *  \return Vector of n indices representing the inverse permutation.
 */
int *perminv(int n, const int *p);

/*! \brief Apply a row permutation given as an inverse permutation to a
 *         (column) vector.
 *
 *    \param n Length of all vectors.
 *    \param b Vector to be permuted.
 *    \param pinv Vector of indices representing the inverse permutation
 *                corresponding the row permutation to be applied.
 *
 *    \return Permuted vector Pb = (b[pinv[1]], ..., b[pinv[n]]).
 *    \brief (Using p instead of pinv would return Pb with Pb[p[i]]=b[i].)
 */
double *vec_permute(const int n, const double *b, const int *pinv);

#endif

/*! \file
 *  \brief GMRes for solving a linear system with CSR matrix.
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

/*! \brief Structure describing an instance of GMRes
 */
struct gmres_s {
  int n;                   /*!< \brief System size. */
  csr *A;                  /*!< \brief System matrix. */
  double *b;               /*!< \brief Right hand side vector. */
  double *x0;              /*!< \brief Initial guess. */
  double *x;               /*!< \brief Current approximation. */
  int iter;                /*!< \brief Current inner iteration number. */
  int iter_max;            /*!< \brief Maximum number of inner iterations. */
  int restart;             /*!< \brief Restart parameter. */
  double tol;              /*!< \brief Tolerance. */
  double *residual_norms;  /*!< \brief Vector of residual norms. */
  int n_ortho;             /*!< \brief Number of orthogonalisations. */
};

/*! \brief GMRes type.
 */
typedef struct gmres_s gmres_t;

/*! \brief Create a new GMRes instance.
 *
 *  \param A System matrix.
 *  \param b Right hand side vector.
 *  \param x0 Start vector/ initial guess.
 */
gmres_t *gmres_new(csr *A, double *b, double *x0);

/*! \brief Set the maximum number of inner iterations to perform.
 *
 *  \param gmres Instance of GMRes.
 *  \param iter_max Maximum number of inner iterations combined across
 *                  all restarts.
 */
void gmres_iter_max(gmres_t *gmres, int iter_max);

/*! \brief Set the restart parameter.
 *
 *  \param gmres Instance of GMRes.
 *  \param restart Number of inner iterations before a restart is performed.
 *                 A zero value means not to restart (full GMRes).
 */
void gmres_restart(gmres_t *gmres, int restart);

/*! \brief Set the relative tolerance to be used as stopping criterion.
 *
 *  \param gmres Instance of GMRes.
 *  \param tol Tolerance interpreted as the ratio of the current residual norm
 *             to the initial residual norm.
 */
void gmres_tol(gmres_t *gmres, double tol);

/*! \brief Run GMRes until either the maximum number of inner iterations or
 *         the desired accuracy has been reached.
 *
 *  \param gmres Instance of GMRes.
 */
void gmres_solve(gmres_t *gmres);

/*! \brief Fetch the current solution after GMRes has been run.
 *
 *  \param gmres Instance of GMRes.
 *  \return Pointer to the current solution internal to the GMRes instance.
 *          Will become invalid when the GMRes instance is destroyed.
 */
double *gmres_solution(gmres_t *gmres);

/*! \brief Query the total number of inner iterations that have been
 *         performed.
 *
 *  \param gmres Instance of GMRes.
 *  \return Total number of inner iterations, at most iter_max.
 */
int gmres_iter(gmres_t *gmres);

/* \brief Fetch vector holding the norm of the initial residual and of the
 *        residual after each step.
 *
 *  \param gmres Instance of GMRes.
 *  \return Vector of length (iter + 1) where iter is the number of inner
 *          iterations that where performed.
 *          The returned pointer will be freed when the GMRes instance is
 *          destroyed.
 */
double *gmres_residual_norms(gmres_t *gmres);

/*! \brief Free all memory associated with a GMRes instance.
 *
 *  \param gmres Instance of GMRes.
 */
void gmres_delete(gmres_t *gmres);

/*! \file
 *  \brief Library header including all individual header files.
 */
#ifndef ITSOLVER_H
#define ITSOLVER_H

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "blas.h"
#include "CSR/csr.h"
#include "CSR/csr_ascii.h"
#include "CSR/csr_mv.h"
#include "CSR/csr_norm.h"
#include "CSR/csr_permute.h"
#include "CSR/csr_spy.h"
#include "CSR/csr_mindeg_startnode.h"
#include "CSR/csr_test_ordering.h"
#include "CSR/csr_transpose.h"
#include "CSR/csr_vm.h"
#include "CSR/vec.h"
#include "LinearSolver/csr_gauss_seidel.h"
#include "LinearSolver/csr_gauss_seidel_reverse.h"
#include "LinearSolver/csr_jacobi.h"
#include "LinearSolver/csr_richardson.h"
#include "LinearSolver/csr_sgs.h"
#include "LinearSolver/csr_sor.h"
#include "LinearSolver/csr_ssor.h"
#include "LinearSolver/prec.h"
#include "LinearSolver/stopcrit.h"
#include "Models/csr_laplace.h"
#include "Models/csr_testmatrices.h"
#include "Permutations/sort_degree.h"
#include "Permutations/perm_level_set.h"
#include "Permutations/perm.h"
#include "Plot/conv_plot.h"
#include "Plot/convplot_steps.h"

#endif

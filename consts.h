/* "consts.h" contains several constants used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
#define PI 3.1415926535897932
// base of natural logarithm
#define EBN 2.7182818284590451
// default numerical epsilon
#define EPS 1.0e-6
// maximum nunber of function evaluations in the class of adaptive \
Gauss-Kronrod quadrature
#define MAX_FUN_EVAL 8192
// tolerance of the adaptive Gauss-Kronrod quadrature
#define ABS_QUAD_TOL 1.0e-9
// IR scale of QCD in units of GeV
#define LAMBDA_QCD 0.234
// Select if the UV term in the Maris-Tandy model is used.
#define IS_UV_MT true
// Select if type-I Gauss-Chebyshev quadrature rule is applied for the grid \
of the Euclidea-space solution.
#define IS_SF_TYPE_I true
// relative and absolute tolerances of the iterative solver
#define TOL_TEST_REL_TOL 1.0e-6
#define TOL_TEST_ABS_TOL 1.0e-6
// name for binary file of the solution from the iterative solver of the SDE
#define SF_LIGHT_FILENAME "sde_solution_light_cpl.bin"
// name for binary file of paramters in the SDE
#define SDE_PMTS_FILENAME "sde_pmts.bin"

/*
"main_sde_iter.cpp" is the source for the main function that runs the
iterative solver of the SDE in the Euclidean space. The solution is to
be saved in a binary file named by the sf_filename member member of
SDE_PMTS_IO.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark.

SDE_quark is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SDE_quark is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SDE_quark. If not, see <https://www.gnu.org/licenses/>.
*/
// C++ header files
#include<cstdlib>// size_t
#include<stdexcept>// runtime_error()
#include<string>// string
#include<cmath>// abs(), pow(), exp(), log()
#include<array>// array<>
#include<vector>// vector<>
#include<iostream>// cout
#include<fstream>// ifstream, ofstream
#include"omp.h"// OpenMP
// definitions of constants and classes
#include"consts.h"
#include"source_gker/g_mt.h"// class Gfun
#include"source_utility/gk_rw.h"// class GK
#include"source_utility/quad_gk.h"// class QUAD_GK<typename T>: private GK
#include"source_sde_quark/sde_ker.h"// class SDE_ker<typename T>: public Gfun
#include"source_sde_quark/sde_ker_z_itg.h"// class Zker<typename T>: \
public SDE_ker<T>, public QUAD_GK<T>
#include"source_utility/gc_rw.h"// class GC
#include"source_utility/gc_grid.h"// class GC_GRID: private GC
#include"source_sde_quark/rvt.h"// class RVT
#include"source_sde_quark/sf.h"// class SF: public RVT, public GC_GRID
#include"source_sde_quark/sde_quark_iter.h"// class SDE_quark_iter ( Zker<double> ): \
public SF
#include"source_utility/tol_test.h"// class TOL_TEST
#include"source_sde_quark/sde_solution.h"// struct sde_solution; \
struct sde_solution_array<unsigned num_grid> class SDE_Solution_IO ( sde_solution )
#include"source_sde_quark/sde_quark.h"// class SDE_quark ( sde_solution ): \
public SDE_quark_iter, private TOL_TEST
#include"sde_pmts_io.h"// class SDE_PMTS_IO
// The main function for the executable of the iterative solver. Modify the \
variable num_processes to the fit the hardware that runs the executable.
using namespace std;
int main( int argc, char* argv[] )
{
  // Load model paramters from the binary file.
  string pmts_filename = SDE_PMTS_FILENAME;// name of the binary file for parameters
  SDE_PMTS_IO sde_pmts( pmts_filename );// Create an instance of the parameter class.
  sde_pmts.print_pmts();// Display parameters.
  // other physical parameters
  double m_rnm_light = 3.6964e-3;// renormalized mass for the quark in GeV
  unsigned quark_flavor = 0;// 0 for light quarks
  double mu2 = pow( 19.0, 2 );// renormalization scale in (GeV)^2
  double uv_cutoff = 1.0e+6;// UV cutoff \Lambda^2_{UV} in (GeV)^2
  // running paramters for the iterative solver of the SDE
  bool is_pv_mt = false;// Do not apply the Pauli-Villars regularization.
  bool is_linear_ren_cond = true;// Applies the regular renormalization condition.
  unsigned max_n_iter = 75;// maximum number of iterations in the iterative solver
  unsigned num_processes = 36;// number of processes in OpenMP
  // iterative solver
  SDE_quark sde_quark( sde_pmts.d_scale, sde_pmts.omega, is_pv_mt,
    sde_pmts.lambda_pv, sde_pmts.c_scale_sde, sde_pmts.n_u_grid, m_rnm_light, mu2,
    uv_cutoff, TOL_TEST_REL_TOL, TOL_TEST_ABS_TOL, max_n_iter, quark_flavor,
    is_linear_ren_cond );
  sde_quark.iter_solve_sde( num_processes );// Run the iterative solver.
  sde_quark.save_solution( sde_pmts.sf_filename );// Save the solution.
  return 0;
}

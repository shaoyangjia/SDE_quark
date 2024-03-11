/*
"bind_sde.cpp" contains the instructions for pybind11 to build an
interface to access the iterative solver of the quark propagator SDE
in Python.

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
#include<pybind11/pybind11.h>// header file of pybind11
#include<pybind11/stl.h>// support of standard libraries in pybind11
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
#include"consts.h"// constants
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
using namespace std;
namespace py = pybind11;
// definitions of pybind11 module named "sde_iter_solver"
PYBIND11_MODULE( sde_iter_solver, m )
{
  m.doc() = "Module compiled in C++ for the class of iterative solver of the "
    "Schwinger-Dyson equation (SDE) for quark propagators in the Euclidean space.";
  // class of the iterative solver of the SDE
  py::class_<SDE_quark>( m, "SDE_quark" )
  .def(py::init<double, double, bool, double, double, unsigned,
    double, double, double, double, double, unsigned, unsigned, bool>())
  .def("iter_solve_sde", &SDE_quark::iter_solve_sde)
  .def("print_solutions", &SDE_quark::print_solutions)
  .def("get_solution_struct", &SDE_quark::get_solution_struct)
  .def("save_solution", &SDE_quark::save_solution)
  .def("load_solution", &SDE_quark::load_solution)
  .def("n_grid", &SDE_quark::get_n_grid)
  .def("c_scale", &SDE_quark::get_c_scale)
  .def("u_grid", &SDE_quark::get_u_grid)
  .def("a_grid", &SDE_quark::get_a_grid)
  .def("b_grid", &SDE_quark::get_b_grid)
  .def("z_2", &SDE_quark::get_z_2)
  .def("z_m", &SDE_quark::get_z_m)
  .def("uv_cutoff", &SDE_quark::get_uv_cutoff)
  .def("get_d_scale", &SDE_quark::SDE_quark_iter::get_d_scale)
  .def("get_omega", &SDE_quark::SDE_quark_iter::get_omega)
  .def("get_is_pv_mt", &SDE_quark::SDE_quark_iter::get_is_pv_mt)
  .def("get_lambda_pv", &SDE_quark::SDE_quark_iter::get_lambda_pv)
  .def("get_c_scale", &SDE_quark::get_c_scale)
  .def("get_m_rnm", &SDE_quark::get_m_rnm)
  .def("get_mu2", &SDE_quark::get_mu2)
  .def("get_iter_rel_tol", &SDE_quark::get_iter_rel_tol)
  .def("get_iter_abs_tol", &SDE_quark::get_iter_abs_tol)
  .def("get_max_n_iter", &SDE_quark::get_max_n_iter)
  .def("get_quark_flavor", &SDE_quark::get_quark_flavor)
  .def("get_is_linear_ren_cond", &SDE_quark::SDE_quark_iter::get_is_linear_ren_cond)
  // pickle suport for python multiprocessing
  .def(py::pickle
    (
      []( SDE_quark &p )
      { return py::make_tuple( p.get_d_scale(), p.get_omega(), p.get_is_pv_mt(),
          p.get_lambda_pv(), p.get_c_scale(), p.get_n_grid(), p.get_m_rnm(),
          p.get_mu2(), p.get_uv_cutoff(), p.get_iter_rel_tol(), p.get_iter_abs_tol(),
          p.get_max_n_iter(), p.get_quark_flavor(), p.get_is_linear_ren_cond() );
      },
      []( py::tuple t )
      {
        if( t.size() != 14 )
        { throw runtime_error("Invalid pickle state of SDE_quark"); }
        SDE_quark p( t[0].cast<double>(), t[1].cast<double>(), t[2].cast<bool>(),
          t[3].cast<double>(), t[4].cast<double>(), t[5].cast<unsigned>(),
          t[6].cast<double>(), t[7].cast<double>(), t[8].cast<double>(),
          t[9].cast<double>(), t[10].cast<double>(), t[11].cast<unsigned>(),
          t[12].cast<unsigned>(), t[13].cast<bool>() );
        return p;
      }
    ));
}

/*
"bind_sde_ker.cpp" contains the instructions for pybind11 to build an
interface for the access of the kernel functions in the quark self-energy
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
#include<pybind11/complex.h>// complex-valued variables in pybind11
#include<cstdlib>// size_t
#include<stdexcept>// runtime_error()
#include<string>// string
#include<complex>// complex-valued variables
#include<cmath>// abs(), pow(), exp(), log()
#include<vector>// vector<>
#include<iostream>// cout
#include<fstream>// ifstream, ofstream
// definitions of constants and classes
#include"consts.h"// constants
#include"source_gker/g_mt.h"// class Gfun
#include"source_sde_quark/sde_ker.h"// class SDE_ker<typename T>: public Gfun
#include"source_sde_quark/sde_solution.h"// struct sde_solution; \
struct sde_solution_array<unsigned num_grid> class SDE_Solution_IO ( sde_solution )
#include"source_sde_quark/quark_prop.h"// class Grid ( NCS, sde_solution, \
  sde_solution_array<> ) ( #include"natural_cubic_spline.h" )
#include"source_sde_quark/sde_quark_ker.h"// class SDE_SE_ker: \
public SDE_ker<T>, public Grid
using namespace std;
namespace py = pybind11;
// definitions of pybind11 module named "sde_ker"
PYBIND11_MODULE( sde_ker, m )
{
  m.doc() = "Module compiled in C++ for the kernel functions of the integrals "
    "in the Schwinger-Dyson equation (SDE) of the quark propagator.";
  // kernel functions of the quark self-energy with real-valued momentum
  py::class_<SDE_SE_ker<double>>( m, "SDE_SE_rker" )
  .def(py::init<double, double, unsigned, double, vector<double>,
    vector<double>, vector<double>, double, bool, double>())
  .def("self_energy_ker", &SDE_SE_ker<double>::self_energy_ker)
  .def("get_omega", &SDE_SE_ker<double>::ker_omega)
  .def("get_d_scale", &SDE_SE_ker<double>::ker_d_scale)
  .def("get_n_grid", &SDE_SE_ker<double>::ker_n_grid)
  .def("get_c_scale", &SDE_SE_ker<double>::ker_c_scale)
  .def("get_u_grid", &SDE_SE_ker<double>::ker_u_grid)
  .def("get_a_grid", &SDE_SE_ker<double>::ker_a_grid)
  .def("get_b_grid", &SDE_SE_ker<double>::ker_b_grid)
  .def("get_uv_cutoff", &SDE_SE_ker<double>::ker_uv_cutoff)
  .def("get_is_pv_mt", &SDE_SE_ker<double>::ker_is_pv_mt)
  .def("get_lambda_pv", &SDE_SE_ker<double>::ker_lambda_pv)
  // pickle suport for python multiprocessing
  .def( py::pickle(
    []( SDE_SE_ker<double> &p )
    {
      return py::make_tuple( p.ker_omega(), p.ker_d_scale(), p.ker_n_grid(),
      p.ker_c_scale(), p.ker_u_grid(), p.ker_a_grid(), p.ker_b_grid(),
      p.ker_uv_cutoff(), p.ker_is_pv_mt(), p.ker_lambda_pv() );
    },
    []( py::tuple t )
    {
      if( t.size() != 10 )
      { throw runtime_error("Invalid pickle state of SDE_SE_ker<double>"); }
      SDE_SE_ker<double> p( t[0].cast<double>(), t[1].cast<double>(),
        t[2].cast<unsigned>(), t[3].cast<double>(), t[4].cast<vector<double>>(),
        t[5].cast<vector<double>>(), t[6].cast<vector<double>>(),
        t[7].cast<double>(), t[8].cast<bool>(), t[9].cast<double>() );
      return p;
    }
  ));
  // kernel functions of the quark self-energy with complex-valued momentum
  py::class_<SDE_SE_ker<complex<double>>>( m, "SDE_SE_cker" )
  .def(py::init<double, double, unsigned, double, vector<double>,
    vector<double>, vector<double>, double, bool, double>())
  .def("self_energy_ker", &SDE_SE_ker<complex<double>>::self_energy_ker)
  .def("get_omega", &SDE_SE_ker<complex<double>>::ker_omega)
  .def("get_d_scale", &SDE_SE_ker<complex<double>>::ker_d_scale)
  .def("get_n_grid", &SDE_SE_ker<complex<double>>::ker_n_grid)
  .def("get_c_scale", &SDE_SE_ker<complex<double>>::ker_c_scale)
  .def("get_u_grid", &SDE_SE_ker<complex<double>>::ker_u_grid)
  .def("get_a_grid", &SDE_SE_ker<complex<double>>::ker_a_grid)
  .def("get_b_grid", &SDE_SE_ker<complex<double>>::ker_b_grid)
  .def("get_uv_cutoff", &SDE_SE_ker<complex<double>>::ker_uv_cutoff)
  .def("get_is_pv_mt", &SDE_SE_ker<complex<double>>::ker_is_pv_mt)
  .def("get_lambda_pv", &SDE_SE_ker<complex<double>>::ker_lambda_pv)
  // pickle suport for python multiprocessing
  .def( py::pickle(
    []( SDE_SE_ker<complex<double>> &p )
    {
      return py::make_tuple( p.ker_omega(), p.ker_d_scale(), p.ker_n_grid(),
      p.ker_c_scale(), p.ker_u_grid(), p.ker_a_grid(), p.ker_b_grid(),
      p.ker_uv_cutoff(), p.ker_is_pv_mt(), p.ker_lambda_pv() );
    },
    []( py::tuple t )
    {
      if( t.size() != 10 )
      { throw runtime_error("Invalid pickle state of SDE_SE_ker<complex<double>>"); }
      SDE_SE_ker<complex<double>> p( t[0].cast<double>(), t[1].cast<double>(),
        t[2].cast<unsigned>(), t[3].cast<double>(), t[4].cast<vector<double>>(),
        t[5].cast<vector<double>>(), t[6].cast<vector<double>>(),
        t[7].cast<double>(), t[8].cast<bool>(), t[9].cast<double>() );
      return p;
    }
  ));
}

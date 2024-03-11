/*
"bind_ker.cpp" contains the instructions for pybind11 to build an interface
for the access of the fixed-grid routines to compute the quark propagator in
Python.

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
// C++ header files
#include<cstdlib>// size_t
#include<stdexcept>// runtime_error()
#include<string>// string
#include<complex>// complex-valued variables
#include<cmath>// abs(), pow(), exp(), log()
#include<array>// array<>
#include<vector>// vector<>
#include<iostream>// cout
#include<fstream>// ifstream, ofstream
// definitions of constants and classes
#include"consts.h"// constants
#include"filenames_complex_prop.h"// classes for the propagator with \
complex-valued momentum
using namespace std;
namespace py = pybind11;
// definitions of pybind11 module named "ker_itg"
PYBIND11_MODULE( ker_itg, m )
{
  m.doc() = "Module compiled in C++ for the computation of quark propagators "
    "with complex-valued momentum applying fixed-grid quadrature rules.";
  // Class of the quark propagator with real-valued momentum.
  py::class_<Complex_prop<double>>( m, "Real_prop" )
    .def(py::init<double, double,
      unsigned, double, vector<double>, vector<double>, vector<double>,
      double, double, double, bool, unsigned, bool, double>())
    .def(py::init<double, double, unsigned, sde_solution>())
    .def(py::init<double, double, unsigned, SDE_Solution_IO>())
    .def("inv_prop_complex", &Complex_prop<double>::inv_prop_complex)
    .def("inv_prop_b", &Complex_prop<double>::inv_prop_b)
    .def("inv_prop_ab", &Complex_prop<double>::inv_prop_ab)
    //.def("prop_fun", &Complex_prop<double>::prop_fun)
    .def("prop_fun", static_cast<double ( Complex_prop<double>::* )( double,
      bool )>( &Complex_prop<double>::prop_fun ))
    .def("prop_fun_vec", static_cast<array<double,2>
      ( Complex_prop<double>::* )( double )>( &Complex_prop<double>::prop_fun ))
    .def("self_energy", &Complex_prop<double>::self_energy)
    .def("reset_grid", &Complex_prop<double>::reset_grid)
    .def("reset_rnm_consts", &Complex_prop<double>::reset_rnm_consts)
    .def("get_omega", &Complex_prop<double>::get_omega)
    .def("get_d_scale", &Complex_prop<double>::get_d_scale)
    .def("get_n_grid", &Complex_prop<double>::get_n_grid)
    .def("get_c_scale", &Complex_prop<double>::get_c_scale)
    .def("get_u_grid", &Complex_prop<double>::get_u_grid)
    .def("get_a_grid", &Complex_prop<double>::get_a_grid)
    .def("get_b_grid", &Complex_prop<double>::get_b_grid)
    .def("get_uv_cutoff", &Complex_prop<double>::get_uv_cutoff)
    .def("get_z_2", &Complex_prop<double>::get_z_2)
    .def("get_m_bare", &Complex_prop<double>::get_m_bare)
    .def("get_is_linear_ren_cond", &Complex_prop<double>::get_is_linear_ren_cond)
    .def("get_n_grid_gc", &Complex_prop<double>::get_n_grid_gc)
    .def("get_is_pv_mt", &Complex_prop<double>::get_is_pv_mt)
    .def("get_lambda_pv", &Complex_prop<double>::get_lambda_pv)
    // pickle suport for python multiprocessing
    .def(py::pickle(
      []( Complex_prop<double> &p )
      {
        return py::make_tuple( p.get_omega(), p.get_d_scale(), p.get_n_grid(),
        p.get_c_scale(), p.get_u_grid(), p.get_a_grid(), p.get_b_grid(),
        p.get_uv_cutoff(), p.get_z_2(), p.get_m_bare(), p.get_is_linear_ren_cond(),
        p.get_n_grid_gc(), p.get_is_pv_mt(), p.get_lambda_pv() );
      },
      []( py::tuple t )
      {
        if( t.size() != 14 )
        { throw runtime_error("Invalid pickle state of Complex_prop<double>"); }
        Complex_prop<double> p( t[0].cast<double>(), t[1].cast<double>(),
          t[2].cast<unsigned>(), t[3].cast<double>(), t[4].cast<vector<double>>(),
          t[5].cast<vector<double>>(), t[6].cast<vector<double>>(),
          t[7].cast<double>(), t[8].cast<double>(), t[9].cast<double>(),
          t[10].cast<bool>(), t[11].cast<unsigned>(), t[12].cast<bool>(),
          t[13].cast<double>() );
        return p;
      }
    ));
  // Class of the quark propagator with complex-valued momentum.
  py::class_<Complex_prop<complex<double>>>( m, "Complex_prop" )
    .def(py::init<double, double,
      unsigned, double, vector<double>, vector<double>, vector<double>,
      double, double, double, bool, unsigned, bool, double>())
    .def(py::init<double, double, unsigned, sde_solution>())
    .def("inv_prop_complex", &Complex_prop<complex<double>>::inv_prop_complex)
    .def("inv_prop_b", &Complex_prop<complex<double>>::inv_prop_b)
    .def("inv_prop_ab", &Complex_prop<complex<double>>::inv_prop_ab)
    //.def("prop_fun", &Complex_prop<complex<double>>::prop_fun)
    .def("prop_fun", static_cast<complex<double>
      ( Complex_prop<complex<double>>::* )( complex<double>, bool )>(
        &Complex_prop<complex<double>>::prop_fun ))
    .def("prop_fun_vec", static_cast<array<complex<double>,2>
      ( Complex_prop<complex<double>>::* )( complex<double> )>(
        &Complex_prop<complex<double>>::prop_fun ))
    .def("self_energy", &Complex_prop<complex<double>>::self_energy)
    .def("reset_grid", &Complex_prop<complex<double>>::reset_grid)
    .def("reset_rnm_consts", &Complex_prop<complex<double>>::reset_rnm_consts)
    .def("get_omega", &Complex_prop<complex<double>>::get_omega)
    .def("get_d_scale", &Complex_prop<complex<double>>::get_d_scale)
    .def("get_n_grid", &Complex_prop<complex<double>>::get_n_grid)
    .def("get_c_scale", &Complex_prop<complex<double>>::get_c_scale)
    .def("get_u_grid", &Complex_prop<complex<double>>::get_u_grid)
    .def("get_a_grid", &Complex_prop<complex<double>>::get_a_grid)
    .def("get_b_grid", &Complex_prop<complex<double>>::get_b_grid)
    .def("get_uv_cutoff", &Complex_prop<complex<double>>::get_uv_cutoff)
    .def("get_z_2", &Complex_prop<complex<double>>::get_z_2)
    .def("get_m_bare", &Complex_prop<complex<double>>::get_m_bare)
    .def("get_is_linear_ren_cond",
      &Complex_prop<complex<double>>::get_is_linear_ren_cond)
    .def("get_n_grid_gc", &Complex_prop<complex<double>>::get_n_grid_gc)
    .def("get_is_pv_mt", &Complex_prop<complex<double>>::get_is_pv_mt)
    .def("get_lambda_pv", &Complex_prop<complex<double>>::get_lambda_pv)
    // pickle suport for python multiprocessing
    .def(py::pickle(
      []( Complex_prop<complex<double>> &p )
      {
        return py::make_tuple( p.get_omega(), p.get_d_scale(), p.get_n_grid(),
        p.get_c_scale(), p.get_u_grid(), p.get_a_grid(), p.get_b_grid(),
        p.get_uv_cutoff(), p.get_z_2(), p.get_m_bare(), p.get_is_linear_ren_cond(),
        p.get_n_grid_gc(), p.get_is_pv_mt(), p.get_lambda_pv() );
      },
      []( py::tuple t )
      {
        if( t.size() != 14 )
        { throw runtime_error("Invalid pickle state of Complex_prop<complex<double>"); }
        Complex_prop<complex<double>> p( t[0].cast<double>(), t[1].cast<double>(),
          t[2].cast<unsigned>(), t[3].cast<double>(), t[4].cast<vector<double>>(),
          t[5].cast<vector<double>>(), t[6].cast<vector<double>>(),
          t[7].cast<double>(), t[8].cast<double>(), t[9].cast<double>(),
          t[10].cast<bool>(), t[11].cast<unsigned>(), t[12].cast<bool>(),
          t[13].cast<double>() );
        return p;
      }
    ));
}

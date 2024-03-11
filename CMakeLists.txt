#[[
"CMakeLists.txt" contains the instructions for CMake to build C++
executables and pybind11 modules for the SDE of the quark propagator.

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
]]
cmake_minimum_required(VERSION 3.9...3.22)
project(SDE_MT)# project name
# choose if pybind11 modules are to be built
set(is_pybind TRUE)
# C++11 or higher versions
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
# dependencies
find_package(GSL REQUIRED)# GNU Scientific Library for natural cubic splines
find_package(OpenMP REQUIRED)# OpenMP for multiprocessing support in C++
message(OpenMP_CXX_version_date = "${OpenMP_CXX_SPEC_DATE}")
message(OpenMP_CXX_version = "${OpenMP_CXX_VERSION}")
message(GSL_location = "${GSL_LIBRARY}")
message(GSL_version = "${GSL_VERSION}")
# build C++ executables
add_executable(main_sde_pmts main_sde_pmts_io.cpp)
add_executable(main_sde_iter main_sde_iter.cpp)
target_link_libraries(main_sde_iter PUBLIC OpenMP::OpenMP_CXX)
add_executable(main_complex_prop main_complex_prop.cpp)
target_link_libraries(main_complex_prop PUBLIC GSL::gsl GSL::gslcblas)
# build pybind11 modules
if( ${is_pybind} )
  find_package(pybind11 REQUIRED)
  pybind11_add_module(ker_itg bind_ker.cpp)
  target_link_libraries(ker_itg PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(ker_itg PUBLIC GSL::gsl GSL::gslcblas)
  pybind11_add_module(sde_ker bind_sde_ker.cpp)
  target_link_libraries(sde_ker PUBLIC GSL::gsl GSL::gslcblas)
  pybind11_add_module(sde_iter_solver bind_sde.cpp)
  target_link_libraries(sde_iter_solver PUBLIC OpenMP::OpenMP_CXX)
endif()

/*
"main_complex_prop.cpp" is the source for the main function that applies
the fixed-grid algorithm to compute the quark propagator with complex-valued
momentum from the SDE and its solution in the Euclidean space.

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
#include"sde_pmts_io.h"// class SDE_PMTS_IO
using namespace std;
// main function as an example of computing the quark propagator with \
complex-valued momentum from the fixed-grid algorithm
int main( int args, char* argv[] )
{
  // load model parameters from the binary file.
  string pmts_filename = SDE_PMTS_FILENAME;// name of the binary file for parameters
  SDE_PMTS_IO sde_pmts( pmts_filename );
  sde_pmts.print_pmts();
  // initialization of the complex propagator by loading the Euclidea-space solution
  Complex_prop<complex<double>> cp( sde_pmts.omega, sde_pmts.d_scale,
    sde_pmts.n_grid_gc, sde_pmts.sf_filename );
  // Test the inverse of the quark propagator with a given complex-valued momentum.
  complex<double> p2( 1.0, 0.25 );
  vector<complex<double>> inv_prop( 2, 0.0 );
  inv_prop = cp.inv_prop_complex( p2 );
  cout<<p2<<endl;// standard value (1.0,0.25)
  cout<<inv_prop[0]<<endl;// standard value (1.41704,-0.0417153)
  cout<<inv_prop[1]<<endl;// standard value (0.171216,-0.0575147)
  return 0;
}

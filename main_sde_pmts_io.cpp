/*
SDE_quark solves the Schwinger-Dyson equation for the quark propagator
in the rainbow-ladder truncation with the application of the Maris-Tandy
model. It also provide fixed-grid routines to compute the quark propagator
with complex-valued momentum.

"main_sde_pmts_io.cpp" is the source of the main function that sets the
model paramters, which is saved to a binary file named by the string
variable sf_filename.

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
#include<cstdlib>// size_t
#include<array>// array<>
#include<cmath>// abs(), pow(), exp(), log()
#include<string>// string
#include<iostream>// cout
#include<fstream>// ifstream, ofstream
#include"consts.h"// constants
#include"sde_pmts_io.h"// class SDE_PMTS_IO
using namespace std;
// The function that create an instance of class SDE_PMTS_IO to be saved in the \
binary file named by pmts_filename.
void init_sde_pmts_io( string pmts_filename )
{
  double omega = 0.4;// momentum scale of the IR term in the Maris-Tandy model in GeV
  double d_scale = 0.859;// strength of the IR term in (GeV)^2
  double c_scale_sde = 10.0;// variable transformation scale in (GeV)^2
  size_t n_u_grid = 2001;// number of grid points in the Euclidean-space solution
  unsigned n_grid_gc = 201;// Number of grid points for the fixed-grid \
  Gauss-Chebyshev quadrature for the angular integral, suggested value is 201.
  string sf_filename = SF_LIGHT_FILENAME;// file name for solution of the iterative solver
  double lambda_pv = 1.0e+4;// Pauli-Villars regularization mass in GeV when applied.
  SDE_PMTS_IO sde_pmts_io( omega, d_scale, c_scale_sde, n_u_grid,
    n_grid_gc, sf_filename, lambda_pv );
  sde_pmts_io.write_pmts_bin( pmts_filename );// Write the binary file.
  cout<<"written parameters for the BSE are"<<endl;
  sde_pmts_io.print_pmts();// Display parameters in sde_pmts_io.
}
// The main function that writes and displays parameters in SDE_PMTS_IO.
int main( int argc, char* argv[] )
{
  string pmts_filename = SDE_PMTS_FILENAME;// Name of the binary file, with \
  default "sde_pmts.bin".
  init_sde_pmts_io( pmts_filename );
  return 0;
}
